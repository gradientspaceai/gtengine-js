// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MinimumSpanningTree.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the minimum spanning tree of a vertex-edge graph. The code is an
// implementation of Prim's algorithm based on the pseudocode in
//   Introduction to Algorithms, 4th edition (April 5, 2022)
//   Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, Clifford Stein,
//   The MIT Press, Cambridge, Massachusetts
// The pseudocode uses a priority queue that is sorted based on a set of keys.
// A plain priority queue cannot be told to update when a key is modified
// outside the class. The MinHeap class provides this capability.
//
// The edges[] input to execute(...) must be unique. An edge is the unordered
// pair of its two vertex indices, so (v0,v1) and (v1,v0) are the same edge
// and only one of them may appear in edges[]. The v0 and v1 must be
// different numbers.
//
// The weights[] input must have the same number of elements as edges[]. Also,
// the weights must be nonnegative.
//
// Set validateInputs to true to have execute() test for valid input. This is
// an expensive operation that the caller might not want if it is known the
// inputs are valid.
//
// The output minimumSpanningTree[] is the minimum spanning tree. The first
// element is a sentinel [nil,vertexIndex], where vertexIndex is the index of
// the starting vertex in the spanning tree. (If the graph is not connected,
// each connected component contributes such a sentinel and the output is a
// minimum spanning forest.)
//
// The output backEdges[] are the graph edges not in the minimum spanning
// tree. The tree has no cycles, but if you were to insert a back edge into
// the tree, the resulting graph has a cycle.
//
// Port notes: the upstream WeightType template argument is 'number'. The
// upstream static Execute(...) with two output parameters becomes the static
// execute(...) that returns { minimumSpanningTree, backEdges }. The upstream
// std::size_t sentinel std::numeric_limits<std::size_t>::max() becomes
// MinimumSpanningTree.nil = Number.MAX_SAFE_INTEGER, and maxWeight is
// Number.MAX_VALUE (the counterpart of std::numeric_limits<double>::max()).
// The iteration order of upstream std::map<Edge,WeightType> determines the
// order of backEdges[]; the port replicates it by sorting the remaining edge
// keys lexicographically. The vertex adjacency lists are built in the same
// order as upstream so that ties among equal weights are broken identically.

import { logAssert } from './Logger.js';
import { MinHeap } from './MinHeap.js';
import type { MinHeapRecord } from './MinHeap.js';

// The port of MinimumSpanningTree::Edge, a pair of vertex indices. The name
// is prefixed because 'Edge' is used by several mesh classes of the library.
export type MSTEdge = [number, number];

// The port of the two output parameters of upstream Execute(...).
export interface MinimumSpanningTreeResult {
    minimumSpanningTree: MSTEdge[];
    backEdges: MSTEdge[];
}

export class MinimumSpanningTree {
    // The invalid-index sentinel. Upstream uses the maximum std::size_t.
    static readonly nil: number = Number.MAX_SAFE_INTEGER;
    static readonly zeroWeight: number = 0;
    static readonly maxWeight: number = Number.MAX_VALUE;

    // The vertices are in {0,...,numVertices-1}. The edges[] need not
    // reference all the vertices.
    static execute(
        edges: readonly MSTEdge[],
        weights: readonly number[],
        validateInputs: boolean): MinimumSpanningTreeResult {
        const minimumSpanningTree: MSTEdge[] = [];
        const backEdges: MSTEdge[] = [];

        if (validateInputs) {
            MinimumSpanningTree.validateInputs(edges, weights);
        }

        // Map the distinct vertex indices to consecutive indices from 0 to
        // maxVertices-1. The map key is the vertex index and the map value is
        // its counterpart in the consecutive indices.
        const vertexMap = new Map<number, number>();
        const inverseVertexMap: number[] = [];
        MinimumSpanningTree.createVertexMap(edges, vertexMap, inverseVertexMap);

        // Create an edge map using the remapped vertex indices. At the same
        // time, create a vertex adjacency map.
        const edgeMap = new Map<string, number>();
        const adjacencyMap = new Map<number, number[]>();
        MinimumSpanningTree.createEdgeAndAdjacencyMaps(edges, weights, vertexMap,
            edgeMap, adjacencyMap);

        // Use a priority queue to extract the minimum spanning tree. The
        // vertex indices are the remapped ones.
        MinimumSpanningTree.extractMinimumSpanningTree(vertexMap.size, edgeMap,
            adjacencyMap, minimumSpanningTree);

        // Remove the minimum spanning tree edges from the edge map. The
        // remaining elements are back edges, but include both (v0,v1) and
        // (v1,v0). The duplicates are omitted by storing only those edges for
        // which v0 < v1.
        MinimumSpanningTree.extractBackEdges(edgeMap, minimumSpanningTree, backEdges);

        // Convert back to the original vertex indices.
        MinimumSpanningTree.convertToOriginalIndices(inverseVertexMap,
            minimumSpanningTree, backEdges);

        return { minimumSpanningTree, backEdges };
    }

    // The key of an edge in the edge map. The upstream key is the
    // std::array<std::size_t,2> itself, ordered lexicographically by
    // std::map.
    private static edgeKey(v0: number, v1: number): string {
        return v0 + ',' + v1;
    }

    private static parseEdgeKey(key: string): MSTEdge {
        const comma = key.indexOf(',');
        return [Number(key.substring(0, comma)), Number(key.substring(comma + 1))];
    }

    private static validateInputs(
        edges: readonly MSTEdge[],
        weights: readonly number[]): void {
        logAssert(
            edges.length === weights.length,
            'The edges.size() and weights.size() must match.');

        const uniqueEdges = new Set<string>();
        for (let e = 0; e < edges.length; ++e) {
            logAssert(
                weights[e] >= MinimumSpanningTree.zeroWeight,
                'Encountered a negative weight.');

            const edge = edges[e];
            logAssert(
                edge[0] !== MinimumSpanningTree.nil &&
                edge[1] !== MinimumSpanningTree.nil &&
                edge[0] !== edge[1],
                'Encountered a degenerate edge.');

            if (edge[0] < edge[1]) {
                uniqueEdges.add(MinimumSpanningTree.edgeKey(edge[0], edge[1]));
            } else {
                uniqueEdges.add(MinimumSpanningTree.edgeKey(edge[1], edge[0]));
            }
        }
        logAssert(
            edges.length === uniqueEdges.size,
            'Encountered a duplicate edge.');
    }

    private static createVertexMap(
        edges: readonly MSTEdge[],
        vertexMap: Map<number, number>,
        inverseVertexMap: number[]): void {
        let numVertices = 0;
        for (const edge of edges) {
            for (let i = 0; i < 2; ++i) {
                if (!vertexMap.has(edge[i])) {
                    vertexMap.set(edge[i], numVertices);
                    ++numVertices;
                }
            }
        }

        inverseVertexMap.length = numVertices;
        for (const [key, value] of vertexMap) {
            inverseVertexMap[value] = key;
        }
    }

    private static createEdgeAndAdjacencyMaps(
        edges: readonly MSTEdge[],
        weights: readonly number[],
        vertexMap: ReadonlyMap<number, number>,
        edgeMap: Map<string, number>,
        adjacencyMap: Map<number, number[]>): void {
        for (let e = 0; e < edges.length; ++e) {
            const edge = edges[e];
            const weight = weights[e];

            const v0 = vertexMap.get(edge[0]) as number;
            const v1 = vertexMap.get(edge[1]) as number;
            const remapped: MSTEdge[] = [[v0, v1], [v1, v0]];

            // The logAssert calls are required if validateInputs is false in
            // the call to execute(...).
            for (let i = 0; i < 2; ++i) {
                const key = MinimumSpanningTree.edgeKey(remapped[i][0], remapped[i][1]);
                logAssert(
                    !edgeMap.has(key),
                    'Unexpected result for validated edges.');
                edgeMap.set(key, weight);

                const adjacents = adjacencyMap.get(remapped[i][0]);
                if (adjacents !== undefined) {
                    adjacents.push(remapped[i][1]);
                } else {
                    adjacencyMap.set(remapped[i][0], [remapped[i][1]]);
                }
            }
        }
    }

    private static extractMinimumSpanningTree(
        numVertices: number,
        edgeMap: ReadonlyMap<string, number>,
        adjacencyMap: ReadonlyMap<number, number[]>,
        minimumSpanningTree: MSTEdge[]): void {
        if (numVertices === 0) {
            // There are no edges, so there is no spanning tree. Upstream
            // writes records[0] into an empty std::vector in this case, which
            // is undefined behavior; see the port notes in the PR.
            return;
        }

        // Initialize the priority queue.
        const heap = new MinHeap<MSTEdge, number>(numVertices);
        const records: (MinHeapRecord<MSTEdge, number> | null)[] = new Array(numVertices);
        records[0] = heap.insert([MinimumSpanningTree.nil, 0],
            MinimumSpanningTree.zeroWeight);
        for (let i = 1; i < numVertices; ++i) {
            records[i] = heap.insert([MinimumSpanningTree.nil, i],
                MinimumSpanningTree.maxWeight);
        }

        const inHeap: number[] = new Array<number>(numVertices).fill(1);
        while (heap.getNumElements() > 0) {
            const removed = heap.remove();
            logAssert(removed !== null, 'Unexpected condition.');
            // C++ copies the key into treeEdge; the copy is required here
            // because the heap record continues to own its key array.
            const treeEdge: MSTEdge = [removed.key[0], removed.key[1]];
            inHeap[treeEdge[1]] = 0;
            minimumSpanningTree.push(treeEdge);

            const adjacents = adjacencyMap.get(treeEdge[1]);
            logAssert(adjacents !== undefined, 'Unexpected condition.');

            for (const a of adjacents) {
                if (inHeap[a]) {
                    const candidateWeight =
                        edgeMap.get(MinimumSpanningTree.edgeKey(treeEdge[1], a));
                    logAssert(candidateWeight !== undefined, 'Unexpected condition.');

                    const record = records[a];
                    logAssert(record !== null, 'Unexpected condition.');
                    const currentWeight = record.value;
                    if (candidateWeight < currentWeight) {
                        record.key[0] = treeEdge[1];
                        heap.update(record, candidateWeight);
                    }
                }
            }
        }
    }

    private static extractBackEdges(
        edgeMap: Map<string, number>,
        minimumSpanningTree: readonly MSTEdge[],
        backEdges: MSTEdge[]): void {
        // Remove the tree edges from the graph.
        for (const treeEdge of minimumSpanningTree) {
            edgeMap.delete(MinimumSpanningTree.edgeKey(treeEdge[0], treeEdge[1]));
            edgeMap.delete(MinimumSpanningTree.edgeKey(treeEdge[1], treeEdge[0]));
        }

        // Locate the back edges. They occur in pairs, so eliminate one of the
        // pair using vertex ordering. The upstream std::map iterates in
        // lexicographically increasing key order, which is replicated here by
        // sorting the remaining keys.
        const remaining: MSTEdge[] = [];
        for (const key of edgeMap.keys()) {
            remaining.push(MinimumSpanningTree.parseEdgeKey(key));
        }
        remaining.sort((a, b) => (a[0] !== b[0] ? a[0] - b[0] : a[1] - b[1]));

        for (const edge of remaining) {
            if (edge[0] < edge[1]) {
                backEdges.push(edge);
            }
        }
    }

    private static convertToOriginalIndices(
        inverseVertexMap: readonly number[],
        minimumSpanningTree: MSTEdge[],
        backEdges: MSTEdge[]): void {
        for (const treeEdge of minimumSpanningTree) {
            for (let i = 0; i < 2; ++i) {
                if (treeEdge[i] !== MinimumSpanningTree.nil) {
                    treeEdge[i] = inverseVertexMap[treeEdge[i]];
                }
            }
        }

        for (const backEdge of backEdges) {
            for (let i = 0; i < 2; ++i) {
                backEdge[i] = inverseVertexMap[backEdge[i]];
            }
        }
    }
}
