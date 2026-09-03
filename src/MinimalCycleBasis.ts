// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MinimalCycleBasis.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Extract the minimal cycle basis for a planar graph. The input vertices and
// edges must form a graph for which edges intersect only at vertices; that
// is, no two edges must intersect at an interior point of one of the edges.
// The algorithm is described in
//   https://www.geometrictools.com/Documentation/MinimalCycleBasis.pdf
// The graph might have isolated vertices (no adjacent vertices via edges).
// These are extracted by the implementation. The graph might have filaments,
// which are subgraphs of polylines that are not shared by a cycle. These are
// also extracted by the implementation.
//
// Port notes:
// * The C++ class template is MinimalCycleBasis<T, IndexType>, where T is a
//   floating-point type and IndexType is an integer type. Both become
//   'number' per PORTING.md, so the port is a single non-generic class.
// * Extract() becomes extract() and the three accessors become
//   getIsolatedVertices(), getFilaments() and getForest().
// * Convexity decisions at a vertex must be exact, so upstream computes the
//   2x2 determinant sign with BSNumber<UIntegerFP32<numWords>>. The port uses
//   the ported BSNumber, which is backed by bigint and therefore needs no
//   worst-case word count; the concrete-instantiation precedent for the
//   arbitrary-precision types is followed (no generic numeric parameter).
// * Upstream stores the adjacency of a vertex in 'std::set<Vertex*>', whose
//   iteration order is the (implementation-defined) address order of the
//   heap-allocated vertices. The port gives each vertex a unique
//   creation-order id and iterates the adjacency in increasing id order, so
//   the traversal is deterministic (the B100+ sorted-iteration precedent).
//   The 'std::set<size_t> detachments', whose increasing-index iteration
//   order upstream does rely on, is likewise iterated in increasing order.
//   'std::map<Vertex*, size_t> duplicates' is only queried by key, so it
//   becomes a Map keyed by vertex id.
// * The 'Vertex::position' pointer into the caller's array becomes a
//   reference to the caller's position tuple (never mutated).
// * The depth-first search 'visited' flags are per-search here rather than
//   persistent members; see the comment on depthFirstSearch for the upstream
//   bug this fixes.

import { BSNumber } from './BSNumber.js';
import { logAssert } from './Logger.js';

// A planar-graph vertex position (x, y).
export type MinimalCycleBasisPosition = readonly [number, number];

// An undirected edge, a pair of indices into the positions array.
export type MinimalCycleBasisEdge = readonly [number, number];

// A node of the cycle forest. The cycle is a closed polyline whose first and
// last indices are the same vertex. The children are the cycles nested
// inside this one.
export interface MinimalCycleBasisTree {
    cycle: number[];
    children: MinimalCycleBasisTree[];
}

// A polyline subgraph that is not shared by a cycle.
export type MinimalCycleBasisFilament = number[];

export type MinimalCycleBasisForest = MinimalCycleBasisTree[];

// The exact rational representation of a position.
type RPosition = readonly [BSNumber, BSNumber];

function rSub(rInput0: RPosition, rInput1: RPosition): RPosition {
    return [rInput0[0].sub(rInput1[0]), rInput0[1].sub(rInput1[1])];
}

function rGetSignDet(rInput0: RPosition, rInput1: RPosition): number {
    const rDet = rInput0[0].mul(rInput1[1]).sub(rInput0[1].mul(rInput1[0]));
    return rDet.getSign();
}

// Lexicographic comparison of positions, the port of the std::array
// 'operator<' that upstream uses to locate the left-most vertex.
function positionLessThan(p0: MinimalCycleBasisPosition,
    p1: MinimalCycleBasisPosition): boolean {
    if (p0[0] < p1[0]) {
        return true;
    }
    if (p0[0] > p1[0]) {
        return false;
    }
    return p0[1] < p1[1];
}

let mcbNextVertexId = 0;

class MCBVertex {
    // A unique identifier that replaces the pointer identity on which the
    // upstream std::set<Vertex*> and std::map<Vertex*, size_t> rely.
    readonly id: number;

    // The index into the 'positions' provided to the call to extract().
    // Vertices cloned during traversal share the index of their original.
    index: number;

    // The position of the vertex, stored as a floating-point tuple.
    position: MinimalCycleBasisPosition;

    // The vertices adjacent to this vertex in the vertex-edge graph, keyed
    // by vertex id.
    adjacents: Map<number, MCBVertex>;

    // The position of the vertex, stored as a rational tuple. Use
    // getRPosition() for memoized conversion of the floating-point tuple to
    // a rational tuple.
    private rPosition: RPosition | null;

    constructor(index: number, position: MinimalCycleBasisPosition) {
        this.id = mcbNextVertexId++;
        this.index = index;
        this.position = position;
        this.adjacents = new Map<number, MCBVertex>();
        this.rPosition = null;
    }

    getRPosition(): RPosition {
        if (this.rPosition === null) {
            this.rPosition = [
                BSNumber.fromNumber(this.position[0]),
                BSNumber.fromNumber(this.position[1])
            ];
        }
        return this.rPosition;
    }

    // The adjacency in increasing-id order, the deterministic stand-in for
    // iterating the upstream std::set<Vertex*>.
    sortedAdjacents(): MCBVertex[] {
        const adjacents = Array.from(this.adjacents.values());
        adjacents.sort((v0, v1) => v0.id - v1.id);
        return adjacents;
    }

    // The port of '*adjacents.begin()' for the cases where upstream has
    // already established that there is exactly one adjacent vertex.
    firstAdjacent(): MCBVertex {
        return this.sortedAdjacents()[0] as MCBVertex;
    }
}

export class MinimalCycleBasis {
    // Storage for vertices of the original graph and for new vertices added
    // during graph traversal.
    private mVertices: MCBVertex[];

    // The output of the extract() function call.
    private mIsolatedVertices: number[];
    private mFilaments: MinimalCycleBasisFilament[];
    private mForest: MinimalCycleBasisForest;

    constructor() {
        this.mVertices = [];
        this.mIsolatedVertices = [];
        this.mFilaments = [];
        this.mForest = [];
    }

    // Extract the cycles, filaments, and isolated vertices. The input
    // positions and edges must form a planar graph for which edges intersect
    // only at vertices; that is, no two edges must intersect at an interior
    // point of one of the edges.
    extract(positions: readonly MinimalCycleBasisPosition[],
        edges: readonly MinimalCycleBasisEdge[],
        verifyInputs: boolean): void {
        this.mVertices = [];
        this.mIsolatedVertices = [];
        this.mFilaments = [];
        this.mForest = [];

        if (positions.length === 0 || edges.length === 0) {
            // The graph is empty, so there are no filaments or cycles.
            return;
        }

        if (verifyInputs) {
            MinimalCycleBasis.verifyInputs(positions, edges);
        }

        this.createGraph(positions, edges);

        // Extract the isolated vertices of the graph.
        this.extractIsolatedVertices();

        // Extract the tree of cycles of the graph. The filaments are
        // extracted during graph traversal.
        this.extractForest();
    }

    // Access to the output of the extract() function.
    getIsolatedVertices(): number[] {
        return this.mIsolatedVertices;
    }

    getFilaments(): MinimalCycleBasisFilament[] {
        return this.mFilaments;
    }

    getForest(): MinimalCycleBasisForest {
        return this.mForest;
    }

    private static verifyInputs(positions: readonly MinimalCycleBasisPosition[],
        edges: readonly MinimalCycleBasisEdge[]): void {
        const uniquePositions = new Set<string>();
        for (const position of positions) {
            uniquePositions.add(position[0] + ',' + position[1]);
        }
        logAssert(uniquePositions.size === positions.length,
            'Input positions must be unique.');

        const numPositions = positions.length;
        for (let i = 0; i < edges.length; ++i) {
            const edge = edges[i] as MinimalCycleBasisEdge;
            logAssert(0 <= edge[0] && edge[0] < numPositions,
                'Input index edge[' + i + '][0] is out of range.');
            logAssert(0 <= edge[1] && edge[1] < numPositions,
                'Input index edge[' + i + '][1] is out of range.');
            logAssert(edge[0] !== edge[1],
                'Input edge[' + i + '] is degenerate.');
        }
    }

    // Create the vertex-edge graph. The edges are undirected.
    private createGraph(positions: readonly MinimalCycleBasisPosition[],
        edges: readonly MinimalCycleBasisEdge[]): void {
        this.mVertices = new Array<MCBVertex>(positions.length);
        for (let index = 0; index < positions.length; ++index) {
            this.mVertices[index] = new MCBVertex(index,
                positions[index] as MinimalCycleBasisPosition);
        }

        for (const edge of edges) {
            MinimalCycleBasis.insert(
                this.mVertices[edge[0]] as MCBVertex,
                this.mVertices[edge[1]] as MCBVertex);
        }
    }

    // Insert an edge into the graph.
    private static insert(vertex0: MCBVertex, vertex1: MCBVertex): void {
        vertex0.adjacents.set(vertex1.id, vertex1);
        vertex1.adjacents.set(vertex0.id, vertex0);
    }

    // Remove an edge from the graph.
    private static remove(vertex0: MCBVertex, vertex1: MCBVertex): void {
        vertex0.adjacents.delete(vertex1.id);
        vertex1.adjacents.delete(vertex0.id);
    }

    // Extract the isolated vertices for the vertex-edge graph.
    private extractIsolatedVertices(): void {
        for (const vertex of this.mVertices) {
            if (vertex.adjacents.size === 0) {
                this.mIsolatedVertices.push(vertex.index);
            }
        }
    }

    // Extract the top-level filaments for the vertex-edge graph. The
    // 'component' array is modified in place, as the upstream reference
    // parameter is.
    private extractFilaments(component: MCBVertex[]): void {
        // Locate all filament endpoints, which are vertices with each having
        // exactly one adjacent vertex.
        const endpoints: MCBVertex[] = [];
        for (const vertex of component) {
            if (vertex.adjacents.size === 1) {
                endpoints.push(vertex);
            }
        }
        if (endpoints.length === 0) {
            // The vertex-edge graph has no filaments.
            return;
        }

        // Remove the filaments from the vertex-edge graph. The greedy removal
        // of vertices allows for removing filaments from a subgraph of
        // filaments that has branch points.
        for (const endpoint of endpoints) {
            let vertex0 = endpoint;
            if (vertex0.adjacents.size === 0) {
                // The endpoint was visited during another filament traversal.
                continue;
            }

            // Traverse the filament and remove the vertices.
            const filament: MinimalCycleBasisFilament = [];
            filament.push(vertex0.index);
            while (vertex0.adjacents.size === 1) {
                const vertex1 = vertex0.firstAdjacent();
                filament.push(vertex1.index);
                MinimalCycleBasis.remove(vertex0, vertex1);
                vertex0 = vertex1;
            }

            // The traversal has terminated because the final vertex is either
            // an endpoint (1 adjacent) or a branch point (at least 3
            // adjacents). When it is an endpoint, the removal in the
            // while-loop reduced the adjacent count to 0. When it is a branch
            // point, the removal in the while-loop reduced the adjacent count
            // to at least 2.
            this.mFilaments.push(filament);
        }

        // At this time the component is either empty because it was an open
        // polyline or it has no filaments and at least one cycle. Identify
        // the remaining vertices and copy to the component, which then has
        // fewer vertices than before the call to extractFilaments.
        MinimalCycleBasis.keepAttached(component);
    }

    // Replace 'component' by the subset of its vertices that still have at
    // least one adjacent vertex.
    private static keepAttached(component: MCBVertex[]): void {
        let numRemaining = 0;
        for (const vertex of component) {
            if (vertex.adjacents.size > 0) {
                component[numRemaining++] = vertex;
            }
        }
        component.length = numRemaining;
    }

    // Extract the minimal cycle basis for the vertex-edge graph, stored as a
    // forest of trees.
    private extractForest(): void {
        const components: MCBVertex[][] = [];
        this.extractConnectedComponents(components);
        for (const component of components) {
            const tree = this.extractBasis(component);
            if (tree.children.length > 0 || tree.cycle.length > 0) {
                this.mForest.push(tree);
            }
        }
    }

    // Extract the connected components of the graph using a depth-first
    // search.
    private extractConnectedComponents(components: MCBVertex[][]): void {
        // The 'visited' set is shared by the searches so that a vertex is
        // assigned to exactly one component. Upstream uses per-vertex
        // 'visited' flags that it resets to zero after component finding.
        const visited = new Set<number>();
        for (const vertex of this.mVertices) {
            if (vertex.adjacents.size >= 2 && !visited.has(vertex.id)) {
                const component: MCBVertex[] = [];
                MinimalCycleBasis.depthFirstSearch(vertex, component, visited);
                components.push(component);
            }
        }
    }

    // Upstream bug (MinimalCycleBasis.h, DepthFirstSearch): the 'visited'
    // flags are members of Vertex that are zeroed only at the end of
    // ExtractConnectedComponents. The searches that collect a detached
    // subgraph (in ExtractCycleFromClosedWalk) therefore see the flags left
    // behind by earlier searches, and a subgraph whose vertices were already
    // visited collapses to the single clone vertex. The nesting of the cycle
    // forest is then wrong: a doubly nested cycle is reported as a sibling of
    // the cycle that contains it instead of as its child. The port fixes this
    // by scoping the visited set to a search (or, for the component finding
    // above, to the loop over components), which reproduces upstream whenever
    // the flags happen to be zero.
    private static depthFirstSearch(vInitial: MCBVertex,
        component: MCBVertex[], visited: Set<number>): void {
        const vStack: MCBVertex[] = [vInitial];
        while (vStack.length > 0) {
            const vertex = vStack[vStack.length - 1] as MCBVertex;
            visited.add(vertex.id);
            const adjacents = vertex.sortedAdjacents();
            let i = 0;
            for (; i < adjacents.length; ++i) {
                const adjacent = adjacents[i] as MCBVertex;
                if (!visited.has(adjacent.id)) {
                    vStack.push(adjacent);
                    break;
                }
            }

            if (i === adjacents.length) {
                component.push(vertex);
                vStack.pop();
            }
        }
    }

    // Extract the minimal cycle basis for a connected component.
    private extractBasis(component: MCBVertex[]): MinimalCycleBasisTree {
        // The top-level tree will not have its cycle member set. The children
        // are the cycle trees extracted from the component.
        const tree: MinimalCycleBasisTree = { cycle: [], children: [] };

        while (component.length > 0) {
            this.extractFilaments(component);
            if (component.length > 0) {
                tree.children.push(this.extractCycleFromComponent(component));
            }
        }

        if (tree.cycle.length === 0 && tree.children.length === 1) {
            // Replace the parent by the child to avoid having two empty
            // cycles in parent/child.
            const child = tree.children[tree.children.length - 1] as MinimalCycleBasisTree;
            tree.cycle = child.cycle;
            tree.children = child.children;
        }

        return tree;
    }

    private extractCycleFromComponent(component: MCBVertex[]): MinimalCycleBasisTree {
        // Search for the left-most vertex of the component. If two or more
        // vertices attain minimum x-value, select the one that has minimum
        // y-value.
        let minVertex = component[0] as MCBVertex;
        for (const vertex of component) {
            if (positionLessThan(vertex.position, minVertex.position)) {
                minVertex = vertex;
            }
        }

        // Traverse the closed walk, duplicating the starting vertex as the
        // last vertex.
        const closedWalk: MCBVertex[] = [];
        let vCurr = minVertex;
        const vStart = vCurr;
        closedWalk.push(vStart);
        let vAdj = MinimalCycleBasis.getClockwiseMost(null, vStart) as MCBVertex;
        while (vAdj !== vStart) {
            closedWalk.push(vAdj);
            const vNext = MinimalCycleBasis.getCounterclockwiseMost(vCurr, vAdj) as MCBVertex;
            vCurr = vAdj;
            vAdj = vNext;
        }
        closedWalk.push(vStart);

        // Recursively process the closed walk to extract cycles.
        const tree = this.extractCycleFromClosedWalk(closedWalk);

        // The isolated vertices generated by cycle removal are also removed
        // from the component.
        MinimalCycleBasis.keepAttached(component);

        return tree;
    }

    private extractCycleFromClosedWalk(closedWalk: MCBVertex[]): MinimalCycleBasisTree {
        const tree: MinimalCycleBasisTree = { cycle: [], children: [] };

        const duplicates = new Map<number, number>();
        const detachments = new Set<number>();
        let numClosedWalk = closedWalk.length;
        for (let i = 1; i < numClosedWalk - 1; ++i) {
            const walkVertex = closedWalk[i] as MCBVertex;
            const diter = duplicates.get(walkVertex.id);
            if (diter === undefined) {
                // We have not yet visited this vertex.
                duplicates.set(walkVertex.id, i);
                continue;
            }

            // The vertex has been visited previously. Collapse the closed
            // walk by removing the subwalk sharing this vertex. Note that the
            // vertex is closedWalk[diter] and closedWalk[i].
            const iMin = diter;
            const iMax = i;
            detachments.add(iMin);
            for (let j = iMin + 1; j < iMax; ++j) {
                const vertex = closedWalk[j] as MCBVertex;
                duplicates.delete(vertex.id);
                detachments.delete(j);
            }
            closedWalk.splice(iMin + 1, iMax - iMin);
            numClosedWalk = closedWalk.length;
            i = iMin;
        }

        if (numClosedWalk > 3) {
            // It is not known whether closedWalk[0] is a detachment point. To
            // determine this, test for any edges strictly contained in the
            // wedge formed by the edges <closedWalk[0],closedWalk[N-1]> and
            // <closedWalk[0],closedWalk[1]>. However, this test must be
            // executed even for the known detachment points. The ensuing
            // logic is designed to handle this and reduce the amount of code,
            // so insert closedWalk[0] into the detachment set and ignore it
            // later if it actually is not.
            detachments.add(0);

            // Detach subgraphs from the vertices of the cycle. The upstream
            // container is a std::set<size_t>, so the indices are visited in
            // increasing order.
            const sortedDetachments = Array.from(detachments);
            sortedDetachments.sort((a, b) => a - b);
            for (const i of sortedDetachments) {
                const orgVertex = closedWalk[i] as MCBVertex;
                const maxVertex = closedWalk[i + 1] as MCBVertex;
                const minVertex = (i > 0
                    ? closedWalk[i - 1] as MCBVertex
                    : closedWalk[numClosedWalk - 2] as MCBVertex);

                const rOrgPos = orgVertex.getRPosition();
                const rDMax = rSub(maxVertex.getRPosition(), rOrgPos);
                const rDMin = rSub(minVertex.getRPosition(), rOrgPos);

                const isConvex = (rGetSignDet(rDMax, rDMin) >= 0);
                const inWedge: MCBVertex[] = [];
                for (const vertex of orgVertex.sortedAdjacents()) {
                    if (vertex.index === minVertex.index ||
                        vertex.index === maxVertex.index) {
                        continue;
                    }

                    const rDVer = rSub(vertex.getRPosition(), rOrgPos);
                    const signDet0 = rGetSignDet(rDVer, rDMin);
                    const signDet1 = rGetSignDet(rDVer, rDMax);
                    const containsVertex = (isConvex
                        ? (signDet0 > 0 && signDet1 < 0)
                        : (signDet0 > 0 || signDet1 < 0));

                    if (containsVertex) {
                        inWedge.push(vertex);
                    }
                }

                if (inWedge.length > 0) {
                    // The clone will manage the adjacents for orgVertex that
                    // lie inside the wedge defined by the first and last
                    // edges of the subgraph rooted at orgVertex. The sorting
                    // is in the clockwise direction.
                    const clone = new MCBVertex(orgVertex.index, orgVertex.position);
                    this.mVertices.push(clone);

                    // Detach the edges inside the wedge.
                    for (const vertex of inWedge) {
                        MinimalCycleBasis.remove(vertex, orgVertex);
                        MinimalCycleBasis.insert(vertex, clone);
                    }

                    // Get the subgraph (it is a single connected component).
                    const component: MCBVertex[] = [];
                    MinimalCycleBasis.depthFirstSearch(clone, component,
                        new Set<number>());

                    // Extract the cycles of the subgraph.
                    tree.children.push(this.extractBasis(component));
                }
                // else the candidate was closedWalk[0] and it has no subgraph
                // to detach.
            }

            tree.cycle = MinimalCycleBasis.extractCycle(closedWalk);
        }
        else {
            // Detach the subgraph from vertex closedWalk[0]; the subgraph is
            // attached via a filament.
            const current = closedWalk[0] as MCBVertex;
            const next = closedWalk[1] as MCBVertex;

            const clone = new MCBVertex(current.index, current.position);
            this.mVertices.push(clone);

            MinimalCycleBasis.remove(next, current);
            MinimalCycleBasis.insert(next, clone);

            // Get the subgraph (it is a single connected component).
            const component: MCBVertex[] = [];
            MinimalCycleBasis.depthFirstSearch(clone, component,
                new Set<number>());

            // Extract the cycles of the subgraph.
            tree.children.push(this.extractBasis(component));
            if (tree.cycle.length === 0 && tree.children.length === 1) {
                // Replace the parent by the child to avoid having two empty
                // cycles in parent/child.
                const child = tree.children[tree.children.length - 1] as MinimalCycleBasisTree;
                tree.cycle = child.cycle;
                tree.children = child.children;
            }
        }

        return tree;
    }

    private static extractCycle(closedWalk: MCBVertex[]): number[] {
        // The logic of this function was designed not to remove filaments
        // after the cycle deletion is complete. This is an iterative process
        // that removes polylines that occur after a cycle has been removed,
        // causing part or all of a cycle boundary to appear to be a filament
        // for the modified graph.

        // The closed walk is a cycle.
        const cycle = new Array<number>(closedWalk.length);
        for (let i = 0; i < closedWalk.length; ++i) {
            cycle[i] = (closedWalk[i] as MCBVertex).index;
        }

        // The clockwise-most edge is always removable.
        let vertex0 = closedWalk[0] as MCBVertex;
        let vertex1 = closedWalk[1] as MCBVertex;
        let vBranch: MCBVertex | null = (vertex0.adjacents.size > 2 ? vertex0 : null);
        MinimalCycleBasis.remove(vertex0, vertex1);

        // Remove edges while traversing counterclockwise.
        while (vertex1 !== vBranch && vertex1.adjacents.size === 1) {
            const adj = vertex1.firstAdjacent();
            MinimalCycleBasis.remove(adj, vertex1);
            vertex1 = adj;
        }

        if (vertex1 !== vertex0) {
            // If vertex1 had exactly 3 adjacent vertices, removal of the CCW
            // edge that shared vertex1 leads to vertex1 having 2 adjacent
            // vertices. When the CW removal occurs and we reach vertex1, the
            // edge deletion will lead to vertex1 having 1 adjacent vertex,
            // making it a filament endpoint. We must ensure we do not delete
            // vertex1 in this case, allowing the recursive algorithm to
            // handle the filament later.
            vBranch = vertex1;

            // Remove edges while traversing clockwise.
            while (vertex0 !== vBranch && vertex0.adjacents.size === 1) {
                vertex1 = vertex0.firstAdjacent();
                MinimalCycleBasis.remove(vertex0, vertex1);
                vertex0 = vertex1;
            }
        }
        // else the cycle is its own connected component.

        return cycle;
    }

    private static getClockwiseMost(vPrev: MCBVertex | null,
        vCurr: MCBVertex): MCBVertex | null {
        let vNext: MCBVertex | null = null;
        let vCurrConvex = false;
        let rDCurr: RPosition;
        let rDNext: RPosition = [BSNumber.fromNumber(0), BSNumber.fromNumber(0)];
        if (vPrev !== null) {
            rDCurr = rSub(vCurr.getRPosition(), vPrev.getRPosition());
        }
        else {
            rDCurr = [BSNumber.fromNumber(0), BSNumber.fromNumber(-1)];
        }

        for (const vAdj of vCurr.sortedAdjacents()) {
            // vAdj is a vertex adjacent to vCurr. No backtracking is allowed.
            if (vAdj === vPrev) {
                continue;
            }

            // Compute the potential direction to move in.
            const rDAdj = rSub(vAdj.getRPosition(), vCurr.getRPosition());

            // Select the first candidate.
            if (vNext === null) {
                vNext = vAdj;
                rDNext = rDAdj;
                vCurrConvex = (rGetSignDet(rDNext, rDCurr) <= 0);
                continue;
            }

            // Update if the next candidate is clockwise of the current
            // clockwise-most vertex.
            const signDet0 = rGetSignDet(rDCurr, rDAdj);
            const signDet1 = rGetSignDet(rDNext, rDAdj);
            if (vCurrConvex) {
                if (signDet0 < 0 || signDet1 < 0) {
                    vNext = vAdj;
                    rDNext = rDAdj;
                    vCurrConvex = (rGetSignDet(rDNext, rDCurr) <= 0);
                }
            }
            else {
                if (signDet0 < 0 && signDet1 < 0) {
                    vNext = vAdj;
                    rDNext = rDAdj;
                    // Upstream bug suspect: this is the only one of the six
                    // convexity updates in the header that tests '< 0'
                    // instead of '<= 0'. The published pseudocode
                    // (MinimalCycleBasis.pdf) uses '<= 0' throughout. The
                    // port preserves the upstream test, which differs only
                    // when the incoming and outgoing directions are
                    // parallel.
                    vCurrConvex = (rGetSignDet(rDNext, rDCurr) < 0);
                }
            }
        }

        return vNext;
    }

    private static getCounterclockwiseMost(vPrev: MCBVertex | null,
        vCurr: MCBVertex): MCBVertex | null {
        let vNext: MCBVertex | null = null;
        let vCurrConvex = false;
        let rDCurr: RPosition;
        let rDNext: RPosition = [BSNumber.fromNumber(0), BSNumber.fromNumber(0)];
        if (vPrev !== null) {
            rDCurr = rSub(vCurr.getRPosition(), vPrev.getRPosition());
        }
        else {
            rDCurr = [BSNumber.fromNumber(0), BSNumber.fromNumber(-1)];
        }

        for (const vAdj of vCurr.sortedAdjacents()) {
            // vAdj is a vertex adjacent to vCurr. No backtracking is allowed.
            if (vAdj === vPrev) {
                continue;
            }

            // Compute the potential direction to move in.
            const rDAdj = rSub(vAdj.getRPosition(), vCurr.getRPosition());

            // Select the first candidate.
            if (vNext === null) {
                vNext = vAdj;
                rDNext = rDAdj;
                vCurrConvex = (rGetSignDet(rDNext, rDCurr) <= 0);
                continue;
            }

            // Select the next candidate if it is counterclockwise of the
            // current counterclockwise-most vertex.
            const signDet0 = rGetSignDet(rDCurr, rDAdj);
            const signDet1 = rGetSignDet(rDNext, rDAdj);
            if (vCurrConvex) {
                if (signDet0 > 0 && signDet1 > 0) {
                    vNext = vAdj;
                    rDNext = rDAdj;
                    vCurrConvex = (rGetSignDet(rDNext, rDCurr) <= 0);
                }
            }
            else {
                if (signDet0 > 0 || signDet1 > 0) {
                    vNext = vAdj;
                    rDNext = rDAdj;
                    vCurrConvex = (rGetSignDet(rDNext, rDCurr) <= 0);
                }
            }
        }

        return vNext;
    }
}
