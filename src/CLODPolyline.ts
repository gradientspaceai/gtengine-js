// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) CLODPolyline.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The continuous level of detail (CLOD) algorithm implemented here is
// described in
// https://www.geometrictools.com/Documentation/PolylineReduction.pdf
// It is an algorithm to reduce incrementally the number of vertices in a
// polyline (open or closed). The sequence of vertex collapses is based on
// vertex weights associated with distance from vertices to polyline segments.
//
// Port notes: upstream 'template <int32_t N, typename Real>' becomes a
// runtime dimension carried by the input vertices (Vector objects). The
// private helper class VertexCollapse, whose 'operator()' writes several
// output reference parameters, becomes module-private functions; the
// collapse driver returns { indices, numEdges, edges } and permutes the
// vertex array in place, as upstream does.
//
// The vertices stored by the object are a reordered copy of the input: the
// vertex-collapse order determines a permutation, and the object's
// getVertices() returns the permuted array. Setting a level of detail L means
// the active vertices are the first L entries of that array and the active
// edges are the first 2*getNumEdges() entries of getEdges().

import { DistPointSegment } from './DistPointSegment';
import { logAssert } from './Logger';
import { MinHeap } from './MinHeap';
import { Segment } from './Segment';
import { Vector, normalize, sub } from './Vector';

export class CLODPolyline {
    // The polyline vertices.
    private mNumVertices: number;
    private mVertices: Vector[];
    private mClosed: boolean;

    // The polyline edges.
    private mNumEdges: number;
    private mEdges: number[];

    // The level of detail information.
    private mVMin: number;
    private mVMax: number;
    private mIndices: number[];

    // Construction. The number of vertices must be 2 or larger (3 or larger
    // when closed). The vertices are assumed to be ordered. The segments are
    // <V[i],V[i+1]> for 0 <= i <= numVertices-2 for an open polyline. If the
    // polyline is closed, an additional segment is <V[numVertices-1],V[0]>.
    constructor(vertices: readonly Vector[], closed: boolean) {
        this.mNumVertices = vertices.length;
        this.mVertices = vertices.map(v => v.clone());
        this.mClosed = closed;
        this.mNumEdges = 0;
        this.mVMin = (closed ? 3 : 2);
        this.mVMax = this.mNumVertices;

        logAssert(closed ? this.mNumVertices >= 3 : this.mNumVertices >= 2,
            'Invalid inputs.');

        // Compute the sequence of vertex collapses. The polyline starts out
        // at the full level of detail (mNumVertices equals mVMax).
        const result = vertexCollapse(this.mVertices, this.mClosed);
        this.mIndices = result.indices;
        this.mNumEdges = result.numEdges;
        this.mEdges = result.edges;
    }

    // Member access.
    getNumVertices(): number {
        return this.mNumVertices;
    }

    // The returned array aliases internal storage.
    getVertices(): readonly Vector[] {
        return this.mVertices;
    }

    getClosed(): boolean {
        return this.mClosed;
    }

    getNumEdges(): number {
        return this.mNumEdges;
    }

    getEdges(): readonly number[] {
        return this.mEdges;
    }

    // Accessors to level of detail (minLOD <= LOD <= maxLOD is required).
    getMinLevelOfDetail(): number {
        return this.mVMin;
    }

    getMaxLevelOfDetail(): number {
        return this.mVMax;
    }

    getLevelOfDetail(): number {
        return this.mNumVertices;
    }

    setLevelOfDetail(numVertices: number): void {
        if (numVertices < this.mVMin || numVertices > this.mVMax) {
            return;
        }

        // Decrease the level of detail.
        while (this.mNumVertices > numVertices) {
            --this.mNumVertices;
            this.mEdges[this.mIndices[this.mNumVertices]] =
                this.mEdges[2 * this.mNumEdges - 1];
            --this.mNumEdges;
        }

        // Increase the level of detail.
        while (this.mNumVertices < numVertices) {
            ++this.mNumEdges;
            this.mEdges[this.mIndices[this.mNumVertices]] = this.mNumVertices;
            ++this.mNumVertices;
        }
    }
}

// Support for vertex collapses in a polyline (upstream's private class
// VertexCollapse). The 'vertices' array is permuted in place.
function vertexCollapse(vertices: Vector[], closed: boolean):
    { indices: number[], numEdges: number, edges: number[] } {
    const numVertices = vertices.length;
    const indices = new Array<number>(numVertices).fill(0);
    let numEdges: number;
    let edges: number[];

    if (closed) {
        numEdges = numVertices;
        edges = new Array<number>(2 * numEdges).fill(0);

        if (numVertices === 3) {
            indices[0] = 0;
            indices[1] = 1;
            indices[2] = 3;
            edges[0] = 0; edges[1] = 1;
            edges[2] = 1; edges[3] = 2;
            edges[4] = 2; edges[5] = 0;
            return { indices, numEdges, edges };
        }
    }
    else {
        numEdges = numVertices - 1;
        edges = new Array<number>(2 * numEdges).fill(0);

        if (numVertices === 2) {
            indices[0] = 0;
            indices[1] = 1;
            edges[0] = 0; edges[1] = 1;
            return { indices, numEdges, edges };
        }
    }

    // Create the heap of weights.
    const heap = new MinHeap<number, number>(numVertices);
    const qm1 = numVertices - 1;
    if (closed) {
        const qm2 = numVertices - 2;
        heap.insert(0, getWeight(qm1, 0, 1, vertices));
        heap.insert(qm1, getWeight(qm2, qm1, 0, vertices));
    }
    else {
        heap.insert(0, Number.MAX_VALUE);
        heap.insert(qm1, Number.MAX_VALUE);
    }
    for (let m = 0, z = 1, p = 2; z < qm1; ++m, ++z, ++p) {
        heap.insert(z, getWeight(m, z, p, vertices));
    }

    // Create the level of detail information for the polyline.
    const collapses = collapseVertices(heap, numVertices);
    computeEdges(numVertices, closed, collapses, indices, numEdges, edges);
    reorderVertices(numVertices, vertices, collapses, numEdges, edges);
    return { indices, numEdges, edges };
}

// Weight calculation for vertex triple <V[m],V[z],V[p]>.
function getWeight(m: number, z: number, p: number,
    vertices: readonly Vector[]): number {
    const direction = sub(vertices[p], vertices[m]);
    const len = normalize(direction);
    if (len > 0) {
        const segment = Segment.fromEndpoints(vertices[m], vertices[p]);
        const query = new DistPointSegment();
        const distance = query.compute(vertices[z], segment).distance;
        return distance / len;
    }
    else {
        return Number.MAX_VALUE;
    }
}

// Create data structures for the polyline. The vertices are removed from the
// min-heap in increasing order of weight, and stored from the end of the
// 'collapses' array backward, so the first vertex to collapse is last in the
// array.
function collapseVertices(heap: MinHeap<number, number>,
    numVertices: number): number[] {
    const collapses = new Array<number>(numVertices).fill(0);
    for (let i = numVertices - 1; i >= 0; --i) {
        const removed = heap.remove();
        collapses[i] = (removed !== null ? removed.key : 0);
    }
    return collapses;
}

function computeEdges(numVertices: number, closed: boolean,
    collapses: readonly number[], indices: number[], numEdges: number,
    edges: number[]): void {
    // Compute the edges (first to collapse is last in array). Do not collapse
    // the last line segment of an open polyline. Do not collapse further when
    // a closed polyline becomes a triangle.
    let i: number;
    let vIndex: number;
    let eIndex = 2 * numEdges - 1;
    if (closed) {
        for (i = numVertices - 1; i >= 0; --i) {
            vIndex = collapses[i];
            edges[eIndex--] = (vIndex + 1) % numVertices;
            edges[eIndex--] = vIndex;
        }
    }
    else {
        for (i = numVertices - 1; i >= 2; --i) {
            vIndex = collapses[i];
            edges[eIndex--] = vIndex + 1;
            edges[eIndex--] = vIndex;
        }

        // The two endpoints of an open polyline have infinite weight, so they
        // are the last two removed from the min-heap and occupy collapses[0]
        // and collapses[1]. The one remaining edge is <v,v+1> for the
        // endpoint v that is not the last vertex.
        //
        // Upstream bug (fixed here): upstream unconditionally uses
        // vIndex = collapses[0]. The min-heap breaks the tie between the two
        // infinite weights so that collapses[0] is the last vertex, hence
        // upstream writes edges[1] = numVertices, an index that is out of
        // range for the vertex array. C++ then reads permute[numVertices] out
        // of bounds in ReorderVertices (and TypeScript would store undefined).
        vIndex = (collapses[0] !== numVertices - 1
            ? collapses[0] : collapses[1]);
        edges[0] = vIndex;
        edges[1] = vIndex + 1;
    }

    // In the given edge order, find the index in the edge array that
    // corresponds to a collapse vertex and save the index for the dynamic
    // change in level of detail. This relies on the assumption that a vertex
    // is shared by at most two edges.
    eIndex = 2 * numEdges - 1;
    for (i = numVertices - 1; i >= 0; --i) {
        vIndex = collapses[i];
        for (let e = 0; e < 2 * numEdges; ++e) {
            if (vIndex === edges[e]) {
                indices[i] = e;
                edges[e] = edges[eIndex];
                break;
            }
        }
        eIndex -= 2;

        if (closed) {
            if (eIndex === 5) {
                break;
            }
        }
        else {
            if (eIndex === 1) {
                break;
            }
        }
    }

    // Restore the edge array to full level of detail.
    if (closed) {
        for (i = 3; i < numVertices; ++i) {
            edges[indices[i]] = collapses[i];
        }
    }
    else {
        for (i = 2; i < numVertices; ++i) {
            edges[indices[i]] = collapses[i];
        }
    }
}

function reorderVertices(numVertices: number, vertices: Vector[],
    collapses: readonly number[], numEdges: number, edges: number[]): void {
    const permute = new Array<number>(numVertices).fill(0);
    const permutedVertex = new Array<Vector>(numVertices);

    for (let i = 0; i < numVertices; ++i) {
        const vIndex = collapses[i];
        permute[vIndex] = i;
        permutedVertex[i] = vertices[vIndex];
    }

    for (let i = 0; i < 2 * numEdges; ++i) {
        edges[i] = permute[edges[i]];
    }

    for (let i = 0; i < numVertices; ++i) {
        vertices[i] = permutedVertex[i];
    }
}
