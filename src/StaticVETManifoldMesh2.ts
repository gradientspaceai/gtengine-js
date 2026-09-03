// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) StaticVETManifoldMesh2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// StaticVETManifoldMesh2 represents a vertex-edge-triangle manifold mesh for
// which the triangles are provided as a single batch and no mesh modification
// operations are performed on the mesh. The class minimizes the memory
// management costs by storing the vertex adjacency information in a single
// block of memory.
//
// Port notes: the upstream 'invalid' sentinel std::numeric_limits<size_t>::max()
// becomes StaticVETManifoldMesh2.invalid = Number.MAX_SAFE_INTEGER, which is
// representable exactly by a float64 and compares greater than any valid
// index. The upstream contiguous std::vector<size_t> storage block, whose
// subblocks are handed to the vertices as raw pointers, becomes a single
// Float64Array with the vertices holding offsets into it; the
// std::array<size_t,3>* edge adjacency records become triples of consecutive
// elements. The upstream index-range tests 'v < mVertices.size()' rely on
// size_t wraparound to reject negative inputs, so the port tests
// '0 <= v && v < numVertices' to obtain the same behavior. The multithreaded
// construction path is omitted (JavaScript has no shared-memory threads
// here); the numThreads argument is accepted for call-site parity and the
// construction always runs the single-threaded path, which computes exactly
// the same result. The out-parameters of GetAdjacentTriangles, GetComponents
// and GetBoundaryPolygons become return values. getAdjacentTriangles also
// fixes an upstream bug; see the comment in that method.

import { logAssert } from './Logger.js';

// The port of the nested class StaticVETManifoldMesh2::Vertex. The name is
// prefixed because the library exports one flat namespace.
export class StaticVETManifoldMesh2Vertex {
    // Only StaticVETManifoldMesh2 writes these members.
    private mStorage: Float64Array;
    private mNumTAdjacents: number;
    private mNumVAdjacents: number;  // <= 2 * mNumTAdjacents
    private mVAdjacents: number;     // offset of [2 * mNumTAdjacents] block
    private mNumEAdjacents: number;  // = mNumTAdjacents after construction
    private mEAdjacents: number;     // offset of [3 * mNumTAdjacents] block

    constructor(storage: Float64Array) {
        this.mStorage = storage;
        this.mNumTAdjacents = 0;
        this.mNumVAdjacents = 0;
        this.mVAdjacents = 0;
        this.mNumEAdjacents = 0;
        this.mEAdjacents = 0;
    }

    // The members are read-only.
    getNumTAdjacents(): number {
        return this.mNumTAdjacents;
    }

    // The number of adjacent vertices is bounded by twice the number of
    // triangles sharing the vertex.
    getNumVAdjacents(): number {
        return this.mNumVAdjacents;
    }

    getVAdjacent(i: number): number {
        return this.mStorage[this.mVAdjacents + i];
    }

    getVAdjacents(): number[] {
        const adjacents: number[] = new Array<number>(this.mNumVAdjacents);
        for (let i = 0; i < this.mNumVAdjacents; ++i) {
            adjacents[i] = this.mStorage[this.mVAdjacents + i];
        }
        return adjacents;
    }

    // The number of adjacent (outgoing) edges is the same as the number of
    // triangles sharing the vertex.
    getNumEAdjacents(): number {
        return this.mNumEAdjacents;
    }

    // The adjacency triple <AV,LT,RT> of the j-th outgoing edge: AV is the
    // other vertex of the edge, LT is the triangle to the left of the edge
    // and RT is the triangle to the right of the edge (RT is 'invalid' when
    // the edge is a boundary edge).
    getEAdjacent(j: number): [number, number, number] {
        const base = this.mEAdjacents + 3 * j;
        return [this.mStorage[base], this.mStorage[base + 1], this.mStorage[base + 2]];
    }

    getEAdjacents(): [number, number, number][] {
        const adjacents: [number, number, number][] = [];
        for (let j = 0; j < this.mNumEAdjacents; ++j) {
            adjacents.push(this.getEAdjacent(j));
        }
        return adjacents;
    }

    // The remaining members are the port of the private upstream members;
    // they are internal to the mesh construction. TypeScript has no 'friend'
    // declaration, so they are marked internal by convention.

    /** @internal */
    initialize(numTAdjacent: number, storageOffset: number): number {
        this.mNumTAdjacents = numTAdjacent;
        this.mNumVAdjacents = 0;
        this.mVAdjacents = storageOffset;
        storageOffset += 2 * this.mNumTAdjacents;
        this.mNumEAdjacents = 0;
        this.mEAdjacents = storageOffset;
        storageOffset += 3 * this.mNumTAdjacents;
        return storageOffset;
    }

    /** @internal */
    insertVAdjacent(v: number): void {
        for (let i = 0; i < this.mNumVAdjacents; ++i) {
            if (v === this.mStorage[this.mVAdjacents + i]) {
                // The vertex v is already in the adjacents list.
                return;
            }
        }

        // The vertex v is not in the adjacents list, so append it.
        this.mStorage[this.mVAdjacents + this.mNumVAdjacents] = v;
        ++this.mNumVAdjacents;
    }

    /** @internal */
    insertEAdjacent(v: number, t: number): void {
        const base = this.mEAdjacents + 3 * this.mNumEAdjacents;
        this.mStorage[base] = v;
        this.mStorage[base + 1] = t;
        this.mStorage[base + 2] = StaticVETManifoldMesh2.invalid;
        ++this.mNumEAdjacents;
    }

    // Return the storage offset of the adjacency triple of the outgoing edge
    // <this,v1>, or -1 when there is no such outgoing edge. This is the port
    // of the pointer returned by StaticVETManifoldMesh2::GetOutgoingEdge.
    /** @internal */
    findEAdjacentOffset(v1: number): number {
        for (let j = 0; j < this.mNumEAdjacents; ++j) {
            const base = this.mEAdjacents + 3 * j;
            if (this.mStorage[base] === v1) {
                return base;
            }
        }
        return -1;
    }
}

// The port of the protected upstream struct BoundaryEdge.
interface VETBoundaryEdge {
    t: number;
    a: number;
    visited: boolean;
}

export class StaticVETManifoldMesh2 {
    // Use the maximum safe integer to denote an invalid index. Upstream uses
    // the maximum size_t, effectively representing -1.
    static readonly invalid: number = Number.MAX_SAFE_INTEGER;

    protected mVertices: StaticVETManifoldMesh2Vertex[];
    protected mStorage: Float64Array;
    protected mTriangles: [number, number, number][];
    protected mAdjacents: [number, number, number][];
    protected mMinTrianglesAtVertex: number;
    protected mMaxTrianglesAtVertex: number;

    // Preconditions.
    //   1. The triangles input must have size 1 or larger.
    //   2. The number of vertices must be 3 or larger.
    //   3. The triangles must form a manifold mesh.
    //   4. Each triangle must be nondegenerate; no repeated vertices.
    //   5. The triangles must all be ordered counterclockwise or all ordered
    //      clockwise; no mixed chirality.
    // The numThreads argument is accepted for parity with upstream but the
    // construction always occurs single threaded.
    constructor(
        numVertices: number,
        triangles: readonly (readonly [number, number, number])[] | readonly (readonly number[])[],
        numThreads: number = 0) {
        void numThreads;
        logAssert(numVertices >= 3 && triangles.length > 0, 'invalid input');

        this.mStorage = new Float64Array(15 * triangles.length)
            .fill(StaticVETManifoldMesh2.invalid);
        this.mVertices = new Array<StaticVETManifoldMesh2Vertex>(numVertices);
        for (let v = 0; v < numVertices; ++v) {
            this.mVertices[v] = new StaticVETManifoldMesh2Vertex(this.mStorage);
        }
        this.mTriangles = triangles.map(
            tri => [tri[0], tri[1], tri[2]] as [number, number, number]);
        this.mAdjacents = triangles.map(() => [
            StaticVETManifoldMesh2.invalid,
            StaticVETManifoldMesh2.invalid,
            StaticVETManifoldMesh2.invalid
        ] as [number, number, number]);
        this.mMinTrianglesAtVertex = 0;
        this.mMaxTrianglesAtVertex = 0;

        const numTrianglesAtVertex: number[] = new Array<number>(numVertices).fill(0);
        this.getNumTrianglesAtVertex(numTrianglesAtVertex);
        this.initializeVertexStorage(numTrianglesAtVertex);
        this.populateVertices();
        this.updateAdjacencyForSharedEdges();
    }

    // Member access.
    getVertices(): StaticVETManifoldMesh2Vertex[] {
        return this.mVertices;
    }

    getTriangles(): [number, number, number][] {
        return this.mTriangles;
    }

    getAdjacents(): [number, number, number][] {
        return this.mAdjacents;
    }

    getMinNumTrianglesAtVertex(): number {
        return this.mMinTrianglesAtVertex;
    }

    getMaxNumTrianglesAtVertex(): number {
        return this.mMaxTrianglesAtVertex;
    }

    // Determine whether or not the undirected edge (v0,v1) exists.
    edgeExists(v0: number, v1: number): boolean {
        const numVertices = this.mVertices.length;
        if (0 <= v0 && v0 < numVertices && 0 <= v1 && v1 < numVertices && v0 !== v1) {
            if (this.mVertices[v0].findEAdjacentOffset(v1) >= 0) {
                return true;
            }

            if (this.mVertices[v1].findEAdjacentOffset(v0) >= 0) {
                return true;
            }
        }
        return false;
    }

    // Get the adjacent triangles for the undirected edge (v0,v1). The
    // returned adjacent triangle indices adj0 and adj1 are the following:
    //
    //   1. <v0,v1> and <v1,v0> are both outgoing edges, so the edge is shared
    //      by two triangles and both adj0 and adj1 are valid (not equal to
    //      'invalid'). The index adj0 is the L-triangle for <v0,v1> and the
    //      index adj1 is the R-triangle for <v0,v1>. Equivalently, adj0 is
    //      the R-triangle for <v1,v0> and adj1 is the L-triangle for <v1,v0>.
    //
    //   2. <v0,v1> is outgoing but <v1,v0> is not outgoing. The index adj0 is
    //      the L-triangle for <v0,v1> and the index adj1 is 'invalid' (no
    //      R-triangle).
    //
    //   3. <v1,v0> is outgoing but <v0,v1> is not outgoing. The index adj0 is
    //      'invalid' (no L-triangle) for <v1,v0> and the index adj1 is the
    //      R-triangle for <v1,v0>.
    //
    //   4. The outgoing edge <v0,v1> does not exist.
    //
    // It is possible to distinguish between the 4 cases by examining the
    // returned indices:
    //   (1) returns (valid, valid) and Boolean 'true'
    //   (2) returns (valid, invalid) and Boolean 'true'
    //   (3) returns (invalid, valid) and Boolean 'true'
    //   (4) returns (invalid, invalid) and Boolean 'false'
    getAdjacentTriangles(v0: number, v1: number):
        { exists: boolean, adj0: number, adj1: number } {
        const numVertices = this.mVertices.length;
        if (0 <= v0 && v0 < numVertices && 0 <= v1 && v1 < numVertices && v0 !== v1) {
            // UPSTREAM BUG (StaticVETManifoldMesh2.h, GetAdjacentTriangles).
            // When the outgoing edge <v0,v1> is not found, upstream searches
            // for the outgoing edge <v1,v0> and then copies that triple's
            // components [1] and [2] straight into (adj0,adj1). Those are the
            // L-triangle and R-triangle of <v1,v0>, which are the R-triangle
            // and the L-triangle of the queried direction <v0,v1>, so the
            // pair must be swapped. Without the swap the documented case 3,
            // "returns (invalid, valid)", can never be produced: a boundary
            // edge whose only outgoing direction is <v1,v0> yields
            // (valid, invalid), which is the signature documented for case 2.
            // The port performs the two directed lookups the documented
            // contract requires: adj0 is the L-triangle of <v0,v1> (the
            // triangle containing the directed edge <v0,v1>) and adj1 is the
            // L-triangle of <v1,v0> (its R-triangle counterpart).
            const offset0 = this.mVertices[v0].findEAdjacentOffset(v1);
            const offset1 = this.mVertices[v1].findEAdjacentOffset(v0);
            if (offset0 >= 0 || offset1 >= 0) {
                return {
                    exists: true,
                    adj0: offset0 >= 0
                        ? this.mStorage[offset0 + 1] : StaticVETManifoldMesh2.invalid,
                    adj1: offset1 >= 0
                        ? this.mStorage[offset1 + 1] : StaticVETManifoldMesh2.invalid
                };
            }
        }

        return {
            exists: false,
            adj0: StaticVETManifoldMesh2.invalid,
            adj1: StaticVETManifoldMesh2.invalid
        };
    }

    // The connected components of the mesh. Each component is an array of
    // triangle indices.
    getComponents(): number[][] {
        const components: number[][] = [];

        // The values are 0 (unvisited), 1 (discovered), 2 (finished).
        const visited: number[] = new Array<number>(this.mTriangles.length).fill(0);

        // Share a stack for the depth-first search. This avoids allocating
        // and deallocating a stack for each call to depthFirstSearch.
        const sharedStack: number[] = new Array<number>(this.mTriangles.length).fill(0);

        // The code reserves maximum space for the component in order to avoid
        // allocation costs associated with resizing caused by push_back.
        const sharedComponents: number[] = new Array<number>(this.mTriangles.length).fill(0);

        for (let t = 0; t < this.mTriangles.length; ++t) {
            if (visited[t] === 0) {
                const numInserted = this.depthFirstSearch(t, visited, sharedStack,
                    sharedComponents);
                components.push(sharedComponents.slice(0, numInserted));
            }
        }

        return components;
    }

    // Compute the boundary-edge components of the mesh. These are polygons
    // that are simple for the strict definition of manifold mesh that
    // disallows bow-tie configurations. The GTE mesh implementations do allow
    // bow-tie configurations, in which case some polygons might not be
    // simple. If you select duplicateEndpoints to be false, a component has
    // consecutive vertices (v[0], v[1], ..., v[n-1]) and the polygon has
    // edges (v[0],v[1]), (v[1],v[2]), ..., (v[n-2],v[n-1]), (v[n-1],v[0]). If
    // duplicateEndpoints is set to true, a component has consecutive vertices
    // (v[0], v[1], ..., v[n-1], v[0]), emphasizing that the component is
    // closed.
    getBoundaryPolygons(duplicateEndpoints: boolean): number[][] {
        const polygons: number[][] = [];

        // Get the boundary edges. The upstream container is a
        // std::map<std::array<size_t,2>, BoundaryEdge>; the port stores the
        // directed edge key as a string and iterates the keys in the
        // lexicographic order of the (v0,v1) pairs that std::map uses.
        const boundaryEdges = new Map<string, VETBoundaryEdge>();
        const directedKeys: [number, number][] = [];
        for (let t = 0; t < this.mTriangles.length; ++t) {
            const tri = this.mTriangles[t];
            for (let a = 0; a < 3; ++a) {
                if (this.mAdjacents[t][a] === StaticVETManifoldMesh2.invalid) {
                    const directed: [number, number] = [tri[a], tri[(a + 1) % 3]];
                    const key = directed[0] + ',' + directed[1];
                    // std::map::insert does not overwrite an existing key.
                    if (!boundaryEdges.has(key)) {
                        boundaryEdges.set(key, { t: t, a: a, visited: false });
                        directedKeys.push(directed);
                    }
                }
            }
        }
        directedKeys.sort((e0, e1) => (e0[0] !== e1[0] ? e0[0] - e1[0] : e0[1] - e1[1]));

        // Extract the polygons. Each polygon is the boundary for a connected
        // component of the mesh.
        for (const key of directedKeys) {
            const initialEdge = boundaryEdges.get(key[0] + ',' + key[1]) as VETBoundaryEdge;
            if (!initialEdge.visited) {
                polygons.push(this.getBoundaryPolygon(initialEdge.t, initialEdge.a,
                    boundaryEdges));
            }
        }

        if (!duplicateEndpoints) {
            for (const polygon of polygons) {
                polygon.length = polygon.length - 1;
            }
        }

        return polygons;
    }

    // Count the number of triangles sharing each vertex. The total number of
    // indices for triangles adjacent to vertices is 3 * numTriangles. This is
    // easy to see from the code where an increment occurs 3 times per
    // triangle.
    protected getNumTrianglesAtVertex(counts: number[]): void {
        for (const tri of this.mTriangles) {
            for (let i = 0; i < 3; ++i) {
                ++counts[tri[i]];
            }
        }

        let minCount = counts[0];
        let maxCount = counts[0];
        for (const count of counts) {
            if (count < minCount) {
                minCount = count;
            }
            if (count > maxCount) {
                maxCount = count;
            }
        }
        this.mMinTrianglesAtVertex = minCount;
        this.mMaxTrianglesAtVertex = maxCount;
    }

    // Assign the storage subblocks to the vertices. The number of adjacent
    // vertices is incremented later during a triangle traversal and is used
    // as an index into the vertex-adjacency subblock during the traversal.
    protected initializeVertexStorage(numTrianglesAtVertex: readonly number[]): void {
        let storageOffset = 0;
        for (let v = 0; v < this.mVertices.length; ++v) {
            storageOffset = this.mVertices[v].initialize(numTrianglesAtVertex[v],
                storageOffset);
        }
    }

    // Populate each vertex with its adjacent L-triangle, adjacent vertices
    // and outgoing edges.
    protected populateVertices(): void {
        for (let t = 0; t < this.mTriangles.length; ++t) {
            const tri = this.mTriangles[t];
            const v0 = tri[0];
            const v1 = tri[1];
            const v2 = tri[2];

            // Update the adjacency information at v0.
            const vertex0 = this.mVertices[v0];
            vertex0.insertVAdjacent(v1);
            vertex0.insertVAdjacent(v2);
            vertex0.insertEAdjacent(v1, t);

            // Update the adjacency information at v1.
            const vertex1 = this.mVertices[v1];
            vertex1.insertVAdjacent(v2);
            vertex1.insertVAdjacent(v0);
            vertex1.insertEAdjacent(v2, t);

            // Update the adjacency information at v2.
            const vertex2 = this.mVertices[v2];
            vertex2.insertVAdjacent(v0);
            vertex2.insertVAdjacent(v1);
            vertex2.insertEAdjacent(v0, t);
        }
    }

    // Update triangle adjacency information for edges that are shared by two
    // triangles.
    protected updateAdjacencyForSharedEdges(): void {
        for (let t = 0; t < this.mTriangles.length; ++t) {
            this.updateAdjacencyForTriangle(t);
        }
    }

    protected updateAdjacencyForTriangle(t: number): void {
        const tri = this.mTriangles[t];
        for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
            // Get an outgoing edge <v0,v1>.
            const v0 = tri[i0];
            const v1 = tri[i1];

            // The outgoing edge from v0 is <v0,v1> and has adjacency triple
            // <v1,LT0,invalid>. If <v1,v0> is an outgoing edge from v1 with
            // adjacency triple <v0,LT1,invalid>, update the v0 adjacent
            // triple to <v1,LT0,LT1>; that is, RT0 = LT1. Although it is
            // possible at this time to update the v1 adjacent triple to
            // <v0,LT1,LT0>, where RT1 = LT0, the triple will be processed
            // when the outgoing edge is visited at another time. By not
            // updating <v0,LT1,invalid> now, two writes are avoided to each
            // of RT0 and RT1.
            const edge0 = this.mVertices[v0].findEAdjacentOffset(v1);
            const edge1 = this.mVertices[v1].findEAdjacentOffset(v0);
            if (edge0 >= 0 && edge1 >= 0) {
                this.mStorage[edge0 + 2] = this.mStorage[edge1 + 1];  // RT0 = LT1
                this.mAdjacents[t][i0] = this.mStorage[edge1 + 1];
            }
        }
    }

    protected depthFirstSearch(tInitial: number, visited: number[],
        tStack: number[], component: number[]): number {
        // The initial 'top' value denotes an empty stack. Upstream relies on
        // size_t wraparound for the empty-stack sentinel; the port uses -1.
        let top = -1;
        let numInserted = 0;

        tStack[++top] = tInitial;
        while (top !== -1) {
            const t = tStack[top];
            visited[t] = 1;
            let i = 0;
            for (i = 0; i < 3; ++i) {
                const tAdjacent = this.mAdjacents[t][i];
                if (tAdjacent !== StaticVETManifoldMesh2.invalid && visited[tAdjacent] === 0) {
                    tStack[++top] = tAdjacent;
                    break;
                }
            }
            if (i === 3) {
                visited[t] = 2;
                component[numInserted] = t;
                ++numInserted;
                --top;
            }
        }

        return numInserted;
    }

    protected getBoundaryPolygon(t: number, a: number,
        boundaryEdges: Map<string, VETBoundaryEdge>): number[] {
        const polygon: number[] = [];
        let tri = this.mTriangles[t];
        let i0 = a;
        let i1 = (i0 + 1) % 3;
        const vEdge: [number, number] = [tri[i0], tri[i1]];
        polygon.push(vEdge[0]);

        // The port of the upstream std::map::operator[], which default
        // constructs an element when the key is not in the map.
        const edgeAt = (key: [number, number]): VETBoundaryEdge => {
            const mapKey = key[0] + ',' + key[1];
            let edge = boundaryEdges.get(mapKey);
            if (edge === undefined) {
                edge = {
                    t: StaticVETManifoldMesh2.invalid,
                    a: StaticVETManifoldMesh2.invalid,
                    visited: false
                };
                boundaryEdges.set(mapKey, edge);
            }
            return edge;
        };

        while (!edgeAt(vEdge).visited) {
            polygon.push(vEdge[1]);
            edgeAt(vEdge).visited = true;

            // Traverse the triangle strip with vertex at vEdge[1] until the
            // last triangle is encountered. The final edge of the last
            // triangle is the next boundary edge and starts at vEdge[1].
            a = this.mAdjacents[t][i1];
            while (a !== StaticVETManifoldMesh2.invalid) {
                // Get the next triangle in the strip.
                t = a;
                tri = this.mTriangles[t];
                for (i1 = 0; i1 < 3; ++i1) {
                    if (vEdge[1] === tri[i1]) {
                        // Get the next interior edge in the triangle strip,
                        // namely, <tri[i0], tri[i1]>.
                        i0 = (i1 + 1) % 3;
                        a = this.mAdjacents[t][i1];
                        break;
                    }
                }
                logAssert(i1 < 3, 'Unexpected condition.');
            }

            const i2 = (i1 + 1) % 3;
            vEdge[0] = vEdge[1];
            vEdge[1] = tri[i2];
            i0 = i1;
            i1 = i2;
        }

        return polygon;
    }
}
