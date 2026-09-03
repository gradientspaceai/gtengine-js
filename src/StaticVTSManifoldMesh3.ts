// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) StaticVTSManifoldMesh3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// StaticVTSManifoldMesh3 represents a vertex-triangle-simplex manifold mesh
// for which tetrahedra (the simplices) are provided as a single batch and no
// mesh modification operations are going to be performed on the mesh. The
// class minimizes the memory management costs of a dynamic mesh. It is a
// requirement that the input tetrahedra form a manifold mesh with
// consistently ordered tetrahedra. In most applications, this requirement is
// already satisfied. See the comments for the static class member 'face'
// regarding ordering of tetrahedra.
//
// Port notes: the upstream 'invalid' sentinel std::numeric_limits<size_t>::max()
// becomes StaticVTSManifoldMesh3.invalid = Number.MAX_SAFE_INTEGER, which is
// representable exactly by a float64 and compares greater than any valid
// index. The upstream contiguous std::vector<size_t> storage block, whose
// subblocks are handed to the vertices as raw pointers, becomes a single
// Float64Array with the vertices holding offsets into it; the
// std::array<size_t,4>* face adjacency records become quadruples of
// consecutive elements. The upstream index-range tests 'v < mVertices.size()'
// rely on size_t wraparound to reject negative inputs, so the port tests
// '0 <= v && v < numVertices' to obtain the same behavior. The multithreaded
// construction path is omitted (JavaScript has no shared-memory threads
// here); the numThreads argument is accepted for call-site parity and the
// construction always runs the single-threaded path, which computes exactly
// the same result. The out-parameters of GetAdjacentTetrahedra and SortFace
// become return values. getAdjacentTetrahedra also fixes an upstream bug; see
// the comment in that method.

import { logAssert } from './Logger.js';

// The port of the nested class StaticVTSManifoldMesh3::Vertex. The name is
// prefixed because the library exports one flat namespace.
export class StaticVTSManifoldMesh3Vertex {
    // Only StaticVTSManifoldMesh3 writes these members.
    private mStorage: Float64Array;
    private mNumSAdjacents: number;
    private mNumVAdjacents: number;  // <= 3 * mNumSAdjacents
    private mVAdjacents: number;     // offset of [3 * mNumSAdjacents] block
    private mNumFAdjacents: number;  // <= 3 * mNumSAdjacents after construction
    private mFAdjacents: number;     // offset of [12 * mNumSAdjacents] block

    constructor(storage: Float64Array) {
        this.mStorage = storage;
        this.mNumSAdjacents = 0;
        this.mNumVAdjacents = 0;
        this.mVAdjacents = 0;
        this.mNumFAdjacents = 0;
        this.mFAdjacents = 0;
    }

    // The members are read-only.
    getNumSAdjacents(): number {
        return this.mNumSAdjacents;
    }

    // The number of adjacent vertices is bounded by three times the number of
    // tetrahedra sharing the vertex.
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

    // The number of adjacent (outgoing) faces is bounded by three times the
    // number of tetrahedra sharing the vertex.
    getNumFAdjacents(): number {
        return this.mNumFAdjacents;
    }

    // The adjacency quadruple <AV0,AV1,LS,RS> of the j-th outgoing face: AV0
    // and AV1 are the other two vertices of the face, LS is the tetrahedron
    // to the left of the face and RS is the tetrahedron to the right of the
    // face (RS is 'invalid' when the face is a boundary face).
    getFAdjacent(j: number): [number, number, number, number] {
        const base = this.mFAdjacents + 4 * j;
        return [
            this.mStorage[base],
            this.mStorage[base + 1],
            this.mStorage[base + 2],
            this.mStorage[base + 3]
        ];
    }

    getFAdjacents(): [number, number, number, number][] {
        const adjacents: [number, number, number, number][] = [];
        for (let j = 0; j < this.mNumFAdjacents; ++j) {
            adjacents.push(this.getFAdjacent(j));
        }
        return adjacents;
    }

    // The remaining members are the port of the private upstream members;
    // they are internal to the mesh construction. TypeScript has no 'friend'
    // declaration, so they are marked internal by convention.

    /** @internal */
    initialize(numSAdjacent: number, storageOffset: number): number {
        this.mNumSAdjacents = numSAdjacent;
        this.mNumVAdjacents = 0;
        this.mVAdjacents = storageOffset;
        storageOffset += 3 * this.mNumSAdjacents;
        this.mNumFAdjacents = 0;
        this.mFAdjacents = storageOffset;
        // 3 faces per vertex, 4 indices per element.
        storageOffset += 12 * this.mNumSAdjacents;
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
    insertFAdjacent(v1: number, v2: number, s: number): void {
        const base = this.mFAdjacents + 4 * this.mNumFAdjacents;
        this.mStorage[base] = v1;
        this.mStorage[base + 1] = v2;
        this.mStorage[base + 2] = s;
        this.mStorage[base + 3] = StaticVTSManifoldMesh3.invalid;
        ++this.mNumFAdjacents;
    }

    // Return the storage offset of the adjacency quadruple of the outgoing
    // face <this,w1,w2>, or -1 when there is no such outgoing face. This is
    // the port of the pointer returned by
    // StaticVTSManifoldMesh3::GetOutgoingFace.
    /** @internal */
    findFAdjacentOffset(w1: number, w2: number): number {
        for (let j = 0; j < this.mNumFAdjacents; ++j) {
            const base = this.mFAdjacents + 4 * j;
            if (this.mStorage[base] === w1 && this.mStorage[base + 1] === w2) {
                return base;
            }
        }
        return -1;
    }
}

export class StaticVTSManifoldMesh3 {
    // Use the maximum safe integer to denote an invalid index. Upstream uses
    // the maximum size_t, effectively representing -1.
    static readonly invalid: number = Number.MAX_SAFE_INTEGER;

    // The tetrahedron is represented as an array of four vertices, V[i] for
    // 0 <= i <= 3. The vertices are ordered so that the triangular faces are
    // counterclockwise-ordered triangles when viewed by an observer outside
    // the tetrahedron: face 0 = <V[0],V[2],V[1]>, face 1 = <V[0],V[1],V[3]>,
    // face 2 = <V[0],V[3],V[2]> and face 3 = <V[1],V[2],V[3]>. The canonical
    // tetrahedron has V[0] = (0,0,0), V[1] = (1,0,0), V[2] = (0,1,0) and
    // V[3] = (0,0,1).
    static readonly face: readonly (readonly [number, number, number])[] = [
        [0, 2, 1],
        [0, 1, 3],
        [0, 3, 2],
        [1, 2, 3]
    ];

    protected mVertices: StaticVTSManifoldMesh3Vertex[];
    protected mStorage: Float64Array;
    protected mTetrahedra: [number, number, number, number][];
    protected mAdjacents: [number, number, number, number][];
    protected mMinTetrahedraAtVertex: number;
    protected mMaxTetrahedraAtVertex: number;

    // Preconditions.
    //   1. The tetrahedra input must have size 1 or larger.
    //   2. The number of vertices must be 4 or larger.
    //   3. The tetrahedra must form a manifold mesh.
    //   4. Each tetrahedron must be nondegenerate; no repeated vertices.
    //   5. The tetrahedra must all be ordered counterclockwise or all ordered
    //      clockwise; no mixed chirality.
    // The numThreads argument is accepted for parity with upstream but the
    // construction always occurs single threaded.
    constructor(
        numVertices: number,
        tetrahedra: readonly (readonly number[])[],
        numThreads: number = 0) {
        void numThreads;
        logAssert(numVertices >= 4 && tetrahedra.length > 0, 'invalid input');

        this.mStorage = new Float64Array(60 * tetrahedra.length)
            .fill(StaticVTSManifoldMesh3.invalid);
        this.mVertices = new Array<StaticVTSManifoldMesh3Vertex>(numVertices);
        for (let v = 0; v < numVertices; ++v) {
            this.mVertices[v] = new StaticVTSManifoldMesh3Vertex(this.mStorage);
        }
        this.mTetrahedra = tetrahedra.map(
            tetra => [tetra[0], tetra[1], tetra[2], tetra[3]] as
                [number, number, number, number]);
        this.mAdjacents = tetrahedra.map(() => [
            StaticVTSManifoldMesh3.invalid,
            StaticVTSManifoldMesh3.invalid,
            StaticVTSManifoldMesh3.invalid,
            StaticVTSManifoldMesh3.invalid
        ] as [number, number, number, number]);
        this.mMinTetrahedraAtVertex = 0;
        this.mMaxTetrahedraAtVertex = 0;

        const numTetrahedraAtVertex: number[] = new Array<number>(numVertices).fill(0);
        this.getNumTetrahedraAtVertex(numTetrahedraAtVertex);
        this.initializeVertexStorage(numTetrahedraAtVertex);
        this.populateVertices();
        this.updateAdjacencyForSharedFaces();
    }

    // Member access.
    getVertices(): StaticVTSManifoldMesh3Vertex[] {
        return this.mVertices;
    }

    // Each 4-tuple contains indices into the vertices.
    getTetrahedra(): [number, number, number, number][] {
        return this.mTetrahedra;
    }

    // Each 4-tuple contains indices into the tetrahedra.
    getAdjacents(): [number, number, number, number][] {
        return this.mAdjacents;
    }

    getMinNumTetrahedraAtVertex(): number {
        return this.mMinTetrahedraAtVertex;
    }

    getMaxNumTetrahedraAtVertex(): number {
        return this.mMaxTetrahedraAtVertex;
    }

    // Determine whether or not the unordered face <v0,v1,v2> exists.
    faceExists(v0: number, v1: number, v2: number): boolean {
        const numVertices = this.mVertices.length;
        if (0 <= v0 && v0 < numVertices && 0 <= v1 && v1 < numVertices &&
            0 <= v2 && v2 < numVertices && v0 !== v1 && v0 !== v2 && v1 !== v2) {
            // Sort the face to <u0,u1,u2> where u0 = min(u0,u1,u2) and the
            // face is CCW when viewed from outside the tetrahedron.
            const u = StaticVTSManifoldMesh3.sortFace(v0, v1, v2);

            if (this.mVertices[u[0]].findFAdjacentOffset(u[1], u[2]) >= 0) {
                return true;
            }

            // An outgoing face was not found. Try to find an incoming face.
            if (this.mVertices[u[0]].findFAdjacentOffset(u[2], u[1]) >= 0) {
                return true;
            }
        }
        return false;
    }

    // Get the adjacent tetrahedra for the unordered face <v0,v1,v2>. The
    // returned adjacent tetrahedron indices adj0 and adj1 are the following:
    //
    //   1. <v0,v1,v2> and <v0,v2,v1> are both outgoing faces, so the face is
    //      shared by two tetrahedra and both adj0 and adj1 are valid (not
    //      equal to 'invalid'). The index adj0 is the L-tetrahedron for
    //      <v0,v1,v2> and the index adj1 is the R-tetrahedron for <v0,v1,v2>.
    //      Equivalently, adj0 is the R-tetrahedron for <v0,v2,v1> and adj1 is
    //      the L-tetrahedron for <v0,v2,v1>.
    //
    //   2. <v0,v1,v2> is outgoing but <v0,v2,v1> is not outgoing. The index
    //      adj0 is the L-tetrahedron for <v0,v1,v2> and the index adj1 is
    //      'invalid' (no R-tetrahedron).
    //
    //   3. <v0,v2,v1> is outgoing but <v0,v1,v2> is not outgoing. The index
    //      adj0 is 'invalid' (no L-tetrahedron) for <v0,v2,v1> and the index
    //      adj1 is the R-tetrahedron for <v0,v2,v1>.
    //
    //   4. The unordered face <v0,v1,v2> does not exist.
    //
    // It is possible to distinguish between the 4 cases by examining the
    // returned indices:
    //   (1) returns (valid, valid) and Boolean 'true'
    //   (2) returns (valid, invalid) and Boolean 'true'
    //   (3) returns (invalid, valid) and Boolean 'true'
    //   (4) returns (invalid, invalid) and Boolean 'false'
    getAdjacentTetrahedra(v0: number, v1: number, v2: number):
        { exists: boolean, adj0: number, adj1: number } {
        const numVertices = this.mVertices.length;
        if (0 <= v0 && v0 < numVertices && 0 <= v1 && v1 < numVertices &&
            0 <= v2 && v2 < numVertices && v0 !== v1 && v0 !== v2 && v1 !== v2) {
            // Sort the face to <u0,u1,u2> where u0 = min(u0,u1,u2) and the
            // face is CCW when viewed from outside the tetrahedron.
            const u = StaticVTSManifoldMesh3.sortFace(v0, v1, v2);

            // UPSTREAM BUG (StaticVTSManifoldMesh3.h, GetAdjacentTetrahedra).
            // When the outgoing face <u0,u1,u2> is not found, upstream
            // searches for the opposite outgoing face <u0,u2,u1> and then
            // copies that quad's components [2] and [3] straight into
            // (adj0,adj1). Those are the L-tetrahedron and R-tetrahedron of
            // <u0,u2,u1>, which are the R-tetrahedron and the L-tetrahedron
            // of the queried orientation <u0,u1,u2>, so the pair must be
            // swapped. Without the swap the documented case 3, "returns
            // (invalid, valid)", can never be produced: a boundary face whose
            // only outgoing orientation is <u0,u2,u1> yields (valid, invalid),
            // which is the signature documented for case 2. The port performs
            // the two oriented lookups the documented contract requires: adj0
            // is the L-tetrahedron of <u0,u1,u2> and adj1 is the
            // L-tetrahedron of <u0,u2,u1> (its R-tetrahedron counterpart).
            const vertex = this.mVertices[u[0]];
            const offset0 = vertex.findFAdjacentOffset(u[1], u[2]);
            const offset1 = vertex.findFAdjacentOffset(u[2], u[1]);
            if (offset0 >= 0 || offset1 >= 0) {
                return {
                    exists: true,
                    adj0: offset0 >= 0
                        ? this.mStorage[offset0 + 2] : StaticVTSManifoldMesh3.invalid,
                    adj1: offset1 >= 0
                        ? this.mStorage[offset1 + 2] : StaticVTSManifoldMesh3.invalid
                };
            }
        }

        return {
            exists: false,
            adj0: StaticVTSManifoldMesh3.invalid,
            adj1: StaticVTSManifoldMesh3.invalid
        };
    }

    // Count the number of tetrahedra sharing each vertex. The total number of
    // indices for tetrahedra adjacent to vertices is 4 * numTetrahedra. This
    // is easy to see from the code where an increment occurs 4 times per
    // tetrahedron.
    protected getNumTetrahedraAtVertex(counts: number[]): void {
        for (const tetra of this.mTetrahedra) {
            for (let i = 0; i < 4; ++i) {
                ++counts[tetra[i]];
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
        this.mMinTetrahedraAtVertex = minCount;
        this.mMaxTetrahedraAtVertex = maxCount;
    }

    // Assign the storage subblocks to the vertices. The number of adjacent
    // vertices is incremented later during a tetrahedron traversal and is
    // used as an index into the vertex-adjacency subblock during the
    // traversal.
    protected initializeVertexStorage(numTetrahedraAtVertex: readonly number[]): void {
        let storageOffset = 0;
        for (let v = 0; v < this.mVertices.length; ++v) {
            storageOffset = this.mVertices[v].initialize(numTetrahedraAtVertex[v],
                storageOffset);
        }
    }

    // Populate each vertex with its adjacent vertices.
    protected updateVertexAdjacents(t: number): void {
        // Iterate over all vertex pairs (u0,u1). Update u0's vertex adjacents
        // with u1 and update u1's vertex adjacents with u0.
        const tetra = this.mTetrahedra[t];
        const v0 = tetra[0];
        const v1 = tetra[1];
        const v2 = tetra[2];
        const v3 = tetra[3];
        const vertex0 = this.mVertices[v0];
        const vertex1 = this.mVertices[v1];
        const vertex2 = this.mVertices[v2];
        const vertex3 = this.mVertices[v3];
        vertex0.insertVAdjacent(v1);
        vertex0.insertVAdjacent(v2);
        vertex0.insertVAdjacent(v3);
        vertex1.insertVAdjacent(v2);
        vertex1.insertVAdjacent(v0);
        vertex1.insertVAdjacent(v3);
        vertex2.insertVAdjacent(v0);
        vertex2.insertVAdjacent(v1);
        vertex2.insertVAdjacent(v3);
        vertex3.insertVAdjacent(v1);
        vertex3.insertVAdjacent(v0);
        vertex3.insertVAdjacent(v2);
    }

    // Populate each vertex with its adjacent L-tetrahedra and adjacent faces.
    protected updateFaceAdjacents(t: number): void {
        const tetra = this.mTetrahedra[t];
        for (let i = 0; i < 4; ++i) {
            // Get an outgoing face <v0,v1,v2>, which is CCW when viewed from
            // outside the tetrahedron.
            const faceIndices = StaticVTSManifoldMesh3.face[i];
            const v0 = tetra[faceIndices[0]];
            const v1 = tetra[faceIndices[1]];
            const v2 = tetra[faceIndices[2]];

            // Sort the face to <u0,u1,u2> where u0 = min(u0,u1,u2) and the
            // face is CCW when viewed from outside the tetrahedron.
            const u = StaticVTSManifoldMesh3.sortFace(v0, v1, v2);

            // Update the face adjacency information at u0.
            this.mVertices[u[0]].insertFAdjacent(u[1], u[2], t);
        }
    }

    protected populateVertices(): void {
        for (let t = 0; t < this.mTetrahedra.length; ++t) {
            this.updateVertexAdjacents(t);
            this.updateFaceAdjacents(t);
        }
    }

    // Update tetrahedra adjacency information for faces that are shared by
    // two tetrahedra.
    protected updateAdjacencyForSharedFaces(): void {
        for (let t = 0; t < this.mTetrahedra.length; ++t) {
            this.updateAdjacencyForTetrahedron(t);
        }
    }

    protected updateAdjacencyForTetrahedron(t: number): void {
        const tetra = this.mTetrahedra[t];
        for (let i = 0; i < 4; ++i) {
            // Get an outgoing face <v0,v1,v2>, which is CCW when viewed from
            // outside the tetrahedron.
            const faceIndices = StaticVTSManifoldMesh3.face[i];
            const v0 = tetra[faceIndices[0]];
            const v1 = tetra[faceIndices[1]];
            const v2 = tetra[faceIndices[2]];

            // Sort the face to <u0,u1,u2> where u0 = min(u0,u1,u2) and the
            // face is CCW when viewed from outside the tetrahedron.
            const u = StaticVTSManifoldMesh3.sortFace(v0, v1, v2);

            // The outgoing face from u0 is <u0,u1,u2> and has adjacency quad
            // <u1,u2,LT0,invalid>. If <u0,u2,u1> is an outgoing face from u0
            // with adjacency quad <u2,u1,LT1,invalid>, update the u0 adjacent
            // quad to <u1,u2,LT0,LT1>; that is, RT0 = LT1. Although it is
            // possible at this time to update the u0 adjacent quad to
            // <u2,u1,LT1,LT0>, where RT1 = LT0, this quad will be processed
            // when the outgoing face is visited at another time. By not
            // updating <u2,u1,LT1,invalid> now, two writes are avoided to
            // each of RT0 and RT1.
            const vertex = this.mVertices[u[0]];
            const face0 = vertex.findFAdjacentOffset(u[1], u[2]);
            const face1 = vertex.findFAdjacentOffset(u[2], u[1]);
            if (face0 >= 0 && face1 >= 0) {
                this.mStorage[face0 + 3] = this.mStorage[face1 + 2];  // RT0 = LT1
                this.mAdjacents[t][i] = this.mStorage[face1 + 2];
            }
        }
    }

    // Sort the face <v0,v1,v2> to <u0,u1,u2> where u0 = min(v0,v1,v2) and the
    // cyclic order of the vertices is preserved, so the face keeps its
    // orientation.
    static sortFace(v0: number, v1: number, v2: number): [number, number, number] {
        if (v0 < v1) {
            if (v0 < v2) {
                // v0 is minimum
                return [v0, v1, v2];
            } else {
                // v2 is minimum
                return [v2, v0, v1];
            }
        } else {
            if (v1 < v2) {
                // v1 is minimum
                return [v1, v2, v0];
            } else {
                // v2 is minimum
                return [v2, v0, v1];
            }
        }
    }
}
