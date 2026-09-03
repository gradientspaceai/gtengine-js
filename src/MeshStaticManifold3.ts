// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MeshStaticManifold3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// MeshStaticManifold3 represents a vertex-face-tetrahedron manifold mesh for
// which tetrahedra (the simplices) are provided as a single batch and no mesh
// modification operations are going to be performed on the mesh. It is a
// requirement that the input tetrahedra form a manifold mesh with
// consistently ordered tetrahedra. In most applications, this requirement is
// already satisfied.
//
// Port notes:
//   - Upstream stores the vertex adjacency 5-tuples in one contiguous block
//     of 20 * numTetrahedra size_t values and hands each vertex a pointer
//     into that block. The port gives each vertex its own array of 5-tuples;
//     the insertion order, and therefore every observable result, is the
//     same.
//   - Upstream optionally distributes the adjacency update over threads. The
//     port always runs the single-threaded loop, which computes the same
//     values. The numThreads constructor parameter is retained and ignored.
//   - 'invalid' is Number.MAX_SAFE_INTEGER rather than SIZE_MAX, so the
//     upstream 'index < numVertices' guards still reject it.
//   - getAdjacentTetrahedra returns { adj0, adj1, exists } instead of using
//     output reference parameters.
//   - getAdjacentTetrahedra fixes an upstream bug; see the comment in that
//     method.

import { logAssert } from './Logger.js';

// The vertices are stored as an array of MeshStaticManifold3Vertex objects,
// one per vertex index. If tetrahedron[t0] = <v0,v1,v2,v3>, then vertex[v0]
// contains a 5-tuple {v1,v2,v3,t0,a0}. The unordered face (v1,v2,v3) is
// opposite v0. If there is no adjacent tetrahedron sharing (v1,v2,v3), then
// a0 is invalid. If there is an adjacent tetrahedron, then a0 is the index
// for that tetrahedron. Let tetrahedron[a0] = <v1,v3,v2,v4>; then vertex[v4]
// contains a 5-tuple {v1,v3,v2,a0,t0}.
export class MeshStaticManifold3Vertex {
    // If tetrahedron t0 is <v0,v1,v2,v3> in counterclockwise order, then the
    // corresponding adjacents element is {v1,v2,v3,t0,a0}, where a0 is
    // invalid when <v1,v2,v3> is contained by a single tetrahedron or a0 is
    // the index for the adjacent tetrahedron when <v1,v2,v3> is contained by
    // two tetrahedra.
    //
    // Upstream declares MeshStaticManifold3 a friend so that only it may
    // write the members. TypeScript has no equivalent; treat 'adjacents' and
    // insert() as internal to MeshStaticManifold3.
    readonly adjacents: number[][] = [];

    getNumAdjacents(): number {
        return this.adjacents.length;
    }

    getAdjacents(): number[][] {
        return this.adjacents;
    }

    insert(v1: number, v2: number, v3: number, s: number, location: number): void {
        this.adjacents.push([v1, v2, v3, s, location]);
    }
}

export class MeshStaticManifold3 {
    // Use the maximum safe integer to denote an invalid index.
    static readonly invalid = Number.MAX_SAFE_INTEGER;

    // The tetrahedron is represented as an array of four vertices, V[i] for
    // 0 <= i <= 3. The vertices are ordered so that the triangle faces are
    // counterclockwise ordered when viewed by an observer outside the
    // tetrahedron: face[0] = <V[1],V[2],V[3]>, face[1] = <V[0],V[3],V[2]>,
    // face[2] = <V[0],V[1],V[3]> and face[3] = <V[0],V[2],V[1]>. Observe that
    // for face[i], the vertex opposite the face is V[i]. The canonical
    // tetrahedron has V[0] = (0,0,0), V[1] = (1,0,0), V[2] = (0,1,0) and
    // V[3] = (0,0,1).
    static readonly face: readonly (readonly number[])[] = [
        [1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]
    ];

    protected mVertices: MeshStaticManifold3Vertex[];
    protected mTetrahedra: number[][];
    protected mAdjacents: number[][];
    protected mMinTetrahedraAtVertex: number;
    protected mMaxTetrahedraAtVertex: number;

    // Preconditions.
    //   1. The tetrahedra input must have size 1 or larger.
    //   2. The number of vertices must be 4 or larger.
    //   3. The tetrahedra must form a manifold mesh.
    //   4. Each tetrahedron must be nondegenerate; no repeated vertices.
    //   5. The tetrahedra must all be ordered counterclockwise.
    // The numThreads parameter of the upstream constructor is accepted for
    // signature compatibility but ignored; see the port notes.
    constructor(numVertices: number, tetrahedra: number[][], numThreads: number = 0) {
        void numThreads;
        logAssert(numVertices >= 4 && tetrahedra.length > 0, 'Invalid input.');

        this.mVertices = new Array<MeshStaticManifold3Vertex>(numVertices);
        for (let v = 0; v < numVertices; ++v) {
            this.mVertices[v] = new MeshStaticManifold3Vertex();
        }
        this.mTetrahedra = tetrahedra.map(tet => [tet[0], tet[1], tet[2], tet[3]]);
        this.mAdjacents = tetrahedra.map(() => [
            MeshStaticManifold3.invalid,
            MeshStaticManifold3.invalid,
            MeshStaticManifold3.invalid,
            MeshStaticManifold3.invalid
        ]);
        this.mMinTetrahedraAtVertex = 0;
        this.mMaxTetrahedraAtVertex = 0;

        this.getNumTetrahedraAtVertex(numVertices);
        this.populateVertices();
        this.updateAdjacencyForSharedFaces();
    }

    // Member access.
    getVertices(): MeshStaticManifold3Vertex[] {
        return this.mVertices;
    }

    // Each 4-tuple contains indices into the vertices.
    getTetrahedra(): number[][] {
        return this.mTetrahedra;
    }

    // Each 4-tuple contains indices into the tetrahedra.
    getAdjacents(): number[][] {
        return this.mAdjacents;
    }

    getMinNumTetrahedraAtVertex(): number {
        return this.mMinTetrahedraAtVertex;
    }

    getMaxNumTetrahedraAtVertex(): number {
        return this.mMaxTetrahedraAtVertex;
    }

    // Determine whether or not the unordered face (v0,v1,v2) exists.
    faceExists(v0: number, v1: number, v2: number): boolean {
        const numVertices = this.mVertices.length;
        if (v0 < numVertices && v1 < numVertices && v2 < numVertices &&
            v0 !== v1 && v0 !== v2 && v1 !== v2) {
            return this.getOrderedFace(v0, v1, v2) !== null
                || this.getOrderedFace(v0, v2, v1) !== null;
        }
        return false;
    }

    // Get the adjacent tetrahedra for the unordered face (v0,v1,v2). The
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
    //   4. Neither <v0,v1,v2> nor <v0,v2,v1> exist; that is, the face does
    //      not occur for any tetrahedron.
    //
    // It is possible to distinguish among the 4 cases by examining the
    // returned indices:
    //   (1) returns (valid, valid) and exists 'true'
    //   (2) returns (valid, invalid) and exists 'true'
    //   (3) returns (invalid, valid) and exists 'true'
    //   (4) returns (invalid, invalid) and exists 'false'
    getAdjacentTetrahedra(v0: number, v1: number, v2: number): {
        adj0: number;
        adj1: number;
        exists: boolean;
    } {
        const numVertices = this.mVertices.length;
        if (v0 < numVertices && v1 < numVertices && v2 < numVertices &&
            v0 !== v1 && v0 !== v2 && v1 !== v2) {
            // UPSTREAM BUG (MeshStaticManifold3.h, GetAdjacentTetrahedra),
            // the 3D form of the MeshStaticManifold2 bug documented in
            // MeshStaticManifold2.ts. Upstream returns components [2] and [3]
            // of a single 5-tuple. The 5-tuple stored at vertex[v0] for
            // tetrahedron <v0,v1,v2,v3> is {v1, v2, v3, tet, adj}, so [2] is
            // the vertex v3 and [3] is the containing tetrahedron: adj0 is
            // handed back a VERTEX index. (The 2D code uses [2] and [3] of a
            // 4-tuple, which suggests this is an unreindexed copy of it.)
            // Even with the indices shifted to [3] and [4], the value at [4]
            // is the neighbor across the face OPPOSITE v0, namely (v1,v2,v3),
            // not across the queried face (v0,v1,v2). And, as in 2D, the
            // documented case 3 return of (invalid, valid) is unreachable in
            // the upstream code.
            //
            // The port computes both ordered-face lookups, as the documented
            // contract requires: adj0 is the tetrahedron containing the
            // ordered face <v0,v1,v2> (the L-tetrahedron) and adj1 is the
            // tetrahedron containing <v0,v2,v1> (the R-tetrahedron).
            const adjacents0 = this.getOrderedFace(v0, v1, v2);
            const adjacents1 = this.getOrderedFace(v0, v2, v1);
            if (adjacents0 !== null || adjacents1 !== null) {
                return {
                    adj0: adjacents0 !== null ? adjacents0[3] : MeshStaticManifold3.invalid,
                    adj1: adjacents1 !== null ? adjacents1[3] : MeshStaticManifold3.invalid,
                    exists: true
                };
            }
        }

        return {
            adj0: MeshStaticManifold3.invalid,
            adj1: MeshStaticManifold3.invalid,
            exists: false
        };
    }

    // Count the number of tetrahedra sharing each vertex. The total number of
    // indices for tetrahedra adjacent to vertices is 4 * numTetrahedra. The
    // minimum and maximum tetrahedra counts are for statistical information.
    protected getNumTetrahedraAtVertex(numVertices: number): void {
        const counts = new Array<number>(numVertices).fill(0);
        for (const tet of this.mTetrahedra) {
            for (let i = 0; i < 4; ++i) {
                ++counts[tet[i]];
            }
        }

        this.mMinTetrahedraAtVertex = counts[0];
        this.mMaxTetrahedraAtVertex = counts[0];
        for (const count of counts) {
            if (count < this.mMinTetrahedraAtVertex) {
                this.mMinTetrahedraAtVertex = count;
            }
            if (count > this.mMaxTetrahedraAtVertex) {
                this.mMaxTetrahedraAtVertex = count;
            }
        }
    }

    // Populate the adjacency information for the vertices.
    protected populateVertices(): void {
        for (let t = 0; t < this.mTetrahedra.length; ++t) {
            const tet = this.mTetrahedra[t];
            const v0 = tet[0];
            const v1 = tet[1];
            const v2 = tet[2];
            const v3 = tet[3];

            // The last arguments (i = 0, 1, 2 or 3) are used to set the
            // correct mAdjacents[][i] indices. These arguments are replaced
            // later by the actual indices for adjacent tetrahedra sharing the
            // face.
            this.mVertices[v0].insert(v1, v2, v3, t, 0);
            this.mVertices[v1].insert(v0, v3, v2, t, 1);
            this.mVertices[v2].insert(v0, v1, v3, t, 2);
            this.mVertices[v3].insert(v0, v2, v1, t, 3);
        }
    }

    // Update tetrahedra adjacency information for faces that are shared by
    // two tetrahedra.
    protected updateAdjacencyForSharedFaces(): void {
        for (let v = 0; v < this.mVertices.length; ++v) {
            this.updateAdjacencyForFace(v);
        }
    }

    protected updateAdjacencyForFace(v0: number): void {
        const vertex0 = this.mVertices[v0];
        for (const adjacents0 of vertex0.adjacents) {
            const v1 = adjacents0[0];
            const v2 = adjacents0[1];
            const v3 = adjacents0[2];

            // The face opposite vertex v0 is (v1,v2,v3). We know that
            // vertex[v0] contains a 5-tuple {v1,v2,v3,tet0,loc0}. Determine
            // whether vertex[v1] contains a 5-tuple that has the face
            // <v1,v3,v2>.
            const adjacents1 = this.getOrderedFace(v1, v3, v2);
            if (adjacents1 !== null) {
                // The face <v1,v2,v3> has a tetrahedron adjacent to
                // tetrahedron tet0. Update the vertex adjacency information
                // for tetrahedron tet0 at that face. Tetrahedron a1 adjacency
                // is not updated. It will be updated when <v1,v3,v2> is
                // visited at another time. This avoids two writes of the
                // adjacent tetrahedron indices.
                const tet0 = adjacents0[3];
                const loc0 = adjacents0[4];
                const adj1 = adjacents1[3];
                adjacents0[4] = adj1;
                this.mAdjacents[tet0][loc0] = adj1;
            } else {
                // Replace the mAdjacents[] location value (0, 1, 2 or 3) by
                // an invalid index because face <v1,v3,v2> does not exist, in
                // which case there is no adjacent tetrahedron to face
                // <v1,v2,v3>.
                adjacents0[4] = MeshStaticManifold3.invalid;
            }
        }
    }

    protected getOrderedFace(v0: number, v1: number, v2: number): number[] | null {
        const vertex0 = this.mVertices[v0];
        for (const adjacents0 of vertex0.adjacents) {
            const tetra = [v0, adjacents0[0], adjacents0[1], adjacents0[2]];

            for (let f = 0; f < 4; ++f) {
                const compareFace = MeshStaticManifold3.face[f];
                if (v0 === tetra[compareFace[0]] &&
                    v1 === tetra[compareFace[1]] &&
                    v2 === tetra[compareFace[2]]) {
                    return adjacents0;
                }
            }
        }
        return null;
    }
}
