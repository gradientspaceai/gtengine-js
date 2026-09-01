// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) TSManifoldMesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The TSManifoldMesh class represents a triangle-tetrahedron mesh. The 'T'
// stands for triangle (face) and the 'S' stands for simplex (tetrahedron). It
// is general purpose, allowing insertion and removal of tetrahedra at any
// time. However, the performance is limited because of the use of the
// container classes. If your application requires a triangle-tetrahedron
// manifold mesh for which no tetrahedra will be removed, a better choice is
// StaticVTSManifoldMesh.
//
// Port notes: this is the tetrahedron analogue of ETManifoldMesh and follows
// the same port conventions. Upstream stores the faces in an
// std::unordered_map keyed by TriangleKey<false> and the tetrahedra in an
// std::unordered_map keyed by TetrahedronKey<true>; the port uses JavaScript
// Maps keyed by FeatureKey.mapKey() and, because the upstream containers are
// unordered (so C++ iteration order is implementation-defined), iterates in
// increasing feature-key order wherever the order is observable. The
// std::unique_ptr map values and the raw back-reference pointers become plain
// object references and null. The TCreator/SCreator function pointers become
// optional constructor callbacks, and the copy constructor and operator=
// become clone() and assign().
//
// Upstream quirk, preserved: when insert() rejects a tetrahedron because the
// mesh would become nonmanifold and throwOnNonmanifoldInsertion(false) has
// been set, upstream returns nullptr but leaves behind the state built by the
// loop iterations that already ran: faces created by this call stay in the
// face map with S[0] referencing a tetrahedron that was never added to the
// tetrahedron map, and an already-present adjacent tetrahedron may have had
// its S[j] pointed at that same phantom tetrahedron. This is the same shape
// as upstream issue #73 for VEManifoldMesh::Insert, and the port reproduces
// it rather than rolling the insertion back. Prefer the default throwing
// behavior.

import { logAssert, logError } from './Logger';
import { TriangleKey } from './TriangleKey';
import { TetrahedronKey } from './TetrahedronKey';
import { FeatureKey } from './FeatureKey';

// The port of TSManifoldMesh::Triangle, a face of the mesh.
export class TSManifoldMeshTriangle {
    // Vertices of the face. These are stored in the order in which the face
    // was first encountered by insert(); the unordered TriangleKey of
    // (V[0],V[1],V[2]) is the map key.
    V: [number, number, number];

    // Tetrahedra sharing the face.
    S: [TSManifoldMeshTetrahedron | null, TSManifoldMeshTetrahedron | null];

    constructor(v0: number, v1: number, v2: number) {
        this.V = [v0, v1, v2];
        this.S = [null, null];
    }
}

// The port of TSManifoldMesh::Tetrahedron.
export class TSManifoldMeshTetrahedron {
    // Vertices, listed in an order so that each face has vertices in
    // counterclockwise order when viewed from outside the tetrahedron.
    V: [number, number, number, number];

    // Adjacent faces. T[i] points to the triangle face opposite V[i].
    //   T[0] points to face (V[1],V[2],V[3])
    //   T[1] points to face (V[0],V[3],V[2])
    //   T[2] points to face (V[0],V[1],V[3])
    //   T[3] points to face (V[0],V[2],V[1])
    T: [TSManifoldMeshTriangle | null, TSManifoldMeshTriangle | null,
        TSManifoldMeshTriangle | null, TSManifoldMeshTriangle | null];

    // Adjacent tetrahedra. S[i] points to the adjacent tetrahedron sharing
    // face T[i].
    S: [TSManifoldMeshTetrahedron | null, TSManifoldMeshTetrahedron | null,
        TSManifoldMeshTetrahedron | null, TSManifoldMeshTetrahedron | null];

    constructor(v0: number, v1: number, v2: number, v3: number) {
        this.V = [v0, v1, v2, v3];
        this.T = [null, null, null, null];
        this.S = [null, null, null, null];
    }
}

export type TSManifoldMeshTCreator =
    (v0: number, v1: number, v2: number) => TSManifoldMeshTriangle;
export type TSManifoldMeshSCreator =
    (v0: number, v1: number, v2: number, v3: number) => TSManifoldMeshTetrahedron;

// The unordered key of a face, TriangleKey<false>(V[0],V[1],V[2]).
function faceKeyOf(face: TSManifoldMeshTriangle): TriangleKey {
    return new TriangleKey(false, face.V[0], face.V[1], face.V[2]);
}

// The ordered key of a tetrahedron, TetrahedronKey<true>(V[0],...,V[3]).
function tetraKeyOf(tetra: TSManifoldMeshTetrahedron): TetrahedronKey {
    return new TetrahedronKey(true, tetra.V[0], tetra.V[1], tetra.V[2], tetra.V[3]);
}

export class TSManifoldMesh {
    protected mTCreator: TSManifoldMeshTCreator;
    protected mTMap: Map<string, TSManifoldMeshTriangle>;
    protected mSCreator: TSManifoldMeshSCreator;
    protected mSMap: Map<string, TSManifoldMeshTetrahedron>;
    protected mThrowOnNonmanifoldInsertion: boolean;  // default: true

    constructor(tCreator?: TSManifoldMeshTCreator, sCreator?: TSManifoldMeshSCreator) {
        this.mTCreator = tCreator ?? TSManifoldMesh.createTriangle;
        this.mTMap = new Map<string, TSManifoldMeshTriangle>();
        this.mSCreator = sCreator ?? TSManifoldMesh.createTetrahedron;
        this.mSMap = new Map<string, TSManifoldMeshTetrahedron>();
        this.mThrowOnNonmanifoldInsertion = true;
    }

    // Support for a deep copy of the mesh; the port of operator=. The face
    // and tetrahedron objects are not shared between the meshes. Note that
    // the tetrahedra are reinserted using the vertices of the tetrahedron
    // keys rather than the vertices of the tetrahedron objects, so a copied
    // tetrahedron may have its vertices permuted relative to the original
    // (the key permutation is even, so the orientation is preserved). This is
    // the upstream behavior.
    assign(mesh: TSManifoldMesh): this {
        this.clear();

        this.mTCreator = mesh.mTCreator;
        this.mSCreator = mesh.mSCreator;
        this.mThrowOnNonmanifoldInsertion = mesh.mThrowOnNonmanifoldInsertion;
        for (const skey of mesh.getTetrahedronKeys()) {
            this.insert(skey.V[0], skey.V[1], skey.V[2], skey.V[3]);
        }

        return this;
    }

    // The port of the upstream copy constructor.
    clone(): TSManifoldMesh {
        return new TSManifoldMesh().assign(this);
    }

    // Member access. Upstream returns the containers themselves; the port
    // returns arrays of the values in increasing feature-key order.
    getTriangles(): TSManifoldMeshTriangle[] {
        const faces = Array.from(this.mTMap.values());
        faces.sort((f0, f1) => FeatureKey.compare(faceKeyOf(f0), faceKeyOf(f1)));
        return faces;
    }

    getTetrahedra(): TSManifoldMeshTetrahedron[] {
        const tetrahedra = Array.from(this.mSMap.values());
        tetrahedra.sort((s0, s1) => FeatureKey.compare(tetraKeyOf(s0), tetraKeyOf(s1)));
        return tetrahedra;
    }

    getTriangleKeys(): TriangleKey[] {
        return this.getTriangles().map(face => faceKeyOf(face));
    }

    getTetrahedronKeys(): TetrahedronKey[] {
        return this.getTetrahedra().map(tetra => tetraKeyOf(tetra));
    }

    // The port of the upstream map find operations.
    getTriangle(v0: number, v1: number, v2: number): TSManifoldMeshTriangle | null {
        return this.mTMap.get(new TriangleKey(false, v0, v1, v2).mapKey()) ?? null;
    }

    getTetrahedron(v0: number, v1: number, v2: number, v3: number):
        TSManifoldMeshTetrahedron | null {
        return this.mSMap.get(new TetrahedronKey(true, v0, v1, v2, v3).mapKey()) ?? null;
    }

    getNumTriangles(): number {
        return this.mTMap.size;
    }

    getNumTetrahedra(): number {
        return this.mSMap.size;
    }

    // If the insertion of a tetrahedron fails because the mesh would become
    // nonmanifold, the default behavior is to throw an exception. You can
    // disable this behavior and continue gracefully without an exception. The
    // return value is the previous value of the internal state.
    throwOnNonmanifoldInsertion(doException: boolean): boolean {
        const previous = this.mThrowOnNonmanifoldInsertion;
        this.mThrowOnNonmanifoldInsertion = doException;
        return previous;
    }

    // If <v0,v1,v2,v3> is not in the mesh, a Tetrahedron object is created
    // and returned; otherwise, <v0,v1,v2,v3> is in the mesh and null is
    // returned. If the insertion leads to a nonmanifold mesh, the call fails
    // with null returned (or an exception, see throwOnNonmanifoldInsertion).
    insert(v0: number, v1: number, v2: number, v3: number): TSManifoldMeshTetrahedron | null {
        const skey = new TetrahedronKey(true, v0, v1, v2, v3).mapKey();
        if (this.mSMap.has(skey)) {
            // The tetrahedron already exists. Return null as a signal to the
            // caller that the insertion failed.
            return null;
        }

        // Add the new tetrahedron.
        const tetra = this.mSCreator(v0, v1, v2, v3);

        // Add the faces to the mesh if they do not already exist.
        for (let i = 0; i < 4; ++i) {
            const opposite = TetrahedronKey.getOppositeFace()[i];
            const u0 = tetra.V[opposite[0]];
            const u1 = tetra.V[opposite[1]];
            const u2 = tetra.V[opposite[2]];
            const tkey = new TriangleKey(false, u0, u1, u2).mapKey();
            let face = this.mTMap.get(tkey);
            if (face === undefined) {
                // This is the first time the face is encountered.
                face = this.mTCreator(u0, u1, u2);
                this.mTMap.set(tkey, face);

                // Update the face and tetrahedron.
                face.S[0] = tetra;
                tetra.T[i] = face;
            } else {
                // This is the second time the face is encountered.
                logAssert(face !== null, 'Unexpected condition.');

                // Update the face.
                if (face.S[1]) {
                    if (this.mThrowOnNonmanifoldInsertion) {
                        logError('Attempt to create nonmanifold mesh.');
                    } else {
                        return null;
                    }
                }
                face.S[1] = tetra;

                // Update the adjacent tetrahedra.
                const adjacent = face.S[0];
                logAssert(adjacent !== null, 'Unexpected condition.');
                for (let j = 0; j < 4; ++j) {
                    if (adjacent.T[j] === face) {
                        adjacent.S[j] = tetra;
                        break;
                    }
                }

                // Update the tetrahedron.
                tetra.T[i] = face;
                tetra.S[i] = adjacent;
            }
        }

        this.mSMap.set(skey, tetra);
        return tetra;
    }

    // If <v0,v1,v2,v3> is in the mesh, it is removed and 'true' is returned;
    // otherwise, <v0,v1,v2,v3> is not in the mesh and 'false' is returned.
    remove(v0: number, v1: number, v2: number, v3: number): boolean {
        const skey = new TetrahedronKey(true, v0, v1, v2, v3).mapKey();
        const tetra = this.mSMap.get(skey);
        if (tetra === undefined) {
            // The tetrahedron does not exist.
            return false;
        }

        // Remove the faces and update adjacent tetrahedra if necessary.
        for (let i = 0; i < 4; ++i) {
            // Inform the faces the tetrahedron is being deleted.
            const face: TSManifoldMeshTriangle | null = tetra.T[i];
            logAssert(face !== null, 'Unexpected condition.');

            if (face.S[0] === tetra) {
                // One-tetrahedron faces always have the reference at index
                // zero.
                face.S[0] = face.S[1];
                face.S[1] = null;
            } else if (face.S[1] === tetra) {
                face.S[1] = null;
            } else {
                logError('Unexpected condition.');
            }

            // Remove the face if you have the last reference to it.
            if (!face.S[0] && !face.S[1]) {
                this.mTMap.delete(faceKeyOf(face).mapKey());
            }

            // Inform adjacent tetrahedra the tetrahedron is being deleted.
            const adjacent: TSManifoldMeshTetrahedron | null = tetra.S[i];
            if (adjacent) {
                for (let j = 0; j < 4; ++j) {
                    if (adjacent.S[j] === tetra) {
                        adjacent.S[j] = null;
                        break;
                    }
                }
            }
        }

        this.mSMap.delete(skey);
        return true;
    }

    // Destroy the triangles and tetrahedra to obtain an empty mesh.
    clear(): void {
        this.mTMap.clear();
        this.mSMap.clear();
    }

    // A manifold mesh is closed if each face is shared twice.
    isClosed(): boolean {
        for (const face of this.mTMap.values()) {
            if (!face.S[0] || !face.S[1]) {
                return false;
            }
        }
        return true;
    }

    // The triangle data and default triangle creation.
    protected static createTriangle(v0: number, v1: number, v2: number): TSManifoldMeshTriangle {
        return new TSManifoldMeshTriangle(v0, v1, v2);
    }

    // The tetrahedron data and default tetrahedron creation.
    protected static createTetrahedron(v0: number, v1: number, v2: number, v3: number):
        TSManifoldMeshTetrahedron {
        return new TSManifoldMeshTetrahedron(v0, v1, v2, v3);
    }
}
