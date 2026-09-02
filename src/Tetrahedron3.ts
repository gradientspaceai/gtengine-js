// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Tetrahedron3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The tetrahedron is represented as an array of four vertices, V[i] for
// 0 <= i <= 3. The vertices are ordered so that the triangular faces are
// counterclockwise-ordered triangles when viewed by an observer outside the
// tetrahedron: face 0 = <V[0],V[2],V[1]>, face 1 = <V[0],V[1],V[3]>,
// face 2 = <V[0],V[3],V[2]> and face 3 = <V[1],V[2],V[3]>. The canonical
// tetrahedron has V[0] = (0,0,0), V[1] = (1,0,0), V[2] = (0,1,0) and
// V[3] = (0,0,1).
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (named static factories that copy their arguments, comparison
// methods). The upstream static-local index tables become module-level frozen
// arrays; the accessors return them directly, as upstream returns const
// references, so callers must not mutate them.

import { logAssert } from './Logger';
import { Hyperplane } from './Hyperplane';
import { Vector, add, dot, mul, negate, sub } from './Vector';
import { unitCross } from './Vector3';

// The vertex indices of the four faces, ordered counterclockwise when viewed
// from outside the tetrahedron.
const sFaceIndices: readonly (readonly number[])[] = Object.freeze([
    Object.freeze([0, 2, 1]),
    Object.freeze([0, 1, 3]),
    Object.freeze([0, 3, 2]),
    Object.freeze([1, 2, 3])
]);

const sAllFaceIndices: readonly number[] = Object.freeze([
    0, 2, 1,
    0, 1, 3,
    0, 3, 2,
    1, 2, 3
]);

// The vertex indices of the six edges.
const sEdgeIndices: readonly (readonly number[])[] = Object.freeze([
    Object.freeze([0, 1]),
    Object.freeze([0, 2]),
    Object.freeze([0, 3]),
    Object.freeze([1, 2]),
    Object.freeze([1, 3]),
    Object.freeze([2, 3])
]);

const sAllEdgeIndices: readonly number[] = Object.freeze([
    0, 1, 0, 2, 0, 3, 1, 2, 1, 3, 2, 3
]);

const sEdgeAugmented: readonly (readonly number[])[] = Object.freeze([
    Object.freeze([0, 1, 2, 3]),
    Object.freeze([0, 2, 3, 1]),
    Object.freeze([0, 3, 1, 2]),
    Object.freeze([1, 2, 0, 3]),
    Object.freeze([1, 3, 2, 0]),
    Object.freeze([2, 3, 0, 1])
]);

const sVertexAugmented: readonly (readonly number[])[] = Object.freeze([
    Object.freeze([0, 1, 3, 2]),
    Object.freeze([1, 3, 0, 2]),
    Object.freeze([2, 1, 0, 3]),
    Object.freeze([3, 2, 0, 1])
]);

export class Tetrahedron3 {
    // Public member access.
    v: Vector[];

    // The port of the default constructor, which sets the vertices to
    // (0,0,0), (1,0,0), (0,1,0) and (0,0,1).
    constructor() {
        this.v = [
            new Vector(3),
            Vector.unit(3, 0),
            Vector.unit(3, 1),
            Vector.unit(3, 2)
        ];
    }

    // The port of 'Tetrahedron3(v0, v1, v2, v3)'. The vectors are copied,
    // matching C++ value semantics.
    static fromVertices(v0: Vector, v1: Vector, v2: Vector,
        v3: Vector): Tetrahedron3 {
        return Tetrahedron3.fromArray([v0, v1, v2, v3]);
    }

    // The port of 'Tetrahedron3(std::array<Vector3<T>, 4> const& inV)'.
    static fromArray(inV: readonly Vector[]): Tetrahedron3 {
        logAssert(inV.length === 4, 'Tetrahedron3: four vertices are required.');
        for (let i = 0; i < 4; ++i) {
            logAssert(inV[i].size === 3, 'Tetrahedron3: mismatched sizes.');
        }
        const tetrahedron = new Tetrahedron3();
        tetrahedron.v = [inV[0].clone(), inV[1].clone(), inV[2].clone(),
            inV[3].clone()];
        return tetrahedron;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Tetrahedron3 {
        return Tetrahedron3.fromArray(this.v);
    }

    // Get the vertex indices for the specified face. The input 'face' must be
    // in {0,1,2,3}. The returned array must not be modified.
    static getFaceIndices(face: number): readonly number[] {
        logAssert(0 <= face && face < 4, 'Invalid face.');
        return sFaceIndices[face];
    }

    static getAllFaceIndices(): readonly number[] {
        return sAllFaceIndices;
    }

    // Get the vertex indices for the specified edge. The input 'edge' must be
    // in {0,1,2,3,4,5}. The returned array must not be modified.
    static getEdgeIndices(edge: number): readonly number[] {
        logAssert(0 <= edge && edge < 6, 'Invalid edge.');
        return sEdgeIndices[edge];
    }

    static getAllEdgeIndices(): readonly number[] {
        return sAllEdgeIndices;
    }

    // Get the vertex indices for the edges with the appropriately ordered
    // adjacent indices. The input 'edge' must be in {0,1,2,3,4,5}. The output
    // is {v0,v1,v2,v3} where the edge is {v0,v1}. The triangles sharing the
    // edge are {v0,v2,v1} and {v0,v1,v3}.
    static getEdgeAugmented(edge: number): readonly number[] {
        logAssert(0 <= edge && edge < 6, 'Invalid edge.');
        return sEdgeAugmented[edge];
    }

    // Get the augmented indices for the vertices with the appropriately
    // ordered adjacent indices. The input 'vertex' must be in {0,1,2,3}. The
    // output is {v0,v1,v2,v3} where the vertex is v0. The triangles sharing
    // the vertex are {v0,v1,v2}, {v0,v2,v3} and {v0,v3,v1}.
    static getVertexAugmented(vertex: number): readonly number[] {
        logAssert(0 <= vertex && vertex < 4, 'Invalid vertex.');
        return sVertexAugmented[vertex];
    }

    // Compute a face normal. The input 'face' must be in {0,1,2,3} and
    // correspond to faces {{0,2,1},{0,1,3},{0,3,2},{1,2,3}}.
    computeFaceNormal(face: number): Vector {
        // Compute the normal for face <v0,v1,v2>.
        const indices = Tetrahedron3.getFaceIndices(face);
        const edge10 = sub(this.v[indices[1]], this.v[indices[0]]);
        const edge20 = sub(this.v[indices[2]], this.v[indices[0]]);
        return unitCross(edge10, edge20);
    }

    // Compute an edge normal, an average of the normals of the 2 faces
    // sharing the edge. The input 'edge' must be in {0,1,2,3,4,5} and
    // correspond to edges {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}}.
    computeEdgeNormal(edge: number): Vector {
        // Compute the weighted average of normals for faces <v0,a0,v1> and
        // <v0,v1,a1> shared by edge <v0,v1>. In the comments,
        // E10 = V[v1]-V[v0], E20 = V[v2]-V[v0], E30 = V[v3]-V[v0] and
        // E23 = V[i2]-V[i3]. The unnormalized vector is
        //   N = E20 x E10 + E10 x E30
        //     = E20 x E10 - E30 x E10
        //     = (E20 - E30) x E10
        //     = E23 x E10
        const indices = Tetrahedron3.getEdgeAugmented(edge);
        const edge23 = sub(this.v[indices[2]], this.v[indices[3]]);
        const edge10 = sub(this.v[indices[1]], this.v[indices[0]]);
        return unitCross(edge23, edge10);
    }

    // Compute a vertex normal, an average of the normals of the 3 faces
    // sharing the vertex. The input 'vertex' must be in {0,1,2,3} and are the
    // indices into the tetrahedron vertex array. The algebra shows that the
    // vertex normal is the negative normal of the face opposite the vertex.
    computeVertexNormal(vertex: number): Vector {
        // Compute the weighted average of normals for faces <v0,v1,v2>,
        // <v0,v2,v3> and <v0,v3,v1>. In the comments, E10 = V[v1]-V[v0],
        // E20 = V[v2]-V[v0], E30 = V[v3]-V[v0], E12 = V[v1]-V[v2],
        // E21 = V[v2]-V[v1] and E31 = V[v3]-V[v1]. The unnormalized vector is
        //   N = E10 x E20 + E20 x E30 + E30 x E10
        //     = E10 x E20 - E30 x E20 + E30 x E10 - E10 x E10
        //     = E13 x E20 + E31 x E10
        //     = E13 x E20 - E13 x E10
        //     = E13 x E21
        const indices = Tetrahedron3.getVertexAugmented(vertex);
        const edge13 = sub(this.v[indices[1]], this.v[indices[3]]);
        const edge21 = sub(this.v[indices[2]], this.v[indices[1]]);
        return unitCross(edge13, edge21);
    }

    // Construct the planes of the faces. The planes have outer pointing
    // normal vectors. The plane indexing is the same as the face indexing
    // mentioned previously. The upstream output parameter becomes the
    // returned array.
    getPlanes(): Hyperplane[] {
        const edge10 = sub(this.v[1], this.v[0]);
        const edge20 = sub(this.v[2], this.v[0]);
        const edge30 = sub(this.v[3], this.v[0]);
        const edge21 = sub(this.v[2], this.v[1]);
        const edge31 = sub(this.v[3], this.v[1]);

        const normals: Vector[] = [
            unitCross(edge20, edge10),  // <v0,v2,v1>
            unitCross(edge10, edge30),  // <v0,v1,v3>
            unitCross(edge30, edge20),  // <v0,v3,v2>
            unitCross(edge21, edge31)   // <v1,v2,v3>
        ];

        const det = dot(edge10, normals[3]);
        if (det < 0) {
            // The normals are inner pointing, reverse their directions.
            for (let i = 0; i < 4; ++i) {
                normals[i] = negate(normals[i]);
            }
        }

        // Upstream assigns only 'normal' and 'constant', leaving the 'origin'
        // member of each Plane3 at its default (0,0,0), which is generally not
        // a point of the plane. The port uses Hyperplane.fromNormalConstant so
        // that all three members are mutually consistent; 'normal' and
        // 'constant' match upstream exactly.
        const plane: Hyperplane[] = new Array<Hyperplane>(4);
        for (let i = 0; i < 4; ++i) {
            plane[i] = Hyperplane.fromNormalConstant(normals[i],
                dot(this.v[i], normals[i]));
        }
        return plane;
    }

    computeCentroid(): Vector {
        return mul(add(add(this.v[0], this.v[1]), add(this.v[2], this.v[3])),
            0.25);
    }

    // Comparisons to support sorted containers. These are the ports of the
    // lexicographic std::array comparisons of the vertices.
    equals(tetrahedron: Tetrahedron3): boolean {
        for (let i = 0; i < 4; ++i) {
            if (!this.v[i].equals(tetrahedron.v[i])) {
                return false;
            }
        }
        return true;
    }

    notEquals(tetrahedron: Tetrahedron3): boolean {
        return !this.equals(tetrahedron);
    }

    lessThan(tetrahedron: Tetrahedron3): boolean {
        for (let i = 0; i < 4; ++i) {
            if (this.v[i].lessThan(tetrahedron.v[i])) {
                return true;
            }
            if (this.v[i].greaterThan(tetrahedron.v[i])) {
                return false;
            }
        }
        return false;
    }

    lessThanOrEqual(tetrahedron: Tetrahedron3): boolean {
        return !tetrahedron.lessThan(this);
    }

    greaterThan(tetrahedron: Tetrahedron3): boolean {
        return tetrahedron.lessThan(this);
    }

    greaterThanOrEqual(tetrahedron: Tetrahedron3): boolean {
        return !this.lessThan(tetrahedron);
    }
}
