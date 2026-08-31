// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MeshSmoother.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Laplacian-style smoothing of a triangle mesh. Each vertex is moved toward
// the average of its neighbors, with the displacement split into the part
// tangent to the estimated surface normal and the part along the normal.
// The tangential motion smooths the mesh (a discrete approximation to mean
// curvature flow when the normal weight is zero) and the normal motion can
// be used to inflate or deflate the surface. The default weights are 1/2 for
// the tangent and 0 for the normal.
//
// Derive from this class and override vertexInfluenced, getTangentWeight and
// getNormalWeight for time-dependent or vertex-dependent evolutions.
//
// Port notes:
//   - The two upstream operator() overloads (pointer/count form and
//     std::vector form) collapse into the single method
//     initialize(vertices, indices). The vertices are Vector objects of size
//     3; they are aliased, not copied, and update() mutates them in place
//     exactly as the upstream pointer does.
//   - The upstream mean is 'mMeans[i] /= neighborCounts[i]'. A vertex that
//     is not referenced by any triangle has neighbor count 0, so upstream
//     divides 0/0 and gets NaN; the port reproduces that (the vertex is
//     still moved by update()). Callers should pass meshes in which every
//     vertex belongs to a triangle.

import { logAssert } from './Logger';
import { Vector, add, dot, mul, normalize, sub } from './Vector';
import { cross } from './Vector3';

export class MeshSmoother {
    protected mNumVertices: number;
    protected mVertices: Vector[];
    protected mNumTriangles: number;
    protected mIndices: ArrayLike<number>;

    protected mNormals: Vector[];
    protected mMeans: Vector[];
    protected mNeighborCounts: number[];

    constructor() {
        this.mNumVertices = 0;
        this.mVertices = [];
        this.mNumTriangles = 0;
        this.mIndices = [];
        this.mNormals = [];
        this.mMeans = [];
        this.mNeighborCounts = [];
    }

    // The input is a triangle mesh with the specified vertex buffer and
    // index buffer. The number of elements of 'indices' must be a multiple
    // of 3, each triple of indices (3*t, 3*t+1, 3*t+2) representing a
    // triangle. The vertices are modified in place by update().
    initialize(vertices: Vector[], indices: ArrayLike<number>): void {
        const numVertices = vertices.length;
        const numTriangles = Math.floor(indices.length / 3);
        logAssert(numVertices >= 3 && numTriangles >= 1, 'Invalid input.');

        this.mNumVertices = numVertices;
        this.mVertices = vertices;
        this.mNumTriangles = numTriangles;
        this.mIndices = indices;

        this.mNormals = new Array<Vector>(numVertices);
        this.mMeans = new Array<Vector>(numVertices);
        this.mNeighborCounts = new Array<number>(numVertices).fill(0);
        for (let i = 0; i < numVertices; ++i) {
            this.mNormals[i] = new Vector(3);
            this.mMeans[i] = new Vector(3);
        }

        // Count the number of vertex neighbors.
        let current = 0;
        for (let i = 0; i < numTriangles; ++i) {
            this.mNeighborCounts[indices[current++]] += 2;
            this.mNeighborCounts[indices[current++]] += 2;
            this.mNeighborCounts[indices[current++]] += 2;
        }
    }

    getNumVertices(): number {
        return this.mNumVertices;
    }

    getVertices(): Vector[] {
        return this.mVertices;
    }

    getNumTriangles(): number {
        return this.mNumTriangles;
    }

    getIndices(): ArrayLike<number> {
        return this.mIndices;
    }

    getNormals(): Vector[] {
        return this.mNormals;
    }

    getMeans(): Vector[] {
        return this.mMeans;
    }

    getNeighborCounts(): number[] {
        return this.mNeighborCounts;
    }

    // Apply one iteration of the smoother. The input time is supported for
    // applications where the surface evolution is time-dependent.
    update(t: number = 0): void {
        for (let i = 0; i < this.mNumVertices; ++i) {
            this.mNormals[i].values.fill(0);
            this.mMeans[i].values.fill(0);
        }

        let current = 0;
        for (let i = 0; i < this.mNumTriangles; ++i) {
            const v0 = this.mIndices[current++];
            const v1 = this.mIndices[current++];
            const v2 = this.mIndices[current++];

            const V0 = this.mVertices[v0];
            const V1 = this.mVertices[v1];
            const V2 = this.mVertices[v2];

            const edge1 = sub(V1, V0);
            const edge2 = sub(V2, V0);
            const normal = cross(edge1, edge2);

            addInPlace(this.mNormals[v0], normal);
            addInPlace(this.mNormals[v1], normal);
            addInPlace(this.mNormals[v2], normal);

            addInPlace(this.mMeans[v0], add(V1, V2));
            addInPlace(this.mMeans[v1], add(V2, V0));
            addInPlace(this.mMeans[v2], add(V0, V1));
        }

        for (let i = 0; i < this.mNumVertices; ++i) {
            normalize(this.mNormals[i]);
            const denom = this.mNeighborCounts[i];
            const mean = this.mMeans[i].values;
            for (let j = 0; j < 3; ++j) {
                mean[j] /= denom;
            }
        }

        for (let i = 0; i < this.mNumVertices; ++i) {
            if (this.vertexInfluenced(i, t)) {
                const diff = sub(this.mMeans[i], this.mVertices[i]);
                const dotDifNor = dot(diff, this.mNormals[i]);
                const surfaceNormal = mul(this.mNormals[i], dotDifNor);
                const tangent = sub(diff, surfaceNormal);

                const tanWeight = this.getTangentWeight(i, t);
                const norWeight = this.getNormalWeight(i, t);
                addInPlace(this.mVertices[i],
                    add(mul(tangent, tanWeight), mul(this.mNormals[i], norWeight)));
            }
        }
    }

    // The input parameters are the vertex index and the time. They are
    // unused by the default implementations.
    protected vertexInfluenced(_i: number, _t: number): boolean {
        return true;
    }

    protected getTangentWeight(_i: number, _t: number): number {
        return 0.5;
    }

    protected getNormalWeight(_i: number, _t: number): number {
        return 0.0;
    }
}

// The port of C++ 'v0 += v1' for vectors, which has no in-place module
// function in src/Vector.ts.
function addInPlace(v0: Vector, v1: Vector): void {
    for (let i = 0; i < v0.values.length; ++i) {
        v0.values[i] += v1.values[i];
    }
}
