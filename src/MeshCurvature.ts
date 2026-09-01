// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MeshCurvature.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The MeshCurvature class estimates principal curvatures and principal
// directions at the vertices of a manifold triangle mesh. The algorithm is
// described in
// https://www.geometrictools.com/Documentation/MeshDifferentialGeometry.pdf
//
// Port notes: the two upstream operator() overloads (raw pointers and
// std::vector) collapse into a single compute() that takes arrays, per
// PORTING.md. The results are read back through the get* accessors.

import { Matrix, mulMatrix, multiplyAB } from './Matrix';
import { inverse3x3 } from './Matrix3x3';
import { Vector, add, dot, mul, normalize, sub } from './Vector';
import { cross, computeOrthogonalComplement3 } from './Vector3';

export class MeshCurvature {
    private mNormals: Vector[] = [];
    private mMinCurvatures: number[] = [];
    private mMaxCurvatures: number[] = [];
    private mMinDirections: Vector[] = [];
    private mMaxDirections: Vector[] = [];

    // The input is a triangle mesh with the specified vertex array and index
    // array. The number of elements of 'indices' must be a multiple of 3,
    // each triple of indices (3*t, 3*t+1, 3*t+2) representing the triangle
    // with vertices (vertices[indices[3*t]], vertices[indices[3*t+1]],
    // vertices[indices[3*t+2]]). The singularity threshold is a small
    // nonnegative number. It is used to characterize whether the DWTrn matrix
    // is singular. In theory, set the threshold to zero. In practice you
    // might have to set this to a small positive number.
    compute(vertices: readonly Vector[], indices: ArrayLike<number>,
        singularityThreshold: number): void {
        const numVertices = vertices.length;
        const numTriangles = Math.floor(indices.length / 3);

        this.mNormals = new Array<Vector>(numVertices);
        this.mMinCurvatures = new Array<number>(numVertices).fill(0);
        this.mMaxCurvatures = new Array<number>(numVertices).fill(0);
        this.mMinDirections = new Array<Vector>(numVertices);
        this.mMaxDirections = new Array<Vector>(numVertices);

        // Compute the normal vectors for the vertices as an area-weighted sum
        // of the triangles sharing a vertex.
        for (let i = 0; i < numVertices; ++i) {
            this.mNormals[i] = new Vector(3);
            this.mMinDirections[i] = new Vector(3);
            this.mMaxDirections[i] = new Vector(3);
        }
        let currentIndex = 0;
        for (let i = 0; i < numTriangles; ++i) {
            // Get vertex indices.
            const v0 = indices[currentIndex++];
            const v1 = indices[currentIndex++];
            const v2 = indices[currentIndex++];

            // Compute the normal (length provides a weighted sum).
            const edge1 = sub(vertices[v1], vertices[v0]);
            const edge2 = sub(vertices[v2], vertices[v0]);
            const normal = cross(edge1, edge2);

            this.mNormals[v0] = add(this.mNormals[v0], normal);
            this.mNormals[v1] = add(this.mNormals[v1], normal);
            this.mNormals[v2] = add(this.mNormals[v2], normal);
        }
        for (let i = 0; i < numVertices; ++i) {
            normalize(this.mNormals[i]);
        }

        // Compute the matrix of normal derivatives.
        const DNormal = new Array<Matrix>(numVertices);
        const WWTrn = new Array<Matrix>(numVertices);
        const DWTrn = new Array<Matrix>(numVertices);
        const DWTrnZero = new Array<boolean>(numVertices).fill(false);
        for (let i = 0; i < numVertices; ++i) {
            DNormal[i] = new Matrix(3, 3);
            WWTrn[i] = new Matrix(3, 3);
            DWTrn[i] = new Matrix(3, 3);
        }

        currentIndex = 0;
        for (let i = 0; i < numTriangles; ++i) {
            // Get vertex indices.
            const v: [number, number, number] = [
                indices[currentIndex++],
                indices[currentIndex++],
                indices[currentIndex++]
            ];

            for (let j = 0; j < 3; ++j) {
                const v0 = v[j];
                const v1 = v[(j + 1) % 3];
                const v2 = v[(j + 2) % 3];

                // Compute the edge direction from vertex v0 to vertex v1,
                // project it to the tangent plane of vertex v0 and compute
                // the difference of adjacent normals.
                let E = sub(vertices[v1], vertices[v0]);
                let W = sub(E, mul(this.mNormals[v0], dot(E, this.mNormals[v0])));
                let D = sub(this.mNormals[v1], this.mNormals[v0]);
                for (let row = 0; row < 3; ++row) {
                    for (let col = 0; col < 3; ++col) {
                        WWTrn[v0].set(row, col, WWTrn[v0].get(row, col) +
                            W.values[row] * W.values[col]);
                        DWTrn[v0].set(row, col, DWTrn[v0].get(row, col) +
                            D.values[row] * W.values[col]);
                    }
                }

                // Compute the edge direction from vertex v0 to vertex v2,
                // project it to the tangent plane of vertex v0 and compute
                // the difference of adjacent normals.
                E = sub(vertices[v2], vertices[v0]);
                W = sub(E, mul(this.mNormals[v0], dot(E, this.mNormals[v0])));
                D = sub(this.mNormals[v2], this.mNormals[v0]);
                for (let row = 0; row < 3; ++row) {
                    for (let col = 0; col < 3; ++col) {
                        WWTrn[v0].set(row, col, WWTrn[v0].get(row, col) +
                            W.values[row] * W.values[col]);
                        DWTrn[v0].set(row, col, DWTrn[v0].get(row, col) +
                            D.values[row] * W.values[col]);
                    }
                }
            }
        }

        // Add in N*N^T to W*W^T for numerical stability. In theory 0*0^T is
        // added to D*W^T, but of course no update is needed in the
        // implementation. Compute the matrix of normal derivatives.
        for (let i = 0; i < numVertices; ++i) {
            for (let row = 0; row < 3; ++row) {
                for (let col = 0; col < 3; ++col) {
                    WWTrn[i].set(row, col, 0.5 * WWTrn[i].get(row, col) +
                        this.mNormals[i].values[row] * this.mNormals[i].values[col]);
                    DWTrn[i].set(row, col, 0.5 * DWTrn[i].get(row, col));
                }
            }

            // Compute the max-abs entry of D*W^T. If this entry is (nearly)
            // zero, flag the DNormal matrix as singular.
            let maxAbs = 0;
            for (let row = 0; row < 3; ++row) {
                for (let col = 0; col < 3; ++col) {
                    const absEntry = Math.abs(DWTrn[i].get(row, col));
                    if (absEntry > maxAbs) {
                        maxAbs = absEntry;
                    }
                }
            }
            if (maxAbs < singularityThreshold) {
                DWTrnZero[i] = true;
            }

            DNormal[i] = multiplyAB(DWTrn[i], inverse3x3(WWTrn[i]).inverse);
        }

        // If N is a unit-length normal at a vertex, let U and V be unit-length
        // tangents so that {U, V, N} is an orthonormal set. Define the matrix
        // J = [U | V], a 3-by-2 matrix whose columns are U and V. Define J^T
        // to be the transpose of J, a 2-by-3 matrix. Let dN/dX denote the
        // matrix of first-order derivatives of the normal vector field. The
        // shape matrix is
        //   S = (J^T * J)^{-1} * J^T * dN/dX * J = J^T * dN/dX * J
        // where the superscript of -1 denotes the inverse; the formula allows
        // for J to be created from non-perpendicular vectors. The matrix S is
        // 2-by-2. The principal curvatures are the eigenvalues of S. If k is a
        // principal curvature and W is the 2-by-1 eigenvector corresponding to
        // it, then S*W = k*W (by definition). The corresponding 3-by-1 tangent
        // vector at the vertex is a principal direction for k and is J*W.
        for (let i = 0; i < numVertices; ++i) {
            // Compute U and V given N.
            // Upstream copies the normal into basis[0] by value; the port clones
            // it because computeOrthogonalComplement3 normalizes in place.
            const basis: Vector[] = [this.mNormals[i].clone(), new Vector(3), new Vector(3)];
            computeOrthogonalComplement3(1, basis);
            const U = basis[1];
            const V = basis[2];

            if (DWTrnZero[i]) {
                // At a locally planar point.
                this.mMinCurvatures[i] = 0;
                this.mMaxCurvatures[i] = 0;
                this.mMinDirections[i] = U;
                this.mMaxDirections[i] = V;
                continue;
            }

            // Compute S = J^T * dN/dX * J. In theory S is symmetric, but
            // because dN/dX is estimated, we must ensure that the computed S
            // is symmetric.
            const s00 = dot(U, mulMatrix(DNormal[i], U) as Vector);
            const s01 = dot(U, mulMatrix(DNormal[i], V) as Vector);
            const s10 = dot(V, mulMatrix(DNormal[i], U) as Vector);
            const s11 = dot(V, mulMatrix(DNormal[i], V) as Vector);
            const avr = 0.5 * (s01 + s10);
            const S = Matrix.fromArray(2, 2, [s00, avr, avr, s11]);

            // Compute the eigenvalues of S (min and max curvatures).
            const trace = S.get(0, 0) + S.get(1, 1);
            const det = S.get(0, 0) * S.get(1, 1) - S.get(0, 1) * S.get(1, 0);
            const discr = trace * trace - 4.0 * det;
            const rootDiscr = Math.sqrt(Math.max(discr, 0));
            this.mMinCurvatures[i] = 0.5 * (trace - rootDiscr);
            this.mMaxCurvatures[i] = 0.5 * (trace + rootDiscr);

            // Compute the eigenvectors of S.
            let W0 = Vector.fromArray([S.get(0, 1), this.mMinCurvatures[i] - S.get(0, 0)]);
            let W1 = Vector.fromArray([this.mMinCurvatures[i] - S.get(1, 1), S.get(1, 0)]);
            if (dot(W0, W0) >= dot(W1, W1)) {
                normalize(W0);
                this.mMinDirections[i] = add(mul(U, W0.values[0]), mul(V, W0.values[1]));
            }
            else {
                normalize(W1);
                this.mMinDirections[i] = add(mul(U, W1.values[0]), mul(V, W1.values[1]));
            }

            W0 = Vector.fromArray([S.get(0, 1), this.mMaxCurvatures[i] - S.get(0, 0)]);
            W1 = Vector.fromArray([this.mMaxCurvatures[i] - S.get(1, 1), S.get(1, 0)]);
            if (dot(W0, W0) >= dot(W1, W1)) {
                normalize(W0);
                this.mMaxDirections[i] = add(mul(U, W0.values[0]), mul(V, W0.values[1]));
            }
            else {
                normalize(W1);
                this.mMaxDirections[i] = add(mul(U, W1.values[0]), mul(V, W1.values[1]));
            }
        }
    }

    getNormals(): readonly Vector[] {
        return this.mNormals;
    }

    getMinCurvatures(): readonly number[] {
        return this.mMinCurvatures;
    }

    getMaxCurvatures(): readonly number[] {
        return this.mMaxCurvatures;
    }

    getMinDirections(): readonly Vector[] {
        return this.mMinDirections;
    }

    getMaxDirections(): readonly Vector[] {
        return this.mMaxDirections;
    }
}
