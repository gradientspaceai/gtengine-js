// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpLinearNonuniform2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Linear interpolation of a network of triangles whose vertices are of the
// form (x,y,f(x,y)). The function samples are F[i] and represent
// f(x[i],y[i]), where i is the index of the input vertex (x[i],y[i]) to
// Delaunay2.
//
// Port notes: upstream the mesh is a template parameter constrained only by
// a duck-typed interface, so the port declares that interface explicitly as
// IntpLinearNonuniform2TriangleMesh. Any object satisfying it (a Delaunay2
// triangulation, a planar mesh, or a hand-built mesh) can be used. The C++
// methods that return 'bool' and write to a reference parameter become
// methods that return the value or null:
//   bool GetIndices(int32_t, std::array<int32_t, 3>&) -> getIndices
//   bool GetBarycentrics(int32_t, Vector2 const&, std::array<Real, 3>&)
//       -> getBarycentrics
// 'operator()(P, F&) -> bool' becomes evaluate(P) returning
// { valid, F }; when valid is false the F value is meaningless (upstream
// leaves the caller's F untouched in that case). Upstream ignores the
// 'bool' returned by GetIndices and would then use the value-initialized
// index triple (0,0,0); the port reports the failure as an invalid result
// instead of interpolating from garbage indices.

import { logAssert } from './Logger.js';
import { Vector } from './Vector.js';

export interface IntpLinearNonuniform2TriangleMesh {
    // The index of the triangle containing P, or -1 when P is outside the
    // mesh.
    getContainingTriangle(P: Vector): number;

    // The three vertex indices of triangle t, or null on failure.
    getIndices(t: number): readonly number[] | null;

    // The barycentric coordinates of P with respect to triangle t, or null
    // when the triangle is degenerate.
    getBarycentrics(t: number, P: Vector): readonly number[] | null;
}

export interface IntpLinearNonuniform2Result {
    // Valid is true if and only if the input point P is in the convex hull
    // of the input vertices, in which case the interpolation is valid.
    valid: boolean;

    // The interpolated function value; meaningful only when valid is true.
    F: number;
}

export class IntpLinearNonuniform2 {
    private mMesh: IntpLinearNonuniform2TriangleMesh;
    private mF: readonly number[];

    // Construction. The mesh and the samples F are aliased, not copied.
    constructor(mesh: IntpLinearNonuniform2TriangleMesh, F: readonly number[]) {
        logAssert(F.length > 0, 'Invalid input.');
        this.mMesh = mesh;
        this.mF = F;
    }

    // Linear interpolation.
    evaluate(P: Vector): IntpLinearNonuniform2Result {
        const t = this.mMesh.getContainingTriangle(P);
        if (t === -1) {
            // The point is outside the triangulation.
            return { valid: false, F: 0 };
        }

        // Get the barycentric coordinates of P with respect to the triangle,
        // P = b0*V0 + b1*V1 + b2*V2, where b0 + b1 + b2 = 1.
        const bary = this.mMesh.getBarycentrics(t, P);
        if (bary === null) {
            // P is in a needle-like or degenerate triangle.
            return { valid: false, F: 0 };
        }

        // The result is a barycentric combination of function values.
        const indices = this.mMesh.getIndices(t);
        if (indices === null) {
            return { valid: false, F: 0 };
        }

        const F = bary[0] * this.mF[indices[0]] + bary[1] * this.mF[indices[1]]
            + bary[2] * this.mF[indices[2]];
        return { valid: true, F };
    }
}
