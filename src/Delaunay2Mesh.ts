// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Delaunay2Mesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A lightweight mesh interface to a Delaunay2 triangulation. The class stores
// a reference to the triangulation; the triangulation must outlive the mesh
// object and must not be recomputed while the mesh object is in use.
//
// Port notes:
// * Upstream Delaunay2Mesh.h declares a variadic class template with two
//   specializations: the deprecated
//   Delaunay2Mesh<InputType, ComputeType, RationalType> and the replacement
//   Delaunay2Mesh<T>. Only the replacement is ported, matching the Delaunay2
//   port.
// * The upstream overloads GetVertices(), GetVertices(t, out),
//   GetIndices(), GetIndices(t, out), GetAdjacencies() and
//   GetAdjacencies(t, out) cannot be distinguished by name in TypeScript, so
//   the per-triangle forms are named getTriangleVertices(t),
//   getTriangleIndices(t) and getTriangleAdjacencies(t), matching the naming
//   already used by Delaunay2. The bool/out-parameter forms return the value
//   or null.
// * The upstream barycentric computation uses BSRational<UIntegerAP32> for
//   exactness; the port does the same with BSRational, then converts the
//   result to double.

import { logAssert } from './Logger';
import { Vector } from './Vector';
import { BSRational } from './BSRational';
import { Delaunay2, Delaunay2SearchInfo } from './Delaunay2';

// The port of Vector2<Rational>.
type RationalPoint2 = [BSRational, BSRational];

// The port of DotPerp(v0, v1) = v0[0]*v1[1] - v0[1]*v1[0] for rational
// vectors.
function rationalDotPerp(v0: RationalPoint2, v1: RationalPoint2): BSRational {
    return v0[0].mul(v1[1]).sub(v0[1].mul(v1[0]));
}

function rationalSub(v0: RationalPoint2, v1: RationalPoint2): RationalPoint2 {
    return [v0[0].sub(v1[0]), v0[1].sub(v1[1])];
}

export class Delaunay2Mesh {
    private mDelaunay: Delaunay2;

    // Construction. The Delaunay triangulation must have dimension 2.
    constructor(delaunay: Delaunay2) {
        logAssert(delaunay.getDimension() === 2, 'Invalid Delaunay dimension.');
        this.mDelaunay = delaunay;
    }

    // Mesh information.
    getNumVertices(): number {
        return this.mDelaunay.getNumVertices();
    }

    getNumTriangles(): number {
        return this.mDelaunay.getNumTriangles();
    }

    getVertices(): readonly Vector[] {
        return this.mDelaunay.getVertices();
    }

    getIndices(): readonly number[] {
        return this.mDelaunay.getIndices();
    }

    getAdjacencies(): readonly number[] {
        return this.mDelaunay.getAdjacencies();
    }

    // The port of Delaunay2::negOne, the "no such triangle" sentinel.
    getInvalidIndex(): number {
        return -1;
    }

    // Containment queries. The returned value is the index of the triangle
    // containing P, or getInvalidIndex() when there is no such triangle.
    getContainingTriangle(p: Vector): number {
        const info = new Delaunay2SearchInfo();
        return this.mDelaunay.getContainingTriangle(p, info);
    }

    // The vertices of triangle t, or null when t is not a valid triangle
    // index.
    getTriangleVertices(t: number): [Vector, Vector, Vector] | null {
        const indices = this.mDelaunay.getTriangleIndices(t);
        if (indices === null) {
            return null;
        }
        const vertices = this.mDelaunay.getVertices();
        return [vertices[indices[0]].clone(), vertices[indices[1]].clone(),
            vertices[indices[2]].clone()];
    }

    // The vertex indices of triangle t, or null when t is not a valid
    // triangle index.
    getTriangleIndices(t: number): [number, number, number] | null {
        return this.mDelaunay.getTriangleIndices(t);
    }

    // The adjacent triangles of triangle t, or null when t is not a valid
    // triangle index. An adjacency of -1 indicates a boundary edge.
    getTriangleAdjacencies(t: number): [number, number, number] | null {
        return this.mDelaunay.getTriangleAdjacencies(t);
    }

    // The barycentric coordinates of P with respect to triangle t. The result
    // is null when t is not a valid triangle index or when the triangle is
    // degenerate. The computation is exact (rational arithmetic) and the
    // result is rounded to double at the end.
    getBarycentrics(t: number, p: Vector): [number, number, number] | null {
        logAssert(p.size === 2, 'Delaunay2Mesh requires 2D points.');
        const indices = this.mDelaunay.getTriangleIndices(t);
        if (indices === null) {
            return null;
        }

        const vertices = this.mDelaunay.getVertices();
        const rtV: RationalPoint2[] = [];
        for (let i = 0; i < 3; ++i) {
            const v = vertices[indices[i]].values;
            rtV.push([BSRational.fromNumber(v[0]), BSRational.fromNumber(v[1])]);
        }
        const rtP: RationalPoint2 =
            [BSRational.fromNumber(p.values[0]), BSRational.fromNumber(p.values[1])];

        // The port of ComputeBarycentrics(P, V0, V1, V2, bary) with the
        // default epsilon of zero.
        const diff0 = rationalSub(rtV[0], rtV[2]);
        const diff1 = rationalSub(rtV[1], rtV[2]);
        const diff2 = rationalSub(rtP, rtV[2]);
        const det = rationalDotPerp(diff0, diff1);
        if (det.getSign() === 0) {
            return null;
        }

        const b0 = rationalDotPerp(diff2, diff1).div(det);
        const b1 = rationalDotPerp(diff0, diff2).div(det);
        const b2 = BSRational.fromNumber(1).sub(b0).sub(b1);
        return [b0.toNumber(), b1.toNumber(), b2.toNumber()];
    }
}
