// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) PrimalQuery3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Queries about the relation of a point to various geometric objects in 3D.
//
// Port notes:
// * Upstream is templated on Real and is intended to be instantiated with an
//   exact arithmetic type (BSNumber or BSRational) so that the signs of the
//   computed determinants are exact and the classifications are therefore
//   robust. This port computes with 'number' (IEEE double), which makes the
//   queries exact only when the inputs and all intermediate products and sums
//   are exactly representable in double precision (for example, small
//   integer-valued coordinates). For inexact inputs, a determinant that is
//   theoretically zero can be computed as a small nonzero value and the query
//   then misclassifies a degenerate configuration. The order of the
//   arithmetic operations is preserved exactly as upstream so that a future
//   exact-arithmetic instantiation reproduces upstream results term by term.
// * The upstream per-query comment tables "Choice of N for UIntegerFP32<N>"
//   are not reproduced. They describe the precision requirements of the
//   BSNumber/BSRational instantiations, and some of them are stale in the
//   upstream headers (see gtengine-js issue #43); they are also irrelevant to
//   the 'number' instantiation ported here.
// * Upstream stores a raw pointer to the vertex array; the port stores the
//   array reference. As upstream, the class does no range checking.
// * The two overloads of each query, one taking a vertex index i and one
//   taking a point 'test', are merged into a single method whose first
//   parameter is 'number | Vector'.

import { Vector } from './Vector.js';

export class PrimalQuery3 {
    private mNumVertices: number;
    private mVertices: readonly Vector[];

    // The caller is responsible for ensuring that the array is not empty
    // before calling queries and that the indices passed to the queries are
    // valid. The class does no range checking. The port of the upstream
    // default constructor is 'new PrimalQuery3()'.
    constructor(numVertices: number = 0, vertices: readonly Vector[] = []) {
        this.mNumVertices = numVertices;
        this.mVertices = vertices;
    }

    // Member access.
    set(numVertices: number, vertices: readonly Vector[]): void {
        this.mNumVertices = numVertices;
        this.mVertices = vertices;
    }

    getNumVertices(): number {
        return this.mNumVertices;
    }

    getVertices(): readonly Vector[] {
        return this.mVertices;
    }

    // In the following, point P refers to vertices[i] or 'test' and Vi refers
    // to vertices[vi].

    // For a plane with origin V0 and normal N = Cross(V1-V0,V2-V0), toPlane
    // returns
    //   +1, P on positive side of plane (side to which N points)
    //   -1, P on negative side of plane (side to which -N points)
    //    0, P on the plane
    toPlane(test: number | Vector, v0: number, v1: number, v2: number): number {
        const p = this.point(test);
        const vec0 = this.mVertices[v0];
        const vec1 = this.mVertices[v1];
        const vec2 = this.mVertices[v2];

        const x0 = p.values[0] - vec0.values[0];
        const y0 = p.values[1] - vec0.values[1];
        const z0 = p.values[2] - vec0.values[2];
        const x1 = vec1.values[0] - vec0.values[0];
        const y1 = vec1.values[1] - vec0.values[1];
        const z1 = vec1.values[2] - vec0.values[2];
        const x2 = vec2.values[0] - vec0.values[0];
        const y2 = vec2.values[1] - vec0.values[1];
        const z2 = vec2.values[2] - vec0.values[2];
        const y1z2 = y1 * z2;
        const y2z1 = y2 * z1;
        const y2z0 = y2 * z0;
        const y0z2 = y0 * z2;
        const y0z1 = y0 * z1;
        const y1z0 = y1 * z0;
        const c0 = y1z2 - y2z1;
        const c1 = y2z0 - y0z2;
        const c2 = y0z1 - y1z0;
        const x0c0 = x0 * c0;
        const x1c1 = x1 * c1;
        const x2c2 = x2 * c2;
        const term = x0c0 + x1c1;
        const det = term + x2c2;
        const zero = 0;

        return (det > zero ? +1 : (det < zero ? -1 : 0));
    }

    // For a tetrahedron with vertices ordered as described in the file
    // TetrahedronKey.h, the function returns
    //   +1, P outside tetrahedron
    //   -1, P inside tetrahedron
    //    0, P on tetrahedron
    toTetrahedron(test: number | Vector, v0: number, v1: number, v2: number,
        v3: number): number {
        const p = this.point(test);

        const sign0 = this.toPlane(p, v1, v2, v3);
        if (sign0 > 0) {
            return +1;
        }

        const sign1 = this.toPlane(p, v0, v2, v3);
        if (sign1 < 0) {
            return +1;
        }

        const sign2 = this.toPlane(p, v0, v1, v3);
        if (sign2 > 0) {
            return +1;
        }

        const sign3 = this.toPlane(p, v0, v1, v2);
        if (sign3 < 0) {
            return +1;
        }

        return ((sign0 !== 0 && sign1 !== 0 && sign2 !== 0 && sign3 !== 0) ? -1 : 0);
    }

    // For a tetrahedron with vertices ordered as described in the file
    // TetrahedronKey.h, the function returns
    //   +1, P outside circumsphere of tetrahedron
    //   -1, P inside circumsphere of tetrahedron
    //    0, P on circumsphere of tetrahedron
    toCircumsphere(test: number | Vector, v0: number, v1: number, v2: number,
        v3: number): number {
        const p = this.point(test);
        const vec0 = this.mVertices[v0];
        const vec1 = this.mVertices[v1];
        const vec2 = this.mVertices[v2];
        const vec3 = this.mVertices[v3];

        const x0 = vec0.values[0] - p.values[0];
        const y0 = vec0.values[1] - p.values[1];
        const z0 = vec0.values[2] - p.values[2];
        const s00 = vec0.values[0] + p.values[0];
        const s01 = vec0.values[1] + p.values[1];
        const s02 = vec0.values[2] + p.values[2];
        const t00 = s00 * x0;
        const t01 = s01 * y0;
        const t02 = s02 * z0;
        const t00pt01 = t00 + t01;
        const w0 = t00pt01 + t02;

        const x1 = vec1.values[0] - p.values[0];
        const y1 = vec1.values[1] - p.values[1];
        const z1 = vec1.values[2] - p.values[2];
        const s10 = vec1.values[0] + p.values[0];
        const s11 = vec1.values[1] + p.values[1];
        const s12 = vec1.values[2] + p.values[2];
        const t10 = s10 * x1;
        const t11 = s11 * y1;
        const t12 = s12 * z1;
        const t10pt11 = t10 + t11;
        const w1 = t10pt11 + t12;

        const x2 = vec2.values[0] - p.values[0];
        const y2 = vec2.values[1] - p.values[1];
        const z2 = vec2.values[2] - p.values[2];
        const s20 = vec2.values[0] + p.values[0];
        const s21 = vec2.values[1] + p.values[1];
        const s22 = vec2.values[2] + p.values[2];
        const t20 = s20 * x2;
        const t21 = s21 * y2;
        const t22 = s22 * z2;
        const t20pt21 = t20 + t21;
        const w2 = t20pt21 + t22;

        const x3 = vec3.values[0] - p.values[0];
        const y3 = vec3.values[1] - p.values[1];
        const z3 = vec3.values[2] - p.values[2];
        const s30 = vec3.values[0] + p.values[0];
        const s31 = vec3.values[1] + p.values[1];
        const s32 = vec3.values[2] + p.values[2];
        const t30 = s30 * x3;
        const t31 = s31 * y3;
        const t32 = s32 * z3;
        const t30pt31 = t30 + t31;
        const w3 = t30pt31 + t32;

        const x0y1 = x0 * y1;
        const x0y2 = x0 * y2;
        const x0y3 = x0 * y3;
        const x1y0 = x1 * y0;
        const x1y2 = x1 * y2;
        const x1y3 = x1 * y3;
        const x2y0 = x2 * y0;
        const x2y1 = x2 * y1;
        const x2y3 = x2 * y3;
        const x3y0 = x3 * y0;
        const x3y1 = x3 * y1;
        const x3y2 = x3 * y2;
        const a0 = x0y1 - x1y0;
        const a1 = x0y2 - x2y0;
        const a2 = x0y3 - x3y0;
        const a3 = x1y2 - x2y1;
        const a4 = x1y3 - x3y1;
        const a5 = x2y3 - x3y2;

        const z0w1 = z0 * w1;
        const z0w2 = z0 * w2;
        const z0w3 = z0 * w3;
        const z1w0 = z1 * w0;
        const z1w2 = z1 * w2;
        const z1w3 = z1 * w3;
        const z2w0 = z2 * w0;
        const z2w1 = z2 * w1;
        const z2w3 = z2 * w3;
        const z3w0 = z3 * w0;
        const z3w1 = z3 * w1;
        const z3w2 = z3 * w2;
        const b0 = z0w1 - z1w0;
        const b1 = z0w2 - z2w0;
        const b2 = z0w3 - z3w0;
        const b3 = z1w2 - z2w1;
        const b4 = z1w3 - z3w1;
        const b5 = z2w3 - z3w2;
        const a0b5 = a0 * b5;
        const a1b4 = a1 * b4;
        const a2b3 = a2 * b3;
        const a3b2 = a3 * b2;
        const a4b1 = a4 * b1;
        const a5b0 = a5 * b0;
        const term0 = a0b5 - a1b4;
        const term1 = term0 + a2b3;
        const term2 = term1 + a3b2;
        const term3 = term2 - a4b1;
        const det = term3 + a5b0;
        const zero = 0;

        return (det > zero ? 1 : (det < zero ? -1 : 0));
    }

    // The port of the upstream overload pairs: a query taking an index i uses
    // mVertices[i] as the test point.
    private point(test: number | Vector): Vector {
        return (typeof test === 'number' ? this.mVertices[test] : test);
    }
}
