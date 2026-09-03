// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) PrimalQuery2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Queries about the relation of a point to various geometric objects in 2D.
//
// Port notes:
// * Upstream is templated on Real and is intended to be instantiated with an
//   exact arithmetic type (BSNumber or BSRational) so that the signs of the
//   computed determinants are exact and the classifications are therefore
//   robust. This port computes with 'number' (IEEE double), which makes the
//   queries exact only when the inputs and all intermediate products and sums
//   are exactly representable in double precision (for example, small
//   integer-valued coordinates). For inexact inputs, a determinant that is
//   theoretically zero can be computed as a small nonzero value and the
//   query then misclassifies a degenerate configuration. The order of the
//   arithmetic operations is preserved exactly as upstream so that a future
//   exact-arithmetic instantiation reproduces upstream results term by term.
// * The upstream per-query comment tables "Choice of N for UIntegerFP32<N>"
//   are not reproduced. They describe the precision requirements of the
//   BSNumber/BSRational instantiations, and some of them are stale in the
//   upstream headers (see gtengine-js issue #43); they are also irrelevant to
//   the 'number' instantiation ported here.
// * Upstream stores a raw pointer to the vertex array; the port stores the
//   array reference. As upstream, the class does no range checking.
// * The ToLine overload with the 'int32_t& order' output parameter becomes
//   toLineWithOrder, which returns an object with the fields 'sign' (the
//   upstream return value) and 'order'. See the comment on
//   PrimalQuery2ToLineOrderResult for an upstream bug in the collinear
//   ordering that the port preserves.
// * The two overloads of each query, one taking a vertex index i and one
//   taking a point 'test', are merged into a single method whose first
//   parameter is 'number | Vector'.
// * The nested enum class OrderType becomes the exported enum
//   PrimalQuery2OrderType (the flat library export requires unique names).

import { Vector } from './Vector.js';

// An extended classification of the relationship of a point to a line
// segment. For noncollinear points, the value is
//   POSITIVE when <P,Q0,Q1> is a counterclockwise triangle
//   NEGATIVE when <P,Q0,Q1> is a clockwise triangle
// For collinear points, the line direction is Q1-Q0. The value is
//   COLLINEAR_LEFT when the line ordering is <P,Q0,Q1>
//   COLLINEAR_RIGHT when the line ordering is <Q0,Q1,P>
//   COLLINEAR_CONTAIN when the line ordering is <Q0,P,Q1>
export enum PrimalQuery2OrderType {
    Q0_EQUALS_Q1,
    P_EQUALS_Q0,
    P_EQUALS_Q1,
    POSITIVE,
    NEGATIVE,
    COLLINEAR_LEFT,
    COLLINEAR_RIGHT,
    COLLINEAR_CONTAIN
}

// The result of toLineWithOrder.
export interface PrimalQuery2ToLineOrderResult {
    // +1 if P is on the right of the line, -1 if P is on the left of the
    // line, 0 if P is on the line.
    sign: number;

    // The relative order of P with respect to the line of <V0,V1>:
    //   +3, points not collinear, P on right of line
    //   -3, points not collinear, P on left of line
    //   -2, P strictly left of V0 on the line
    //   -1, P = V0
    //    0, P interior to line segment [V0,V1]
    //   +1, P = V1
    //   +2, P strictly right of V0 on the line
    //
    // Upstream bug, preserved: for the collinear cases the upstream code
    // compares dot = Dot(P-V0,V1-V0) against sqrLength = |P-V0|^2 rather than
    // against |V1-V0|^2 (it squares x0,y0 = P-V0 instead of x1,y1 = V1-V0).
    // Writing P-V0 = c*(V1-V0) with c > 0, the comparison dot > sqrLength
    // becomes c > c^2, which holds for 0 < c < 1, so 'order' is set to +2
    // (documented as P beyond the segment) when P is interior to the segment,
    // and to 0 (documented as interior) when c > 1, i.e. P is beyond V1. The
    // c = 1 case (P = V1, order +1) and the dot <= 0 cases are unaffected.
    // The sibling query toLineExtended does the analogous test correctly. No
    // upstream code calls this overload, so the port preserves the behavior
    // rather than silently changing it; use toLineExtended when the collinear
    // ordering must be correct.
    order: number;
}

export class PrimalQuery2 {
    private mNumVertices: number;
    private mVertices: readonly Vector[];

    // The caller is responsible for ensuring that the array is not empty
    // before calling queries and that the indices passed to the queries are
    // valid. The class does no range checking. The port of the upstream
    // default constructor is 'new PrimalQuery2()'.
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

    // For a line with origin V0 and direction <V0,V1>, toLine returns
    //   +1, P on right of line
    //   -1, P on left of line
    //    0, P on the line
    toLine(test: number | Vector, v0: number, v1: number): number {
        const p = this.point(test);
        const vec0 = this.mVertices[v0];
        const vec1 = this.mVertices[v1];

        const x0 = p.values[0] - vec0.values[0];
        const y0 = p.values[1] - vec0.values[1];
        const x1 = vec1.values[0] - vec0.values[0];
        const y1 = vec1.values[1] - vec0.values[1];
        const x0y1 = x0 * y1;
        const x1y0 = x1 * y0;
        const det = x0y1 - x1y0;
        const zero = 0;

        return (det > zero ? +1 : (det < zero ? -1 : 0));
    }

    // For a line with origin V0 and direction <V0,V1>, toLineWithOrder
    // returns the same sign as toLine together with the 'order' of P
    // described by PrimalQuery2ToLineOrderResult.
    toLineWithOrder(test: number | Vector, v0: number, v1: number):
        PrimalQuery2ToLineOrderResult {
        const p = this.point(test);
        const vec0 = this.mVertices[v0];
        const vec1 = this.mVertices[v1];

        const x0 = p.values[0] - vec0.values[0];
        const y0 = p.values[1] - vec0.values[1];
        const x1 = vec1.values[0] - vec0.values[0];
        const y1 = vec1.values[1] - vec0.values[1];
        const x0y1 = x0 * y1;
        const x1y0 = x1 * y0;
        const det = x0y1 - x1y0;
        const zero = 0;

        if (det > zero) {
            return { sign: +1, order: +3 };
        }

        if (det < zero) {
            return { sign: -1, order: -3 };
        }

        let order: number;
        const x0x1 = x0 * x1;
        const y0y1 = y0 * y1;
        const dot = x0x1 + y0y1;
        if (dot === zero) {
            order = -1;
        }
        else if (dot < zero) {
            order = -2;
        }
        else {
            // Upstream compares against |P-V0|^2 here rather than
            // |V1-V0|^2; see PrimalQuery2ToLineOrderResult. The quirk is
            // preserved.
            const x0x0 = x0 * x0;
            const y0y0 = y0 * y0;
            const sqrLength = x0x0 + y0y0;
            if (dot === sqrLength) {
                order = +1;
            }
            else if (dot > sqrLength) {
                order = +2;
            }
            else {
                order = 0;
            }
        }

        return { sign: 0, order };
    }

    // For a triangle with counterclockwise vertices V0, V1, and V2,
    // toTriangle returns
    //   +1, P outside triangle
    //   -1, P inside triangle
    //    0, P on triangle
    toTriangle(test: number | Vector, v0: number, v1: number, v2: number): number {
        const p = this.point(test);

        const sign0 = this.toLine(p, v1, v2);
        if (sign0 > 0) {
            return +1;
        }

        const sign1 = this.toLine(p, v0, v2);
        if (sign1 < 0) {
            return +1;
        }

        const sign2 = this.toLine(p, v0, v1);
        if (sign2 > 0) {
            return +1;
        }

        return ((sign0 !== 0 && sign1 !== 0 && sign2 !== 0) ? -1 : 0);
    }

    // For a triangle with counterclockwise vertices V0, V1, and V2,
    // toCircumcircle returns
    //   +1, P outside circumcircle of triangle
    //   -1, P inside circumcircle of triangle
    //    0, P on circumcircle of triangle
    toCircumcircle(test: number | Vector, v0: number, v1: number, v2: number): number {
        const p = this.point(test);
        const vec0 = this.mVertices[v0];
        const vec1 = this.mVertices[v1];
        const vec2 = this.mVertices[v2];

        const x0 = vec0.values[0] - p.values[0];
        const y0 = vec0.values[1] - p.values[1];
        const s00 = vec0.values[0] + p.values[0];
        const s01 = vec0.values[1] + p.values[1];
        const t00 = s00 * x0;
        const t01 = s01 * y0;
        const z0 = t00 + t01;

        const x1 = vec1.values[0] - p.values[0];
        const y1 = vec1.values[1] - p.values[1];
        const s10 = vec1.values[0] + p.values[0];
        const s11 = vec1.values[1] + p.values[1];
        const t10 = s10 * x1;
        const t11 = s11 * y1;
        const z1 = t10 + t11;

        const x2 = vec2.values[0] - p.values[0];
        const y2 = vec2.values[1] - p.values[1];
        const s20 = vec2.values[0] + p.values[0];
        const s21 = vec2.values[1] + p.values[1];
        const t20 = s20 * x2;
        const t21 = s21 * y2;
        const z2 = t20 + t21;

        const y0z1 = y0 * z1;
        const y0z2 = y0 * z2;
        const y1z0 = y1 * z0;
        const y1z2 = y1 * z2;
        const y2z0 = y2 * z0;
        const y2z1 = y2 * z1;
        const c0 = y1z2 - y2z1;
        const c1 = y2z0 - y0z2;
        const c2 = y0z1 - y1z0;
        const x0c0 = x0 * c0;
        const x1c1 = x1 * c1;
        const x2c2 = x2 * c2;
        const term = x0c0 + x1c1;
        const det = term + x2c2;
        const zero = 0;

        return (det < zero ? 1 : (det > zero ? -1 : 0));
    }

    // An extended classification of the relationship of the point P to the
    // line segment <Q0,Q1>. See PrimalQuery2OrderType for the meaning of the
    // returned value. Unlike the other queries, this one takes the three
    // points directly rather than indices into the vertex array.
    toLineExtended(P: Vector, Q0: Vector, Q1: Vector): PrimalQuery2OrderType {
        const zero = 0;

        const x0 = Q1.values[0] - Q0.values[0];
        const y0 = Q1.values[1] - Q0.values[1];
        if (x0 === zero && y0 === zero) {
            return PrimalQuery2OrderType.Q0_EQUALS_Q1;
        }

        const x1 = P.values[0] - Q0.values[0];
        const y1 = P.values[1] - Q0.values[1];
        if (x1 === zero && y1 === zero) {
            return PrimalQuery2OrderType.P_EQUALS_Q0;
        }

        const x2 = P.values[0] - Q1.values[0];
        const y2 = P.values[1] - Q1.values[1];
        if (x2 === zero && y2 === zero) {
            return PrimalQuery2OrderType.P_EQUALS_Q1;
        }

        // The theoretical classification relies on computing exactly the sign
        // of the determinant. Numerical roundoff errors can cause
        // misclassification.
        const x0y1 = x0 * y1;
        const x1y0 = x1 * y0;
        const det = x0y1 - x1y0;

        if (det !== zero) {
            if (det > zero) {
                // The points form a counterclockwise triangle <P,Q0,Q1>.
                return PrimalQuery2OrderType.POSITIVE;
            }
            else {
                // The points form a clockwise triangle <P,Q1,Q0>.
                return PrimalQuery2OrderType.NEGATIVE;
            }
        }
        else {
            // The points are collinear; P is on the line through Q0 and Q1.
            const x0x1 = x0 * x1;
            const y0y1 = y0 * y1;
            const dot = x0x1 + y0y1;
            if (dot < zero) {
                // The line ordering is <P,Q0,Q1>.
                return PrimalQuery2OrderType.COLLINEAR_LEFT;
            }

            const x0x0 = x0 * x0;
            const y0y0 = y0 * y0;
            const sqrLength = x0x0 + y0y0;
            if (dot > sqrLength) {
                // The line ordering is <Q0,Q1,P>.
                return PrimalQuery2OrderType.COLLINEAR_RIGHT;
            }

            // The line ordering is <Q0,P,Q1> with P strictly between Q0 and
            // Q1.
            return PrimalQuery2OrderType.COLLINEAR_CONTAIN;
        }
    }

    // The port of the upstream overload pairs: a query taking an index i uses
    // mVertices[i] as the test point.
    private point(test: number | Vector): Vector {
        return (typeof test === 'number' ? this.mVertices[test] : test);
    }
}
