// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine3Torus3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The line is parameterized by L(t) = P + t * D, where P is a point on the
// line and D is a nonzero direction vector that is not necessarily unit
// length.
//
// The standard torus has center (0,0,0), plane of symmetry z = 0, axis of
// symmetry containing (0,0,0) in the direction (0,0,1), outer radius r0 and
// inner radius r1 with r0 > r1 (a "ring torus"). It is defined implicitly by
//   (x^2 + y^2 + z^2 + r0^2 - r1^2)^2 - 4 * r0^2 * (x^2 + y^2) = 0
// where (x,y,z) is a point on the torus. A parameterization is
//   x(u,v) = (r0 + r1 * cos(v)) * cos(u)
//   y(u,v) = (r0 + r1 * cos(v)) * sin(u)
//   z(u,v) = r1 * sin(v)
// for u in [0,2*pi) and v in [0,2*pi).
//
// Generally, the torus has center C with plane of symmetry containing C and
// having unit-length normal N. The axis of symmetry is the normal line to the
// plane at C. If X is a point on the torus, the implicit formulation is
//   (|X-C|^2 + r0^2 - r1^2)^2 - 4 * r0^2 * (|X-C|^2 - (Dot(N,X-C))^2) = 0
// Let D0 and D1 be unit-length vectors that span the symmetry plane where
// {D0,D1,N} is a right-handed orthonormal basis. A parameterization for the
// torus is
//   X(u,v) = C + (r0 + r1*cos(v))*(cos(u)*D0 + sin(u)*D1) + r1*sin(v)*N
// for u in [0,2*pi) and v in [0,2*pi).
//
// Compute the intersections of a line with a torus. The number of
// intersections is between 0 and 4. As noted, the line direction D does not
// have to be unit length. The normal vector N must be unit length, but notice
// that the implicit formulation has a term
//   (Dot(N,X-C))^2 = (X-C)^T * (N * N^T) * (X - C)
// If the normal were chosen to be nonzero but not unit length, say M, then
// N = M/|M|. The term can be modified to
//   (Dot(N,X-C))^2 = (X-C)^T * ((M * M^T)/|M|^2) * (X - C)
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream provides
// only an FIQuery specialization for this pair of primitives, which becomes
// IntrLine3Torus3FI. The upstream std::map<T,int32_t> of roots-with-
// multiplicities becomes the RootMultiplicity[] returned by
// RootsPolynomial.solveQuartic, which is already sorted ascending by root, so
// the ordering of the reported intersections matches upstream.

import type { Line } from './Line.js';
import type { Torus3 } from './Torus3.js';
import type { FIQuery } from './FIQuery.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import { Polynomial1 } from './Polynomial1.js';
import { RootsPolynomial } from './RootsPolynomial.js';

// The result of IntrLine3Torus3FI queries. The number of intersections is
// between 0 and 4. Only the first 'numIntersections' entries of
// lineParameter, torusParameter and point are valid.
export interface IntrLine3Torus3FIResult {
    intersect: boolean;
    numIntersections: number;

    // The line parameters t of the intersection points.
    lineParameter: [number, number, number, number];

    // The torus surface parameters (u,v) of the intersection points.
    torusParameter: [
        [number, number], [number, number],
        [number, number], [number, number]
    ];

    // The intersection points.
    point: [Vector, Vector, Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrLine3Torus3FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        lineParameter: [0, 0, 0, 0],
        torusParameter: [[0, 0], [0, 0], [0, 0], [0, 0]],
        point: [Vector.zero(3), Vector.zero(3), Vector.zero(3), Vector.zero(3)]
    };
}

// Read a polynomial coefficient, returning zero beyond the degree. Upstream
// indexes quartic[0] through quartic[4] directly; Polynomial1 drops leading
// zero coefficients, so a degenerate input that lowers the degree would be an
// out-of-range access in C++ and 'undefined' here. The helper keeps the
// mathematically correct value of zero for the missing high-order terms.
function coefficient(p: Polynomial1, i: number): number {
    return i <= p.getDegree() ? p.get(i) : 0;
}

// Find-intersection query for a line and a torus in 3D.
export class IntrLine3Torus3FI implements
    FIQuery<Line, Torus3, IntrLine3Torus3FIResult> {

    find(line: Line, torus: Torus3): IntrLine3Torus3FIResult {
        const result = defaultFIResult();

        // Short names for readability.
        const P = line.origin;
        const D = line.direction;
        const C = torus.center;
        const N = torus.normal;
        const r0 = torus.radius0;  // outer radius
        const r1 = torus.radius1;  // inner radius

        // Common intermediate terms.
        const r0Sqr = r0 * r0;
        const r1Sqr = r1 * r1;
        const PmC = sub(P, C);
        const sqrLenPmC = dot(PmC, PmC);
        const dotDPmC = dot(D, PmC);
        const sqrLenD = dot(D, D);
        const sqrLenN = dot(N, N);
        const dotND = dot(N, D);
        const dotNPmC = dot(N, PmC);

        // |X-C|^2
        const quad0 = new Polynomial1(2);
        quad0.set(0, sqrLenPmC);
        quad0.set(1, 2 * dotDPmC);
        quad0.set(2, sqrLenD);

        // |X-C|^2 + r0^2 - r1^2
        const quad1 = quad0.clone();
        quad1.set(0, quad1.get(0) + r0Sqr - r1Sqr);

        // Dot(N,X-C)
        const linear = Polynomial1.fromCoefficients([dotNPmC, dotND]);

        // Dot(N,X-C)^2 with adjustment for non-unit N
        const quad2 = linear.mul(linear).div(sqrLenN);

        // |X-C|^2 - (Dot(N,X-C))^2
        const quad3 = quad0.sub(quad2);

        // (|X-C|^2 + r0^2-r1^2)^2 - 4*r0^2 * (|X-C|^2 - (Dot(N,X-C))^2)
        const quartic = quad1.mul(quad1).sub(quad3.mul(4 * r0Sqr));

        // Solve the quartic.
        const rmMap = RootsPolynomial.solveQuartic(
            coefficient(quartic, 0), coefficient(quartic, 1),
            coefficient(quartic, 2), coefficient(quartic, 3),
            coefficient(quartic, 4));

        // Get the intersection parameters and points.
        result.numIntersections = rmMap.length;
        result.intersect = (result.numIntersections > 0);
        for (let i = 0; i < rmMap.length; ++i) {
            const t = rmMap[i].root;
            result.lineParameter[i] = t;
            result.point[i] = add(line.origin, mul(t, line.direction));
            const { u, v } = torus.getParameters(result.point[i]);
            result.torusParameter[i] = [u, v];
        }

        return result;
    }
}
