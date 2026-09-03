// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Projection.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The algorithm for the perspective projection of an ellipsoid onto a plane
// is described in
// https://www.geometrictools.com/Documentation/PerspectiveProjectionEllipsoid.pdf
//
// Port notes:
// - Upstream overloads 'Project' on the primitive dimension. The name
//   'project' is already taken by Vector.ts (the port of Vector's Project),
//   so the ports are named projectEllipse2 and projectEllipsoid3, following
//   the Cont* precedent of suffixing with the primitive type.
// - Output reference parameters become returned object literals:
//   'Project(..., smin, smax)' returns { smin, smax } and
//   'PerspectiveProject(..., ellipse)' returns { ellipse, isEllipse }, where
//   isEllipse is the value that Hyperellipsoid.fromCoefficientsABC returns.
//   Upstream discards that value; the port surfaces it because the
//   preconditions of the query (the ellipsoid strictly between the eyepoint
//   and the view plane) are the caller's responsibility and a violation makes
//   the returned ellipse meaningless.
// - The Ellipse2/Ellipsoid3/Line2/Line3/Plane3 aliases are the runtime-
//   dimension Hyperellipsoid/Line/Hyperplane; the dimensions are checked at
//   runtime.

import type { Ellipse2, Ellipsoid3 } from './Hyperellipsoid.js';
import { Hyperellipsoid } from './Hyperellipsoid.js';
import type { Plane3 } from './Hyperplane.js';
import type { Line2, Line3 } from './Line.js';
import { logAssert } from './Logger.js';
import { Matrix, mulMatrix, outerProduct, subMatrix } from './Matrix.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import { computeOrthogonalComplement3 } from './Vector3.js';

// Orthogonally project an ellipse onto a line. The projection interval is
// [smin,smax] and corresponds to the line segment P + s * D, where
// smin <= s <= smax.
export function projectEllipse2(ellipse: Ellipse2, line: Line2):
    { smin: number, smax: number } {
    logAssert(ellipse.dimension === 2 && line.dimension === 2,
        'Projection: expecting dimension 2.');

    // Center of projection interval.
    const center = dot(line.direction, sub(ellipse.center, line.origin));

    // Radius of projection interval.
    const tmp = [
        ellipse.extent.values[0] * dot(line.direction, ellipse.axis[0]),
        ellipse.extent.values[1] * dot(line.direction, ellipse.axis[1])
    ];
    const rSqr = tmp[0] * tmp[0] + tmp[1] * tmp[1];
    const radius = Math.sqrt(rSqr);

    return { smin: center - radius, smax: center + radius };
}

// Orthogonally project an ellipsoid onto a line. The projection interval is
// [smin,smax] and corresponds to the line segment P + s * D, where
// smin <= s <= smax.
export function projectEllipsoid3(ellipsoid: Ellipsoid3, line: Line3):
    { smin: number, smax: number } {
    logAssert(ellipsoid.dimension === 3 && line.dimension === 3,
        'Projection: expecting dimension 3.');

    // Center of projection interval.
    const center = dot(line.direction, sub(ellipsoid.center, line.origin));

    // Radius of projection interval.
    const tmp = [
        ellipsoid.extent.values[0] * dot(line.direction, ellipsoid.axis[0]),
        ellipsoid.extent.values[1] * dot(line.direction, ellipsoid.axis[1]),
        ellipsoid.extent.values[2] * dot(line.direction, ellipsoid.axis[2])
    ];
    const rSqr = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
    const radius = Math.sqrt(rSqr);

    return { smin: center - radius, smax: center + radius };
}

// Perspectively project an ellipsoid onto a plane.
//
// The ellipsoid has center C, axes A[i] and extents e[i] for 0 <= i <= 2.
//
// The eyepoint is E.
//
// The view plane is Dot(N,X) = d, where N is a unit-length normal vector.
// Choose U and V so that {U,V,N} is a right-handed orthonormal set; that is,
// the vectors are unit length, mutually perpendicular and N = Cross(U,V). N
// must be directed away from E in the sense that the point K on the plane
// closest to E is K = E + n * N with n > 0. When using a view frustum, n is
// the 'near' distance (from the eyepoint to the view plane). The plane
// equation is then
//   0 = Dot(N,X-K) = Dot(N,X) - Dot(N,E) - n = d - Dot(N,E) - n
// so that n = d - Dot(N,E).
//
// The ellipsoid must be between the eyepoint and the view plane in the sense
// that all rays from the eyepoint that intersect the ellipsoid must also
// intersect the view plane. The precondition test is to project the ellipsoid
// onto the line E + s * N to obtain interval [smin,smax] where smin > 0. The
// function projectEllipsoid3(ellipsoid, line) defined previously in this file
// can be used to verify the precondition. If the precondition is satisfied,
// the projection is an ellipse in the plane. If the precondition is not
// satisfied, the projection is a conic section that is not an ellipse or it
// is the empty set.
//
// The output is the equation of the ellipse in 2D. The projected ellipse
// coordinates Y = (y0,y1) are the view plane coordinates of the actual 3D
// ellipse points X = K + y0 * U + y1 * V = K + J * Y, where J is a 3x2 matrix
// whose columns are U and V.

// Use this query when you have a single plane and a single ellipsoid to
// project onto the plane.
export function perspectiveProject(ellipsoid: Ellipsoid3, E: Vector,
    plane: Plane3): { ellipse: Ellipse2, isEllipse: boolean };

// Use this query when you have a single plane and multiple ellipsoids to
// project onto the plane. The vectors U and V and the near value n are
// precomputed.
export function perspectiveProject(ellipsoid: Ellipsoid3, E: Vector,
    N: Vector, U: Vector, V: Vector, n: number):
    { ellipse: Ellipse2, isEllipse: boolean };

export function perspectiveProject(ellipsoid: Ellipsoid3, E: Vector,
    arg2: Plane3 | Vector, U?: Vector, V?: Vector, n?: number):
    { ellipse: Ellipse2, isEllipse: boolean } {
    if (U === undefined) {
        const plane = arg2 as Plane3;
        logAssert(plane.dimension === 3, 'Projection: expecting dimension 3.');
        const basis: Vector[] = [plane.normal.clone(), new Vector(3),
            new Vector(3)];
        computeOrthogonalComplement3(1, basis);
        const planeN = plane.normal;
        const planeU = basis[1];
        const planeV = basis[2];
        const near = plane.constant - dot(planeN, E);
        return perspectiveProject(ellipsoid, E, planeN, planeU, planeV, near);
    }

    const N = arg2 as Vector;
    logAssert(ellipsoid.dimension === 3 && E.size === 3 && N.size === 3
        && U.size === 3 && (V as Vector).size === 3,
        'Projection: expecting dimension 3.');
    const VV = V as Vector;
    const near = n as number;

    const two = 2;
    const four = 4;

    // Compute the coefficients for the ellipsoid represented by the quadratic
    // equation X^T*A*X + B^T*X + C = 0.
    const { A, B, C } = ellipsoid.toCoefficientsABC();

    // Compute the matrix M; see PerspectiveProjectionEllipsoid.pdf for the
    // mathematical details.
    const AE = mulMatrix(A, E) as Vector;
    const qformEAE = dot(E, AE);
    const dotBE = dot(B, E);
    const quadE = four * (qformEAE + dotBE + C);
    const Bp2AE = add(B, mul(AE, two));
    const M = subMatrix(outerProduct(Bp2AE, Bp2AE),
        mulMatrix(quadE, A) as Matrix);

    // Compute the coefficients for the projected ellipse.
    const MU = mulMatrix(M, U) as Vector;
    const MV = mulMatrix(M, VV) as Vector;
    const MN = mulMatrix(M, N) as Vector;
    const twoN = two * near;
    const AOut = new Matrix(2, 2);
    const BOut = new Vector(2);
    AOut.set(0, 0, dot(U, MU));
    AOut.set(0, 1, dot(U, MV));
    AOut.set(1, 0, AOut.get(0, 1));
    AOut.set(1, 1, dot(VV, MV));
    BOut.values[0] = twoN * dot(U, MN);
    BOut.values[1] = twoN * dot(VV, MN);
    const COut = near * near * dot(N, MN);

    // Extract the ellipse center, axis directions and extents.
    const ellipse = new Hyperellipsoid(2);
    const isEllipse = ellipse.fromCoefficientsABC(AOut, BOut, COut);
    return { ellipse, isEllipse };
}
