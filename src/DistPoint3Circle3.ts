// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPoint3Circle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The 3D point-circle distance algorithm is described in
// https://www.geometrictools.com/Documentation/DistanceToCircle3.pdf
// The notation used in the code matches that of the document.
//
// The input point is stored in the member closest[0]. If a single point on
// the circle is closest to the input point, the member closest[1] is set to
// that point and the equidistant member is set to false. If the entire circle
// is equidistant to the point, the member closest[1] is set to C+r*U, where C
// is the circle center, r is the circle radius and U is a vector
// perpendicular to the normal N for the plane of the circle. Moreover, the
// equidistant member is set to true.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector3<T>, Circle3<T>>' becomes the
// class DistPoint3Circle3 with the result type DistPoint3Circle3Result. The
// private member function ComputeOrthogonalBasis becomes the module-private
// computeOrthogonalBasis; only the numInputs == 1 path is reachable from this
// file, but the port keeps the other branch for fidelity.

import type { Circle3 } from './Circle3.js';
import type { DCPQuery } from './DCPQuery.js';
import { logAssert } from './Logger.js';
import { Vector, add, div, dot, getOrthogonal, length, mul, sub } from './Vector.js';
import { cross } from './Vector3.js';

export interface DistPoint3Circle3Result {
    distance: number;
    sqrDistance: number;

    // closest[0] is the input point, closest[1] is the closest circle point.
    closest: [Vector, Vector];

    // True when every circle point is equidistant from the input point, that
    // is, when the point is on the normal line through the circle center.
    equidistant: boolean;
}

function isZero3(v: Vector): boolean {
    return v.values[0] === 0 && v.values[1] === 0 && v.values[2] === 0;
}

// The upstream 'ComputeOrthogonalBasis(numInputs, v0, v1, v2)'. The vectors
// v1 and v2 are outputs (v1 is also an input when numInputs is 2 or 3). The
// boolean return value reports whether v2 is nonzero.
function computeOrthogonalBasis(numInputs: number, v0: Vector,
    v1In: Vector): { v1: Vector, v2: Vector, valid: boolean } {
    logAssert(1 <= numInputs && numInputs <= 3,
        'Invalid number of inputs.');

    let v1: Vector;
    if (numInputs === 1) {
        if (Math.abs(v0.values[0]) > Math.abs(v0.values[1])) {
            v1 = Vector.fromArray([-v0.values[2], 0, v0.values[0]]);
        }
        else {
            v1 = Vector.fromArray([0, v0.values[2], -v0.values[1]]);
        }
    }
    else {
        // numInputs is 2 or 3.
        v1 = sub(mul(dot(v0, v0), v1In), mul(dot(v1In, v0), v0));
    }

    if (isZero3(v1)) {
        return { v1, v2: new Vector(3), valid: false };
    }

    const v2 = cross(v0, v1);
    return { v1, v2, valid: !isZero3(v2) };
}

export class DistPoint3Circle3
    implements DCPQuery<Vector, Circle3, DistPoint3Circle3Result> {
    compute(point: Vector, circle: Circle3): DistPoint3Circle3Result {
        // The projection of P-C onto the plane of the circle is
        // Q - C = (P - C) - Dot(N, P - C) * N. When P is nearly on the normal
        // line C + t * N, Q - C is nearly the zero vector. In this case,
        // floating-point rounding errors are a problem when the closest point
        // is computed as C + r * (Q - C) / Length(Q - C). The rounding errors
        // in Q - C are magnified by the division by length, leading to an
        // inaccurate result. Experiments indicate it is better to compute an
        // orthonormal basis {U, V, N}, where the vectors are unit length and
        // mutually perpendicular. The point is P = C + x * U + y * V + z * N,
        // with x = Dot(U, P - C), y = Dot(V, Q - C) and z = Dot(N, Q - C).
        // The projection is Q = C + x * U + y * V. The computation of U and V
        // involves normalizations (divisions by square roots) which can be
        // avoided by instead computing an orthogonal basis {U, V, N}, where
        // the vectors are mutually perpendicular but not required to be unit
        // length. U is computed by swapping two components of N with at least
        // one component nonzero and then negating a component. V is computed
        // as Cross(N, U). For example, if N = (n0, n1, n2) with n0 != 0 or
        // n1 != 0, then U = (-n1, n0, 0) and V = (-n0*n2, -n1*n2, n0^2+n1^2).
        // Observe that the length of V is |V| = |N|*|U|. In this case the
        // projection is
        //   Q - C = x * U + y * V,
        //   x = Dot(U, Q - C) / Dot(U, U)
        //   y = Dot(V, Q - C) / (Dot(U, U) * Dot(N, N))
        // It is sufficient to process the scaled
        //   Dot(N, N) * Dot(U, U) * (Q - C)
        // to avoid the divisions before normalization.
        const PmC = sub(point, circle.center);
        const N = circle.normal;
        const { v1: U, v2: V } = computeOrthogonalBasis(1, N, new Vector(3));
        const scaledQmC = add(mul(dot(N, N) * dot(U, PmC), U),
            mul(dot(V, PmC), V));
        const lengthScaledQmC = length(scaledQmC);
        if (lengthScaledQmC > 0) {
            const height = dot(N, PmC);
            const radial = length(cross(N, PmC)) - circle.radius;
            const sqrDistance = height * height + radial * radial;
            return {
                distance: Math.sqrt(sqrDistance),
                sqrDistance,
                closest: [
                    point.clone(),
                    add(circle.center,
                        mul(circle.radius, div(scaledQmC, lengthScaledQmC)))
                ],
                equidistant: false
            };
        }

        // All circle points are equidistant from P. Return one of them.
        const sqrDistance = dot(PmC, PmC) + circle.radius * circle.radius;
        return {
            distance: Math.sqrt(sqrDistance),
            sqrDistance,
            closest: [
                point.clone(),
                add(circle.center,
                    mul(circle.radius, getOrthogonal(N, true)))
            ],
            equidistant: true
        };
    }
}
