// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Slerp.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The spherical linear interpolation (slerp) of unit-length vectors q0 and
// q1 for t in [0,1] and theta in (0,pi) is
//   slerp(t,q0,q1) = [sin((1-t)*theta)*q0 + sin(t*theta)*q1]/sin(theta)
// where theta is the angle between q0 and q1 [cos(theta) = Dot(q0,q1)]. This
// function is a parameterization of the great spherical arc between q0 and
// q1 on the unit hypersphere. Moreover, the parameterization has the
// property that a particle traveling along the arc does so with constant
// speed, where t is time.
//
// When applying slerp to unit-length quaternions (N = 4) that represent
// rotations, q and -q represent the same rotation. It is typical that a
// quaternion sequence is preprocessed by
//   const q: number[][] = ...;  // assuming initialized
//   for (let i0 = 0, i1 = 1; i1 < q.length; i0 = i1++) {
//       const cosA = dot(q[i0], q[i1]);
//       if (cosA < 0) {
//           q[i1] = negate(q[i1]);  // now dot(q[i0], q[i1]) >= 0
//       }
//   }
// so that the angle between consecutive quaternions is in [0,pi/2].
//
// The cosines might also be precomputed and passed to slerpUsingCosAngle.
//
// For numerical robustness of slerp, the quaternions can be preprocessed so
// that a quaternion is inserted between each pair of original quaternions.
// Given q0 and q1, the midpoint of the arc connecting them is qh so that
// A = Dot(q0,q1) and A/2 = Dot(q0,qh) = Dot(qh,q1). The midpoint is
// qh = slerp(1/2,q0,q1) = (q0 + q1)/|q0 + q1|. The preprocessing computes
// cosAH = sqrt((1 + cosA)/2) and qh = (q0 + q1)/(2*cosAH) for each pair, to
// be passed to slerpUsingMidpoint.

import { chebyshevRatiosUsingCosAngle } from './ChebyshevRatio';
import { logAssert } from './Logger';

function assertDimensions(q0: readonly number[], q1: readonly number[]): void {
    logAssert(q0.length >= 2, 'Invalid dimension.');
    logAssert(q0.length === q1.length, 'Mismatched dimensions.');
}

// The angle between q0 and q1 is in [0,pi).
export function slerp(t: number, q0: readonly number[], q1: readonly number[]): number[] {
    assertDimensions(q0, q1);

    const n = q0.length;
    let cosA = 0;
    for (let i = 0; i < n; ++i) {
        cosA += q0[i] * q1[i];
    }

    const f = chebyshevRatiosUsingCosAngle(t, cosA);
    const result = new Array<number>(n);
    for (let i = 0; i < n; ++i) {
        result[i] = f[0] * q0[i] + f[1] * q1[i];
    }
    return result;
}

// The angle between q0 and q1 must be in [0,pi) and cosA = Dot(q0,q1).
export function slerpUsingCosAngle(t: number, q0: readonly number[],
    q1: readonly number[], cosA: number): number[]
{
    assertDimensions(q0, q1);

    const n = q0.length;
    const f = chebyshevRatiosUsingCosAngle(t, cosA);
    const result = new Array<number>(n);
    for (let i = 0; i < n; ++i) {
        result[i] = f[0] * q0[i] + f[1] * q1[i];
    }
    return result;
}

// The angle between q0 and q1 is in [0,pi). The input qh is halfway between
// q0 and q1 along a hyperspherical arc. If cosA = Dot(q0,q1), then
// cosAH = sqrt((1+cosA)/2) and qh = (q0+q1)/(2*cosAH).
export function slerpUsingMidpoint(t: number, q0: readonly number[],
    q1: readonly number[], qh: readonly number[], cosAH: number): number[]
{
    assertDimensions(q0, q1);
    logAssert(qh.length === q0.length, 'Mismatched dimensions.');

    const n = q0.length;
    const result = new Array<number>(n);
    const twoT = 2 * t;
    if (twoT <= 1) {
        const f = chebyshevRatiosUsingCosAngle(twoT, cosAH);
        for (let i = 0; i < n; ++i) {
            result[i] = f[0] * q0[i] + f[1] * qh[i];
        }
    } else {
        const f = chebyshevRatiosUsingCosAngle(twoT - 1, cosAH);
        for (let i = 0; i < n; ++i) {
            result[i] = f[0] * qh[i] + f[1] * q1[i];
        }
    }
    return result;
}
