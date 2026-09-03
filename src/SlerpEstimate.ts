// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SlerpEstimate.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Read the comments in Slerp.ts about the slerp function. In particular, if
// you are using quaternions to represent rotations, read the comments about
// preprocessing the quaternions before calling slerp. The slerp functions in
// Slerp.ts require angles in [0,pi). The first two slerp estimates
// implemented in this file require angles in [0,pi/2], because the estimates
// are based on Chebyshev ratio estimates that have the same angle
// requirement. The third estimate that uses the qh inputs allows for angles
// in [0,pi).
//
// The upstream C++ selects the degree with a template parameter checked by
// static_assert. This port takes the degree as a runtime argument validated
// by logAssert.

import { chebyshevRatioEstimate } from './ChebyshevRatioEstimate.js';
import { logAssert } from './Logger.js';

function assertInputs(q0: readonly number[], q1: readonly number[], degree: number): void {
    logAssert(q0.length >= 2, 'Invalid dimension.');
    logAssert(q0.length === q1.length, 'Mismatched dimensions.');
    logAssert(1 <= degree && degree <= 16 && degree === Math.trunc(degree),
        'Invalid degree.');
}

// The angle between q0 and q1 is in [0,pi/2].
export function slerpEstimate(t: number, q0: readonly number[],
    q1: readonly number[], degree: number): number[]
{
    assertInputs(q0, q1, degree);

    const n = q0.length;
    let cosA = 0;
    for (let i = 0; i < n; ++i) {
        cosA += q0[i] * q1[i];
    }

    const f = chebyshevRatioEstimate(t, cosA, degree);
    const result = new Array<number>(n);
    for (let i = 0; i < n; ++i) {
        result[i] = f[0] * q0[i] + f[1] * q1[i];
    }
    return result;
}

// The angle between q0 and q1 must be in [0,pi/2] and cosA = Dot(q0,q1).
export function slerpEstimateUsingCosAngle(t: number, q0: readonly number[],
    q1: readonly number[], cosA: number, degree: number): number[]
{
    assertInputs(q0, q1, degree);

    const n = q0.length;
    const f = chebyshevRatioEstimate(t, cosA, degree);
    const result = new Array<number>(n);
    for (let i = 0; i < n; ++i) {
        result[i] = f[0] * q0[i] + f[1] * q1[i];
    }
    return result;
}

// The angle between q0 and q1 is in [0,pi). The input qh is halfway between
// q0 and q1 along a hyperspherical arc. If cosA = Dot(q0,q1), then
// cosAH = sqrt((1+cosA)/2) and qh = (q0+q1)/(2*cosAH).
export function slerpEstimateUsingMidpoint(t: number, q0: readonly number[],
    q1: readonly number[], qh: readonly number[], cosAH: number,
    degree: number): number[]
{
    assertInputs(q0, q1, degree);
    logAssert(qh.length === q0.length, 'Mismatched dimensions.');

    const n = q0.length;
    const result = new Array<number>(n);
    const twoT = 2 * t;
    if (twoT <= 1) {
        const f = chebyshevRatioEstimate(twoT, cosAH, degree);
        for (let i = 0; i < n; ++i) {
            result[i] = f[0] * q0[i] + f[1] * qh[i];
        }
    } else {
        const f = chebyshevRatioEstimate(twoT - 1, cosAH, degree);
        for (let i = 0; i < n; ++i) {
            result[i] = f[0] * qh[i] + f[1] * q1[i];
        }
    }
    return result;
}
