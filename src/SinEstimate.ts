// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SinEstimate.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Minimax polynomial approximations to sin(x). The polynomial p(x) of
// degree D has only odd-power terms, is required to have linear term x,
// and p(pi/2) = sin(pi/2) = 1. It minimizes the quantity
// maximum{|sin(x) - p(x)| : x in [-pi/2,pi/2]} over all polynomials of
// degree D subject to the constraints mentioned.
//
// The upstream C++ selects the degree with a template parameter checked by
// static_assert. This port takes the degree as a runtime argument validated
// by logAssert; see sinEstimate below.

import { GTE_C_HALF_PI, GTE_C_PI, GTE_C_TWO_PI } from './Constants';
import { logAssert } from './Logger';

const C_SIN_EST_COEFF: readonly (readonly number[])[] = [
    [   // degree 3
        +1.0,
        -1.4727245910375519e-1
    ],
    [   // degree 5
        +1.0,
        -1.6600599923812209e-1,
        +7.5924178409012000e-3
    ],
    [   // degree 7
        +1.0,
        -1.6665578084732124e-1,
        +8.3109378830028557e-3,
        -1.8447486103462252e-4
    ],
    [   // degree 9
        +1.0,
        -1.6666656235308897e-1,
        +8.3329962509886002e-3,
        -1.9805100675274190e-4,
        +2.5967200279475300e-6
    ],
    [   // degree 11
        +1.0,
        -1.6666666601721269e-1,
        +8.3333303183525942e-3,
        -1.9840782426250314e-4,
        +2.7521557770526783e-6,
        -2.3828544692960918e-8
    ]
];

const C_SIN_EST_MAX_ERROR: readonly number[] = [
    1.3481903639146e-2,  // degree 3
    1.4001209384651e-4,  // degree 5
    1.0205878939740e-6,  // degree 7
    5.2010783457846e-9,  // degree 9
    1.9323431743601e-11  // degree 11
];

function validateDegree(degree: number): void {
    logAssert(degree === Math.trunc(degree) && (degree & 1) === 1
        && 1 <= (degree - 1) / 2 && (degree - 1) / 2 <= 5,
        'Invalid degree.');
}

// Port of std::remainder(x, y): x - n*y where n is the integer nearest to
// x/y, with ties rounded to the nearest even integer. JavaScript has no
// equivalent (% truncates), so it is computed explicitly.
function ieeeRemainder(x: number, y: number): number {
    const q = x / y;
    // Math.round rounds halfway cases toward +infinity; adjust ties to the
    // even neighbor as IEEE remainder requires.
    let n = Math.round(q);
    if (n - q === 0.5 && n % 2 !== 0) {
        n -= 1;
    }
    return x - n * y;
}

// The input constraint is x in [-pi/2,pi/2]. The degree must be odd in
// {3..11}. For example a degree-3 estimate is
//   const x = ...;  // in [-pi/2,pi/2]
//   const result = sinEstimate(x, 3);
export function sinEstimate(x: number, degree: number): number {
    validateDegree(degree);

    const coeff = C_SIN_EST_COEFF[(degree - 3) / 2];
    const last = (degree - 1) / 2;
    const xsqr = x * x;
    let poly = coeff[last];
    for (let i = 0, index = last - 1; i < last; ++i, --index) {
        poly = coeff[index] + poly * xsqr;
    }
    poly = poly * x;
    return poly;
}

// The input x can be any real number. Range reduction is used to
// generate a value y in [-pi/2,pi/2] for which sin(y) = sin(x).
// For example a degree-3 estimate is
//   const x = ...;  // any real number
//   const result = sinEstimateRR(x, 3);
export function sinEstimateRR(x: number, degree: number): number {
    validateDegree(degree);

    // Map x to r in [-pi,pi].
    const r = ieeeRemainder(x, GTE_C_TWO_PI);

    // Map r to y in [-pi/2,pi/2] with sin(y) = sin(x).
    if (r > GTE_C_HALF_PI) {
        // r is in (pi/2,pi], so y = pi - r is in (-pi/2,0]
        return sinEstimate(GTE_C_PI - r, degree);
    } else if (r < -GTE_C_HALF_PI) {
        // r is in [-pi,-pi/2), so y = -pi - r is in [0,pi/2)
        return sinEstimate(-GTE_C_PI - r, degree);
    } else {
        // r is in [-pi/2,pi/2], y = r
        return sinEstimate(r, degree);
    }
}

export function getSinEstimateMaxError(degree: number): number {
    validateDegree(degree);

    return C_SIN_EST_MAX_ERROR[(degree - 3) / 2];
}
