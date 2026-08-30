// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) CosEstimate.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Minimax polynomial approximations to cos(x). The polynomial p(x) of
// degree D has only even-power terms, is required to have constant term 1,
// and p(pi/2) = cos(pi/2) = 0. It minimizes the quantity
// maximum{|cos(x) - p(x)| : x in [-pi/2,pi/2]} over all polynomials of
// degree D subject to the constraints mentioned.
//
// The upstream C++ selects the degree with a template parameter checked by
// static_assert. This port takes the degree as a runtime argument validated
// by logAssert; see cosEstimate below.

import { GTE_C_HALF_PI, GTE_C_PI, GTE_C_TWO_PI } from './Constants';
import { logAssert } from './Logger';

const C_COS_EST_COEFF: readonly (readonly number[])[] = [
    [   // degree 2
        +1.0,
        -4.0528473456935105e-1
    ],
    [   // degree 4
        +1.0,
        -4.9607181958647262e-1,
        +3.6794619653489236e-2
    ],
    [   // degree 6
        +1.0,
        -4.9992746217057404e-1,
        +4.1493920348353308e-2,
        -1.2712435011987822e-3
    ],
    [   // degree 8
        +1.0,
        -4.9999925121358291e-1,
        +4.1663780117805693e-2,
        -1.3854239405310942e-3,
        +2.3154171575501259e-5
    ],
    [   // degree 10
        +1.0,
        -4.9999999508695869e-1,
        +4.1666638865338612e-2,
        -1.3888377661039897e-3,
        +2.4760495088926859e-5,
        -2.6051615464872668e-7
    ]
];

const C_COS_EST_MAX_ERROR: readonly number[] = [
    5.6009595954128e-2,  // degree 2
    9.1879932449727e-4,  // degree 4
    9.2028470144446e-6,  // degree 6
    5.9804535233743e-8,  // degree 8
    2.7008567604626e-10  // degree 10
];

function validateDegree(degree: number): void {
    logAssert(degree === Math.trunc(degree) && (degree & 1) === 0
        && 1 <= degree / 2 && degree / 2 <= 5,
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

// The input constraint is x in [-pi/2,pi/2]. The degree must be even in
// {2..10}. For example a degree-4 estimate is
//   const x = ...;  // in [-pi/2,pi/2]
//   const result = cosEstimate(x, 4);
export function cosEstimate(x: number, degree: number): number {
    validateDegree(degree);

    const coeff = C_COS_EST_COEFF[(degree - 2) / 2];
    const last = degree / 2;
    const xsqr = x * x;
    let poly = coeff[last];
    for (let i = 0, index = last - 1; i < last; ++i, --index) {
        poly = coeff[index] + poly * xsqr;
    }
    return poly;
}

// The input x can be any real number. Range reduction is used to generate
// a value y in [-pi/2,pi/2] and a sign s for which cos(y) = s * cos(x).
// For example a degree-4 estimate is
//   const x = ...;  // any real number
//   const result = cosEstimateRR(x, 4);
export function cosEstimateRR(x: number, degree: number): number {
    validateDegree(degree);

    // Map x to r in [-pi,pi].
    const r = ieeeRemainder(x, GTE_C_TWO_PI);

    // Map r to y in [-pi/2,pi/2] with cos(y) = sign * cos(x).
    if (r > GTE_C_HALF_PI) {
        // r is in (pi/2,pi], so y = pi - r is in (-pi/2,0], sign = -1
        return -cosEstimate(GTE_C_PI - r, degree);
    } else if (r < -GTE_C_HALF_PI) {
        // r is in [-pi,-pi/2), so y = -pi - r is in [0,pi/2), sign = -1
        return -cosEstimate(-GTE_C_PI - r, degree);
    } else {
        // r is in [-pi/2,pi/2], y = r, sign = +1
        return cosEstimate(r, degree);
    }
}

export function getCosEstimateMaxError(degree: number): number {
    validateDegree(degree);

    return C_COS_EST_MAX_ERROR[(degree - 2) / 2];
}
