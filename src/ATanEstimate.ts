// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ATanEstimate.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Minimax polynomial approximations to atan(x). The polynomial p(x) of
// degree D has only odd-power terms, is required to have linear term x,
// and p(1) = atan(1) = pi/4. It minimizes the quantity
// maximum{|atan(x) - p(x)| : x in [-1,1]} over all polynomials of
// degree D subject to the constraints mentioned.
//
// The upstream C++ selects the degree with a template parameter checked by
// static_assert. This port takes the degree as a runtime argument validated
// by logAssert; see atanEstimate below.

import { GTE_C_HALF_PI } from './Constants.js';
import { logAssert } from './Logger.js';

const C_ATAN_EST_COEFF: readonly (readonly number[])[] = [
    [   // degree 3
        +1.0,
        -2.1460183660255172e-1
    ],
    [   // degree 5
        +1.0,
        -3.0189478312144946e-1,
        +8.7292946518897740e-2
    ],
    [   // degree 7
        +1.0,
        -3.2570157599356531e-1,
        +1.5342994884206673e-1,
        -4.2330209451053591e-2
    ],
    [   // degree 9
        +1.0,
        -3.3157878236439586e-1,
        +1.8383034738018011e-1,
        -8.9253037587244677e-2,
        +2.2399635968909593e-2
    ],
    [   // degree 11
        +1.0,
        -3.3294527685374087e-1,
        +1.9498657165383548e-1,
        -1.1921576270475498e-1,
        +5.5063351366968050e-2,
        -1.2490720064867844e-2
    ],
    [   // degree 13
        +1.0,
        -3.3324998579202170e-1,
        +1.9856563505717162e-1,
        -1.3374657325451267e-1,
        +8.1675882859940430e-2,
        -3.5059680836411644e-2,
        +7.2128853633444123e-3
    ]
];

const C_ATAN_EST_MAX_ERROR: readonly number[] = [
    1.5970326392625e-2,  // degree 3
    1.3509832247375e-3,  // degree 5
    1.5051227215525e-4,  // degree 7
    1.8921598624725e-5,  // degree 9
    2.5477725020825e-6,  // degree 11
    3.5859106295450e-7   // degree 13
];

function validateDegree(degree: number): void {
    logAssert(degree === Math.trunc(degree) && (degree & 1) === 1
        && 1 <= (degree - 1) / 2 && (degree - 1) / 2 <= 6,
        'Invalid degree.');
}

// The input constraint is x in [-1,1]. The degree must be odd in {3..13}.
// For example a degree-3 estimate is
//   const x = ...;  // in [-1,1]
//   const result = atanEstimate(x, 3);
export function atanEstimate(x: number, degree: number): number {
    validateDegree(degree);

    const coeff = C_ATAN_EST_COEFF[(degree - 3) / 2];
    const last = (degree - 1) / 2;
    const xsqr = x * x;
    let poly = coeff[last];
    for (let i = 0, index = last - 1; i < last; ++i, --index) {
        poly = coeff[index] + poly * xsqr;
    }
    poly = poly * x;
    return poly;
}

// The input x can be any real number. Range reduction is used
// via the identities atan(x) = pi/2 - atan(1/x) for x > 0 and
// atan(x) = -pi/2 - atan(1/x) for x < 0. For example,
//   const x = ...;  // any real number
//   const result = atanEstimateRR(x, 3);
export function atanEstimateRR(x: number, degree: number): number {
    validateDegree(degree);

    if (Math.abs(x) <= 1) {
        return atanEstimate(x, degree);
    } else if (x > 1) {
        return GTE_C_HALF_PI - atanEstimate(1 / x, degree);
    } else {
        return -GTE_C_HALF_PI - atanEstimate(1 / x, degree);
    }
}

export function getATanEstimateMaxError(degree: number): number {
    validateDegree(degree);

    return C_ATAN_EST_MAX_ERROR[(degree - 3) / 2];
}
