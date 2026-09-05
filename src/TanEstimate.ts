// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) TanEstimate.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Minimax polynomial approximations to tan(x). The polynomial p(x) of
// degree D has only odd-power terms, is required to have linear term x,
// and p(pi/4) = tan(pi/4) = 1. It minimizes the quantity
// maximum{|tan(x) - p(x)| : x in [-pi/4,pi/4]} over all polynomials of
// degree D subject to the constraints mentioned.
//
// The upstream C++ selects the degree with a template parameter checked by
// static_assert. This port takes the degree as a runtime argument validated
// by logAssert; see tanEstimate below.

import { GTE_C_HALF_PI, GTE_C_PI, GTE_C_QUARTER_PI } from './Constants.js';
import { logAssert } from './Logger.js';

const C_TAN_EST_COEFF: readonly (readonly number[])[] = [
    [   // degree 3
        1.0,
        4.4295926544736286e-1
    ],
    [   // degree 5
        1.0,
        3.1401320403542421e-1,
        2.0903948109240345e-1
    ],
    [   // degree 7
        1.0,
        3.3607213284422555e-1,
        1.1261037305184907e-1,
        9.8352099470524479e-2
    ],
    [   // degree 9
        1.0,
        3.3299232843941784e-1,
        1.3747843432474838e-1,
        3.7696344813028304e-2,
        4.6097377279281204e-2
    ],
    [   // degree 11
        1.0,
        3.3337224456224224e-1,
        1.3264516053824593e-1,
        5.8145237645931047e-2,
        1.0732193237572574e-2,
        2.1558456793513869e-2
    ],
    [   // degree 13
        1.0,
        3.3332916426394554e-1,
        1.3343404625112498e-1,
        5.3104565343119248e-2,
        2.5355038312682154e-2,
        1.8253255966556026e-3,
        1.0069407176615641e-2
    ]
];

const C_TAN_EST_MAX_ERROR: readonly number[] = [
    1.1661892256205e-2,  // degree 3
    5.8431854390146e-4,  // degree 5
    3.5418688397793e-5,  // degree 7
    2.2988173248307e-6,  // degree 9
    1.5426258070939e-7,  // degree 11
    1.0550265105991e-8   // degree 13
];

function validateDegree(degree: number): void {
    logAssert(degree === Math.trunc(degree) && (degree & 1) === 1
        && 1 <= (degree - 1) / 2 && (degree - 1) / 2 <= 6,
        'Invalid degree.');
}

// Port of std::remainder(x, y): x - n*y where n is the integer nearest to
// x/y, with ties rounded to the nearest even integer. JavaScript has no
// equivalent (% truncates), and the obvious x - Math.round(x/y)*y is not a
// substitute: std::remainder is computed as if in infinite precision, while
// the quotient x/y and the product n*y each round. That approximation is
// already wrong by 1e-6 at |x| = 1e10 and returns values outside
// [-|y|/2, |y|/2] beyond |x| = 1e14 (for example remainder(-9.77e15, 2*pi)
// evaluates to 4 instead of -3.141005...), which drives the range-reduced
// estimates below far outside the domain their polynomials approximate.
//
// The reduction here is exact: the loop is binary long division of |x| by
// |y| in which every d is |y| scaled by a power of two (exact) and every
// subtraction satisfies d <= r < 2*d, hence is exact by Sterbenz's lemma.
// The remaining r is |x| mod |y| in [0,|y|), and quotientOdd is the low bit
// of the quotient, which decides the halfway case.
function ieeeRemainder(x: number, y: number): number {
    if (!Number.isFinite(x) || Number.isNaN(y) || y === 0) {
        // std::remainder(+-inf, y), remainder(NaN, y) and remainder(x, 0)
        // are all NaN.
        return Number.NaN;
    }
    if (!Number.isFinite(y)) {
        // std::remainder(x, +-inf) is x for finite x.
        return x;
    }

    const ay = Math.abs(y);
    let r = Math.abs(x);
    let quotientOdd = false;
    if (r >= ay) {
        let d = ay;
        while (r >= 2 * d) {
            d *= 2;
        }
        for (; d >= ay; d *= 0.5) {
            if (r >= d) {
                r -= d;
                quotientOdd = (d === ay);
            }
        }
    }

    // Round the quotient to nearest with ties to even: subtract one more |y|
    // when the fractional part exceeds 1/2, or equals 1/2 with an odd
    // quotient. Doubling r is exact because r < |y|; in the one case where
    // 2*r would overflow, |y| is normal and halving it is exact instead.
    if (r <= Number.MAX_VALUE * 0.5) {
        const twiceR = 2 * r;
        if (twiceR > ay || (twiceR === ay && quotientOdd)) {
            r -= ay;
        }
    } else {
        const half = ay * 0.5;
        if (r > half || (r === half && quotientOdd)) {
            r -= ay;
        }
    }

    // std::remainder(-x, y) = -remainder(x, y), including the sign of zero.
    return (x < 0 || Object.is(x, -0)) ? -r : r;
}

// The input constraint is x in [-pi/4,pi/4]. The degree must be odd in
// {3..13}. For example a degree-3 estimate is
//   const x = ...;  // in [-pi/4,pi/4]
//   const result = tanEstimate(x, 3);
export function tanEstimate(x: number, degree: number): number {
    validateDegree(degree);

    const coeff = C_TAN_EST_COEFF[(degree - 3) / 2];
    const last = (degree - 1) / 2;
    const xsqr = x * x;
    let poly = coeff[last];
    for (let i = 0, index = last - 1; i < last; ++i, --index) {
        poly = coeff[index] + poly * xsqr;
    }
    poly = poly * x;
    return poly;
}

// The input x can be any real number. Range reduction is used to generate
// a value y in [-pi/2,pi/2]. If |y| <= pi/4, then the polynomial is
// evaluated. If y in (pi/4,pi/2), set z = y - pi/4 and use the identity
//   tan(y) = tan(z + pi/4) = [1 + tan(z)]/[1 - tan(z)]
// If y in (-pi/2,-pi/4), set z = y + pi/4 and use the identity
//   tan(y) = tan(z - pi/4) = -[1 - tan(z)]/[1 + tan(z)]
// Be careful when evaluating at y nearly pi/2, because tan(y) becomes
// infinite. For example a degree-3 estimate is
//   const x = ...;  // any real number
//   const result = tanEstimateRR(x, 3);
export function tanEstimateRR(x: number, degree: number): number {
    validateDegree(degree);

    // Map x to r in [-pi/2,pi/2].
    const r = ieeeRemainder(x, GTE_C_PI);

    // Map r to y in [-pi/2,pi/2] with tan(y) = tan(r).
    let y: number;
    if (r > GTE_C_HALF_PI) {
        y = r - GTE_C_PI;
    } else if (r < -GTE_C_HALF_PI) {
        y = r + GTE_C_PI;
    } else {
        y = r;
    }

    if (Math.abs(y) <= GTE_C_QUARTER_PI) {
        return tanEstimate(y, degree);
    } else if (y > GTE_C_QUARTER_PI) {
        const poly = tanEstimate(y - GTE_C_QUARTER_PI, degree);
        return (1 + poly) / (1 - poly);
    } else {
        const poly = tanEstimate(y + GTE_C_QUARTER_PI, degree);
        return (-1 + poly) / (1 + poly);
    }
}

export function getTanEstimateMaxError(degree: number): number {
    validateDegree(degree);

    return C_TAN_EST_MAX_ERROR[(degree - 3) / 2];
}
