// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) InvSqrtEstimate.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Minimax polynomial approximations to 1/sqrt(x). The polynomial p(x) of
// degree D minimizes the quantity maximum{|1/sqrt(x) - p(x)| : x in [1,2]}
// over all polynomials of degree D.
//
// The upstream C++ selects the degree with a template parameter checked by
// static_assert. This port takes the degree as a runtime argument validated
// by logAssert; see invSqrtEstimate below.

import { GTE_C_INV_SQRT_2 } from './Constants';
import { logAssert } from './Logger';

const C_INVSQRT_EST_COEFF: readonly (readonly number[])[] = [
    [   // degree 1
        +1.0,
        -2.9289321881345254e-1
    ],
    [   // degree 2
        +1.0,
        -4.4539812104566801e-1,
        +1.5250490223221547e-1
    ],
    [   // degree 3
        +1.0,
        -4.8703230993068791e-1,
        +2.8163710486669835e-1,
        -8.7498013749463421e-2
    ],
    [   // degree 4
        +1.0,
        -4.9710061558048779e-1,
        +3.4266247597676802e-1,
        -1.9106356536293490e-1,
        +5.2608486153198797e-2
    ],
    [   // degree 5
        +1.0,
        -4.9937760586004143e-1,
        +3.6508741295133973e-1,
        -2.5884890281853501e-1,
        +1.3275782221320753e-1,
        -3.2511945299404488e-2
    ],
    [   // degree 6
        +1.0,
        -4.9987029229547453e-1,
        +3.7220923604495226e-1,
        -2.9193067713256937e-1,
        +1.9937605991094642e-1,
        -9.3135712130901993e-2,
        +2.0458166789566690e-2
    ],
    [   // degree 7
        +1.0,
        -4.9997357250704977e-1,
        +3.7426216884998809e-1,
        -3.0539882498248971e-1,
        +2.3976005607005391e-1,
        -1.5410326351684489e-1,
        +6.5598809723041995e-2,
        -1.3038592450470787e-2
    ],
    [   // degree 8
        +1.0,
        -4.9999471066120371e-1,
        +3.7481415745794067e-1,
        -3.1023804387422160e-1,
        +2.5977002682930106e-1,
        -1.9818790717727097e-1,
        +1.1882414252613671e-1,
        -4.6270038088550791e-2,
        +8.3891541755747312e-3
    ]
];

const C_INVSQRT_EST_MAX_ERROR: readonly number[] = [
    3.7814314552702e-2,  // degree 1
    4.1953446330581e-3,  // degree 2
    5.6307702007275e-4,  // degree 3
    8.1513919990229e-5,  // degree 4
    1.2289367490981e-5,  // degree 5
    1.9001451476708e-6,  // degree 6
    2.9887737629242e-7,  // degree 7
    4.7597402907940e-8   // degree 8
];

// Equivalent of std::ldexp(value, power): value * 2^power with correct
// handling of overflow to infinity and underflow toward the subnormals.
// See Exp2Estimate.ts for details of the two-step scaling.
function ldexp(value: number, power: number): number {
    if (power >= -1022) {
        if (power <= 1023) {
            return value * Math.pow(2, power);
        }
        return value * Math.pow(2, 1023) * Math.pow(2, Math.min(power - 1023, 1023));
    }
    return value * Math.pow(2, -1022) * Math.pow(2, Math.max(power + 1022, -1074));
}

// The exponent of the std::frexp decomposition x = f * 2^p with f in
// [1/2,1), for finite x > 0. The IEEE exponent field is read directly so
// the result is exact, including for subnormal inputs.
function frexpExponent(x: number): number {
    if (x === 0) {
        // std::frexp(0) returns 0 with exponent 0.
        return 0;
    }
    const dv = new DataView(new ArrayBuffer(8));
    dv.setFloat64(0, x);
    let biasedExponent = (dv.getUint32(0) >>> 20) & 0x7ff;
    if (biasedExponent === 0) {
        // Subnormal input: normalize by scaling with 2^64 (exact) and
        // compensate in the exponent.
        dv.setFloat64(0, x * 18446744073709551616.0);
        biasedExponent = ((dv.getUint32(0) >>> 20) & 0x7ff) - 64;
    }
    return biasedExponent - 1022;
}

// The input constraint is x in [1,2]. For example a degree-3 estimate is
//   const x = ...;  // in [1,2]
//   const result = invSqrtEstimate(x, 3);
export function invSqrtEstimate(x: number, degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    const coeff = C_INVSQRT_EST_COEFF[degree - 1];
    const last = degree;
    const t = x - 1;  // t in [0,1]
    let poly = coeff[last];
    for (let i = 0, index = last - 1; i < last; ++i, --index) {
        poly = coeff[index] + poly * t;
    }
    return poly;
}

// The input constraint is x > 0. Range reduction is used to generate
// a value y in [1,2], call invSqrtEstimate(y) and then combine the output
// with the proper exponent to obtain the approximation. For example a
// degree-3 estimate is
//   const x = ...;  // x > 0
//   const result = invSqrtEstimateRR(x, 3);
export function invSqrtEstimateRR(x: number, degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    // Apply the reduction.
    let p = frexpExponent(x);
    const y = 2 * ldexp(x, -p);  // y in [1,2)
    --p;
    const adj = (1 & p) * GTE_C_INV_SQRT_2 + (1 & ~p) * 1;
    p = -(p >> 1);

    // Evaluate the polynomial on the reduced range.
    const poly = invSqrtEstimate(y, degree);

    // Combine the outputs.
    return adj * ldexp(poly, p);
}

export function getInvSqrtEstimateMaxError(degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    return C_INVSQRT_EST_MAX_ERROR[degree - 1];
}
