// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SqrtEstimate.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Minimax polynomial approximations to sqrt(x). The polynomial p(x) of
// degree D minimizes the quantity maximum{|sqrt(x) - p(x)| : x in [1,2]}
// over all polynomials of degree D.
//
// The upstream C++ selects the degree with a template parameter checked by
// static_assert. This port takes the degree as a runtime argument validated
// by logAssert; see sqrtEstimate below.

import { GTE_C_SQRT_2 } from './Constants';
import { logAssert } from './Logger';

const C_SQRT_EST_COEFF: readonly (readonly number[])[] = [
    [   // degree 1
        +1.0,
        +4.1421356237309505e-1
    ],
    [   // degree 2
        +1.0,
        +4.8563183076125260e-1,
        -7.1418268388157458e-2
    ],
    [   // degree 3
        +1.0,
        +4.9750045320242231e-1,
        -1.0787308044477850e-1,
        +2.4586189615451115e-2
    ],
    [   // degree 4
        +1.0,
        +4.9955939832918816e-1,
        -1.2024066151943025e-1,
        +4.5461507257698486e-2,
        -1.0566681694362146e-2
    ],
    [   // degree 5
        +1.0,
        +4.9992197660031912e-1,
        -1.2378506719245053e-1,
        +5.6122776972699739e-2,
        -2.3128836281145482e-2,
        +5.0827122737047148e-3
    ],
    [   // degree 6
        +1.0,
        +4.9998616695784914e-1,
        -1.2470733323278438e-1,
        +6.0388587356982271e-2,
        -3.1692053551807930e-2,
        +1.2856590305148075e-2,
        -2.6183954624343642e-3
    ],
    [   // degree 7
        +1.0,
        +4.9999754817809228e-1,
        -1.2493243476353655e-1,
        +6.1859954146370910e-2,
        -3.6091595023208356e-2,
        +1.9483946523450868e-2,
        -7.5166134568007692e-3,
        +1.4127567687864939e-3
    ],
    [   // degree 8
        +1.0,
        +4.9999956583056759e-1,
        -1.2498490369914350e-1,
        +6.2318494667579216e-2,
        -3.7982961896432244e-2,
        +2.3642612312869460e-2,
        -1.2529377587270574e-2,
        +4.5382426960713929e-3,
        -7.8810995273670414e-4
    ]
];

const C_SQRT_EST_MAX_ERROR: readonly number[] = [
    1.7766952966369e-2,  // degree 1
    1.1795695163111e-3,  // degree 2
    1.1309620116485e-4,  // degree 3
    1.2741170151820e-5,  // degree 4
    1.5725569051384e-6,  // degree 5
    2.0584162152560e-7,  // degree 6
    2.8072338675856e-8,  // degree 7
    3.9468401880072e-9   // degree 8
];

// Equivalent of std::ldexp(value, power); see Exp2Estimate.ts for details
// of the two-step scaling.
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
//   const result = sqrtEstimate(x, 3);
export function sqrtEstimate(x: number, degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    const coeff = C_SQRT_EST_COEFF[degree - 1];
    const last = degree;
    const t = x - 1;  // t in [0,1]
    let poly = coeff[last];
    for (let i = 0, index = last - 1; i < last; ++i, --index) {
        poly = coeff[index] + poly * t;
    }
    return poly;
}

// The input constraint is x >= 0. Range reduction is used to generate a
// value y in [1,2], call sqrtEstimate(y) and then combine the output with
// the proper exponent to obtain the approximation. For example a degree-3
// estimate is
//   const x = ...;  // x >= 0
//   const result = sqrtEstimateRR(x, 3);
export function sqrtEstimateRR(x: number, degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    // Apply the reduction.
    let p = frexpExponent(x);
    const y = 2 * ldexp(x, -p);  // y in [1,2)
    --p;
    const adj = (1 & p) * GTE_C_SQRT_2 + (1 & ~p) * 1;
    p >>= 1;

    // Evaluate the polynomial on the reduced range.
    const poly = sqrtEstimate(y, degree);

    // Combine the outputs.
    return adj * ldexp(poly, p);
}

export function getSqrtEstimateMaxError(degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    return C_SQRT_EST_MAX_ERROR[degree - 1];
}
