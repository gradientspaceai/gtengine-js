// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Log2Estimate.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Minimax polynomial approximations to log2(x). The polynomial p(x) of
// degree D minimizes the quantity maximum{|log2(x) - p(x)| : x in [1,2]}
// over all polynomials of degree D.
//
// The upstream C++ selects the degree with a template parameter checked by
// static_assert. This port takes the degree as a runtime argument validated
// by logAssert; see log2Estimate below.

import { logAssert } from './Logger';

const C_LOG2_EST_COEFF: readonly (readonly number[])[] = [
    [   // degree 1
        +1.0
    ],
    [   // degree 2
        +1.3465553856377803,
        -3.4655538563778032e-1
    ],
    [   // degree 3
        +1.4228653756681227,
        -5.8208556916449616e-1,
        +1.5922019349637218e-1
    ],
    [   // degree 4
        +1.4387257478171547,
        -6.7778401359918661e-1,
        +3.2118898377713379e-1,
        -8.2130717995088531e-2
    ],
    [   // degree 5
        +1.4419170408633741,
        -7.0909645927612530e-1,
        +4.1560609399164150e-1,
        -1.9357573729558908e-1,
        +4.5149061716699634e-2
    ],
    [   // degree 6
        +1.4425449435950917,
        -7.1814525675038965e-1,
        +4.5754919692564044e-1,
        -2.7790534462849337e-1,
        +1.2179791068763279e-1,
        -2.5841449829670182e-2
    ],
    [   // degree 7
        +1.4426664401536078,
        -7.2055423726162360e-1,
        +4.7332419162501083e-1,
        -3.2514018752954144e-1,
        +1.9302965529095673e-1,
        -7.8534970641157997e-2,
        +1.5209108363023915e-2
    ],
    [   // degree 8
        +1.4426896453621882,
        -7.2115893912535967e-1,
        +4.7861716616785088e-1,
        -3.4699935395019565e-1,
        +2.4114048765477492e-1,
        -1.3657398692885181e-1,
        +5.1421382871922106e-2,
        -9.1364020499895560e-3
    ]
];

const C_LOG2_EST_MAX_ERROR: readonly number[] = [
    8.6071332055935e-2,  // degree 1
    7.6362868906659e-3,  // degree 2
    8.7902902652948e-4,  // degree 3
    1.1318551356388e-4,  // degree 4
    1.5521274483455e-5,  // degree 5
    2.2162052037978e-6,  // degree 6
    3.2546558681457e-7,  // degree 7
    4.8798286744756e-8   // degree 8
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
//   const result = log2Estimate(x, 3);
export function log2Estimate(x: number, degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    const coeff = C_LOG2_EST_COEFF[degree - 1];
    const last = degree - 1;
    const t = x - 1;  // t in [0,1]
    let poly = coeff[last];
    for (let i = 0, index = last - 1; i < last; ++i, --index) {
        poly = coeff[index] + poly * t;
    }
    poly = poly * t;
    return poly;
}

// The input constraint is x > 0. Range reduction is used to generate
// a value y in [1,2], call log2Estimate(y) and then add the exponent
// for the power of two in the binary scientific representation of x.
// For example a degree-3 estimate is
//   const x = ...;  // x > 0
//   const result = log2EstimateRR(x, 3);
export function log2EstimateRR(x: number, degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    let p = frexpExponent(x);
    const y = 2 * ldexp(x, -p);  // y in [1,2)
    --p;
    const poly = log2Estimate(y, degree);
    return poly + p;
}

export function getLog2EstimateMaxError(degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    return C_LOG2_EST_MAX_ERROR[degree - 1];
}
