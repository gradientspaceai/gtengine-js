// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Exp2Estimate.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Minimax polynomial approximations to 2^x. The polynomial p(x) of
// degree D minimizes the quantity maximum{|2^x - p(x)| : x in [0,1]}
// over all polynomials of degree D.
//
// The upstream C++ selects the degree with a template parameter checked by
// static_assert. This port takes the degree as a runtime argument validated
// by logAssert; see exp2Estimate below.

import { logAssert } from './Logger';

const C_EXP2_EST_COEFF: readonly (readonly number[])[] = [
    [   // degree 1
        1.0,
        1.0
    ],
    [   // degree 2
        1.0,
        6.5571332605741528e-1,
        3.4428667394258472e-1
    ],
    [   // degree 3
        1.0,
        6.9589012084456225e-1,
        2.2486494900110188e-1,
        7.9244930154334980e-2
    ],
    [   // degree 4
        1.0,
        6.9300392358459195e-1,
        2.4154981722455560e-1,
        5.1744260331489045e-2,
        1.3701998859367848e-2
    ],
    [   // degree 5
        1.0,
        6.9315298010274962e-1,
        2.4014712313022102e-1,
        5.5855296413199085e-2,
        8.9477503096873079e-3,
        1.8968500441332026e-3
    ],
    [   // degree 6
        1.0,
        6.9314698914837525e-1,
        2.4023013440952923e-1,
        5.5481276898206033e-2,
        9.6838443037086108e-3,
        1.2388324048515642e-3,
        2.1892283501756538e-4
    ],
    [   // degree 7
        1.0,
        6.9314718588750690e-1,
        2.4022637363165700e-1,
        5.5505235570535660e-2,
        9.6136265387940512e-3,
        1.3429234504656051e-3,
        1.4299202757683815e-4,
        2.1662892777385423e-5
    ]
];

const C_EXP2_EST_MAX_ERROR: readonly number[] = [
    8.6071332055935e-2,  // degree 1
    3.8132476831059e-3,  // degree 2
    1.4694877755229e-4,  // degree 3
    4.7617792662269e-6,  // degree 4
    1.3162098788655e-7,  // degree 5
    3.1590552396211e-9,  // degree 6
    6.7157390759576e-11  // degree 7
];

// Equivalent of std::ldexp(value, power): value * 2^power with correct
// handling of overflow to infinity and underflow toward the subnormals.
// Math.pow(2, power) is exact for integer power in [-1074, 1023], so a
// single multiplication is exact rounding there; the extreme ranges are
// handled by scaling in two steps.
function ldexp(value: number, power: number): number {
    if (power >= -1022) {
        if (power <= 1023) {
            return value * Math.pow(2, power);
        }
        // Overflow territory unless value is tiny.
        return value * Math.pow(2, 1023) * Math.pow(2, Math.min(power - 1023, 1023));
    }
    // Subnormal territory unless value is huge. The first factor keeps the
    // intermediate product exact for normal inputs.
    return value * Math.pow(2, -1022) * Math.pow(2, Math.max(power + 1022, -1074));
}

// The input constraint is x in [0,1]. For example a degree-3 estimate is
//   const x = ...;  // in [0,1]
//   const result = exp2Estimate(x, 3);
export function exp2Estimate(x: number, degree: number): number {
    logAssert(1 <= degree && degree <= 7 && degree === Math.trunc(degree),
        'Invalid degree.');

    const coeff = C_EXP2_EST_COEFF[degree - 1];
    const last = degree;
    let poly = coeff[last];
    for (let i = 0, index = last - 1; i < last; ++i, --index) {
        poly = coeff[index] + poly * x;
    }
    return poly;
}

// The input x can be any real number. Range reduction is used to generate
// a value y in [0,1], call exp2Estimate(y) and then combine the output
// with the proper exponent to obtain the approximation. For example a
// degree-3 estimate is
//   const x = ...;  // any real number
//   const result = exp2EstimateRR(x, 3);
export function exp2EstimateRR(x: number, degree: number): number {
    const p = Math.floor(x);
    const y = x - p;
    const poly = exp2Estimate(y, degree);
    return ldexp(poly, p);
}

export function getExp2EstimateMaxError(degree: number): number {
    logAssert(1 <= degree && degree <= 7 && degree === Math.trunc(degree),
        'Invalid degree.');

    return C_EXP2_EST_MAX_ERROR[degree - 1];
}
