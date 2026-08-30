// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ExpEstimate.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Minimax polynomial approximations to exp(x) of the form f(x) =
// p(x/log(2)), where log(2) is the natural logarithm of 2 and the polynomial
// p(y) of degree D minimizes the quantity maximum{|2^y - p(y)| : y in [0,1]}
// over all polynomials of degree D. The identity exp(x) = 2^{x/log(2)} is
// used.
//
// The upstream C++ selects the degree with a template parameter checked by
// static_assert. This port takes the degree as a runtime argument validated
// by logAssert; see expEstimate below.

import { GTE_C_INV_LN_2 } from './Constants';
import { exp2Estimate, exp2EstimateRR } from './Exp2Estimate';
import { logAssert } from './Logger';

const C_EXP_EST_MAX_ERROR: readonly number[] = [
    8.6071332055935e-2,  // degree 1
    3.8132476831059e-3,  // degree 2
    1.4694877755229e-4,  // degree 3
    4.7617792662269e-6,  // degree 4
    1.3162098766451e-7,  // degree 5
    3.1590550175765e-9,  // degree 6
    6.7157168714971e-11  // degree 7
];

// The input constraint is x in [0,log(2)], where log(2) is the natural
// logarithm of 2. For example a degree-3 estimate is
//   const x = ...;  // in [0,log(2)]
//   const result = expEstimate(x, 3);
export function expEstimate(x: number, degree: number): number {
    logAssert(1 <= degree && degree <= 7 && degree === Math.trunc(degree),
        'Invalid degree.');

    return exp2Estimate(x * GTE_C_INV_LN_2, degree);
}

// The input x can be any real number. Range reduction is used to generate
// a value y in [0,log(2)], call expEstimate(y) and then combine the
// output with the proper exponent to obtain the approximation. For
// example,
//   const x = ...;  // any real number
//   const result = expEstimateRR(x, 3);
export function expEstimateRR(x: number, degree: number): number {
    logAssert(1 <= degree && degree <= 7 && degree === Math.trunc(degree),
        'Invalid degree.');

    return exp2EstimateRR(x * GTE_C_INV_LN_2, degree);
}

export function getExpEstimateMaxError(degree: number): number {
    logAssert(1 <= degree && degree <= 7 && degree === Math.trunc(degree),
        'Invalid degree.');

    return C_EXP_EST_MAX_ERROR[degree - 1];
}
