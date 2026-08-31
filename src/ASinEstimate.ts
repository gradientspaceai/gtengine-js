// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ASinEstimate.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Approximations to asin(x) of the form f(x) = pi/2 - sqrt(1-x)*p(x), where
// the polynomial p(x) of degree D minimizes the quantity
// maximum{|acos(x)/sqrt(1-x) - p(x)| : x in [0,1]} over all polynomials of
// degree D. The identity asin(x) = pi/2 - acos(x) is used.
//
// The upstream C++ selects the degree with a template parameter checked by
// static_assert. This port takes the degree as a runtime argument validated
// by logAssert; see asinEstimate below.

import { acosEstimate } from './ACosEstimate';
import { GTE_C_HALF_PI } from './Constants';
import { logAssert } from './Logger';

// These are the upstream C_ASIN_EST_MAX_ERROR values. Because
// asin(x) - asinEstimate(x) = acosEstimate(x) - acos(x), the bounds are
// mathematically identical to those in ACosEstimate; the two upstream tables
// differ only in the last one or two decimal digits, which is noise in the
// minimax solver rather than a difference in the estimates. The published
// values are preserved.
const C_ASIN_EST_MAX_ERROR: readonly number[] = [
    9.0128265558586e-3,  // degree 1
    8.1851275863202e-4,  // degree 2
    8.8200141836567e-5,  // degree 3
    1.0563052499871e-5,  // degree 4
    1.3535063235066e-6,  // degree 5
    1.8169471743823e-7,  // degree 6
    2.5231622315797e-8,  // degree 7
    3.5952707963527e-9   // degree 8
];

// The input constraint is x in [0,1]. For example a degree-3 estimate is
//   const x = ...;  // in [0,1]
//   const result = asinEstimate(x, 3);
export function asinEstimate(x: number, degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    return GTE_C_HALF_PI - acosEstimate(x, degree);
}

export function getASinEstimateMaxError(degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    return C_ASIN_EST_MAX_ERROR[degree - 1];
}
