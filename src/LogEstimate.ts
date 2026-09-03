// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) LogEstimate.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Minimax polynomial approximations to log(x) of the form f(x) =
// p(x)*log(2), where log(2) is the natural logarithm of 2 and the polynomial
// p(x) of degree D minimizes the quantity
// maximum{|log2(x) - p(x)| : x in [1,2]} over all polynomials of degree D.
// The identity log(x) = log2(x)*log(2) is used.
//
// The upstream C++ selects the degree with a template parameter checked by
// static_assert. This port takes the degree as a runtime argument validated
// by logAssert; see logEstimate below.

import { GTE_C_LN_2 } from './Constants.js';
import { getLog2EstimateMaxError, log2Estimate, log2EstimateRR } from './Log2Estimate.js';
import { logAssert } from './Logger.js';

// The input constraint is x in [1,2]. For example a degree-3 estimate is
//   const x = ...;  // in [1,2]
//   const result = logEstimate(x, 3);
export function logEstimate(x: number, degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    return log2Estimate(x, degree) * GTE_C_LN_2;
}

// The input constraint is x > 0. Range reduction is used to generate a value
// y in [1,2], call logEstimate(y) and then add the exponent for the power of
// two in the binary scientific representation of x. For example a degree-3
// estimate is
//   const x = ...;  // x > 0
//   const result = logEstimateRR(x, 3);
export function logEstimateRR(x: number, degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    return log2EstimateRR(x, degree) * GTE_C_LN_2;
}

// Upstream GetLogEstimateMaxError forwards to GetLog2EstimateMaxError, so the
// reported bound is the log2 error rather than the log error. Since
// log(x) - logEstimate(x) = log(2) * (log2(x) - log2Estimate(x)), the tight
// bound is log(2) times smaller. The upstream value is a valid but loose
// bound and is preserved here (compare ExpEstimate.h, which carries its own
// scaled table rather than forwarding to Exp2Estimate).
export function getLogEstimateMaxError(degree: number): number {
    logAssert(1 <= degree && degree <= 8 && degree === Math.trunc(degree),
        'Invalid degree.');

    return getLog2EstimateMaxError(degree);
}
