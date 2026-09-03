// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ChebyshevRatioEstimate.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Chebyshev ratio is f(t,A) = sin(t*A)/sin(A) for t in [0,1] and A in
// [0,pi/2]. Let x = cos(A) and y = 1 - cos(A), both in [0,1]. As a function
// of y, a series representation for f(t,y) is
//   f(t,y) = sum_{i=0}^{infinity} c_{i}(t) y^{i}
// where c_0(t) = t, c_{i}(t) = c_{i-1}(t)*(i^2 - t^2)/(i*(2*i+1)) for i >= 1.
// The c_{i}(t) are polynomials in t of degree 2*i+1. The document
// https://www.geometrictools.com/Documentation/FastAndAccurateSlerp.pdf
// derives an approximation
//   g(t,y) = sum_{i=0}^{n-1} c_{i}(t) y^{i} + u_n c_{n}(t) y^n
// which has degree 2*n+1 in t and degree n in y. The constants u_n are chosen
// for balanced error bounds. chebyshevRatioEstimate implements this
// algorithm. If the angle A is restricted to [0,pi/4], then better estimates
// are obtained for the same computational cost. All that differs are the
// u_n-values. chebyshevRatioEstimateR implements this algorithm. The
// functions return pairs {f(1-t,A), f(t,A)}, which is useful for spherical
// linear interpolation.
//
// The upstream C++ selects the degree with a template parameter checked by
// static_assert. This port takes the degree as a runtime argument validated
// by logAssert.

import { logAssert } from './Logger.js';

// Constants for chebyshevRatioEstimate.
const C_CHBRAT_EST_U: readonly number[] = [
    1.5149656562200644050,
    1.6410179946672027729,
    1.7124880779005808851,
    1.7593545031636841358,
    1.7927054757060019163,
    1.8177479632959470113,
    1.8372872973294931409,
    1.8529805143706497006,
    1.8658739107798316681,
    1.8766626700393858052,
    1.8858276947289707159,
    1.8937127486228939599,
    1.9005703533887863266,
    1.9065903281211855624,
    1.9119182032942771965,
    1.9166674811124804201
];

const C_CHBRAT_EST_MAX_ERROR: readonly number[] = [
    1.8249897492955e-2,
    5.2760601519022e-3,
    1.8055057987877e-3,
    6.7244299646175e-4,
    2.6386437427495e-4,
    1.0731422197408e-4,
    4.4805894183764e-5,
    1.9088088593749e-5,
    8.2629028074211e-6,
    3.6237273527418e-6,
    1.6064797200289e-6,
    7.1872518425665e-7,
    3.2407757655229e-7,
    1.4712279927665e-7,
    6.7187475472075e-8,
    3.0844086507110e-8
];

// Constants for chebyshevRatioEstimateR.
const C_CHBRAT_ESTR_U: readonly number[] = [
    1.1021472152138613865,
    1.1239349540626744073,
    1.1351870374370363059,
    1.1421060160698368602,
    1.1468020192623136211,
    1.1502017494201659531,
    1.1527782928466798751,
    1.1547990001678465344,
    1.1564265502929687024,
    1.1577657226562501069,
    1.1588859375000000185,
    1.1598375000000000767
];

const C_CHBRAT_ESTR_MAX_ERROR: readonly number[] = [
    8.6832275204274e-4,
    6.6040175097815e-5,
    6.1949661303018e-6,
    6.4578503422564e-7,
    7.1792162659179e-8,
    8.3364721792379e-9,
    9.9903230132981e-10,
    1.2262002524466e-10,
    1.5335510639148e-11,
    1.9472201628901e-12,
    2.5046631435544e-13,
    3.2751579226443e-14
];

// Shared evaluator for the two estimate variants. The upstream C++ builds
// the a[] and b[] coefficient arrays with constexpr generators
// C_CHBRAT_ACOEFF/C_CHBRAT_BCOEFF; here the coefficients are computed in the
// loop. Only the last term (i = degree - 1) is scaled by the u-value.
function evaluate(t: number, x: number, degree: number, u: readonly number[]): [number, number] {
    const y = 1 - x;
    let term0 = 1 - t;
    let term1 = t;
    const sqr0 = term0 * term0;
    const sqr1 = term1 * term1;
    const f: [number, number] = [term0, term1];
    for (let i = 0; i < degree; ++i) {
        const uFactor = (degree !== i + 1 ? 1 : u[i]);
        const a = uFactor / ((i + 1) * (2 * (i + 1) + 1));
        const b = uFactor * (i + 1) / (2 * (i + 1) + 1);
        term0 *= (b - a * sqr0) * y;
        term1 *= (b - a * sqr1) * y;
        f[0] += term0;
        f[1] += term1;
    }
    return f;
}

// Compute estimates for f(t,x) = sin(t*A)/sin(A), where t in [0,1],
// A in [0,pi/2], x = cos(A) in [0,1]; the returned pair contains the
// estimate for f(1-t,x) and the estimate for f(t,x). The approximating
// function is a polynomial of two variables. The degree must be in {1..16}.
// The degree in t is 2*degree+1 and the degree in x is degree.
export function chebyshevRatioEstimate(t: number, x: number, degree: number): [number, number] {
    logAssert(1 <= degree && degree <= 16 && degree === Math.trunc(degree),
        'Invalid degree.');

    return evaluate(t, x, degree, C_CHBRAT_EST_U);
}

export function getChebyshevRatioEstimateMaxError(degree: number): number {
    logAssert(1 <= degree && degree <= 16 && degree === Math.trunc(degree),
        'Invalid degree.');

    return C_CHBRAT_EST_MAX_ERROR[degree - 1];
}

// Compute estimates for f(t,x) = sin(t*A)/sin(A), where t in [0,1],
// A in [0,pi/4], x = cos(A) in [sqrt(2)/2,1]; the returned pair contains
// the estimate for f(1-t,x) and the estimate for f(t,x). The approximating
// function is a polynomial of two variables. The degree must be in {1..12}.
// The degree in t is 2*degree+1 and the degree in x is degree.
export function chebyshevRatioEstimateR(t: number, x: number, degree: number): [number, number] {
    logAssert(1 <= degree && degree <= 12 && degree === Math.trunc(degree),
        'Invalid degree.');

    return evaluate(t, x, degree, C_CHBRAT_ESTR_U);
}

export function getChebyshevRatioEstimateRMaxError(degree: number): number {
    logAssert(1 <= degree && degree <= 12 && degree === Math.trunc(degree),
        'Invalid degree.');

    return C_CHBRAT_ESTR_MAX_ERROR[degree - 1];
}
