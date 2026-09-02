// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RootsGeneralPolynomial.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The solve functions return the real-valued roots of the polynomial
//   p(x) = p[0] + p[1] * x + p[2] * x^2 + ... + p[d] * x^d
// The derivative is
//   p'(x) = p[1] + 2 * p[2] * x + ... + d * p[d] * x^{d-1}.
// If r0 and r1 are consecutive roots of p'(x), say r0 < r1, then p(t) is
// monotonic on the open interval (r0,r1). Additionally, if
// p(r0) * p(r1) <= 0, then p(x) has a unique root on the closed interval
// [r0,r1]. Using this observation, one can compute the derivatives through
// order d for p(x), find roots for the derivative of order k+1, and then
// use these to bound roots for the derivative of order k. This is a
// recursive formulation, implemented as recursive function calls.
//
// The old code, now deprecated, is RootsPolynomial.find(...). It uses only
// floating-point arithmetic. The rounding errors in computing the
// coefficients of the polynomial derivatives can be catastrophic, leading to
// extremely inaccurate roots.
//
// The new code uses a mixture of rational arithmetic and floating-point
// arithmetic. The coefficients of the polynomial derivatives are computed
// using rational arithmetic. When it comes time for bisection, intervals are
// located for which the rational polynomial values at the endpoints have
// opposite signs. Rational-valued bisection effectively does not converge
// (the number of bits in a rational is extremely large), so instead the
// polynomial coefficients are rounded to the nearest floating-point numbers
// and the polynomial is evaluated at the endpoints using floating-point
// arithmetic. Special handling is given to the case where the rational
// polynomial values have opposite signs but the floating-point polynomial
// values do not.
//
// A root of even multiplicity is not reported, because p(x) does not change
// sign there; a root of odd multiplicity is reported once, whatever its
// multiplicity. The roots are ordered increasingly.
//
// Port notes:
//   - The upstream 'Rational' is BSRational<UIntegerAP32>; the port uses
//     BSRational (see BSRational.ts), whose arithmetic returns new objects.
//   - The two upstream 'Solve' overloads (floating-point input and rational
//     input) become 'solve' and 'solveRational'; the output std::vector
//     parameters become return values.
//   - 'useThreading' is accepted for API compatibility and ignored: the port
//     is single-threaded. The upstream threaded branch pushes the roots of
//     the subintervals in increasing interval order, exactly as the
//     sequential branch does, so the results are identical.
//   - The upstream Evaluate relies on std::size_t wrap-around as its loop
//     terminator; the port uses a signed index.
//   - Upstream bug fixed (see the PR body): the vector of rational
//     coefficients is allocated with p.size() entries but only degree+1 are
//     assigned, so when high-order zero-valued coefficients are removed the
//     solver runs on a polynomial padded with zeros. The leading coefficient
//     is then zero, which contradicts the monic normalization the Cauchy
//     bound and the recursion assume. The port allocates degree+1 entries.

import { BSRational } from './BSRational';

function ratio(x: number): BSRational {
    return BSRational.fromNumber(x);
}

// The port of Evaluate<T> (Horner's method).
function evaluateNumber(p: readonly number[], x: number): number {
    let i = p.length - 1;
    let result = p[i];
    while (--i >= 0) {
        result = x * result + p[i];
    }
    return result;
}

// The port of Evaluate<Rational>.
function evaluateRational(p: readonly BSRational[], x: BSRational): BSRational {
    let i = p.length - 1;
    let result = p[i].clone();
    while (--i >= 0) {
        result = x.mul(result).add(p[i]);
    }
    return result;
}

// The port of Bisect. The return value is the root, or null when the
// interval has no root to report.
function bisect(tP: readonly number[], rP: readonly BSRational[],
    rXMin: BSRational, rXMax: BSRational): BSRational | null {
    // The first interval is [-cauchyBound,derivativeRoot.first]. It is
    // possible that p'(x) has a root smaller than the minimum root of p(x),
    // in which case the incoming interval endpoints are not correctly
    // ordered. Such an interval cannot produce a root of p(x). The last
    // interval is [derivative.last,+cauchyBound]. It is possible that p'(x)
    // has a root larger than the maximum root of p(x), in which case the
    // incoming interval endpoints are not correctly ordered. Such an interval
    // cannot produce a root of p(x).
    if (rXMin.greaterThanOrEqual(rXMax)) {
        return null;
    }

    const rPMin = evaluateRational(rP, rXMin);
    const signRPMin = rPMin.getSign();
    if (signRPMin === 0) {
        return rXMin.clone();
    }

    const rPMax = evaluateRational(rP, rXMax);
    const signRPMax = rPMax.getSign();
    if (signRPMax === 0) {
        // Do not return the root rXMax. The next interval will be responsible
        // for managing this root.
        return null;
    }

    if (signRPMin * signRPMax > 0) {
        // The polynomial p(x) is monotone on [rXMin,rXMax], so it cannot have
        // a root on the interval.
        return null;
    }

    // At this time rPMin and rPMax have opposite signs. There must be a
    // unique root on [rXMin,rXMax] because the derivative is not zero on the
    // interval, which implies that p(x) is monotone on the interval.

    // Use floating-point arithmetic for speed. Be aware that the conversions
    // from rational numbers to floating-point numbers can affect sign tests.
    // Recompute the endpoint tests for the floating-point numbers.
    const zero = 0;

    let tXMin = rXMin.toNumber();
    const tPMin = evaluateNumber(tP, tXMin);
    const signTPMin = (tPMin > zero ? +1 : (tPMin < zero ? -1 : 0));
    if (signTPMin === 0) {
        return rXMin.clone();
    }

    let tXMax = rXMax.toNumber();
    const tPMax = evaluateNumber(tP, tXMax);
    const signTPMax = (tPMax > zero ? +1 : (tPMax < zero ? -1 : 0));
    if (signTPMax === 0) {
        // Do not return the root rXMax. The next interval will be responsible
        // for managing this root.
        return null;
    }

    if (signTPMin * signTPMax > 0) {
        // We know that rPMin and rPMax have opposite signs. Rounding errors
        // lead to tPMin and tPMax having the same sign. Rather than return
        // 'null' as in the rational arithmetic case, return the root given by
        // the intersection of the x-axis with the line through (rXMin,rPMin)
        // and (rXMax,rPMax). This amounts to approximating a nearly flat
        // polynomial on the interval by a line segment.
        return rXMin.mul(rPMax).sub(rXMax.mul(rPMin)).div(rPMax.sub(rPMin));
    }

    // At this time tPMin and tPMax have opposite signs. Bisect to find a
    // root. In theory the root is unique, but floating-point rounding errors
    // can lead to multiple roots (all approximately the same floating-point
    // number).

    // The maximum number of iterations suffices for convergence when using
    // floating-point numbers.
    const maxIterations = 4096;
    const tHalf = 0.5;
    let tRoot = 0;
    for (let i = 0; i < maxIterations; ++i) {
        tRoot = tHalf * (tXMin + tXMax);

        // The test is designed for floating-point numbers when tXMin and
        // tXMax are consecutive.
        if (tRoot === tXMin || tRoot === tXMax) {
            break;
        }

        const tPAtRoot = evaluateNumber(tP, tRoot);
        const signTPAtRoot = (tPAtRoot > zero ? +1 : (tPAtRoot < zero ? -1 : 0));
        const sign = signTPAtRoot * signTPMin;
        if (sign < 0) {
            tXMax = tRoot;
        } else if (sign > 0) {
            tXMin = tRoot;
        } else {
            // The root is exactly tRoot.
            break;
        }
    }

    return ratio(tRoot);
}

// The port of SolveRecursive. The roots are appended to rRoots.
function solveRecursive(rP: readonly BSRational[], rXMin: BSRational,
    rXMax: BSRational, rRoots: BSRational[]): void {
    // The base of the recursion.
    const degree = rP.length - 1;
    if (degree === 1) {
        if (rP[1].getSign() !== 0) {
            rRoots.push(rP[0].negated().div(rP[1]));
        }
        return;
    }

    // Compute the derivative polynomial p'(x) of p(x) using rational numbers.
    const rPDerivative = new Array<BSRational>(degree);
    for (let i0 = 0, i1 = 1; i1 <= degree; i0 = i1++) {
        rPDerivative[i0] = rP[i1].mul(ratio(i1));
    }

    // Estimate the roots of the derivative polynomial.
    const rRootsDerivative: BSRational[] = [];
    solveRecursive(rPDerivative, rXMin, rXMax, rRootsDerivative);

    // Round the coefficients of rP(x) to floating-point numbers. This is used
    // for fast performance by floating-point-based bisection.
    const tP = new Array<number>(rP.length);
    for (let i = 0; i < rP.length; ++i) {
        tP[i] = rP[i].toNumber();
    }

    // The polynomial is monotonic between consecutive roots of the
    // derivative. This feature and the polynomial values at the derivative
    // roots are used to compute polynomial roots via bisection.
    if (rRootsDerivative.length > 0) {
        // Let rXUpper = rRootsDerivative[0]. Estimate a root, if any, on the
        // interval [rXMin,rXUpper].
        let rRoot = bisect(tP, rP, rXMin, rRootsDerivative[0]);
        if (rRoot !== null) {
            rRoots.push(rRoot);
        }

        // Let rXLower = rRootsDerivative[i0] and let rXUpper =
        // rRootsDerivative[i1]. Estimate a root, if any, on
        // [rXLower,rXUpper].
        for (let i0 = 0, i1 = 1; i1 < rRootsDerivative.length; i0 = i1++) {
            rRoot = bisect(tP, rP, rRootsDerivative[i0], rRootsDerivative[i1]);
            if (rRoot !== null) {
                rRoots.push(rRoot);
            }
        }

        // Let rXLower be the last derivative root. Estimate a root, if any,
        // on the interval [rXLower,rXMax].
        rRoot = bisect(tP, rP,
            rRootsDerivative[rRootsDerivative.length - 1], rXMax);
        if (rRoot !== null) {
            rRoots.push(rRoot);
        }
    } else {
        // The polynomial is monotone on [rXMin,rXMax], so it has at most one
        // root.
        const rRoot = bisect(tP, rP, rXMin, rXMax);
        if (rRoot !== null) {
            rRoots.push(rRoot);
        }
    }
}

// The port of InitiateSolver.
function initiateSolver(rP: readonly BSRational[], rRoots: BSRational[]): void {
    // Compute Cauchy bounds to obtain an interval containing the roots of
    // p(x). At this time the polynomial is monic.
    const degree = rP.length - 1;
    let rCauchyBound = BSRational.fabs(rP[0]);
    for (let i = 1; i + 1 <= degree; ++i) {
        const rCandidate = BSRational.fabs(rP[i]);
        if (rCandidate.greaterThan(rCauchyBound)) {
            rCauchyBound = rCandidate;
        }
    }
    rCauchyBound = rCauchyBound.add(ratio(1));

    // Solve recursively in degree.
    solveRecursive(rP, rCauchyBound.negated(), rCauchyBound, rRoots);
}

export class RootsGeneralPolynomial {
    // The port of Solve(std::vector<T> const&, bool, std::vector<T>&). The
    // input p stores the coefficients of the polynomial in increasing order
    // of power. The returned roots are ordered increasingly and are distinct.
    static solve(p: readonly number[], useThreading: boolean = false): number[] {
        // The port is single-threaded; see the port notes.
        void useThreading;

        // The order of p is p.length. The degree of p is p.length - 1.
        if (p.length <= 1) {
            // The polynomial is identically a constant. Do not report roots
            // even when that constant is 0.
            return [];
        }

        // Remove high-order zero-valued coefficients.
        let degree = p.length - 1;
        while (degree >= 1 && p[degree] === 0) {
            --degree;
        }

        if (degree === 1) {
            return [-p[0] / p[1]];
        }

        if (degree === 0) {
            // The polynomial is identically a constant. Do not report roots
            // even when that constant is 0.
            return [];
        }

        // At this time the degree is at least 2. Create a polynomial for p(x)
        // that has rational coefficients. Upstream allocates p.size() entries
        // rather than degree+1; see the port notes.
        const rP = new Array<BSRational>(degree + 1);
        for (let i = 0; i <= degree; ++i) {
            rP[i] = ratio(p[i]);
        }

        // Make the polynomial monic. Theoretically, this is irrelevant when
        // estimating roots for a polynomial with rational coefficients.
        // However, during the recursion the rational coefficients can be
        // quite large, so using a monic polynomial helps with robustness.
        const rOne = ratio(1);
        const rLast = rP[degree];
        if (rLast.notEquals(rOne)) {
            for (let i = 0; i < degree; ++i) {
                rP[i] = rP[i].div(rLast);
            }
            rP[degree] = rOne;
        }

        // Compute Cauchy bounds and solve for roots using recursion on the
        // polynomial degree.
        const rRoots: BSRational[] = [];
        initiateSolver(rP, rRoots);

        // Convert the rational roots to floating-point.
        return rRoots.map(r => r.toNumber());
    }

    // The port of Solve(std::vector<Rational> const&, bool,
    // std::vector<Rational>&).
    static solveRational(rP: readonly BSRational[],
        useThreading: boolean = false): BSRational[] {
        // The port is single-threaded; see the port notes.
        void useThreading;

        // The order of p is rP.length. The degree of p is rP.length - 1.
        if (rP.length <= 1) {
            // The polynomial is identically a constant. Do not report roots
            // even when that constant is 0.
            return [];
        }

        // Remove high-order zero-valued coefficients.
        let degree = rP.length - 1;
        while (degree >= 1 && rP[degree].getSign() === 0) {
            --degree;
        }

        if (degree === 1) {
            return [rP[0].negated().div(rP[1])];
        }

        if (degree === 0) {
            // The polynomial is identically a constant. Do not report roots
            // even when that constant is 0.
            return [];
        }

        // At this time the degree is at least 2.

        // Make the polynomial monic; see the comment in solve(). Upstream
        // allocates rP.size() entries rather than degree+1; see the port
        // notes.
        const rOne = ratio(1);
        const rPMonic = new Array<BSRational>(degree + 1);
        const rLast = rP[degree];
        if (rLast.notEquals(rOne)) {
            for (let i = 0; i < degree; ++i) {
                rPMonic[i] = rP[i].div(rLast);
            }
            rPMonic[degree] = rOne;
        } else {
            for (let i = 0; i <= degree; ++i) {
                rPMonic[i] = rP[i].clone();
            }
        }

        // Compute Cauchy bounds and solve for roots using recursion on the
        // polynomial degree.
        const rRoots: BSRational[] = [];
        initiateSolver(rPMonic, rRoots);
        return rRoots;
    }
}
