// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RootsQuadratic.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the real-valued roots of a quadratic polynomial with real-valued
// coefficients. The general quadratic polynomial is
//   g(x) = g0 + g1 * x + g2 * x^2
// where g2 is not zero. The monic quadratic polynomial is
//   m(x) = m0 + m1 * x + x^2
// The depressed quadratic polynomial is
//   d(x) = d0 + x^2
// The classification of roots and multiplicities is performed using rational
// arithmetic for exactness. For algorithmic details, see
// https://www.geometrictools.com/Documentation/LowDegreePolynomialRoots.pdf
//
// The code uses bisection on bounding intervals for roots. For a polynomial
// of degree n, Lagrange's bound is
//   b = max(1,|p[0]/p[n]|, |p[1]/p[n]|, ..., |p[n-1]/p[n]|)
// The real roots lie in the interval [-b,b].
//
// Port notes:
//   - See RootsLinear.ts for the RootsScalarOps<T, R> mechanism that stands
//     in for the C++ template parameter T; the 'double' instantiation is
//     RootsQuadratic.solve/solveMonic/solveDepressed and the Rational one is
//     the *Rational variants.
//   - The upstream out-parameter 'PolynomialRoot<T>* roots' plus the size_t
//     count become a returned array of roots.
//   - The three upstream 'Solve' overloads are distinguished by arity in C++;
//     the port names them solve (general), solveMonic and solveDepressed.

import { BSRational } from './BSRational';
import { fma } from './Functions';
import { PolynomialRoot, polynomialRootBisect } from './PolynomialRoot';
import {
    PolynomialRootRational, type RootsScalarOps, rootsNumberOps, rootsRationalOps,
    solveLinearGeneric, solveMonicLinearGeneric, sortRoots
} from './RootsLinear';

const HALF = BSRational.fromNumber(0.5);

// The port of ComputeClassifiers(rG0, rG1, rG2, rD0, rM1Div2).
function computeClassifiersGeneral(rG0: BSRational, rG1: BSRational, rG2: BSRational):
    { rD0: BSRational, rM1Div2: BSRational } {
    const rM0 = rG0.div(rG2);
    const rM1 = rG1.div(rG2);
    return computeClassifiersMonic(rM0, rM1);
}

// The port of ComputeClassifiers(rM0, rM1, rD0, rM1Div2).
function computeClassifiersMonic(rM0: BSRational, rM1: BSRational):
    { rD0: BSRational, rM1Div2: BSRational } {
    const rM1Div2 = HALF.mul(rM1);
    const rD0 = rM0.sub(rM1Div2.mul(rM1Div2));
    return { rD0: rD0, rM1Div2: rM1Div2 };
}

// The port of ComputeDepressedRootsBisection.
function computeDepressedRootsBisection(rD0: BSRational): PolynomialRootRational[] {
    const signD0 = rD0.getSign();
    if (signD0 > 0) {
        // Two non-real roots, each multiplicity 1.
        return [];
    }

    if (signD0 === 0) {
        // One real root, multiplicity 2.
        return [new PolynomialRootRational(BSRational.fromNumber(0), 2)];
    }

    // Two real roots, each multiplicity 1. The Cauchy bound for F(x) is
    // b = max{1,|d0|}. Use bisection on the interval [-b,b] to estimate the
    // roots.
    const d0 = rD0.toNumber();
    const b = Math.max(1.0, Math.abs(d0));
    const F = (x: number) => fma(x, x, d0);

    // Bisect on the interval [0,b]. The polynomial is an even function, so we
    // do not have to bisect on the interval [-b,0].
    const bisected = polynomialRootBisect(F, -1, +1, 0.0, b);
    const average = HALF.mul(BSRational.fromNumber(bisected.xMin).add(
        BSRational.fromNumber(bisected.xMax)));
    return [
        new PolynomialRootRational(average.negated(), 1),
        new PolynomialRootRational(average, 1)
    ];
}

// The port of ComputeDepressedRootsClosedForm.
function computeDepressedRootsClosedForm(rD0: BSRational): PolynomialRootRational[] {
    const signD0 = rD0.getSign();
    if (signD0 > 0) {
        // Two non-real roots, each multiplicity 1.
        return [];
    }

    if (signD0 === 0) {
        // One real root, multiplicity 2.
        return [new PolynomialRootRational(BSRational.fromNumber(0), 2)];
    }

    // Two real roots, each multiplicity 1. Use the closed-form representation
    // of the roots.
    const rSqrtNegD0 = BSRational.sqrt(rD0.negated());
    return [
        new PolynomialRootRational(rSqrtNegD0.negated(), 1),
        new PolynomialRootRational(rSqrtNegD0, 1)
    ];
}

// The port of RootsQuadratic<T>::ComputeDepressedRoots, which does not depend
// on T.
export function computeDepressedQuadraticRoots(useBisection: boolean, rD0: BSRational):
    PolynomialRootRational[] {
    if (useBisection) {
        return computeDepressedRootsBisection(rD0);
    } else {
        return computeDepressedRootsClosedForm(rD0);
    }
}

// The non-negative square root of rV expressed through the depressed-
// quadratic solver, which upstream writes inline as
//   RootsQuadratic<Rational>::ComputeDepressedRoots(useBisection, -rV, rQRoots);
//   Rational rSqrt = rQRoots[1].x;
// Upstream reuses a single rQRoots array across such calls, so when rV is
// zero (one root written, at index 0) or negative (no roots written) the
// value read from index 1 is whatever a previous call left there. The port
// returns 0 in both of those cases: for rV = 0 that is the mathematically
// correct square root (the upstream stale read is a latent bug), and for
// rV < 0 the value has no real square root and every upstream caller
// discards it.
export function rationalSqrtViaQuadratic(useBisection: boolean, rV: BSRational): BSRational {
    if (rV.getSign() <= 0) {
        return BSRational.fromNumber(0);
    }
    const rQRoots = computeDepressedQuadraticRoots(useBisection, rV.negated());
    return rQRoots[1].x;
}

// The port of HasZeroValuedRoots(g1, g2, roots).
function hasZeroValuedRootsGeneral<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, g1: T, g2: T): R[] {
    if (ops.isZero(g1)) {
        return [ops.makeRoot(ops.zero(), 2)];
    } else {
        const roots = solveLinearGeneric(ops, g1, g2);
        roots.push(ops.makeRoot(ops.zero(), 1));
        sortRoots(ops, roots);
        return roots;
    }
}

// The port of HasZeroValuedRoots(m1, roots).
function hasZeroValuedRootsMonic<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, m1: T): R[] {
    if (ops.isZero(m1)) {
        return [ops.makeRoot(ops.zero(), 2)];
    } else {
        const roots = solveMonicLinearGeneric(ops, m1);
        roots.push(ops.makeRoot(ops.zero(), 1));
        sortRoots(ops, roots);
        return roots;
    }
}

// The port of RootsQuadratic<T>::Solve(useBisection, g0, g1, g2, roots).
export function solveQuadraticGeneric<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, useBisection: boolean, g0: T, g1: T, g2: T): R[] {
    // Test whether the degree is smaller than 2.
    if (ops.isZero(g2)) {
        return solveLinearGeneric(ops, g0, g1);
    }

    // Test for zero-valued roots.
    if (ops.isZero(g0)) {
        return hasZeroValuedRootsGeneral(ops, g1, g2);
    }

    // At this time g0 and g2 are not zero. Transform the general quadratic to
    // a depressed quadratic, solve for its roots, and inverse transform them
    // to roots of the general quadratic.
    const c = computeClassifiersGeneral(ops.toRational(g0), ops.toRational(g1),
        ops.toRational(g2));
    const rRoots = computeDepressedQuadraticRoots(useBisection, c.rD0);
    return rRoots.map(r => ops.makeRoot(ops.fromRational(r.x.sub(c.rM1Div2)), r.m));
}

// The port of RootsQuadratic<T>::Solve(useBisection, m0, m1, roots).
export function solveMonicQuadraticGeneric<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, useBisection: boolean, m0: T, m1: T): R[] {
    // Test for zero-valued roots.
    if (ops.isZero(m0)) {
        return hasZeroValuedRootsMonic(ops, m1);
    }

    // At this time m0 is not zero. Transform the monic quadratic to a
    // depressed quadratic, solve for its roots, and inverse transform them to
    // roots of the monic quadratic.
    const c = computeClassifiersMonic(ops.toRational(m0), ops.toRational(m1));
    const rRoots = computeDepressedQuadraticRoots(useBisection, c.rD0);
    return rRoots.map(r => ops.makeRoot(ops.fromRational(r.x.sub(c.rM1Div2)), r.m));
}

// The port of RootsQuadratic<T>::Solve(useBisection, d0, roots).
export function solveDepressedQuadraticGeneric<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, useBisection: boolean, d0: T): R[] {
    // The quadratic is already depressed, so no transforming and inverse
    // transforming are necessary.
    const rRoots = computeDepressedQuadraticRoots(useBisection, ops.toRational(d0));
    return rRoots.map(r => ops.makeRoot(ops.fromRational(r.x), r.m));
}

export class RootsQuadratic {
    // Solve the general quadratic g0 + g1*x + g2*x^2 = 0.
    static solve(useBisection: boolean, g0: number, g1: number, g2: number): PolynomialRoot[] {
        return solveQuadraticGeneric(rootsNumberOps, useBisection, g0, g1, g2);
    }

    // Solve the monic quadratic m0 + m1*x + x^2 = 0.
    static solveMonic(useBisection: boolean, m0: number, m1: number): PolynomialRoot[] {
        return solveMonicQuadraticGeneric(rootsNumberOps, useBisection, m0, m1);
    }

    // Solve the depressed quadratic d0 + x^2 = 0.
    static solveDepressed(useBisection: boolean, d0: number): PolynomialRoot[] {
        return solveDepressedQuadraticGeneric(rootsNumberOps, useBisection, d0);
    }

    // Solve the general quadratic with rational coefficients.
    static solveRational(useBisection: boolean, g0: BSRational, g1: BSRational,
        g2: BSRational): PolynomialRootRational[] {
        return solveQuadraticGeneric(rootsRationalOps, useBisection, g0, g1, g2);
    }

    // Solve the monic quadratic with rational coefficients.
    static solveMonicRational(useBisection: boolean, m0: BSRational,
        m1: BSRational): PolynomialRootRational[] {
        return solveMonicQuadraticGeneric(rootsRationalOps, useBisection, m0, m1);
    }

    // Solve the depressed quadratic with rational coefficients.
    static solveDepressedRational(useBisection: boolean, d0: BSRational): PolynomialRootRational[] {
        return solveDepressedQuadraticGeneric(rootsRationalOps, useBisection, d0);
    }

    // The upstream public ComputeDepressedRoots, which does not depend on T.
    static computeDepressedRoots(useBisection: boolean, rD0: BSRational): PolynomialRootRational[] {
        return computeDepressedQuadraticRoots(useBisection, rD0);
    }
}
