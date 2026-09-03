// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RootsCubic.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the real-valued roots of a cubic polynomial with real-valued
// coefficients. The general cubic polynomial is
//   g(x) = g0 + g1 * x + g2 * x^2 + g3 * x^3
// where g3 is not zero. The monic cubic polynomial is
//   m(x) = m0 + m1 * x + m2 * x^2 + x^3
// The depressed cubic polynomial is
//   d(x) = d0 + d1 * x + x^3
// The classification of roots and multiplicities is performed using rational
// arithmetic for exactness. For algorithmic details, see
// https://www.geometrictools.com/Documentation/LowDegreePolynomialRoots.pdf
//
// The code uses bisection on bounding intervals for roots. For a polynomial
// of degree n with all real roots, Samuelson's inequality
// https://en.wikipedia.org/wiki/Samuelson%27s_inequality
// provides a bounding interval [b0,b1] where
//   b0 = (-p[n-1] - (n-1) * s) / (n * p[n])
//   b1 = (-p[n-1] + (n-1) * s) / (n * p[n])
//   s = sqrt(p[n-1]^2 - 2 * n * p[n] * p[n-2] / (n-1))
// Applied to the general cubic,
//   b0 = (-p2 - 2 * s) / (3 * p3)
//   b1 = (-p2 + 2 * s) / (3 * p3)
//   s = sqrt(p2^2 - 3 * p3 * p1)
// Applied to the depressed cubic,
//   b0 = -sqrt(-4 * p1 / 3)
//   b1 = +sqrt(-4 * p1 / 3)
//
// For a monic polynomial of degree n, Cauchy's bound is
//   b = 1 + max(|p[0]|, |p[1]|, ..., |p[n-1]|)
// The real roots lie in the interval [-b,b].
//
// Port notes:
//   - See RootsLinear.ts for the RootsScalarOps<T, R> mechanism that stands
//     in for the C++ template parameter T.
//   - The upstream out-parameter 'PolynomialRoot<T>* roots' plus the size_t
//     count become a returned array of roots.
//   - Upstream bug fixed (see the PR body): the one-real-root branch of
//     ComputeDepressedRootsBisection (negative discriminant) uses
//     b = max(1,|d0|,|d1|) as a root bound. That is not a valid bound for
//     x^3 + d1*x + d0; for example x^3 - 0.9*x - 0.9 has b = 1 while its only
//     real root is about 1.2695, so F(b) has the wrong sign and
//     PolynomialRootBisect returns the endpoint b as the "root". The port
//     uses Cauchy's bound b = 1 + max(|d0|,|d1|), which is valid for a monic
//     polynomial.

import { BSRational } from './BSRational.js';
import { fma } from './Functions.js';
import { PolynomialRoot, polynomialRootBisect } from './PolynomialRoot.js';
import {
    PolynomialRootRational, type RootsScalarOps, rootsNumberOps, rootsRationalOps,
    solveLinearGeneric, solveMonicLinearGeneric, sortRoots
} from './RootsLinear.js';
import {
    computeDepressedQuadraticRoots, solveMonicQuadraticGeneric, solveQuadraticGeneric
} from './RootsQuadratic.js';

const HALF = BSRational.fromNumber(0.5);
const ONE_DIV_3 = BSRational.fromNumber(1, 3);

function rat(x: number): BSRational {
    return BSRational.fromNumber(x);
}

// The port of ComputeClassifiers(rG0, rG1, rG2, rG3, rD0, rD1, rM2Div3).
function computeClassifiersGeneral(rG0: BSRational, rG1: BSRational, rG2: BSRational,
    rG3: BSRational): { rD0: BSRational, rD1: BSRational, rM2Div3: BSRational } {
    return computeClassifiersMonic(rG0.div(rG3), rG1.div(rG3), rG2.div(rG3));
}

// The port of ComputeClassifiers(rM0, rM1, rM2, rD0, rD1, rM2Div3).
function computeClassifiersMonic(rM0: BSRational, rM1: BSRational, rM2: BSRational):
    { rD0: BSRational, rD1: BSRational, rM2Div3: BSRational } {
    const rM2Div3 = ONE_DIV_3.mul(rM2);
    const rD0 = rM0.sub(rM2Div3.mul(rM1.sub(rat(2).mul(rM2Div3).mul(rM2Div3))));
    const rD1 = rM1.sub(rM2.mul(rM2Div3));
    return { rD0: rD0, rD1: rD1, rM2Div3: rM2Div3 };
}

// The discriminant delta = -27*d0^2 - 4*d1^3 of the depressed cubic.
function computeDelta(rD0: BSRational, rD1: BSRational): BSRational {
    return rat(-27).mul(rD0).mul(rD0).add(rat(-4).mul(rD1).mul(rD1).mul(rD1));
}

// The midpoint of the bisection interval, as a rational number.
function midpoint(xMin: number, xMax: number): BSRational {
    return HALF.mul(rat(xMin).add(rat(xMax)));
}

// The port of ComputeDepressedRootsBisection.
function computeDepressedRootsBisection(rD0: BSRational, rD1: BSRational):
    PolynomialRootRational[] {
    const signD0 = rD0.getSign();
    const signD1 = rD1.getSign();
    if (signD0 === 0) {
        if (signD1 > 0) {
            // One real root, multiplicity 1.
            return [new PolynomialRootRational(rat(0), 1)];
        } else if (signD1 < 0) {
            // Three real roots, each multiplicity 1.
            const roots = computeDepressedQuadraticRoots(true, rD1);
            roots.push(new PolynomialRootRational(rat(0), 1));
            sortRoots(rootsRationalOps, roots);
            return roots;
        } else {
            // One real root, multiplicity 3.
            return [new PolynomialRootRational(rat(0), 3)];
        }
    }

    if (signD1 === 0) {
        // One real root, multiplicity 1. The root of x^3 + d0 has magnitude
        // |d0|^(1/3), so b = max(1,|d0|) bounds it.
        const d0 = rD0.toNumber();
        const b = Math.max(1.0, Math.abs(d0));
        const F = (x: number) => fma(x, x * x, d0);

        // Bisect on the interval [-b,b].
        const bisected = polynomialRootBisect(F, -1, +1, -b, b);
        return [new PolynomialRootRational(midpoint(bisected.xMin, bisected.xMax), 1)];
    }

    const rDelta = computeDelta(rD0, rD1);
    const signDelta = rDelta.getSign();
    if (signDelta > 0) {
        // Three real roots, each multiplicity 1. The derivative of
        // F(x) = x^3 + d1 * x + d0 is F'(x) = 3 * x^2 + d1 and must have two
        // real roots x0 and x1, which means d1 < 0. Let s = sqrt(-d1 / 3).
        // The F'(x) roots are x0 = -s and x1 = s. Using Samuelson's
        // inequality, an interval bounding the roots is [-2 * s, 2 * s].
        // Partition the interval into [-2 * s, -s], [-s, s], and [s, 2 * s].
        // Use bisection on each interval to estimate the roots of F(x).
        const rQRoots = computeDepressedQuadraticRoots(true, ONE_DIV_3.mul(rD1));
        const rS = rQRoots[1].x;
        const rTwoS = rat(2).mul(rS);
        const d0 = rD0.toNumber();
        const d1 = rD1.toNumber();
        const s = rS.toNumber();
        const twoS = rTwoS.toNumber();
        const F = (x: number) => fma(x, fma(x, x, d1), d0);

        // Bisect on the interval [-2 * s, -s].
        const b0 = polynomialRootBisect(F, -1, +1, -twoS, -s);
        // Bisect on the interval [-s, s].
        const b1 = polynomialRootBisect(F, +1, -1, -s, s);
        // Bisect on the interval [s, 2 * s].
        const b2 = polynomialRootBisect(F, -1, +1, s, twoS);
        return [
            new PolynomialRootRational(midpoint(b0.xMin, b0.xMax), 1),
            new PolynomialRootRational(midpoint(b1.xMin, b1.xMax), 1),
            new PolynomialRootRational(midpoint(b2.xMin, b2.xMax), 1)
        ];
    } else if (signDelta < 0) {
        // One real root, multiplicity 1. Cauchy's bound for the monic F(x) is
        // b = 1 + max{|d0|,|d1|}. Use bisection on the interval [-b,b] to
        // estimate the root. (Upstream uses max{1,|d0|,|d1|}, which is not a
        // valid bound; see the port notes at the top of the file.)
        const d0 = rD0.toNumber();
        const d1 = rD1.toNumber();
        const b = 1.0 + Math.max(Math.abs(d0), Math.abs(d1));
        const F = (x: number) => fma(x, fma(x, x, d1), d0);

        // Bisect on the interval [-b,b].
        const bisected = polynomialRootBisect(F, -1, +1, -b, b);
        return [new PolynomialRootRational(midpoint(bisected.xMin, bisected.xMax), 1)];
    } else {
        // One real root, multiplicity 1. One real root, multiplicity 2. The
        // roots are rational numbers, so F(x) = 0 exactly for each root x.
        return computeDeltaZeroRoots(rD0, rD1);
    }
}

// The delta = 0 branch shared by the bisection and closed-form solvers.
function computeDeltaZeroRoots(rD0: BSRational, rD1: BSRational): PolynomialRootRational[] {
    const rX0 = BSRational.fromNumber(-3, 2).mul(rD0).div(rD1);
    const rX1 = rat(-2).mul(rX0);
    if (rX0.lessThan(rX1)) {
        return [
            new PolynomialRootRational(rX0, 2),
            new PolynomialRootRational(rX1, 1)
        ];
    } else {
        return [
            new PolynomialRootRational(rX1, 1),
            new PolynomialRootRational(rX0, 2)
        ];
    }
}

// The port of ComputeDepressedRootsClosedForm.
function computeDepressedRootsClosedForm(rD0: BSRational, rD1: BSRational):
    PolynomialRootRational[] {
    const signD0 = rD0.getSign();
    const signD1 = rD1.getSign();
    if (signD0 === 0) {
        if (signD1 > 0) {
            // One real root, multiplicity 1.
            return [new PolynomialRootRational(rat(0), 1)];
        } else if (signD1 < 0) {
            // Three real roots, each multiplicity 1.
            const rSqrtNegD1 = BSRational.sqrt(rD1.negated());
            return [
                new PolynomialRootRational(rSqrtNegD1.negated(), 1),
                new PolynomialRootRational(rat(0), 1),
                new PolynomialRootRational(rSqrtNegD1, 1)
            ];
        } else {
            // One real root, multiplicity 3.
            return [new PolynomialRootRational(rat(0), 3)];
        }
    }

    if (signD1 === 0) {
        // One real root, multiplicity 1.
        if (signD0 > 0) {
            return [new PolynomialRootRational(
                BSRational.pow(rD0, ONE_DIV_3).negated(), 1)];
        } else {
            return [new PolynomialRootRational(
                BSRational.pow(rD0.negated(), ONE_DIV_3), 1)];
        }
    }

    const rDelta = computeDelta(rD0, rD1);
    const signDelta = rDelta.getSign();
    if (signDelta > 0) {
        // Three real roots, each multiplicity 1.
        const rSqrt3 = BSRational.sqrt(rat(3));
        const r3Div2 = BSRational.fromNumber(3, 2);
        const rD1Div3 = rD1.mul(ONE_DIV_3);
        const rRho = BSRational.pow(BSRational.fabs(rD1Div3), r3Div2);
        const rCbrtRho = BSRational.pow(rRho, ONE_DIV_3);
        const rTheta = rat(Math.atan2(
            BSRational.sqrt(rDelta.div(rat(27))).toNumber(), rD0.negated().toNumber()));
        const rThetaDiv3 = rTheta.mul(ONE_DIV_3);
        const rCosThetaDiv3 = BSRational.cos(rThetaDiv3);
        const rSinThetaDiv3 = BSRational.sin(rThetaDiv3);
        const rTemp0 = rCbrtRho.mul(rCosThetaDiv3);
        const rTemp1 = rSqrt3.mul(rCbrtRho).mul(rSinThetaDiv3);
        const r0 = rat(2).mul(rTemp0);
        const r1 = rTemp0.negated().sub(rTemp1);
        const r2 = rTemp0.negated().add(rTemp1);
        if (rSinThetaDiv3.getSign() > 0) {
            return [
                new PolynomialRootRational(r1, 1),
                new PolynomialRootRational(r2, 1),
                new PolynomialRootRational(r0, 1)
            ];
        } else {
            return [
                new PolynomialRootRational(r2, 1),
                new PolynomialRootRational(r1, 1),
                new PolynomialRootRational(r0, 1)
            ];
        }
    } else if (signDelta < 0) {
        // One real root, multiplicity 1.
        const rSqrtNegDeltaDiv27 = BSRational.sqrt(rDelta.negated().div(rat(27)));
        const rD1Div3 = rD1.mul(ONE_DIV_3);
        if (signD0 < 0) {
            const rW = HALF.mul(rD0.negated().add(rSqrtNegDeltaDiv27));
            const rCbrtW = BSRational.pow(rW, ONE_DIV_3);
            const r0 = rCbrtW.sub(rD1Div3.div(rCbrtW));
            return [new PolynomialRootRational(r0, 1)];
        } else {
            const rNegY = HALF.mul(rD0.add(rSqrtNegDeltaDiv27));
            const rCbrtY = BSRational.pow(rNegY, ONE_DIV_3).negated();
            const r0 = rCbrtY.sub(rD1Div3.div(rCbrtY));
            return [new PolynomialRootRational(r0, 1)];
        }
    } else {
        // One real root, multiplicity 1. One real root, multiplicity 2. The
        // roots are rational numbers, so F(x) = 0 exactly for each root x.
        return computeDeltaZeroRoots(rD0, rD1);
    }
}

// The port of RootsCubic<T>::ComputeDepressedRoots, which does not depend
// on T.
export function computeDepressedCubicRoots(useBisection: boolean, rD0: BSRational,
    rD1: BSRational): PolynomialRootRational[] {
    if (useBisection) {
        return computeDepressedRootsBisection(rD0, rD1);
    } else {
        return computeDepressedRootsClosedForm(rD0, rD1);
    }
}

// The port of HasZeroValuedRoots(useBisection, g1, g2, g3, roots).
function hasZeroValuedRootsGeneral<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, useBisection: boolean, g1: T, g2: T, g3: T): R[] {
    if (ops.isZero(g1)) {
        if (ops.isZero(g2)) {
            return [ops.makeRoot(ops.zero(), 3)];
        } else {
            const roots = solveLinearGeneric(ops, g2, g3);
            roots.push(ops.makeRoot(ops.zero(), 2));
            sortRoots(ops, roots);
            return roots;
        }
    } else {
        const roots = solveQuadraticGeneric(ops, useBisection, g1, g2, g3);
        roots.push(ops.makeRoot(ops.zero(), 1));
        sortRoots(ops, roots);
        return roots;
    }
}

// The port of HasZeroValuedRoots(useBisection, m1, m2, roots).
function hasZeroValuedRootsMonic<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, useBisection: boolean, m1: T, m2: T): R[] {
    if (ops.isZero(m1)) {
        if (ops.isZero(m2)) {
            return [ops.makeRoot(ops.zero(), 3)];
        } else {
            const roots = solveMonicLinearGeneric(ops, m2);
            roots.push(ops.makeRoot(ops.zero(), 2));
            sortRoots(ops, roots);
            return roots;
        }
    } else {
        const roots = solveMonicQuadraticGeneric(ops, useBisection, m1, m2);
        roots.push(ops.makeRoot(ops.zero(), 1));
        sortRoots(ops, roots);
        return roots;
    }
}

// The port of RootsCubic<T>::Solve(useBisection, g0, g1, g2, g3, roots).
export function solveCubicGeneric<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, useBisection: boolean, g0: T, g1: T, g2: T, g3: T): R[] {
    // Test whether the degree is smaller than 3.
    if (ops.isZero(g3)) {
        return solveQuadraticGeneric(ops, useBisection, g0, g1, g2);
    }

    // Test for zero-valued roots.
    if (ops.isZero(g0)) {
        return hasZeroValuedRootsGeneral(ops, useBisection, g1, g2, g3);
    }

    // At this time g0 and g3 are not zero. Transform the general cubic to a
    // depressed cubic, solve for its roots, and inverse transform them to
    // roots of the general cubic.
    const c = computeClassifiersGeneral(ops.toRational(g0), ops.toRational(g1),
        ops.toRational(g2), ops.toRational(g3));
    const rRoots = computeDepressedCubicRoots(useBisection, c.rD0, c.rD1);
    return rRoots.map(r => ops.makeRoot(ops.fromRational(r.x.sub(c.rM2Div3)), r.m));
}

// The port of RootsCubic<T>::Solve(useBisection, m0, m1, m2, roots).
export function solveMonicCubicGeneric<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, useBisection: boolean, m0: T, m1: T, m2: T): R[] {
    // Test for zero-valued roots.
    if (ops.isZero(m0)) {
        return hasZeroValuedRootsMonic(ops, useBisection, m1, m2);
    }

    // At this time m0 is not zero. Transform the monic cubic to a depressed
    // cubic, solve for its roots, and inverse transform them to roots of the
    // monic cubic.
    const c = computeClassifiersMonic(ops.toRational(m0), ops.toRational(m1),
        ops.toRational(m2));
    const rRoots = computeDepressedCubicRoots(useBisection, c.rD0, c.rD1);
    return rRoots.map(r => ops.makeRoot(ops.fromRational(r.x.sub(c.rM2Div3)), r.m));
}

// The port of RootsCubic<T>::Solve(useBisection, d0, d1, roots).
export function solveDepressedCubicGeneric<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, useBisection: boolean, d0: T, d1: T): R[] {
    // The cubic is already depressed, so no transforming and inverse
    // transforming are necessary.
    const rRoots = computeDepressedCubicRoots(useBisection, ops.toRational(d0),
        ops.toRational(d1));
    return rRoots.map(r => ops.makeRoot(ops.fromRational(r.x), r.m));
}

export class RootsCubic {
    // Solve the general cubic g0 + g1*x + g2*x^2 + g3*x^3 = 0.
    static solve(useBisection: boolean, g0: number, g1: number, g2: number,
        g3: number): PolynomialRoot[] {
        return solveCubicGeneric(rootsNumberOps, useBisection, g0, g1, g2, g3);
    }

    // Solve the monic cubic m0 + m1*x + m2*x^2 + x^3 = 0.
    static solveMonic(useBisection: boolean, m0: number, m1: number,
        m2: number): PolynomialRoot[] {
        return solveMonicCubicGeneric(rootsNumberOps, useBisection, m0, m1, m2);
    }

    // Solve the depressed cubic d0 + d1*x + x^3 = 0.
    static solveDepressed(useBisection: boolean, d0: number, d1: number): PolynomialRoot[] {
        return solveDepressedCubicGeneric(rootsNumberOps, useBisection, d0, d1);
    }

    // Solve the general cubic with rational coefficients.
    static solveRational(useBisection: boolean, g0: BSRational, g1: BSRational,
        g2: BSRational, g3: BSRational): PolynomialRootRational[] {
        return solveCubicGeneric(rootsRationalOps, useBisection, g0, g1, g2, g3);
    }

    // Solve the monic cubic with rational coefficients.
    static solveMonicRational(useBisection: boolean, m0: BSRational, m1: BSRational,
        m2: BSRational): PolynomialRootRational[] {
        return solveMonicCubicGeneric(rootsRationalOps, useBisection, m0, m1, m2);
    }

    // Solve the depressed cubic with rational coefficients.
    static solveDepressedRational(useBisection: boolean, d0: BSRational,
        d1: BSRational): PolynomialRootRational[] {
        return solveDepressedCubicGeneric(rootsRationalOps, useBisection, d0, d1);
    }

    // The upstream public ComputeDepressedRoots, which does not depend on T.
    static computeDepressedRoots(useBisection: boolean, rD0: BSRational,
        rD1: BSRational): PolynomialRootRational[] {
        return computeDepressedCubicRoots(useBisection, rD0, rD1);
    }
}
