// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RootsQuartic.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the real-valued roots of a quartic polynomial with real-valued
// coefficients. The general quartic polynomial is
//   g(x) = g0 + g1 * x + g2 * x^2 + g3 * x^3 + g4 * x^4
// where g4 is not zero. The monic quartic polynomial is
//   m(x) = m0 + m1 * x + m2 * x^2 + m3 * x^3 + x^4
// The depressed quartic polynomial is
//   d(x) = d0 + d1 * x + d2 * x^2 + x^4
// The classification of roots and multiplicities is performed using rational
// arithmetic for exactness. For algorithmic details, see
// https://www.geometrictools.com/Documentation/LowDegreePolynomialRoots.pdf
//
// The code uses bounding intervals for roots. For a polynomial of degree n
// with all real roots, Samuelson's inequality
// https://en.wikipedia.org/wiki/Samuelson%27s_inequality
// provides an interval [b0,b1] where
//   b0 = (-p[n-1] - (n-1) * s) / (n * p[n])
//   b1 = (-p[n-1] + (n-1) * s) / (n * p[n])
//   s = sqrt(p[n-1]^2 - 2 * n * p[n] * p[n-2] / (n-1))
// Applied to the general quartic (n = 4),
//   b0 = (-p3 - 3 * s) / (4 * p4)
//   b1 = (-p3 + 3 * s) / (4 * p4)
//   s = sqrt(p3^2 - 8 * p4 * p2 / 3)
// Applied to the depressed quartic when it has all real roots,
//   b0 = -sqrt(-3 * p2 / 2)
//   b1 = +sqrt(-3 * p2 / 2)
//
// Port notes:
//   - See RootsLinear.ts for the RootsScalarOps<T, R> mechanism that stands
//     in for the C++ template parameter T.
//   - The upstream out-parameter 'PolynomialRoot<T>* roots' plus the size_t
//     count become a returned array of roots.
//   - Upstream ComputeDepressedRootsBisection and ComputeDepressedRootsClosedForm
//     are character-for-character identical except that the former passes
//     'true' and the latter 'false' to the cubic and quadratic sub-solvers.
//     The port keeps a single function parameterized by useBisection, which
//     is behavior-identical.
//   - Upstream bug fixed (see the PR body): the 'two complex-conjugate pairs'
//     early-out for a positive discriminant tests only d2 > 0. The full
//     criterion for four real roots of the depressed quartic
//     x^4 + d2*x^2 + d1*x + d0 with positive discriminant is P = 8*d2 < 0 and
//     D = 64*d0 - 16*d2^2 < 0. For example
//     2.3588208672590554 + 1.996267369017005*x - 0.8864723690785468*x^2
//     - 2.481981527991593*x^3 + 1.4308462475892156*x^4 has a positive
//     discriminant, d2 < 0 and D > 0, so upstream falls through and reports
//     four roots for a quartic that has none.
//   - Upstream repeatedly reuses one 'rQRoots' array to extract square roots
//     via RootsQuadratic<Rational>::ComputeDepressedRoots, reading index 1
//     even when the call writes fewer than two entries (a stale read). The
//     port routes every such extraction through rationalSqrtViaQuadratic (see
//     RootsQuadratic.ts), which returns 0 for a zero or negative argument.

import { BSRational } from './BSRational.js';
import { PolynomialRoot } from './PolynomialRoot.js';
import {
    PolynomialRootRational, type RootsScalarOps, rootsNumberOps, rootsRationalOps,
    solveLinearGeneric, solveMonicLinearGeneric, sortRoots
} from './RootsLinear.js';
import {
    rationalSqrtViaQuadratic, solveMonicQuadraticGeneric, solveQuadraticGeneric
} from './RootsQuadratic.js';
import {
    computeDepressedCubicRoots, solveCubicGeneric, solveMonicCubicGeneric
} from './RootsCubic.js';

const HALF = BSRational.fromNumber(0.5);
const ONE_DIV_3 = BSRational.fromNumber(1, 3);
const ONE_DIV_4 = BSRational.fromNumber(1, 4);

function rat(x: number): BSRational {
    return BSRational.fromNumber(x);
}

// The port of ComputeClassifiers(rG0, rG1, rG2, rG3, rG4, rD0, rD1, rD2, rM3Div4).
function computeClassifiersGeneral(rG0: BSRational, rG1: BSRational, rG2: BSRational,
    rG3: BSRational, rG4: BSRational):
    { rD0: BSRational, rD1: BSRational, rD2: BSRational, rM3Div4: BSRational } {
    return computeClassifiersMonic(rG0.div(rG4), rG1.div(rG4), rG2.div(rG4), rG3.div(rG4));
}

// The port of ComputeClassifiers(rM0, rM1, rM2, rM3, rD0, rD1, rD2, rM3Div4).
function computeClassifiersMonic(rM0: BSRational, rM1: BSRational, rM2: BSRational,
    rM3: BSRational):
    { rD0: BSRational, rD1: BSRational, rD2: BSRational, rM3Div4: BSRational } {
    const rM3Div4 = ONE_DIV_4.mul(rM3);
    const rM3Div4Sqr = rM3Div4.mul(rM3Div4);
    const rD0 = rM0.sub(rM3Div4.mul(rM1.sub(rM3Div4.mul(rM2.sub(rat(3).mul(rM3Div4Sqr))))));
    const rD1 = rM1.sub(rat(2).mul(rM3Div4).mul(rM2.sub(rat(4).mul(rM3Div4Sqr))));
    const rD2 = rM2.sub(rat(6).mul(rM3Div4Sqr));
    return { rD0: rD0, rD1: rD1, rD2: rD2, rM3Div4: rM3Div4 };
}

// The port of SolveBiquadratic, which solves d0 + d2*x^2 + x^4 = 0.
function solveBiquadratic(useBisection: boolean, rD0: BSRational, rD2: BSRational):
    PolynomialRootRational[] {
    const rS = rat(-0.5).mul(rD2);
    const rT = rS.mul(rS).sub(rD0);
    const signT = rT.getSign();
    if (signT > 0) {
        const rSqrtT = rationalSqrtViaQuadratic(useBisection, rT);
        const rSPsqrtT = rS.add(rSqrtT);
        const rSMsqrtT = rD0.div(rSPsqrtT);
        const signSPsqrtT = rSPsqrtT.getSign();
        const signSMsqrtT = rSMsqrtT.getSign();
        if (signSMsqrtT > 0) {
            // Four real roots.
            const r0 = rationalSqrtViaQuadratic(useBisection, rSMsqrtT);
            const r1 = rationalSqrtViaQuadratic(useBisection, rSPsqrtT);
            const roots = [
                new PolynomialRootRational(r0, 1),
                new PolynomialRootRational(r0.negated(), 1),
                new PolynomialRootRational(r1, 1),
                new PolynomialRootRational(r1.negated(), 1)
            ];
            sortRoots(rootsRationalOps, roots);
            return roots;
        } else if (signSPsqrtT < 0) {
            // Two complex-conjugate pairs.
            return [];
        } else {
            // signSMsqrtT < 0 and signSPsqrtT > 0: two real roots, one
            // complex-conjugate pair.
            const r0 = rationalSqrtViaQuadratic(useBisection, rSPsqrtT);
            if (r0.getSign() > 0) {
                return [
                    new PolynomialRootRational(r0.negated(), 1),
                    new PolynomialRootRational(r0, 1)
                ];
            } else {
                return [
                    new PolynomialRootRational(r0, 1),
                    new PolynomialRootRational(r0.negated(), 1)
                ];
            }
        }
    } else if (signT < 0) {
        // Two complex-conjugate pairs.
        return [];
    } else {
        if (rS.getSign() > 0) {
            // Two real roots, each of multiplicity 2.
            const r0 = rationalSqrtViaQuadratic(useBisection, rS);
            if (r0.getSign() > 0) {
                return [
                    new PolynomialRootRational(r0.negated(), 2),
                    new PolynomialRootRational(r0, 2)
                ];
            } else {
                return [
                    new PolynomialRootRational(r0, 2),
                    new PolynomialRootRational(r0.negated(), 2)
                ];
            }
        } else {
            // One complex-conjugate pair of multiplicity 2.
            return [];
        }
    }
}

// The port of RootsQuartic<T>::ComputeDepressedRoots, which does not depend
// on T. Upstream has two copies of this body, one passing 'true' to the
// sub-solvers (ComputeDepressedRootsBisection) and one passing 'false'
// (ComputeDepressedRootsClosedForm).
export function computeDepressedQuarticRoots(useBisection: boolean, rD0: BSRational,
    rD1: BSRational, rD2: BSRational): PolynomialRootRational[] {
    const signD0 = rD0.getSign();
    const signD1 = rD1.getSign();
    const signD2 = rD2.getSign();
    if (signD0 === 0) {
        if (signD1 === 0) {
            if (signD2 > 0) {
                // One real root, multiplicity 2.
                return [new PolynomialRootRational(rat(0), 2)];
            } else if (signD2 < 0) {
                // Three real roots, one with multiplicity 2, two with
                // multiplicity 1.
                const rSqrtNegD2 = rationalSqrtViaQuadratic(useBisection, rD2.negated());
                return [
                    new PolynomialRootRational(rSqrtNegD2.negated(), 1),
                    new PolynomialRootRational(rat(0), 2),
                    new PolynomialRootRational(rSqrtNegD2, 1)
                ];
            } else {
                // One real root, multiplicity 4.
                return [new PolynomialRootRational(rat(0), 4)];
            }
        } else {
            // Zero is a root of multiplicity 1. The cubic solver computes the
            // other roots.
            const roots = computeDepressedCubicRoots(useBisection, rD1, rD2);
            roots.push(new PolynomialRootRational(rat(0), 1));
            sortRoots(rootsRationalOps, roots);
            return roots;
        }
    }

    // At this time d0 != 0.
    if (signD1 === 0) {
        return solveBiquadratic(useBisection, rD0, rD2);
    }

    // At this time, d0 != 0 and d1 != 0.
    const rD0sqr = rD0.mul(rD0);
    const rD1sqr = rD1.mul(rD1);
    const rD2sqr = rD2.mul(rD2);
    const rDelta = rD1sqr.mul(rat(-27).mul(rD1sqr).add(
        rat(4).mul(rD2).mul(rat(36).mul(rD0).sub(rD2sqr)))).add(
        rat(16).mul(rD0).mul(rD2sqr.mul(rD2sqr.sub(rat(8).mul(rD0))).add(
            rat(16).mul(rD0sqr))));

    const signDelta = rDelta.getSign();
    if (signDelta === 0) {
        // Process the repeated roots.
        const rA0 = rat(12).mul(rD0).add(rD2sqr);
        if (rA0.getSign() === 0) {
            // Case (x-r0)^3 (x-r1), d2 < 0 guaranteed.
            const r0 = rat(-0.75).mul(rD1).div(rD2);
            const r1 = rat(-3).mul(r0);
            if (r0.lessThan(r1)) {
                return [
                    new PolynomialRootRational(r0, 3),
                    new PolynomialRootRational(r1, 1)
                ];
            } else {
                return [
                    new PolynomialRootRational(r1, 1),
                    new PolynomialRootRational(r0, 3)
                ];
            }
        }

        // Non-zero denominator guaranteed.
        const rA1 = rat(4).mul(rD0).sub(rD2sqr);
        const r0 = rD1.negated().mul(rA0).div(
            rat(9).mul(rD1sqr).sub(rat(2).mul(rD2).mul(rA1)));
        const rQDiscriminant = rD2.add(rat(2).mul(r0).mul(r0)).negated();
        if (rQDiscriminant.getSign() > 0) {
            // Case (x-r0)^2 (x-r1) (x-r2).
            const rSqrtQDiscriminant = rationalSqrtViaQuadratic(useBisection, rQDiscriminant);
            const r1 = r0.negated().sub(rSqrtQDiscriminant);
            const r2 = r0.negated().add(rSqrtQDiscriminant);
            const roots = [
                new PolynomialRootRational(r0, 2),
                new PolynomialRootRational(r1, 1),
                new PolynomialRootRational(r2, 1)
            ];
            sortRoots(rootsRationalOps, roots);
            return roots;
        }

        // Case (x-r0)^2 (x-z0) (x-z0c).
        return [new PolynomialRootRational(r0, 2)];
    }

    if (signDelta > 0) {
        // When delta > 0 the quartic has either four real roots or two
        // complex-conjugate pairs. It has four real roots if and only if
        // P = 8*d2 < 0 and D = 64*d0 - 16*d2^2 < 0; otherwise there are no
        // real roots. Upstream tests only d2 > 0, so it falls through to the
        // root extraction below for quartics with d2 <= 0 and D >= 0 and
        // reports four spurious roots (see the port notes). The factor 16 of
        // D is dropped here because only its sign matters.
        const rD = rat(4).mul(rD0).sub(rD2sqr);
        if (signD2 >= 0 || rD.getSign() >= 0) {
            // Two complex-conjugate pairs.
            return [];
        }
    }

    // Transform the discriminant (monic cubic) to a depressed cubic.
    const rM0 = rat(0.125).mul(rat(4).mul(rD0).mul(rD2).sub(rD1sqr));
    const rM1 = rD0.negated();
    const rM2 = rat(-0.5).mul(rD2);
    const rM2Div3 = ONE_DIV_3.mul(rM2);
    const rC0 = rM0.sub(rM2Div3.mul(rM1.sub(rat(2).mul(rM2Div3).mul(rM2Div3))));
    const rC1 = rM1.sub(rM2.mul(rM2Div3));

    // Compute the roots of the depressed cubic. The minimum root of the
    // depressed polynomial corresponds to the maximum root of the monic
    // polynomial. Also inverse-transform the root.
    const rCRoots = computeDepressedCubicRoots(useBisection, rC0, rC1);
    const rT = rCRoots[0].x.sub(rM2Div3);

    const rAlphaSqr = rat(2).mul(rT).sub(rD2);
    const rAlpha = rationalSqrtViaQuadratic(useBisection, rAlphaSqr);
    const rSignD1 = (signD1 > 0 ? rat(1) : rat(-1));
    const rArg = rT.mul(rT).sub(rD0);
    const rBeta = rSignD1.mul(rationalSqrtViaQuadratic(useBisection, rArg));
    const rDiscr0 = rAlphaSqr.sub(rat(4).mul(rT.add(rBeta)));
    const rSqrtDiscr0 = rationalSqrtViaQuadratic(useBisection, rDiscr0);
    const rDiscr1 = rAlphaSqr.sub(rat(4).mul(rT.sub(rBeta)));
    const rSqrtDiscr1 = rationalSqrtViaQuadratic(useBisection, rDiscr1);

    if (signDelta > 0) {
        // Case (x-r0)(x-r1)(x-r2)(x-r3).
        const r0 = HALF.mul(rAlpha.sub(rSqrtDiscr0));
        const r1 = HALF.mul(rAlpha.add(rSqrtDiscr0));
        const r2 = HALF.mul(rAlpha.negated().sub(rSqrtDiscr1));
        const r3 = HALF.mul(rAlpha.negated().add(rSqrtDiscr1));
        const roots = [
            new PolynomialRootRational(r0, 1),
            new PolynomialRootRational(r1, 1),
            new PolynomialRootRational(r2, 1),
            new PolynomialRootRational(r3, 1)
        ];
        sortRoots(rootsRationalOps, roots);
        return roots;
    } else {
        // signDelta < 0. Case (x-r0)(x-r1)(x-z0)(x-z0c).
        let r0: BSRational;
        let r1: BSRational;
        if (signD1 > 0) {
            r0 = HALF.mul(rAlpha.negated().sub(rSqrtDiscr1));
            r1 = HALF.mul(rAlpha.negated().add(rSqrtDiscr1));
        } else {
            r0 = HALF.mul(rAlpha.sub(rSqrtDiscr0));
            r1 = HALF.mul(rAlpha.add(rSqrtDiscr0));
        }
        return [
            new PolynomialRootRational(r0, 1),
            new PolynomialRootRational(r1, 1)
        ];
    }
}

// The port of HasZeroValuedRoots(useBisection, g1, g2, g3, g4, roots).
function hasZeroValuedRootsGeneral<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, useBisection: boolean, g1: T, g2: T, g3: T, g4: T): R[] {
    if (ops.isZero(g1)) {
        if (ops.isZero(g2)) {
            if (ops.isZero(g3)) {
                return [ops.makeRoot(ops.zero(), 4)];
            } else {
                const roots = solveLinearGeneric(ops, g3, g4);
                roots.push(ops.makeRoot(ops.zero(), 3));
                sortRoots(ops, roots);
                return roots;
            }
        } else {
            const roots = solveQuadraticGeneric(ops, useBisection, g2, g3, g4);
            roots.push(ops.makeRoot(ops.zero(), 2));
            sortRoots(ops, roots);
            return roots;
        }
    } else {
        const roots = solveCubicGeneric(ops, useBisection, g1, g2, g3, g4);
        roots.push(ops.makeRoot(ops.zero(), 1));
        sortRoots(ops, roots);
        return roots;
    }
}

// The port of HasZeroValuedRoots(useBisection, m1, m2, m3, roots).
function hasZeroValuedRootsMonic<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, useBisection: boolean, m1: T, m2: T, m3: T): R[] {
    if (ops.isZero(m1)) {
        if (ops.isZero(m2)) {
            if (ops.isZero(m3)) {
                return [ops.makeRoot(ops.zero(), 4)];
            } else {
                const roots = solveMonicLinearGeneric(ops, m3);
                roots.push(ops.makeRoot(ops.zero(), 3));
                sortRoots(ops, roots);
                return roots;
            }
        } else {
            const roots = solveMonicQuadraticGeneric(ops, useBisection, m2, m3);
            roots.push(ops.makeRoot(ops.zero(), 2));
            sortRoots(ops, roots);
            return roots;
        }
    } else {
        const roots = solveMonicCubicGeneric(ops, useBisection, m1, m2, m3);
        roots.push(ops.makeRoot(ops.zero(), 1));
        sortRoots(ops, roots);
        return roots;
    }
}

// The port of RootsQuartic<T>::Solve(useBisection, g0, g1, g2, g3, g4, roots).
export function solveQuarticGeneric<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, useBisection: boolean, g0: T, g1: T, g2: T, g3: T,
    g4: T): R[] {
    // Test whether the degree is smaller than 4.
    if (ops.isZero(g4)) {
        return solveCubicGeneric(ops, useBisection, g0, g1, g2, g3);
    }

    // Test for zero-valued roots.
    if (ops.isZero(g0)) {
        return hasZeroValuedRootsGeneral(ops, useBisection, g1, g2, g3, g4);
    }

    // At this time g0 and g4 are not zero. Transform the general quartic to a
    // depressed quartic, solve for its roots, and inverse transform them to
    // roots of the general quartic.
    const c = computeClassifiersGeneral(ops.toRational(g0), ops.toRational(g1),
        ops.toRational(g2), ops.toRational(g3), ops.toRational(g4));
    const rRoots = computeDepressedQuarticRoots(useBisection, c.rD0, c.rD1, c.rD2);
    return rRoots.map(r => ops.makeRoot(ops.fromRational(r.x.sub(c.rM3Div4)), r.m));
}

// The port of RootsQuartic<T>::Solve(useBisection, m0, m1, m2, m3, roots).
export function solveMonicQuarticGeneric<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, useBisection: boolean, m0: T, m1: T, m2: T, m3: T): R[] {
    // Test for zero-valued roots.
    if (ops.isZero(m0)) {
        return hasZeroValuedRootsMonic(ops, useBisection, m1, m2, m3);
    }

    // At this time m0 is not zero. Transform the monic quartic to a depressed
    // quartic, solve for its roots, and inverse transform them to roots of
    // the monic quartic.
    const c = computeClassifiersMonic(ops.toRational(m0), ops.toRational(m1),
        ops.toRational(m2), ops.toRational(m3));
    const rRoots = computeDepressedQuarticRoots(useBisection, c.rD0, c.rD1, c.rD2);
    return rRoots.map(r => ops.makeRoot(ops.fromRational(r.x.sub(c.rM3Div4)), r.m));
}

// The port of RootsQuartic<T>::Solve(useBisection, d0, d1, d2, roots).
export function solveDepressedQuarticGeneric<T, R extends { x: T, m: number }>(
    ops: RootsScalarOps<T, R>, useBisection: boolean, d0: T, d1: T, d2: T): R[] {
    // The quartic is already depressed, so no transforming and inverse
    // transforming are necessary.
    const rRoots = computeDepressedQuarticRoots(useBisection, ops.toRational(d0),
        ops.toRational(d1), ops.toRational(d2));
    return rRoots.map(r => ops.makeRoot(ops.fromRational(r.x), r.m));
}

export class RootsQuartic {
    // Solve the general quartic g0 + g1*x + g2*x^2 + g3*x^3 + g4*x^4 = 0.
    static solve(useBisection: boolean, g0: number, g1: number, g2: number, g3: number,
        g4: number): PolynomialRoot[] {
        return solveQuarticGeneric(rootsNumberOps, useBisection, g0, g1, g2, g3, g4);
    }

    // Solve the monic quartic m0 + m1*x + m2*x^2 + m3*x^3 + x^4 = 0.
    static solveMonic(useBisection: boolean, m0: number, m1: number, m2: number,
        m3: number): PolynomialRoot[] {
        return solveMonicQuarticGeneric(rootsNumberOps, useBisection, m0, m1, m2, m3);
    }

    // Solve the depressed quartic d0 + d1*x + d2*x^2 + x^4 = 0.
    static solveDepressed(useBisection: boolean, d0: number, d1: number,
        d2: number): PolynomialRoot[] {
        return solveDepressedQuarticGeneric(rootsNumberOps, useBisection, d0, d1, d2);
    }

    // Solve the general quartic with rational coefficients.
    static solveRational(useBisection: boolean, g0: BSRational, g1: BSRational,
        g2: BSRational, g3: BSRational, g4: BSRational): PolynomialRootRational[] {
        return solveQuarticGeneric(rootsRationalOps, useBisection, g0, g1, g2, g3, g4);
    }

    // Solve the monic quartic with rational coefficients.
    static solveMonicRational(useBisection: boolean, m0: BSRational, m1: BSRational,
        m2: BSRational, m3: BSRational): PolynomialRootRational[] {
        return solveMonicQuarticGeneric(rootsRationalOps, useBisection, m0, m1, m2, m3);
    }

    // Solve the depressed quartic with rational coefficients.
    static solveDepressedRational(useBisection: boolean, d0: BSRational, d1: BSRational,
        d2: BSRational): PolynomialRootRational[] {
        return solveDepressedQuarticGeneric(rootsRationalOps, useBisection, d0, d1, d2);
    }

    // The upstream public ComputeDepressedRoots, which does not depend on T.
    static computeDepressedRoots(useBisection: boolean, rD0: BSRational, rD1: BSRational,
        rD2: BSRational): PolynomialRootRational[] {
        return computeDepressedQuarticRoots(useBisection, rD0, rD1, rD2);
    }
}
