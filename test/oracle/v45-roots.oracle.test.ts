// Replays oracle/cpp/cases/v45-roots.cpp. Keep the two files in the same
// order.
//
// RootsCubic.h and RootsQuartic.h classify roots with exact rational
// arithmetic. The "bisection" entry points (useBisection = true) use only
// +, *, /, comparisons and fma, so they are compared bit for bit; the
// "closedForm" ones go through std::pow / std::atan2 / std::cos / std::sin,
// where MSVC and V8 differ by an ulp, so they carry a tolerance.
import { describe } from 'vitest';
import { BSRational } from '../../src/BSRational.js';
import { PolynomialRoot } from '../../src/PolynomialRoot.js';
import { RootsCubic } from '../../src/RootsCubic.js';
import { RootsQuartic } from '../../src/RootsQuartic.js';
import { PolynomialRootRational } from '../../src/RootsLinear.js';
import { OracleFamily, type OracleIO } from './harness.js';

function emitRoots(io: OracleIO, roots: PolynomialRoot[]): void {
    io.outInt(roots.length);
    for (const r of roots) {
        io.outReal(r.x);
        io.outInt(r.m);
    }
}

function emitRationalRoots(io: OracleIO, roots: PolynomialRootRational[]): void {
    io.outInt(roots.length);
    for (const r of roots) {
        io.outReal(r.x.toNumber());
        io.outInt(r.m);
    }
}

// The tolerance of the closed-form cases. The depressed cubic's closed form
// evaluates pow(rho, 1/3) * cos(atan2(...)/3) and the quartic feeds that root
// through alpha = sqrt(2*T - d2) and beta = sqrt(T^2 - d0), so a 1 ulp libm
// difference is amplified by the conditioning of those square roots.
const CLOSED_FORM_TOL = 1e-9;

describe('oracle: v45-roots', () => {
    const family = new OracleFamily('v45-roots');

    family.case('RootsCubic.solve.bisection', (io) => {
        const c = io.reals(4);
        emitRoots(io, RootsCubic.solve(true, c[0], c[1], c[2], c[3]));
    }, { exact: true });

    family.case('RootsCubic.solve.closedForm', (io) => {
        const c = io.reals(4);
        emitRoots(io, RootsCubic.solve(false, c[0], c[1], c[2], c[3]));
    }, { tol: CLOSED_FORM_TOL });

    family.case('RootsCubic.solveMonic.bisection', (io) => {
        const c = io.reals(3);
        emitRoots(io, RootsCubic.solveMonic(true, c[0], c[1], c[2]));
    }, { exact: true });

    family.case('RootsCubic.solveMonic.closedForm', (io) => {
        const c = io.reals(3);
        emitRoots(io, RootsCubic.solveMonic(false, c[0], c[1], c[2]));
    }, { tol: CLOSED_FORM_TOL });

    family.case('RootsCubic.solveDepressed.bisection', (io) => {
        const c = io.reals(2);
        emitRoots(io, RootsCubic.solveDepressed(true, c[0], c[1]));
    }, { exact: true });

    family.case('RootsCubic.solveDepressed.closedForm', (io) => {
        const c = io.reals(2);
        emitRoots(io, RootsCubic.solveDepressed(false, c[0], c[1]));
    }, { tol: CLOSED_FORM_TOL });

    family.case('RootsCubic.computeDepressedRoots.bisection', (io) => {
        const c = io.reals(2);
        emitRationalRoots(io, RootsCubic.computeDepressedRoots(true,
            BSRational.fromNumber(c[0]), BSRational.fromNumber(c[1])));
    }, { exact: true });

    family.case('RootsCubic.computeDepressedRoots.closedForm', (io) => {
        const c = io.reals(2);
        emitRationalRoots(io, RootsCubic.computeDepressedRoots(false,
            BSRational.fromNumber(c[0]), BSRational.fromNumber(c[1])));
    }, { tol: CLOSED_FORM_TOL });

    family.case('RootsCubic.solveRational.bisection', (io) => {
        const c = io.reals(4);
        emitRationalRoots(io, RootsCubic.solveRational(true,
            BSRational.fromNumber(c[0]), BSRational.fromNumber(c[1]),
            BSRational.fromNumber(c[2]), BSRational.fromNumber(c[3])));
    }, { exact: true });

    family.case('RootsCubic.solveMonicRational.closedForm', (io) => {
        const c = io.reals(3);
        emitRationalRoots(io, RootsCubic.solveMonicRational(false,
            BSRational.fromNumber(c[0]), BSRational.fromNumber(c[1]),
            BSRational.fromNumber(c[2])));
    }, { tol: CLOSED_FORM_TOL });

    family.case('RootsCubic.solveDepressedRational.bisection', (io) => {
        const c = io.reals(2);
        emitRationalRoots(io, RootsCubic.solveDepressedRational(true,
            BSRational.fromNumber(c[0]), BSRational.fromNumber(c[1])));
    }, { exact: true });

    // The ill-conditioned closed-form regime the cases above reject. Only the
    // root count and the multiplicities are emitted; both come from exact
    // rational arithmetic and are comparable even there.
    family.case('RootsCubic.solveDepressed.closedForm.illConditioned', (io) => {
        const c = io.reals(2);
        const roots = RootsCubic.solveDepressed(false, c[0], c[1]);
        io.outInt(roots.length);
        for (const r of roots) { io.outInt(r.m); }
    }, { exact: true });

    // Deliberate deviation: upstream's 'signDelta < 0' bisection branch uses
    // the invalid root bound max(1,|d0|,|d1|) and returns an interval endpoint
    // instead of the root. See docs/UPSTREAM-FINDINGS.md, RootsCubic.h item 1.
    family.case('RootsCubic.solveDepressed.deviation.bisectionBound', (io) => {
        const c = io.reals(2);
        emitRoots(io, RootsCubic.solveDepressed(true, c[0], c[1]));
    }, { exact: true, deviation: 'issue #340, RootsCubic.h item 1' });

    family.case('RootsQuartic.solve.bisection', (io) => {
        const c = io.reals(5);
        emitRoots(io, RootsQuartic.solve(true, c[0], c[1], c[2], c[3], c[4]));
    }, { exact: true });

    family.case('RootsQuartic.solveMonic.bisection', (io) => {
        const c = io.reals(4);
        emitRoots(io, RootsQuartic.solveMonic(true, c[0], c[1], c[2], c[3]));
    }, { exact: true });

    family.case('RootsQuartic.solveDepressed.bisection', (io) => {
        const c = io.reals(3);
        emitRoots(io, RootsQuartic.solveDepressed(true, c[0], c[1], c[2]));
    }, { exact: true });

    family.case('RootsQuartic.computeDepressedRoots.bisection', (io) => {
        const c = io.reals(3);
        emitRationalRoots(io, RootsQuartic.computeDepressedRoots(true,
            BSRational.fromNumber(c[0]), BSRational.fromNumber(c[1]),
            BSRational.fromNumber(c[2])));
    }, { exact: true });

    family.case('RootsQuartic.solveRational.bisection', (io) => {
        const c = io.reals(5);
        emitRationalRoots(io, RootsQuartic.solveRational(true,
            BSRational.fromNumber(c[0]), BSRational.fromNumber(c[1]),
            BSRational.fromNumber(c[2]), BSRational.fromNumber(c[3]),
            BSRational.fromNumber(c[4])));
    }, { exact: true });

    family.case('RootsQuartic.solveDepressedRational.bisection', (io) => {
        const c = io.reals(3);
        emitRationalRoots(io, RootsQuartic.solveDepressedRational(true,
            BSRational.fromNumber(c[0]), BSRational.fromNumber(c[1]),
            BSRational.fromNumber(c[2])));
    }, { exact: true });

    family.case('RootsQuartic.solveMonicRational.closedForm', (io) => {
        const c = io.reals(4);
        emitRationalRoots(io, RootsQuartic.solveMonicRational(false,
            BSRational.fromNumber(c[0]), BSRational.fromNumber(c[1]),
            BSRational.fromNumber(c[2]), BSRational.fromNumber(c[3])));
    }, { tol: CLOSED_FORM_TOL });

    family.case('RootsQuartic.solve.closedForm', (io) => {
        const c = io.reals(5);
        emitRoots(io, RootsQuartic.solve(false, c[0], c[1], c[2], c[3], c[4]));
    }, { tol: CLOSED_FORM_TOL });

    family.case('RootsQuartic.solveMonic.closedForm', (io) => {
        const c = io.reals(4);
        emitRoots(io, RootsQuartic.solveMonic(false, c[0], c[1], c[2], c[3]));
    }, { tol: CLOSED_FORM_TOL });

    family.case('RootsQuartic.solveDepressed.closedForm', (io) => {
        const c = io.reals(3);
        emitRoots(io, RootsQuartic.solveDepressed(false, c[0], c[1], c[2]));
    }, { tol: CLOSED_FORM_TOL });

    family.case('RootsQuartic.computeDepressedRoots.closedForm', (io) => {
        const c = io.reals(3);
        emitRationalRoots(io, RootsQuartic.computeDepressedRoots(false,
            BSRational.fromNumber(c[0]), BSRational.fromNumber(c[1]),
            BSRational.fromNumber(c[2])));
    }, { tol: CLOSED_FORM_TOL });

    // d1 is exactly zero, so SolveBiquadratic runs: two nested square roots
    // and no libm call at all, hence exact in both variants.
    family.case('RootsQuartic.solveBiquadratic.bisection', (io) => {
        const c = io.reals(3);
        emitRoots(io, RootsQuartic.solveDepressed(true, c[0], c[1], c[2]));
    }, { exact: true });

    family.case('RootsQuartic.solveBiquadratic.closedForm', (io) => {
        const c = io.reals(3);
        emitRoots(io, RootsQuartic.solveDepressed(false, c[0], c[1], c[2]));
    }, { exact: true });

    // Deliberate deviation: upstream's positive-discriminant early-out tests
    // only d2 > 0 and reports four roots for a quartic with two
    // complex-conjugate pairs. docs/UPSTREAM-FINDINGS.md, RootsQuartic.h
    // item 2.
    family.case('RootsQuartic.solveDepressed.deviation.complexPairs', (io) => {
        const c = io.reals(3);
        emitRoots(io, RootsQuartic.solveDepressed(true, c[0], c[1], c[2]));
    }, { exact: true, deviation: 'issue #340, RootsQuartic.h item 2' });

    // Deliberate deviation: upstream reads index 1 of a reused root array that
    // the call did not write. docs/UPSTREAM-FINDINGS.md, RootsQuartic.h item 3.
    // The case runs the closed form, so it is compared at 1e-6 rather than
    // bit-for-bit: that keeps the libm disagreement of the resolvent cubic
    // (measured maximum 6.7e-9 over 2000 records) from being counted as a
    // deviation, while a stale read substitutes alpha = sqrt(2*T - d2), of
    // order 1, for a beta that is of order 1e-200, which moves the roots by
    // order 1.
    family.case('RootsQuartic.solveDepressed.deviation.staleSqrt', (io) => {
        const c = io.reals(3);
        emitRoots(io, RootsQuartic.solveDepressed(false, c[0], c[1], c[2]));
    }, { tol: 1e-6, deviation: 'issue #340, RootsQuartic.h item 3' });

    family.finish();
});
