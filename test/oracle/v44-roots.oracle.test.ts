// Replays oracle/cpp/cases/v44-roots.cpp. Keep the two files in the same
// order.
//
// Everything in this family is arithmetic only (+ - * /, comparisons,
// Math.sqrt, Math.abs, stdMax) except RootsPolynomial's low-degree closed
// forms, which call Math.pow / Math.atan2 / Math.cos / Math.sin; those three
// cases carry a libm tolerance and everything else is compared bit for bit.
import { describe } from 'vitest';
import { BSRational } from '../../src/BSRational.js';
import { CubicRootsQR } from '../../src/CubicRootsQR.js';
import { Polynomial1, greatestCommonDivisor, squareFreeFactorization }
    from '../../src/Polynomial1.js';
import { PolynomialCurve } from '../../src/PolynomialCurve.js';
import { PolynomialRoot, polynomialRootBisect } from '../../src/PolynomialRoot.js';
import { QuarticRootsQR } from '../../src/QuarticRootsQR.js';
import { RootsBisection } from '../../src/RootsBisection.js';
import { RootsBisection1 } from '../../src/RootsBisection1.js';
import { RootsBisection2 } from '../../src/RootsBisection2.js';
import { RootsBrentsMethod } from '../../src/RootsBrentsMethod.js';
import { RootsGeneralPolynomial } from '../../src/RootsGeneralPolynomial.js';
import { RootsLinear, type PolynomialRootRational } from '../../src/RootsLinear.js';
import { RootsPolynomial, type RootMultiplicity } from '../../src/RootsPolynomial.js';
import { RootsQuadratic } from '../../src/RootsQuadratic.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// The C++ side records the degree first, then degree+1 coefficients, and
// builds the polynomial with the degree constructor (no leading-zero
// elimination).
function readPoly(io: OracleIO): Polynomial1 {
    const degree = io.integer();
    const p = new Polynomial1(degree);
    for (let i = 0; i <= degree; ++i) { p.set(i, io.real()); }
    return p;
}

function emitPoly(io: OracleIO, p: Polynomial1): void {
    const degree = p.getDegree();
    io.outInt(degree);
    for (let i = 0; i <= degree; ++i) { io.outReal(p.get(i)); }
}

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

function emitRootMap(io: OracleIO, rmMap: RootMultiplicity[]): void {
    io.outInt(rmMap.length);
    for (const rm of rmMap) {
        io.outReal(rm.root);
        io.outInt(rm.multiplicity);
    }
}

// Horner, written identically to HornerAt in the C++ case file.
function hornerAt(c: readonly number[], t: number): number {
    let i = c.length - 1;
    let result = c[i];
    while (i-- > 0) { result = t * result + c[i]; }
    return result;
}

// MSVC and V8 round pow/atan2/cos/sin differently in the last bit, and
// RootsPolynomial's closed forms build roots out of those values.
const LIBM_TOL = 1e-9;

describe('oracle: v44-roots', () => {
    const family = new OracleFamily('v44-roots');

    family.case('Polynomial1.arithmetic', (io) => {
        const p0 = readPoly(io);
        const p1 = readPoly(io);
        const s = io.real();
        const nonzero = io.real();
        emitPoly(io, p0.clone());
        emitPoly(io, p0.negate());
        emitPoly(io, p0.add(p1));
        emitPoly(io, p0.sub(p1));
        emitPoly(io, p0.mul(p1));
        emitPoly(io, p0.add(s));
        emitPoly(io, p0.add(s));
        emitPoly(io, p0.sub(s));
        emitPoly(io, p0.subFrom(s));
        emitPoly(io, p0.mul(s));
        emitPoly(io, p0.mul(s));
        emitPoly(io, p0.div(nonzero));
        io.outBool(p0.equals(p1));
        io.outBool(p0.notEquals(p1));
        io.outBool(p0.lessThan(p1));
        io.outBool(p0.lessThanOrEqual(p1));
        io.outBool(p0.greaterThan(p1));
        io.outBool(p0.greaterThanOrEqual(p1));
    }, { exact: true });

    family.case('Polynomial1.evaluateAndTransforms', (io) => {
        const p = readPoly(io);
        const t0 = io.real();
        const t1 = io.real();
        io.outReal(p.evaluate(t0));
        io.outReal(p.evaluate(t1));
        const d1 = p.getDerivative();
        emitPoly(io, d1);
        emitPoly(io, d1.getDerivative());
        emitPoly(io, p.getInversion());
        emitPoly(io, p.getTranslation(t0));
        emitPoly(io, p.getTranslation(t1));
    }, { exact: true });

    family.case('Polynomial1.mutators', (io) => {
        const p = readPoly(io);
        const newDegree = io.integer();
        const q = p.clone();
        q.eliminateLeadingZeros();
        emitPoly(io, q);
        const r = p.clone();
        r.makeMonic();
        emitPoly(io, r);
        const t = p.clone();
        t.setDegree(newDegree);
        emitPoly(io, t);
        const u = p.clone();
        u.setCoefficients(2.5);
        emitPoly(io, u);
    }, { exact: true });

    family.case('Polynomial1.divide', (io) => {
        const p = readPoly(io);
        const d = readPoly(io);
        const r = p.divide(d);
        emitPoly(io, r.quotient);
        emitPoly(io, r.remainder);
    }, { exact: true });

    family.case('Polynomial1.greatestCommonDivisor', (io) => {
        const p0 = readPoly(io);
        const p1 = readPoly(io);
        emitPoly(io, greatestCommonDivisor(p0, p1));
        emitPoly(io, greatestCommonDivisor(p1, p0));
    }, { exact: true });

    family.case('Polynomial1.squareFreeFactorization', (io) => {
        const degree = io.integer();
        const f = new Polynomial1(degree);
        for (let i = 0; i <= degree; ++i) { f.set(i, io.real()); }
        const factors = squareFreeFactorization(f);
        io.outInt(factors.length);
        for (const factor of factors) { emitPoly(io, factor); }
    }, { exact: true });

    family.case('PolynomialCurve.evaluate', (io) => {
        const tmin = io.real();
        const tmax = io.real();
        const components = [readPoly(io), readPoly(io), readPoly(io)];
        const t = io.real();
        const order = io.integer();
        const curve = new PolynomialCurve(3, tmin, tmax, components);
        for (let i = 0; i < 3; ++i) {
            emitPoly(io, curve.getPolynomial(i));
            emitPoly(io, curve.getDer1Polynomial(i));
            emitPoly(io, curve.getDer2Polynomial(i));
            emitPoly(io, curve.getDer3Polynomial(i));
        }
        const jet = [new Vector(3), new Vector(3), new Vector(3), new Vector(3)];
        curve.evaluate(t, order, jet);
        for (let k = 0; k <= order; ++k) { io.outVec(jet[k]); }
    }, { exact: true });

    // Deliberate deviation: upstream never sets mConstructed, so its
    // operator bool is false for every valid curve. See
    // docs/UPSTREAM-FINDINGS.md, PolynomialCurve.h item 3, issue #319.
    family.case('PolynomialCurve.deviation.constructed', (io) => {
        const tmin = io.real();
        const tmax = io.real();
        const components = [readPoly(io), readPoly(io), readPoly(io)];
        const curve = new PolynomialCurve(3, tmin, tmax, components);
        io.outBool(curve.isConstructed());
    }, { exact: true, deviation: 'issue #319, PolynomialCurve mConstructed' });

    function readBracketPoly(io: OracleIO, degree: number): number[] {
        return io.reals(degree + 1);
    }

    family.case('PolynomialRoot.bisect', (io) => {
        const c = readBracketPoly(io, 3);
        const xMin = io.real();
        const xMax = io.real();
        const signFMin = io.integer() === 0 ? -1 : +1;
        const signFMax = -signFMin;
        const r = polynomialRootBisect((x) => hornerAt(c, x), signFMin, signFMax,
            xMin, xMax);
        io.outReal(r.xMin);
        io.outReal(r.xMax);
    }, { exact: true });

    family.case('PolynomialRoot.compare', (io) => {
        const x0 = io.real();
        const x1 = io.real();
        const m0 = io.integer();
        const m1 = io.integer();
        const r0 = new PolynomialRoot(x0, m0);
        const r1 = new PolynomialRoot(x1, m1);
        io.outBool(r0.equals(r1));
        io.outBool(r0.lessThan(r1));
        io.outBool(r1.lessThan(r0));
        const def = new PolynomialRoot();
        io.outReal(def.x);
        io.outInt(def.m);
    }, { exact: true });

    family.case('RootsBisection.find', (io) => {
        const c = readBracketPoly(io, 4);
        const t0 = io.real();
        const t1 = io.real();
        const maxIterations = io.integer();
        const r = RootsBisection.find((t) => hornerAt(c, t), t0, t1, maxIterations);
        io.outInt(r.iterations);
        io.outReal(r.root);
    }, { exact: true });

    family.case('RootsBisection.findWithValues', (io) => {
        const c = readBracketPoly(io, 4);
        const t0 = io.real();
        const t1 = io.real();
        const maxIterations = io.integer();
        const useSigns = io.boolean();
        const F = (t: number) => hornerAt(c, t);
        let f0 = F(t0), f1 = F(t1);
        if (useSigns) {
            f0 = (f0 > 0 ? 1 : (f0 < 0 ? -1 : 0));
            f1 = (f1 > 0 ? 1 : (f1 < 0 ? -1 : 0));
        }
        const r = RootsBisection.find(F, t0, t1, f0, f1, maxIterations);
        io.outInt(r.iterations);
        io.outReal(r.root);
    }, { exact: true });

    family.case('RootsBisection1.find', (io) => {
        const c = readBracketPoly(io, 4);
        const tMin = io.real();
        const tMax = io.real();
        const maxIterations = io.integer();
        const passValues = io.boolean();
        const F = (t: number) => hornerAt(c, t);
        const bisector = new RootsBisection1(maxIterations);
        const r = passValues
            ? bisector.find(F, tMin, tMax, F(tMin), F(tMax))
            : bisector.find(F, tMin, tMax);
        io.outInt(r.iterations);
        io.outReal(r.root);
        io.outReal(r.fAtRoot);
    }, { exact: true });

    // Deliberate deviation: with maxIterations == 1 upstream returns its
    // output references unassigned, so the caller reads back the value it
    // passed in (a sentinel here); the port returns explicit zeros. See
    // docs/UPSTREAM-FINDINGS.md, RootsBisection1.h item 1, issue #84.
    family.case('RootsBisection1.deviation.maxIterationsOne', (io) => {
        const c = readBracketPoly(io, 3);
        const tMin = io.real();
        const tMax = io.real();
        io.real();  // the C++ sentinel, unused by the port
        const bisector = new RootsBisection1(1);
        const r = bisector.find((t) => hornerAt(c, t), tMin, tMax);
        io.outInt(r.iterations);
        io.outReal(r.root);
        io.outReal(r.fAtRoot);
    }, { exact: true, deviation: 'issue #84, RootsBisection1 maxIterations == 1' });

    family.case('RootsBisection1.throwParity', (io) => {
        const c = readBracketPoly(io, 2);
        const tMin = io.real();
        const tMax = io.real();
        const bisector = new RootsBisection1(20);
        const r = bisector.find((t) => hornerAt(c, t), tMin, tMax);
        io.outInt(r.iterations);
        io.outReal(r.root);
        io.outReal(r.fAtRoot);
    }, { exact: true });

    family.case('RootsBisection2.find', (io) => {
        const a = io.real(), b = io.real(), c = io.real();
        const d = io.real(), e = io.real(), f = io.real();
        const xMaxIterations = io.integer();
        const yMaxIterations = io.integer();
        const F = (x: number, y: number) => x * x * x + a * x + b * y + c;
        const G = (x: number, y: number) => y * y * y + d * y + e * x + f;
        const bisector = new RootsBisection2(xMaxIterations, yMaxIterations);
        const r = bisector.find(F, G, -4, 4, -4, 4);
        io.outInt(r.iterations);
        io.outReal(r.xRoot);
        io.outReal(r.yRoot);
        io.outReal(r.fAtRoot);
        io.outReal(r.gAtRoot);
        io.outBool(bisector.noGuaranteeForRootBound());
    }, { exact: true });

    // Deliberate deviation: upstream's members keep the previous call's
    // values when the x-bisector does not write them. See
    // docs/UPSTREAM-FINDINGS.md, RootsBisection2.h item 2, issues #84, #152.
    family.case('RootsBisection2.deviation.staleOutputs', (io) => {
        const a = io.real();
        const c2 = io.real();
        const d = io.real(), e = io.real(), f = io.real();
        const c1 = io.real();
        const G = (x: number, y: number) => y * y * y + d * y + e * x + f;
        const bisector = new RootsBisection2(1, 20);

        const F1 = (x: number) => x * x * x + a * x + c1;
        const r1 = bisector.find((x) => F1(x), G, -4, 4, -4, 4);
        io.outInt(r1.iterations);
        io.outReal(r1.xRoot);
        io.outReal(r1.fAtRoot);

        const F2 = (x: number) => x * x * x + a * x + c2;
        const r2 = bisector.find((x) => F2(x), G, -4, 4, -4, 4);
        io.outInt(r2.iterations);
        io.outReal(r2.xRoot);
        io.outReal(r2.fAtRoot);
    }, { exact: true, deviation: 'issues #84 and #152, RootsBisection2 stale outputs' });

    family.case('RootsBrentsMethod.find', (io) => {
        const c = readBracketPoly(io, 4);
        const t0 = io.real();
        const t1 = io.real();
        const maxIterations = io.integer();
        const mode = io.integer();
        const scale = io.real();
        let negFTolerance = 0, posFTolerance = 0;
        let stepTTolerance = 0, convTTolerance = 0;
        switch (mode) {
            case 0: break;
            case 1: negFTolerance = -scale; posFTolerance = scale; break;
            case 2: stepTTolerance = scale; break;
            case 3: convTTolerance = scale; break;
            case 4: negFTolerance = scale; break;
            default: posFTolerance = -scale; break;
        }
        const r = RootsBrentsMethod.find((t) => hornerAt(c, t), t0, t1, maxIterations,
            negFTolerance, posFTolerance, stepTTolerance, convTTolerance);
        io.outBool(r.found);
        io.outReal(r.root);
    }, { exact: true });

    family.case('CubicRootsQR.solve', (io) => {
        const c = io.reals(3);
        const maxIterations = io.integer();
        const r = new CubicRootsQR().solve(maxIterations, c[0], c[1], c[2]);
        io.outInt(r.iterations);
        io.outInt(r.numRoots);
        io.outReals(r.roots);
    }, { exact: true });

    family.case('CubicRootsQR.solveMatrix', (io) => {
        const A: number[][] = [];
        for (let r = 0; r < 3; ++r) { A.push(io.reals(3)); }
        const maxIterations = io.integer();
        const r = new CubicRootsQR().solveMatrix(maxIterations, A);
        io.outInt(r.iterations);
        io.outInt(r.numRoots);
        io.outReals(r.roots);
        for (let i = 0; i < 3; ++i) { io.outReals(A[i]); }
    }, { exact: true });

    family.case('QuarticRootsQR.solve', (io) => {
        const c = io.reals(4);
        const maxIterations = io.integer();
        const r = new QuarticRootsQR().solve(maxIterations, c[0], c[1], c[2], c[3]);
        io.outInt(r.iterations);
        io.outInt(r.numRoots);
        io.outReals(r.roots);
    }, { exact: true });

    family.case('QuarticRootsQR.solveMatrix', (io) => {
        const A: number[][] = [];
        for (let r = 0; r < 4; ++r) { A.push(io.reals(4)); }
        const maxIterations = io.integer();
        const r = new QuarticRootsQR().solveMatrix(maxIterations, A);
        io.outInt(r.iterations);
        io.outInt(r.numRoots);
        io.outReals(r.roots);
        for (let i = 0; i < 4; ++i) { io.outReals(A[i]); }
    }, { exact: true });

    family.case('RootsLinear.solve', (io) => {
        const c = io.reals(2);
        emitRoots(io, RootsLinear.solve(c[0], c[1]));
    }, { exact: true });

    family.case('RootsLinear.solveMonic', (io) => {
        const c = io.reals(1);
        emitRoots(io, RootsLinear.solveMonic(c[0]));
    }, { exact: true });

    family.case('RootsLinear.solveRational', (io) => {
        const c = io.reals(2);
        emitRationalRoots(io, RootsLinear.solveRational(
            BSRational.fromNumber(c[0]), BSRational.fromNumber(c[1])));
        emitRationalRoots(io, RootsLinear.solveMonicRational(
            BSRational.fromNumber(c[0])));
    }, { exact: true });

    family.case('RootsQuadratic.solve.bisection', (io) => {
        const c = io.reals(3);
        emitRoots(io, RootsQuadratic.solve(true, c[0], c[1], c[2]));
    }, { exact: true });

    family.case('RootsQuadratic.solve.closedForm', (io) => {
        const c = io.reals(3);
        emitRoots(io, RootsQuadratic.solve(false, c[0], c[1], c[2]));
    }, { exact: true });

    family.case('RootsQuadratic.solveMonic', (io) => {
        const c = io.reals(2);
        const useBisection = io.boolean();
        emitRoots(io, RootsQuadratic.solveMonic(useBisection, c[0], c[1]));
    }, { exact: true });

    family.case('RootsQuadratic.solveDepressed', (io) => {
        const d0 = io.real();
        const useBisection = io.boolean();
        emitRoots(io, RootsQuadratic.solveDepressed(useBisection, d0));
    }, { exact: true });

    family.case('RootsQuadratic.computeDepressedRoots', (io) => {
        const d0 = io.real();
        const useBisection = io.boolean();
        emitRationalRoots(io, RootsQuadratic.computeDepressedRoots(useBisection,
            BSRational.fromNumber(d0)));
    }, { exact: true });

    family.case('RootsQuadratic.solveRational', (io) => {
        const c = io.reals(3);
        const useBisection = io.boolean();
        emitRationalRoots(io, RootsQuadratic.solveRational(useBisection,
            BSRational.fromNumber(c[0]), BSRational.fromNumber(c[1]),
            BSRational.fromNumber(c[2])));
        emitRationalRoots(io, RootsQuadratic.solveMonicRational(useBisection,
            BSRational.fromNumber(c[0]), BSRational.fromNumber(c[1])));
        emitRationalRoots(io, RootsQuadratic.solveDepressedRational(useBisection,
            BSRational.fromNumber(c[0])));
    }, { exact: true });

    family.case('RootsPolynomial.solveQuadratic', (io) => {
        const c = io.reals(3);
        emitRootMap(io, RootsPolynomial.solveQuadratic(c[0], c[1], c[2]));
    }, { exact: true });

    family.case('RootsPolynomial.solveCubic', (io) => {
        const c = io.reals(4);
        emitRootMap(io, RootsPolynomial.solveCubic(c[0], c[1], c[2], c[3]));
    }, { tol: LIBM_TOL });

    family.case('RootsPolynomial.solveQuartic', (io) => {
        const c = io.reals(5);
        emitRootMap(io, RootsPolynomial.solveQuartic(c[0], c[1], c[2], c[3], c[4]));
    }, { tol: LIBM_TOL });

    family.case('RootsPolynomial.getRootInfo', (io) => {
        const c = io.reals(5);
        for (const info of [
            RootsPolynomial.getRootInfoQuadratic(c[0], c[1], c[2]),
            RootsPolynomial.getRootInfoCubic(c[0], c[1], c[2], c[3]),
            RootsPolynomial.getRootInfoQuartic(c[0], c[1], c[2], c[3], c[4])
        ]) {
            io.outInt(info.length);
            for (const m of info) { io.outInt(m); }
        }
    }, { exact: true });

    family.case('RootsPolynomial.find', (io) => {
        const degree = io.integer();
        const c = io.reals(degree + 1);
        const maxIterations = io.integer();
        const roots = RootsPolynomial.find(degree, c, maxIterations);
        io.outInt(roots.length);
        io.outReals(roots);
    }, { exact: true });

    family.case('RootsPolynomial.findBounded', (io) => {
        const degree = io.integer();
        const c = io.reals(degree + 1);
        const tmin = io.real();
        const tmax = io.real();
        const maxIterations = io.integer();
        const r = RootsPolynomial.find(degree, c, tmin, tmax, maxIterations);
        io.outBool(r.found);
        io.outReal(r.root);
    }, { exact: true });

    family.case('RootsGeneralPolynomial.solve', (io) => {
        const degree = io.integer();
        const p = io.reals(degree + 1);
        const useThreading = io.boolean();
        const roots = RootsGeneralPolynomial.solve(p, useThreading);
        io.outInt(roots.length);
        io.outReals(roots);
    }, { exact: true });

    family.case('RootsGeneralPolynomial.solveRational', (io) => {
        const degree = io.integer();
        const rP = io.reals(degree + 1).map((x) => BSRational.fromNumber(x));
        const rRoots = RootsGeneralPolynomial.solveRational(rP, false);
        io.outInt(rRoots.length);
        for (const r of rRoots) { io.outReal(r.toNumber()); }
    }, { exact: true });

    // Deliberate deviation: with high-order zero coefficients upstream runs
    // the solver on a zero-padded rational polynomial. See
    // docs/UPSTREAM-FINDINGS.md, RootsGeneralPolynomial.h item 4, issue #340.
    family.case('RootsGeneralPolynomial.deviation.zeroPadding', (io) => {
        const degree = io.integer();
        const numZeros = io.integer();
        const p = io.reals(degree + 1);
        for (let i = 0; i < numZeros; ++i) { p.push(0); }
        const roots = RootsGeneralPolynomial.solve(p, false);
        io.outInt(roots.length);
        io.outReals(roots);
    }, { exact: true, deviation: 'issue #340, RootsGeneralPolynomial zero padding' });

    family.finish();
});
