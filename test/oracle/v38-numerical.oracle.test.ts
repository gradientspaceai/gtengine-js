// Replays oracle/cpp/cases/v38-numerical.cpp. Keep the two files in the same
// order.
//
// Almost everything in group 38 is + - * / sqrt on doubles, so the cases are
// declared exact. The exceptions are stated at each case: the non-iterative
// NISymmetricEigensolver3x3 (std::acos, std::cos), one Integration case with
// a deliberately libm integrand, and the FPInterval cases, where the port
// emulates the directed rounding that upstream obtains from std::fesetround.
import { describe } from 'vitest';
import { BandedMatrix } from '../../src/BandedMatrix.js';
import { FPInterval } from '../../src/FPInterval.js';
import { GaussianElimination } from '../../src/GaussianElimination.js';
import { Integration } from '../../src/Integration.js';
import { LCPSolver } from '../../src/LCPSolver.js';
import { SymmetricEigensolver2x2 } from '../../src/SymmetricEigensolver2x2.js';
import {
    NISymmetricEigensolver3x3, SortEigenstuff, SymmetricEigensolver3x3,
    type EigenBasis3, type EigenTriple
} from '../../src/SymmetricEigensolver3x3.js';
import { UnsymmetricEigenvalues } from '../../src/UnsymmetricEigenvalues.js';
import { OracleFamily, type OracleIO } from './harness.js';

// The replay of RecordSym3: the (ignored) generator mode, then the six
// unique entries a00, a01, a02, a11, a12, a22.
function sym3(io: OracleIO): [number, number, number, number, number, number] {
    io.integer();
    return [io.real(), io.real(), io.real(), io.real(), io.real(), io.real()];
}

function outSym3Result(io: OracleIO, evals: EigenTriple, evecs: EigenBasis3,
    tol?: number): void {
    for (let i = 0; i < 3; ++i) { io.outReal(evals[i], tol); }
    for (let i = 0; i < 3; ++i) {
        for (let j = 0; j < 3; ++j) { io.outReal(evecs[i][j], tol); }
    }
}

describe('oracle: v38-numerical', () => {
    const family = new OracleFamily('v38-numerical');

    // ------------------------------------------ SymmetricEigensolver2x2.h

    family.case('SymmetricEigensolver2x2.solve', (io) => {
        io.integer();
        const a00 = io.real();
        const a01 = io.real();
        const a11 = io.real();
        const sortType = io.integer();
        const r = new SymmetricEigensolver2x2().solve(a00, a01, a11, sortType);
        io.outReal(r.evals[0]);
        io.outReal(r.evals[1]);
        io.outReal(r.evecs[0][0]);
        io.outReal(r.evecs[0][1]);
        io.outReal(r.evecs[1][0]);
        io.outReal(r.evecs[1][1]);
    }, { exact: true });

    // ------------------------------------------ SymmetricEigensolver3x3.h

    family.case('SortEigenstuff.sort', (io) => {
        const evals: EigenTriple = [io.real(), io.real(), io.real()];
        const evecs: EigenBasis3 = [
            [io.real(), io.real(), io.real()],
            [io.real(), io.real(), io.real()],
            [io.real(), io.real(), io.real()]
        ];
        const isRotation = io.boolean();
        const sortType = io.integer();
        new SortEigenstuff().sort(sortType, isRotation, evals, evecs);
        outSym3Result(io, evals, evecs);
    }, { exact: true });

    family.case('SymmetricEigensolver3x3.solve', (io) => {
        const a = sym3(io);
        const aggressive = io.boolean();
        const sortType = io.integer();
        const r = new SymmetricEigensolver3x3().solve(a[0], a[1], a[2], a[3],
            a[4], a[5], aggressive, sortType);
        io.outInt(r.iterations);
        outSym3Result(io, r.evals, r.evecs);
    }, { exact: true });

    family.case('SymmetricEigensolver3x3.solve.outOfBandScale', (io) => {
        const a: number[] = [io.real(), io.real(), io.real(), io.real(),
            io.real(), io.real()];
        const aggressive = io.boolean();
        const sortType = io.integer();
        const r = new SymmetricEigensolver3x3().solve(a[0], a[1], a[2], a[3],
            a[4], a[5], aggressive, sortType);
        io.outInt(r.iterations);
        outSym3Result(io, r.evals, r.evecs);
    }, { exact: true, deviation: 'issue #379: GetCosSin lacks the maxAbsComp rescaling' });

    // The only libm-dependent eigensolver: std::acos and std::cos decide the
    // eigenvalues, and the eigenvectors are built from them. The generator
    // keeps the eigenvalue gaps above 1e-3 of the spectral radius so that
    // the libm-decided control flow is the same on both sides; the residual
    // difference is pure value error.
    family.case('NISymmetricEigensolver3x3.solve', (io) => {
        const a = sym3(io);
        const sortType = io.integer();
        const r = new NISymmetricEigensolver3x3().solve(a[0], a[1], a[2], a[3],
            a[4], a[5], sortType);
        outSym3Result(io, r.evals, r.evecs);
    }, { tol: 1e-12 });

    family.case('NISymmetricEigensolver3x3.solve.degenerate', (io) => {
        const a = sym3(io);
        const sortType = io.integer();
        const r = new NISymmetricEigensolver3x3().solve(a[0], a[1], a[2], a[3],
            a[4], a[5], sortType);
        outSym3Result(io, r.evals, r.evecs);
    }, { exact: true });

    // -------------------------------------------- UnsymmetricEigenvalues.h

    function unsymmetric(io: OracleIO): void {
        const n = io.integer();
        const maxIterations = io.integer();
        const m = io.reals(n * n);
        const sortType = io.integer();
        const solver = new UnsymmetricEigenvalues(n, maxIterations);
        const numIterations = solver.solve(m, sortType);
        io.outInt(numIterations);
        const r = solver.getEigenvalues();
        io.outInt(r.numEigenvalues);
        for (let i = 0; i < r.numEigenvalues; ++i) { io.outReal(r.eigenvalues[i]); }
    }

    family.case('UnsymmetricEigenvalues.solve', unsymmetric, { exact: true });

    family.case('UnsymmetricEigenvalues.solve.trailingBlock', unsymmetric,
        { exact: true, deviation: 'issue #42: the packing loop drops A(N-1,N-1)' });

    family.case('UnsymmetricEigenvalues.solve.invalidSize', (io) => {
        const n = io.integer();
        const maxIterations = io.integer();
        const storage = Math.max(n, 1);
        const m = io.reals(storage * storage);
        const sortType = io.integer();
        const solver = new UnsymmetricEigenvalues(n, maxIterations);
        io.outInt(solver.solve(m, sortType));
        const r = solver.getEigenvalues();
        io.outInt(r.numEigenvalues);
        for (let i = 0; i < r.numEigenvalues; ++i) { io.outReal(r.eigenvalues[i]); }
    }, { exact: true });

    family.case('UnsymmetricEigenvalues.solve.cycling', (io) => {
        const t = io.real();
        const maxIterations = io.integer();
        const sortType = io.integer();
        const m = [0, t, 0, t, 0, t, 0, t, 0];
        const solver = new UnsymmetricEigenvalues(3, maxIterations);
        io.outInt(solver.solve(m, sortType));
        const r = solver.getEigenvalues();
        io.outInt(r.numEigenvalues);
        for (let i = 0; i < r.numEigenvalues; ++i) { io.outReal(r.eigenvalues[i]); }
    }, { exact: true });

    // Two exact reproductions of the preserved non-convergence defect
    // (issue #476): FrancisQRStep has no exceptional shift.
    family.case('UnsymmetricEigenvalues.solve.nonConvergence', (io) => {
        io.integer();
        const maxIterations = io.integer();
        const m = io.reals(9);
        const sortType = io.integer();
        const solver = new UnsymmetricEigenvalues(3, maxIterations);
        io.outInt(solver.solve(m, sortType));
        const r = solver.getEigenvalues();
        io.outInt(r.numEigenvalues);
        for (let i = 0; i < r.numEigenvalues; ++i) { io.outReal(r.eigenvalues[i]); }
    }, { exact: true });

    // ------------------------------------------------------ BandedMatrix.h

    // The replay of RecordBanded: the shape, then the diagonal, the lower
    // bands and the upper bands.
    function banded(io: OracleIO): BandedMatrix {
        const size = io.integer();
        const numLBands = io.integer();
        const numUBands = io.integer();
        io.integer();
        const matrix = new BandedMatrix(size, numLBands, numUBands);
        const d = matrix.getDBand();
        for (let i = 0; i < size; ++i) { d[i] = io.real(); }
        const l = matrix.getLBands();
        for (let k = 0; k < numLBands; ++k) {
            for (let i = 0; i < size - 1 - k; ++i) { l[k][i] = io.real(); }
        }
        const u = matrix.getUBands();
        for (let k = 0; k < numUBands; ++k) {
            for (let i = 0; i < size - 1 - k; ++i) { u[k][i] = io.real(); }
        }
        return matrix;
    }

    function outBanded(io: OracleIO, matrix: BandedMatrix): void {
        for (const e of matrix.getDBand()) { io.outReal(e); }
        for (const band of matrix.getLBands()) {
            for (const e of band) { io.outReal(e); }
        }
        for (const band of matrix.getUBands()) {
            for (const e of band) { io.outReal(e); }
        }
    }

    family.case('BandedMatrix.accessor', (io) => {
        const matrix = banded(io);
        const size = matrix.getSize();
        io.outInt(size);
        for (let r = -1; r <= size; ++r) {
            for (let c = -1; c <= size; ++c) { io.outReal(matrix.get(r, c)); }
        }
        for (let r = -1; r <= size; ++r) {
            for (let c = -1; c <= size; ++c) {
                matrix.set(r, c, 1 + r * (size + 2) + c);
            }
        }
        outBanded(io, matrix);
    }, { exact: true });

    family.case('BandedMatrix.choleskyFactor', (io) => {
        const matrix = banded(io);
        io.outBool(matrix.choleskyFactor());
        outBanded(io, matrix);
    }, { exact: true });

    family.case('BandedMatrix.solveSystem.vector', (io) => {
        const matrix = banded(io);
        const rhs = io.reals(matrix.getSize());
        io.outBool(matrix.solveSystem(rhs));
        io.outReals(rhs);
        outBanded(io, matrix);
    }, { exact: true });

    family.case('BandedMatrix.solveSystem.matrix', (io) => {
        const matrix = banded(io);
        const numBColumns = io.integer();
        const rowMajor = io.boolean();
        const rhs = io.reals(matrix.getSize() * numBColumns);
        io.outBool(matrix.solveSystemMatrix(rhs, numBColumns, rowMajor));
        io.outReals(rhs);
        outBanded(io, matrix);
    }, { exact: true });

    family.case('BandedMatrix.computeInverse', (io) => {
        const matrix = banded(io);
        const rowMajor = io.boolean();
        const size = matrix.getSize();
        const inverse = new Array<number>(size * size).fill(0);
        const success = matrix.computeInverse(inverse, rowMajor);
        io.outBool(success);
        if (success) { io.outReals(inverse); }
        outBanded(io, matrix);
    }, { exact: true });

    // Upstream defect found by this oracle: ComputeInverse eliminates without
    // pivoting and only rejects an exactly zero pivot, so a matrix with a
    // vanishing leading principal minor gets a wrong inverse and a 'true'
    // return. The port preserves that, so the two sides agree bit for bit;
    // see oracle/reports/v38-numerical.md.
    family.case('BandedMatrix.computeInverse.zeroLeadingMinor', (io) => {
        const matrix = banded(io);
        const rowMajor = io.boolean();
        const inverse = new Array<number>(36).fill(0);
        io.outBool(matrix.computeInverse(inverse, rowMajor));
        io.outReals(inverse);
    }, { exact: true });

    family.case('BandedMatrix.invalidShape', (io) => {
        const size = io.integer();
        const numLBands = io.integer();
        const numUBands = io.integer();
        const matrix = new BandedMatrix(size, numLBands, numUBands);
        io.outInt(matrix.getSize());
        io.outInt(matrix.getDBand().length);
        io.outInt(matrix.getLBands().length);
        io.outInt(matrix.getUBands().length);
        io.outBool(matrix.choleskyFactor());
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 3; ++c) { io.outReal(matrix.get(r, c)); }
        }
    }, { exact: true });

    // ----------------------------------------------- GaussianElimination.h

    family.case('GaussianElimination.compute', (io) => {
        const n = io.integer();
        io.integer();
        const m = io.reals(n * n);
        const numCols = io.integer();
        const B = io.reals(n);
        const C = io.reals(n * numCols);
        const r = new GaussianElimination().compute(n, m,
            { wantInverse: true, B, C, numCols });
        io.outBool(r.invertible);
        io.outReal(r.determinant);
        io.outReals(r.inverseM!);
        io.outReals(r.X!);
        io.outReals(r.Y!);
    }, { exact: true });

    family.case('GaussianElimination.compute.partial', (io) => {
        const n = io.integer();
        io.integer();
        const m = io.reals(n * n);
        const selector = io.integer();
        const numCols = io.integer();
        const B = io.reals(n);
        const C = io.reals(n * numCols);
        const solver = new GaussianElimination();
        let r;
        if (selector === 0) {
            r = solver.compute(n, m);
        } else if (selector === 1) {
            r = solver.compute(n, m, { B });
        } else if (selector === 2) {
            r = solver.compute(n, m, { C, numCols });
        } else {
            r = solver.compute(n, m, { wantInverse: true });
        }
        io.outBool(r.invertible);
        io.outReal(r.determinant);
        if (selector === 3) { io.outReals(r.inverseM!); }
        if (selector === 1) { io.outReals(r.X!); }
        if (selector === 2) { io.outReals(r.Y!); }
    }, { exact: true });

    // Preserved upstream defect (issue #375): a subnormal-magnitude matrix
    // passes the pivot test but 1 / pivot overflows, so the inverse comes
    // back NaN with invertible = true. Compared bit for bit.
    family.case('GaussianElimination.compute.subnormal', (io) => {
        const n = io.integer();
        const m = io.reals(n * n);
        const r = new GaussianElimination().compute(n, m, { wantInverse: true });
        io.outBool(r.invertible);
        io.outReal(r.determinant);
        io.outReals(r.inverseM!);
    }, { exact: true });

    family.case('GaussianElimination.compute.invalidInput', (io) => {
        const which = io.integer();
        const n = io.integer();
        const storage = Math.max(n, 1);
        const m = io.reals(storage * storage);
        const C = io.reals(storage);
        const solver = new GaussianElimination();
        const r = which === 0
            ? solver.compute(n, m)
            : solver.compute(n, m, { C, numCols: 0 });
        io.outBool(r.invertible);
        io.outReal(r.determinant);
    }, { exact: true });

    // -------------------------------------------------------- Integration.h

    // The replay of DrawIntegrand and Evaluate. The two sides must write the
    // integrand identically, so that the case compares the quadrature rule
    // and not the integrand.
    interface Integrand { kind: number; degree: number; coefficient: number[]; }

    function integrand(io: OracleIO, kind: number): Integrand {
        const degree = io.integer();
        return { kind, degree, coefficient: io.reals(degree + 1) };
    }

    function evaluate(f: Integrand): (t: number) => number {
        return (t: number): number => {
            let result = f.coefficient[f.degree];
            for (let i = f.degree - 1; i >= 0; --i) {
                result = result * t + f.coefficient[i];
            }
            if (f.kind === 1) {
                result = result / (t * t + 1);
            } else if (f.kind === 2) {
                result = Math.exp(-t * t) * result;
            }
            return result;
        };
    }

    family.case('Integration.trapezoidRule', (io) => {
        const numSamples = io.integer();
        const a = io.real();
        const b = io.real();
        const kind = io.integer();
        const f = integrand(io, kind);
        io.outReal(Integration.trapezoidRule(numSamples, a, b, evaluate(f)));
    }, { exact: true });

    family.case('Integration.romberg', (io) => {
        const order = io.integer();
        const a = io.real();
        const b = io.real();
        const kind = io.integer();
        const f = integrand(io, kind);
        io.outReal(Integration.romberg(order, a, b, evaluate(f)));
    }, { exact: true });

    family.case('Integration.computeQuadratureInfo', (io) => {
        const degree = io.integer();
        const r = Integration.computeQuadratureInfo(degree);
        io.outInt(r.roots.length);
        io.outReals(r.roots);
        io.outInt(r.coefficients.length);
        io.outReals(r.coefficients);
    }, { exact: true });

    family.case('Integration.gaussianQuadrature', (io) => {
        const degree = io.integer();
        const a = io.real();
        const b = io.real();
        const kind = io.integer();
        const f = integrand(io, kind);
        const q = Integration.computeQuadratureInfo(degree);
        io.outReal(Integration.gaussianQuadrature(q.roots, q.coefficients,
            a, b, evaluate(f)));
        io.outInt(q.roots.length);
        io.outReals(q.roots);
        io.outReals(q.coefficients);
    }, { exact: true });

    // The one libm integrand of the group: Math.exp and std::exp are both
    // allowed to be an ulp off, and the quadrature sums amplify that.
    family.case('Integration.libmIntegrand', (io) => {
        const order = io.integer();
        const numSamples = io.integer();
        const a = io.real();
        const b = io.real();
        const f = integrand(io, 2);
        io.outReal(Integration.trapezoidRule(numSamples, a, b, evaluate(f)));
        io.outReal(Integration.romberg(order, a, b, evaluate(f)));
    }, { tol: 1e-12 });

    // --------------------------------------------------------- LCPSolver.h

    // The replay of RecordLcp: the dimension, the (ignored) mode, q, then M.
    function lcp(io: OracleIO): { n: number; q: number[]; m: number[] } {
        const n = io.integer();
        io.integer();
        const q = io.reals(n);
        const m = io.reals(n * n);
        return { n, q, m };
    }

    function runLcp(io: OracleIO, problem: { n: number; q: number[]; m: number[] },
        maxIterations: number): void {
        const solver = new LCPSolver(problem.n);
        if (maxIterations !== 0) { solver.setMaxIterations(maxIterations); }
        io.outInt(solver.getMaxIterations());
        const r = solver.solve(problem.q, problem.m);
        io.outBool(r.success);
        io.outInt(r.result);
        io.outInt(solver.getNumIterations());
        io.outReals(r.w);
        io.outReals(r.z);
    }

    family.case('LCPSolver.solve', (io) => {
        const problem = lcp(io);
        const maxIterations = io.integer();
        runLcp(io, problem, maxIterations);
    }, { exact: true });

    family.case('LCPSolver.solve.positiveDefinite', (io) => {
        runLcp(io, lcp(io), 0);
    }, { exact: true });

    family.case('LCPSolver.solve.trivial', (io) => {
        runLcp(io, lcp(io), 0);
    }, { exact: true });

    family.case('LCPSolver.solve.noSolution', (io) => {
        runLcp(io, lcp(io), 0);
    }, { exact: true });

    family.case('LCPSolver.solve.maxIterations', (io) => {
        const problem = lcp(io);
        const maxIterations = io.integer();
        runLcp(io, problem, maxIterations);
    }, { exact: true });

    // Upstream defect found by this oracle: a round-off-sized pivot accepted
    // by the ratio test turns a provably infeasible LCP into a reported
    // solution. The port preserves it, so the two sides agree bit for bit;
    // see oracle/reports/v38-numerical.md.
    family.case('LCPSolver.solve.roundoffPivot', (io) => {
        const problem = lcp(io);
        const maxIterations = io.integer();
        runLcp(io, problem, maxIterations);
    }, { exact: true });

    family.case('LCPSolver.solve.invalidInput', (io) => {
        const n = io.integer();
        const shortfall = io.integer();
        const shortenQ = io.boolean();
        const qSize = shortenQ ? Math.max(0, n - shortfall) : n;
        const mSize = shortenQ ? n * n : Math.max(0, n * n - shortfall);
        const q = io.reals(qSize);
        const m = io.reals(mSize);
        const r = new LCPSolver(n).solve(q, m);
        io.outBool(r.success);
        io.outInt(r.result);
        io.outReals(r.w);
        io.outReals(r.z);
    }, { exact: true });

    // -------------------------------------------------------- FPInterval.h

    function outInterval(io: OracleIO, w: FPInterval): void {
        io.outReal(w.get(0));
        io.outReal(w.get(1));
    }

    family.case('FPInterval.leaf', (io) => {
        const u = io.real();
        io.integer();
        const v = io.real();
        outInterval(io, FPInterval.add(u, v));
        outInterval(io, FPInterval.sub(u, v));
        outInterval(io, FPInterval.mul(u, v));
        outInterval(io, FPInterval.div(u, v));
    }, { exact: true });

    family.case('FPInterval.internal', (io) => {
        io.integer();
        const u0 = io.real();
        const u1 = io.real();
        io.integer();
        const v0 = io.real();
        const v1 = io.real();
        outInterval(io, FPInterval.add(u0, u1, v0, v1));
        outInterval(io, FPInterval.sub(u0, u1, v0, v1));
        outInterval(io, FPInterval.mul(u0, u1, v0, v1));
        outInterval(io, FPInterval.mul2(u0, u1, v0, v1));
        outInterval(io, FPInterval.div(u0, u1, v0, v1));
        outInterval(io, FPInterval.reciprocal(v0, v1));
        outInterval(io, FPInterval.reciprocalDown(v1));
        outInterval(io, FPInterval.reciprocalUp(v0));
        outInterval(io, FPInterval.reals());
        outInterval(io, FPInterval.fromEndpoints([v0, v1]));
        outInterval(io, new FPInterval(u0));
        outInterval(io, new FPInterval());
    }, { exact: true });

    family.case('FPInterval.operators', (io) => {
        io.integer();
        const u = new FPInterval(io.real(), io.real());
        io.integer();
        const v = new FPInterval(io.real(), io.real());
        io.integer();
        const d = new FPInterval(io.real(), io.real());
        io.integer();
        const scalar = io.real();

        // The port has no unary plus; upstream's returns its argument.
        outInterval(io, u);
        outInterval(io, u.negate());
        outInterval(io, u.add(v));
        outInterval(io, u.add(scalar));
        outInterval(io, u.add(scalar));
        outInterval(io, u.sub(v));
        outInterval(io, u.sub(scalar));
        outInterval(io, FPInterval.scalarSub(scalar, u));
        outInterval(io, u.mul(v));
        outInterval(io, u.mul(scalar));
        outInterval(io, u.mul(scalar));
        outInterval(io, u.div(d));
        outInterval(io, u.div(scalar));
        outInterval(io, FPInterval.scalarDiv(scalar, d));
    }, { exact: true });

    family.case('FPInterval.productBounds', (io) => {
        io.integer();
        const u = [io.real(), io.real()];
        io.integer();
        const v = [io.real(), io.real()];
        io.outReal(FPInterval.productLowerBound(u, v));
        io.outReal(FPInterval.productUpperBound(u, v));
    }, { exact: true });

    family.case('FPInterval.directedRounding', (io) => {
        const u0 = io.real();
        const u1 = io.real();
        const v0 = io.real();
        const v1 = io.real();
        const u = new FPInterval(u0, u1);
        const v = new FPInterval(v0, v1);
        outInterval(io, FPInterval.add(u0, v0));
        outInterval(io, FPInterval.sub(u0, v0));
        outInterval(io, FPInterval.mul(u0, v0));
        outInterval(io, FPInterval.div(u0, v0));
        outInterval(io, u.add(v));
        outInterval(io, u.sub(v));
        outInterval(io, u.mul(v));
        outInterval(io, u.div(v));
        outInterval(io, FPInterval.reciprocal(v0, v1));
        outInterval(io, FPInterval.reciprocal(-v1, 0));
        outInterval(io, FPInterval.reciprocalDown(0));
        outInterval(io, FPInterval.reciprocalUp(0));
    }, {
        exact: true,
        deviation: 'the port emulates the directed rounding that upstream takes '
            + 'from std::fesetround; see the PORT DEVIATION note in src/FPInterval.ts'
    });

    family.finish();
});
