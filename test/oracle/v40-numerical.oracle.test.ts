// Replays oracle/cpp/cases/v40-numerical.cpp (verify group 40, numerical).
// Keep the two files in the same order.
//
// Every case is arithmetic only (+ - * / and sqrt), so every case is
// declared { exact: true }. The residual functions F and the Jacobians J
// below are the same expressions in the same association as the C++ structs
// Rosenbrock, CircleFit and LinearModel.
import { describe } from 'vitest';
import { GaussNewtonMinimizer } from '../../src/GaussNewtonMinimizer.js';
import type { GaussNewtonMinimizerResult } from '../../src/GaussNewtonMinimizer.js';
import {
    LevenbergMarquardtMinimizer
} from '../../src/LevenbergMarquardtMinimizer.js';
import type {
    LevenbergMarquardtMinimizerResult
} from '../../src/LevenbergMarquardtMinimizer.js';
import { LinearSystem } from '../../src/LinearSystem.js';
import type { LinearSystemSparseEntry } from '../../src/LinearSystem.js';
import { Matrix } from '../../src/Matrix.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// Reads an n-by-n matrix recorded row-major by FillMatrix/FillSymmetric.
function squareMatrix(io: OracleIO, n: number): Matrix {
    const m = new Matrix(n, n);
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) { m.set(r, c, io.real()); }
    }
    return m;
}

function emitVectorExact(io: OracleIO, x: readonly number[]): void {
    for (const value of x) { io.outRealExact(value); }
}

// ---- the minimizer problems, mirroring the C++ structs --------------------

interface Problem {
    numP: number;
    numF: number;
    p0: Vector;
    f: (p: Vector, f: Vector) => void;
    j: (p: Vector, j: Matrix) => void;
}

// Rosenbrock: F = (a*(p1 - p0^2), b - p0). Reads a, b, p0[0], p0[1].
function makeRosenbrock(io: OracleIO): Problem {
    const a = io.real();
    const b = io.real();
    const p0 = io.vec(2);
    return {
        numP: 2,
        numF: 2,
        p0,
        f: (p, out) => {
            out.set(0, a * (p.get(1) - p.get(0) * p.get(0)));
            out.set(1, b - p.get(0));
        },
        j: (p, out) => {
            out.set(0, 0, -2.0 * a * p.get(0));
            out.set(0, 1, a);
            out.set(1, 0, -1.0);
            out.set(1, 1, 0.0);
        }
    };
}

interface CircleProblem extends Problem {
    normal: (p: Vector, jtj: Matrix, negJTF: Vector) => void;
}

// Circle fit by algebraic residuals. Reads n, then n points, then p0.
function makeCircleFit(io: OracleIO): CircleProblem {
    const n = io.integer();
    const points: number[][] = [];
    for (let i = 0; i < n; ++i) { points.push([io.real(), io.real()]); }
    const p0 = io.vec(3);
    return {
        numP: 3,
        numF: n,
        p0,
        f: (p, out) => {
            for (let i = 0; i < points.length; ++i) {
                const dx = points[i][0] - p.get(0);
                const dy = points[i][1] - p.get(1);
                out.set(i, dx * dx + dy * dy - p.get(2) * p.get(2));
            }
        },
        j: (p, out) => {
            for (let i = 0; i < points.length; ++i) {
                const dx = points[i][0] - p.get(0);
                const dy = points[i][1] - p.get(1);
                out.set(i, 0, -2.0 * dx);
                out.set(i, 1, -2.0 * dy);
                out.set(i, 2, -2.0 * p.get(2));
            }
        },
        normal: (p, jtj, negJTF) => {
            for (let r = 0; r < 3; ++r) {
                negJTF.set(r, 0.0);
                for (let c = 0; c < 3; ++c) { jtj.set(r, c, 0.0); }
            }
            for (let i = 0; i < points.length; ++i) {
                const dx = points[i][0] - p.get(0);
                const dy = points[i][1] - p.get(1);
                const fi = dx * dx + dy * dy - p.get(2) * p.get(2);
                const d = [-2.0 * dx, -2.0 * dy, -2.0 * p.get(2)];
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 3; ++c) {
                        jtj.set(r, c, jtj.get(r, c) + d[r] * d[c]);
                    }
                    negJTF.set(r, negJTF.get(r) - fi * d[r]);
                }
            }
        }
    };
}

// The linear model whose Jacobian is rank deficient in modes 0 and 1.
// Reads n, then n abscissas, then n ordinates, then p0.
function makeLinearModel(io: OracleIO, mode: number): Problem {
    const n = io.integer();
    const t: number[] = [];
    for (let i = 0; i < n; ++i) { t.push(io.real()); }
    const y: number[] = [];
    for (let i = 0; i < n; ++i) { y.push(io.real()); }
    const p0 = io.vec(2);
    return {
        numP: 2,
        numF: n,
        p0,
        f: (p, out) => {
            for (let i = 0; i < t.length; ++i) {
                const model = (mode === 0 ? p.get(0) * t[i]
                    : mode === 1 ? (p.get(0) + p.get(1)) * t[i]
                    : p.get(0) * t[i] + p.get(1));
                out.set(i, model - y[i]);
            }
        },
        j: (_p, out) => {
            for (let i = 0; i < t.length; ++i) {
                out.set(i, 0, t[i]);
                out.set(i, 1, mode === 0 ? 0.0 : mode === 1 ? t[i] : 1.0);
            }
        }
    };
}

function emitGN(io: OracleIO, r: GaussNewtonMinimizerResult): void {
    io.outInt(r.numIterations);
    io.outBool(r.converged);
    io.outRealExact(r.minError);
    io.outRealExact(r.minErrorDifference);
    io.outRealExact(r.minUpdateLength);
    io.outVecExact(r.minLocation);
}

function emitLM(io: OracleIO, r: LevenbergMarquardtMinimizerResult): void {
    io.outInt(r.numIterations);
    io.outInt(r.numAdjustments);
    io.outBool(r.converged);
    io.outRealExact(r.minError);
    io.outRealExact(r.minErrorDifference);
    io.outRealExact(r.minUpdateLength);
    io.outVecExact(r.minLocation);
}

describe('oracle: v40-numerical', () => {
    const family = new OracleFamily('v40-numerical');

    // ---- LinearSystem: the closed-form fixed-size solvers ------------------

    for (const n of [2, 3, 4] as const) {
        family.case(`LinearSystem.solve.${n}x${n}`, (io) => {
            const A = squareMatrix(io, n);
            const B = io.vec(n);
            const r = n === 2 ? LinearSystem.solve2x2(A, B)
                : n === 3 ? LinearSystem.solve3x3(A, B)
                : LinearSystem.solve4x4(A, B);
            io.outBool(r.invertible);
            io.outVecExact(r.X);
        }, { exact: true });
    }

    // ---- LinearSystem: Gaussian elimination --------------------------------

    family.case('LinearSystem.solve.nxn', (io) => {
        const n = io.integer();
        const A = io.reals(n * n);
        const B = io.reals(n);
        const r = LinearSystem.solve(n, A, B);
        io.outBool(r.invertible);
        emitVectorExact(io, r.X);
    }, { exact: true });

    family.case('LinearSystem.solve.nxm', (io) => {
        const n = io.integer();
        const m = io.integer();
        const A = io.reals(n * n);
        const B = io.reals(n * m);
        const r = LinearSystem.solveMultiple(n, m, A, B);
        io.outBool(r.invertible);
        emitVectorExact(io, r.X);
    }, { exact: true });

    // ---- LinearSystem: the tridiagonal solvers -----------------------------

    family.case('LinearSystem.solveTridiagonal', (io) => {
        const n = io.integer();
        const sub = io.reals(n - 1);
        const diag = io.reals(n);
        const sup = io.reals(n - 1);
        const B = io.reals(n);
        const r = LinearSystem.solveTridiagonal(n, sub, diag, sup, B);
        io.outBool(r.solved);
        emitVectorExact(io, r.X);
    }, { exact: true });

    family.case('LinearSystem.solveConstantTridiagonal', (io) => {
        const n = io.integer();
        const sub = io.real();
        const diag = io.real();
        const sup = io.real();
        const B = io.reals(n);
        const r = LinearSystem.solveConstantTridiagonal(n, sub, diag, sup, B);
        io.outBool(r.solved);
        emitVectorExact(io, r.X);
    }, { exact: true });

    // ---- LinearSystem: conjugate gradient ----------------------------------

    family.case('LinearSystem.solveSymmetricCG.dense', (io) => {
        const n = io.integer();
        const A = io.reals(n * n);
        const B = io.reals(n);
        const maxIterations = io.integer();
        const tolerance = io.real();
        const r = LinearSystem.solveSymmetricCG(n, A, B, maxIterations,
            tolerance);
        io.outInt(r.iterations);
        emitVectorExact(io, r.X);
    }, { exact: true });

    family.case('LinearSystem.solveSymmetricCG.sparse', (io) => {
        const n = io.integer();
        const entries: LinearSystemSparseEntry[] = [];
        for (let r = 0; r < n; ++r) {
            entries.push({ row: r, col: r, value: io.real() });
        }
        for (let r = 0; r < n; ++r) {
            for (let c = r + 1; c < n; ++c) {
                const present = io.boolean();
                const value = io.real();
                if (present) { entries.push({ row: r, col: c, value }); }
            }
        }
        const B = io.reals(n);
        const maxIterations = io.integer();
        const tolerance = io.real();
        const r = LinearSystem.solveSymmetricCGSparse(n, entries, B,
            maxIterations, tolerance);
        io.outInt(r.iterations);
        emitVectorExact(io, r.X);
    }, { exact: true });

    // A preserved upstream defect (issue #261): B = 0 makes alpha = 0/0 and
    // overwrites the exact solution with NaN. The NaN pattern is emitted as
    // booleans because the harness treats any NaN as equal to any NaN.
    family.case('LinearSystem.solveSymmetricCG.zeroRHS', (io) => {
        const n = io.integer();
        const A = io.reals(n * n);
        const B = io.reals(n);
        const maxIterations = io.integer();
        const r = LinearSystem.solveSymmetricCG(n, A, B, maxIterations, 1e-8);
        io.outInt(r.iterations);
        for (let i = 0; i < n; ++i) { io.outBool(Number.isNaN(r.X[i])); }
    }, { exact: true });

    // ---- GaussNewtonMinimizer ---------------------------------------------

    family.case('GaussNewtonMinimizer.minimize.rosenbrock', (io) => {
        const problem = makeRosenbrock(io);
        const maxIterations = io.integer();
        const updateLengthTolerance = io.real();
        const errorDifferenceTolerance = io.real();
        const minimizer = GaussNewtonMinimizer.fromJFunction(problem.numP,
            problem.numF, problem.f, problem.j);
        emitGN(io, minimizer.minimize(problem.p0, maxIterations,
            updateLengthTolerance, errorDifferenceTolerance));
    }, { exact: true });

    family.case('GaussNewtonMinimizer.minimize.circleFit', (io) => {
        const problem = makeCircleFit(io);
        const maxIterations = io.integer();
        const updateLengthTolerance = io.real();
        const errorDifferenceTolerance = io.real();
        const minimizer = GaussNewtonMinimizer.fromJFunction(problem.numP,
            problem.numF, problem.f, problem.j);
        emitGN(io, minimizer.minimize(problem.p0, maxIterations,
            updateLengthTolerance, errorDifferenceTolerance));
    }, { exact: true });

    family.case('GaussNewtonMinimizer.minimize.jPlus', (io) => {
        const problem = makeCircleFit(io);
        const maxIterations = io.integer();
        const minimizer = GaussNewtonMinimizer.fromJPlusFunction(problem.numP,
            problem.numF, problem.f, problem.normal);
        emitGN(io, minimizer.minimize(problem.p0, maxIterations, 0, 0));
    }, { exact: true });

    family.case('GaussNewtonMinimizer.minimize.rankDeficient', (io) => {
        const problem = makeLinearModel(io, io.index % 3);
        const maxIterations = io.integer();
        const minimizer = GaussNewtonMinimizer.fromJFunction(problem.numP,
            problem.numF, problem.f, problem.j);
        emitGN(io, minimizer.minimize(problem.p0, maxIterations, 0, 0));
    }, { exact: true });

    family.case('GaussNewtonMinimizer.minimize.zeroResidual', (io) => {
        const problem = makeLinearModel(io, 2);
        const maxIterations = io.integer();
        const minimizer = GaussNewtonMinimizer.fromJFunction(problem.numP,
            problem.numF, problem.f, problem.j);
        emitGN(io, minimizer.minimize(problem.p0, maxIterations, 0, 0));
    }, { exact: true });

    family.case('GaussNewtonMinimizer.construct.invalidDimensions', (io) => {
        const numP = io.integer();
        const numF = io.integer();
        const minimizer = GaussNewtonMinimizer.fromJFunction(numP, numF,
            (_p, f) => { f.set(0, 0); }, (_p, j) => { j.set(0, 0, 0); });
        io.outInt(minimizer.getNumPDimensions());
        io.outInt(minimizer.getNumFDimensions());
    }, { exact: true });

    // ---- LevenbergMarquardtMinimizer --------------------------------------
    // maxIterations is an aimed input: the C++ side chose the longest prefix
    // on which upstream never rebuilds -J^T*F from a stale residual.

    family.case('LevenbergMarquardtMinimizer.minimize.rosenbrock', (io) => {
        const problem = makeRosenbrock(io);
        const maxAdjustments = io.integer();
        const lambdaFactor = io.real();
        const lambdaAdjust = io.real();
        const maxIterations = io.integer();
        const minimizer = LevenbergMarquardtMinimizer.fromJFunction(
            problem.numP, problem.numF, problem.f, problem.j);
        emitLM(io, minimizer.minimize(problem.p0, maxIterations, 0, 0,
            lambdaFactor, lambdaAdjust, maxAdjustments));
    }, { exact: true });

    family.case('LevenbergMarquardtMinimizer.minimize.circleFit', (io) => {
        const problem = makeCircleFit(io);
        const maxAdjustments = io.integer();
        const lambdaFactor = io.real();
        const lambdaAdjust = io.real();
        const updateLengthTolerance = io.real();
        const maxIterations = io.integer();
        const minimizer = LevenbergMarquardtMinimizer.fromJFunction(
            problem.numP, problem.numF, problem.f, problem.j);
        emitLM(io, minimizer.minimize(problem.p0, maxIterations,
            updateLengthTolerance, 0, lambdaFactor, lambdaAdjust,
            maxAdjustments));
    }, { exact: true });

    family.case('LevenbergMarquardtMinimizer.minimize.staleResidual.deviation',
        (io) => {
            const problem = makeRosenbrock(io);
            const maxAdjustments = io.integer();
            const lambdaFactor = io.real();
            const lambdaAdjust = io.real();
            const maxIterations = io.integer();
            const minimizer = LevenbergMarquardtMinimizer.fromJFunction(
                problem.numP, problem.numF, problem.f, problem.j);
            emitLM(io, minimizer.minimize(problem.p0, maxIterations, 0, 0,
                lambdaFactor, lambdaAdjust, maxAdjustments));
        },
        { deviation: 'issue #261: LevenbergMarquardtMinimizer::DoIteration '
            + 'builds -J^T*F from the residual at the previously rejected '
            + 'candidate' });

    family.case('LevenbergMarquardtMinimizer.minimize.jPlus', (io) => {
        const problem = makeCircleFit(io);
        const maxIterations = io.integer();
        const maxAdjustments = io.integer();
        const lambdaFactor = io.real();
        const lambdaAdjust = io.real();
        const minimizer = LevenbergMarquardtMinimizer.fromJPlusFunction(
            problem.numP, problem.numF, problem.f, problem.normal);
        emitLM(io, minimizer.minimize(problem.p0, maxIterations, 0, 0,
            lambdaFactor, lambdaAdjust, maxAdjustments));
    }, { exact: true });

    family.case('LevenbergMarquardtMinimizer.minimize.rankDeficient', (io) => {
        const problem = makeLinearModel(io, io.index % 3);
        const maxAdjustments = io.integer();
        const lambdaFactor = io.real();
        const lambdaAdjust = io.real();
        const maxIterations = io.integer();
        const minimizer = LevenbergMarquardtMinimizer.fromJFunction(
            problem.numP, problem.numF, problem.f, problem.j);
        emitLM(io, minimizer.minimize(problem.p0, maxIterations, 0, 0,
            lambdaFactor, lambdaAdjust, maxAdjustments));
    }, { exact: true });

    family.case('LevenbergMarquardtMinimizer.construct.invalidDimensions',
        (io) => {
            const numP = io.integer();
            const numF = io.integer();
            const minimizer = LevenbergMarquardtMinimizer.fromJFunction(
                numP, numF, (_p, f) => { f.set(0, 0); },
                (_p, j) => { j.set(0, 0, 0); });
            io.outInt(minimizer.getNumPDimensions());
            io.outInt(minimizer.getNumFDimensions());
        }, { exact: true });

    family.finish();
});
