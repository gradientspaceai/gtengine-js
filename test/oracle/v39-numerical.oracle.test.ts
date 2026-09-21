// Replay of the group 39 (numerical) golden records produced by the C++
// oracle. See ORACLE.md and oracle/cpp/cases/v39-numerical.cpp: each case
// body here mirrors the ORACLE_CASE of the same name, reading the inputs in
// the order the C++ side generated them and pushing the outputs in the order
// the C++ side recorded them.

import { Matrix } from '../../src/Matrix.js';
import { Vector } from '../../src/Vector.js';
import { Minimize1 } from '../../src/Minimize1.js';
import { MinimizeN } from '../../src/MinimizeN.js';
import { OdeEuler } from '../../src/OdeEuler.js';
import { OdeMidpoint } from '../../src/OdeMidpoint.js';
import { OdeRungeKutta4 } from '../../src/OdeRungeKutta4.js';
import { OdeImplicitEuler } from '../../src/OdeImplicitEuler.js';
import type { OdeSolver } from '../../src/OdeSolver.js';
import { inverse } from '../../src/Matrix.js';
import { SingularValueDecomposition } from '../../src/SingularValueDecomposition.js';
import { SymmetricEigensolver } from '../../src/SymmetricEigensolver.js';
import { RemezAlgorithm, REMEZ_FAILURE } from '../../src/RemezAlgorithm.js';
import {
    BlockCholeskyDecomposition, CholeskyDecomposition
} from '../../src/CholeskyDecomposition.js';
import {
    BlockLDLTDecomposition, LDLTDecomposition
} from '../../src/LDLTDecomposition.js';
import { OracleFamily, type OracleIO } from './harness.js';

const family = new OracleFamily('v39-numerical');

// ----------------------------------------------------------- 1-D objective
//
// The mirror of Objective1 / Evaluate1 in the case file. Every mode records
// one integer and five doubles.
interface Objective1 {
    kind: number;
    c: number[];
}

function readObjective1(io: OracleIO): Objective1 {
    const kind = io.integer();
    return { kind, c: io.reals(5) };
}

function evaluate1(f: Objective1, t: number): number {
    if (f.kind === 3) {
        const d = t - f.c[2];
        return f.c[0] + f.c[1] / (d * d + f.c[3]);
    }
    let result = f.c[4];
    result = result * t + f.c[3];
    result = result * t + f.c[2];
    result = result * t + f.c[1];
    result = result * t + f.c[0];
    return result;
}

interface Min1Draw {
    f: Objective1;
    t0: number; t1: number; tInitial: number;
    maxSubdivisions: number; maxBisections: number;
    epsilon: number; tolerance: number;
}

function readMin1(io: OracleIO): Min1Draw {
    const f = readObjective1(io);
    const t0 = io.real();
    const t1 = io.real();
    const tInitial = io.real();
    const maxSubdivisions = io.integer();
    const maxBisections = io.integer();
    const epsilon = io.real();
    const tolerance = io.real();
    return { f, t0, t1, tInitial, maxSubdivisions, maxBisections, epsilon, tolerance };
}

// ----------------------------------------------------------- N-D objective
interface ObjectiveN {
    kind: number;
    dimensions: number;
    w: number[];
    a: number[];
    g: number;
}

function evaluateN(f: ObjectiveN, x: readonly number[]): number {
    let result = 0;
    if (f.kind === 0) {
        for (let i = 0; i < f.dimensions; ++i) {
            const t = x[i] - f.a[i];
            result += f.w[i] * t * t;
        }
        for (let i = 0; i + 1 < f.dimensions; ++i) {
            result += f.g * x[i] * x[i + 1];
        }
    }
    else {
        for (let i = 0; i + 1 < f.dimensions; ++i) {
            const u = x[i + 1] - x[i] * x[i];
            const v = f.a[i] - x[i];
            result += f.w[i] * u * u + v * v;
        }
    }
    return result;
}

interface MinNDraw {
    f: ObjectiveN;
    t0: number[]; t1: number[]; tInitial: number[];
    maxLevel: number; maxBracket: number; maxIterations: number;
    epsilon: number;
}

function readMinN(io: OracleIO): MinNDraw {
    const dimensions = io.integer();
    const kind = io.integer();
    const w = io.reals(dimensions);
    const a = io.reals(dimensions);
    const g = io.real();
    const t0 = io.reals(dimensions);
    const t1 = io.reals(dimensions);
    const tInitial = io.reals(dimensions);
    const maxLevel = io.integer();
    const maxBracket = io.integer();
    const maxIterations = io.integer();
    const epsilon = io.real();
    return {
        f: { kind, dimensions, w, a, g },
        t0, t1, tInitial, maxLevel, maxBracket, maxIterations, epsilon
    };
}

// --------------------------------------------------------------- Ode system
//
// The mirror of OdeSystem<N> in the case file.
interface OdeSystem {
    kind: number;
    n: number;
    A: number[][];
    b: number[];
}

function odeF(sys: OdeSystem, t: number, x: Vector): Vector {
    const f = new Vector(sys.n);
    for (let i = 0; i < sys.n; ++i) {
        let s = sys.b[i] * t;
        for (let j = 0; j < sys.n; ++j) {
            s += (sys.kind === 0 ? sys.A[i][j] * x.get(j)
                : sys.A[i][j] * x.get(j) * x.get(j));
        }
        f.set(i, s);
    }
    return f;
}

function odeDF(sys: OdeSystem, _t: number, x: Vector): Matrix {
    const df = new Matrix(sys.n, sys.n);
    for (let r = 0; r < sys.n; ++r) {
        for (let c = 0; c < sys.n; ++c) {
            df.set(r, c, sys.kind === 0 ? sys.A[r][c] : 2 * sys.A[r][c] * x.get(c));
        }
    }
    return df;
}

function readOdeSystem(io: OracleIO, n: number, kind: number): OdeSystem {
    const A: number[][] = [];
    for (let r = 0; r < n; ++r) { A.push(io.reals(n)); }
    const b = io.reals(n);
    return { kind, n, A, b };
}

function runOde(io: OracleIO, solverId: number): void {
    const n = io.integer();
    io.integer();  // the generator mode, which only shapes the draws
    const kind = io.integer();
    const sys = readOdeSystem(io, n, kind);
    const tDelta = io.real();
    const newTDelta = io.real();
    let t = io.real();
    let x = io.vec(n);
    const numSteps = io.integer();

    const F = (tIn: number, xIn: Vector): Vector => odeF(sys, tIn, xIn);
    const DF = (tIn: number, xIn: Vector): Matrix => odeDF(sys, tIn, xIn);
    let solver: OdeSolver<Vector>;
    if (solverId === 0) { solver = new OdeEuler(tDelta, F); }
    else if (solverId === 1) { solver = new OdeMidpoint(tDelta, F); }
    else if (solverId === 2) { solver = new OdeRungeKutta4(tDelta, F); }
    else { solver = new OdeImplicitEuler(tDelta, F, DF); }

    io.outReal(solver.getTDelta());
    for (let step = 0; step < numSteps; ++step) {
        const { tOut, xOut } = solver.update(t, x);
        io.outReal(tOut);
        io.outVec(xOut);
        t = tOut;
        x = xOut;
        if (step === 0) {
            solver.setTDelta(newTDelta);
            io.outReal(solver.getTDelta());
        }
    }
}

// --------------------------------------------------------------- utilities

function readMatrix(io: OracleIO, numRows: number, numCols: number): Matrix {
    const M = new Matrix(numRows, numCols);
    for (let r = 0; r < numRows; ++r) {
        for (let c = 0; c < numCols; ++c) { M.set(r, c, io.real()); }
    }
    return M;
}

function outMatrixEntries(io: OracleIO, M: Matrix): void {
    for (let r = 0; r < M.numRows; ++r) {
        for (let c = 0; c < M.numCols; ++c) { io.outReal(M.get(r, c)); }
    }
}

// The C++ side emits -1 for the "did not converge" sentinel, which exceeds
// 2^53 in both libraries.
function outIterations(io: OracleIO, iterations: number, invalid: number): void {
    io.outInt(iterations === invalid ? -1 : iterations);
}

// =========================================================== Minimize1

family.case('Minimize1.getMinimum', (io) => {
    const d = readMin1(io);
    const minimizer = new Minimize1((t) => evaluate1(d.f, t),
        d.maxSubdivisions, d.maxBisections, d.epsilon, d.tolerance);
    const { tMin, fMin } = minimizer.getMinimum(d.t0, d.t1, d.tInitial);
    io.outReal(tMin);
    io.outReal(fMin);
    io.outReal(minimizer.getEpsilon());
    io.outReal(minimizer.getTolerance());
}, { exact: true });

family.case('Minimize1.getMinimum.defaultGuess', (io) => {
    const d = readMin1(io);
    const minimizer = new Minimize1((t) => evaluate1(d.f, t),
        d.maxSubdivisions, d.maxBisections, d.epsilon, d.tolerance);
    // The port merges upstream's two overloads with a defaulted tInitial.
    const { tMin, fMin } = minimizer.getMinimum(d.t0, d.t1);
    io.outReal(tMin);
    io.outReal(fMin);
}, { exact: true });

family.case('Minimize1.setEpsilonTolerance', (io) => {
    const d = readMin1(io);
    const newEpsilon = io.real();
    const newTolerance = io.real();
    const minimizer = new Minimize1((t) => evaluate1(d.f, t),
        d.maxSubdivisions, d.maxBisections);
    io.outReal(minimizer.getEpsilon());
    io.outReal(minimizer.getTolerance());
    const { tMin, fMin } = minimizer.getMinimum(d.t0, d.t1, d.tInitial);
    io.outReal(tMin);
    io.outReal(fMin);
    minimizer.setEpsilon(newEpsilon);
    minimizer.setTolerance(newTolerance);
    io.outReal(minimizer.getEpsilon());
    io.outReal(minimizer.getTolerance());
}, { exact: true });

family.case('Minimize1.getMinimum.deviation', (io) => {
    const d = readMin1(io);
    const minimizer = new Minimize1((t) => evaluate1(d.f, t),
        d.maxSubdivisions, d.maxBisections, d.epsilon, d.tolerance);
    const { tMin, fMin } = minimizer.getMinimum(d.t0, d.t1, d.tInitial);
    io.outReal(tMin);
    io.outReal(fMin);
}, { exact: true, deviation: 'issue #298 (Minimize1 bracket collapse)' });

family.case('Minimize1.getMinimum.endpointGuessDeviation', (io) => {
    const t0 = io.real();
    const t1 = io.real();
    const tInitial = io.real();
    const c2 = io.real();
    const maxSubdivisions = io.integer();
    const maxBisections = io.integer();
    const minimizer = new Minimize1((t) => c2 * t * t, maxSubdivisions, maxBisections);
    const { tMin, fMin } = minimizer.getMinimum(t0, t1, tInitial);
    io.outReal(tMin);
    io.outReal(fMin);
}, { exact: true, deviation: 'issue #298 (degenerate initial bracket)' });

family.case('Minimize1.invalidInput', (io) => {
    const maxSubdivisions = io.integer();
    const maxBisections = io.integer();
    const t0 = io.real();
    const t1 = io.real();
    const tInitial = io.real();
    const c2 = io.real();
    const minimizer = new Minimize1((t) => c2 * t * t, maxSubdivisions, maxBisections);
    const { tMin, fMin } = minimizer.getMinimum(t0, t1, tInitial);
    io.outReal(tMin);
    io.outReal(fMin);
}, { exact: true });

// =========================================================== MinimizeN

family.case('MinimizeN.getMinimum', (io) => {
    const d = readMinN(io);
    const minimizer = new MinimizeN(d.f.dimensions, (x) => evaluateN(d.f, x),
        d.maxLevel, d.maxBracket, d.maxIterations, d.epsilon);
    const { tMin, fMin } = minimizer.getMinimum(d.t0, d.t1, d.tInitial);
    for (const value of tMin) { io.outReal(value); }
    io.outReal(fMin);
    io.outReal(minimizer.getEpsilon());
}, { exact: true });

family.case('MinimizeN.getMinimum.deviation', (io) => {
    const d = readMinN(io);
    const minimizer = new MinimizeN(d.f.dimensions, (x) => evaluateN(d.f, x),
        d.maxLevel, d.maxBracket, d.maxIterations, d.epsilon);
    const { tMin, fMin } = minimizer.getMinimum(d.t0, d.t1, d.tInitial);
    for (const value of tMin) { io.outReal(value); }
    io.outReal(fMin);
}, { exact: true, deviation: 'issue #146 (Powell direction-set update)' });

family.case('MinimizeN.getMinimum.reuseDeviation', (io) => {
    const d = readMinN(io);
    const shift = io.real();
    const minimizer = new MinimizeN(d.f.dimensions, (x) => evaluateN(d.f, x),
        d.maxLevel, d.maxBracket, d.maxIterations, d.epsilon);
    minimizer.getMinimum(d.t0, d.t1, d.tInitial);

    const second: number[] = [];
    for (let i = 0; i < d.f.dimensions; ++i) {
        second.push(d.tInitial[i] + shift * (d.t1[i] - d.tInitial[i]));
    }
    const { tMin, fMin } = minimizer.getMinimum(d.t0, d.t1, second);
    for (const value of tMin) { io.outReal(value); }
    io.outReal(fMin);
}, { exact: true, deviation: 'issue #146 (mDConjIndex is never reset)' });

family.case('MinimizeN.setEpsilon', (io) => {
    const d = readMinN(io);
    const minimizer = new MinimizeN(d.f.dimensions, (x) => evaluateN(d.f, x),
        d.maxLevel, d.maxBracket, d.maxIterations, -1);
    io.outReal(minimizer.getEpsilon());
    minimizer.setEpsilon(d.epsilon);
    io.outReal(minimizer.getEpsilon());
    const { tMin, fMin } = minimizer.getMinimum(d.t0, d.t1, d.tInitial);
    for (const value of tMin) { io.outReal(value); }
    io.outReal(fMin);
}, { exact: true });

// ================================================================ Ode

family.case('OdeEuler.update', (io) => { runOde(io, 0); }, { exact: true });
family.case('OdeMidpoint.update', (io) => { runOde(io, 1); }, { exact: true });
family.case('OdeRungeKutta4.update', (io) => { runOde(io, 2); }, { exact: true });
family.case('OdeImplicitEuler.update', (io) => { runOde(io, 3); }, { exact: true });

family.case('OdeImplicitEuler.update.singularJacobian', (io) => {
    const a = io.real();
    const b0 = io.real();
    const b1 = io.real();
    let t = io.real();
    const x0 = io.real();
    const x1 = io.real();
    const tDelta = io.real();
    const numSteps = io.integer();

    const F = (tIn: number, xIn: Vector): Vector => {
        const f = new Vector(2);
        f.set(0, b0 * tIn);
        f.set(1, b1 * tIn + a * xIn.get(1));
        return f;
    };
    const DF = (): Matrix => {
        const df = new Matrix(2, 2);
        df.set(0, 0, 0); df.set(0, 1, 0);
        df.set(1, 0, 0); df.set(1, 1, a);
        return df;
    };

    const solver = new OdeImplicitEuler(tDelta, F, DF);
    let x = new Vector(2);
    x.set(0, x0);
    x.set(1, x1);
    for (let step = 0; step < numSteps; ++step) {
        const result = solver.update(t, x);
        io.outReal(result.tOut);
        io.outVec(result.xOut);
        t = result.tOut;
        x = result.xOut;
    }
}, { exact: true });

// The port's Matrix.inverse is the Gaussian elimination of
// GaussianElimination.h, which is also what the C++ case's Inverse(...)
// resolves to in that translation unit; the first boolean records that the
// C++ side really compared equal to an explicit GaussianElimination call.
family.case('Matrix.inverse.isGaussianElimination', (io) => {
    const M3 = readMatrix(io, 3, 3);
    const M2 = readMatrix(io, 2, 2);
    const r3 = inverse(M3);
    const r2 = inverse(M2);
    io.outBool(true);
    io.outBool(r3.invertible);
    outMatrixEntries(io, r3.inverse);
    io.outBool(r2.invertible);
    outMatrixEntries(io, r2.inverse);
}, { exact: true });

// ========================================= SingularValueDecomposition

function runSvd(io: OracleIO, readMode: boolean): void {
    if (readMode) { io.integer(); }
    const numCols = io.integer();
    const numRows = io.integer();
    const maxIterations = io.integer();
    const multiplier = io.real();
    const A = io.reals(numRows * numCols);
    const uIndex = io.integer();
    const vIndex = io.integer();

    const svd = new SingularValueDecomposition(numRows, numCols, maxIterations);
    const iterations = svd.solve(A, multiplier);
    outIterations(io, iterations, SingularValueDecomposition.invalid);
    io.outReals(svd.getSingularValues());
    io.outReals(svd.getU());
    io.outReals(svd.getV());
    io.outReals(svd.getS());
    io.outReals(svd.getUColumn(uIndex));
    io.outReals(svd.getVColumn(vIndex));
    io.outReal(svd.getSingularValue(vIndex));
}

family.case('SingularValueDecomposition.solve',
    (io) => { runSvd(io, true); }, { exact: true });
family.case('SingularValueDecomposition.solve.rankDeficient',
    (io) => { runSvd(io, false); }, { exact: true });
family.case('SingularValueDecomposition.solve.tiedSingularValues',
    (io) => { runSvd(io, false); }, { exact: true });
family.case('SingularValueDecomposition.solve.zeroMatrix',
    (io) => { runSvd(io, false); }, { exact: true });
family.case('SingularValueDecomposition.solve.extremeScale',
    (io) => { runSvd(io, false); }, { exact: true });

family.case('SingularValueDecomposition.invalidInput', (io) => {
    const numRows = io.integer();
    const numCols = io.integer();
    const maxIterations = io.integer();
    const multiplier = io.real();
    const A = io.reals(numRows * numCols);

    const svd = new SingularValueDecomposition(numRows, numCols, maxIterations);
    const iterations = svd.solve(A, multiplier);
    outIterations(io, iterations, SingularValueDecomposition.invalid);
}, { exact: true });

family.case('SingularValueDecomposition.invalidIndex', (io) => {
    const numCols = io.integer();
    const numRows = io.integer();
    const maxIterations = io.integer();
    const accessor = io.integer();
    const index = io.integer();
    const A = io.reals(numRows * numCols);

    const svd = new SingularValueDecomposition(numRows, numCols, maxIterations);
    svd.solve(A, 8);
    if (accessor === 0) { io.outReals(svd.getUColumn(index)); }
    else if (accessor === 1) { io.outReals(svd.getVColumn(index)); }
    else { io.outReal(svd.getSingularValue(index)); }
}, { exact: true });

// ============================================== SymmetricEigensolver (NxN)

function runSymmetricEigensolver(io: OracleIO): void {
    const size = io.integer();
    const A = io.reals(size * size);
    const maxIterations = io.integer();
    const sortType = io.integer();
    const index = io.integer();

    const solver = new SymmetricEigensolver(size, maxIterations);
    const iterations = solver.solve(A, sortType);
    outIterations(io, iterations, SymmetricEigensolver.noConvergence);
    io.outInt(solver.getEigenvectorMatrixType());
    io.outReals(solver.getEigenvalues());
    io.outReals(solver.getEigenvectors());
    io.outInt(solver.getEigenvectorMatrixType());
    io.outReals(solver.getEigenvector(index));
    io.outReal(solver.getEigenvalue(index));
}

family.case('SymmetricEigensolver.solve',
    (io) => { runSymmetricEigensolver(io); }, { exact: true });

family.case('SymmetricEigensolver.solve.decoupledDeviation',
    (io) => { runSymmetricEigensolver(io); },
    { exact: true, deviation: 'issue #80 (degenerate Householder step)' });

// GetEigenvectors is deliberately not called: after a Solve that did not
// converge it loops forever, on both sides. See the group report.
family.case('SymmetricEigensolver.solve.nonConvergence', (io) => {
    const size = io.integer();
    const A = io.reals(size * size);
    const maxIterations = io.integer();
    const sortType = io.integer();
    const index = io.integer();

    const solver = new SymmetricEigensolver(size, maxIterations);
    const iterations = solver.solve(A, sortType);
    outIterations(io, iterations, SymmetricEigensolver.noConvergence);
    io.outInt(solver.getEigenvectorMatrixType());
    io.outReals(solver.getEigenvalues());
    io.outReals(solver.getEigenvector(index));
    io.outReal(solver.getEigenvalue(index));
}, { exact: true });

family.case('SymmetricEigensolver.solve.invalidSize', (io) => {
    const size = io.integer();
    const maxIterations = io.integer();
    const index = io.integer();
    const A = io.reals(size * size);

    const solver = new SymmetricEigensolver(size, maxIterations);
    const iterations = solver.solve(A, 0);
    outIterations(io, iterations, SymmetricEigensolver.noConvergence);
    io.outInt(solver.getEigenvectorMatrixType());
    io.outReal(solver.getEigenvalue(index));
}, { exact: true });

// ======================================================= RemezAlgorithm
//
// The mirror of RemezFunction in the case file.
function remezF(kind: number, c: readonly number[], x: number): number {
    if (kind === 0) {
        const n = c[0] + c[1] * x;
        const d = c[2] + x * x;
        return n / d;
    }
    let result = c[5];
    result = result * x + c[4];
    result = result * x + c[3];
    result = result * x + c[2];
    result = result * x + c[1];
    result = result * x + c[0];
    return result;
}

function remezFDer(kind: number, c: readonly number[], x: number): number {
    if (kind === 0) {
        const n = c[0] + c[1] * x;
        const d = c[2] + x * x;
        return (c[1] * d - n * (2 * x)) / (d * d);
    }
    let result = 5 * c[5];
    result = result * x + 4 * c[4];
    result = result * x + 3 * c[3];
    result = result * x + 2 * c[2];
    result = result * x + c[1];
    return result;
}

function runRemez(io: OracleIO): void {
    const kind = io.integer();
    const c = io.reals(6);
    const xMin = io.real();
    const xMax = io.real();
    const degree = io.integer();
    const maxRemezIterations = io.integer();
    const maxBisectionIterations = io.integer();
    const maxBracketIterations = io.integer();

    const remez = new RemezAlgorithm();
    const iterations = remez.execute(
        (x) => remezF(kind, c, x), (x) => remezFDer(kind, c, x),
        xMin, xMax, degree, maxRemezIterations, maxBisectionIterations,
        maxBracketIterations);
    outIterations(io, iterations, REMEZ_FAILURE);
    io.outReal(remez.getEstimatedMaxError());
    io.outReals(remez.getCoefficients());
    io.outReals(remez.getXNodes());
    io.outReals(remez.getErrors());
}

family.case('RemezAlgorithm.execute.degreeOne',
    (io) => { runRemez(io); }, { exact: true });
family.case('RemezAlgorithm.execute.rational',
    (io) => { runRemez(io); }, { exact: true });
family.case('RemezAlgorithm.execute.polynomial',
    (io) => { runRemez(io); }, { exact: true });
family.case('RemezAlgorithm.execute.singleIteration',
    (io) => { runRemez(io); }, { exact: true });

family.case('RemezAlgorithm.execute.invalidInput', (io) => {
    const xMin = io.real();
    const xMax = io.real();
    const degree = io.integer();
    const maxRemezIterations = io.integer();
    const maxBisectionIterations = io.integer();
    const maxBracketIterations = io.integer();
    const c = io.reals(3);

    const remez = new RemezAlgorithm();
    const iterations = remez.execute(
        (x) => remezF(0, c, x), (x) => remezFDer(0, c, x),
        xMin, xMax, degree, maxRemezIterations, maxBisectionIterations,
        maxBracketIterations);
    outIterations(io, iterations, REMEZ_FAILURE);
    io.outReal(remez.getEstimatedMaxError());
    io.outReals(remez.getCoefficients());
}, { exact: true });

// ================================================= CholeskyDecomposition

family.case('CholeskyDecomposition.factorAndSolve', (io) => {
    const n = io.integer();
    io.integer();  // the generator mode, which only shapes the draws
    const M = readMatrix(io, n, n);
    const B = io.reals(n);

    const decomposer = new CholeskyDecomposition(n);
    const success = decomposer.factor(M);
    io.outBool(success);
    outMatrixEntries(io, M);
    if (success) {
        const Y = Vector.fromArray(B);
        decomposer.solveLower(M, Y);
        io.outVec(Y);
        decomposer.solveUpper(M, Y);
        io.outVec(Y);
    }
}, { exact: true });

family.case('CholeskyDecomposition.factorAndSolve.fixedSize', (io) => {
    io.integer();  // the generator mode
    const M = readMatrix(io, 3, 3);
    const B = io.reals(3);

    const decomposer = new CholeskyDecomposition(3);
    const success = decomposer.factor(M);
    io.outBool(success);
    outMatrixEntries(io, M);
    if (success) {
        const Y = Vector.fromArray(B);
        decomposer.solveLower(M, Y);
        io.outVec(Y);
        decomposer.solveUpper(M, Y);
        io.outVec(Y);
    }
}, { exact: true });

family.case('CholeskyDecomposition.invalidSize', (io) => {
    const n = io.integer();
    io.boolean();  // whether the generator used the valid shape
    const rows = io.integer();
    const cols = io.integer();
    const vectorSize = io.integer();
    const M = readMatrix(io, rows, cols);
    const B = io.reals(vectorSize);

    const decomposer = new CholeskyDecomposition(n);
    const success = decomposer.factor(M);
    const Y = Vector.fromArray(B);
    decomposer.solveLower(M, Y);
    const afterLower = Y.values.slice();
    decomposer.solveUpper(M, Y);

    io.outBool(success);
    io.outReals(afterLower);
    io.outVec(Y);
}, { exact: true });

// ============================================ BlockCholeskyDecomposition

// The block matrix layout both sides use: numBlocks*numBlocks blocks of
// blockSize-by-blockSize, in row-major order over the blocks.
function newBlocks(blockSize: number, numBlocks: number): Matrix[] {
    const blocks: Matrix[] = [];
    for (let i = 0; i < numBlocks * numBlocks; ++i) {
        blocks.push(new Matrix(blockSize, blockSize));
    }
    return blocks;
}

function newBlockVector(blockSize: number, numBlocks: number,
    rhs: readonly number[]): Vector[] {
    const Y: Vector[] = [];
    for (let r = 0, k = 0; r < numBlocks; ++r) {
        const block = new Vector(blockSize);
        for (let i = 0; i < blockSize; ++i, ++k) { block.set(i, rhs[k]); }
        Y.push(block);
    }
    return Y;
}

family.case('BlockCholeskyDecomposition.factorAndSolve', (io) => {
    const shapes = [[1, 2], [1, 3], [1, 4], [2, 1], [2, 2], [2, 3],
        [3, 1], [3, 2], [1, 1]];
    const shape = io.integer();
    const B = shapes[shape][0], Nb = shapes[shape][1];
    const dim = B * Nb;
    io.integer();  // the generator mode
    const A = io.reals(dim * dim);
    const rhs = io.reals(dim);

    const decomposer = new BlockCholeskyDecomposition(B, Nb);
    const M = newBlocks(B, Nb);
    for (let r = 0; r < dim; ++r) {
        for (let c = 0; c < dim; ++c) { decomposer.set(M, r, c, A[c + dim * r]); }
    }
    for (let r = 0; r < dim; ++r) {
        for (let c = 0; c < dim; ++c) { io.outReal(decomposer.get(M, r, c)); }
    }

    const success = decomposer.factor(M);
    io.outBool(success);
    for (let r = 0; r < dim; ++r) {
        for (let c = 0; c < dim; ++c) { io.outReal(decomposer.get(M, r, c)); }
    }

    if (success) {
        const Y = newBlockVector(B, Nb, rhs);
        decomposer.solveLower(M, Y);
        for (const block of Y) { io.outVec(block); }
        decomposer.solveUpper(M, Y);
        for (const block of Y) { io.outVec(block); }
    }
}, { exact: true });

family.case('BlockCholeskyDecomposition.runtimeStrideDeviation', (io) => {
    const shapes = [[3, 2], [4, 2]];
    const shape = io.integer();
    const B = shapes[shape][0], Nb = shapes[shape][1];
    const dim = B * Nb;
    io.integer();  // the generator mode
    const A = io.reals(dim * dim);
    const rhs = io.reals(dim);

    const decomposer = new BlockCholeskyDecomposition(B, Nb);
    const M = newBlocks(B, Nb);
    for (let r = 0; r < dim; ++r) {
        for (let c = 0; c < dim; ++c) { decomposer.set(M, r, c, A[c + dim * r]); }
    }

    const success = decomposer.factor(M);
    io.outBool(success);
    for (const block of M) { outMatrixEntries(io, block); }

    if (success) {
        const Y = newBlockVector(B, Nb, rhs);
        decomposer.solveLower(M, Y);
        for (const block of Y) { io.outVec(block); }
        decomposer.solveUpper(M, Y);
        for (const block of Y) { io.outVec(block); }
    }
}, { exact: true, deviation: 'issue #209 (in-block stride of the run-time class)' });

family.case('BlockCholeskyDecomposition.invalidSize', (io) => {
    const B = io.integer();
    const Nb = io.integer();
    const value = io.real();
    const decomposer = new BlockCholeskyDecomposition(B, Nb);
    io.outInt(decomposer.numDimensions);
    io.outReal(value);
}, { exact: true });

// ==================================================== LDLTDecomposition

family.case('LDLTDecomposition.factorAndSolve', (io) => {
    const n = io.integer();
    io.integer();  // the generator mode
    const M = readMatrix(io, n, n);
    const B = io.reals(n);
    const rhs = Vector.fromArray(B);

    const decomposer = new LDLTDecomposition(n);
    const { success, L, D } = decomposer.factor(M);
    io.outBool(success);
    outMatrixEntries(io, L);
    outMatrixEntries(io, D);
    if (success) {
        io.outVec(decomposer.solveFactored(L, D, rhs));
    }

    const result = decomposer.solve(M, rhs);
    io.outBool(result.success);
    if (result.success) { io.outVec(result.X); }
}, { exact: true });

family.case('LDLTDecomposition.factorAndSolve.fixedSize', (io) => {
    io.integer();  // the generator mode
    const M = readMatrix(io, 4, 4);
    const B = io.reals(4);

    const decomposer = new LDLTDecomposition(4);
    const { success, L, D } = decomposer.factor(M);
    io.outBool(success);
    outMatrixEntries(io, L);
    outMatrixEntries(io, D);
    const result = decomposer.solve(M, Vector.fromArray(B));
    io.outBool(result.success);
    if (result.success) { io.outVec(result.X); }
}, { exact: true });

family.case('LDLTDecomposition.invalidSize', (io) => {
    const n = io.integer();
    io.boolean();  // whether the generator used the valid shape
    const rows = io.integer();
    const vectorSize = io.integer();
    const M = readMatrix(io, rows, rows);
    const B = io.reals(vectorSize);

    const decomposer = new LDLTDecomposition(n);
    const result = decomposer.solve(M, Vector.fromArray(B));
    io.outBool(result.success);
    if (result.success) { io.outVec(result.X); }
}, { exact: true });

// =============================================== BlockLDLTDecomposition

function runBlockLdlt(io: OracleIO, B: number, Nb: number,
    A: readonly number[], rhs: readonly number[]): void {
    const dim = B * Nb;
    const decomposer = new BlockLDLTDecomposition(B, Nb);
    const full = new Matrix(dim, dim);
    for (let r = 0; r < dim; ++r) {
        for (let c = 0; c < dim; ++c) { full.set(r, c, A[c + dim * r]); }
    }
    const fullRhs = Vector.fromArray(rhs);

    const MBlock = decomposer.convertMatrixToBlock(full);
    const BBlock = decomposer.convertVectorToBlock(fullRhs);

    outMatrixEntries(io, decomposer.convertBlockToMatrix(MBlock));
    for (let r = 0; r < dim; ++r) {
        for (let c = 0; c < dim; ++c) { io.outReal(decomposer.get(MBlock, r, c)); }
    }

    const { success, L, D } = decomposer.factor(MBlock);
    io.outBool(success);
    for (const block of L) { outMatrixEntries(io, block); }
    for (const block of D) { outMatrixEntries(io, block); }

    const result = decomposer.solve(MBlock, BBlock);
    io.outBool(result.success);
    if (result.success) {
        for (const block of result.X) { io.outVec(block); }
    }
}

family.case('BlockLDLTDecomposition.factorAndSolve', (io) => {
    const shapes = [[1, 2], [1, 3], [2, 2], [2, 3], [3, 2], [1, 1], [2, 1], [3, 3]];
    const shape = io.integer();
    const B = shapes[shape][0], Nb = shapes[shape][1];
    const dim = B * Nb;
    io.integer();  // the generator mode
    const A = io.reals(dim * dim);
    const rhs = io.reals(dim);
    runBlockLdlt(io, B, Nb, A, rhs);
}, { exact: true });

family.case('BlockLDLTDecomposition.getSetConvert', (io) => {
    const B = io.integer();
    const Nb = io.integer();
    const dim = B * Nb;
    const A = io.reals(dim * dim);
    const v = io.reals(dim);

    const decomposer = new BlockLDLTDecomposition(B, Nb);
    const MBlock = newBlocks(B, Nb);
    for (let r = 0; r < dim; ++r) {
        for (let c = 0; c < dim; ++c) { decomposer.set(MBlock, r, c, A[c + dim * r]); }
    }
    for (let r = 0; r < dim; ++r) {
        for (let c = 0; c < dim; ++c) { io.outReal(decomposer.get(MBlock, r, c)); }
    }

    const VBlock = decomposer.convertVectorToBlock(Vector.fromArray(v));
    io.outVec(decomposer.convertBlockToVector(VBlock));
}, { exact: true });

family.case('BlockLDLTDecomposition.convertBlockToVector.deviation', (io) => {
    const B = io.integer();
    const Nb = io.integer();
    const dim = B * Nb;
    const v = io.reals(dim);

    const decomposer = new BlockLDLTDecomposition(B, Nb);
    const VBlock = decomposer.convertVectorToBlock(Vector.fromArray(v));
    const back = decomposer.convertBlockToVector(VBlock);
    io.outInt(back.size);
    io.outVec(back);
}, { exact: true, deviation: 'issue #209 (Convert checks NumBlocks, not BlockSize)' });

family.case('BlockLDLTDecomposition.factorAndSolve.fixedSize', (io) => {
    io.integer();  // the generator mode
    const A = io.reals(16);
    const rhs = io.reals(4);

    const decomposer = new BlockLDLTDecomposition(2, 2);
    const full = new Matrix(4, 4);
    for (let r = 0; r < 4; ++r) {
        for (let c = 0; c < 4; ++c) { full.set(r, c, A[c + 4 * r]); }
    }
    const MBlock = decomposer.convertMatrixToBlock(full);
    const BBlock = decomposer.convertVectorToBlock(Vector.fromArray(rhs));

    outMatrixEntries(io, decomposer.convertBlockToMatrix(MBlock));

    const { success, L, D } = decomposer.factor(MBlock);
    io.outBool(success);
    for (const block of L) { outMatrixEntries(io, block); }
    for (const block of D) { outMatrixEntries(io, block); }

    const result = decomposer.solve(MBlock, BBlock);
    io.outBool(result.success);
    if (result.success) {
        io.outVec(decomposer.convertBlockToVector(result.X));
    }
}, { exact: true });

family.case('BlockLDLTDecomposition.invalidSize', (io) => {
    const B = io.integer();
    const Nb = io.integer();
    const value = io.real();
    const decomposer = new BlockLDLTDecomposition(B, Nb);
    io.outInt(decomposer.numDimensions);
    io.outReal(value);
}, { exact: true });

family.finish();
