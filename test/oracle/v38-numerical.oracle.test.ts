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
import { LCPSolver, LCPSolverResult } from '../../src/LCPSolver.js';
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

    family.finish();
});
