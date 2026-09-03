import { describe, it, expect } from 'vitest';
import {
    SymmetricEigensolver3x3,
    NISymmetricEigensolver3x3,
    SortEigenstuff,
    type EigenTriple,
    type EigenBasis3
} from '../src/SymmetricEigensolver3x3.js';

// A simple deterministic pseudorandom generator so test runs are repeatable.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

interface Sym3 {
    a00: number; a01: number; a02: number;
    a11: number; a12: number; a22: number;
}

function randomSym3(random: () => number): Sym3 {
    return {
        a00: 4 * random() - 2,
        a01: 4 * random() - 2,
        a02: 4 * random() - 2,
        a11: 4 * random() - 2,
        a12: 4 * random() - 2,
        a22: 4 * random() - 2
    };
}

function matVec(m: Sym3, v: EigenTriple): EigenTriple {
    return [
        m.a00 * v[0] + m.a01 * v[1] + m.a02 * v[2],
        m.a01 * v[0] + m.a11 * v[1] + m.a12 * v[2],
        m.a02 * v[0] + m.a12 * v[1] + m.a22 * v[2]
    ];
}

function dot(u: EigenTriple, v: EigenTriple): number {
    return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

function det3(evecs: EigenBasis3): number {
    return evecs[0][0] * (evecs[1][1] * evecs[2][2] - evecs[1][2] * evecs[2][1])
        - evecs[0][1] * (evecs[1][0] * evecs[2][2] - evecs[1][2] * evecs[2][0])
        + evecs[0][2] * (evecs[1][0] * evecs[2][1] - evecs[1][1] * evecs[2][0]);
}

// Verify A*v = lambda*v for each pair, orthonormality of the basis and,
// when sorted, the eigenvalue ordering and right-handedness.
function checkEigensystem(m: Sym3, evals: EigenTriple, evecs: EigenBasis3,
    sortType: number, tol: number): void {
    for (let i = 0; i < 3; ++i) {
        const av = matVec(m, evecs[i]);
        for (let j = 0; j < 3; ++j) {
            expect(Math.abs(av[j] - evals[i] * evecs[i][j])).toBeLessThanOrEqual(tol);
        }
    }

    for (let i = 0; i < 3; ++i) {
        for (let j = i; j < 3; ++j) {
            const expected = (i === j ? 1 : 0);
            expect(Math.abs(dot(evecs[i], evecs[j]) - expected)).toBeLessThanOrEqual(1e-12);
        }
    }

    if (sortType === 1) {
        expect(evals[0]).toBeLessThanOrEqual(evals[1]);
        expect(evals[1]).toBeLessThanOrEqual(evals[2]);
    }
    else if (sortType === -1) {
        expect(evals[0]).toBeGreaterThanOrEqual(evals[1]);
        expect(evals[1]).toBeGreaterThanOrEqual(evals[2]);
    }

    if (sortType !== 0) {
        // {evecs[0], evecs[1], evecs[2]} must be right-handed when sorted.
        expect(Math.abs(det3(evecs) - 1)).toBeLessThanOrEqual(1e-12);
    }
}

// Build A = R * diag(d) * R^T from a rotation R specified by axis and angle.
function fromEigensystem(d: EigenTriple, axis: EigenTriple, angle: number): Sym3 {
    const len = Math.sqrt(dot(axis, axis));
    const u: EigenTriple = [axis[0] / len, axis[1] / len, axis[2] / len];
    const c = Math.cos(angle), s = Math.sin(angle), omc = 1 - c;
    // Rotation matrix R (rows).
    const r: number[][] = [
        [c + u[0] * u[0] * omc, u[0] * u[1] * omc - u[2] * s, u[0] * u[2] * omc + u[1] * s],
        [u[1] * u[0] * omc + u[2] * s, c + u[1] * u[1] * omc, u[1] * u[2] * omc - u[0] * s],
        [u[2] * u[0] * omc - u[1] * s, u[2] * u[1] * omc + u[0] * s, c + u[2] * u[2] * omc]
    ];
    // A[i][j] = sum_k R[i][k] * d[k] * R[j][k].
    const a = (i: number, j: number): number =>
        r[i][0] * d[0] * r[j][0] + r[i][1] * d[1] * r[j][1] + r[i][2] * d[2] * r[j][2];
    return { a00: a(0, 0), a01: a(0, 1), a02: a(0, 2), a11: a(1, 1), a12: a(1, 2), a22: a(2, 2) };
}

describe('SymmetricEigensolver3x3 (iterative)', () => {
    const solver = new SymmetricEigensolver3x3();

    it('solves randomized symmetric matrices for all sort types', () => {
        const random = makeRandom(3030);
        for (let trial = 0; trial < 100; ++trial) {
            const m = randomSym3(random);
            for (const sortType of [-1, 0, 1]) {
                for (const aggressive of [false, true]) {
                    const { iterations, evals, evecs } =
                        solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, aggressive, sortType);
                    expect(iterations).toBeGreaterThanOrEqual(0);
                    checkEigensystem(m, evals, evecs, sortType, 1e-12);
                }
            }
        }
    });

    it('reproduces known eigenvalues of a constructed matrix', () => {
        const d: EigenTriple = [-2, 1, 4];
        const m = fromEigensystem(d, [1, 2, 3], 0.7);
        const { evals, evecs } = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, false, +1);
        expect(evals[0]).toBeCloseTo(-2, 12);
        expect(evals[1]).toBeCloseTo(1, 12);
        expect(evals[2]).toBeCloseTo(4, 12);
        checkEigensystem(m, evals, evecs, +1, 1e-12);
    });

    it('handles a diagonal matrix', () => {
        const m: Sym3 = { a00: 3, a01: 0, a02: 0, a11: -5, a12: 0, a22: 1 };
        const { evals, evecs } = solver.solve(3, 0, 0, -5, 0, 1, false, +1);
        expect(evals[0]).toBeCloseTo(-5, 14);
        expect(evals[1]).toBeCloseTo(1, 14);
        expect(evals[2]).toBeCloseTo(3, 14);
        checkEigensystem(m, evals, evecs, +1, 1e-13);
    });

    it('handles the zero matrix', () => {
        const m: Sym3 = { a00: 0, a01: 0, a02: 0, a11: 0, a12: 0, a22: 0 };
        const { evals, evecs } = solver.solve(0, 0, 0, 0, 0, 0, false, +1);
        expect(evals).toEqual([0, 0, 0]);
        checkEigensystem(m, evals, evecs, +1, 1e-14);
    });

    it('handles repeated eigenvalues', () => {
        // A with eigenvalues {2, 2, 5}.
        const m = fromEigensystem([2, 2, 5], [1, -1, 2], 1.1);
        const { evals, evecs } = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, false, +1);
        expect(evals[0]).toBeCloseTo(2, 12);
        expect(evals[1]).toBeCloseTo(2, 12);
        expect(evals[2]).toBeCloseTo(5, 12);
        checkEigensystem(m, evals, evecs, +1, 1e-12);
    });

    it('handles a triple eigenvalue (multiple of identity)', () => {
        const m: Sym3 = { a00: 7, a01: 0, a02: 0, a11: 7, a12: 0, a22: 7 };
        const { evals, evecs } = solver.solve(7, 0, 0, 7, 0, 7, false, +1);
        expect(evals[0]).toBeCloseTo(7, 14);
        expect(evals[1]).toBeCloseTo(7, 14);
        expect(evals[2]).toBeCloseTo(7, 14);
        checkEigensystem(m, evals, evecs, +1, 1e-13);
    });
});

describe('NISymmetricEigensolver3x3 (noniterative)', () => {
    const solver = new NISymmetricEigensolver3x3();

    it('solves randomized symmetric matrices for all sort types', () => {
        const random = makeRandom(4040);
        for (let trial = 0; trial < 100; ++trial) {
            const m = randomSym3(random);
            for (const sortType of [-1, 0, 1]) {
                const { evals, evecs } =
                    solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, sortType);
                checkEigensystem(m, evals, evecs, sortType, 1e-10);
            }
        }
    });

    it('agrees with the iterative solver on random matrices', () => {
        const random = makeRandom(5050);
        const iterative = new SymmetricEigensolver3x3();
        for (let trial = 0; trial < 50; ++trial) {
            const m = randomSym3(random);
            const ni = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, +1);
            const it2 = iterative.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, false, +1);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(ni.evals[i] - it2.evals[i])).toBeLessThanOrEqual(1e-10);
            }
        }
    });

    it('reproduces known eigenvalues of a constructed matrix', () => {
        const m = fromEigensystem([-3, 0.5, 6], [2, 1, -1], -0.4);
        const { evals, evecs } = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, +1);
        expect(evals[0]).toBeCloseTo(-3, 11);
        expect(evals[1]).toBeCloseTo(0.5, 11);
        expect(evals[2]).toBeCloseTo(6, 11);
        checkEigensystem(m, evals, evecs, +1, 1e-11);
    });

    it('handles a diagonal matrix', () => {
        const m: Sym3 = { a00: 3, a01: 0, a02: 0, a11: -5, a12: 0, a22: 1 };
        const { evals, evecs } = solver.solve(3, 0, 0, -5, 0, 1, +1);
        expect(evals[0]).toBeCloseTo(-5, 14);
        expect(evals[1]).toBeCloseTo(1, 14);
        expect(evals[2]).toBeCloseTo(3, 14);
        checkEigensystem(m, evals, evecs, +1, 1e-13);
    });

    it('handles the zero matrix', () => {
        const m: Sym3 = { a00: 0, a01: 0, a02: 0, a11: 0, a12: 0, a22: 0 };
        const { evals, evecs } = solver.solve(0, 0, 0, 0, 0, 0, +1);
        expect(evals).toEqual([0, 0, 0]);
        checkEigensystem(m, evals, evecs, +1, 1e-14);
    });

    it('handles repeated eigenvalues', () => {
        // Repeated eigenvalues are computed to about sqrt(epsilon) accuracy
        // by the noniterative trigonometric formulas; the eigen equation
        // residual is of the same order.
        const m = fromEigensystem([2, 2, 5], [1, -1, 2], 1.1);
        const { evals, evecs } = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, +1);
        expect(evals[0]).toBeCloseTo(2, 7);
        expect(evals[1]).toBeCloseTo(2, 7);
        expect(evals[2]).toBeCloseTo(5, 11);
        checkEigensystem(m, evals, evecs, +1, 1e-7);
    });

    it('handles decreasing sort order', () => {
        const m = fromEigensystem([-1, 2, 9], [0, 1, 1], 0.3);
        const { evals, evecs } = solver.solve(m.a00, m.a01, m.a02, m.a11, m.a12, m.a22, -1);
        expect(evals[0]).toBeCloseTo(9, 11);
        expect(evals[1]).toBeCloseTo(2, 11);
        expect(evals[2]).toBeCloseTo(-1, 11);
        checkEigensystem(m, evals, evecs, -1, 1e-11);
    });
});

describe('SortEigenstuff', () => {
    it('sorts eigenvalues increasingly and keeps eigenvector pairing', () => {
        const evals: EigenTriple = [3, 1, 2];
        const evecs: EigenBasis3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        new SortEigenstuff().sort(+1, true, evals, evecs);
        expect(evals).toEqual([1, 2, 3]);
        // Eigenvector rows follow their eigenvalues (up to sign of row 2,
        // which may be negated to keep the basis right-handed).
        expect(evecs[0]).toEqual([0, 1, 0]);
        expect(evecs[1]).toEqual([0, 0, 1]);
        expect(Math.abs(evecs[2][0])).toBe(1);
        expect(Math.abs(det3(evecs) - 1)).toBeLessThanOrEqual(1e-15);
    });

    it('sorts eigenvalues decreasingly and preserves right-handedness', () => {
        const evals: EigenTriple = [3, 1, 2];
        const evecs: EigenBasis3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        new SortEigenstuff().sort(-1, true, evals, evecs);
        expect(evals).toEqual([3, 2, 1]);
        expect(Math.abs(det3(evecs) - 1)).toBeLessThanOrEqual(1e-15);
    });

    it('with sortType 0 only enforces right-handedness', () => {
        const evals: EigenTriple = [3, 1, 2];
        const evecs: EigenBasis3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        new SortEigenstuff().sort(0, false, evals, evecs);
        expect(evals).toEqual([3, 1, 2]);
        expect(evecs[2]).toEqual([-0, -0, -1]);
    });
});
