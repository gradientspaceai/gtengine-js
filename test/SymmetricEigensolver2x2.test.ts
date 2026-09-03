import { describe, it, expect } from 'vitest';
import { SymmetricEigensolver2x2 } from '../src/SymmetricEigensolver2x2.js';

// A simple deterministic pseudorandom generator so test runs are repeatable.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        // Numerical Recipes LCG constants.
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function checkEigenpair(a00: number, a01: number, a11: number,
    lambda: number, v: [number, number], tol: number): void {
    // A*v = lambda*v
    const r0 = a00 * v[0] + a01 * v[1] - lambda * v[0];
    const r1 = a01 * v[0] + a11 * v[1] - lambda * v[1];
    expect(Math.abs(r0)).toBeLessThanOrEqual(tol);
    expect(Math.abs(r1)).toBeLessThanOrEqual(tol);
}

describe('SymmetricEigensolver2x2', () => {
    const solver = new SymmetricEigensolver2x2();

    it('solves randomized symmetric matrices (eigen equation, orthonormality)', () => {
        const random = makeRandom(1010);
        for (let trial = 0; trial < 200; ++trial) {
            const a00 = 4 * random() - 2;
            const a01 = 4 * random() - 2;
            const a11 = 4 * random() - 2;
            const scale = Math.max(Math.abs(a00), Math.abs(a01), Math.abs(a11), 1);
            const tol = 1e-13 * scale;

            for (const sortType of [-1, 0, 1]) {
                const { evals, evecs } = solver.solve(a00, a01, a11, sortType);

                checkEigenpair(a00, a01, a11, evals[0], evecs[0], tol);
                checkEigenpair(a00, a01, a11, evals[1], evecs[1], tol);

                // Orthonormality of the eigenvectors.
                const len0 = Math.hypot(evecs[0][0], evecs[0][1]);
                const len1 = Math.hypot(evecs[1][0], evecs[1][1]);
                const dot01 = evecs[0][0] * evecs[1][0] + evecs[0][1] * evecs[1][1];
                expect(Math.abs(len0 - 1)).toBeLessThanOrEqual(1e-14);
                expect(Math.abs(len1 - 1)).toBeLessThanOrEqual(1e-14);
                expect(Math.abs(dot01)).toBeLessThanOrEqual(1e-14);

                // Right-handedness: det[evecs] = +1.
                const det = evecs[0][0] * evecs[1][1] - evecs[0][1] * evecs[1][0];
                expect(Math.abs(det - 1)).toBeLessThanOrEqual(1e-14);
            }
        }
    });

    it('orders the eigenvalues according to sortType', () => {
        const random = makeRandom(2020);
        for (let trial = 0; trial < 100; ++trial) {
            const a00 = 4 * random() - 2;
            const a01 = 4 * random() - 2;
            const a11 = 4 * random() - 2;

            const inc = solver.solve(a00, a01, a11, +1);
            expect(inc.evals[0]).toBeLessThanOrEqual(inc.evals[1]);

            const dec = solver.solve(a00, a01, a11, -1);
            expect(dec.evals[0]).toBeGreaterThanOrEqual(dec.evals[1]);
        }
    });

    it('matches the analytic eigenvalues of a known matrix', () => {
        // A = [[2, 1], [1, 2]] has eigenvalues 1 and 3.
        const { evals, evecs } = solver.solve(2, 1, 2, +1);
        expect(evals[0]).toBeCloseTo(1, 14);
        expect(evals[1]).toBeCloseTo(3, 14);
        // Eigenvector for eigenvalue 1 is parallel to (1,-1).
        expect(Math.abs(evecs[0][0] + evecs[0][1])).toBeLessThanOrEqual(1e-14);
        // Eigenvector for eigenvalue 3 is parallel to (1,1).
        expect(Math.abs(evecs[1][0] - evecs[1][1])).toBeLessThanOrEqual(1e-14);
    });

    it('handles a diagonal matrix', () => {
        const { evals, evecs } = solver.solve(5, 0, -7, +1);
        expect(evals[0]).toBeCloseTo(-7, 14);
        expect(evals[1]).toBeCloseTo(5, 14);
        checkEigenpair(5, 0, -7, evals[0], evecs[0], 1e-13);
        checkEigenpair(5, 0, -7, evals[1], evecs[1], 1e-13);
    });

    it('handles the zero matrix', () => {
        const { evals, evecs } = solver.solve(0, 0, 0, +1);
        expect(evals[0]).toBe(0);
        expect(evals[1]).toBe(0);
        const det = evecs[0][0] * evecs[1][1] - evecs[0][1] * evecs[1][0];
        expect(Math.abs(det - 1)).toBeLessThanOrEqual(1e-14);
    });

    it('handles repeated eigenvalues (multiple of identity)', () => {
        const { evals, evecs } = solver.solve(3, 0, 3, +1);
        expect(evals[0]).toBeCloseTo(3, 14);
        expect(evals[1]).toBeCloseTo(3, 14);
        const dot01 = evecs[0][0] * evecs[1][0] + evecs[0][1] * evecs[1][1];
        expect(Math.abs(dot01)).toBeLessThanOrEqual(1e-14);
    });
});
