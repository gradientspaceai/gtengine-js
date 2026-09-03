import { describe, it, expect } from 'vitest';
import { LCPSolver, LCPSolverResult } from '../src/LCPSolver.js';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Verify the LCP conditions w = q + M*z, w >= 0, z >= 0 and w^T*z = 0.
function expectLCPSolution(n: number, q: readonly number[],
    M: readonly number[], w: readonly number[], z: readonly number[],
    tolerance = 1e-9): void {
    for (let r = 0; r < n; ++r) {
        let sum = q[r];
        for (let c = 0; c < n; ++c) {
            sum += M[c + n * r] * z[c];
        }
        expect(Math.abs(w[r] - sum)).toBeLessThan(tolerance);
        expect(w[r]).toBeGreaterThan(-tolerance);
        expect(z[r]).toBeGreaterThan(-tolerance);
        expect(Math.abs(w[r] * z[r])).toBeLessThan(tolerance);
    }
}

describe('LCPSolver trivial and degenerate cases', () => {
    it('reports the trivial solution when q >= 0', () => {
        const solver = new LCPSolver(3);
        const q = [1, 2, 0];
        const M = [2, 0, 0, 0, 3, 0, 0, 0, 4];
        const out = solver.solve(q, M);
        expect(out.success).toBe(true);
        expect(out.result).toBe(LCPSolverResult.HAS_TRIVIAL_SOLUTION);
        expect(out.w).toEqual([1, 2, 0]);
        expect(out.z).toEqual([0, 0, 0]);
        expectLCPSolution(3, q, M, out.w, out.z);
    });

    it('solves the 1x1 problem with a negative q', () => {
        // w = -1 + z, complementarity forces z = 1 and w = 0.
        const solver = new LCPSolver(1);
        // The default iteration budget of n*n is 1 here, which is not enough
        // for the two pivots this problem needs; upstream documents raising
        // the budget and solving again.
        expect(solver.solve([-1], [1]).result)
            .toBe(LCPSolverResult.FAILED_TO_CONVERGE);
        solver.setMaxIterations(10);
        const out = solver.solve([-1], [1]);
        expect(out.success).toBe(true);
        expect(out.result).toBe(LCPSolverResult.HAS_NONTRIVIAL_SOLUTION);
        expect(out.w[0]).toBeCloseTo(0, 12);
        expect(out.z[0]).toBeCloseTo(1, 12);
    });

    it('reports no solution when the driving column cannot leave', () => {
        // w = -1 - z is negative for every z >= 0.
        const solver = new LCPSolver(1);
        const out = solver.solve([-1], [-1]);
        expect(out.success).toBe(false);
        expect(out.result).toBe(LCPSolverResult.NO_SOLUTION);
        expect(out.w).toEqual([0]);
        expect(out.z).toEqual([0]);
    });

    it('reports invalid input', () => {
        const solver = new LCPSolver(3);
        expect(solver.solve([1, 2], new Array<number>(9).fill(1)).result)
            .toBe(LCPSolverResult.INVALID_INPUT);
        expect(solver.solve([1, 2, 3], new Array<number>(8).fill(1)).result)
            .toBe(LCPSolverResult.INVALID_INPUT);

        // Upstream's dynamic solver dereferences null pointers when it is
        // constructed with n <= 0; the port reports INVALID_INPUT instead.
        const empty = new LCPSolver(0);
        expect(empty.getDimension()).toBe(0);
        const out = empty.solve([], []);
        expect(out.success).toBe(false);
        expect(out.result).toBe(LCPSolverResult.INVALID_INPUT);
        expect(out.w).toEqual([]);
        expect(out.z).toEqual([]);
    });

    it('defaults the iteration limit to n*n and allows overriding it', () => {
        const solver = new LCPSolver(4);
        expect(solver.getMaxIterations()).toBe(16);
        solver.setMaxIterations(100);
        expect(solver.getMaxIterations()).toBe(100);
        // A nonpositive value restores the default.
        solver.setMaxIterations(0);
        expect(solver.getMaxIterations()).toBe(16);
        solver.setMaxIterations(-5);
        expect(solver.getMaxIterations()).toBe(16);
    });

    it('records the number of iterations used', () => {
        const solver = new LCPSolver(2);
        const q = [-1, -1];
        const M = [2, 1, 1, 2];
        const out = solver.solve(q, M);
        expect(out.success).toBe(true);
        expect(solver.getNumIterations()).toBeGreaterThan(0);
        expect(solver.getNumIterations())
            .toBeLessThanOrEqual(solver.getMaxIterations());
    });
});

describe('LCPSolver nontrivial solutions', () => {
    it('solves a 2x2 problem with a known solution', () => {
        // M = [[2, 1], [1, 2]], q = (-3, -3). The solution of the linear
        // system M*z = -q is z = (1, 1) with w = 0.
        const solver = new LCPSolver(2);
        const q = [-3, -3];
        const M = [2, 1, 1, 2];
        const out = solver.solve(q, M);
        expect(out.success).toBe(true);
        expect(out.result).toBe(LCPSolverResult.HAS_NONTRIVIAL_SOLUTION);
        expect(out.z[0]).toBeCloseTo(1, 10);
        expect(out.z[1]).toBeCloseTo(1, 10);
        expect(out.w[0]).toBeCloseTo(0, 10);
        expect(out.w[1]).toBeCloseTo(0, 10);
        expectLCPSolution(2, q, M, out.w, out.z);
    });

    it('solves a problem whose solution is active in only one component', () => {
        // The identity M with q = (-2, 5) forces z = (2, 0), w = (0, 5).
        const solver = new LCPSolver(2);
        const q = [-2, 5];
        const M = [1, 0, 0, 1];
        const out = solver.solve(q, M);
        expect(out.success).toBe(true);
        expect(out.z[0]).toBeCloseTo(2, 10);
        expect(out.z[1]).toBeCloseTo(0, 10);
        expect(out.w[0]).toBeCloseTo(0, 10);
        expect(out.w[1]).toBeCloseTo(5, 10);
        expectLCPSolution(2, q, M, out.w, out.z);
    });

    it('solves the classic 3x3 Lemke example', () => {
        // A copositive M with a mixed active set.
        const solver = new LCPSolver(3);
        const q = [-3, 6, -1];
        const M = [
            1, 0, 2,
            0, 1, 1,
            2, 1, 4
        ];
        const out = solver.solve(q, M);
        expect(out.success).toBe(true);
        expectLCPSolution(3, q, M, out.w, out.z, 1e-8);
    });

    it('reuses the solver object across calls', () => {
        const solver = new LCPSolver(2);
        const M = [2, 1, 1, 2];
        const first = solver.solve([-3, -3], M);
        const second = solver.solve([-3, -3], M);
        expect(second.w).toEqual(first.w);
        expect(second.z).toEqual(first.z);
        // A trivial problem after a nontrivial one.
        const third = solver.solve([1, 1], M);
        expect(third.result).toBe(LCPSolverResult.HAS_TRIVIAL_SOLUTION);
    });

    it('solves random positive definite problems (randomized cross-check)', () => {
        const rand = makeRandom(13579);
        for (let n = 1; n <= 6; ++n) {
            for (let trial = 0; trial < 25; ++trial) {
                // M = A^T*A + I is symmetric positive definite, so the LCP
                // has a unique solution and Lemke's algorithm converges.
                const A: number[] = [];
                for (let i = 0; i < n * n; ++i) {
                    A.push(2 * rand() - 1);
                }
                const M = new Array<number>(n * n).fill(0);
                for (let r = 0; r < n; ++r) {
                    for (let c = 0; c < n; ++c) {
                        let sum = r === c ? 1 : 0;
                        for (let k = 0; k < n; ++k) {
                            sum += A[r + n * k] * A[c + n * k];
                        }
                        M[c + n * r] = sum;
                    }
                }

                const q: number[] = [];
                for (let i = 0; i < n; ++i) {
                    q.push(4 * rand() - 3);
                }

                const solver = new LCPSolver(n);
                solver.setMaxIterations(200);
                const out = solver.solve(q, M);
                expect(out.success).toBe(true);
                expectLCPSolution(n, q, M, out.w, out.z, 1e-7);
            }
        }
    });

    it('finds the nearest point of a box (a quadratic program as an LCP)', () => {
        // Minimize |x - p|^2 subject to x >= 0. The KKT conditions are the
        // LCP w = -2*p + 2*I*x, so z = max(p, 0) componentwise.
        const p = [1.5, -2.25, 0.5, -0.125];
        const n = p.length;
        const q = p.map(value => -2 * value);
        const M = new Array<number>(n * n).fill(0);
        for (let i = 0; i < n; ++i) {
            M[i + n * i] = 2;
        }
        const solver = new LCPSolver(n);
        const out = solver.solve(q, M);
        expect(out.success).toBe(true);
        for (let i = 0; i < n; ++i) {
            expect(out.z[i]).toBeCloseTo(Math.max(p[i], 0), 10);
        }
        expectLCPSolution(n, q, M, out.w, out.z);
    });
});

describe('LCPSolver failure reporting', () => {
    it('reports a failure to converge when the iteration budget is too small', () => {
        const solver = new LCPSolver(3);
        solver.setMaxIterations(1);
        const out = solver.solve([-3, -4, -5], [2, 1, 0, 1, 2, 1, 0, 1, 2]);
        expect(out.success).toBe(false);
        expect(out.result).toBe(LCPSolverResult.FAILED_TO_CONVERGE);
    });

    it('optionally throws on failure (GTE_THROW_ON_LCPSOLVER_ERRORS)', () => {
        const solver = new LCPSolver(3);
        solver.setMaxIterations(1);
        LCPSolver.throwOnErrors = true;
        try {
            expect(() => solver.solve([-3, -4, -5], [2, 1, 0, 1, 2, 1, 0, 1, 2]))
                .toThrow(/failed to converge/);
        } finally {
            LCPSolver.throwOnErrors = false;
        }
    });
});
