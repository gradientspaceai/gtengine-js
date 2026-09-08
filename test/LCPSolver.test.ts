import { describe, it, expect } from 'vitest';
import { LCPSolver, LCPSolverResult } from '../src/LCPSolver.js';
import { check, fc, scaled, seededRandom } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (group V38).
// ---------------------------------------------------------------------------
describe('LCPSolver verification', () => {
    // Row-major M times z.
    function apply(n: number, M: readonly number[], z: readonly number[]):
        number[] {
        const out = new Array<number>(n).fill(0);
        for (let r = 0; r < n; ++r) {
            let sum = 0;
            for (let c = 0; c < n; ++c) { sum += M[r * n + c] * z[c]; }
            out[r] = sum;
        }
        return out;
    }

    // The three LCP conditions, with a tolerance scaled by the data.
    function expectSolves(n: number, q: readonly number[],
        M: readonly number[], w: readonly number[], z: readonly number[],
        tol: number): void {
        const mz = apply(n, M, z);
        for (let i = 0; i < n; ++i) {
            expect(w[i]).toBeGreaterThanOrEqual(-tol);
            expect(z[i]).toBeGreaterThanOrEqual(-tol);
            expect(Math.abs(w[i] - (q[i] + mz[i]))).toBeLessThanOrEqual(tol);
            // Complementarity: at most one of w[i], z[i] is nonzero. Lemke's
            // method keeps one of each pair nonbasic, hence exactly zero.
            expect(w[i] === 0 || z[i] === 0).toBe(true);
        }
    }

    // M = A^T*A + s*I with A a random n x n matrix is symmetric positive
    // definite, hence a P-matrix, for which Lemke's method always terminates
    // with a solution.
    const pMatrixProblem = fc.integer({ min: 1, max: 5 }).chain(n =>
        fc.record({
            n: fc.constant(n),
            entries: fc.array(scaled(-2, 2),
                { minLength: n * n, maxLength: n * n }),
            shift: fc.double({ min: 0.5, max: 3, noNaN: true,
                noDefaultInfinity: true }),
            q: fc.array(scaled(-4, 4), { minLength: n, maxLength: n })
        })).map(({ n, entries, shift, q }) => {
            const M = new Array<number>(n * n).fill(0);
            for (let r = 0; r < n; ++r) {
                for (let c = 0; c < n; ++c) {
                    let sum = (r === c ? shift : 0);
                    for (let k = 0; k < n; ++k) {
                        sum += entries[k * n + r] * entries[k * n + c];
                    }
                    M[r * n + c] = sum;
                }
            }
            return { n, M, q };
        });

    it('solves every positive definite problem and the solution satisfies the '
        + 'LCP conditions', () => {
            check(pMatrixProblem, ({ n, M, q }) => {
                const solver = new LCPSolver(n);
                // The default budget of n*n is 1 for n = 1, which is one
                // pivot short of the two a nontrivial 1x1 problem needs (see
                // the 1x1 test above). Raise it so this property is about the
                // algorithm and not about the arbitrary default.
                solver.setMaxIterations(4 * n * n + 8);
                const out = solver.solve(q, M);
                expect(out.success).toBe(true);
                expect([LCPSolverResult.HAS_TRIVIAL_SOLUTION,
                    LCPSolverResult.HAS_NONTRIVIAL_SOLUTION])
                    .toContain(out.result);
                const scale = Math.max(...M.map(Math.abs),
                    ...q.map(Math.abs), 1)
                    * Math.max(...out.z.map(Math.abs), 1);
                expectSolves(n, q, M, out.w, out.z, 1e-9 * scale);
                expect(solver.getNumIterations())
                    .toBeLessThanOrEqual(solver.getMaxIterations());
            }, 150);
        });

    it('reports the trivial solution exactly when q >= 0', () => {
        check(pMatrixProblem, ({ n, M, q }) => {
            const nonNegativeQ = q.map(Math.abs);
            const out = new LCPSolver(n).solve(nonNegativeQ, M);
            expect(out.result).toBe(LCPSolverResult.HAS_TRIVIAL_SOLUTION);
            expect(out.w).toEqual(nonNegativeQ);
            expect(out.z).toEqual(new Array<number>(n).fill(0));

            // And it is not the trivial one when some q[r] is strictly
            // negative: the perturbation polynomial q[r] + t^(r+1) is then
            // negative for small t.
            const out2 = new LCPSolver(n).solve(q, M);
            const anyNegative = q.some(x => x < 0);
            expect(out2.result === LCPSolverResult.HAS_TRIVIAL_SOLUTION)
                .toBe(!anyNegative);
        }, 150);
    });

    it('rejects invalid input the way the port documents', () => {
        check(fc.tuple(fc.integer({ min: -2, max: 4 }),
            fc.integer({ min: 0, max: 5 }), fc.integer({ min: 0, max: 30 })),
            ([n, qLen, mLen]) => {
                const solver = new LCPSolver(n);
                expect(solver.getDimension()).toBe(Math.max(n, 0));
                const out = solver.solve(new Array<number>(qLen).fill(1),
                    new Array<number>(mLen).fill(1));
                const valid = n > 0 && n <= qLen && n * n <= mLen;
                if (!valid) {
                    expect(out.result).toBe(LCPSolverResult.INVALID_INPUT);
                    expect(out.success).toBe(false);
                    expect(out.w.length).toBe(Math.max(n, 0));
                    expect(out.z.length).toBe(Math.max(n, 0));
                } else {
                    expect(out.result).not.toBe(LCPSolverResult.INVALID_INPUT);
                }
            });
    });

    it('setMaxIterations falls back to n*n for a nonpositive request', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 6 }),
            fc.integer({ min: -5, max: 40 })), ([n, requested]) => {
                const solver = new LCPSolver(n);
                expect(solver.getMaxIterations()).toBe(n * n);
                solver.setMaxIterations(requested);
                expect(solver.getMaxIterations())
                    .toBe(requested > 0 ? requested : n * n);
            });
    });

    it('a solver object is reusable and its results do not alias its state',
        () => {
            check(fc.tuple(pMatrixProblem, pMatrixProblem),
                ([first, second]) => {
                    if (first.n !== second.n) { return; }
                    const solver = new LCPSolver(first.n);
                    solver.setMaxIterations(4 * first.n * first.n + 8);
                    const a = solver.solve(first.q, first.M);
                    const aw = a.w.slice();
                    const az = a.z.slice();
                    solver.solve(second.q, second.M);
                    // The first result must be untouched by the second call.
                    expect(a.w).toEqual(aw);
                    expect(a.z).toEqual(az);
                    // And re-solving the first problem reproduces it.
                    const again = solver.solve(first.q, first.M);
                    expect(again.w).toEqual(aw);
                    expect(again.z).toEqual(az);
                }, 100);
        });

    it('upstream (preserved, #431): a degenerate PSD problem reports a '
        + 'nontrivial solution that does not satisfy w = q + M*z', () => {
            // The box-quadrilateral Hessian of DistOrientedBox3Cone3 is the
            // Gram matrix of five vectors in R^3, so its 10-D LCP is singular
            // and the Lemke pivots divide by numerically zero denominators.
            // Reproduce that shape: M is the Gram matrix of ten vectors in
            // R^3, which is positive semidefinite of rank 3.
            const rnd = seededRandom(12345);
            const n = 10;
            const A: number[][] = [];
            for (let i = 0; i < n; ++i) {
                A.push([2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1]);
            }
            const M: number[] = [];
            for (let i = 0; i < n; ++i) {
                for (let j = 0; j < n; ++j) {
                    M.push(A[i][0] * A[j][0] + A[i][1] * A[j][1]
                        + A[i][2] * A[j][2]);
                }
            }
            const q: number[] = [];
            for (let i = 0; i < n; ++i) { q.push(2 * rnd() - 1); }

            const out = new LCPSolver(n).solve(q, M);
            expect(out.result).toBe(LCPSolverResult.HAS_NONTRIVIAL_SOLUTION);
            expect(out.success).toBe(true);

            // w >= 0, z >= 0 and complementarity still hold ...
            for (let i = 0; i < n; ++i) {
                expect(out.w[i]).toBeGreaterThanOrEqual(0);
                expect(out.z[i]).toBeGreaterThanOrEqual(0);
                expect(out.w[i] === 0 || out.z[i] === 0).toBe(true);
            }
            // ... but the reported w is not q + M*z: the pivots ran through
            // denominators that are zero in exact arithmetic, and z blew up
            // to about 1e15. Upstream produces the same failure; the port
            // matches it, and this pins the signature.
            const mz = apply(n, M, out.z);
            let maxResidual = 0;
            for (let i = 0; i < n; ++i) {
                maxResidual = Math.max(maxResidual,
                    Math.abs(out.w[i] - (q[i] + mz[i])));
            }
            expect(maxResidual).toBeGreaterThan(1);
            expect(Math.max(...out.z)).toBeGreaterThan(1e12);
        });

    it('upstream: degenerate PSD problems break the LCP identity for a small '
        + 'fraction of inputs', () => {
            // The same construction over many draws: measure the rate rather
            // than assert per-instance correctness, so that a regression that
            // makes it dramatically worse (or a fix that removes it) is
            // visible.
            const rnd = seededRandom(20260907);
            const n = 8;
            let solved = 0;
            let broken = 0;
            for (let trial = 0; trial < 200; ++trial) {
                const A: number[][] = [];
                for (let i = 0; i < n; ++i) {
                    A.push([2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1]);
                }
                const M: number[] = [];
                for (let i = 0; i < n; ++i) {
                    for (let j = 0; j < n; ++j) {
                        M.push(A[i][0] * A[j][0] + A[i][1] * A[j][1]
                            + A[i][2] * A[j][2]);
                    }
                }
                const q: number[] = [];
                for (let i = 0; i < n; ++i) { q.push(2 * rnd() - 1); }
                const out = new LCPSolver(n).solve(q, M);
                if (out.result !== LCPSolverResult.HAS_NONTRIVIAL_SOLUTION) {
                    continue;
                }
                ++solved;
                const mz = apply(n, M, out.z);
                let maxResidual = 0;
                for (let i = 0; i < n; ++i) {
                    maxResidual = Math.max(maxResidual,
                        Math.abs(out.w[i] - (q[i] + mz[i])));
                }
                if (maxResidual > 1e-6) { ++broken; }
            }
            expect(solved).toBeGreaterThan(20);
            // Deterministic given the seed; the assertion is a band, so an
            // unrelated change that shifts the rate slightly does not fail
            // but a qualitative change does.
            expect(broken).toBeGreaterThan(0);
            expect(broken).toBeLessThan(solved);
        }, 30000);

    it('upstream (preserved): a subnormal pivot overflows and turns an '
        + 'infeasible problem into a reported solution', () => {
            // M <= 0 entrywise with q < 0 is infeasible (w = q + M*z <= q < 0
            // for every z >= 0), so the only correct answer is NO_SOLUTION.
            // When the pivot Augmented(basic, driving) is subnormal, the
            // reciprocal 1/pivot overflows to infinity and the pivoting
            // produces an "answer" of z = infinity with w = 0, which the
            // convergence test accepts. Upstream divides identically; this is
            // the same family as GaussianElimination's #375. Preserved and
            // pinned.
            const q = [-2.1357582139478715, -0.8588343517869644];
            const M = [-5.4e-323, -9.821122365124614e-284,
                -8.4e-323, -7.600258093787983e-155];
            const solver = new LCPSolver(2);
            solver.setMaxIterations(32);
            const out = solver.solve(q, M);
            expect(out.result).toBe(LCPSolverResult.HAS_NONTRIVIAL_SOLUTION);
            expect(out.z.every(Number.isFinite)).toBe(false);
            // With no subnormal entries the same shape of problem is reported
            // correctly.
            const solver2 = new LCPSolver(2);
            solver2.setMaxIterations(32);
            expect(solver2.solve(q, [-1, -2, -3, -4]).result)
                .toBe(LCPSolverResult.NO_SOLUTION);
        });

    it('a nonpositive M with a negative q has no solution', () => {
        // With M <= 0 entrywise and z >= 0, w = q + M*z <= q < 0 for every
        // feasible z, so the LCP is infeasible and the solver must say so
        // rather than return a wrong answer.
        check(fc.integer({ min: 1, max: 5 }).chain(n => fc.record({
            n: fc.constant(n),
            // scaled() draws from a uniform grid, so no entry is subnormal:
            // see the pin below for what a subnormal pivot does here.
            M: fc.array(scaled(-4, 0), { minLength: n * n, maxLength: n * n }),
            q: fc.array(scaled(-4, -0.5), { minLength: n, maxLength: n })
        })), ({ n, M, q }) => {
            const solver = new LCPSolver(n);
            solver.setMaxIterations(4 * n * n + 8);
            const out = solver.solve(q, M);
            expect(out.result).toBe(LCPSolverResult.NO_SOLUTION);
            expect(out.success).toBe(false);
            expect(out.w).toEqual(new Array<number>(n).fill(0));
            expect(out.z).toEqual(new Array<number>(n).fill(0));
        }, 100);
    });
});
