import { describe, expect, it } from 'vitest';
import { MinimizeN } from '../src/MinimizeN.js';

// Powell's direction set method with maxLevel/maxBracket handed to Minimize1.
const MAX_LEVEL = 8;
const MAX_BRACKET = 64;
const MAX_ITERATIONS = 32;

function makeMinimizer(dimensions: number, F: (t: readonly number[]) => number,
    maxIterations: number = MAX_ITERATIONS, epsilon?: number): MinimizeN {
    return new MinimizeN(dimensions, F, MAX_LEVEL, MAX_BRACKET, maxIterations, epsilon);
}

describe('MinimizeN', () => {
    describe('member access', () => {
        it('clamps a negative epsilon to zero', () => {
            const minimizer = makeMinimizer(2, (t) => t[0] * t[0] + t[1] * t[1]);
            expect(minimizer.getEpsilon()).toBe(1e-06);
            minimizer.setEpsilon(1e-03);
            expect(minimizer.getEpsilon()).toBe(1e-03);
            minimizer.setEpsilon(-1);
            expect(minimizer.getEpsilon()).toBe(0);
        });

        it('accepts the epsilon passed to the constructor', () => {
            const minimizer = makeMinimizer(2, (t) => t[0] * t[0], MAX_ITERATIONS, 1e-09);
            expect(minimizer.getEpsilon()).toBe(1e-09);
        });
    });

    describe('axis-aligned quadratic bowls', () => {
        // F(t) = (t0-c0)^2 + (t1-c1)^2 + fOffset has its minimum at c. The
        // Euclidean basis directions are the principal axes, so a single
        // direction-set pass locates the minimum.
        const c = [0.3, -0.4];
        const fOffset = 1.25;
        const F = (t: readonly number[]): number => {
            const d0 = t[0] - c[0];
            const d1 = t[1] - c[1];
            return d0 * d0 + d1 * d1 + fOffset;
        };

        const starts = [
            [0, 0],
            [0.9, 0.9],
            [-0.9, -0.9],
            [-0.75, 0.8],
            [0.3, -0.4],
        ];

        for (const start of starts) {
            it(`recovers the minimum from (${start[0]}, ${start[1]})`, () => {
                const minimizer = makeMinimizer(2, F);
                const result = minimizer.getMinimum([-1, -1], [1, 1], start);
                expect(result.tMin[0]).toBeCloseTo(c[0], 6);
                expect(result.tMin[1]).toBeCloseTo(c[1], 6);
                expect(result.fMin).toBeCloseTo(fOffset, 10);
                expect(result.fMin).toBeLessThanOrEqual(F(start));
            });
        }

        it('works in three dimensions', () => {
            const center = [0.25, -0.5, 0.75];
            const G = (t: readonly number[]): number => {
                let sum = 0;
                for (let i = 0; i < 3; ++i) {
                    const d = t[i] - center[i];
                    sum += (i + 1) * d * d;
                }
                return sum;
            };
            const minimizer = makeMinimizer(3, G);
            const result = minimizer.getMinimum([-1, -1, -1], [1, 1, 1], [0, 0, 0]);
            for (let i = 0; i < 3; ++i) {
                expect(result.tMin[i]).toBeCloseTo(center[i], 6);
            }
            expect(result.fMin).toBeCloseTo(0, 10);
            expect(result.tMin).toHaveLength(3);
        });

        it('returns a copy of the location that the caller may modify', () => {
            const minimizer = makeMinimizer(2, F);
            const first = minimizer.getMinimum([-1, -1], [1, 1], [0, 0]);
            first.tMin[0] = 100;
            const second = minimizer.getMinimum([-1, -1], [1, 1], [0.5, 0.5]);
            expect(second.tMin[0]).toBeCloseTo(c[0], 6);
        });
    });

    describe('coupled quadratic (non-axis-aligned valley)', () => {
        // F(t) = 5*u^2 + v^2 with u = (t0+t1)/sqrt(2), v = (t0-t1)/sqrt(2)
        // shifted so the minimum is at (0.2, -0.1). The principal axes are
        // rotated 45 degrees from the Euclidean basis, so the conjugate
        // direction is what drives the convergence.
        const center = [0.2, -0.1];
        const F = (t: readonly number[]): number => {
            const d0 = t[0] - center[0];
            const d1 = t[1] - center[1];
            const u = (d0 + d1) / Math.SQRT2;
            const v = (d0 - d1) / Math.SQRT2;
            return 5 * u * u + v * v;
        };

        it('locates the minimum of the rotated bowl', () => {
            const minimizer = makeMinimizer(2, F);
            const result = minimizer.getMinimum([-1, -1], [1, 1], [-0.8, 0.7]);
            expect(result.tMin[0]).toBeCloseTo(center[0], 5);
            expect(result.tMin[1]).toBeCloseTo(center[1], 5);
            expect(result.fMin).toBeCloseTo(0, 8);
        });

        it('never returns a value larger than the value at the initial guess', () => {
            const starts = [[-1, -1], [1, 1], [0, 0], [0.9, -0.9], [-0.4, 0.6]];
            for (const start of starts) {
                const minimizer = makeMinimizer(2, F);
                const result = minimizer.getMinimum([-1, -1], [1, 1], start);
                expect(result.fMin).toBeLessThanOrEqual(F(start));
                // The reported value is consistent with the reported location
                // to within the tolerance of the line searches.
                expect(Math.abs(result.fMin - F(result.tMin))).toBeLessThan(1e-08);
            }
        });
    });

    describe('Rosenbrock-style valley', () => {
        // F(t) = (1-t0)^2 + 100*(t1-t0^2)^2 has its minimum 0 at (1,1). The
        // curved valley is the classic hard case for direction-set methods;
        // the test verifies substantial progress from the standard start.
        const F = (t: readonly number[]): number => {
            const a = 1 - t[0];
            const b = t[1] - t[0] * t[0];
            return a * a + 100 * b * b;
        };

        it('converges toward (1,1) from (-1.2, 1)', () => {
            const start = [-1.2, 1];
            const fStart = F(start);
            expect(fStart).toBeCloseTo(24.2, 10);

            const minimizer = makeMinimizer(2, F, 64);
            const result = minimizer.getMinimum([-2, -2], [2, 2], start);

            expect(result.fMin).toBeLessThan(1e-08);
            // The iterate lies in the valley t1 = t0^2 near (1,1).
            expect(result.tMin[1] - result.tMin[0] * result.tMin[0]).toBeCloseTo(0, 5);
            expect(result.tMin[0]).toBeCloseTo(1, 4);
            expect(result.tMin[1]).toBeCloseTo(1, 4);
            expect(result.fMin).toBe(F(result.tMin));
        });

        it('makes progress from several starting points', () => {
            for (const start of [[0, 0], [-1.5, 1.5], [1.8, -1.8], [0.5, 0.5]]) {
                const minimizer = makeMinimizer(2, F, 64);
                const result = minimizer.getMinimum([-2, -2], [2, 2], start);
                expect(result.fMin).toBeLessThanOrEqual(F(start));
                expect(result.fMin).toBeLessThan(1e-06);
                expect(result.tMin[0]).toBeCloseTo(1, 3);
                expect(result.tMin[1]).toBeCloseTo(1, 3);
            }
        });

        it('reuses a minimizer object across calls', () => {
            // The conjugate-direction slot index must be restored to a usable
            // state between calls; a second call from the same start must
            // reproduce the first result exactly.
            const minimizer = makeMinimizer(2, F, 64);
            const first = minimizer.getMinimum([-2, -2], [2, 2], [-1.2, 1]);
            minimizer.getMinimum([-2, -2], [2, 2], [0.4, -0.4]);
            const third = minimizer.getMinimum([-2, -2], [2, 2], [-1.2, 1]);
            expect(third.tMin).toEqual(first.tMin);
            expect(third.fMin).toBe(first.fMin);
        });

        it('improves monotonically as the iteration budget grows', () => {
            const start = [-1.2, 1];
            let previous = F(start);
            for (const maxIterations of [1, 2, 4, 8, 16]) {
                const minimizer = makeMinimizer(2, F, maxIterations);
                const result = minimizer.getMinimum([-2, -2], [2, 2], start);
                expect(result.fMin).toBeLessThanOrEqual(previous + 1e-12);
                previous = result.fMin;
            }
        });
    });

    describe('degenerate inputs', () => {
        it('handles a function with a flat direction', () => {
            // F ignores t1 entirely, so every line search along the t1 axis
            // finds no improvement and the t1 coordinate is left in the
            // domain. The t0 coordinate must still converge.
            const F = (t: readonly number[]): number => {
                const d = t[0] - 0.5;
                return d * d;
            };
            const minimizer = makeMinimizer(2, F);
            const result = minimizer.getMinimum([-1, -1], [1, 1], [-0.9, 0.4]);
            expect(result.tMin[0]).toBeCloseTo(0.5, 6);
            expect(result.fMin).toBeCloseTo(0, 10);
            expect(result.tMin[1]).toBeGreaterThanOrEqual(-1);
            expect(result.tMin[1]).toBeLessThanOrEqual(1);
        });

        it('reports the constant value for a constant function', () => {
            // Minimize1 keeps the first sampled point whose value is below the
            // running minimum, so for a constant function every line search
            // walks to the low end of the clipped s-interval. The iteration
            // then breaks because the conjugate direction is degenerate. The
            // value is exact and the location stays inside the domain.
            const minimizer = makeMinimizer(2, () => 7);
            const result = minimizer.getMinimum([-1, -1], [1, 1], [0.25, -0.125]);
            expect(result.fMin).toBe(7);
            expect(result.tMin[0]).toBeGreaterThanOrEqual(-1);
            expect(result.tMin[0]).toBeLessThanOrEqual(1);
            expect(result.tMin[1]).toBeGreaterThanOrEqual(-1);
            expect(result.tMin[1]).toBeLessThanOrEqual(1);
        });

        it('respects a degenerate domain where t0 equals t1', () => {
            // The clipped s-interval is [0,0] in every direction, so no motion
            // is possible and the answer is the (only) domain point.
            const F = (t: readonly number[]): number => t[0] * t[0] + t[1] * t[1];
            const minimizer = makeMinimizer(2, F);
            const result = minimizer.getMinimum([0.5, 0.5], [0.5, 0.5], [0.5, 0.5]);
            expect(result.tMin[0]).toBeCloseTo(0.5, 12);
            expect(result.tMin[1]).toBeCloseTo(0.5, 12);
            expect(result.fMin).toBeCloseTo(0.5, 12);
        });

        it('finds a minimum located on the domain boundary', () => {
            // F is monotone increasing in each variable, so the minimum of the
            // restriction to the box is the corner (t0min, t1min).
            const F = (t: readonly number[]): number => t[0] + 2 * t[1];
            const minimizer = makeMinimizer(2, F);
            const result = minimizer.getMinimum([-1, -3], [1, 3], [0, 0]);
            expect(result.tMin[0]).toBeCloseTo(-1, 10);
            expect(result.tMin[1]).toBeCloseTo(-3, 10);
            expect(result.fMin).toBeCloseTo(-7, 10);
        });

        it('never leaves the Cartesian-product domain', () => {
            // A function whose unconstrained minimum is outside the box.
            const F = (t: readonly number[]): number => {
                const d0 = t[0] - 5;
                const d1 = t[1] + 5;
                return d0 * d0 + d1 * d1;
            };
            const t0 = [-1, -2];
            const t1 = [1, 2];
            const minimizer = makeMinimizer(2, F);
            const result = minimizer.getMinimum(t0, t1, [0, 0]);
            for (let i = 0; i < 2; ++i) {
                expect(result.tMin[i]).toBeGreaterThanOrEqual(t0[i] - 1e-12);
                expect(result.tMin[i]).toBeLessThanOrEqual(t1[i] + 1e-12);
            }
            expect(result.tMin[0]).toBeCloseTo(1, 8);
            expect(result.tMin[1]).toBeCloseTo(-2, 8);
        });

        it('works with dimension 1', () => {
            const F = (t: readonly number[]): number => {
                const d = t[0] - 0.375;
                return d * d * d * d + d * d;
            };
            const minimizer = makeMinimizer(1, F);
            const result = minimizer.getMinimum([-2], [2], [1.5]);
            expect(result.tMin).toHaveLength(1);
            expect(result.tMin[0]).toBeCloseTo(0.375, 5);
            expect(result.fMin).toBeCloseTo(0, 8);
        });

        it('does nothing when maxIterations is zero', () => {
            const F = (t: readonly number[]): number => t[0] * t[0] + t[1] * t[1];
            const minimizer = makeMinimizer(2, F, 0);
            const result = minimizer.getMinimum([-1, -1], [1, 1], [0.6, -0.7]);
            expect(result.tMin[0]).toBe(0.6);
            expect(result.tMin[1]).toBe(-0.7);
            expect(result.fMin).toBe(F([0.6, -0.7]));
        });
    });

    describe('randomized cross-check', () => {
        it('matches the analytic minimizer of separable quadratics', () => {
            // F(t) = sum_i w[i]*(t[i]-c[i])^2 is separable, so the minimum on
            // the box is the componentwise clamp of c into [t0,t1].
            let seed = 987654321;
            const nextRandom = (): number => {
                // A small LCG keeps the test deterministic.
                seed = (1103515245 * seed + 12345) % 2147483648;
                return seed / 2147483648;
            };

            for (let trial = 0; trial < 20; ++trial) {
                const dimensions = 2 + (trial % 3);
                const c: number[] = [];
                const w: number[] = [];
                const lo: number[] = [];
                const hi: number[] = [];
                const start: number[] = [];
                for (let i = 0; i < dimensions; ++i) {
                    c.push(4 * nextRandom() - 2);
                    w.push(0.5 + 2 * nextRandom());
                    lo.push(-1);
                    hi.push(1);
                    start.push(2 * nextRandom() - 1);
                }
                const F = (t: readonly number[]): number => {
                    let sum = 0;
                    for (let i = 0; i < dimensions; ++i) {
                        const d = t[i] - c[i];
                        sum += w[i] * d * d;
                    }
                    return sum;
                };

                const expected = c.map((ci, i) => Math.max(lo[i], Math.min(ci, hi[i])));
                const minimizer = makeMinimizer(dimensions, F);
                const result = minimizer.getMinimum(lo, hi, start);
                for (let i = 0; i < dimensions; ++i) {
                    expect(result.tMin[i]).toBeCloseTo(expected[i], 5);
                }
                expect(result.fMin).toBeCloseTo(F(expected), 8);
            }
        });
    });
});
