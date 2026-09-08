import { describe, it, expect } from 'vitest';
import { Minimize1 } from '../src/Minimize1.js';
import {
    check, expectClose, fc, positive, scaled
} from './helpers/arbitraries.js';

// A brute-force scan of F on [t0,t1] used as an independent estimate of the
// global minimum.
function scanMinimum(F: (t: number) => number, t0: number, t1: number,
    samples: number): { tMin: number; fMin: number } {
    let tMin = t0;
    let fMin = F(t0);
    for (let i = 1; i <= samples; ++i) {
        const t = t0 + (t1 - t0) * (i / samples);
        const f = F(t);
        if (f < fMin) {
            fMin = f;
            tMin = t;
        }
    }
    return { tMin, fMin };
}

describe('Minimize1', () => {
    it('rejects invalid subdivision and bisection counts', () => {
        expect(() => new Minimize1((t) => t, 0, 8)).toThrow('Invalid argument.');
        expect(() => new Minimize1((t) => t, 8, 0)).toThrow('Invalid argument.');
        expect(() => new Minimize1((t) => t, -1, -1)).toThrow('Invalid argument.');
    });

    it('rejects an initial guess outside the interval', () => {
        const minimizer = new Minimize1((t) => t * t, 8, 64);
        expect(() => minimizer.getMinimum(0, 1, -0.5)).toThrow('Invalid initial t value.');
        expect(() => minimizer.getMinimum(0, 1, 1.5)).toThrow('Invalid initial t value.');
        expect(() => minimizer.getMinimum(0, 1, 0)).not.toThrow();
        expect(() => minimizer.getMinimum(0, 1, 1)).not.toThrow();
    });

    it('clamps epsilon and tolerance to nonnegative values', () => {
        const minimizer = new Minimize1((t) => t * t, 8, 64, -1, -2);
        expect(minimizer.getEpsilon()).toBe(0);
        expect(minimizer.getTolerance()).toBe(0);

        minimizer.setEpsilon(1e-10);
        minimizer.setTolerance(1e-6);
        expect(minimizer.getEpsilon()).toBe(1e-10);
        expect(minimizer.getTolerance()).toBe(1e-6);

        minimizer.setEpsilon(-1e-10);
        minimizer.setTolerance(-1e-6);
        expect(minimizer.getEpsilon()).toBe(0);
        expect(minimizer.getTolerance()).toBe(0);
    });

    it('uses the midpoint of the interval as the default initial guess', () => {
        const samples: number[] = [];
        const minimizer = new Minimize1((t) => {
            samples.push(t);
            return t * t;
        }, 8, 64);
        minimizer.getMinimum(-2, 6);
        // The first three evaluations are at t0, the initial guess and t1.
        expect(samples[0]).toBe(-2);
        expect(samples[1]).toBe(2);
        expect(samples[2]).toBe(6);
    });

    it('finds the vertex of a parabola essentially exactly', () => {
        // F(t) = (t - 1)^2 + 3. The parabolic interpolation of three samples
        // of a quadratic is exact, so the vertex is located immediately.
        const F = (t: number) => (t - 1) * (t - 1) + 3;
        const minimizer = new Minimize1(F, 8, 64, 1e-12, 1e-12);
        const { tMin, fMin } = minimizer.getMinimum(-1, 4, 0.5);
        expect(tMin).toBeCloseTo(1, 10);
        expect(fMin).toBeCloseTo(3, 10);
    });

    it('finds the minimum of a convex quartic', () => {
        // F(t) = (t - 0.3)^4 + 2 has a single (flat) minimum at t = 0.3.
        const F = (t: number) => {
            const d = t - 0.3;
            return d * d * d * d + 2;
        };
        const minimizer = new Minimize1(F, 16, 128, 1e-14, 1e-8);
        const { tMin, fMin } = minimizer.getMinimum(-2, 5);
        expect(tMin).toBeCloseTo(0.3, 4);
        expect(fMin).toBeCloseTo(2, 10);
    });

    it('finds the minimum of a convex function whose minimum is an endpoint', () => {
        // F(t) = t is monotone, so no polyline is V-shaped and the search is
        // a pure subdivision. The smallest sample is the left endpoint.
        const minimizer = new Minimize1((t) => t, 6, 32);
        const { tMin, fMin } = minimizer.getMinimum(0, 1);
        expect(tMin).toBe(0);
        expect(fMin).toBe(0);
    });

    it('handles the symmetric case where the parabola vertex is the midpoint', () => {
        // F(t) = t^4 - 3*t^2 has minima at t = +/-sqrt(1.5) with value
        // -2.25. The initial polyline on [-3,3] is symmetric about t = 0, so
        // the parabola vertex is the middle sample and the neighborhood
        // subdivision branch of the bracketed search is exercised.
        const F = (t: number) => t * t * t * t - 3 * t * t;
        const minimizer = new Minimize1(F, 16, 128, 1e-14, 1e-10);
        const { tMin, fMin } = minimizer.getMinimum(-3, 3);
        expect(Math.abs(tMin)).toBeCloseTo(Math.sqrt(1.5), 6);
        expect(fMin).toBeCloseTo(-2.25, 10);
    });

    it('finds the global minimum of a multimodal quartic', () => {
        // F(t) = t^4 + t^3 - 3*t^2 has a local minimum near t = 0.906 with
        // value about -1.045 and a global minimum near t = -1.656 with value
        // about -5.249.
        const F = (t: number) => t * t * (t * t + t - 3);
        const minimizer = new Minimize1(F, 16, 128, 1e-14, 1e-10);
        const { tMin, fMin } = minimizer.getMinimum(-3, 3);
        const scan = scanMinimum(F, -3, 3, 1000000);
        expect(fMin).toBeLessThanOrEqual(scan.fMin + 1e-9);
        expect(tMin).toBeCloseTo(scan.tMin, 4);
        expect(fMin).toBeCloseTo(F(tMin), 12);
    });

    it('finds the global minimum of a multimodal trigonometric function', () => {
        // F(t) = sin(t) + 0.1*t on [0, 4*pi] has two local minima; the
        // smaller is the one nearer the left end because of the linear term.
        const F = (t: number) => Math.sin(t) + 0.1 * t;
        const minimizer = new Minimize1(F, 16, 128, 1e-14, 1e-10);
        const { tMin, fMin } = minimizer.getMinimum(0, 4 * Math.PI);
        const scan = scanMinimum(F, 0, 4 * Math.PI, 1000000);
        expect(fMin).toBeLessThanOrEqual(scan.fMin + 1e-9);
        expect(tMin).toBeCloseTo(scan.tMin, 4);
        // The derivative cos(t) + 0.1 vanishes at the minimum.
        expect(Math.cos(tMin) + 0.1).toBeCloseTo(0, 6);
    });

    it('honors the bisection budget for a slowly converging search', () => {
        // With a single bisection allowed, the search returns the best of the
        // few samples it has taken rather than a converged minimum.
        const F = (t: number) => (t - 0.3) * (t - 0.3);
        const coarse = new Minimize1(F, 1, 1, 0, 0);
        const fine = new Minimize1(F, 16, 128, 1e-14, 1e-12);
        const coarseResult = coarse.getMinimum(-2, 5);
        const fineResult = fine.getMinimum(-2, 5);
        expect(fineResult.fMin).toBeLessThanOrEqual(coarseResult.fMin);
        expect(fineResult.tMin).toBeCloseTo(0.3, 8);
    });

    it('respects the tolerance-based convergence test', () => {
        // A larger tolerance stops the bisection sooner, so the reported
        // minimum is no better than the one from a tighter tolerance.
        const F = (t: number) => Math.cosh(t - 1.25);
        const loose = new Minimize1(F, 8, 64, 1e-2, 1e-1);
        const tight = new Minimize1(F, 8, 64, 1e-14, 1e-12);
        const looseResult = loose.getMinimum(-4, 4);
        const tightResult = tight.getMinimum(-4, 4);
        expect(tightResult.fMin).toBeLessThanOrEqual(looseResult.fMin + 1e-15);
        expect(tightResult.tMin).toBeCloseTo(1.25, 6);
        expect(tightResult.fMin).toBeCloseTo(1, 12);
    });

    it('locates a minimum at an interior sample of a piecewise function', () => {
        // A V-shaped function with a corner at t = 0.5. The parabolic fit is
        // not exact, so the search relies on the bracketing logic.
        const F = (t: number) => Math.abs(t - 0.5);
        const minimizer = new Minimize1(F, 16, 200, 1e-16, 1e-14);
        const { tMin, fMin } = minimizer.getMinimum(-1, 2);
        expect(tMin).toBeCloseTo(0.5, 6);
        expect(fMin).toBeLessThan(1e-6);
    });
});

// ---------------------------------------------------------------------------
// Verification (V39): independent review against upstream Minimize1.h.
//
// The bracket-collapse defect described in the port notes of src/Minimize1.ts
// (upstream #298) is fixed here; the first property below is its regression
// test and fails on the unfixed algorithm.
// ---------------------------------------------------------------------------

describe('Minimize1 verification', () => {
    // F(t) = (t^2 - A^2) * (1 + skew*t) + c*t on [-A, A].
    //
    // The first factor vanishes at both endpoints and is negative inside, so
    // F(-A) and F(A) are equal up to the c*t term: exactly the situation that
    // defeats the upstream bracket update, where the two endpoints of the
    // search interval describe the same configuration and their function
    // values differ only by round-off. The skew moves the minimum away from
    // the midpoint, so a search that collapses onto the midpoint reports a
    // value that is visibly too large.
    const A = Math.PI / 2;

    function nearSymmetric(skew: number, c: number) {
        return (t: number) => (t * t - A * A) * (1 + skew * t) + c * t;
    }

    // The minimum of (t^2 - A^2) * (1 + skew*t) inside (-A, A): the root of
    // the derivative 3*skew*t^2 + 2*t - skew*A^2 that lies in the interval.
    function trueMinimum(skew: number): number {
        if (skew === 0) { return 0; }
        const a = 3 * skew, b = 2, d = -skew * A * A;
        const disc = Math.sqrt(b * b - 4 * a * d);
        const r0 = (-b + disc) / (2 * a);
        const r1 = (-b - disc) / (2 * a);
        return (r0 > -A && r0 < A) ? r0 : r1;
    }

    it('finds the minimum when the two endpoint values agree to round-off '
        + '(upstream #298)', () => {
        check(fc.tuple(fc.oneof(scaled(0.1, 0.6), scaled(-0.6, -0.1)),
            fc.constantFrom(0, 1e-16, 5e-16, 1e-15, 1e-14, 1e-13, 1e-12,
                1e-11, 1e-10, 1e-9)),
        ([skew, c]) => {
            for (const sign of [1, -1]) {
                const F = nearSymmetric(skew, sign * c);
                const tTrue = trueMinimum(skew);
                const fTrue = nearSymmetric(skew, sign * c)(tTrue);
                const { tMin, fMin } = new Minimize1(F, 8, 128)
                    .getMinimum(-A, A);
                // The search converges to the tolerance-limited neighborhood
                // of the true minimizer, where F is flat to second order, so
                // the value is far more accurate than the location.
                expect(fMin).toBeLessThanOrEqual(fTrue + 1e-9);
                expectClose(tMin, tTrue, 1e-3, 1e-3);
                // The upstream collapse stops at the midpoint value.
                expect(fMin).toBeLessThan(F(0));
            }
        });
    });

    it('reports a pair (tMin, fMin) with fMin = F(tMin) inside [t0,t1]',
        () => {
            check(fc.tuple(scaled(-0.6, 0.6), scaled(-1e-9, 1e-9),
                fc.integer({ min: 1, max: 16 }),
                fc.integer({ min: 1, max: 64 })),
            ([skew, c, maxSubdivisions, maxBisections]) => {
                let last = NaN;
                const F = (t: number) => {
                    last = nearSymmetric(skew, c)(t);
                    return last;
                };
                const { tMin, fMin } = new Minimize1(F, maxSubdivisions,
                    maxBisections).getMinimum(-A, A);
                expect(tMin).toBeGreaterThanOrEqual(-A);
                expect(tMin).toBeLessThanOrEqual(A);
                // Every reported minimum is an evaluated sample, so
                // re-evaluating F at tMin reproduces fMin exactly.
                expect(F(tMin)).toBe(fMin);
                expect(last).toBe(fMin);
            });
        });

    it('never reports a worse minimum when the bisection budget grows', () => {
        // The iteration for a budget of k bisections is a prefix of the one
        // for k+1, and the recorded minimum only ever decreases along it.
        check(fc.tuple(scaled(-0.6, 0.6), scaled(-1e-10, 1e-10)),
            ([skew, c]) => {
                const F = nearSymmetric(skew, c);
                let previous = Number.POSITIVE_INFINITY;
                for (const maxBisections of [1, 2, 3, 5, 8, 13, 21, 34]) {
                    const { fMin } = new Minimize1(F, 8, maxBisections)
                        .getMinimum(-A, A);
                    expect(fMin).toBeLessThanOrEqual(previous);
                    previous = fMin;
                }
            });
    });

    it('finds the minimum of a smooth unimodal function', () => {
        check(fc.tuple(scaled(-3, 3), positive(4, 0.1), positive(2, 0.05)),
            ([c, w, v]) => {
                // F(t) = w*(t-c)^2 + v*(t-c)^4 is strictly convex with its
                // only minimum at c, which lies inside [-4, 4].
                const F = (t: number) => {
                    const d = t - c;
                    return w * d * d + v * d * d * d * d;
                };
                const { tMin, fMin } = new Minimize1(F, 16, 200, 1e-15, 1e-13)
                    .getMinimum(-4, 4);
                expect(fMin).toBeGreaterThanOrEqual(0);
                // Parabolic interpolation on a smooth convex function
                // converges superlinearly, so the value is essentially exact
                // and the location is good to the square root of that.
                expect(fMin).toBeLessThan(1e-12);
                expectClose(tMin, c, 1e-5, 1e-5);
            });
    });

    it('reports the smallest value it evaluated on a multimodal function',
        () => {
            check(fc.tuple(scaled(0.2, 1.5), scaled(-1, 1)), ([w, phase]) => {
                // Several local minima on [-6, 6]. The search is not
                // guaranteed to find the deepest one - the subdivision phase
                // only refines brackets it happens to detect - but what it
                // reports must be the smallest value it ever evaluated, and
                // it must improve on the three samples it starts from.
                const base = (t: number) => Math.sin(2 * t + phase)
                    + (w * t * t) / 10;
                const values: number[] = [];
                const F = (t: number) => {
                    const f = base(t);
                    values.push(f);
                    return f;
                };
                const { tMin, fMin } = new Minimize1(F, 12, 128, 1e-14, 1e-12)
                    .getMinimum(-6, 6);
                expect(fMin).toBe(Math.min(...values));
                expect(base(tMin)).toBe(fMin);
                expect(fMin).toBeLessThanOrEqual(
                    Math.min(base(-6), base(0), base(6)));
            });
        });

    it('returns an endpoint when the minimum is at an endpoint', () => {
        check(fc.tuple(scaled(0.2, 3), fc.boolean()), ([slope, atLeft]) => {
            // A strictly monotone function has its minimum at whichever
            // endpoint the slope points to.
            const F = (t: number) => (atLeft ? slope : -slope) * t;
            const { tMin, fMin } = new Minimize1(F, 8, 64).getMinimum(-2, 3);
            expect(tMin).toBe(atLeft ? -2 : 3);
            expect(fMin).toBe(F(tMin));
        });
    });

    it('is unaffected by the fix when the parabola vertex is resolvable',
        () => {
            // The two changes made for #298 only alter the iteration when
            // the vertex of the interpolating parabola is within a few ulps
            // of the middle sample, or when the parabola is degenerate. For
            // an asymmetric bracket neither happens, and the search takes
            // exactly the upstream path: the classical parabolic
            // interpolation of a quadratic converges in one step.
            check(fc.tuple(scaled(-1.5, 1.5), positive(4, 0.5)),
                ([c, w]) => {
                    const F = (t: number) => w * (t - c) * (t - c);
                    let evaluations = 0;
                    const counted = (t: number) => { ++evaluations; return F(t); };
                    const { tMin, fMin } = new Minimize1(counted, 8, 64,
                        1e-14, 1e-12).getMinimum(-2.25, 2.5);
                    // Parabolic interpolation is exact for a quadratic, so
                    // the very first vertex is the minimizer.
                    expectClose(tMin, c, 1e-9, 1e-9);
                    expect(fMin).toBeLessThan(1e-16);
                    expect(evaluations).toBeGreaterThanOrEqual(4);
                });
        });
});
