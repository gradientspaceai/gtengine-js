import { describe, it, expect } from 'vitest';
import { RootsBisection } from '../src/RootsBisection.js';
import {
    check, expectClose, fc, positive, wellScaled
} from './helpers/arbitraries.js';

describe('RootsBisection', () => {
    it('finds the root of cos on [0,2] to machine precision', () => {
        const { iterations, root } = RootsBisection.find(Math.cos, 0, 2, 2048);
        expect(iterations).toBeGreaterThanOrEqual(2);
        expect(iterations).toBeLessThanOrEqual(2048);
        expect(Math.abs(root - Math.PI / 2)).toBeLessThanOrEqual(Number.EPSILON * 2);
    });

    it('converges to the floating-point neighbor of 1/3 with large iteration budget', () => {
        const F = (t: number) => t - 1 / 3;
        const { iterations, root } = RootsBisection.find(F, 0, 1, 200);
        // Bisection halts when the midpoint rounds to an interval endpoint
        // or F evaluates to exactly zero.
        expect(iterations).toBeLessThanOrEqual(60);
        expect(Math.abs(root - 1 / 3)).toBeLessThanOrEqual(Number.EPSILON);
    });

    it('returns 1 when the root is at the left bracket endpoint', () => {
        const F = (t: number) => t;
        const { iterations, root } = RootsBisection.find(F, 0, 1, 100);
        expect(iterations).toBe(1);
        expect(root).toBe(0);
    });

    it('returns 1 when the root is at the right bracket endpoint', () => {
        const F = (t: number) => t - 1;
        const { iterations, root } = RootsBisection.find(F, -3, 1, 100);
        expect(iterations).toBe(1);
        expect(root).toBe(1);
    });

    it('returns 0 when there is no sign change', () => {
        const F = (t: number) => t * t + 1;
        const { iterations } = RootsBisection.find(F, -1, 1, 100);
        expect(iterations).toBe(0);
    });

    it('returns 0 for an invalid interval (t0 >= t1)', () => {
        const F = (t: number) => t;
        expect(RootsBisection.find(F, 1, 1, 100).iterations).toBe(0);
        expect(RootsBisection.find(F, 2, 1, 100).iterations).toBe(0);
    });

    it('breaks early when the midpoint evaluates to exactly zero', () => {
        // F(t) = t - 0.25 on [0,1]: midpoints 0.5, then 0.25 where F = 0.
        const F = (t: number) => t - 0.25;
        const { iterations, root } = RootsBisection.find(F, 0, 1, 100);
        expect(root).toBe(0.25);
        expect(iterations).toBe(3);
    });

    it('returns maxIterations+1 when the budget is exhausted', () => {
        const F = (t: number) => t - 1 / 3;
        const { iterations, root } = RootsBisection.find(F, 0, 1, 5);
        expect(iterations).toBe(6);
        // The root estimate is still the best midpoint so far.
        expect(Math.abs(root - 1 / 3)).toBeLessThanOrEqual(1 / 32);
    });

    it('handles a flat region around the root (cubic contact)', () => {
        const F = (t: number) => Math.pow(t - 0.5, 3);
        const { root } = RootsBisection.find(F, 0.1, 1.3, 2048);
        expect(Math.abs(root - 0.5)).toBeLessThanOrEqual(1e-15);
    });

    it('handles a discontinuous sign function', () => {
        const x0 = 0.7234;
        const F = (t: number) => Math.sign(t - x0);
        const { root } = RootsBisection.find(F, 0, 1, 2048);
        expect(Math.abs(root - x0)).toBeLessThanOrEqual(Number.EPSILON);
    });

    describe('overload with known endpoint values', () => {
        it('accepts sign surrogates for infinite endpoint values', () => {
            // F(t) = 1/(t - 0.5) - 1/(t + 0.5) - like function with a root;
            // use F(t) = tan(t) on an interval straddling pi/2 asymmetric?
            // Instead: F(t) = 1/t - 4 on (0, 1]; F(0+) = +infinity. Pass
            // sign surrogates f0 = 1, f1 = F(1) = -3.
            const F = (t: number) => 1 / t - 4;
            const { iterations, root } = RootsBisection.find(F, 0, 1, 1, -3, 2048);
            expect(iterations).toBeGreaterThanOrEqual(2);
            expect(Math.abs(root - 0.25)).toBeLessThanOrEqual(1e-15);
        });

        it('returns 1 immediately when a supplied endpoint value is zero', () => {
            const F = (t: number) => t - 2;
            expect(RootsBisection.find(F, 2, 3, 0, 1, 100)).toEqual({ iterations: 1, root: 2 });
            expect(RootsBisection.find(F, 1, 2, -1, 0, 100)).toEqual({ iterations: 1, root: 2 });
        });

        it('returns 0 when supplied endpoint values have the same sign', () => {
            const F = (t: number) => t * t + 1;
            expect(RootsBisection.find(F, -1, 1, 2, 2, 100).iterations).toBe(0);
        });

        it('returns 0 for an invalid interval', () => {
            const F = (t: number) => t;
            expect(RootsBisection.find(F, 1, 0, -1, 1, 100).iterations).toBe(0);
        });

        it('matches the two-argument overload on a regular bracket', () => {
            const F = (t: number) => Math.exp(t) - 2;
            const a = RootsBisection.find(F, 0, 1, 2048);
            const b = RootsBisection.find(F, 0, 1, F(0), F(1), 2048);
            expect(b.root).toBe(a.root);
            expect(b.iterations).toBe(a.iterations);
            expect(Math.abs(a.root - Math.LN2)).toBeLessThanOrEqual(Number.EPSILON * 2);
        });
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V44): property-based comparison against upstream
// RootsBisection.h.
// ---------------------------------------------------------------------------

// A strictly increasing function whose only real root is r.
function monotoneF(r: number, c: number): (t: number) => number {
    return (t: number) => {
        const d = t - r;
        return d * d * d + c * d;
    };
}

const bracket = fc.tuple(wellScaled(-5, 5), positive(4, 0.5), positive(3, 0.1),
    positive(3, 0.1));

describe('RootsBisection verification', () => {
    it('rejects an interval whose endpoints are not ordered', () => {
        check(fc.tuple(wellScaled(-5, 5), wellScaled(-5, 5)), ([a, b]) => {
            const t0 = Math.max(a, b);
            const t1 = Math.min(a, b);
            const out = RootsBisection.find(Math.cos, t0, t1, 64);
            expect(out.iterations).toBe(0);
            // Upstream initializes root to t0 before the ordering test.
            expect(out.root === t0).toBe(true);
        });
    });

    it('reports failure when the endpoint values have the same sign', () => {
        check(fc.tuple(wellScaled(-5, 5), positive(4, 0.5), positive(3, 0.1),
            positive(3, 0.1)),
            ([r, c, gap, width]) => {
                const F = monotoneF(r, c);
                // Both endpoints lie above the root, so F > 0 at both.
                const t0 = r + gap;
                const out = RootsBisection.find(F, t0, t0 + width, 64);
                expect(out.iterations).toBe(0);
                expect(out.root === t0).toBe(true);
            });
    });

    it('detects an exact root at either endpoint in one iteration', () => {
        check(fc.tuple(wellScaled(-5, 5), positive(4, 0.5), positive(3, 0.1)),
            ([r, c, width]) => {
                const F = monotoneF(r, c);
                const low = RootsBisection.find(F, r, r + width, 64);
                expect(low.iterations).toBe(1);
                expect(low.root === r).toBe(true);
                const high = RootsBisection.find(F, r - width, r, 64);
                expect(high.iterations).toBe(1);
                expect(high.root === r).toBe(true);
            });
    });

    it('converges to the root with a large iteration budget', () => {
        check(bracket, ([r, c, below, above]) => {
            const F = monotoneF(r, c);
            const t0 = r - below, t1 = r + above;
            const out = RootsBisection.find(F, t0, t1, 4096);
            expect(out.iterations).toBeGreaterThanOrEqual(2);
            expect(out.iterations).toBeLessThanOrEqual(4097);
            // Full bisection to consecutive doubles, so the root is exact to
            // a couple of ulps of the interval scale.
            expectClose(out.root, r, 1e-14, 1e-14);
        });
    });

    it('halves the bracket once per iteration when the budget is small', () => {
        check(fc.tuple(bracket, fc.integer({ min: 2, max: 20 })),
            ([[r, c, below, above], maxIterations]) => {
                const F = monotoneF(r, c);
                const t0 = r - below, t1 = r + above;
                const out = RootsBisection.find(F, t0, t1, maxIterations);
                expect(out.iterations).toBeLessThanOrEqual(maxIterations + 1);
                // After k executed bisections the bracket has width
                // (t1 - t0) / 2^k and the reported midpoint is inside it.
                const k = Math.min(out.iterations, maxIterations) - 1;
                expect(Math.abs(out.root - r))
                    .toBeLessThanOrEqual((t1 - t0) / Math.pow(2, k) + 1e-15);
            });
    });

    it('gives the same answer from precomputed endpoint values and from signs', () => {
        check(bracket, ([r, c, below, above]) => {
            const F = monotoneF(r, c);
            const t0 = r - below, t1 = r + above;
            const plain = RootsBisection.find(F, t0, t1, 4096);
            const withValues = RootsBisection.find(F, t0, t1, F(t0), F(t1), 4096);
            expect(withValues.iterations).toBe(plain.iterations);
            expect(withValues.root === plain.root).toBe(true);
            // Upstream documents that only the signs of f0 and f1 matter, so
            // passing +-1 (or an infinity's sign) must not change the result.
            const withSigns = RootsBisection.find(F, t0, t1, -1, +1, 4096);
            expect(withSigns.iterations).toBe(plain.iterations);
            expect(withSigns.root === plain.root).toBe(true);
        });
    });

    it('brackets a root of any sign-changing function it accepts', () => {
        check(fc.tuple(wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
            wellScaled(-3, 3)),
            ([a, b, c, d]) => {
                const t0 = Math.min(a, b), t1 = Math.max(a, b);
                const F = (t: number) => (t - c) * (t - d) * (t * t + 1);
                const out = RootsBisection.find(F, t0, t1, 4096);
                if (out.iterations === 0) { return; }
                expect(out.root).toBeGreaterThanOrEqual(t0);
                expect(out.root).toBeLessThanOrEqual(t1);
                // A returned root is either c or d (the only real roots).
                const err = Math.min(Math.abs(out.root - c), Math.abs(out.root - d));
                expect(err).toBeLessThanOrEqual(1e-12 * (1 + Math.abs(out.root)));
            });
    });
});
