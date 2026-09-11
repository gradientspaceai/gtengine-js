import { describe, it, expect } from 'vitest';
import { RootsBisection1 } from '../src/RootsBisection1.js';
import {
    check, expectClose, fc, positive, wellScaled
} from './helpers/arbitraries.js';

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

// Enough iterations that the bisector always terminates by exhausting the
// 53-bit mantissa rather than by the iteration cap.
const MAX_ITERATIONS = 4096;

describe('RootsBisection1 construction', () => {
    it('rejects a nonpositive maximum iteration count', () => {
        expect(() => new RootsBisection1(0)).toThrow('Invalid maximum iterations.');
        expect(() => new RootsBisection1(-3)).toThrow('Invalid maximum iterations.');
    });

    it('rejects a misordered interval', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        expect(() => bisector.find((t) => t, 1, 0)).toThrow(
            'Invalid ordering of t-interval endpoints.');
        expect(() => bisector.find((t) => t, 1, 1)).toThrow(
            'Invalid ordering of t-interval endpoints.');
        expect(() => bisector.find((t) => t, 1, 0, -1, 1)).toThrow(
            'Invalid ordering of t-interval endpoints.');
    });
});

describe('RootsBisection1 endpoint handling', () => {
    it('returns tMin with 1 iteration when F(tMin) is exactly zero', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        const { iterations, root, fAtRoot } = bisector.find((t) => t - 2, 2, 5);
        expect(iterations).toBe(1);
        expect(root).toBe(2);
        expect(fAtRoot).toBe(0);
    });

    it('returns tMax with 1 iteration when F(tMax) is exactly zero', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        const { iterations, root, fAtRoot } = bisector.find((t) => t - 5, 2, 5);
        expect(iterations).toBe(1);
        expect(root).toBe(5);
        expect(fAtRoot).toBe(0);
    });

    it('prefers tMin when both endpoints are roots', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        // F(t) = (t + 1) * (t - 1) is zero at both endpoints of [-1, 1].
        const { iterations, root } = bisector.find((t) => (t + 1) * (t - 1), -1, 1);
        expect(iterations).toBe(1);
        expect(root).toBe(-1);
    });

    it('returns 0 iterations when the endpoint signs agree', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        // F(t) = t^2 + 1 is positive on [-1, 1].
        const { iterations, root, fAtRoot } = bisector.find((t) => t * t + 1, -1, 1);
        expect(iterations).toBe(0);
        expect(root).toBe(0);
        expect(fAtRoot).toBe(0);
    });

    it('returns 0 iterations when a root exists but is not bracketed', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        // F(t) = t^2 has a root at 0 but does not change sign on [-1, 1].
        const { iterations } = bisector.find((t) => t * t, -1, 1);
        expect(iterations).toBe(0);
    });
});

describe('RootsBisection1 bisection', () => {
    it('finds sqrt(2) to full double precision', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        const { iterations, root, fAtRoot } = bisector.find((t) => t * t - 2, 0, 2);
        expect(iterations).toBeGreaterThan(2);
        expect(iterations).toBeLessThanOrEqual(MAX_ITERATIONS);
        // Bisection on [0,2] resolves the root to the last bit.
        expect(Math.abs(root - Math.SQRT2)).toBeLessThanOrEqual(
            2 * Number.EPSILON * Math.SQRT2);
        expect(Math.abs(fAtRoot)).toBeLessThanOrEqual(1e-15);
    });

    it('stops immediately when the midpoint is an exact root', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        // The first midpoint of [-1, 3] is 1, where F(1) = 0 exactly.
        const { iterations, root, fAtRoot } = bisector.find((t) => t - 1, -1, 3);
        expect(iterations).toBe(2);
        expect(root).toBe(1);
        expect(fAtRoot).toBe(0);
    });

    it('handles a decreasing function (negative endpoint sign at tMax)', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        const { root } = bisector.find((t) => 2 - t * t, 0, 2);
        expect(Math.abs(root - Math.SQRT2)).toBeLessThanOrEqual(
            2 * Number.EPSILON * Math.SQRT2);
    });

    it('finds a root on an interval with negative endpoints', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        const { root } = bisector.find((t) => t * t * t + 8, -5, 0);
        expect(Math.abs(root + 2)).toBeLessThanOrEqual(4 * Number.EPSILON * 2);
    });

    it('honors a small iteration budget', () => {
        // Ten bisections of [0,2] narrow the bracket to 2/2^9, and the
        // midpoint is within half of that of the true root.
        const bisector = new RootsBisection1(10);
        const { iterations, root } = bisector.find((t) => t * t - 2, 0, 2);
        expect(iterations).toBe(11);
        expect(Math.abs(root - Math.SQRT2)).toBeLessThanOrEqual(2 / Math.pow(2, 9));
        expect(Math.abs(root - Math.SQRT2)).toBeGreaterThan(0);
    });

    it('returns zero-valued outputs for a one-iteration budget', () => {
        // Upstream leaves tRoot/fAtTRoot unassigned in this case (the loop
        // body never executes); the port initializes them to zero.
        const bisector = new RootsBisection1(1);
        const { iterations, root, fAtRoot } = bisector.find((t) => t * t - 2, 0, 2);
        expect(iterations).toBe(2);
        expect(root).toBe(0);
        expect(fAtRoot).toBe(0);
    });

    it('uses caller-supplied endpoint values, including infinite ones', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        // F(t) = 1/(t - 1) - 1 has a pole at t = 1 and a root at t = 2. On
        // [1, 4] the value at t = 1 is +infinity; pass its sign instead.
        const F = (t: number) => 1 / (t - 1) - 1;
        const { iterations, root } = bisector.find(F, 1, 4, +1, F(4));
        expect(iterations).toBeGreaterThan(2);
        expect(Math.abs(root - 2)).toBeLessThanOrEqual(1e-14);
    });

    it('agrees with the no-value overload when the values are supplied', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        const F = (t: number) => Math.cos(t) - t;
        const a = bisector.find(F, 0, 1);
        const b = bisector.find(F, 0, 1, F(0), F(1));
        expect(b.iterations).toBe(a.iterations);
        expect(b.root).toBe(a.root);
        expect(b.fAtRoot).toBe(a.fAtRoot);
    });

    it('reports zero-valued outputs when a supplied endpoint value is zero', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        const { iterations, root } = bisector.find((t) => t, -1, 1, 0, 1);
        expect(iterations).toBe(1);
        expect(root).toBe(-1);
    });

    it('brackets randomized roots of a cubic to near machine precision', () => {
        const rng = makeRng(0x815ec);
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        for (let trial = 0; trial < 100; ++trial) {
            // A strictly increasing cubic F(t) = (t - r)^3 + (t - r) has the
            // single root r.
            const r = 8 * rng() - 4;
            const F = (t: number) => {
                const u = t - r;
                return u * u * u + u;
            };
            const { iterations, root, fAtRoot } = bisector.find(F, r - 3, r + 5);
            expect(iterations).toBeGreaterThanOrEqual(2);
            expect(Math.abs(root - r)).toBeLessThanOrEqual(1e-14 * (1 + Math.abs(r)));
            expect(Math.abs(fAtRoot)).toBeLessThanOrEqual(1e-13);
        }
    });

    it('converges when the root sits at a tiny distance from tMin', () => {
        const bisector = new RootsBisection1(MAX_ITERATIONS);
        const root0 = 1 + 4 * Number.EPSILON;
        const { iterations, root } = bisector.find((t) => t - root0, 1, 3);
        expect(iterations).toBeGreaterThan(2);
        expect(Math.abs(root - root0)).toBeLessThanOrEqual(Number.EPSILON);
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V44): property-based comparison against upstream
// RootsBisection1.h. Only the floating-point instantiation is ported; the
// arbitrary-precision constructor (precision + maxIterations, with Convert
// rounding in RoundInitial/RoundAverage) is deliberately not present, so no
// property here exercises it.
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

describe('RootsBisection1 verification', () => {
    it('rejects a non-positive iteration budget and an unordered interval', () => {
        check(fc.tuple(fc.integer({ min: -4, max: 0 }), wellScaled(-5, 5)),
            ([maxIterations, t]) => {
                expect(() => new RootsBisection1(maxIterations)).toThrow();
                const bisector = new RootsBisection1(16);
                expect(() => bisector.find(Math.cos, t, t)).toThrow();
                expect(() => bisector.find(Math.cos, t + 1, t)).toThrow();
            }, 50);
    });

    it('reports 0 with zeroed outputs when the endpoint signs agree', () => {
        check(fc.tuple(wellScaled(-5, 5), positive(4, 0.5), positive(3, 0.1),
            positive(3, 0.1)),
            ([r, c, gap, width]) => {
                const F = monotoneF(r, c);
                const out = new RootsBisection1(64).find(F, r + gap, r + gap + width);
                expect(out.iterations).toBe(0);
                // Upstream writes tRoot = fAtTRoot = 0 on this path.
                expect(out.root === 0).toBe(true);
                expect(out.fAtRoot === 0).toBe(true);
            });
    });

    it('reports 1 with fAtRoot = 0 for an exact endpoint root', () => {
        check(fc.tuple(wellScaled(-5, 5), positive(4, 0.5), positive(3, 0.1)),
            ([r, c, width]) => {
                const F = monotoneF(r, c);
                const bisector = new RootsBisection1(64);
                const low = bisector.find(F, r, r + width);
                expect(low.iterations).toBe(1);
                expect(low.root === r).toBe(true);
                // Upstream assigns the literal zero, not F(tRoot).
                expect(low.fAtRoot === 0).toBe(true);
                const high = bisector.find(F, r - width, r);
                expect(high.iterations).toBe(1);
                expect(high.root === r).toBe(true);
                expect(high.fAtRoot === 0).toBe(true);
            });
    });

    it('converges with fAtRoot equal to F(root)', () => {
        check(bracket, ([r, c, below, above]) => {
            const F = monotoneF(r, c);
            const out = new RootsBisection1(4096).find(F, r - below, r + above);
            expect(out.iterations).toBeGreaterThanOrEqual(2);
            expect(out.iterations).toBeLessThanOrEqual(4097);
            expectClose(out.root, r, 1e-14, 1e-14);
            // Unlike the endpoint paths, this path reports the true value.
            expect(out.fAtRoot === F(out.root)).toBe(true);
        });
    });

    it('accepts precomputed endpoint values and only uses their signs', () => {
        check(bracket, ([r, c, below, above]) => {
            const F = monotoneF(r, c);
            const t0 = r - below, t1 = r + above;
            const bisector = new RootsBisection1(4096);
            const plain = bisector.find(F, t0, t1);
            const withValues = bisector.find(F, t0, t1, F(t0), F(t1));
            expect(withValues.iterations).toBe(plain.iterations);
            expect(withValues.root === plain.root).toBe(true);
            // The documented use: pass sign(f) when |f| is infinite.
            const withSigns = bisector.find(F, t0, t1, -1, +1);
            expect(withSigns.iterations).toBe(plain.iterations);
            expect(withSigns.root === plain.root).toBe(true);
            const withInfinities = bisector.find(F, t0, t1, -Infinity, Infinity);
            expect(withInfinities.root === plain.root).toBe(true);
        });
    });

    it('respects the iteration budget and halves the bracket each step', () => {
        check(fc.tuple(bracket, fc.integer({ min: 2, max: 20 })),
            ([[r, c, below, above], maxIterations]) => {
                const F = monotoneF(r, c);
                const t0 = r - below, t1 = r + above;
                const out = new RootsBisection1(maxIterations).find(F, t0, t1);
                expect(out.iterations).toBeLessThanOrEqual(maxIterations + 1);
                const k = Math.min(out.iterations, maxIterations) - 1;
                expect(Math.abs(out.root - r))
                    .toBeLessThanOrEqual((t1 - t0) / Math.pow(2, k) + 1e-15);
            });
    });

    it('returns the documented zeros when the budget allows no bisection step', () => {
        // Upstream leaves tRoot and fAtTRoot unwritten when maxIterations is
        // 1 and the endpoints bracket a root (gtengine-js issue #84); the
        // port returns zeros instead of uninitialized memory, and still
        // reports iteration = 2 as upstream's loop counter does.
        check(bracket, ([r, c, below, above]) => {
            const F = monotoneF(r, c);
            const out = new RootsBisection1(1).find(F, r - below, r + above);
            expect(out.iterations).toBe(2);
            expect(out.root === 0).toBe(true);
            expect(out.fAtRoot === 0).toBe(true);
        });
    });

    it('is reusable: a second call does not inherit state from the first', () => {
        check(fc.tuple(bracket, bracket), ([a, b]) => {
            const bisector = new RootsBisection1(4096);
            const F0 = monotoneF(a[0], a[1]);
            const F1 = monotoneF(b[0], b[1]);
            const first = bisector.find(F0, a[0] - a[2], a[0] + a[3]);
            const second = bisector.find(F1, b[0] - b[2], b[0] + b[3]);
            const fresh = new RootsBisection1(4096)
                .find(F1, b[0] - b[2], b[0] + b[3]);
            expect(second.iterations).toBe(fresh.iterations);
            expect(second.root === fresh.root).toBe(true);
            expectClose(first.root, a[0], 1e-14, 1e-14);
        }, 100);
    });
});
