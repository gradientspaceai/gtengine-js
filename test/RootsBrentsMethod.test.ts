import { describe, it, expect } from 'vitest';
import { RootsBisection } from '../src/RootsBisection.js';
import { RootsBrentsMethod } from '../src/RootsBrentsMethod.js';
import {
    check, expectClose, fc, positive, wellScaled
} from './helpers/arbitraries.js';

// Zero tolerances: iterate until F(t) is exactly zero or the bracket
// collapses to consecutive floating-point numbers.
const tight = { maxIter: 200, negF: 0, posF: 0, stepT: 0, convT: 0 } as const;

function findTight(F: (t: number) => number, t0: number, t1: number) {
    return RootsBrentsMethod.find(F, t0, t1, tight.maxIter, tight.negF,
        tight.posF, tight.stepT, tight.convT);
}

describe('RootsBrentsMethod', () => {
    it('finds the root of cos on [0,3] to machine precision', () => {
        const { found, root } = findTight(Math.cos, 0, 3);
        expect(found).toBe(true);
        expect(Math.abs(root - Math.PI / 2)).toBeLessThanOrEqual(Number.EPSILON * 4);
    });

    it('finds the root of a cubic polynomial with tight tolerances', () => {
        // p(t) = (t - 1/3)(t^2 + 1), single real root 1/3.
        const F = (t: number) => (t - 1 / 3) * (t * t + 1);
        const { found, root } = findTight(F, -2, 2);
        expect(found).toBe(true);
        expect(Math.abs(root - 1 / 3)).toBeLessThanOrEqual(Number.EPSILON);
    });

    it('accepts a bracket endpoint that satisfies the function tolerance', () => {
        const F = (t: number) => t * (t - 2);
        // Root exactly at t0 = 0.
        let r = findTight(F, 0, 1);
        expect(r).toEqual({ found: true, root: 0 });
        // Root exactly at t1 = 2.
        r = findTight(F, 1, 2);
        expect(r).toEqual({ found: true, root: 2 });
    });

    it('uses the f-tolerances to accept near-root endpoints', () => {
        const F = (t: number) => t - 1;
        // F(0.9999) = -1e-4 is within [-1e-3, 1e-3].
        const r = RootsBrentsMethod.find(F, 0.9999, 2, 100, -1e-3, 1e-3, 0, 0);
        expect(r.found).toBe(true);
        expect(r.root).toBe(0.9999);
    });

    it('rejects an interval without a sign change', () => {
        const F = (t: number) => t * t + 1;
        expect(findTight(F, -1, 1).found).toBe(false);
    });

    it('rejects invalid parameters', () => {
        const F = (t: number) => t;
        // t1 <= t0.
        expect(RootsBrentsMethod.find(F, 1, 1, 100, 0, 0, 0, 0).found).toBe(false);
        expect(RootsBrentsMethod.find(F, 2, 1, 100, 0, 0, 0, 0).found).toBe(false);
        // maxIterations == 0.
        expect(RootsBrentsMethod.find(F, -1, 1, 0, 0, 0, 0, 0).found).toBe(false);
        // negFTolerance > 0.
        expect(RootsBrentsMethod.find(F, -1, 1, 100, 1e-8, 0, 0, 0).found).toBe(false);
        // posFTolerance < 0.
        expect(RootsBrentsMethod.find(F, -1, 1, 100, 0, -1e-8, 0, 0).found).toBe(false);
        // stepTTolerance < 0.
        expect(RootsBrentsMethod.find(F, -1, 1, 100, 0, 0, -1e-8, 0).found).toBe(false);
        // convTTolerance < 0.
        expect(RootsBrentsMethod.find(F, -1, 1, 100, 0, 0, 0, -1e-8).found).toBe(false);
    });

    it('terminates via the interval-size tolerance', () => {
        const F = (t: number) => t - 1 / 3;
        const r = RootsBrentsMethod.find(F, 0, 1, 200, 0, 0, 0, 1e-6);
        expect(r.found).toBe(true);
        expect(Math.abs(r.root - 1 / 3)).toBeLessThanOrEqual(1e-6);
    });

    it('handles a pathological flat region (ninth-order contact)', () => {
        const F = (t: number) => Math.pow(t - 1, 9);
        const { found, root } = findTight(F, 0, 1.75);
        expect(found).toBe(true);
        expect(Math.abs(root - 1)).toBeLessThanOrEqual(1e-12);
    });

    it('handles a discontinuous sign function', () => {
        const x0 = 0.7234;
        const F = (t: number) => Math.sign(t - x0);
        const { found, root } = findTight(F, 0, 1);
        expect(found).toBe(true);
        expect(Math.abs(root - x0)).toBeLessThanOrEqual(1e-15);
    });

    it('converges much faster than bisection on a smooth function', () => {
        // Count function evaluations; Brent should need far fewer than the
        // ~50+ bisections required for machine-precision convergence.
        let count = 0;
        const F = (t: number) => { ++count; return Math.exp(t) - 2; };
        const { found, root } = findTight(F, 0, 1);
        expect(found).toBe(true);
        expect(Math.abs(root - Math.LN2)).toBeLessThanOrEqual(Number.EPSILON * 2);
        expect(count).toBeLessThan(30);
    });

    it('fails to converge when the iteration budget is too small', () => {
        // One iteration is not enough for exact-zero tolerances on a
        // transcendental function; the algorithm reports failure rather
        // than an insufficiently accurate root.
        const r = RootsBrentsMethod.find(Math.cos, 0, 3, 1, 0, 0, 0, 0);
        expect(r.found).toBe(false);
    });

    it('solves a linear function in a single secant step', () => {
        // The secant iterate is exact for linear F, so even a budget of one
        // iteration finds the root with exact-zero tolerances.
        const F = (t: number) => t - 1 / 3;
        const r = RootsBrentsMethod.find(F, 0, 1, 1, 0, 0, 0, 0);
        expect(r.found).toBe(true);
        expect(r.root).toBe(1 / 3);
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V44): property-based comparison against upstream
// RootsBrentsMethod.h.
// ---------------------------------------------------------------------------

// A strictly increasing function whose only real root is r.
function monotoneF(r: number, c: number): (t: number) => number {
    return (t: number) => {
        const d = t - r;
        return d * d * d + c * d;
    };
}

const brentBracket = fc.tuple(wellScaled(-5, 5), positive(4, 0.5),
    positive(3, 0.1), positive(3, 0.1));

describe('RootsBrentsMethod verification', () => {
    it('rejects every invalid parameter combination', () => {
        check(fc.tuple(wellScaled(-5, 5), positive(3, 0.1), positive(1, 1e-3)),
            ([r, width, eps]) => {
                const F = monotoneF(r, 1);
                const t0 = r - width, t1 = r + width;
                const ok = (a: number, b: number, mi: number, nf: number,
                    pf: number, st: number, ct: number) =>
                    RootsBrentsMethod.find(F, a, b, mi, nf, pf, st, ct).found;
                // t1 <= t0
                expect(ok(t1, t0, 64, 0, 0, 0, 0)).toBe(false);
                expect(ok(t0, t0, 64, 0, 0, 0, 0)).toBe(false);
                // maxIterations == 0
                expect(ok(t0, t1, 0, 0, 0, 0, 0)).toBe(false);
                // negFTolerance > 0, posFTolerance < 0
                expect(ok(t0, t1, 64, eps, 0, 0, 0)).toBe(false);
                expect(ok(t0, t1, 64, 0, -eps, 0, 0)).toBe(false);
                // stepTTolerance < 0, convTTolerance < 0
                expect(ok(t0, t1, 64, 0, 0, -eps, 0)).toBe(false);
                expect(ok(t0, t1, 64, 0, 0, 0, -eps)).toBe(false);
            });
    });

    it('reports failure when the endpoint values have the same sign', () => {
        check(fc.tuple(wellScaled(-5, 5), positive(4, 0.5), positive(3, 0.1),
            positive(3, 0.1)),
            ([r, c, gap, width]) => {
                const F = monotoneF(r, c);
                const out = RootsBrentsMethod.find(F, r + gap, r + gap + width,
                    64, 0, 0, 0, 0);
                expect(out.found).toBe(false);
            });
    });

    it('accepts an endpoint that already satisfies the function tolerance', () => {
        check(fc.tuple(wellScaled(-5, 5), positive(4, 0.5), positive(3, 0.1)),
            ([r, c, width]) => {
                const F = monotoneF(r, c);
                const low = RootsBrentsMethod.find(F, r, r + width, 64, 0, 0, 0, 0);
                expect(low.found).toBe(true);
                expect(low.root === r).toBe(true);
                const high = RootsBrentsMethod.find(F, r - width, r, 64, 0, 0, 0, 0);
                expect(high.found).toBe(true);
                // The low endpoint is tested first, so a root at the high
                // endpoint is only reported after F(t0) fails the tolerance.
                expect(high.root === r).toBe(true);
            });
    });

    it('converges to the root with zero tolerances', () => {
        check(brentBracket, ([r, c, below, above]) => {
            const F = monotoneF(r, c);
            const out = RootsBrentsMethod.find(F, r - below, r + above, 500,
                0, 0, 0, 0);
            expect(out.found).toBe(true);
            expectClose(out.root, r, 1e-13, 1e-13);
        });
    });

    it('agrees with RootsBisection on the same bracket', () => {
        check(brentBracket, ([r, c, below, above]) => {
            const F = monotoneF(r, c);
            const t0 = r - below, t1 = r + above;
            const brent = RootsBrentsMethod.find(F, t0, t1, 500, 0, 0, 0, 0);
            const bisect = RootsBisection.find(F, t0, t1, 4096);
            expect(brent.found).toBe(true);
            expect(bisect.iterations).toBeGreaterThanOrEqual(1);
            // Two independent algorithms on the same bracket must agree to
            // the resolution of the doubles near the root.
            expectClose(brent.root, bisect.root, 1e-13, 1e-13);
        });
    });

    it('honors the function tolerance band', () => {
        check(fc.tuple(brentBracket, positive(1e-2, 1e-6)),
            ([[r, c, below, above], tol]) => {
                const F = monotoneF(r, c);
                const out = RootsBrentsMethod.find(F, r - below, r + above, 500,
                    -tol, tol, 0, 0);
                expect(out.found).toBe(true);
                const f = F(out.root);
                // Either the accepted estimate is inside the band, or the
                // bracket collapsed to consecutive doubles / the convergence
                // tolerance, both of which give a root accurate to an ulp.
                const inBand = -tol <= f && f <= tol;
                expect(inBand || Math.abs(out.root - r)
                    <= 1e-13 * (1 + Math.abs(r))).toBe(true);
            });
    });

    it('honors the subinterval convergence tolerance', () => {
        check(fc.tuple(brentBracket, positive(1e-2, 1e-6)),
            ([[r, c, below, above], convT]) => {
                const F = monotoneF(r, c);
                const out = RootsBrentsMethod.find(F, r - below, r + above, 500,
                    0, 0, 0, convT);
                expect(out.found).toBe(true);
                // Terminating on |t1 - t0| <= convT leaves the estimate within
                // convT of the true root; the other exits are far tighter.
                expect(Math.abs(out.root - r))
                    .toBeLessThanOrEqual(convT + 1e-13 * (1 + Math.abs(r)));
            });
    });

    it('never returns a root outside the input bracket', () => {
        check(fc.tuple(wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
            wellScaled(-3, 3)),
            ([a, b, c, d]) => {
                const t0 = Math.min(a, b), t1 = Math.max(a, b);
                const F = (t: number) => (t - c) * (t - d) * (t * t + 1);
                const out = RootsBrentsMethod.find(F, t0, t1, 500, 0, 0, 0, 0);
                if (!out.found) { return; }
                expect(out.root).toBeGreaterThanOrEqual(t0);
                expect(out.root).toBeLessThanOrEqual(t1);
                const err = Math.min(Math.abs(out.root - c), Math.abs(out.root - d));
                expect(err).toBeLessThanOrEqual(1e-11 * (1 + Math.abs(out.root)));
            });
    });
});
