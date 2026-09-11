import { describe, it, expect } from 'vitest';
import { PolynomialRoot, polynomialRootBisect } from '../src/PolynomialRoot.js';
import {
    check, expectClose, expectStrictWeakOrder, fc, positive, wellScaled
} from './helpers/arbitraries.js';

describe('PolynomialRoot', () => {
    it('default-constructs an invalid root (x = 0, m = 0)', () => {
        const root = new PolynomialRoot();
        expect(root.x).toBe(0);
        expect(root.m).toBe(0);
    });

    it('stores the root estimate and its multiplicity', () => {
        const root = new PolynomialRoot(-2.5, 3);
        expect(root.x).toBe(-2.5);
        expect(root.m).toBe(3);
    });

    it('compares only the root estimate, ignoring the multiplicity', () => {
        const a = new PolynomialRoot(1.5, 1);
        const b = new PolynomialRoot(1.5, 4);
        const c = new PolynomialRoot(2.5, 1);
        expect(a.equals(b)).toBe(true);
        expect(b.equals(a)).toBe(true);
        expect(a.equals(c)).toBe(false);
        expect(a.lessThan(c)).toBe(true);
        expect(c.lessThan(a)).toBe(false);
        expect(a.lessThan(b)).toBe(false);
        expect(b.lessThan(a)).toBe(false);
    });

    it('sorts a list of roots by the estimate', () => {
        const roots = [new PolynomialRoot(3, 1), new PolynomialRoot(-1, 2),
            new PolynomialRoot(0.5, 1)];
        roots.sort((r0, r1) => (r0.lessThan(r1) ? -1 : (r1.lessThan(r0) ? 1 : 0)));
        expect(roots.map(r => r.x)).toEqual([-1, 0.5, 3]);
    });
});

describe('polynomialRootBisect', () => {
    it('brackets sqrt(2) by consecutive floating-point numbers', () => {
        const F = (x: number): number => x * x - 2;
        const result = polynomialRootBisect(F, -1, +1, 1, 2);
        expect(result.xMin).toBeLessThanOrEqual(Math.SQRT2);
        expect(result.xMax).toBeGreaterThanOrEqual(Math.SQRT2);
        // Either the function is exactly zero at a common endpoint, or the
        // endpoints are consecutive floating-point numbers.
        if (result.xMin !== result.xMax) {
            expect(0.5 * (result.xMin + result.xMax) === result.xMin ||
                0.5 * (result.xMin + result.xMax) === result.xMax).toBe(true);
            expect(F(result.xMin)).toBeLessThan(0);
            expect(F(result.xMax)).toBeGreaterThan(0);
        }
        expect(Math.abs(result.xMin - Math.SQRT2)).toBeLessThan(1e-15);
    });

    it('collapses the interval to the exact root when F is zero there', () => {
        // The midpoint of [0,1] is 0.5, where F is exactly zero.
        const F = (x: number): number => x - 0.5;
        const result = polynomialRootBisect(F, -1, +1, 0, 1);
        expect(result.xMin).toBe(0.5);
        expect(result.xMax).toBe(0.5);
    });

    it('handles a decreasing function (signFMin = +1, signFMax = -1)', () => {
        const F = (x: number): number => 2 - x * x;
        const result = polynomialRootBisect(F, +1, -1, 1, 2);
        expect(result.xMin).toBeCloseTo(Math.SQRT2, 15);
        expect(result.xMax).toBeCloseTo(Math.SQRT2, 15);
        expect(F(result.xMin)).toBeGreaterThanOrEqual(0);
        expect(F(result.xMax)).toBeLessThanOrEqual(0);
    });

    it('collapses to xMin when the sign at xMin is not the claimed sign', () => {
        // Upstream: rounding errors prevent the correct classification of the
        // multiplicity of roots, so the interval degenerates.
        const F = (x: number): number => x * x + 1;
        const result = polynomialRootBisect(F, -1, +1, 1, 2);
        expect(result.xMin).toBe(1);
        expect(result.xMax).toBe(1);
    });

    it('collapses to xMax when the sign at xMax is not the claimed sign', () => {
        // F(1) = -1 matches signFMin, but F(2) = -1 does not match signFMax.
        const F = (x: number): number => -1;
        const result = polynomialRootBisect(F, -1, +1, 1, 2);
        expect(result.xMin).toBe(2);
        expect(result.xMax).toBe(2);
    });

    it('returns immediately when F is zero at xMin (sign mismatch)', () => {
        // A zero at an endpoint has true sign 0, which never equals the
        // required +1 or -1, so the interval collapses.
        const F = (x: number): number => x;
        const result = polynomialRootBisect(F, -1, +1, 0, 1);
        expect(result.xMin).toBe(0);
        expect(result.xMax).toBe(0);
    });

    it('bisects transcendental functions to machine precision', () => {
        const cases: { F: (x: number) => number, xMin: number, xMax: number,
            signMin: number, signMax: number, expected: number }[] = [
            { F: (x) => Math.cos(x), xMin: 1, xMax: 2, signMin: +1, signMax: -1,
                expected: Math.PI / 2 },
            { F: (x) => Math.exp(x) - 3, xMin: 0, xMax: 2, signMin: -1,
                signMax: +1, expected: Math.log(3) },
            { F: (x) => x * x * x - 7, xMin: 1, xMax: 3, signMin: -1,
                signMax: +1, expected: Math.cbrt(7) }
        ];
        for (const c of cases) {
            const result = polynomialRootBisect(c.F, c.signMin, c.signMax,
                c.xMin, c.xMax);
            expect(result.xMin).toBeLessThanOrEqual(result.xMax);
            expect(Math.abs(result.xMin - c.expected))
                .toBeLessThan(4 * Number.EPSILON * Math.abs(c.expected));
            expect(Math.abs(result.xMax - c.expected))
                .toBeLessThan(4 * Number.EPSILON * Math.abs(c.expected));
        }
    });

    it('keeps the root inside the returned interval for random line functions', () => {
        // F(x) = a * (x - r) has the exact root r; bisection must bracket it.
        let s = 4242 >>> 0;
        const rand = (): number => {
            s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
            return s / 4294967296;
        };
        for (let trial = 0; trial < 200; ++trial) {
            const r = 10 * (2 * rand() - 1);
            const a = 1 + 9 * rand();
            const F = (x: number): number => a * (x - r);
            const xMin = r - (0.5 + rand());
            const xMax = r + (0.5 + rand());
            const result = polynomialRootBisect(F, -1, +1, xMin, xMax);
            expect(result.xMin).toBeLessThanOrEqual(result.xMax);
            expect(Math.abs(result.xMin - r))
                .toBeLessThan(1e-14 * (1 + Math.abs(r)));
            expect(Math.abs(result.xMax - r))
                .toBeLessThan(1e-14 * (1 + Math.abs(r)));
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V44): property-based comparison against upstream
// PolynomialRoot.h.
// ---------------------------------------------------------------------------

// A strictly increasing function with the single real root r:
// F(x) = (x - r)^3 + c * (x - r), c >= 0.
function monotoneF(r: number, c: number): (x: number) => number {
    return (x: number) => {
        const d = x - r;
        return d * d * d + c * d;
    };
}

const sign = (v: number): number => (v > 0 ? +1 : (v < 0 ? -1 : 0));

describe('PolynomialRoot verification', () => {
    it('lessThan is a strict weak ordering that ignores the multiplicity', () => {
        check(fc.array(fc.tuple(fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: 1, max: 3 })), { minLength: 2, maxLength: 6 }),
            pairs => {
                const roots = pairs.map(([x, m]) => new PolynomialRoot(x, m));
                expectStrictWeakOrder(roots, (a, b) => a.lessThan(b));
                for (const a of roots) {
                    for (const b of roots) {
                        // operator== and operator< compare only x upstream.
                        expect(a.equals(b)).toBe(a.x === b.x);
                    }
                }
            }, 50);
    });

    it('brackets the true root and stays inside the initial interval', () => {
        // c and the half-widths are bounded away from zero: for a tiny c and
        // a tiny interval the products c*(x-r) underflow to exactly zero near
        // the root and bisection stops at a subnormal, which says nothing
        // about the bracketing contract.
        check(fc.tuple(wellScaled(-5, 5), positive(4, 0.5), positive(3, 0.1),
            positive(3, 0.1)),
            ([r, c, below, above]) => {
                const F = monotoneF(r, c);
                const xMin = r - below;
                const xMax = r + above;
                // F is increasing, so the theoretical signs are -1 and +1.
                const out = polynomialRootBisect(F, -1, +1, xMin, xMax);
                expect(out.xMin).toBeGreaterThanOrEqual(xMin);
                expect(out.xMax).toBeLessThanOrEqual(xMax);
                expect(out.xMin).toBeLessThanOrEqual(out.xMax);
                if (out.xMin === out.xMax) {
                    // Upstream terminates here only with F exactly zero.
                    // Compare with === so that a -0 result does not fail the
                    // Object.is comparison of toBe.
                    expect(F(out.xMin) === 0).toBe(true);
                } else {
                    // Otherwise the endpoints keep the theoretical signs and
                    // are consecutive doubles (subdivision cannot separate
                    // them any further).
                    expect(F(out.xMin)).toBeLessThan(0);
                    expect(F(out.xMax)).toBeGreaterThan(0);
                    const mid = 0.5 * (out.xMin + out.xMax);
                    expect(mid === out.xMin || mid === out.xMax).toBe(true);
                }
                // Either way the bracket contains the true root r, up to the
                // rounding of the last F evaluation.
                if (out.xMin !== out.xMax) {
                    expect(out.xMin).toBeLessThanOrEqual(r);
                    expect(out.xMax).toBeGreaterThanOrEqual(r);
                } else {
                    expectClose(out.xMin, r, 1e-14, 1e-14);
                }
            });
    });

    it('handles the decreasing orientation symmetrically', () => {
        check(fc.tuple(wellScaled(-5, 5), positive(4, 0.5), positive(3, 0.1),
            positive(3, 0.1)),
            ([r, c, below, above]) => {
                const G = monotoneF(r, c);
                const F = (x: number) => -G(x);
                const out = polynomialRootBisect(F, +1, -1, r - below, r + above);
                expect(out.xMin).toBeLessThanOrEqual(r);
                expect(out.xMax).toBeGreaterThanOrEqual(r);
                expect(sign(F(out.xMin))).toBeGreaterThanOrEqual(0);
                expect(sign(F(out.xMax))).toBeLessThanOrEqual(0);
            });
    });

    it('collapses the interval when the claimed endpoint signs are wrong', () => {
        // Upstream's documented failure mode: floating-point rounding errors
        // prevent the correct classification of the multiplicity of roots, so
        // the interval degenerates to one endpoint rather than reporting an
        // error. A root exactly at an endpoint has true sign 0, which never
        // matches the required +-1, so it collapses too (gtengine-js #319).
        check(fc.tuple(wellScaled(-5, 5), positive(4, 0.5), positive(3, 0.1)),
            ([r, c, width]) => {
                const F = monotoneF(r, c);
                // Both endpoints above the root: F > 0 at both.
                let out = polynomialRootBisect(F, -1, +1, r + width, r + 2 * width);
                expect(out.xMin).toBe(r + width);
                expect(out.xMax).toBe(r + width);
                // Correct sign at xMin, wrong sign at xMax.
                out = polynomialRootBisect(F, -1, -1, r - width, r + width);
                expect(out.xMin).toBe(r + width);
                expect(out.xMax).toBe(r + width);
                // The root sits exactly on the low endpoint.
                out = polynomialRootBisect(F, -1, +1, r, r + width);
                expect(out.xMin).toBe(r);
                expect(out.xMax).toBe(r);
            });
    });

    it('terminates and returns a finite bracket for every ordered interval', () => {
        check(fc.tuple(wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
            fc.integer({ min: -1, max: 1 }), fc.integer({ min: -1, max: 1 })),
            ([a, b, c, s0, s1]) => {
                const lo = Math.min(a, b);
                const hi = Math.max(a, b);
                const F = (x: number) => (x - c) * (x * x + 1);
                const out = polynomialRootBisect(F, s0, s1, lo, hi);
                expect(Number.isFinite(out.xMin)).toBe(true);
                expect(Number.isFinite(out.xMax)).toBe(true);
                expect(out.xMin).toBeGreaterThanOrEqual(lo);
                expect(out.xMax).toBeLessThanOrEqual(hi);
                expect(out.xMin).toBeLessThanOrEqual(out.xMax);
            });
    });
});
