import { describe, it, expect } from 'vitest';
import { PolynomialRoot, polynomialRootBisect } from '../src/PolynomialRoot';

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
