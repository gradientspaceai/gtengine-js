import { describe, it, expect } from 'vitest';
import { HermiteQuintic, HermiteQuinticSample } from '../src/HermiteQuintic.js';
import { check, expectClose, fc, finite, scaled } from './helpers/arbitraries.js';

describe('HermiteQuintic', () => {
    it('default constructor creates the identically zero polynomial', () => {
        const h = new HermiteQuintic();
        for (let order = 0; order <= 5; ++order) {
            for (const x of [0, 0.25, 0.5, 0.75, 1]) {
                expect(h.evaluate(order, x)).toBe(0);
            }
        }
    });

    it('sample default constructor zero-fills', () => {
        const s = new HermiteQuinticSample();
        expect(s.f).toBe(0);
        expect(s.fx).toBe(0);
        expect(s.fxx).toBe(0);
    });

    it('reproduces function, first and second derivatives at the endpoints', () => {
        const s0 = new HermiteQuinticSample(2, -1, 4);
        const s1 = new HermiteQuinticSample(1, 3, -2);
        const h = new HermiteQuintic([s0, s1]);

        expect(h.evaluate(0, 0)).toBeCloseTo(s0.f, 14);
        expect(h.evaluate(1, 0)).toBeCloseTo(s0.fx, 14);
        expect(h.evaluate(2, 0)).toBeCloseTo(s0.fxx, 13);
        expect(h.evaluate(0, 1)).toBeCloseTo(s1.f, 14);
        expect(h.evaluate(1, 1)).toBeCloseTo(s1.fx, 14);
        expect(h.evaluate(2, 1)).toBeCloseTo(s1.fxx, 13);
    });

    it('matches a hand-computed value at x = 1/2', () => {
        // Coefficients generated from the samples are
        //   c0 = F0 = 2
        //   c1 = 5*F0 + Fx0 = 9
        //   c2 = 10*F0 + 4*Fx0 + Fxx0/2 = 18
        //   c3 = 10*F1 - 4*Fx1 + Fxx1/2 = -3
        //   c4 = 5*F1 - Fx1 = 2
        //   c5 = F1 = 1
        // At x = 1/2 every basis value P(i,0,1/2) = (1/2)^5 = 1/32, so
        //   H(1/2) = (2 + 9 + 18 - 3 + 2 + 1) / 32 = 29/32 = 0.90625.
        const h = new HermiteQuintic([
            new HermiteQuinticSample(2, -1, 4),
            new HermiteQuinticSample(1, 3, -2)
        ]);
        expect(h.evaluate(0, 0.5)).toBe(0.90625);
    });

    it('reproduces an arbitrary quintic polynomial on [0,1]', () => {
        // g(x) = x^5 - x^4 + 2 x^3 - x + 1; the interpolator matches g, g'
        // and g'' at both endpoints, and the quintic space has dimension 6,
        // so H = g.
        const g = (x: number): number =>
            (((((x - 1) * x + 2) * x) * x - 1) * x) + 1;
        const dg = (x: number): number =>
            (((5 * x - 4) * x + 6) * x) * x - 1;
        const ddg = (x: number): number =>
            ((20 * x - 12) * x + 12) * x;
        const h = new HermiteQuintic([
            new HermiteQuinticSample(g(0), dg(0), ddg(0)),
            new HermiteQuinticSample(g(1), dg(1), ddg(1))
        ]);
        for (let i = 0; i <= 10; ++i) {
            const x = i / 10;
            expect(h.evaluate(0, x)).toBeCloseTo(g(x), 13);
            expect(h.evaluate(1, x)).toBeCloseTo(dg(x), 12);
            expect(h.evaluate(2, x)).toBeCloseTo(ddg(x), 12);
        }
    });

    it('basis functions match P(i,t) = (1-t)^{5-i} t^i', () => {
        for (let i = 0; i <= 5; ++i) {
            for (let k = 0; k <= 8; ++k) {
                const t = k / 8;
                const expected = Math.pow(1 - t, 5 - i) * Math.pow(t, i);
                expect(HermiteQuintic.p(i, 0, t)).toBeCloseTo(expected, 14);
            }
        }
    });

    it('basis first derivatives match the analytic differentiation', () => {
        // d/dt [(1-t)^{5-i} t^i]
        //   = -(5-i)(1-t)^{4-i} t^i + i (1-t)^{5-i} t^{i-1}, 0 < t < 1.
        for (let i = 0; i <= 5; ++i) {
            for (let k = 1; k <= 7; ++k) {
                const t = k / 8;
                const expected =
                    -(5 - i) * Math.pow(1 - t, 4 - i) * Math.pow(t, i) +
                    i * Math.pow(1 - t, 5 - i) * Math.pow(t, i - 1);
                expect(HermiteQuintic.p(i, 1, t)).toBeCloseTo(expected, 12);
            }
        }
    });

    it('derivative orders agree with central finite differences', () => {
        const h = new HermiteQuintic([
            new HermiteQuinticSample(2, -1, 4),
            new HermiteQuinticSample(1, 3, -2)
        ]);
        const step = 1e-5;
        for (let order = 0; order <= 4; ++order) {
            for (const x of [0.1, 0.35, 0.6, 0.85]) {
                const fd =
                    (h.evaluate(order, x + step) - h.evaluate(order, x - step)) /
                    (2 * step);
                expect(h.evaluate(order + 1, x)).toBeCloseTo(fd, 4);
            }
        }
    });

    it('returns zero for derivative orders beyond the degree', () => {
        const h = new HermiteQuintic([
            new HermiteQuinticSample(2, -1, 4),
            new HermiteQuinticSample(1, 3, -2)
        ]);
        expect(h.evaluate(6, 0.5)).toBe(0);
        expect(h.evaluate(9, 0.25)).toBe(0);
        expect(HermiteQuintic.p(3, 6, 0.5)).toBe(0);
    });

    it('generate replaces the coefficients of an existing polynomial', () => {
        const h = new HermiteQuintic([
            new HermiteQuinticSample(2, -1, 4),
            new HermiteQuinticSample(1, 3, -2)
        ]);
        h.generate([
            new HermiteQuinticSample(7, 0, 0),
            new HermiteQuinticSample(7, 0, 0)
        ]);
        for (const x of [0, 0.3, 0.5, 1]) {
            expect(h.evaluate(0, x)).toBeCloseTo(7, 13);
        }
        expect(h.evaluate(1, 0.5)).toBeCloseTo(0, 13);
        expect(h.evaluate(2, 0.5)).toBeCloseTo(0, 12);
    });
});

// ---------------------------------------------------------------------------
// Verification (V28): property-based cross-checks against HermiteQuintic.h.
// ---------------------------------------------------------------------------
describe('HermiteQuintic verification', () => {
    const basis = (i: number, t: number): number =>
        Math.pow(1 - t, 5 - i) * Math.pow(t, i);

    it('p(select, 0, t) equals (1-t)^{5-i} t^i', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 5 }), scaled(-1, 2)), ([i, t]) => {
            expectClose(HermiteQuintic.p(i, 0, t), basis(i, t), 1e-14, 1e-13);
        });
    });

    it('p(select, order, t) matches a central difference of p(select, order-1)', () => {
        const h = 1e-4;
        check(fc.tuple(fc.integer({ min: 0, max: 5 }),
            fc.integer({ min: 1, max: 5 }), scaled(0.1, 0.9)), ([i, order, t]) => {
            const fd = (HermiteQuintic.p(i, order - 1, t + h)
                - HermiteQuintic.p(i, order - 1, t - h)) / (2 * h);
            // Truncation is (h^2/6) * P^{(order+2)}, at most O(120*1e-8/6);
            // roundoff is O(eps * 120 / h) = O(1e-11).
            expectClose(HermiteQuintic.p(i, order, t), fd, 1e-5, 1e-5);
        });
    });

    it('interpolates value, slope and curvature at the cell corners', () => {
        check(fc.array(finite(-5, 5), { minLength: 6, maxLength: 6 }), s => {
            const h = new HermiteQuintic([new HermiteQuinticSample(s[0], s[1], s[2]),
                new HermiteQuinticSample(s[3], s[4], s[5])]);
            expectClose(h.evaluate(0, 0), s[0], 1e-13, 1e-13);
            expectClose(h.evaluate(1, 0), s[1], 1e-11, 1e-11);
            expectClose(h.evaluate(2, 0), s[2], 1e-10, 1e-10);
            expectClose(h.evaluate(0, 1), s[3], 1e-13, 1e-13);
            expectClose(h.evaluate(1, 1), s[4], 1e-11, 1e-11);
            expectClose(h.evaluate(2, 1), s[5], 1e-10, 1e-10);
        });
    });

    it('reproduces every quintic from exact samples of f, f\' and f\'\'', () => {
        check(fc.tuple(fc.array(finite(-3, 3), { minLength: 6, maxLength: 6 }),
            scaled(0, 1)), ([a, x]) => {
            const p = (t: number) => a.reduceRight((acc, ai) => acc * t + ai, 0);
            const d1 = (t: number) => a.slice(1).map((ai, i) => (i + 1) * ai)
                .reduceRight((acc, ai) => acc * t + ai, 0);
            const d2 = (t: number) => a.slice(2).map((ai, i) => (i + 1) * (i + 2) * ai)
                .reduceRight((acc, ai) => acc * t + ai, 0);
            const h = new HermiteQuintic([
                new HermiteQuinticSample(p(0), d1(0), d2(0)),
                new HermiteQuinticSample(p(1), d1(1), d2(1))]);
            // The largest coefficient multiplier in generateSingle is 10, and
            // the basis functions are O(1) on [0,1], so the reconstruction
            // keeps about 13 digits.
            expectClose(h.evaluate(0, x), p(x), 1e-11, 1e-11);
            expectClose(h.evaluate(1, x), d1(x), 1e-10, 1e-10);
            expectClose(h.evaluate(2, x), d2(x), 1e-9, 1e-9);
        });
    });

    it('is linear in the sample data', () => {
        check(fc.tuple(fc.array(finite(-3, 3), { minLength: 6, maxLength: 6 }),
            fc.array(finite(-3, 3), { minLength: 6, maxLength: 6 }),
            scaled(0, 1)), ([u, v, x]) => {
            const mk = (s: number[]) => new HermiteQuintic([
                new HermiteQuinticSample(s[0], s[1], s[2]),
                new HermiteQuinticSample(s[3], s[4], s[5])]);
            const sum = u.map((ui, i) => ui + v[i]);
            for (let order = 0; order <= 5; ++order) {
                expectClose(mk(sum).evaluate(order, x),
                    mk(u).evaluate(order, x) + mk(v).evaluate(order, x), 1e-10, 1e-10);
            }
        });
    });

    it('returns zero for derivative orders above the degree', () => {
        check(fc.tuple(fc.integer({ min: 6, max: 40 }), scaled(0, 1)), ([order, x]) => {
            const h = new HermiteQuintic([new HermiteQuinticSample(1, 2, 3),
                new HermiteQuinticSample(-3, 4, -5)]);
            expect(h.evaluate(order, x)).toBe(0);
            expect(HermiteQuintic.p(3, order, x)).toBe(0);
        });
    });
});
