import { describe, it, expect } from 'vitest';
import { HermiteCubic, HermiteCubicSample } from '../src/HermiteCubic';

describe('HermiteCubic', () => {
    it('default constructor creates the identically zero polynomial', () => {
        const h = new HermiteCubic();
        for (let order = 0; order <= 3; ++order) {
            for (const x of [0, 0.25, 0.5, 0.75, 1]) {
                expect(h.evaluate(order, x)).toBe(0);
            }
        }
    });

    it('sample default constructor zero-fills', () => {
        const s = new HermiteCubicSample();
        expect(s.f).toBe(0);
        expect(s.fx).toBe(0);
    });

    it('reproduces function and derivative values at the endpoints', () => {
        const s0 = new HermiteCubicSample(1, 2);
        const s1 = new HermiteCubicSample(3, -1);
        const h = new HermiteCubic([s0, s1]);

        expect(h.evaluate(0, 0)).toBeCloseTo(s0.f, 14);
        expect(h.evaluate(1, 0)).toBeCloseTo(s0.fx, 14);
        expect(h.evaluate(0, 1)).toBeCloseTo(s1.f, 14);
        expect(h.evaluate(1, 1)).toBeCloseTo(s1.fx, 14);
    });

    it('matches hand-computed values at x = 1/2', () => {
        // Coefficients generated from the samples are
        //   c0 = F0 = 1, c1 = 3*F0 + Fx0 = 5,
        //   c2 = 3*F1 - Fx1 = 10, c3 = F1 = 3.
        // At x = 1/2 every basis value P(i,0,1/2) = (1/2)^3 = 1/8, so
        //   H(1/2) = (1 + 5 + 10 + 3) / 8 = 19/8 = 2.375.
        // The basis derivatives at 1/2 are -3/4, -1/4, 1/4, 3/4, so
        //   H'(1/2) = -3/4 - 5/4 + 10/4 + 9/4 = 11/4 = 2.75.
        const h = new HermiteCubic([
            new HermiteCubicSample(1, 2),
            new HermiteCubicSample(3, -1)
        ]);
        expect(h.evaluate(0, 0.5)).toBe(2.375);
        expect(h.evaluate(1, 0.5)).toBe(2.75);
    });

    it('reproduces an arbitrary cubic polynomial on [0,1]', () => {
        // g(x) = x^3 - 2 x^2 + 3 x - 1; the interpolator matches g and g'
        // at both endpoints, and the cubic space has dimension 4, so H = g.
        const g = (x: number): number => ((x - 2) * x + 3) * x - 1;
        const dg = (x: number): number => (3 * x - 4) * x + 3;
        const h = new HermiteCubic([
            new HermiteCubicSample(g(0), dg(0)),
            new HermiteCubicSample(g(1), dg(1))
        ]);
        for (let i = 0; i <= 10; ++i) {
            const x = i / 10;
            expect(h.evaluate(0, x)).toBeCloseTo(g(x), 13);
            expect(h.evaluate(1, x)).toBeCloseTo(dg(x), 13);
        }
    });

    it('basis functions match P(i,t) = (1-t)^{3-i} t^i', () => {
        for (let i = 0; i <= 3; ++i) {
            for (let k = 0; k <= 8; ++k) {
                const t = k / 8;
                const expected = Math.pow(1 - t, 3 - i) * Math.pow(t, i);
                expect(HermiteCubic.p(i, 0, t)).toBeCloseTo(expected, 14);
            }
        }
    });

    it('basis first derivatives match the analytic differentiation', () => {
        // d/dt [(1-t)^{3-i} t^i]
        //   = -(3-i)(1-t)^{2-i} t^i + i (1-t)^{3-i} t^{i-1}, 0 < t < 1.
        for (let i = 0; i <= 3; ++i) {
            for (let k = 1; k <= 7; ++k) {
                const t = k / 8;
                const expected =
                    -(3 - i) * Math.pow(1 - t, 2 - i) * Math.pow(t, i) +
                    i * Math.pow(1 - t, 3 - i) * Math.pow(t, i - 1);
                expect(HermiteCubic.p(i, 1, t)).toBeCloseTo(expected, 12);
            }
        }
    });

    it('derivative orders agree with central finite differences', () => {
        const h = new HermiteCubic([
            new HermiteCubicSample(1, 2),
            new HermiteCubicSample(3, -1)
        ]);
        const step = 1e-5;
        for (let order = 0; order <= 2; ++order) {
            for (const x of [0.1, 0.35, 0.6, 0.85]) {
                const fd =
                    (h.evaluate(order, x + step) - h.evaluate(order, x - step)) /
                    (2 * step);
                expect(h.evaluate(order + 1, x)).toBeCloseTo(fd, 5);
            }
        }
    });

    it('returns zero for derivative orders beyond the degree', () => {
        const h = new HermiteCubic([
            new HermiteCubicSample(1, 2),
            new HermiteCubicSample(3, -1)
        ]);
        expect(h.evaluate(4, 0.5)).toBe(0);
        expect(h.evaluate(7, 0.25)).toBe(0);
        expect(HermiteCubic.p(2, 4, 0.5)).toBe(0);
    });

    it('generate replaces the coefficients of an existing polynomial', () => {
        const h = new HermiteCubic([
            new HermiteCubicSample(1, 2),
            new HermiteCubicSample(3, -1)
        ]);
        h.generate([new HermiteCubicSample(5, 0), new HermiteCubicSample(5, 0)]);
        for (const x of [0, 0.3, 0.5, 1]) {
            expect(h.evaluate(0, x)).toBeCloseTo(5, 13);
        }
        expect(h.evaluate(1, 0)).toBeCloseTo(0, 13);
        expect(h.evaluate(1, 1)).toBeCloseTo(0, 13);
    });
});
