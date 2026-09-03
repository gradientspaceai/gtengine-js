import { describe, it, expect } from 'vitest';
import { IntpAkimaUniform1 } from '../src/IntpAkimaUniform1.js';

describe('IntpAkimaUniform1', () => {
    it('throws for invalid inputs', () => {
        expect(() => new IntpAkimaUniform1(2, 0, 1, [0, 1]))
            .toThrow('Invalid input to IntpAkima1 constructor.');
        expect(() => new IntpAkimaUniform1(3, 0, 0, [0, 1, 2]))
            .toThrow('Spacing must be positive.');
        expect(() => new IntpAkimaUniform1(3, 0, -1, [0, 1, 2]))
            .toThrow('Spacing must be positive.');
    });

    it('provides member access', () => {
        const F = [0, 1, 4, 9, 16];
        const interp = new IntpAkimaUniform1(5, -1, 0.5, F);
        expect(interp.getQuantity()).toBe(5);
        expect(interp.getF()).toBe(F);
        expect(interp.getXMin()).toBe(-1);
        expect(interp.getXSpacing()).toBe(0.5);
        expect(interp.getXMax()).toBeCloseTo(1, 15);
    });

    it('interpolates the samples exactly', () => {
        const xMin = -2;
        const dx = 0.75;
        const F = [1, -3, 2.5, 7, 0, -1.25, 4, 4];
        const interp = new IntpAkimaUniform1(F.length, xMin, dx, F);
        for (let i = 0; i < F.length; ++i) {
            expect(interp.evaluate(xMin + i * dx)).toBeCloseTo(F[i], 12);
            expect(interp.evaluate(0, xMin + i * dx)).toBeCloseTo(F[i], 12);
        }
    });

    it('reproduces linear data exactly', () => {
        // All slopes are equal, so the estimated derivatives are the common
        // slope and the quadratic and cubic coefficients vanish.
        const xMin = 1;
        const dx = 0.25;
        const a = -0.5, b = 3;
        const F: number[] = [];
        for (let i = 0; i < 6; ++i) {
            F.push(a + b * (xMin + i * dx));
        }
        const interp = new IntpAkimaUniform1(F.length, xMin, dx, F);
        for (let k = 0; k <= 40; ++k) {
            const x = xMin + (k / 40) * (F.length - 1) * dx;
            expect(interp.evaluate(x)).toBeCloseTo(a + b * x, 12);
            expect(interp.evaluate(1, x)).toBeCloseTo(b, 12);
            expect(interp.evaluate(2, x)).toBeCloseTo(0, 12);
            expect(interp.evaluate(3, x)).toBeCloseTo(0, 12);
        }
    });

    it('reproduces constant data exactly', () => {
        const F = [2, 2, 2, 2, 2];
        const interp = new IntpAkimaUniform1(5, 0, 1, F);
        for (let k = 0; k <= 20; ++k) {
            const x = (k / 20) * 4;
            expect(interp.evaluate(x)).toBeCloseTo(2, 14);
            expect(interp.evaluate(1, x)).toBeCloseTo(0, 14);
        }
    });

    it('averages the two slopes when both slope pairs are equal', () => {
        // The slopes are (1,1,1,2,2,2). At sample index 3 the four-slope
        // window is (1,1,2,2), so slope[0] == slope[1] and
        // slope[2] == slope[3] and the estimated derivative is
        // 0.5*(1 + 2) = 1.5.
        const F = [0, 1, 2, 3, 5, 7, 9];
        const interp = new IntpAkimaUniform1(F.length, 0, 1, F);
        expect(interp.evaluate(1, 3)).toBeCloseTo(1.5, 14);

        // Away from the corner the data is linear, so the derivative is the
        // local slope.
        expect(interp.evaluate(1, 1)).toBeCloseTo(1, 14);
        expect(interp.evaluate(1, 5)).toBeCloseTo(2, 14);
    });

    it('uses the weighted-average derivative formula', () => {
        // The slopes are (0,1,3,6). Prepending the extrapolated boundary
        // slopes, the window centered at sample index 2 is
        // (0, 1, 3, 6), so ad0 = |6 - 3| = 3, ad1 = |0 - 1| = 1 and the
        // derivative is (3*1 + 1*3)/(3 + 1) = 1.5.
        const F = [0, 0, 1, 4, 10];
        const interp = new IntpAkimaUniform1(F.length, 0, 1, F);
        expect(interp.evaluate(1, 2)).toBeCloseTo(1.5, 14);
    });

    it('is C1 continuous at the interior knots', () => {
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const interp = new IntpAkimaUniform1(F.length, 0, 1, F);
        const h = 1e-9;
        for (let i = 1; i < F.length - 1; ++i) {
            const left = interp.evaluate(1, i - h);
            const right = interp.evaluate(1, i + h);
            expect(left).toBeCloseTo(right, 6);
            expect(interp.evaluate(i - h)).toBeCloseTo(interp.evaluate(i + h), 6);
        }
    });

    it('matches finite differences of the function', () => {
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const interp = new IntpAkimaUniform1(F.length, 0, 1, F);
        const h = 1e-5;
        for (const x of [0.3, 1.7, 2.5, 3.1, 4.8, 5.4]) {
            const d1 = (interp.evaluate(x + h) - interp.evaluate(x - h)) / (2 * h);
            expect(interp.evaluate(1, x)).toBeCloseTo(d1, 6);
            const d2 = (interp.evaluate(1, x + h) - interp.evaluate(1, x - h)) / (2 * h);
            expect(interp.evaluate(2, x)).toBeCloseTo(d2, 5);
        }
    });

    it('clamps inputs outside the domain', () => {
        const F = [1, -3, 2.5, 7, 0];
        const xMin = -1, dx = 0.5;
        const interp = new IntpAkimaUniform1(F.length, xMin, dx, F);
        expect(interp.evaluate(-100)).toBeCloseTo(F[0], 14);
        expect(interp.evaluate(100)).toBeCloseTo(F[F.length - 1], 12);
        expect(interp.evaluate(1, -100)).toBe(interp.evaluate(1, xMin));
        expect(interp.evaluate(1, 100)).toBe(interp.evaluate(1, interp.getXMax()));
    });

    it('returns zero for derivative orders larger than three', () => {
        const F = [1, -3, 2.5, 7, 0];
        const interp = new IntpAkimaUniform1(5, 0, 1, F);
        expect(interp.evaluate(4, 1.5)).toBe(0);
        expect(interp.evaluate(7, 1.5)).toBe(0);
    });

    it('is cubic on each interval, matching the third derivative', () => {
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const interp = new IntpAkimaUniform1(F.length, 0, 1, F);
        // Within a single interval the polynomial has a constant third
        // derivative.
        const c3 = interp.evaluate(3, 2.25);
        expect(interp.evaluate(3, 2.75)).toBeCloseTo(c3, 14);
        // The second derivative varies linearly with x on the interval.
        const s0 = interp.evaluate(2, 2.25);
        const s1 = interp.evaluate(2, 2.75);
        expect(s1 - s0).toBeCloseTo(c3 * 0.5, 12);
    });

    it('handles the minimum sample count', () => {
        const F = [0, 2, 1];
        const interp = new IntpAkimaUniform1(3, 0, 1, F);
        expect(interp.evaluate(0)).toBeCloseTo(0, 14);
        expect(interp.evaluate(1)).toBeCloseTo(2, 14);
        expect(interp.evaluate(2)).toBeCloseTo(1, 14);
    });
});
