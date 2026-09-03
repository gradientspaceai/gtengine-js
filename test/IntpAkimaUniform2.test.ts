import { describe, it, expect } from 'vitest';
import { IntpAkimaUniform2 } from '../src/IntpAkimaUniform2.js';

// Build the row-major sample array F[c + xBound*r] = f(xMin + c*dx, yMin + r*dy).
function makeSamples(xBound: number, yBound: number, xMin: number, dx: number,
    yMin: number, dy: number, f: (x: number, y: number) => number): number[] {
    const F: number[] = [];
    for (let r = 0; r < yBound; ++r) {
        for (let c = 0; c < xBound; ++c) {
            F.push(f(xMin + c * dx, yMin + r * dy));
        }
    }
    return F;
}

describe('IntpAkimaUniform2', () => {
    it('throws for invalid inputs', () => {
        const F9 = new Array<number>(9).fill(0);
        expect(() => new IntpAkimaUniform2(2, 3, 0, 1, 0, 1, F9)).toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform2(3, 2, 0, 1, 0, 1, F9)).toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform2(3, 3, 0, 0, 0, 1, F9)).toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform2(3, 3, 0, 1, 0, -1, F9)).toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform2(4, 3, 0, 1, 0, 1, F9)).toThrow('Invalid input.');
    });

    it('provides member access', () => {
        const F = makeSamples(4, 5, -1, 0.5, 2, 0.25, (x, y) => x + y);
        const interp = new IntpAkimaUniform2(4, 5, -1, 0.5, 2, 0.25, F);
        expect(interp.getXBound()).toBe(4);
        expect(interp.getYBound()).toBe(5);
        expect(interp.getQuantity()).toBe(20);
        expect(interp.getF()).toBe(F);
        expect(interp.getXMin()).toBe(-1);
        expect(interp.getXMax()).toBeCloseTo(0.5, 14);
        expect(interp.getXSpacing()).toBe(0.5);
        expect(interp.getYMin()).toBe(2);
        expect(interp.getYMax()).toBeCloseTo(3, 14);
        expect(interp.getYSpacing()).toBe(0.25);
    });

    it('interpolates a hand-computed 3x3 grid of f(x,y) = x*y', () => {
        // For f = x*y all x-slopes in a row equal y (and y-slopes in a
        // column equal x), so the Akima estimates fx = y and fy = x are
        // exact everywhere, and the interior/min-boundary fxy estimates are
        // exactly 1. In the cell [0,1]^2 the polynomial is exactly x*y.
        const F = makeSamples(3, 3, 0, 1, 0, 1, (x, y) => x * y);
        const interp = new IntpAkimaUniform2(3, 3, 0, 1, 0, 1, F);

        for (const [x, y] of [[0.5, 0.5], [0.25, 0.75], [0, 0.3], [1, 1]]) {
            expect(interp.evaluate(x, y)).toBeCloseTo(x * y, 12);
            expect(interp.evaluate(0, 0, x, y)).toBeCloseTo(x * y, 12);
            expect(interp.evaluate(1, 0, x, y)).toBeCloseTo(y, 12);
            expect(interp.evaluate(0, 1, x, y)).toBeCloseTo(x, 12);
            expect(interp.evaluate(1, 1, x, y)).toBeCloseTo(1, 12);
            expect(interp.evaluate(2, 0, x, y)).toBeCloseTo(0, 12);
            expect(interp.evaluate(0, 2, x, y)).toBeCloseTo(0, 12);
        }

        // Upstream-behavior regression pin: GetFXY reuses the min-boundary
        // one-sided difference coefficients at the max boundaries without
        // negating them for the reversed direction, so the fxy estimates on
        // the top row and right column (except the top-right corner, where
        // the two sign flips cancel) come out as -1 instead of +1 for
        // f = x*y. Cells touching those points therefore deviate from x*y;
        // the value below is hand-computed from the upstream algorithm.
        expect(interp.evaluate(1.5, 0.75)).toBeCloseTo(1.1015625, 12);
        // The samples themselves are still interpolated exactly.
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 3; ++c) {
                expect(interp.evaluate(c, r)).toBeCloseTo(c * r, 12);
            }
        }
    });

    it('passes through the samples of a non-polynomial function', () => {
        const f = (x: number, y: number): number => Math.sin(x) * Math.exp(0.5 * y);
        const xBound = 6, yBound = 5;
        const xMin = -1, dx = 0.5, yMin = 0, dy = 0.4;
        const F = makeSamples(xBound, yBound, xMin, dx, yMin, dy, f);
        const interp = new IntpAkimaUniform2(xBound, yBound, xMin, dx, yMin, dy, F);

        for (let r = 0; r < yBound; ++r) {
            for (let c = 0; c < xBound; ++c) {
                const x = xMin + c * dx;
                const y = yMin + r * dy;
                expect(interp.evaluate(x, y)).toBeCloseTo(F[c + xBound * r], 10);
                expect(interp.evaluate(0, 0, x, y)).toBeCloseTo(F[c + xBound * r], 10);
            }
        }
    });

    it('exactly reproduces a per-variable quadratic away from the max boundaries', () => {
        // f(x,y) = (1 + 2x - x^2) * (3 - y + 0.5 y^2). Uniform slopes of a
        // quadratic are linear in the index, so the Akima weighted average
        // reduces to the exact central difference and the boundary
        // extrapolations are exact as well; the interpolant reproduces f in
        // every cell whose corners avoid the max-boundary fxy estimates
        // (see the upstream GetFXY sign note in the previous test).
        const g = (x: number): number => 1 + 2 * x - x * x;
        const dg = (x: number): number => 2 - 2 * x;
        const ddg = (): number => -2;
        const h = (y: number): number => 3 - y + 0.5 * y * y;
        const dh = (y: number): number => -1 + y;
        const ddh = (): number => 1;

        const xBound = 5, yBound = 6;
        const xMin = -0.5, dx = 0.5, yMin = 1, dy = 0.25;
        const F = makeSamples(xBound, yBound, xMin, dx, yMin, dy,
            (x, y) => g(x) * h(y));
        const interp = new IntpAkimaUniform2(xBound, yBound, xMin, dx, yMin, dy, F);

        for (const [x, y] of [[-0.3, 1.1], [0.42, 1.87], [0.9, 1.6], [0.75, 1.95]]) {
            expect(interp.evaluate(x, y)).toBeCloseTo(g(x) * h(y), 9);
            expect(interp.evaluate(1, 0, x, y)).toBeCloseTo(dg(x) * h(y), 9);
            expect(interp.evaluate(0, 1, x, y)).toBeCloseTo(g(x) * dh(y), 9);
            expect(interp.evaluate(1, 1, x, y)).toBeCloseTo(dg(x) * dh(y), 9);
            expect(interp.evaluate(2, 0, x, y)).toBeCloseTo(ddg() * h(y), 9);
            expect(interp.evaluate(0, 2, x, y)).toBeCloseTo(g(x) * ddh(), 9);
            expect(interp.evaluate(2, 2, x, y)).toBeCloseTo(ddg() * ddh(), 8);
        }
    });

    it('is C1 across interior cell boundaries', () => {
        const f = (x: number, y: number): number => Math.sin(x) * Math.exp(0.5 * y);
        const xBound = 6, yBound = 6;
        const xMin = 0, dx = 0.5, yMin = 0, dy = 0.5;
        const F = makeSamples(xBound, yBound, xMin, dx, yMin, dy, f);
        const interp = new IntpAkimaUniform2(xBound, yBound, xMin, dx, yMin, dy, F);

        const eps = 1e-7;
        // Cross the vertical grid line x = 1.0 and the horizontal grid line
        // y = 1.5 at several positions; value and first derivatives must
        // agree in the limit from both sides.
        for (const y of [0.3, 1.2, 2.1]) {
            for (const [xo, yo] of [[0, 0], [1, 0], [0, 1]]) {
                const left = interp.evaluate(xo, yo, 1 - eps, y);
                const right = interp.evaluate(xo, yo, 1 + eps, y);
                expect(right - left).toBeCloseTo(0, 5);
            }
        }
        for (const x of [0.4, 1.3, 2.2]) {
            for (const [xo, yo] of [[0, 0], [1, 0], [0, 1]]) {
                const below = interp.evaluate(xo, yo, x, 1.5 - eps);
                const above = interp.evaluate(xo, yo, x, 1.5 + eps);
                expect(above - below).toBeCloseTo(0, 5);
            }
        }
    });

    it('clamps evaluations to the domain', () => {
        const F = makeSamples(4, 4, 0, 1, 0, 1, (x, y) => x * x + y);
        const interp = new IntpAkimaUniform2(4, 4, 0, 1, 0, 1, F);
        expect(interp.evaluate(-10, 1.5)).toBeCloseTo(interp.evaluate(0, 1.5), 12);
        expect(interp.evaluate(10, 1.5)).toBeCloseTo(interp.evaluate(3, 1.5), 12);
        expect(interp.evaluate(1.5, -10)).toBeCloseTo(interp.evaluate(1.5, 0), 12);
        expect(interp.evaluate(1.5, 10)).toBeCloseTo(interp.evaluate(1.5, 3), 12);
    });

    it('returns zero for derivative orders beyond the degree', () => {
        const F = makeSamples(3, 3, 0, 1, 0, 1, (x, y) => x + y);
        const interp = new IntpAkimaUniform2(3, 3, 0, 1, 0, 1, F);
        expect(interp.evaluate(4, 0, 0.5, 0.5)).toBe(0);
        expect(interp.evaluate(0, 4, 0.5, 0.5)).toBe(0);
        expect(interp.evaluate(3, 0, 0.5, 0.5)).toBeCloseTo(0, 12);
    });
});
