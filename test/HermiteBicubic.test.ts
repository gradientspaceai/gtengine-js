import { describe, it, expect } from 'vitest';
import { HermiteBicubic, HermiteBicubicSample } from '../src/HermiteBicubic.js';

// Polynomial helpers: value and derivative of sum_i c[i] x^i.
function polyval(c: readonly number[], x: number): number {
    let result = 0;
    for (let i = c.length - 1; i >= 0; --i) {
        result = result * x + c[i];
    }
    return result;
}

function polyder(c: readonly number[]): number[] {
    const d: number[] = [];
    for (let i = 1; i < c.length; ++i) {
        d.push(i * c[i]);
    }
    return d;
}

// All derivatives X^{(a)} for a = 0..maxOrder as coefficient arrays.
function derivatives(c: readonly number[], maxOrder: number): number[][] {
    const all: number[][] = [c.slice()];
    for (let a = 1; a <= maxOrder; ++a) {
        all.push(polyder(all[a - 1]));
    }
    return all;
}

describe('HermiteBicubic', () => {
    it('default constructor creates the identically zero polynomial', () => {
        const h = new HermiteBicubic();
        for (const x of [0, 0.25, 0.75, 1]) {
            for (const y of [0, 0.5, 1]) {
                for (let a = 0; a <= 3; ++a) {
                    for (let b = 0; b <= 3; ++b) {
                        expect(h.evaluate(a, b, x, y)).toBe(0);
                    }
                }
            }
        }
    });

    it('sample default constructor zero-fills', () => {
        const s = new HermiteBicubicSample();
        expect(s.f).toBe(0);
        expect(s.fx).toBe(0);
        expect(s.fy).toBe(0);
        expect(s.fxy).toBe(0);
    });

    it('reproduces the samples (values and mixed partials) at the corners', () => {
        // Arbitrary distinct sample data at the four corners.
        const blocks = [
            [new HermiteBicubicSample(1, 2, -1, 0.5),
             new HermiteBicubicSample(-2, 0.25, 3, -1.5)],
            [new HermiteBicubicSample(4, -3, 1.25, 2),
             new HermiteBicubicSample(0.5, 1, -0.75, -2.5)]
        ] as const;
        const h = new HermiteBicubic(blocks);

        for (let b0 = 0; b0 <= 1; ++b0) {
            for (let b1 = 0; b1 <= 1; ++b1) {
                const s = blocks[b0][b1];
                expect(h.evaluate(0, 0, b0, b1)).toBeCloseTo(s.f, 12);
                expect(h.evaluate(1, 0, b0, b1)).toBeCloseTo(s.fx, 12);
                expect(h.evaluate(0, 1, b0, b1)).toBeCloseTo(s.fy, 12);
                expect(h.evaluate(1, 1, b0, b1)).toBeCloseTo(s.fxy, 12);
            }
        }
    });

    it('exactly reproduces a tensor-product cubic polynomial', () => {
        // g(x,y) = X(x) * Y(y) with cubic factors. The bicubic Hermite
        // interpolant matches g and its relevant partials at the corners, and
        // the tensor-product cubic space is reproduced exactly, so H = g on
        // [0,1]^2 including all partial derivatives.
        const X = derivatives([1, -2, 3, -1], 3);
        const Y = derivatives([2, 1, -1, 0.5], 3);

        const corner = (x: number, y: number): HermiteBicubicSample =>
            new HermiteBicubicSample(
                polyval(X[0], x) * polyval(Y[0], y),
                polyval(X[1], x) * polyval(Y[0], y),
                polyval(X[0], x) * polyval(Y[1], y),
                polyval(X[1], x) * polyval(Y[1], y));

        const h = new HermiteBicubic([
            [corner(0, 0), corner(0, 1)],
            [corner(1, 0), corner(1, 1)]
        ]);

        for (let a = 0; a <= 3; ++a) {
            for (let b = 0; b <= 3; ++b) {
                for (let i = 0; i <= 4; ++i) {
                    for (let j = 0; j <= 4; ++j) {
                        const x = i / 4;
                        const y = j / 4;
                        const expected = polyval(X[a], x) * polyval(Y[b], y);
                        expect(h.evaluate(a, b, x, y)).toBeCloseTo(expected, 10);
                    }
                }
            }
        }
    });

    it('mixed partials agree with central finite differences', () => {
        const h = new HermiteBicubic([
            [new HermiteBicubicSample(1, 2, -1, 0.5),
             new HermiteBicubicSample(-2, 0.25, 3, -1.5)],
            [new HermiteBicubicSample(4, -3, 1.25, 2),
             new HermiteBicubicSample(0.5, 1, -0.75, -2.5)]
        ]);
        const step = 1e-5;
        for (const [x, y] of [[0.3, 0.6], [0.75, 0.2]]) {
            const fdX = (h.evaluate(0, 0, x + step, y) - h.evaluate(0, 0, x - step, y)) / (2 * step);
            expect(h.evaluate(1, 0, x, y)).toBeCloseTo(fdX, 5);
            const fdY = (h.evaluate(0, 0, x, y + step) - h.evaluate(0, 0, x, y - step)) / (2 * step);
            expect(h.evaluate(0, 1, x, y)).toBeCloseTo(fdY, 5);
            const fdXY = (h.evaluate(1, 0, x, y + step) - h.evaluate(1, 0, x, y - step)) / (2 * step);
            expect(h.evaluate(1, 1, x, y)).toBeCloseTo(fdXY, 5);
        }
    });

    it('returns zero for derivative orders beyond the degree', () => {
        const h = new HermiteBicubic([
            [new HermiteBicubicSample(1, 2, -1, 0.5),
             new HermiteBicubicSample(-2, 0.25, 3, -1.5)],
            [new HermiteBicubicSample(4, -3, 1.25, 2),
             new HermiteBicubicSample(0.5, 1, -0.75, -2.5)]
        ]);
        expect(h.evaluate(4, 0, 0.5, 0.5)).toBe(0);
        expect(h.evaluate(0, 4, 0.5, 0.5)).toBe(0);
        expect(h.evaluate(5, 5, 0.25, 0.75)).toBe(0);
    });
});
