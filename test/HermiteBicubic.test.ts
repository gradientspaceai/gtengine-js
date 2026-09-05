import { describe, it, expect } from 'vitest';
import { HermiteBicubic, HermiteBicubicSample } from '../src/HermiteBicubic.js';
import { HermiteCubic, HermiteCubicSample } from '../src/HermiteCubic.js';
import { check, expectClose, fc, finite, scaled } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V28): property-based cross-checks against HermiteBicubic.h.
// ---------------------------------------------------------------------------
describe('HermiteBicubic verification', () => {
    // A random polynomial of degree <= 3 in each variable, with exact
    // derivatives. deriv(a, k) differentiates the 1D coefficient list k times.
    const deriv = (a: number[], k: number): number[] => {
        let c = a.slice();
        for (let n = 0; n < k; ++n) {
            c = c.slice(1).map((ci, i) => (i + 1) * ci);
        }
        return c.length > 0 ? c : [0];
    };
    const horner = (c: number[], t: number): number =>
        c.reduceRight((acc, ci) => acc * t + ci, 0);

    // p(x,y) = sum_{i,j} a[i][j] x^i y^j; d(p)/dx^m dy^n at (x,y).
    const polyArb = (deg: number) =>
        fc.array(finite(-2, 2),
            { minLength: (deg + 1) * (deg + 1), maxLength: (deg + 1) * (deg + 1) });
    const evalPoly = (a: number[], deg: number, mx: number, my: number,
        x: number, y: number): number => {
        let sum = 0;
        for (let i = 0; i <= deg; ++i) {
            // Row i of the coefficient table as a polynomial in y.
            const row: number[] = [];
            for (let j = 0; j <= deg; ++j) { row.push(a[i * (deg + 1) + j]); }
            sum += horner(deriv(row, my), y)
                * horner(deriv(unit(i, deg), mx), x);
        }
        return sum;
    };
    // The 1D monomial x^i as a coefficient list of length deg+1.
    const unit = (i: number, deg: number): number[] => {
        const c = new Array<number>(deg + 1).fill(0);
        c[i] = 1;
        return c;
    };

    const sampleAt = (a: number[], bx: number, by: number): HermiteBicubicSample =>
        new HermiteBicubicSample(
            evalPoly(a, 3, 0, 0, bx, by),
            evalPoly(a, 3, 1, 0, bx, by),
            evalPoly(a, 3, 0, 1, bx, by),
            evalPoly(a, 3, 1, 1, bx, by));

    const blocksOf = (a: number[]) => [
        [sampleAt(a, 0, 0), sampleAt(a, 0, 1)],
        [sampleAt(a, 1, 0), sampleAt(a, 1, 1)]
    ] as const;

    it('reproduces every bicubic from exact samples of f, fx, fy, fxy', () => {
        check(fc.tuple(polyArb(3), scaled(0, 1), scaled(0, 1)), ([a, x, y]) => {
            const h = new HermiteBicubic(blocksOf(a));
            // generateSingle mixes the data with weights up to 9, so the
            // reconstruction keeps roughly 13 significant digits for
            // coefficients bounded by 2.
            expectClose(h.evaluate(0, 0, x, y), evalPoly(a, 3, 0, 0, x, y),
                1e-11, 1e-11);
            expectClose(h.evaluate(1, 0, x, y), evalPoly(a, 3, 1, 0, x, y),
                1e-10, 1e-10);
            expectClose(h.evaluate(0, 1, x, y), evalPoly(a, 3, 0, 1, x, y),
                1e-10, 1e-10);
            expectClose(h.evaluate(1, 1, x, y), evalPoly(a, 3, 1, 1, x, y),
                1e-10, 1e-10);
            expectClose(h.evaluate(2, 0, x, y), evalPoly(a, 3, 2, 0, x, y),
                1e-9, 1e-9);
            expectClose(h.evaluate(0, 2, x, y), evalPoly(a, 3, 0, 2, x, y),
                1e-9, 1e-9);
        });
    });

    it('interpolates the prescribed data at all four corners', () => {
        check(fc.array(finite(-4, 4), { minLength: 16, maxLength: 16 }), s => {
            const mk = (k: number) =>
                new HermiteBicubicSample(s[4 * k], s[4 * k + 1], s[4 * k + 2], s[4 * k + 3]);
            const blocks = [[mk(0), mk(1)], [mk(2), mk(3)]] as const;
            const h = new HermiteBicubic(blocks);
            for (let b0 = 0; b0 <= 1; ++b0) {
                for (let b1 = 0; b1 <= 1; ++b1) {
                    const b = blocks[b0][b1];
                    expectClose(h.evaluate(0, 0, b0, b1), b.f, 1e-12, 1e-12);
                    expectClose(h.evaluate(1, 0, b0, b1), b.fx, 1e-11, 1e-11);
                    expectClose(h.evaluate(0, 1, b0, b1), b.fy, 1e-11, 1e-11);
                    expectClose(h.evaluate(1, 1, b0, b1), b.fxy, 1e-10, 1e-10);
                }
            }
        });
    });

    it('agrees with the tensor product of two HermiteCubics on separable data', () => {
        check(fc.tuple(fc.array(finite(-3, 3), { minLength: 4, maxLength: 4 }),
            scaled(0, 1), scaled(0, 1)), ([s, x, y]) => {
            // f(x,y) = g(x) h(y) with g, h Hermite cubics; then
            // F = g h, Fx = g' h, Fy = g h', Fxy = g' h'.
            const g = new HermiteCubic([new HermiteCubicSample(s[0], s[1]),
                new HermiteCubicSample(s[2], s[3])]);
            const hh = new HermiteCubic([new HermiteCubicSample(s[3], s[0]),
                new HermiteCubicSample(s[1], s[2])]);
            const mk = (bx: number, by: number) => new HermiteBicubicSample(
                g.evaluate(0, bx) * hh.evaluate(0, by),
                g.evaluate(1, bx) * hh.evaluate(0, by),
                g.evaluate(0, bx) * hh.evaluate(1, by),
                g.evaluate(1, bx) * hh.evaluate(1, by));
            const b = new HermiteBicubic([[mk(0, 0), mk(0, 1)], [mk(1, 0), mk(1, 1)]]);
            expectClose(b.evaluate(0, 0, x, y),
                g.evaluate(0, x) * hh.evaluate(0, y), 1e-11, 1e-11);
            expectClose(b.evaluate(1, 1, x, y),
                g.evaluate(1, x) * hh.evaluate(1, y), 1e-10, 1e-10);
        });
    });

    it('exposes c publicly for manual coefficient assignment (upstream API)', () => {
        // Upstream declares 'c' public with the comment "Set the coefficients
        // manually as desired"; this is the regression guard for that.
        const h = new HermiteBicubic();
        h.c[2][1] = 3;
        // P(2,x) P(1,y) = (1-x) x^2 * (1-y)^2 y, so H(1/2,1/2) = 3/8 * 1/8.
        expectClose(h.evaluate(0, 0, 0.5, 0.5), 3 * 0.125 * 0.125, 1e-15, 1e-15);
    });

    it('returns zero when either order exceeds the degree', () => {
        check(fc.tuple(fc.integer({ min: 4, max: 20 }), fc.integer({ min: 0, max: 3 }),
            scaled(0, 1), scaled(0, 1)), ([big, small, x, y]) => {
            const h = new HermiteBicubic([
                [new HermiteBicubicSample(1, 2, 3, 4), new HermiteBicubicSample(5, 6, 7, 8)],
                [new HermiteBicubicSample(-1, -2, -3, -4), new HermiteBicubicSample(9, 8, 7, 6)]
            ]);
            expect(h.evaluate(big, small, x, y)).toBe(0);
            expect(h.evaluate(small, big, x, y)).toBe(0);
        });
    });
});
