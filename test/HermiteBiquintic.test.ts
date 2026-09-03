import { describe, it, expect } from 'vitest';
import { HermiteBiquintic, HermiteBiquinticSample } from '../src/HermiteBiquintic.js';

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

function derivatives(c: readonly number[], maxOrder: number): number[][] {
    const all: number[][] = [c.slice()];
    for (let a = 1; a <= maxOrder; ++a) {
        all.push(polyder(all[a - 1]));
    }
    return all;
}

// The sample field holding the partial derivative of orders (a,b) is named
// 'f' + 'x' repeated a times + 'y' repeated b times, e.g. (2,1) -> 'fxxy'.
function fieldName(a: number, b: number): string {
    return 'f' + 'x'.repeat(a) + 'y'.repeat(b);
}

function getField(s: HermiteBiquinticSample, a: number, b: number): number {
    return (s as unknown as Record<string, number>)[fieldName(a, b)];
}

function setField(s: HermiteBiquinticSample, a: number, b: number, value: number): void {
    (s as unknown as Record<string, number>)[fieldName(a, b)] = value;
}

describe('HermiteBiquintic', () => {
    it('default constructor creates the identically zero polynomial', () => {
        const h = new HermiteBiquintic();
        for (const x of [0, 0.25, 0.75, 1]) {
            for (const y of [0, 0.5, 1]) {
                for (let a = 0; a <= 5; ++a) {
                    for (let b = 0; b <= 5; ++b) {
                        expect(h.evaluate(a, b, x, y)).toBe(0);
                    }
                }
            }
        }
    });

    it('sample default constructor zero-fills', () => {
        const s = new HermiteBiquinticSample();
        for (let a = 0; a <= 2; ++a) {
            for (let b = 0; b <= 2; ++b) {
                expect(getField(s, a, b)).toBe(0);
            }
        }
    });

    it('reproduces the samples (values and mixed partials) at the corners', () => {
        // Arbitrary distinct sample data at the four corners: field for
        // orders (a,b) at corner (b0,b1) gets a deterministic unique value.
        const value = (b0: number, b1: number, a: number, b: number): number =>
            ((7 * b0 + 13 * b1 + 3 * a + 5 * b + 1) % 11) - 5 + 0.5 * a * b;

        const makeCorner = (b0: number, b1: number): HermiteBiquinticSample => {
            const s = new HermiteBiquinticSample();
            for (let a = 0; a <= 2; ++a) {
                for (let b = 0; b <= 2; ++b) {
                    setField(s, a, b, value(b0, b1, a, b));
                }
            }
            return s;
        };

        const blocks = [
            [makeCorner(0, 0), makeCorner(0, 1)],
            [makeCorner(1, 0), makeCorner(1, 1)]
        ] as const;
        const h = new HermiteBiquintic(blocks);

        for (let b0 = 0; b0 <= 1; ++b0) {
            for (let b1 = 0; b1 <= 1; ++b1) {
                for (let a = 0; a <= 2; ++a) {
                    for (let b = 0; b <= 2; ++b) {
                        expect(h.evaluate(a, b, b0, b1))
                            .toBeCloseTo(value(b0, b1, a, b), 10);
                    }
                }
            }
        }
    });

    it('exactly reproduces a tensor-product quintic polynomial', () => {
        // g(x,y) = X(x) * Y(y) with quintic factors. The biquintic Hermite
        // interpolant matches g and all partials through order (2,2) at the
        // corners, so H = g on [0,1]^2 including partial derivatives.
        const X = derivatives([1, -2, 3, -1, 0.5, 2], 5);
        const Y = derivatives([2, 1, -1, 0.5, -0.25, 1], 5);

        const corner = (x: number, y: number): HermiteBiquinticSample => {
            const s = new HermiteBiquinticSample();
            for (let a = 0; a <= 2; ++a) {
                for (let b = 0; b <= 2; ++b) {
                    setField(s, a, b, polyval(X[a], x) * polyval(Y[b], y));
                }
            }
            return s;
        };

        const h = new HermiteBiquintic([
            [corner(0, 0), corner(0, 1)],
            [corner(1, 0), corner(1, 1)]
        ]);

        for (let a = 0; a <= 5; ++a) {
            for (let b = 0; b <= 5; ++b) {
                for (let i = 0; i <= 4; ++i) {
                    for (let j = 0; j <= 4; ++j) {
                        const x = i / 4;
                        const y = j / 4;
                        const expected = polyval(X[a], x) * polyval(Y[b], y);
                        expect(h.evaluate(a, b, x, y)).toBeCloseTo(expected, 8);
                    }
                }
            }
        }
    });

    it('returns zero for derivative orders beyond the degree', () => {
        const h = new HermiteBiquintic([
            [new HermiteBiquinticSample(1, 2, 3, 4, 5, 6, 7, 8, 9),
             new HermiteBiquinticSample(9, 8, 7, 6, 5, 4, 3, 2, 1)],
            [new HermiteBiquinticSample(-1, 2, -3, 4, -5, 6, -7, 8, -9),
             new HermiteBiquinticSample(2, -2, 2, -2, 2, -2, 2, -2, 2)]
        ]);
        expect(h.evaluate(6, 0, 0.5, 0.5)).toBe(0);
        expect(h.evaluate(0, 6, 0.5, 0.5)).toBe(0);
        expect(h.evaluate(7, 7, 0.25, 0.75)).toBe(0);
    });
});
