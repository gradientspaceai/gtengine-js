import { describe, it, expect } from 'vitest';
import { HermiteTricubic, HermiteTricubicSample } from '../src/HermiteTricubic.js';

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

// The sample field holding the partial derivative of orders (a,b,c) with
// a,b,c in {0,1} is named 'f' + 'x'*a + 'y'*b + 'z'*c.
function fieldName(a: number, b: number, c: number): string {
    return 'f' + 'x'.repeat(a) + 'y'.repeat(b) + 'z'.repeat(c);
}

function getField(s: HermiteTricubicSample, a: number, b: number, c: number): number {
    return (s as unknown as Record<string, number>)[fieldName(a, b, c)];
}

function setField(s: HermiteTricubicSample, a: number, b: number, c: number,
    value: number): void {
    (s as unknown as Record<string, number>)[fieldName(a, b, c)] = value;
}

type Pair = readonly [HermiteTricubicSample, HermiteTricubicSample];
type Blocks = readonly [readonly [Pair, Pair], readonly [Pair, Pair]];

function makeBlocks(
    corner: (b0: number, b1: number, b2: number) => HermiteTricubicSample): Blocks {
    return [
        [[corner(0, 0, 0), corner(0, 0, 1)], [corner(0, 1, 0), corner(0, 1, 1)]],
        [[corner(1, 0, 0), corner(1, 0, 1)], [corner(1, 1, 0), corner(1, 1, 1)]]
    ] as const;
}

describe('HermiteTricubic', () => {
    it('default constructor creates the identically zero polynomial', () => {
        const h = new HermiteTricubic();
        for (const t of [0, 0.4, 1]) {
            for (let a = 0; a <= 3; ++a) {
                for (let b = 0; b <= 3; ++b) {
                    for (let c = 0; c <= 3; ++c) {
                        expect(h.evaluate(a, b, c, t, 1 - t, 0.5 * t)).toBe(0);
                    }
                }
            }
        }
    });

    it('sample default constructor zero-fills', () => {
        const s = new HermiteTricubicSample();
        for (let a = 0; a <= 1; ++a) {
            for (let b = 0; b <= 1; ++b) {
                for (let c = 0; c <= 1; ++c) {
                    expect(getField(s, a, b, c)).toBe(0);
                }
            }
        }
    });

    it('reproduces the samples (values and mixed partials) at the corners', () => {
        const value = (b0: number, b1: number, b2: number,
            a: number, b: number, c: number): number =>
            ((7 * b0 + 13 * b1 + 17 * b2 + 3 * a + 5 * b + 11 * c + 1) % 13) - 6;

        const corner = (b0: number, b1: number, b2: number): HermiteTricubicSample => {
            const s = new HermiteTricubicSample();
            for (let a = 0; a <= 1; ++a) {
                for (let b = 0; b <= 1; ++b) {
                    for (let c = 0; c <= 1; ++c) {
                        setField(s, a, b, c, value(b0, b1, b2, a, b, c));
                    }
                }
            }
            return s;
        };

        const h = new HermiteTricubic(makeBlocks(corner));

        for (let b0 = 0; b0 <= 1; ++b0) {
            for (let b1 = 0; b1 <= 1; ++b1) {
                for (let b2 = 0; b2 <= 1; ++b2) {
                    for (let a = 0; a <= 1; ++a) {
                        for (let b = 0; b <= 1; ++b) {
                            for (let c = 0; c <= 1; ++c) {
                                expect(h.evaluate(a, b, c, b0, b1, b2))
                                    .toBeCloseTo(value(b0, b1, b2, a, b, c), 10);
                            }
                        }
                    }
                }
            }
        }
    });

    it('exactly reproduces a tensor-product cubic polynomial', () => {
        // g(x,y,z) = X(x) * Y(y) * Z(z) with cubic factors; the tricubic
        // Hermite interpolant reproduces g on [0,1]^3 including partials.
        const X = derivatives([1, -2, 3, -1], 3);
        const Y = derivatives([2, 1, -1, 0.5], 3);
        const Z = derivatives([-1, 0.5, 2, 1], 3);

        const corner = (x: number, y: number, z: number): HermiteTricubicSample => {
            const s = new HermiteTricubicSample();
            for (let a = 0; a <= 1; ++a) {
                for (let b = 0; b <= 1; ++b) {
                    for (let c = 0; c <= 1; ++c) {
                        setField(s, a, b, c,
                            polyval(X[a], x) * polyval(Y[b], y) * polyval(Z[c], z));
                    }
                }
            }
            return s;
        };

        const h = new HermiteTricubic(makeBlocks(corner));

        for (let a = 0; a <= 3; ++a) {
            for (let b = 0; b <= 3; ++b) {
                for (let c = 0; c <= 3; ++c) {
                    for (const [x, y, z] of
                        [[0, 0, 0], [0.25, 0.5, 0.75], [0.5, 0.5, 0.5],
                         [1, 0.25, 0.75], [1, 1, 1]]) {
                        const expected =
                            polyval(X[a], x) * polyval(Y[b], y) * polyval(Z[c], z);
                        expect(h.evaluate(a, b, c, x, y, z)).toBeCloseTo(expected, 9);
                    }
                }
            }
        }
    });

    it('returns zero for derivative orders beyond the degree', () => {
        const corner = (): HermiteTricubicSample =>
            new HermiteTricubicSample(1, 2, 3, 4, 5, 6, 7, 8);
        const h = new HermiteTricubic(makeBlocks(corner));
        expect(h.evaluate(4, 0, 0, 0.5, 0.5, 0.5)).toBe(0);
        expect(h.evaluate(0, 4, 0, 0.5, 0.5, 0.5)).toBe(0);
        expect(h.evaluate(0, 0, 4, 0.5, 0.5, 0.5)).toBe(0);
        expect(h.evaluate(5, 5, 5, 0.25, 0.75, 0.5)).toBe(0);
    });
});
