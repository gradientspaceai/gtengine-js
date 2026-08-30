import { describe, it, expect } from 'vitest';
import { HermiteTriquintic, HermiteTriquinticSample } from '../src/HermiteTriquintic';

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
// a,b,c in {0,1,2} is named 'f' + 'x'*a + 'y'*b + 'z'*c.
function fieldName(a: number, b: number, c: number): string {
    return 'f' + 'x'.repeat(a) + 'y'.repeat(b) + 'z'.repeat(c);
}

function getField(s: HermiteTriquinticSample, a: number, b: number, c: number): number {
    return (s as unknown as Record<string, number>)[fieldName(a, b, c)];
}

function setField(s: HermiteTriquinticSample, a: number, b: number, c: number,
    value: number): void {
    (s as unknown as Record<string, number>)[fieldName(a, b, c)] = value;
}

type Pair = readonly [HermiteTriquinticSample, HermiteTriquinticSample];
type Blocks = readonly [readonly [Pair, Pair], readonly [Pair, Pair]];

function makeBlocks(
    corner: (b0: number, b1: number, b2: number) => HermiteTriquinticSample): Blocks {
    return [
        [[corner(0, 0, 0), corner(0, 0, 1)], [corner(0, 1, 0), corner(0, 1, 1)]],
        [[corner(1, 0, 0), corner(1, 0, 1)], [corner(1, 1, 0), corner(1, 1, 1)]]
    ] as const;
}

describe('HermiteTriquintic', () => {
    it('default constructor creates the identically zero polynomial', () => {
        const h = new HermiteTriquintic();
        for (const t of [0, 0.4, 1]) {
            for (let a = 0; a <= 5; ++a) {
                for (let b = 0; b <= 5; ++b) {
                    for (let c = 0; c <= 5; ++c) {
                        expect(h.evaluate(a, b, c, t, 1 - t, 0.5 * t)).toBe(0);
                    }
                }
            }
        }
    });

    it('sample default constructor zero-fills all 27 fields', () => {
        const s = new HermiteTriquinticSample();
        for (let a = 0; a <= 2; ++a) {
            for (let b = 0; b <= 2; ++b) {
                for (let c = 0; c <= 2; ++c) {
                    expect(getField(s, a, b, c)).toBe(0);
                }
            }
        }
    });

    it('sample constructor assigns all 27 fields in upstream order', () => {
        // Argument order must match the upstream Sample constructor.
        const s = new HermiteTriquinticSample(
            1,                              // f
            2, 3, 4,                        // fx, fy, fz
            5, 6, 7, 8, 9, 10,              // fxx, fxy, fxz, fyy, fyz, fzz
            11, 12, 13, 14, 15, 16, 17,     // fxxy, fxxz, fxyy, fxyz, fxzz, fyyz, fyzz
            18, 19, 20, 21, 22, 23,         // fxxyy, fxxyz, fxxzz, fxyyz, fxyzz, fyyzz
            24, 25, 26,                     // fxxyyz, fxxyzz, fxyyzz
            27);                            // fxxyyzz
        expect(s.f).toBe(1);
        expect(s.fx).toBe(2);
        expect(s.fy).toBe(3);
        expect(s.fz).toBe(4);
        expect(s.fxx).toBe(5);
        expect(s.fxy).toBe(6);
        expect(s.fxz).toBe(7);
        expect(s.fyy).toBe(8);
        expect(s.fyz).toBe(9);
        expect(s.fzz).toBe(10);
        expect(s.fxxy).toBe(11);
        expect(s.fxxz).toBe(12);
        expect(s.fxyy).toBe(13);
        expect(s.fxyz).toBe(14);
        expect(s.fxzz).toBe(15);
        expect(s.fyyz).toBe(16);
        expect(s.fyzz).toBe(17);
        expect(s.fxxyy).toBe(18);
        expect(s.fxxyz).toBe(19);
        expect(s.fxxzz).toBe(20);
        expect(s.fxyyz).toBe(21);
        expect(s.fxyzz).toBe(22);
        expect(s.fyyzz).toBe(23);
        expect(s.fxxyyz).toBe(24);
        expect(s.fxxyzz).toBe(25);
        expect(s.fxyyzz).toBe(26);
        expect(s.fxxyyzz).toBe(27);
    });

    it('reproduces the samples (values and mixed partials) at the corners', () => {
        const value = (b0: number, b1: number, b2: number,
            a: number, b: number, c: number): number =>
            ((7 * b0 + 13 * b1 + 17 * b2 + 3 * a + 5 * b + 11 * c + 1) % 13) - 6
            + 0.25 * a * b * c;

        const corner = (b0: number, b1: number, b2: number): HermiteTriquinticSample => {
            const s = new HermiteTriquinticSample();
            for (let a = 0; a <= 2; ++a) {
                for (let b = 0; b <= 2; ++b) {
                    for (let c = 0; c <= 2; ++c) {
                        setField(s, a, b, c, value(b0, b1, b2, a, b, c));
                    }
                }
            }
            return s;
        };

        const h = new HermiteTriquintic(makeBlocks(corner));

        for (let b0 = 0; b0 <= 1; ++b0) {
            for (let b1 = 0; b1 <= 1; ++b1) {
                for (let b2 = 0; b2 <= 1; ++b2) {
                    for (let a = 0; a <= 2; ++a) {
                        for (let b = 0; b <= 2; ++b) {
                            for (let c = 0; c <= 2; ++c) {
                                expect(h.evaluate(a, b, c, b0, b1, b2))
                                    .toBeCloseTo(value(b0, b1, b2, a, b, c), 9);
                            }
                        }
                    }
                }
            }
        }
    });

    it('exactly reproduces a tensor-product quintic polynomial', () => {
        // g(x,y,z) = X(x) * Y(y) * Z(z) with quintic factors; the
        // triquintic Hermite interpolant reproduces g on [0,1]^3 including
        // partial derivatives.
        const X = derivatives([1, -2, 3, -1, 0.5, 2], 5);
        const Y = derivatives([2, 1, -1, 0.5, -0.25, 1], 5);
        const Z = derivatives([-1, 0.5, 2, 1, -0.5, -1], 5);

        const corner = (x: number, y: number, z: number): HermiteTriquinticSample => {
            const s = new HermiteTriquinticSample();
            for (let a = 0; a <= 2; ++a) {
                for (let b = 0; b <= 2; ++b) {
                    for (let c = 0; c <= 2; ++c) {
                        setField(s, a, b, c,
                            polyval(X[a], x) * polyval(Y[b], y) * polyval(Z[c], z));
                    }
                }
            }
            return s;
        };

        const h = new HermiteTriquintic(makeBlocks(corner));

        for (let a = 0; a <= 5; ++a) {
            for (let b = 0; b <= 5; ++b) {
                for (let c = 0; c <= 5; ++c) {
                    for (const [x, y, z] of
                        [[0, 0, 0], [0.25, 0.5, 0.75], [0.5, 0.5, 0.5], [1, 0.25, 0.75]]) {
                        const expected =
                            polyval(X[a], x) * polyval(Y[b], y) * polyval(Z[c], z);
                        expect(h.evaluate(a, b, c, x, y, z)).toBeCloseTo(expected, 7);
                    }
                }
            }
        }
    });

    it('returns zero for derivative orders beyond the degree', () => {
        const corner = (): HermiteTriquinticSample =>
            new HermiteTriquinticSample(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
                13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27);
        const h = new HermiteTriquintic(makeBlocks(corner));
        expect(h.evaluate(6, 0, 0, 0.5, 0.5, 0.5)).toBe(0);
        expect(h.evaluate(0, 6, 0, 0.5, 0.5, 0.5)).toBe(0);
        expect(h.evaluate(0, 0, 6, 0.5, 0.5, 0.5)).toBe(0);
        expect(h.evaluate(7, 7, 7, 0.25, 0.75, 0.5)).toBe(0);
    });
});
