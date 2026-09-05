import { describe, it, expect } from 'vitest';
import { HermiteTricubic, HermiteTricubicSample } from '../src/HermiteTricubic.js';
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

// ---------------------------------------------------------------------------
// Verification (V28): property-based cross-checks against HermiteTricubic.h.
// ---------------------------------------------------------------------------
describe('HermiteTricubic verification', () => {
    const mono = (k: number, i: number, t: number): number => {
        if (i < k) { return 0; }
        let c = 1;
        for (let m = 0; m < k; ++m) { c *= i - m; }
        return c * Math.pow(t, i - k);
    };
    // p(x,y,z) = sum_{i,j,k<=3} a[16i+4j+k] x^i y^j z^k, differentiated
    // (mx, my, mz) times.
    const evalPoly = (a: number[], mx: number, my: number, mz: number,
        x: number, y: number, z: number): number => {
        let sum = 0;
        for (let i = 0; i <= 3; ++i) {
            const mi = mono(mx, i, x);
            if (mi === 0) { continue; }
            for (let j = 0; j <= 3; ++j) {
                const mj = mono(my, j, y);
                if (mj === 0) { continue; }
                for (let k = 0; k <= 3; ++k) {
                    sum += a[16 * i + 4 * j + k] * mi * mj * mono(mz, k, z);
                }
            }
        }
        return sum;
    };
    const sampleAt = (a: number[], bx: number, by: number, bz: number) =>
        new HermiteTricubicSample(
            evalPoly(a, 0, 0, 0, bx, by, bz),
            evalPoly(a, 1, 0, 0, bx, by, bz),
            evalPoly(a, 0, 1, 0, bx, by, bz),
            evalPoly(a, 0, 0, 1, bx, by, bz),
            evalPoly(a, 1, 1, 0, bx, by, bz),
            evalPoly(a, 1, 0, 1, bx, by, bz),
            evalPoly(a, 0, 1, 1, bx, by, bz),
            evalPoly(a, 1, 1, 1, bx, by, bz));
    const blocksOf = (a: number[]) => [
        [[sampleAt(a, 0, 0, 0), sampleAt(a, 0, 0, 1)],
            [sampleAt(a, 0, 1, 0), sampleAt(a, 0, 1, 1)]],
        [[sampleAt(a, 1, 0, 0), sampleAt(a, 1, 0, 1)],
            [sampleAt(a, 1, 1, 0), sampleAt(a, 1, 1, 1)]]
    ] as const;

    it('reproduces every tricubic from exact samples of its eight derivatives', () => {
        check(fc.tuple(fc.array(finite(-1, 1), { minLength: 64, maxLength: 64 }),
            scaled(0, 1), scaled(0, 1), scaled(0, 1)), ([a, x, y, z]) => {
            const h = new HermiteTricubic(blocksOf(a));
            // The coefficient weights reach 27 and the corner data reach ~64
            // for unit coefficients, so about three digits are lost.
            for (const [mx, my, mz] of [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
                [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1], [2, 0, 0]]) {
                expectClose(h.evaluate(mx, my, mz, x, y, z),
                    evalPoly(a, mx, my, mz, x, y, z), 1e-10, 1e-10);
            }
        });
    });

    it('interpolates the prescribed data at all eight corners', () => {
        check(fc.array(finite(-4, 4), { minLength: 64, maxLength: 64 }), s => {
            const mk = (k: number) => new HermiteTricubicSample(
                s[8 * k], s[8 * k + 1], s[8 * k + 2], s[8 * k + 3],
                s[8 * k + 4], s[8 * k + 5], s[8 * k + 6], s[8 * k + 7]);
            const blocks = [[[mk(0), mk(1)], [mk(2), mk(3)]],
                [[mk(4), mk(5)], [mk(6), mk(7)]]] as const;
            const h = new HermiteTricubic(blocks);
            for (let b0 = 0; b0 <= 1; ++b0) {
                for (let b1 = 0; b1 <= 1; ++b1) {
                    for (let b2 = 0; b2 <= 1; ++b2) {
                        const b = blocks[b0][b1][b2];
                        expectClose(h.evaluate(0, 0, 0, b0, b1, b2), b.f, 1e-11, 1e-11);
                        expectClose(h.evaluate(1, 0, 0, b0, b1, b2), b.fx, 1e-10, 1e-10);
                        expectClose(h.evaluate(0, 1, 0, b0, b1, b2), b.fy, 1e-10, 1e-10);
                        expectClose(h.evaluate(0, 0, 1, b0, b1, b2), b.fz, 1e-10, 1e-10);
                        expectClose(h.evaluate(1, 1, 1, b0, b1, b2), b.fxyz, 1e-9, 1e-9);
                    }
                }
            }
        });
    });

    it('agrees with the tensor product of three HermiteCubics', () => {
        check(fc.tuple(fc.array(finite(-2, 2), { minLength: 12, maxLength: 12 }),
            scaled(0, 1), scaled(0, 1), scaled(0, 1)), ([s, x, y, z]) => {
            const mk1 = (o: number) => new HermiteCubic([
                new HermiteCubicSample(s[o], s[o + 1]),
                new HermiteCubicSample(s[o + 2], s[o + 3])]);
            const g = mk1(0), k = mk1(4), l = mk1(8);
            const mk = (bx: number, by: number, bz: number) =>
                new HermiteTricubicSample(
                    g.evaluate(0, bx) * k.evaluate(0, by) * l.evaluate(0, bz),
                    g.evaluate(1, bx) * k.evaluate(0, by) * l.evaluate(0, bz),
                    g.evaluate(0, bx) * k.evaluate(1, by) * l.evaluate(0, bz),
                    g.evaluate(0, bx) * k.evaluate(0, by) * l.evaluate(1, bz),
                    g.evaluate(1, bx) * k.evaluate(1, by) * l.evaluate(0, bz),
                    g.evaluate(1, bx) * k.evaluate(0, by) * l.evaluate(1, bz),
                    g.evaluate(0, bx) * k.evaluate(1, by) * l.evaluate(1, bz),
                    g.evaluate(1, bx) * k.evaluate(1, by) * l.evaluate(1, bz));
            const h = new HermiteTricubic([
                [[mk(0, 0, 0), mk(0, 0, 1)], [mk(0, 1, 0), mk(0, 1, 1)]],
                [[mk(1, 0, 0), mk(1, 0, 1)], [mk(1, 1, 0), mk(1, 1, 1)]]]);
            expectClose(h.evaluate(0, 0, 0, x, y, z),
                g.evaluate(0, x) * k.evaluate(0, y) * l.evaluate(0, z), 1e-10, 1e-10);
            expectClose(h.evaluate(1, 1, 1, x, y, z),
                g.evaluate(1, x) * k.evaluate(1, y) * l.evaluate(1, z), 1e-9, 1e-9);
        });
    });

    it('exposes c publicly for manual coefficient assignment (upstream API)', () => {
        const h = new HermiteTricubic();
        h.c[3][3][3] = 5;
        // P(3,t) = t^3, so H(x,y,z) = 5 x^3 y^3 z^3.
        expectClose(h.evaluate(0, 0, 0, 1, 1, 1), 5, 1e-15, 1e-15);
        expectClose(h.evaluate(0, 0, 0, 0.5, 0.5, 0.5), 5 / 512, 1e-15, 1e-15);
    });

    it('returns zero when any order exceeds the degree', () => {
        check(fc.tuple(fc.integer({ min: 4, max: 20 }), scaled(0, 1)), ([big, t]) => {
            const h = new HermiteTricubic();
            h.c[1][2][0] = 1;
            expect(h.evaluate(big, 0, 0, t, t, t)).toBe(0);
            expect(h.evaluate(0, big, 0, t, t, t)).toBe(0);
            expect(h.evaluate(0, 0, big, t, t, t)).toBe(0);
        });
    });
});
