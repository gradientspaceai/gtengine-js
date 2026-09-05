import { describe, it, expect } from 'vitest';
import { HermiteTriquintic, HermiteTriquinticSample } from '../src/HermiteTriquintic.js';
import { HermiteQuintic, HermiteQuinticSample } from '../src/HermiteQuintic.js';
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

// ---------------------------------------------------------------------------
// Verification (V28): property-based cross-checks against HermiteTriquintic.h.
// ---------------------------------------------------------------------------
type Tuple27 = [number, number, number, number, number, number, number, number,
    number, number, number, number, number, number, number, number, number,
    number, number, number, number, number, number, number, number, number,
    number];

describe('HermiteTriquintic verification', () => {
    // The 27 sample fields in constructor order, with the (x,y,z) derivative
    // orders each one denotes. This list is also the vXYZ naming used by
    // upstream's GenerateSingle, so it doubles as a check of that mapping.
    const ORDERS: ReadonlyArray<readonly [number, number, number]> = [
        [0, 0, 0],
        [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2],
        [2, 1, 0], [2, 0, 1], [1, 2, 0], [1, 1, 1], [1, 0, 2], [0, 2, 1], [0, 1, 2],
        [2, 2, 0], [2, 1, 1], [2, 0, 2], [1, 2, 1], [1, 1, 2], [0, 2, 2],
        [2, 2, 1], [2, 1, 2], [1, 2, 2],
        [2, 2, 2]
    ];

    const FIELDS: ReadonlyArray<keyof HermiteTriquinticSample> = [
        'f', 'fx', 'fy', 'fz', 'fxx', 'fxy', 'fxz', 'fyy', 'fyz', 'fzz',
        'fxxy', 'fxxz', 'fxyy', 'fxyz', 'fxzz', 'fyyz', 'fyzz',
        'fxxyy', 'fxxyz', 'fxxzz', 'fxyyz', 'fxyzz', 'fyyzz',
        'fxxyyz', 'fxxyzz', 'fxyyzz', 'fxxyyzz'];

    const mono = (k: number, i: number, t: number): number => {
        if (i < k) { return 0; }
        let c = 1;
        for (let m = 0; m < k; ++m) { c *= i - m; }
        return c * Math.pow(t, i - k);
    };

    // p(x,y,z) = sum_{i,j,k<=5} a[36i+6j+k] x^i y^j z^k, differentiated m times.
    const evalPoly = (a: number[], m: readonly [number, number, number],
        x: number, y: number, z: number): number => {
        let sum = 0;
        for (let i = 0; i <= 5; ++i) {
            const mi = mono(m[0], i, x);
            if (mi === 0) { continue; }
            for (let j = 0; j <= 5; ++j) {
                const mj = mono(m[1], j, y);
                if (mj === 0) { continue; }
                for (let k = 0; k <= 5; ++k) {
                    sum += a[36 * i + 6 * j + k] * mi * mj * mono(m[2], k, z);
                }
            }
        }
        return sum;
    };

    const sampleOf = (values: number[]) =>
        new HermiteTriquinticSample(...(values as Tuple27));

    const sampleAt = (a: number[], bx: number, by: number, bz: number) =>
        sampleOf(ORDERS.map(m => evalPoly(a, m, bx, by, bz)));

    const blocksOf = (a: number[]) => [
        [[sampleAt(a, 0, 0, 0), sampleAt(a, 0, 0, 1)],
            [sampleAt(a, 0, 1, 0), sampleAt(a, 0, 1, 1)]],
        [[sampleAt(a, 1, 0, 0), sampleAt(a, 1, 0, 1)],
            [sampleAt(a, 1, 1, 0), sampleAt(a, 1, 1, 1)]]
    ] as const;

    it('reproduces every triquintic from exact samples of its 27 derivatives', () => {
        check(fc.tuple(fc.array(finite(-1, 1), { minLength: 216, maxLength: 216 }),
            scaled(0, 1), scaled(0, 1), scaled(0, 1)), ([a, x, y, z]) => {
            const h = new HermiteTriquintic(blocksOf(a));
            // The mixed high-order corner data of a unit-coefficient
            // triquintic reach 8000 per term and generateSingle weights the
            // value by up to 1000, so intermediates are of order 1e5 while the
            // reconstructed value is of order 10. The absolute tolerance is
            // that intermediate scale times a few hundred ulps.
            const orders: ReadonlyArray<readonly [number, number, number]> = [
                [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1], [2, 0, 0]];
            for (const m of orders) {
                expectClose(h.evaluate(m[0], m[1], m[2], x, y, z),
                    evalPoly(a, m, x, y, z), 1e-7, 1e-9);
            }
        }, 50);
    }, 30000);

    it('interpolates the prescribed data at all eight corners', () => {
        check(fc.array(finite(-2, 2), { minLength: 27 * 8, maxLength: 27 * 8 }), s => {
            const mk = (n: number) => sampleOf(s.slice(27 * n, 27 * n + 27));
            const blocks = [[[mk(0), mk(1)], [mk(2), mk(3)]],
                [[mk(4), mk(5)], [mk(6), mk(7)]]] as const;
            const h = new HermiteTriquintic(blocks);
            for (let b0 = 0; b0 <= 1; ++b0) {
                for (let b1 = 0; b1 <= 1; ++b1) {
                    for (let b2 = 0; b2 <= 1; ++b2) {
                        const b = blocks[b0][b1][b2];
                        for (let n = 0; n < ORDERS.length; ++n) {
                            const m = ORDERS[n];
                            // Each derivative order multiplies the roundoff by
                            // up to 120 (the largest basis derivative), so the
                            // sixth-order mixed term keeps about eight digits.
                            expectClose(h.evaluate(m[0], m[1], m[2], b0, b1, b2),
                                b[FIELDS[n]], 1e-6, 1e-7);
                        }
                    }
                }
            }
        }, 50);
    }, 30000);

    it('agrees with the tensor product of three HermiteQuintics', () => {
        check(fc.tuple(fc.array(finite(-2, 2), { minLength: 18, maxLength: 18 }),
            scaled(0, 1), scaled(0, 1), scaled(0, 1)), ([s, x, y, z]) => {
            const mk1 = (o: number) => new HermiteQuintic([
                new HermiteQuinticSample(s[o], s[o + 1], s[o + 2]),
                new HermiteQuinticSample(s[o + 3], s[o + 4], s[o + 5])]);
            const g = mk1(0), k = mk1(6), l = mk1(12);
            const mk = (bx: number, by: number, bz: number) =>
                sampleOf(ORDERS.map(m => g.evaluate(m[0], bx) * k.evaluate(m[1], by)
                    * l.evaluate(m[2], bz)));
            const h = new HermiteTriquintic([
                [[mk(0, 0, 0), mk(0, 0, 1)], [mk(0, 1, 0), mk(0, 1, 1)]],
                [[mk(1, 0, 0), mk(1, 0, 1)], [mk(1, 1, 0), mk(1, 1, 1)]]]);
            expectClose(h.evaluate(0, 0, 0, x, y, z),
                g.evaluate(0, x) * k.evaluate(0, y) * l.evaluate(0, z), 1e-8, 1e-9);
            expectClose(h.evaluate(1, 1, 1, x, y, z),
                g.evaluate(1, x) * k.evaluate(1, y) * l.evaluate(1, z), 1e-7, 1e-8);
        }, 100);
    }, 30000);

    it('exposes c publicly for manual coefficient assignment (upstream API)', () => {
        const h = new HermiteTriquintic();
        h.c[5][5][5] = 7;
        // P(5,t) = t^5, so H(x,y,z) = 7 x^5 y^5 z^5.
        expectClose(h.evaluate(0, 0, 0, 1, 1, 1), 7, 1e-15, 1e-15);
        expect(h.evaluate(0, 0, 0, 0, 1, 1)).toBe(0);
    });

    it('returns zero when any order exceeds the degree', () => {
        check(fc.tuple(fc.integer({ min: 6, max: 20 }), scaled(0, 1)), ([big, t]) => {
            const h = new HermiteTriquintic();
            h.c[2][3][4] = 1;
            expect(h.evaluate(big, 0, 0, t, t, t)).toBe(0);
            expect(h.evaluate(0, big, 0, t, t, t)).toBe(0);
            expect(h.evaluate(0, 0, big, t, t, t)).toBe(0);
        });
    });
});
