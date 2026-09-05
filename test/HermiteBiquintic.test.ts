import { describe, it, expect } from 'vitest';
import { HermiteBiquintic, HermiteBiquinticSample } from '../src/HermiteBiquintic.js';
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

// ---------------------------------------------------------------------------
// Verification (V28): property-based cross-checks against HermiteBiquintic.h.
// ---------------------------------------------------------------------------
describe('HermiteBiquintic verification', () => {
    // d^k/dt^k t^i evaluated at t (zero once k exceeds i).
    const mono = (k: number, i: number, t: number): number => {
        if (i < k) { return 0; }
        let c = 1;
        for (let m = 0; m < k; ++m) { c *= i - m; }
        return c * Math.pow(t, i - k);
    };
    // p(x,y) = sum_{i,j<=5} a[6i+j] x^i y^j, differentiated mx times in x and
    // my times in y.
    const evalPoly = (a: number[], mx: number, my: number,
        x: number, y: number): number => {
        let sum = 0;
        for (let i = 0; i <= 5; ++i) {
            const mi = mono(mx, i, x);
            if (mi === 0) { continue; }
            for (let j = 0; j <= 5; ++j) {
                sum += a[6 * i + j] * mi * mono(my, j, y);
            }
        }
        return sum;
    };
    const sampleAt = (a: number[], bx: number, by: number) =>
        new HermiteBiquinticSample(
            evalPoly(a, 0, 0, bx, by), evalPoly(a, 1, 0, bx, by),
            evalPoly(a, 0, 1, bx, by), evalPoly(a, 2, 0, bx, by),
            evalPoly(a, 1, 1, bx, by), evalPoly(a, 0, 2, bx, by),
            evalPoly(a, 2, 1, bx, by), evalPoly(a, 1, 2, bx, by),
            evalPoly(a, 2, 2, bx, by));

    it('reproduces every biquintic from exact samples of its nine derivatives', () => {
        check(fc.tuple(fc.array(finite(-1, 1), { minLength: 36, maxLength: 36 }),
            scaled(0, 1), scaled(0, 1)), ([a, x, y]) => {
            const h = new HermiteBiquintic([
                [sampleAt(a, 0, 0), sampleAt(a, 0, 1)],
                [sampleAt(a, 1, 0), sampleAt(a, 1, 1)]]);
            // generateSingle weights the data by up to 100 and the corner data
            // themselves reach ~400 for unit coefficients, so the worst-case
            // cancellation costs about four digits of the 16 available.
            for (const [mx, my] of [[0, 0], [1, 0], [0, 1], [2, 0], [1, 1], [0, 2]]) {
                expectClose(h.evaluate(mx, my, x, y), evalPoly(a, mx, my, x, y),
                    1e-8, 1e-9);
            }
        });
    });

    it('interpolates the prescribed data at all four corners', () => {
        check(fc.array(finite(-4, 4), { minLength: 36, maxLength: 36 }), s => {
            const mk = (k: number) => new HermiteBiquinticSample(
                s[9 * k], s[9 * k + 1], s[9 * k + 2], s[9 * k + 3], s[9 * k + 4],
                s[9 * k + 5], s[9 * k + 6], s[9 * k + 7], s[9 * k + 8]);
            const blocks = [[mk(0), mk(1)], [mk(2), mk(3)]] as const;
            const h = new HermiteBiquintic(blocks);
            for (let b0 = 0; b0 <= 1; ++b0) {
                for (let b1 = 0; b1 <= 1; ++b1) {
                    const b = blocks[b0][b1];
                    expectClose(h.evaluate(0, 0, b0, b1), b.f, 1e-11, 1e-11);
                    expectClose(h.evaluate(1, 0, b0, b1), b.fx, 1e-9, 1e-9);
                    expectClose(h.evaluate(0, 1, b0, b1), b.fy, 1e-9, 1e-9);
                    expectClose(h.evaluate(2, 0, b0, b1), b.fxx, 1e-8, 1e-8);
                    expectClose(h.evaluate(1, 1, b0, b1), b.fxy, 1e-8, 1e-8);
                    expectClose(h.evaluate(0, 2, b0, b1), b.fyy, 1e-8, 1e-8);
                    expectClose(h.evaluate(2, 2, b0, b1), b.fxxyy, 1e-6, 1e-7);
                }
            }
        });
    });

    it('agrees with the tensor product of two HermiteQuintics on separable data', () => {
        check(fc.tuple(fc.array(finite(-3, 3), { minLength: 6, maxLength: 6 }),
            scaled(0, 1), scaled(0, 1)), ([s, x, y]) => {
            const g = new HermiteQuintic([new HermiteQuinticSample(s[0], s[1], s[2]),
                new HermiteQuinticSample(s[3], s[4], s[5])]);
            const k = new HermiteQuintic([new HermiteQuinticSample(s[5], s[3], s[1]),
                new HermiteQuinticSample(s[0], s[2], s[4])]);
            const mk = (bx: number, by: number) => new HermiteBiquinticSample(
                g.evaluate(0, bx) * k.evaluate(0, by),
                g.evaluate(1, bx) * k.evaluate(0, by),
                g.evaluate(0, bx) * k.evaluate(1, by),
                g.evaluate(2, bx) * k.evaluate(0, by),
                g.evaluate(1, bx) * k.evaluate(1, by),
                g.evaluate(0, bx) * k.evaluate(2, by),
                g.evaluate(2, bx) * k.evaluate(1, by),
                g.evaluate(1, bx) * k.evaluate(2, by),
                g.evaluate(2, bx) * k.evaluate(2, by));
            const b = new HermiteBiquintic([[mk(0, 0), mk(0, 1)], [mk(1, 0), mk(1, 1)]]);
            expectClose(b.evaluate(0, 0, x, y),
                g.evaluate(0, x) * k.evaluate(0, y), 1e-9, 1e-10);
            expectClose(b.evaluate(2, 2, x, y),
                g.evaluate(2, x) * k.evaluate(2, y), 1e-7, 1e-8);
        });
    });

    it('exposes c publicly for manual coefficient assignment (upstream API)', () => {
        const h = new HermiteBiquintic();
        h.c[5][0] = 2;
        // P(5,x) P(0,y) = x^5 (1-y)^5, so H(1,0) = 2.
        expectClose(h.evaluate(0, 0, 1, 0), 2, 1e-15, 1e-15);
        expect(h.evaluate(0, 0, 0, 1)).toBe(0);
    });

    it('returns zero when either order exceeds the degree', () => {
        check(fc.tuple(fc.integer({ min: 6, max: 20 }), fc.integer({ min: 0, max: 5 }),
            scaled(0, 1), scaled(0, 1)), ([big, small, x, y]) => {
            const h = new HermiteBiquintic();
            h.c[3][2] = 1;
            expect(h.evaluate(big, small, x, y)).toBe(0);
            expect(h.evaluate(small, big, x, y)).toBe(0);
        });
    });
});
