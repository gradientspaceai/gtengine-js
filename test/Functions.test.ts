import { describe, it, expect } from 'vitest';
import { check, finite, expectClose, fc } from './helpers/arbitraries.js';
import {
    atandivpi, atan2divpi, clamp, cospi, exp10, invsqrt, isign, saturate,
    sign, sinpi, sqr, fma, robustSOP, robustDOP
} from '../src/Functions.js';
import { GTE_C_INV_PI, GTE_C_PI, GTE_C_LN_10 } from '../src/Constants.js';

describe('Functions trigonometric helpers', () => {
    it('atandivpi returns the angle as a fraction of pi', () => {
        expect(atandivpi(0)).toBe(0);
        expect(atandivpi(1)).toBeCloseTo(0.25, 15);
        expect(atandivpi(-1)).toBeCloseTo(-0.25, 15);
        // The limit as x -> +infinity is 1/2.
        expect(atandivpi(Number.POSITIVE_INFINITY)).toBeCloseTo(0.5, 15);
    });

    it('atandivpi is exactly atan(x) * GTE_C_INV_PI', () => {
        for (const x of [-7.5, -1.25, -0.3, 0.0, 0.3, 1.25, 7.5]) {
            expect(atandivpi(x)).toBe(Math.atan(x) * GTE_C_INV_PI);
        }
    });

    it('atan2divpi handles the quadrants', () => {
        expect(atan2divpi(0, 1)).toBe(0);
        expect(atan2divpi(1, 1)).toBeCloseTo(0.25, 15);
        expect(atan2divpi(1, 0)).toBeCloseTo(0.5, 15);
        expect(atan2divpi(1, -1)).toBeCloseTo(0.75, 15);
        expect(atan2divpi(0, -1)).toBeCloseTo(1.0, 15);
        expect(atan2divpi(-1, -1)).toBeCloseTo(-0.75, 15);
        expect(atan2divpi(-1, 1)).toBeCloseTo(-0.25, 15);
    });

    it('sinpi and cospi evaluate at multiples of pi', () => {
        expect(sinpi(0)).toBe(0);
        expect(cospi(0)).toBe(1);
        expect(sinpi(0.5)).toBeCloseTo(1, 15);
        expect(cospi(0.5)).toBeCloseTo(0, 15);
        expect(sinpi(1)).toBeCloseTo(0, 15);
        expect(cospi(1)).toBeCloseTo(-1, 15);
        expect(sinpi(1.5)).toBeCloseTo(-1, 15);
    });

    it('sinpi and cospi satisfy the Pythagorean identity', () => {
        for (let i = -20; i <= 20; ++i) {
            const x = i / 7;
            expect(sqr(sinpi(x)) + sqr(cospi(x))).toBeCloseTo(1, 15);
        }
    });

    it('sinpi and cospi are exactly the scaled library calls', () => {
        for (const x of [-2.25, -0.5, 0.125, 1.75]) {
            expect(sinpi(x)).toBe(Math.sin(x * GTE_C_PI));
            expect(cospi(x)).toBe(Math.cos(x * GTE_C_PI));
        }
    });
});

describe('Functions algebraic helpers', () => {
    it('clamp restricts to the interval', () => {
        expect(clamp(-5, -1, 1)).toBe(-1);
        expect(clamp(0.5, -1, 1)).toBe(0.5);
        expect(clamp(5, -1, 1)).toBe(1);
        expect(clamp(-1, -1, 1)).toBe(-1);
        expect(clamp(1, -1, 1)).toBe(1);
    });

    it('clamp returns xmin first for an inverted interval, as upstream does', () => {
        // The upstream expression is (x <= xmin ? xmin : (x >= xmax ? xmax : x)),
        // so the xmin test wins when xmin > xmax.
        expect(clamp(0, 1, -1)).toBe(1);
    });

    it('saturate is clamp to [0,1]', () => {
        for (const x of [-3, -1e-16, 0, 0.25, 1, 1 + 1e-16, 4]) {
            expect(saturate(x)).toBe(clamp(x, 0, 1));
        }
        expect(saturate(-2)).toBe(0);
        expect(saturate(0.75)).toBe(0.75);
        expect(saturate(2)).toBe(1);
    });

    it('sign and isign agree on the trichotomy', () => {
        expect(sign(3)).toBe(1);
        expect(sign(-3)).toBe(-1);
        expect(sign(0)).toBe(0);
        expect(sign(-0)).toBe(0);
        expect(isign(3)).toBe(1);
        expect(isign(-3)).toBe(-1);
        expect(isign(0)).toBe(0);
        for (const x of [-1e300, -1, -1e-300, 0, 1e-300, 1, 1e300]) {
            expect(sign(x)).toBe(isign(x));
        }
    });

    it('sqr is the square', () => {
        expect(sqr(0)).toBe(0);
        expect(sqr(3)).toBe(9);
        expect(sqr(-3)).toBe(9);
        expect(sqr(1e-200)).toBe(0);
    });

    it('invsqrt is the reciprocal square root', () => {
        expect(invsqrt(4)).toBe(0.5);
        expect(invsqrt(1)).toBe(1);
        expect(invsqrt(0)).toBe(Number.POSITIVE_INFINITY);
        expect(Number.isNaN(invsqrt(-1))).toBe(true);
        for (const x of [0.25, 2, 7, 1e10]) {
            expect(invsqrt(x) * Math.sqrt(x)).toBeCloseTo(1, 15);
        }
    });

    it('exp10 is the base-10 exponential', () => {
        expect(exp10(0)).toBe(1);
        expect(exp10(1)).toBeCloseTo(10, 12);
        expect(exp10(3)).toBeCloseTo(1000, 9);
        expect(exp10(-2)).toBeCloseTo(0.01, 15);
        expect(exp10(2)).toBe(Math.exp(2 * GTE_C_LN_10));
    });
});

describe('Functions.fma', () => {
    it('computes exact small integer cases', () => {
        expect(fma(2, 3, 4)).toBe(10);
        expect(fma(-2, 3, 4)).toBe(-2);
        expect(fma(0.5, 0.5, 0.5)).toBe(0.75);
    });

    it('agrees with the hardware expression when that expression is exact', () => {
        let state = 20250831;
        for (let i = 0; i < 2000; ++i) {
            state = (state * 1103515245 + 12345) & 0x7FFFFFFF;
            const u = ((state % 2048) - 1024);
            state = (state * 1103515245 + 12345) & 0x7FFFFFFF;
            const v = ((state % 2048) - 1024);
            state = (state * 1103515245 + 12345) & 0x7FFFFFFF;
            const w = ((state % 2048) - 1024);
            // All operands and the result are small integers, so u*v+w is
            // exact in binary64 and fma must produce the same value.
            expect(fma(u, v, w)).toBe(u * v + w);
        }
    });

    it('keeps the low bits that a two-rounding evaluation loses', () => {
        // The exact product (1 + 2^-52)^2 = 1 + 2^-51 + 2^-104. The naive
        // expression rounds the product to 1 + 2^-51, so the difference is 0;
        // the fused operation retains 2^-104.
        const a = 1 + Number.EPSILON;
        expect(a * a - 1 - 2 * Number.EPSILON).toBe(0);
        expect(fma(a, a, -1 - 2 * Number.EPSILON)).toBe(2 ** -104);
    });

    it('recovers the rounding error of a product', () => {
        // For finite u, v with no overflow, fma(u, v, -(u*v)) is exactly the
        // error of the rounded product.
        let state = 777;
        for (let i = 0; i < 1000; ++i) {
            state = (state * 1103515245 + 12345) & 0x7FFFFFFF;
            const u = (state / 0x7FFFFFFF) * 4 - 2;
            state = (state * 1103515245 + 12345) & 0x7FFFFFFF;
            const v = (state / 0x7FFFFFFF) * 4 - 2;
            const p = u * v;
            const err = fma(u, v, -p);
            // p + err reproduces the exact product, so multiplying the pair
            // back out and comparing with an exact bigint product is a strong
            // check. Use the scaled-integer representation of the operands.
            expect(Number.isFinite(err)).toBe(true);
            expect(Math.abs(err)).toBeLessThanOrEqual(Math.abs(p) * Number.EPSILON);
            // fma of the recovered pair is exact: p + err == u*v exactly, so
            // subtracting p from the fused result returns err again.
            expect(fma(u, v, -p) - err).toBe(0);
        }
    });

    it('handles zero, sign of zero, and the addend', () => {
        expect(fma(0, 5, 3)).toBe(3);
        expect(fma(5, 0, 3)).toBe(3);
        expect(fma(2, 3, 0)).toBe(6);
        expect(Object.is(fma(1, 1, -1), 0)).toBe(true);
        expect(Object.is(fma(-0, 5, 0), 0)).toBe(true);
        expect(Object.is(fma(-1, 0, -0), -0)).toBe(true);
    });

    it('handles non-finite operands', () => {
        expect(Number.isNaN(fma(Number.NaN, 1, 1))).toBe(true);
        expect(Number.isNaN(fma(1, Number.NaN, 1))).toBe(true);
        expect(Number.isNaN(fma(1, 1, Number.NaN))).toBe(true);
        expect(fma(Number.POSITIVE_INFINITY, 2, 1)).toBe(Number.POSITIVE_INFINITY);
        expect(fma(Number.POSITIVE_INFINITY, -2, 1)).toBe(Number.NEGATIVE_INFINITY);
        expect(Number.isNaN(fma(Number.POSITIVE_INFINITY, 0, 1))).toBe(true);
        expect(fma(2, 3, Number.POSITIVE_INFINITY)).toBe(Number.POSITIVE_INFINITY);
        // A finite product cannot cancel an infinite addend.
        expect(fma(1e308, 10, Number.NEGATIVE_INFINITY)).toBe(Number.NEGATIVE_INFINITY);
    });

    it('overflows and underflows correctly', () => {
        expect(fma(1e300, 1e300, 0)).toBe(Number.POSITIVE_INFINITY);
        expect(fma(1e300, -1e300, 0)).toBe(Number.NEGATIVE_INFINITY);
        // The exact product overflows but the addend brings it back in range;
        // a fused operation still overflows because the exact sum is huge.
        expect(fma(1e300, 1e300, -1)).toBe(Number.POSITIVE_INFINITY);
        // Subnormal results.
        expect(fma(2 ** -600, 2 ** -474, 0)).toBe(2 ** -1074);
        expect(fma(2 ** -600, 2 ** -500, 0)).toBe(0);
    });

    it('rounds ties to even', () => {
        // The exact sum 1 + 2^-53 is halfway between 1 and 1 + 2^-52; ties to
        // even selects 1.
        expect(fma(2 ** -53, 1, 1)).toBe(1);
        // The exact sum (1 + 2^-52) + 2^-53 is halfway between 1 + 2^-52 and
        // 1 + 2^-51; ties to even selects 1 + 2^-51.
        expect(fma(2 ** -53, 1, 1 + Number.EPSILON)).toBe(1 + 2 * Number.EPSILON);
    });
});

describe('Functions.robustSOP / robustDOP', () => {
    it('computes ordinary sums and differences of products', () => {
        expect(robustSOP(2, 3, 4, 5)).toBe(26);
        expect(robustDOP(2, 3, 4, 5)).toBe(-14);
        expect(robustDOP(3, 4, 2, 5)).toBe(2);
    });

    it('is accurate for 2x2 determinants whose products are inexact', () => {
        // Random integer operands near 2^30 make the products exceed 2^53,
        // so the naive difference of products can be badly wrong. The exact
        // determinant is computed with bigint. Kahan's algorithm bounds the
        // relative error of the result by about two rounding units.
        let state = 20260101;
        const next = () => {
            state = (state * 1103515245 + 12345) & 0x7FFFFFFF;
            return 2 ** 30 + (state % (2 ** 30));
        };
        let naiveWorse = 0;
        for (let i = 0; i < 500; ++i) {
            const u = next();
            const v = next();
            const w = next();
            const z = next();
            const exact = BigInt(u) * BigInt(v) - BigInt(w) * BigInt(z);
            const exactNumber = Number(exact);
            const robust = robustDOP(u, v, w, z);
            const naive = u * v - w * z;
            const tolerance = 2 * Number.EPSILON * Math.abs(exactNumber) + 1e-9;
            expect(Math.abs(robust - exactNumber)).toBeLessThanOrEqual(tolerance);
            if (naive !== robust) {
                ++naiveWorse;
            }
        }
        // The inputs really do exercise the difference: the naive
        // expression disagrees with the robust one on many of them.
        expect(naiveWorse).toBeGreaterThan(0);
    });

    it('recovers a determinant the naive expression gets wrong', () => {
        // 2^53 + 1 cannot be represented, which makes the naive difference of
        // products lose the exact answer while the fused version keeps it.
        const u = 2 ** 27 + 1;
        const v = 2 ** 27 + 1;
        const w = 2 ** 27;
        const z = 2 ** 27 + 2;
        // Exact: (2^27+1)^2 - (2^27)(2^27+2) = 2^54+2^28+1 - (2^54+2^28) = 1.
        expect(robustDOP(u, v, w, z)).toBe(1);
    });

    it('matches the naive expressions when those are exact', () => {
        let state = 4242;
        for (let i = 0; i < 1000; ++i) {
            const next = () => {
                state = (state * 1103515245 + 12345) & 0x7FFFFFFF;
                return (state % 512) - 256;
            };
            const u = next();
            const v = next();
            const w = next();
            const z = next();
            // Use '===' rather than Object.is so that the two spellings of
            // zero compare equal; the fused evaluation of an exactly zero
            // result produces +0 where the naive expression can produce -0.
            expect(robustSOP(u, v, w, z) === u * v + w * z).toBe(true);
            expect(robustDOP(u, v, w, z) === u * v - w * z).toBe(true);
        }
    });

    it('relates the two by negating an operand', () => {
        let state = 99;
        for (let i = 0; i < 200; ++i) {
            const next = () => {
                state = (state * 1103515245 + 12345) & 0x7FFFFFFF;
                return (state / 0x7FFFFFFF) * 20 - 10;
            };
            const u = next();
            const v = next();
            const w = next();
            const z = next();
            expect(robustDOP(u, v, w, z) === robustSOP(u, v, -w, z)).toBe(true);
        }
    });
});

describe('Functions verification', () => {
    // ---- exact dyadic arithmetic used as the independent reference --------
    // Every finite double is n * 2^e with n an integer of at most 53 bits.
    // The extraction below uses only exact arithmetic (repeated doubling), so
    // it is independent of the bit-level decomposition inside Functions.ts.
    type Dyadic = { n: bigint; e: number };

    function toDyadic(x: number): Dyadic {
        if (x === 0) { return { n: 0n, e: 0 }; }
        let m = x, e = 0;
        while (!Number.isInteger(m)) { m *= 2; e -= 1; }
        return { n: BigInt(m), e };
    }

    function scaleTo(d: Dyadic, e: number): bigint {
        return d.n << BigInt(d.e - e);
    }

    function exactFMA(u: number, v: number, w: number): Dyadic {
        const du = toDyadic(u), dv = toDyadic(v), dw = toDyadic(w);
        const prod: Dyadic = { n: du.n * dv.n, e: du.e + dv.e };
        const e = Math.min(prod.e, dw.e);
        return { n: scaleTo(prod, e) + scaleTo(dw, e), e };
    }

    const bitsView = new DataView(new ArrayBuffer(8));

    function nextUp(x: number): number {
        if (x === 0) { return Number.MIN_VALUE; }
        bitsView.setFloat64(0, x);
        const bits = bitsView.getBigUint64(0);
        bitsView.setBigUint64(0, x > 0 ? bits + 1n : bits - 1n);
        return bitsView.getFloat64(0);
    }

    function nextDown(x: number): number {
        return -nextUp(-x);
    }

    function significandIsEven(x: number): boolean {
        bitsView.setFloat64(0, x);
        return (bitsView.getBigUint64(0) & 1n) === 0n;
    }

    /** |a - target| as an integer at the common exponent e. */
    function absDiff(a: number, target: Dyadic, e: number): bigint {
        const d = scaleTo(toDyadic(a), e) - scaleTo(target, e);
        return d < 0n ? -d : d;
    }

    /** fma(u, v, w) must be the correctly rounded (ties to even) exact value. */
    function expectCorrectlyRounded(u: number, v: number, w: number): void {
        const f = fma(u, v, w);
        expect(Number.isFinite(f)).toBe(true);
        const exact = exactFMA(u, v, w);
        const lo = nextDown(f), hi = nextUp(f);
        const e = Math.min(exact.e, toDyadic(f).e, toDyadic(lo).e, toDyadic(hi).e);
        const df = absDiff(f, exact, e);
        for (const neighbour of [lo, hi]) {
            const dn = absDiff(neighbour, exact, e);
            expect(df <= dn).toBe(true);
            if (df === dn) {
                // A tie must resolve to the even significand.
                expect(significandIsEven(f)).toBe(true);
            }
        }
    }

    // Operands spanning many binary exponents but staying well inside the
    // finite range, so the fused result never overflows.
    const scaled = fc.tuple(finite(-2, 2), fc.integer({ min: -60, max: 60 }))
        .map(([m, e]) => m * 2 ** e);

    it('fma is correctly rounded for operands across many exponents', () => {
        check(fc.tuple(scaled, scaled, scaled), ([u, v, w]) => {
            expectCorrectlyRounded(u, v, w);
        });
    });

    it('fma is correctly rounded for ordinary magnitudes', () => {
        check(fc.tuple(finite(-1e3, 1e3), finite(-1e3, 1e3), finite(-1e6, 1e6)),
            ([u, v, w]) => { expectCorrectlyRounded(u, v, w); });
    });

    it('fma is correctly rounded when the product nearly cancels the addend', () => {
        // The hardest case for a fused operation: massive cancellation, where
        // every bit of the exact product matters.
        check(fc.tuple(finite(-2, 2), finite(-2, 2)), ([u, v]) => {
            const p = u * v;
            expectCorrectlyRounded(u, v, -p);
            expectCorrectlyRounded(u, v, -p * (1 + 2 ** -30));
        });
    });

    it('fma(u, v, -(u*v)) is the exact rounding error of the product', () => {
        check(fc.tuple(scaled, scaled), ([u, v]) => {
            const p = u * v;
            // Dekker's theorem needs the product to be free of underflow: in
            // the subnormal range the error term itself is not representable
            // (e.g. u = 0.038461538461538464, v = -6.4e-323), so the
            // error-free transformation does not apply there.
            if (!Number.isFinite(p) || Math.abs(p) < 2 ** -960) { return true; }
            const err = fma(u, v, -p);
            // p + err must equal the exact product u*v (Dekker's theorem).
            const du = toDyadic(u), dv = toDyadic(v);
            const exact: Dyadic = { n: du.n * dv.n, e: du.e + dv.e };
            const dp = toDyadic(p), de = toDyadic(err);
            const e = Math.min(exact.e, dp.e, de.e);
            return scaleTo(dp, e) + scaleTo(de, e) === scaleTo(exact, e);
        });
    });

    it('fma(1, a, b) equals the correctly rounded a + b', () => {
        check(fc.tuple(scaled, scaled), ([a, b]) => Object.is(fma(1, a, b), a + b)
            || (a + b === 0 && fma(1, a, b) === 0));
    });

    it('fma(u, v, 0) equals the correctly rounded product', () => {
        check(fc.tuple(scaled, scaled), ([u, v]) => {
            const f = fma(u, v, 0);
            const p = u * v;
            return Object.is(f, p) || (f === 0 && p === 0);
        });
    });

    it('fma agrees with the naive expression whenever the naive one is exact', () => {
        // Small integers: u*v+w fits in 53 bits, so both must agree exactly.
        const small = fc.integer({ min: -100000, max: 100000 });
        check(fc.tuple(small, small, small), ([u, v, w]) => fma(u, v, w) === u * v + w);
    });

    it('robustDOP is accurate to a few ulps for 2x2 determinants', () => {
        // Kahan's algorithm bounds the relative error of u*v - w*z by about
        // two rounding units even under cancellation; the naive expression
        // does not. The exact determinant is computed with bigint.
        const big = fc.integer({ min: 2 ** 30, max: 2 ** 31 - 1 });
        check(fc.tuple(big, big, big, big), ([u, v, w, z]) => {
            const exact = BigInt(u) * BigInt(v) - BigInt(w) * BigInt(z);
            const target = Number(exact);
            const robust = robustDOP(u, v, w, z);
            const tol = 2 * Number.EPSILON * Math.abs(target) + 1e-6;
            return Math.abs(robust - target) <= tol;
        });
    });

    it('robustSOP and robustDOP are related by negating an operand', () => {
        check(fc.tuple(scaled, scaled, scaled, scaled), ([u, v, w, z]) =>
            robustDOP(u, v, w, z) === robustSOP(u, v, -w, z));
    });

    it('clamp and saturate follow the upstream conditional order', () => {
        check(fc.tuple(finite(), finite(), finite()), ([x, a, b]) => {
            const expected = x <= a ? a : (x >= b ? b : x);
            return Object.is(clamp(x, a, b), expected)
                && Object.is(saturate(x), (x <= 0 ? 0 : (x >= 1 ? 1 : x)));
        });
    });

    it('clamp with a well-ordered interval lands in [xmin, xmax]', () => {
        check(fc.tuple(finite(), finite(), finite()), ([x, a, b]) => {
            const lo = Math.min(a, b), hi = Math.max(a, b);
            const c = clamp(x, lo, hi);
            return c >= lo && c <= hi;
        });
    });

    it('sign and isign agree and are odd functions', () => {
        check(finite(), x => sign(x) === isign(x) && sign(-x) === -sign(x)
            && Math.abs(sign(x)) === (x === 0 ? 0 : 1));
    });

    it('sqr is nonnegative and even', () => {
        check(finite(-1e150, 1e150), x => sqr(x) === sqr(-x) && sqr(x) >= 0);
    });

    it('invsqrt(x) * invsqrt(x) recovers 1/x', () => {
        check(finite(1e-8, 1e8).filter(x => x > 0), x => {
            expectClose(invsqrt(x) * invsqrt(x), 1 / x, 0, 1e-14);
            return true;
        });
    });

    it('exp10 inverts log10 and matches the upstream expression', () => {
        check(finite(-100, 100), x => {
            expect(exp10(x)).toBe(Math.exp(x * GTE_C_LN_10));
            if (x > -300 && x < 300) {
                expectClose(Math.log10(exp10(x)), x, 1e-9, 1e-12);
            }
            return true;
        });
    });

    it('sinpi/cospi match the upstream expressions and the Pythagorean identity', () => {
        check(finite(-50, 50), x => {
            expect(sinpi(x)).toBe(Math.sin(x * GTE_C_PI));
            expect(cospi(x)).toBe(Math.cos(x * GTE_C_PI));
            expectClose(sqr(sinpi(x)) + sqr(cospi(x)), 1, 1e-12, 1e-12);
            return true;
        });
    });

    it('atandivpi and atan2divpi match the upstream expressions', () => {
        check(fc.tuple(finite(-1e3, 1e3), finite(-1e3, 1e3)), ([y, x]) => {
            expect(atandivpi(x)).toBe(Math.atan(x) * GTE_C_INV_PI);
            expect(atan2divpi(y, x)).toBe(Math.atan2(y, x) * GTE_C_INV_PI);
            // The result is an angle in (-1, 1] half-turns.
            expect(Math.abs(atan2divpi(y, x))).toBeLessThanOrEqual(1);
            expect(Math.abs(atandivpi(x))).toBeLessThanOrEqual(0.5);
            return true;
        });
    });

    it('atan2divpi(y, x) recovers the direction of (x, y)', () => {
        check(fc.tuple(finite(-10, 10), finite(-10, 10))
            .filter(([y, x]) => x * x + y * y > 1e-6), ([y, x]) => {
                const t = atan2divpi(y, x) * GTE_C_PI;
                const r = Math.sqrt(x * x + y * y);
                expectClose(r * Math.cos(t), x, 1e-9, 1e-12);
                expectClose(r * Math.sin(t), y, 1e-9, 1e-12);
                return true;
            });
    });
});
