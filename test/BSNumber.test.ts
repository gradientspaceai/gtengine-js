import { describe, it, expect } from 'vitest';
import { check, fc } from './helpers/arbitraries.js';
import { BSNumber, BSNumberRoundingMode, convertBSNumber } from '../src/BSNumber.js';
import { isArbitraryPrecision, hasDivisionOperator } from '../src/TypeTraits.js';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The exact value of a BSNumber as a rational with a power-of-two
// denominator: value = num * 2^exp. Used for independent cross-checks.
function exact(x: BSNumber): { num: bigint, exp: number } {
    const magnitude = x.getUInteger();
    const num = (x.getSign() < 0 ? -magnitude : magnitude);
    return { num, exp: x.getBiasedExponent() };
}

// Exact comparison of two BSNumber values via a common exponent.
function exactCompare(a: BSNumber, b: BSNumber): number {
    const ea = exact(a), eb = exact(b);
    const minExp = Math.min(ea.exp, eb.exp);
    const na = ea.num << BigInt(ea.exp - minExp);
    const nb = eb.num << BigInt(eb.exp - minExp);
    return na < nb ? -1 : (na > nb ? 1 : 0);
}

// The exact double value 'value' as a (numerator, power-of-two) pair, taken
// from the IEEE-754 bits rather than from BSNumber, for an independent check.
function doubleToExact(value: number): { num: bigint, exp: number } {
    const buffer = new DataView(new ArrayBuffer(8));
    buffer.setFloat64(0, value);
    const bits = buffer.getBigUint64(0);
    const s = (bits >> 63n) === 1n ? -1n : 1n;
    const e = Number((bits >> 52n) & 0x7FFn);
    const t = bits & 0xFFFFFFFFFFFFFn;
    if (e === 0) {
        return { num: s * t, exp: -1074 };
    }
    return { num: s * (t | 0x10000000000000n), exp: e - 1075 };
}

describe('BSNumber construction and representation', () => {
    it('default constructs zero', () => {
        const zero = new BSNumber();
        expect(zero.getSign()).toBe(0);
        expect(zero.getBiasedExponent()).toBe(0);
        expect(zero.getUInteger()).toBe(0n);
        expect(zero.isValid()).toBe(true);
        expect(zero.toNumber()).toBe(0);
    });

    it('normalizes the unsigned integer to an odd number', () => {
        // 1.5 = 3 * 2^{-1}
        const x = BSNumber.fromNumber(1.5);
        expect(x.getSign()).toBe(1);
        expect(x.getUInteger()).toBe(3n);
        expect(x.getBiasedExponent()).toBe(-1);
        expect(x.getExponent()).toBe(0);

        // -12 = -(3 * 2^2)
        const y = BSNumber.fromNumber(-12);
        expect(y.getSign()).toBe(-1);
        expect(y.getUInteger()).toBe(3n);
        expect(y.getBiasedExponent()).toBe(2);
        expect(y.getExponent()).toBe(3);

        // 1 = 1 * 2^0
        const one = BSNumber.fromNumber(1);
        expect(one.getUInteger()).toBe(1n);
        expect(one.getBiasedExponent()).toBe(0);
    });

    it('converts the minimum subnormal double exactly', () => {
        const x = BSNumber.fromNumber(Number.MIN_VALUE);
        expect(x.getSign()).toBe(1);
        expect(x.getUInteger()).toBe(1n);
        expect(x.getBiasedExponent()).toBe(-1074);
        expect(x.toNumber()).toBe(Number.MIN_VALUE);
    });

    it('converts a general subnormal double exactly', () => {
        // 5 * 2^{-1074}
        const value = 5 * Number.MIN_VALUE;
        const x = BSNumber.fromNumber(value);
        expect(x.getUInteger()).toBe(5n);
        expect(x.getBiasedExponent()).toBe(-1074);
        expect(x.toNumber()).toBe(value);
    });

    it('converts the maximum normal double exactly', () => {
        const x = BSNumber.fromNumber(Number.MAX_VALUE);
        const expected = doubleToExact(Number.MAX_VALUE);
        expect(exact(x).num << BigInt(exact(x).exp - expected.exp)).toBe(expected.num);
        expect(x.toNumber()).toBe(Number.MAX_VALUE);
    });

    it('has a graceful exit for infinities and NaN', () => {
        // Infinity has no representation; upstream returns
        // (-1)^s * 2^{1+EXPONENT_BIAS}, which converts back to infinity.
        const pos = BSNumber.fromNumber(Number.POSITIVE_INFINITY);
        expect(pos.getSign()).toBe(1);
        expect(pos.getUInteger()).toBe(1n);
        expect(pos.getBiasedExponent()).toBe(1024);
        expect(pos.toNumber()).toBe(Number.POSITIVE_INFINITY);

        const neg = BSNumber.fromNumber(Number.NEGATIVE_INFINITY);
        expect(neg.getSign()).toBe(-1);
        expect(neg.toNumber()).toBe(Number.NEGATIVE_INFINITY);

        // NaN has no representation; upstream returns zero.
        const nan = BSNumber.fromNumber(Number.NaN);
        expect(nan.getSign()).toBe(0);
        expect(nan.toNumber()).toBe(0);
    });

    it('treats negative zero as zero', () => {
        const x = BSNumber.fromNumber(-0);
        expect(x.getSign()).toBe(0);
        expect(x.isValid()).toBe(true);
    });

    it('converts big integers exactly', () => {
        const value = 123456789012345678901234567890n;
        const x = BSNumber.fromBigInt(value);
        expect(x.getSign()).toBe(1);
        expect(x.isValid()).toBe(true);
        expect(exact(x).num << BigInt(exact(x).exp)).toBe(value);

        const negative = BSNumber.fromBigInt(-value);
        expect(negative.getSign()).toBe(-1);
        expect(negative.getUInteger()).toBe(x.getUInteger());
        expect(negative.getBiasedExponent()).toBe(x.getBiasedExponent());

        expect(BSNumber.fromBigInt(0n).getSign()).toBe(0);
    });

    it('converts decimal strings exactly', () => {
        const digits = '123456789012345678901234567890';
        const x = BSNumber.fromString(digits);
        const y = BSNumber.fromBigInt(BigInt(digits));
        expect(exactCompare(x, y)).toBe(0);

        expect(exactCompare(BSNumber.fromString('+7'), BSNumber.fromNumber(7))).toBe(0);
        expect(exactCompare(BSNumber.fromString('-7'), BSNumber.fromNumber(-7))).toBe(0);
        expect(BSNumber.fromString('0').getSign()).toBe(0);
        // The sign is not applied to zero.
        expect(BSNumber.fromString('-0').getSign()).toBe(0);

        expect(() => BSNumber.fromString('')).toThrow();
        expect(() => BSNumber.fromString('-')).toThrow();
        expect(() => BSNumber.fromString('012')).toThrow();
        expect(() => BSNumber.fromString('1x2')).toThrow();
    });

    it('implements the arbitrary-precision marker interface', () => {
        const x = BSNumber.fromNumber(1);
        expect(isArbitraryPrecision(x)).toBe(true);
        expect(hasDivisionOperator(x)).toBe(false);
        expect(hasDivisionOperator(1)).toBe(true);
    });
});

describe('BSNumber double round trips', () => {
    it('round trips a set of known values', () => {
        const values = [
            0, 1, -1, 0.5, -0.5, 2, 1024, 1 / 3, -1 / 7, Math.PI, Math.E,
            1e-300, 1e300, -1e300, Number.MIN_VALUE, -Number.MIN_VALUE,
            Number.MAX_VALUE, -Number.MAX_VALUE, Number.EPSILON,
            2 ** 52 + 1, 2 ** 53, -(2 ** 53)
        ];
        for (const value of values) {
            expect(BSNumber.fromNumber(value).toNumber()).toBe(value);
        }
    });

    it('round trips random doubles over a wide exponent range', () => {
        const rand = makeRandom(12345);
        for (let i = 0; i < 500; ++i) {
            const mantissa = 2 * rand() - 1;
            const exponent = Math.floor(600 * rand()) - 300;
            const value = mantissa * Math.pow(2, exponent);
            expect(BSNumber.fromNumber(value).toNumber()).toBe(value);
        }
    });

    it('matches the IEEE bit pattern of the input', () => {
        const rand = makeRandom(999);
        for (let i = 0; i < 200; ++i) {
            const value = (2 * rand() - 1) * Math.pow(2, Math.floor(60 * rand()) - 30);
            const fromBits = doubleToExact(value);
            const fromBSN = exact(BSNumber.fromNumber(value));
            // Both represent the same rational; compare with a common
            // exponent.
            const minExp = Math.min(fromBits.exp, fromBSN.exp);
            expect(fromBSN.num << BigInt(fromBSN.exp - minExp))
                .toBe(fromBits.num << BigInt(fromBits.exp - minExp));
        }
    });

    it('rounds to nearest with ties to even when the value needs more than 53 bits', () => {
        const twoTo53 = 2n ** 53n;
        // 2^53 + 1 is exactly halfway between 2^53 and 2^53 + 2; the even
        // significand wins, so the result is 2^53.
        expect(BSNumber.fromBigInt(twoTo53 + 1n).toNumber()).toBe(2 ** 53);
        // 2^53 + 3 is exactly halfway between 2^53 + 2 and 2^53 + 4; the even
        // significand wins, so the result is 2^53 + 4.
        expect(BSNumber.fromBigInt(twoTo53 + 3n).toNumber()).toBe(2 ** 53 + 4);
        // 2^53 + 5 is closer to 2^53 + 4 than to 2^53 + 6.
        expect(BSNumber.fromBigInt(twoTo53 + 5n).toNumber()).toBe(2 ** 53 + 4);
        // Not a tie: the remainder is larger than 1/2.
        expect(BSNumber.fromBigInt(4n * twoTo53 + 5n).toNumber()).toBe(4 * 2 ** 53 + 8);
    });

    it('overflows to infinity and underflows to zero', () => {
        // 2^1024 exceeds the maximum normal.
        const big = BSNumber.ldexp(BSNumber.fromNumber(1), 1024);
        expect(big.toNumber()).toBe(Number.POSITIVE_INFINITY);
        expect(big.negated().toNumber()).toBe(Number.NEGATIVE_INFINITY);

        // 2^{-1075} is exactly halfway between 0 and the minimum subnormal,
        // so ties-to-even rounds it to zero.
        const tiny = BSNumber.ldexp(BSNumber.fromNumber(1), -1075);
        expect(tiny.toNumber()).toBe(0);

        // 1.5 * 2^{-1075} is closer to the minimum subnormal than to zero.
        const tiny2 = BSNumber.ldexp(BSNumber.fromNumber(3), -1076);
        expect(tiny2.toNumber()).toBe(Number.MIN_VALUE);
    });

    it('converts to binary32 the same way Math.fround does', () => {
        const rand = makeRandom(4242);
        const values = [0, 1, -1, 0.1, Math.PI, 1e-40, 1e-45, 3.4e38, -3.4e38];
        for (let i = 0; i < 200; ++i) {
            values.push((2 * rand() - 1) * Math.pow(2, Math.floor(60 * rand()) - 30));
        }
        for (const value of values) {
            expect(BSNumber.fromNumber(value).toFloat32()).toBe(Math.fround(value));
        }
    });

    it('converts from binary32 inputs', () => {
        const x = BSNumber.fromFloat32(0.1);
        expect(x.toNumber()).toBe(Math.fround(0.1));
        expect(x.isValid()).toBe(true);
    });
});

describe('BSNumber arithmetic', () => {
    it('adds and subtracts exactly across huge exponent gaps', () => {
        const a = BSNumber.fromNumber(1e300);
        const b = BSNumber.fromNumber(1e-300);
        const sum = a.add(b);
        // The double sum loses b entirely; the exact sum does not.
        expect(sum.toNumber()).toBe(1e300);
        expect(exactCompare(sum.sub(b), a)).toBe(0);
        expect(exactCompare(sum.sub(a), b)).toBe(0);
        expect(sum.isValid()).toBe(true);
    });

    it('satisfies (a + b) - b == a exactly for random inputs', () => {
        const rand = makeRandom(7);
        for (let i = 0; i < 300; ++i) {
            const a = BSNumber.fromNumber(
                (2 * rand() - 1) * Math.pow(2, Math.floor(200 * rand()) - 100));
            const b = BSNumber.fromNumber(
                (2 * rand() - 1) * Math.pow(2, Math.floor(200 * rand()) - 100));
            expect(exactCompare(a.add(b).sub(b), a)).toBe(0);
            expect(exactCompare(a.sub(b).add(b), a)).toBe(0);
            expect(a.add(b).isValid()).toBe(true);
            expect(a.sub(b).isValid()).toBe(true);
        }
    });

    it('cancels exactly', () => {
        const a = BSNumber.fromNumber(Math.PI);
        expect(a.sub(a).getSign()).toBe(0);
        expect(a.add(a.negated()).getSign()).toBe(0);
        expect(a.sub(a).isValid()).toBe(true);
    });

    it('multiplies exactly and cross-checks against bigint arithmetic', () => {
        const rand = makeRandom(31337);
        for (let i = 0; i < 300; ++i) {
            const av = (2 * rand() - 1) * Math.pow(2, Math.floor(100 * rand()) - 50);
            const bv = (2 * rand() - 1) * Math.pow(2, Math.floor(100 * rand()) - 50);
            const a = BSNumber.fromNumber(av);
            const b = BSNumber.fromNumber(bv);
            const product = a.mul(b);
            const ea = exact(a), eb = exact(b), ep = exact(product);
            expect(ep.num).toBe(ea.num * eb.num);
            expect(ep.exp).toBe(ea.exp + eb.exp);
            expect(product.isValid()).toBe(true);
        }
    });

    it('has 0 as the additive identity and annihilator of multiplication', () => {
        const zero = new BSNumber();
        const a = BSNumber.fromNumber(-3.25);
        expect(exactCompare(a.add(zero), a)).toBe(0);
        expect(exactCompare(zero.add(a), a)).toBe(0);
        expect(exactCompare(zero.sub(a), a.negated())).toBe(0);
        expect(exactCompare(a.sub(zero), a)).toBe(0);
        expect(a.mul(zero).getSign()).toBe(0);
        expect(zero.mul(a).getSign()).toBe(0);
    });

    it('computes exact sums and differences of products', () => {
        // A determinant that vanishes exactly but not in double arithmetic.
        const u = BSNumber.fromNumber(1 + 2 ** -52);
        const v = BSNumber.fromNumber(1 - 2 ** -52);
        const w = BSNumber.fromNumber(1);
        const z = BSNumber.fromNumber(1);
        const dop = BSNumber.robustDOP(u, v, w, z);
        // u*v - 1 = -2^{-104} exactly, whereas the double computation gives 0.
        expect(u.toNumber() * v.toNumber() - 1).toBe(0);
        expect(dop.getSign()).toBe(-1);
        expect(dop.toNumber()).toBe(-(2 ** -104));

        const sop = BSNumber.robustSOP(u, v, w, z);
        expect(exactCompare(sop, u.mul(v).add(w.mul(z)))).toBe(0);
        expect(exactCompare(BSNumber.fma(u, v, w), u.mul(v).add(w))).toBe(0);
    });

    it('negates in place and by copy', () => {
        const a = BSNumber.fromNumber(2.5);
        const b = a.negated();
        expect(a.getSign()).toBe(1);
        expect(b.getSign()).toBe(-1);
        b.negate();
        expect(b.getSign()).toBe(1);
        expect(exactCompare(a, b)).toBe(0);
    });

    it('clones without aliasing', () => {
        const a = BSNumber.fromNumber(3);
        const b = a.clone();
        b.negate();
        expect(a.getSign()).toBe(1);
        expect(b.getSign()).toBe(-1);
    });
});

describe('BSNumber comparisons', () => {
    it('orders known values', () => {
        const neg = BSNumber.fromNumber(-2);
        const zero = new BSNumber();
        const small = BSNumber.fromNumber(0.5);
        const big = BSNumber.fromNumber(1e100);

        expect(neg.lessThan(zero)).toBe(true);
        expect(zero.lessThan(small)).toBe(true);
        expect(small.lessThan(big)).toBe(true);
        expect(big.lessThan(neg)).toBe(false);
        expect(zero.lessThan(neg)).toBe(false);
        expect(neg.greaterThan(zero)).toBe(false);
        expect(zero.greaterThanOrEqual(zero)).toBe(true);
        expect(zero.lessThanOrEqual(zero)).toBe(true);
        expect(zero.equals(new BSNumber())).toBe(true);
        expect(zero.notEquals(small)).toBe(true);
    });

    it('agrees with double comparisons on random values', () => {
        const rand = makeRandom(24680);
        for (let i = 0; i < 500; ++i) {
            const av = (2 * rand() - 1) * Math.pow(2, Math.floor(80 * rand()) - 40);
            const bv = (2 * rand() - 1) * Math.pow(2, Math.floor(80 * rand()) - 40);
            const a = BSNumber.fromNumber(av);
            const b = BSNumber.fromNumber(bv);
            expect(a.lessThan(b)).toBe(av < bv);
            expect(a.greaterThan(b)).toBe(av > bv);
            expect(a.equals(b)).toBe(av === bv);
            expect(a.lessThanOrEqual(b)).toBe(av <= bv);
            expect(a.greaterThanOrEqual(b)).toBe(av >= bv);
        }
    });

    it('distinguishes values that are equal as doubles', () => {
        const a = BSNumber.fromBigInt(2n ** 53n);
        const b = BSNumber.fromBigInt(2n ** 53n + 1n);
        expect(a.toNumber()).toBe(b.toNumber());
        expect(a.equals(b)).toBe(false);
        expect(a.lessThan(b)).toBe(true);
    });

    it('compares numbers with the same exponent but different bit counts', () => {
        // 5 = 101b * 2^0 and 6 = 11b * 2^1 have the same exponent 2 but the
        // stored magnitudes compare the other way round (5 > 3). The unsigned
        // integers must be compared left-aligned, as fractions 1.u and 1.v.
        const a = BSNumber.fromNumber(5);
        const b = BSNumber.fromNumber(6);
        expect(a.getExponent()).toBe(b.getExponent());
        expect(a.getUInteger()).toBe(5n);
        expect(b.getUInteger()).toBe(3n);
        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(false);
        // The same case reached through subtraction, which uses the
        // comparison to order the operands.
        expect(a.sub(b).toNumber()).toBe(-1);
        expect(b.sub(a).toNumber()).toBe(1);
        expect(a.sub(b).isValid()).toBe(true);

        // 3 * 2^0 = 3 and 7 * 2^{-1} = 3.5 also share the exponent 1.
        const c = BSNumber.fromNumber(3);
        const d = BSNumber.fromNumber(3.5);
        expect(c.getExponent()).toBe(d.getExponent());
        expect(c.lessThan(d)).toBe(true);
        expect(d.lessThan(c)).toBe(false);
    });
});

describe('BSNumber math helpers', () => {
    it('frexp and ldexp are exact inverses', () => {
        const rand = makeRandom(555);
        for (let i = 0; i < 100; ++i) {
            const value = (2 * rand() - 1) * Math.pow(2, Math.floor(100 * rand()) - 50);
            const x = BSNumber.fromNumber(value);
            const { result, exponent } = BSNumber.frexp(x);
            // result is in [1/2, 1) in magnitude.
            const magnitude = Math.abs(result.toNumber());
            expect(magnitude).toBeGreaterThanOrEqual(0.5);
            expect(magnitude).toBeLessThan(1);
            expect(exactCompare(BSNumber.ldexp(result, exponent), x)).toBe(0);
        }

        const zero = BSNumber.frexp(new BSNumber());
        expect(zero.exponent).toBe(0);
        expect(zero.result.getSign()).toBe(0);
    });

    it('fabs is exact', () => {
        const a = BSNumber.fromNumber(-1.25);
        expect(BSNumber.fabs(a).toNumber()).toBe(1.25);
        expect(BSNumber.fabs(a.negated()).toNumber()).toBe(1.25);
        expect(BSNumber.fabs(new BSNumber()).getSign()).toBe(0);
    });

    it('routes the transcendental functions through double', () => {
        const half = BSNumber.fromNumber(0.5);
        const two = BSNumber.fromNumber(2);
        expect(BSNumber.sqrt(two).toNumber()).toBe(Math.SQRT2);
        expect(BSNumber.sin(half).toNumber()).toBe(Math.sin(0.5));
        expect(BSNumber.cos(half).toNumber()).toBe(Math.cos(0.5));
        expect(BSNumber.atan2(half, two).toNumber()).toBe(Math.atan2(0.5, 2));
        expect(BSNumber.floor(BSNumber.fromNumber(1.75)).toNumber()).toBe(1);
        expect(BSNumber.ceil(BSNumber.fromNumber(1.25)).toNumber()).toBe(2);
        expect(BSNumber.exp2(BSNumber.fromNumber(10)).toNumber()).toBe(1024);
        expect(BSNumber.log2(BSNumber.fromNumber(1024)).toNumber()).toBe(10);
        expect(BSNumber.sqr(BSNumber.fromNumber(-3)).toNumber()).toBe(9);
        expect(BSNumber.isign(BSNumber.fromNumber(-3))).toBe(-1);
        expect(BSNumber.sign(BSNumber.fromNumber(3)).toNumber()).toBe(1);
        expect(BSNumber.saturate(BSNumber.fromNumber(2)).toNumber()).toBe(1);
        expect(BSNumber.clamp(BSNumber.fromNumber(5), BSNumber.fromNumber(0),
            BSNumber.fromNumber(2)).toNumber()).toBe(2);
        expect(BSNumber.fmod(BSNumber.fromNumber(7), BSNumber.fromNumber(3)).toNumber()).toBe(1);
    });

    it('remainder rounds the quotient to the nearest even integer', () => {
        // 5/2 = 2.5 is a tie; the even quotient 2 wins, so 5 - 2*2 = 1.
        expect(BSNumber.remainder(BSNumber.fromNumber(5), BSNumber.fromNumber(2))
            .toNumber()).toBe(1);
        // 7/2 = 3.5 is a tie; the even quotient 4 wins, so 7 - 4*2 = -1.
        expect(BSNumber.remainder(BSNumber.fromNumber(7), BSNumber.fromNumber(2))
            .toNumber()).toBe(-1);
        expect(BSNumber.remainder(BSNumber.fromNumber(7), BSNumber.fromNumber(3))
            .toNumber()).toBe(1);
    });
});

describe('convertBSNumber', () => {
    it('returns the input when the precision is already satisfied', () => {
        const x = BSNumber.fromNumber(3);  // 3 = 11b, two bits
        for (const mode of [BSNumberRoundingMode.FE_TONEAREST,
            BSNumberRoundingMode.FE_UPWARD, BSNumberRoundingMode.FE_DOWNWARD,
            BSNumberRoundingMode.FE_TOWARDZERO]) {
            expect(exactCompare(convertBSNumber(x, 2, mode), x)).toBe(0);
            expect(exactCompare(convertBSNumber(x, 10, mode), x)).toBe(0);
        }
        expect(convertBSNumber(new BSNumber(), 5,
            BSNumberRoundingMode.FE_TONEAREST).getSign()).toBe(0);
    });

    it('rounds to nearest with ties to even', () => {
        // 11b = 3: rounding to one bit is a tie (remainder exactly 1/2) and
        // the retained bit is 1 (odd), so round up to 100b = 4.
        const three = BSNumber.fromNumber(3);
        expect(convertBSNumber(three, 1, BSNumberRoundingMode.FE_TONEAREST)
            .toNumber()).toBe(4);
        // 101b = 5: rounding to two bits is a tie and the retained low bit is
        // 0 (even), so round down to 100b = 4.
        const five = BSNumber.fromNumber(5);
        expect(convertBSNumber(five, 2, BSNumberRoundingMode.FE_TONEAREST)
            .toNumber()).toBe(4);
        // 111b = 7: rounding to two bits has a remainder larger than 1/2, so
        // round up to 1000b = 8.
        const seven = BSNumber.fromNumber(7);
        expect(convertBSNumber(seven, 2, BSNumberRoundingMode.FE_TONEAREST)
            .toNumber()).toBe(8);
        // 1101b = 13: rounding to two bits truncates to 1100b = 12 (the
        // remainder 01b is smaller than 1/2).
        const thirteen = BSNumber.fromNumber(13);
        expect(convertBSNumber(thirteen, 2, BSNumberRoundingMode.FE_TONEAREST)
            .toNumber()).toBe(12);
    });

    it('rounds toward zero, up and down', () => {
        const seven = BSNumber.fromNumber(7);
        const negSeven = BSNumber.fromNumber(-7);
        expect(convertBSNumber(seven, 2, BSNumberRoundingMode.FE_TOWARDZERO)
            .toNumber()).toBe(6);
        expect(convertBSNumber(negSeven, 2, BSNumberRoundingMode.FE_TOWARDZERO)
            .toNumber()).toBe(-6);
        expect(convertBSNumber(seven, 2, BSNumberRoundingMode.FE_UPWARD)
            .toNumber()).toBe(8);
        expect(convertBSNumber(negSeven, 2, BSNumberRoundingMode.FE_UPWARD)
            .toNumber()).toBe(-6);
        expect(convertBSNumber(seven, 2, BSNumberRoundingMode.FE_DOWNWARD)
            .toNumber()).toBe(6);
        expect(convertBSNumber(negSeven, 2, BSNumberRoundingMode.FE_DOWNWARD)
            .toNumber()).toBe(-8);
    });

    it('keeps the odd-integer invariant and the requested precision', () => {
        const rand = makeRandom(864);
        for (let i = 0; i < 200; ++i) {
            const value = (2 * rand() - 1) * Math.pow(2, Math.floor(60 * rand()) - 30);
            const x = BSNumber.fromNumber(value);
            if (x.getSign() === 0) {
                continue;
            }
            const precision = 1 + Math.floor(40 * rand());
            const y = convertBSNumber(x, precision, BSNumberRoundingMode.FE_TONEAREST);
            expect(y.isValid()).toBe(true);
            expect(y.getNumBits()).toBeLessThanOrEqual(Math.max(precision, x.getNumBits()));
            // The rounded value is within one unit in the last place.
            const ulp = Math.pow(2, x.getExponent() - precision + 1);
            expect(Math.abs(y.toNumber() - value)).toBeLessThanOrEqual(ulp);
        }
    });

    it('rounds up all-ones bit patterns into the next power of two', () => {
        // 2^8 - 1 = 11111111b rounded to four bits rounds up to 2^8.
        const x = BSNumber.fromNumber(255);
        const y = convertBSNumber(x, 4, BSNumberRoundingMode.FE_TONEAREST);
        expect(y.toNumber()).toBe(256);
        expect(y.getUInteger()).toBe(1n);
        expect(y.getBiasedExponent()).toBe(8);
        expect(y.isValid()).toBe(true);
    });

    it('rejects a nonpositive precision and unsupported modes', () => {
        const x = BSNumber.fromNumber(3);
        expect(() => convertBSNumber(x, 0, BSNumberRoundingMode.FE_TONEAREST)).toThrow();
        expect(() => convertBSNumber(x, -1, BSNumberRoundingMode.FE_TONEAREST)).toThrow();
        expect(() => convertBSNumber(x, 1, 17 as BSNumberRoundingMode)).toThrow();
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V05). The port replaces the upstream 32-bit-block
// UInteger with bigint, so these properties check the *semantics* -- sign,
// biased exponent, the odd-integer invariant, round-to-nearest-ties-to-even
// in both directions of conversion, and the rounding modes of Convert --
// against independent bigint computations rather than the storage layout.
// ---------------------------------------------------------------------------
describe('BSNumber verification', () => {
    // Any finite double, including subnormals and the extremes.
    const anyDouble = fc.double({ noNaN: true, noDefaultInfinity: true });

    // A BSNumber whose unsigned integer typically needs far more than 53
    // bits, built by exact arithmetic on ordinary doubles.
    const wideBSNumber: fc.Arbitrary<BSNumber> = fc.tuple(anyDouble, anyDouble, anyDouble)
        .map(([a, b, c]) => BSNumber.fromNumber(a)
            .mul(BSNumber.fromNumber(b))
            .add(BSNumber.fromNumber(c)));

    const anyBSNumber: fc.Arbitrary<BSNumber> =
        fc.oneof(anyDouble.map(d => BSNumber.fromNumber(d)), wideBSNumber);

    function nextUp(v: number): number {
        if (!Number.isFinite(v)) { return v; }
        if (v === 0) { return Number.MIN_VALUE; }
        const buffer = new DataView(new ArrayBuffer(8));
        buffer.setFloat64(0, v);
        buffer.setBigUint64(0, buffer.getBigUint64(0) + (v > 0 ? 1n : -1n));
        return buffer.getFloat64(0);
    }

    function nextDown(v: number): number {
        return -nextUp(-v);
    }

    /** True when the significand of the double is even (the ties-to-even winner). */
    function significandIsEven(v: number): boolean {
        const buffer = new DataView(new ArrayBuffer(8));
        buffer.setFloat64(0, v);
        return (buffer.getBigUint64(0) & 1n) === 0n;
    }

    /** |a - b| as an exact BSNumber. */
    function absDiff(a: BSNumber, b: BSNumber): BSNumber {
        return BSNumber.fabs(a.sub(b));
    }

    // ---- representation invariants ----------------------------------------

    it('every constructed and computed number satisfies the odd-integer invariant', () => {
        check(fc.tuple(anyBSNumber, anyBSNumber), ([a, b]) => {
            const results = [a, b, a.add(b), a.sub(b), a.mul(b), a.negated(),
                BSNumber.fabs(a), BSNumber.frexp(a).result];
            if (a.getSign() !== 0) {
                // ldexp of zero is the one upstream operation that can break
                // the invariant; see the dedicated test below.
                results.push(BSNumber.ldexp(a, 37));
            }
            for (const x of results) {
                expect(x.isValid()).toBe(true);
            }
            return true;
        });
    });

    it('zero has the canonical representation and is absorbing/neutral', () => {
        const zero = new BSNumber();
        check(anyBSNumber, a => {
            expect(a.add(zero).equals(a)).toBe(true);
            expect(zero.add(a).equals(a)).toBe(true);
            expect(a.sub(zero).equals(a)).toBe(true);
            expect(zero.sub(a).equals(a.negated())).toBe(true);
            const p = a.mul(zero);
            expect(p.getSign()).toBe(0);
            expect(p.getBiasedExponent()).toBe(0);
            expect(p.getUInteger()).toBe(0n);
            expect(a.sub(a).getSign()).toBe(0);
            // Regression: the sign is an int32 upstream, so negating zero
            // must yield +0 and not JavaScript's -0.
            expect(Object.is(zero.negated().getSign(), 0)).toBe(true);
            expect(Object.is(zero.sub(a).getSign(), -a.getSign() | 0)).toBe(true);
            return true;
        });
    });

    // ---- exact arithmetic against independent bigint rationals -------------

    it('add, sub and mul agree with exact bigint arithmetic', () => {
        check(fc.tuple(anyBSNumber, anyBSNumber), ([a, b]) => {
            const ea = exact(a);
            const eb = exact(b);
            const minExp = Math.min(ea.exp, eb.exp);
            const na = ea.num << BigInt(ea.exp - minExp);
            const nb = eb.num << BigInt(eb.exp - minExp);

            for (const [result, expectedNum] of [
                [a.add(b), na + nb], [a.sub(b), na - nb]] as Array<[BSNumber, bigint]>) {
                const er = exact(result);
                if (er.num === 0n) {
                    expect(expectedNum).toBe(0n);
                    continue;
                }
                // Compare on the common exponent minExp; a nonzero result's
                // exponent is never smaller because its bits are a subset of
                // the operands' bits.
                expect(er.exp).toBeGreaterThanOrEqual(minExp);
                expect(er.num << BigInt(er.exp - minExp)).toBe(expectedNum);
            }

            const prod = exact(a.mul(b));
            if (ea.num === 0n || eb.num === 0n) {
                expect(prod.num).toBe(0n);
            } else {
                expect(prod.num).toBe(ea.num * eb.num);
                expect(prod.exp).toBe(ea.exp + eb.exp);
            }
            return true;
        });
    });

    it('is a commutative ring on the generated values', () => {
        check(fc.tuple(anyBSNumber, anyBSNumber, anyBSNumber), ([a, b, c]) => {
            expect(a.add(b).equals(b.add(a))).toBe(true);
            expect(a.mul(b).equals(b.mul(a))).toBe(true);
            expect(a.add(b).add(c).equals(a.add(b.add(c)))).toBe(true);
            expect(a.mul(b).mul(c).equals(a.mul(b.mul(c)))).toBe(true);
            // Distributivity is exact only because there is no rounding.
            expect(a.add(b).mul(c).equals(a.mul(c).add(b.mul(c)))).toBe(true);
            expect(a.sub(b).equals(a.add(b.negated()))).toBe(true);
            return true;
        });
    });

    it('robustSOP, robustDOP and fma are the exact expressions', () => {
        check(fc.tuple(anyBSNumber, anyBSNumber, anyBSNumber, anyBSNumber),
            ([u, v, w, z]) => {
                expect(BSNumber.fma(u, v, w).equals(u.mul(v).add(w))).toBe(true);
                expect(BSNumber.robustSOP(u, v, w, z).equals(u.mul(v).add(w.mul(z)))).toBe(true);
                expect(BSNumber.robustDOP(u, v, w, z).equals(u.mul(v).sub(w.mul(z)))).toBe(true);
                return true;
            });
    });

    // ---- ordering ----------------------------------------------------------

    it('the comparison operators realise the exact total order', () => {
        check(fc.tuple(anyBSNumber, anyBSNumber), ([a, b]) => {
            const c = exactCompare(a, b);
            expect(a.lessThan(b)).toBe(c < 0);
            expect(a.greaterThan(b)).toBe(c > 0);
            expect(a.equals(b)).toBe(c === 0);
            expect(a.notEquals(b)).toBe(c !== 0);
            expect(a.lessThanOrEqual(b)).toBe(c <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(c >= 0);
            // Order-reversal under negation.
            expect(a.negated().lessThan(b.negated())).toBe(c > 0);
            return true;
        });
    });

    it('a - b has the sign implied by the comparison', () => {
        check(fc.tuple(anyBSNumber, anyBSNumber), ([a, b]) => {
            expect(a.sub(b).getSign()).toBe(exactCompare(a, b));
            return true;
        });
    });

    // ---- conversions to double --------------------------------------------

    it('fromNumber/toNumber round trips every finite double exactly', () => {
        check(anyDouble, d => {
            const x = BSNumber.fromNumber(d);
            expect(x.toNumber()).toBe(d === 0 ? 0 : d);   // -0 maps to +0
            const e = exact(x);
            const ed = doubleToExact(d);
            const minExp = Math.min(e.exp, ed.exp);
            expect(e.num << BigInt(e.exp - minExp)).toBe(ed.num << BigInt(ed.exp - minExp));
            return true;
        });
    });

    // toNumber must produce the double nearest the exact value, ties to even.
    // This characterises round-to-nearest-ties-to-even without reimplementing
    // the conversion: the returned double must be at least as close as both
    // of its neighbours, and a tie must be broken toward the even significand.
    it('toNumber returns the nearest double with ties to even', () => {
        check(wideBSNumber, x => {
            const d = x.toNumber();
            if (!Number.isFinite(d)) { return true; }   // overflow, checked below
            const bd = BSNumber.fromNumber(d);
            const err = absDiff(x, bd);

            for (const neighbour of [nextUp(d), nextDown(d)]) {
                if (!Number.isFinite(neighbour)) { continue; }
                const errN = absDiff(x, BSNumber.fromNumber(neighbour));
                expect(err.lessThanOrEqual(errN)).toBe(true);
                if (err.equals(errN)) {
                    expect(significandIsEven(d)).toBe(true);
                }
            }
            return true;
        });
    });

    it('toFloat32 agrees with Math.fround on every finite double', () => {
        // BSNumber has no signed zero (upstream too: the sign of the zero
        // representation is 0), so a result that underflows to zero is
        // compared with == rather than Object.is.
        const sameDouble = (a: number, b: number): boolean =>
            Object.is(a, b) || (a === 0 && b === 0);
        check(anyDouble, d => {
            expect(sameDouble(BSNumber.fromNumber(d).toFloat32(), Math.fround(d))).toBe(true);
            if (Number.isFinite(Math.fround(d))) {
                expect(sameDouble(BSNumber.fromFloat32(d).toNumber(), Math.fround(d))).toBe(true);
            } else {
                // A binary32 overflow becomes an infinity encoding, which
                // BSNumber cannot represent; upstream's graceful exit stores
                // (-1)^s * 2^{1 + EXPONENT_BIAS} instead.
                expect(BSNumber.fromFloat32(d).toNumber()).toBe(Math.sign(d) * 2 ** 128);
            }
            return true;
        });
    });

    it('the zero representation is unsigned', () => {
        expect(Object.is(BSNumber.fromNumber(-0).toNumber(), 0)).toBe(true);
        expect(BSNumber.fromNumber(-0).getSign()).toBe(0);
        // A nonzero negative that underflows keeps its sign, as in C++.
        expect(Object.is(BSNumber.fromNumber(-1e-50).toFloat32(), -0)).toBe(true);
    });

    it('toFloat32 of a wide exact value rounds to nearest with ties to even', () => {
        check(wideBSNumber, x => {
            const f = x.toFloat32();
            if (!Number.isFinite(f)) { return true; }
            const bf = BSNumber.fromNumber(f);
            const err = absDiff(x, bf);
            for (const step of [1, -1]) {
                // Neighbouring binary32 values, obtained through the binary32
                // bit pattern.
                const buffer = new DataView(new ArrayBuffer(4));
                buffer.setFloat32(0, f);
                const bits = buffer.getUint32(0);
                const magnitudeBits = f === 0
                    ? (step > 0 ? 1 : 0x80000001)
                    : ((f > 0) === (step > 0) ? bits + 1 : bits - 1);
                buffer.setUint32(0, magnitudeBits >>> 0);
                const neighbour = buffer.getFloat32(0);
                if (!Number.isFinite(neighbour)) { continue; }
                const errN = absDiff(x, BSNumber.fromNumber(neighbour));
                expect(err.lessThanOrEqual(errN)).toBe(true);
            }
            return true;
        });
    });

    it('overflow saturates to infinity and underflow rounds to zero or the min subnormal', () => {
        // 2^1024 is the first power of two above the double range.
        const overflow = BSNumber.fromNumber(2 ** 1023).mul(BSNumber.fromNumber(2));
        expect(overflow.toNumber()).toBe(Infinity);
        expect(overflow.negated().toNumber()).toBe(-Infinity);

        const minSub = BSNumber.fromNumber(Number.MIN_VALUE);
        // Exactly half the min subnormal is a tie; zero has the even
        // significand, so it wins.
        expect(minSub.mul(BSNumber.fromNumber(0.5)).toNumber()).toBe(0);
        // Three quarters is above the tie, so it rounds up.
        expect(minSub.mul(BSNumber.fromNumber(0.75)).toNumber()).toBe(Number.MIN_VALUE);
        // A quarter is below the tie, so it rounds down to zero.
        expect(minSub.mul(BSNumber.fromNumber(0.25)).toNumber()).toBe(0);
        // One and a half is a tie between one and two ulps; two is even.
        expect(minSub.mul(BSNumber.fromNumber(1.5)).toNumber()).toBe(2 * Number.MIN_VALUE);
    });

    // ---- frexp / ldexp -----------------------------------------------------

    it('frexp splits into a significand in [1/2, 1) and an exponent', () => {
        check(anyBSNumber, x => {
            const { result, exponent } = BSNumber.frexp(x);
            if (x.getSign() === 0) {
                expect(exponent).toBe(0);
                expect(result.getSign()).toBe(0);
                return true;
            }
            expect(result.getExponent()).toBe(-1);
            expect(BSNumber.ldexp(result, exponent).equals(x)).toBe(true);
            const magnitude = BSNumber.fabs(result);
            expect(magnitude.greaterThanOrEqual(BSNumber.fromNumber(0.5))).toBe(true);
            expect(magnitude.lessThan(BSNumber.fromNumber(1))).toBe(true);
            return true;
        });
    });

    it('ldexp scales by a power of two exactly', () => {
        check(fc.tuple(anyBSNumber, fc.integer({ min: -200, max: 200 })), ([x, k]) => {
            const scaled = BSNumber.ldexp(x, k);
            if (x.getSign() === 0) { return true; }
            expect(scaled.getUInteger()).toBe(x.getUInteger());
            expect(scaled.getBiasedExponent()).toBe(x.getBiasedExponent() + k);
            expect(BSNumber.ldexp(scaled, -k).equals(x)).toBe(true);
            return true;
        });
    });

    // Upstream quirk, preserved: ldexp of zero shifts the biased exponent of
    // a number whose sign is 0, producing a representation that violates the
    // BSNumber invariant. The value is still zero on conversion.
    it('ldexp of zero produces an invalid (but harmless) representation', () => {
        const z = BSNumber.ldexp(new BSNumber(), 5);
        expect(z.getSign()).toBe(0);
        expect(z.getBiasedExponent()).toBe(5);
        expect(z.isValid()).toBe(false);
        expect(z.toNumber()).toBe(0);
    });

    // ---- strings -----------------------------------------------------------

    it('fromString parses signed decimal integers exactly', () => {
        check(fc.bigInt({ min: -(10n ** 30n), max: 10n ** 30n }), v => {
            const s = (v >= 0n ? v.toString() : v.toString());
            const x = BSNumber.fromString(s);
            expect(x.equals(BSNumber.fromBigInt(v))).toBe(true);
            if (v > 0n) {
                expect(BSNumber.fromString('+' + s).equals(x)).toBe(true);
            }
            return true;
        }, 60);
    });

    // ---- remainder ---------------------------------------------------------

    // Regression: forming the quotient as a double loses exactness above
    // 2^53. std::remainder is exact, so the port computes the quotient with
    // integer arithmetic.
    it('remainder is exact for quotients larger than 2^53', () => {
        const r = BSNumber.remainder(BSNumber.fromNumber(1e17), BSNumber.fromNumber(3));
        expect(r.toNumber()).toBe(1);
        const r2 = BSNumber.remainder(BSNumber.fromNumber(12345678901234568),
            BSNumber.fromNumber(7));
        expect(r2.toNumber()).toBe(1);
        const r3 = BSNumber.remainder(BSNumber.fromNumber(1e22), BSNumber.fromNumber(1e-3));
        expect(r3.toNumber()).toBe(-0.00011721684699608857);
    });

    it('remainder satisfies the IEEE definition against exact bigint arithmetic', () => {
        check(fc.tuple(anyDouble.filter(d => d !== 0), anyDouble.filter(d => d !== 0)),
            ([dx, dy]) => {
                const x = BSNumber.fromNumber(dx);
                const y = BSNumber.fromNumber(dy);
                const r = BSNumber.remainder(x, y);

                // r = x - n*y for some integer n, and 2|r| <= |y| with ties
                // resolved to the even n.
                const ex = exact(x), ey = exact(y), er = exact(r);
                const minExp = Math.min(ex.exp, ey.exp, er.exp);
                const nx = ex.num << BigInt(ex.exp - minExp);
                const ny = ey.num << BigInt(ey.exp - minExp);
                const nr = er.num << BigInt(er.exp - minExp);
                const diff = nx - nr;
                expect(diff % ny).toBe(0n);          // n is an integer
                const n = diff / ny;
                const absY = ny < 0n ? -ny : ny;
                const absR = nr < 0n ? -nr : nr;
                expect(2n * absR).toBeLessThanOrEqual(absY);
                if (2n * absR === absY) {
                    expect(n % 2n).toBe(0n);          // ties to even
                }
                return true;
            }, 100);
    });

    it('remainder degenerates gracefully like the upstream double round trip', () => {
        // std::remainder(x, 0) is NaN and BSNumber(NaN) is zero.
        expect(BSNumber.remainder(BSNumber.fromNumber(5), new BSNumber()).getSign()).toBe(0);
        // std::remainder(0, y) is 0.
        expect(BSNumber.remainder(new BSNumber(), BSNumber.fromNumber(3)).getSign()).toBe(0);
    });

    // ---- convertBSNumber ---------------------------------------------------

    // Independent model of the upstream rounding: keep the leading p bits of
    // the odd magnitude, then apply the mode to the discarded remainder.
    function referenceConvert(input: BSNumber, precision: number,
        mode: BSNumberRoundingMode): BSNumber {
        if (input.getSign() === 0) { return new BSNumber(); }
        const k = input.getNumBits() - precision;
        if (k <= 0) { return input.clone(); }

        const magnitude = input.getUInteger();
        const hi = magnitude >> BigInt(k);
        const rem = magnitude - (hi << BigInt(k));
        const half = 1n << BigInt(k - 1);
        const signValue = input.getSign();

        let rounded: bigint;
        if (mode === BSNumberRoundingMode.FE_TOWARDZERO) {
            rounded = hi;
        } else if (mode === BSNumberRoundingMode.FE_UPWARD) {
            rounded = signValue > 0 ? hi + 1n : hi;
        } else if (mode === BSNumberRoundingMode.FE_DOWNWARD) {
            rounded = signValue < 0 ? hi + 1n : hi;
        } else {
            // The magnitude is odd, so rem == half can only happen for k == 1.
            rounded = rem > half || (rem === half && (hi & 1n) === 1n) ? hi + 1n : hi;
        }

        // Normalise to the odd-integer invariant.
        let shift = 0;
        while (rounded > 0n && (rounded & 1n) === 0n) { rounded >>= 1n; ++shift; }
        return BSNumber.fromParts(signValue, input.getBiasedExponent() + k + shift, rounded);
    }

    const allModes = [
        BSNumberRoundingMode.FE_TONEAREST, BSNumberRoundingMode.FE_DOWNWARD,
        BSNumberRoundingMode.FE_TOWARDZERO, BSNumberRoundingMode.FE_UPWARD
    ];

    it('convertBSNumber matches an independent bigint rounding model', () => {
        check(fc.tuple(anyBSNumber, fc.integer({ min: 1, max: 80 }),
            fc.constantFrom(...allModes)), ([x, precision, mode]) => {
                const got = convertBSNumber(x, precision, mode);
                const want = referenceConvert(x, precision, mode);
                expect(got.getSign()).toBe(want.getSign());
                expect(got.getUInteger()).toBe(want.getUInteger());
                if (got.getSign() !== 0) {
                    expect(got.getBiasedExponent()).toBe(want.getBiasedExponent());
                }
                expect(got.isValid()).toBe(true);
                return true;
            });
    });

    it('convertBSNumber respects the direction of each rounding mode', () => {
        check(fc.tuple(anyBSNumber.filter(x => x.getSign() !== 0),
            fc.integer({ min: 1, max: 60 })), ([x, precision]) => {
                const down = convertBSNumber(x, precision, BSNumberRoundingMode.FE_DOWNWARD);
                const up = convertBSNumber(x, precision, BSNumberRoundingMode.FE_UPWARD);
                const zero = convertBSNumber(x, precision, BSNumberRoundingMode.FE_TOWARDZERO);
                const near = convertBSNumber(x, precision, BSNumberRoundingMode.FE_TONEAREST);

                expect(down.lessThanOrEqual(x)).toBe(true);
                expect(up.greaterThanOrEqual(x)).toBe(true);
                expect(BSNumber.fabs(zero).lessThanOrEqual(BSNumber.fabs(x))).toBe(true);
                // Nearest is between the two directed roundings.
                expect(near.greaterThanOrEqual(down)).toBe(true);
                expect(near.lessThanOrEqual(up)).toBe(true);
                // Nearest is at least as close as either directed rounding.
                expect(absDiff(near, x).lessThanOrEqual(absDiff(down, x))).toBe(true);
                expect(absDiff(near, x).lessThanOrEqual(absDiff(up, x))).toBe(true);
                // The result never uses more than the requested precision.
                expect(near.getNumBits()).toBeLessThanOrEqual(precision);
                expect(down.getNumBits()).toBeLessThanOrEqual(precision);
                expect(up.getNumBits()).toBeLessThanOrEqual(precision);
                expect(zero.getNumBits()).toBeLessThanOrEqual(precision);
                return true;
            });
    });

    it('convertBSNumber is the identity when the precision is already satisfied', () => {
        check(fc.tuple(anyBSNumber, fc.constantFrom(...allModes)), ([x, mode]) => {
            const precision = Math.max(1, x.getNumBits());
            const got = convertBSNumber(x, precision, mode);
            expect(got.equals(x)).toBe(true);
            return true;
        });
    });
});
