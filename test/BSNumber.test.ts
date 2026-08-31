import { describe, it, expect } from 'vitest';
import { BSNumber, BSNumberRoundingMode, convertBSNumber } from '../src/BSNumber';
import { isArbitraryPrecision, hasDivisionOperator } from '../src/TypeTraits';

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
