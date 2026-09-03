import { describe, it, expect } from 'vitest';
import { check, finite, fc } from './helpers/arbitraries.js';
import { IEEEBinary32, IEEEBinary64, IEEEClassification } from '../src/IEEEBinary.js';

describe('IEEEBinary32', () => {
    it('has the upstream-documented constants', () => {
        expect(IEEEBinary32.NUM_ENCODING_BITS).toBe(32);
        expect(IEEEBinary32.NUM_EXPONENT_BITS).toBe(8);
        expect(IEEEBinary32.NUM_SIGNIFICAND_BITS).toBe(24);
        expect(IEEEBinary32.NUM_TRAILING_BITS).toBe(23);
        expect(IEEEBinary32.EXPONENT_BIAS).toBe(127);
        expect(IEEEBinary32.MAX_BIASED_EXPONENT).toBe(255);
        expect(IEEEBinary32.MIN_SUB_EXPONENT).toBe(-126);
        expect(IEEEBinary32.MIN_EXPONENT).toBe(-149);
        expect(IEEEBinary32.SIGN_SHIFT).toBe(31);
        expect(IEEEBinary32.SIGN_MASK).toBe(0x80000000);
        expect(IEEEBinary32.NOT_SIGN_MASK).toBe(0x7FFFFFFF);
        expect(IEEEBinary32.TRAILING_MASK).toBe(0x007FFFFF);
        expect(IEEEBinary32.EXPONENT_MASK).toBe(0x7F800000);
        expect(IEEEBinary32.NAN_QUIET_MASK).toBe(0x00400000);
        expect(IEEEBinary32.NAN_PAYLOAD_MASK).toBe(0x003FFFFF);
        expect(IEEEBinary32.SUP_TRAILING).toBe(0x00800000);
        expect(IEEEBinary32.MIN_NORMAL).toBe(0x00800000);
        expect(IEEEBinary32.MAX_NORMAL).toBe(0x7F7FFFFF);
        expect(IEEEBinary32.POS_INFINITY).toBe(0x7F800000);
        expect(IEEEBinary32.NEG_INFINITY).toBe(0xFF800000);
    });

    it('encodes known bit patterns', () => {
        expect(IEEEBinary32.fromNumber(1).encoding).toBe(0x3F800000);
        expect(IEEEBinary32.fromNumber(-2).encoding).toBe(0xC0000000);
        expect(IEEEBinary32.fromNumber(0.5).encoding).toBe(0x3F000000);
        expect(IEEEBinary32.fromNumber(0.75).encoding).toBe(0x3F400000);
        expect(IEEEBinary32.fromNumber(0).encoding).toBe(IEEEBinary32.POS_ZERO);
        expect(IEEEBinary32.fromNumber(-0).encoding).toBe(IEEEBinary32.NEG_ZERO);
        expect(IEEEBinary32.fromNumber(Infinity).encoding).toBe(IEEEBinary32.POS_INFINITY);
        expect(IEEEBinary32.fromNumber(-Infinity).encoding).toBe(IEEEBinary32.NEG_INFINITY);
    });

    it('decodes known bit patterns', () => {
        expect(IEEEBinary32.fromEncoding(0x40490FDB).number).toBe(Math.fround(Math.PI));
        expect(IEEEBinary32.fromEncoding(IEEEBinary32.MAX_NORMAL).number)
            .toBe(3.4028234663852886e38);
        expect(IEEEBinary32.fromEncoding(IEEEBinary32.MIN_SUBNORMAL).number)
            .toBe(2 ** -149);
        expect(IEEEBinary32.fromEncoding(IEEEBinary32.MIN_NORMAL).number)
            .toBe(2 ** -126);
    });

    it('rounds to binary32 on assignment', () => {
        const x = IEEEBinary32.fromNumber(0.1);
        expect(x.number).toBe(Math.fround(0.1));
        expect(x.number).not.toBe(0.1);
    });

    it('round trips number <-> encoding for random values', () => {
        for (let i = 0; i < 1000; ++i) {
            const value = Math.fround((Math.random() - 0.5) * 2 ** (Math.random() * 60 - 30));
            const x = IEEEBinary32.fromNumber(value);
            expect(IEEEBinary32.fromEncoding(x.encoding).number).toBe(value);
        }
    });

    it('decomposes and recomposes sign/biased/trailing', () => {
        const x = IEEEBinary32.fromNumber(-6.5); // -1.625 * 2^2
        const { sign, biased, trailing } = x.getEncoding();
        expect(sign).toBe(1);
        expect(biased).toBe(129); // exponent 2 + bias 127
        expect(trailing).toBe(0b10100000000000000000000);
        const y = IEEEBinary32.fromParts(sign, biased, trailing);
        expect(y.encoding).toBe(x.encoding);
        expect(y.number).toBe(-6.5);
    });

    it('classifies numbers', () => {
        expect(IEEEBinary32.fromNumber(0).getClassification())
            .toBe(IEEEClassification.POS_ZERO);
        expect(IEEEBinary32.fromNumber(-0).getClassification())
            .toBe(IEEEClassification.NEG_ZERO);
        expect(IEEEBinary32.fromNumber(1.5).getClassification())
            .toBe(IEEEClassification.POS_NORMAL);
        expect(IEEEBinary32.fromNumber(-1.5).getClassification())
            .toBe(IEEEClassification.NEG_NORMAL);
        expect(IEEEBinary32.fromEncoding(IEEEBinary32.MAX_SUBNORMAL).getClassification())
            .toBe(IEEEClassification.POS_SUBNORMAL);
        expect(IEEEBinary32.fromEncoding((IEEEBinary32.SIGN_MASK | 1) >>> 0).getClassification())
            .toBe(IEEEClassification.NEG_SUBNORMAL);
        expect(IEEEBinary32.fromNumber(Infinity).getClassification())
            .toBe(IEEEClassification.POS_INFINITY);
        expect(IEEEBinary32.fromNumber(-Infinity).getClassification())
            .toBe(IEEEClassification.NEG_INFINITY);
        expect(IEEEBinary32.fromNumber(NaN).getClassification())
            .toBe(IEEEClassification.QUIET_NAN);
        // A signaling NaN: max exponent, quiet bit clear, nonzero payload.
        expect(IEEEBinary32.fromEncoding(0x7F800001).getClassification())
            .toBe(IEEEClassification.SIGNALING_NAN);
    });

    it('implements the Is* predicates', () => {
        const one = IEEEBinary32.fromNumber(1);
        expect(one.isNormal()).toBe(true);
        expect(one.isFinite()).toBe(true);
        expect(one.isZero()).toBe(false);
        expect(one.isSubnormal()).toBe(false);
        expect(one.isSignMinus()).toBe(false);
        expect(one.isNaN()).toBe(false);

        expect(IEEEBinary32.fromNumber(-0).isZero()).toBe(true);
        expect(IEEEBinary32.fromNumber(-0).isSignMinus()).toBe(true);
        expect(IEEEBinary32.fromEncoding(1).isSubnormal()).toBe(true);
        expect(IEEEBinary32.fromNumber(Infinity).isInfinite()).toBe(true);
        expect(IEEEBinary32.fromNumber(Infinity).isFinite()).toBe(false);
        expect(IEEEBinary32.fromNumber(NaN).isNaN()).toBe(true);
        expect(IEEEBinary32.fromNumber(NaN).isQuietNaN()).toBe(true);
        expect(IEEEBinary32.fromEncoding(0x7F800001).isSignalingNaN()).toBe(true);
    });

    it('computes next-up encodings', () => {
        expect(IEEEBinary32.fromNumber(1).getNextUp()).toBe(0x3F800001);
        expect(IEEEBinary32.fromNumber(0).getNextUp()).toBe(IEEEBinary32.MIN_SUBNORMAL);
        expect(IEEEBinary32.fromNumber(-0).getNextUp()).toBe(IEEEBinary32.MIN_SUBNORMAL);
        // -MIN_SUBNORMAL -> -0.
        expect(IEEEBinary32.fromEncoding((IEEEBinary32.SIGN_MASK | 1) >>> 0).getNextUp())
            .toBe(IEEEBinary32.NEG_ZERO);
        // MAX_SUBNORMAL -> MIN_NORMAL.
        expect(IEEEBinary32.fromEncoding(IEEEBinary32.MAX_SUBNORMAL).getNextUp())
            .toBe(IEEEBinary32.MIN_NORMAL);
        // -INFINITY -> -MAX_NORMAL.
        expect(IEEEBinary32.fromNumber(-Infinity).getNextUp())
            .toBe((IEEEBinary32.SIGN_MASK | IEEEBinary32.MAX_NORMAL) >>> 0);
        // +INFINITY -> +INFINITY.
        expect(IEEEBinary32.fromNumber(Infinity).getNextUp())
            .toBe(IEEEBinary32.POS_INFINITY);
        // NaN -> itself.
        const nan = IEEEBinary32.fromNumber(NaN);
        expect(nan.getNextUp()).toBe(nan.encoding);
    });

    it('computes next-down encodings', () => {
        expect(IEEEBinary32.fromNumber(1).getNextDown()).toBe(0x3F7FFFFF);
        expect(IEEEBinary32.fromNumber(0).getNextDown())
            .toBe((IEEEBinary32.SIGN_MASK | IEEEBinary32.MIN_SUBNORMAL) >>> 0);
        // MIN_SUBNORMAL -> +0.
        expect(IEEEBinary32.fromEncoding(1).getNextDown()).toBe(IEEEBinary32.POS_ZERO);
        // -MAX_SUBNORMAL -> -MIN_NORMAL.
        expect(IEEEBinary32.fromEncoding(
            (IEEEBinary32.SIGN_MASK | IEEEBinary32.MAX_SUBNORMAL) >>> 0).getNextDown())
            .toBe((IEEEBinary32.SIGN_MASK | IEEEBinary32.MIN_NORMAL) >>> 0);
        // +INFINITY -> MAX_NORMAL, -INFINITY -> -INFINITY.
        expect(IEEEBinary32.fromNumber(Infinity).getNextDown()).toBe(IEEEBinary32.MAX_NORMAL);
        expect(IEEEBinary32.fromNumber(-Infinity).getNextDown()).toBe(IEEEBinary32.NEG_INFINITY);
    });

    it('next-up/next-down are inverses on random finite values', () => {
        for (let i = 0; i < 500; ++i) {
            const x = IEEEBinary32.fromNumber(
                Math.fround((Math.random() - 0.5) * 1e6));
            if (x.isZero()) {
                continue; // next-down(next-up(+/-0)) is +0, losing the sign of -0
            }
            const up = IEEEBinary32.fromEncoding(x.getNextUp());
            expect(up.getNextDown()).toBe(x.encoding);
            expect(up.number).toBeGreaterThan(x.number);
        }
    });
});

describe('IEEEBinary64', () => {
    it('has the upstream-documented constants', () => {
        expect(IEEEBinary64.NUM_ENCODING_BITS).toBe(64);
        expect(IEEEBinary64.NUM_EXPONENT_BITS).toBe(11);
        expect(IEEEBinary64.NUM_SIGNIFICAND_BITS).toBe(53);
        expect(IEEEBinary64.NUM_TRAILING_BITS).toBe(52);
        expect(IEEEBinary64.EXPONENT_BIAS).toBe(1023);
        expect(IEEEBinary64.MAX_BIASED_EXPONENT).toBe(2047);
        expect(IEEEBinary64.MIN_SUB_EXPONENT).toBe(-1022);
        expect(IEEEBinary64.MIN_EXPONENT).toBe(-1074);
        expect(IEEEBinary64.SIGN_SHIFT).toBe(63);
        expect(IEEEBinary64.SIGN_MASK).toBe(0x8000000000000000n);
        expect(IEEEBinary64.NOT_SIGN_MASK).toBe(0x7FFFFFFFFFFFFFFFn);
        expect(IEEEBinary64.TRAILING_MASK).toBe(0x000FFFFFFFFFFFFFn);
        expect(IEEEBinary64.EXPONENT_MASK).toBe(0x7FF0000000000000n);
        expect(IEEEBinary64.NAN_QUIET_MASK).toBe(0x0008000000000000n);
        expect(IEEEBinary64.NAN_PAYLOAD_MASK).toBe(0x0007FFFFFFFFFFFFn);
        expect(IEEEBinary64.SUP_TRAILING).toBe(0x0010000000000000n);
        expect(IEEEBinary64.MIN_NORMAL).toBe(0x0010000000000000n);
        expect(IEEEBinary64.MAX_NORMAL).toBe(0x7FEFFFFFFFFFFFFFn);
        expect(IEEEBinary64.POS_INFINITY).toBe(0x7FF0000000000000n);
        expect(IEEEBinary64.NEG_INFINITY).toBe(0xFFF0000000000000n);
    });

    it('encodes known bit patterns', () => {
        expect(IEEEBinary64.fromNumber(1).encoding).toBe(0x3FF0000000000000n);
        expect(IEEEBinary64.fromNumber(-2).encoding).toBe(0xC000000000000000n);
        expect(IEEEBinary64.fromNumber(0.5).encoding).toBe(0x3FE0000000000000n);
        expect(IEEEBinary64.fromNumber(Math.PI).encoding).toBe(0x400921FB54442D18n);
        expect(IEEEBinary64.fromNumber(0).encoding).toBe(IEEEBinary64.POS_ZERO);
        expect(IEEEBinary64.fromNumber(-0).encoding).toBe(IEEEBinary64.NEG_ZERO);
        expect(IEEEBinary64.fromNumber(Infinity).encoding).toBe(IEEEBinary64.POS_INFINITY);
        expect(IEEEBinary64.fromNumber(-Infinity).encoding).toBe(IEEEBinary64.NEG_INFINITY);
    });

    it('decodes special encodings to the JavaScript limit constants', () => {
        expect(IEEEBinary64.fromEncoding(IEEEBinary64.MAX_NORMAL).number)
            .toBe(Number.MAX_VALUE);
        expect(IEEEBinary64.fromEncoding(IEEEBinary64.MIN_SUBNORMAL).number)
            .toBe(Number.MIN_VALUE); // 5e-324
        expect(IEEEBinary64.fromEncoding(IEEEBinary64.MIN_NORMAL).number)
            .toBe(2 ** -1022);
        // 1 + 2^-52 is the double after 1.
        expect(IEEEBinary64.fromEncoding(0x3FF0000000000001n).number)
            .toBe(1 + Number.EPSILON);
    });

    it('round trips number <-> encoding for random values', () => {
        for (let i = 0; i < 1000; ++i) {
            const value = (Math.random() - 0.5) * 2 ** (Math.random() * 120 - 60);
            const x = IEEEBinary64.fromNumber(value);
            expect(IEEEBinary64.fromEncoding(x.encoding).number).toBe(value);
        }
    });

    it('decomposes and recomposes sign/biased/trailing', () => {
        const x = IEEEBinary64.fromNumber(-6.5); // -1.625 * 2^2
        const { sign, biased, trailing } = x.getEncoding();
        expect(sign).toBe(1);
        expect(biased).toBe(1025); // exponent 2 + bias 1023
        expect(trailing).toBe(0b1010n << 48n);
        const y = IEEEBinary64.fromParts(sign, biased, trailing);
        expect(y.encoding).toBe(x.encoding);
        expect(y.number).toBe(-6.5);
    });

    it('classifies numbers', () => {
        expect(IEEEBinary64.fromNumber(0).getClassification())
            .toBe(IEEEClassification.POS_ZERO);
        expect(IEEEBinary64.fromNumber(-0).getClassification())
            .toBe(IEEEClassification.NEG_ZERO);
        expect(IEEEBinary64.fromNumber(Math.E).getClassification())
            .toBe(IEEEClassification.POS_NORMAL);
        expect(IEEEBinary64.fromNumber(-Math.E).getClassification())
            .toBe(IEEEClassification.NEG_NORMAL);
        expect(IEEEBinary64.fromNumber(Number.MIN_VALUE).getClassification())
            .toBe(IEEEClassification.POS_SUBNORMAL);
        expect(IEEEBinary64.fromNumber(-Number.MIN_VALUE).getClassification())
            .toBe(IEEEClassification.NEG_SUBNORMAL);
        expect(IEEEBinary64.fromNumber(Infinity).getClassification())
            .toBe(IEEEClassification.POS_INFINITY);
        expect(IEEEBinary64.fromNumber(-Infinity).getClassification())
            .toBe(IEEEClassification.NEG_INFINITY);
        // The JavaScript NaN is the quiet NaN 0x7FF8000000000000.
        expect(IEEEBinary64.fromNumber(NaN).getClassification())
            .toBe(IEEEClassification.QUIET_NAN);
        expect(IEEEBinary64.fromEncoding(0x7FF0000000000001n).getClassification())
            .toBe(IEEEClassification.SIGNALING_NAN);
    });

    it('implements the Is* predicates', () => {
        expect(IEEEBinary64.fromNumber(1).isNormal()).toBe(true);
        expect(IEEEBinary64.fromNumber(-0).isZero()).toBe(true);
        expect(IEEEBinary64.fromNumber(-0).isSignMinus()).toBe(true);
        expect(IEEEBinary64.fromNumber(Number.MIN_VALUE).isSubnormal()).toBe(true);
        expect(IEEEBinary64.fromNumber(Number.MIN_VALUE).isNormal()).toBe(false);
        expect(IEEEBinary64.fromNumber(Infinity).isInfinite()).toBe(true);
        expect(IEEEBinary64.fromNumber(Infinity).isFinite()).toBe(false);
        expect(IEEEBinary64.fromNumber(NaN).isNaN()).toBe(true);
        expect(IEEEBinary64.fromNumber(NaN).isQuietNaN()).toBe(true);
        expect(IEEEBinary64.fromNumber(NaN).isSignalingNaN()).toBe(false);
        expect(IEEEBinary64.fromEncoding(0x7FF0000000000001n).isSignalingNaN()).toBe(true);
    });

    it('computes next-up and next-down encodings', () => {
        // next-up(1) is 1 + 2^-52.
        expect(IEEEBinary64.fromNumber(1).getNextUp()).toBe(0x3FF0000000000001n);
        expect(IEEEBinary64.fromEncoding(IEEEBinary64.fromNumber(1).getNextUp()).number)
            .toBe(1 + Number.EPSILON);
        expect(IEEEBinary64.fromNumber(0).getNextUp()).toBe(IEEEBinary64.MIN_SUBNORMAL);
        expect(IEEEBinary64.fromNumber(-0).getNextUp()).toBe(IEEEBinary64.MIN_SUBNORMAL);
        expect(IEEEBinary64.fromNumber(0).getNextDown())
            .toBe(IEEEBinary64.SIGN_MASK | IEEEBinary64.MIN_SUBNORMAL);
        // MAX_SUBNORMAL -> MIN_NORMAL.
        expect(IEEEBinary64.fromEncoding(IEEEBinary64.MAX_SUBNORMAL).getNextUp())
            .toBe(IEEEBinary64.MIN_NORMAL);
        // -MIN_SUBNORMAL -> -0.
        expect(IEEEBinary64.fromEncoding(
            IEEEBinary64.SIGN_MASK | IEEEBinary64.MIN_SUBNORMAL).getNextUp())
            .toBe(IEEEBinary64.NEG_ZERO);
        // Infinities.
        expect(IEEEBinary64.fromNumber(Infinity).getNextUp()).toBe(IEEEBinary64.POS_INFINITY);
        expect(IEEEBinary64.fromNumber(Infinity).getNextDown()).toBe(IEEEBinary64.MAX_NORMAL);
        expect(IEEEBinary64.fromNumber(-Infinity).getNextUp())
            .toBe(IEEEBinary64.SIGN_MASK | IEEEBinary64.MAX_NORMAL);
        expect(IEEEBinary64.fromNumber(-Infinity).getNextDown()).toBe(IEEEBinary64.NEG_INFINITY);
        // NaN -> itself.
        const nan = IEEEBinary64.fromNumber(NaN);
        expect(nan.getNextUp()).toBe(nan.encoding);
        expect(nan.getNextDown()).toBe(nan.encoding);
    });

    it('next-up/next-down are inverses on random finite values', () => {
        for (let i = 0; i < 500; ++i) {
            const x = IEEEBinary64.fromNumber((Math.random() - 0.5) * 1e12);
            if (x.isZero()) {
                continue;
            }
            const up = IEEEBinary64.fromEncoding(x.getNextUp());
            expect(up.getNextDown()).toBe(x.encoding);
            expect(up.number).toBeGreaterThan(x.number);
            const down = IEEEBinary64.fromEncoding(x.getNextDown());
            expect(down.getNextUp()).toBe(x.encoding);
            expect(down.number).toBeLessThan(x.number);
        }
    });
});

describe('IEEEBinary verification', () => {
    // ---- independent references ------------------------------------------
    const dv = new DataView(new ArrayBuffer(8));

    function bits32(x: number): number {
        dv.setFloat32(0, x);
        return dv.getUint32(0);
    }
    function value32(encoding: number): number {
        dv.setUint32(0, encoding >>> 0);
        return dv.getFloat32(0);
    }
    function bits64(x: number): bigint {
        dv.setFloat64(0, x);
        return dv.getBigUint64(0);
    }
    function value64(encoding: bigint): number {
        dv.setBigUint64(0, BigInt.asUintN(64, encoding));
        return dv.getFloat64(0);
    }

    // Doubles spanning normals, subnormals and the extremes.
    const anyDouble = fc.oneof(
        finite(-1e6, 1e6),
        fc.tuple(finite(-2, 2), fc.integer({ min: -1080, max: 1020 }))
            .map(([m, e]) => m * 2 ** e),
        fc.constantFrom(0, -0, 1, -1, Number.MIN_VALUE, -Number.MIN_VALUE,
            Number.MAX_VALUE, -Number.MAX_VALUE, 2 ** -1022, -(2 ** -1022),
            Number.POSITIVE_INFINITY, Number.NEGATIVE_INFINITY, Number.EPSILON));

    const anyEncoding64 = fc.bigInt({ min: 0n, max: (1n << 64n) - 1n });
    const anyEncoding32 = fc.integer({ min: 0, max: 0xFFFFFFFF });

    // ---- constants -------------------------------------------------------
    it('the constants follow the upstream template formulas', () => {
        const derive = (numBits: number, precision: number) => {
            const numExponentBits = numBits - precision;
            const numTrailingBits = precision - 1;
            return {
                NUM_ENCODING_BITS: numBits,
                NUM_EXPONENT_BITS: numExponentBits,
                NUM_SIGNIFICAND_BITS: precision,
                NUM_TRAILING_BITS: numTrailingBits,
                EXPONENT_BIAS: 2 ** (numExponentBits - 1) - 1,
                MAX_BIASED_EXPONENT: 2 ** numExponentBits - 1,
                MIN_SUB_EXPONENT: 1 - (2 ** (numExponentBits - 1) - 1),
                MIN_EXPONENT: 1 - (2 ** (numExponentBits - 1) - 1) - numTrailingBits,
                SIGN_SHIFT: numBits - 1
            };
        };
        for (const [cls, numBits, precision] of
            [[IEEEBinary32, 32, 24], [IEEEBinary64, 64, 53]] as const) {
            const d = derive(numBits, precision);
            expect(cls.NUM_ENCODING_BITS).toBe(d.NUM_ENCODING_BITS);
            expect(cls.NUM_EXPONENT_BITS).toBe(d.NUM_EXPONENT_BITS);
            expect(cls.NUM_SIGNIFICAND_BITS).toBe(d.NUM_SIGNIFICAND_BITS);
            expect(cls.NUM_TRAILING_BITS).toBe(d.NUM_TRAILING_BITS);
            expect(cls.EXPONENT_BIAS).toBe(d.EXPONENT_BIAS);
            expect(cls.MAX_BIASED_EXPONENT).toBe(d.MAX_BIASED_EXPONENT);
            expect(cls.MIN_SUB_EXPONENT).toBe(d.MIN_SUB_EXPONENT);
            expect(cls.MIN_EXPONENT).toBe(d.MIN_EXPONENT);
            expect(cls.SIGN_SHIFT).toBe(d.SIGN_SHIFT);
        }

        // The mask constants, derived the same way as upstream.
        expect(IEEEBinary32.SIGN_MASK).toBe(2 ** 31);
        expect(IEEEBinary32.NOT_SIGN_MASK).toBe(2 ** 31 - 1);
        expect(IEEEBinary32.TRAILING_MASK).toBe(2 ** 23 - 1);
        expect(IEEEBinary32.EXPONENT_MASK).toBe((2 ** 31 - 1) & ~(2 ** 23 - 1));
        expect(IEEEBinary32.NAN_QUIET_MASK).toBe(2 ** 22);
        expect(IEEEBinary32.NAN_PAYLOAD_MASK).toBe((2 ** 23 - 1) >> 1);
        expect(IEEEBinary32.MAX_TRAILING).toBe(IEEEBinary32.TRAILING_MASK);
        expect(IEEEBinary32.SUP_TRAILING).toBe(2 ** 23);
        expect(IEEEBinary32.MAX_SUBNORMAL).toBe(IEEEBinary32.TRAILING_MASK);
        expect(IEEEBinary32.MIN_NORMAL).toBe(IEEEBinary32.SUP_TRAILING);
        expect(IEEEBinary32.MAX_NORMAL)
            .toBe((IEEEBinary32.NOT_SIGN_MASK & ~IEEEBinary32.SUP_TRAILING) >>> 0);
        expect(IEEEBinary32.POS_INFINITY).toBe(IEEEBinary32.EXPONENT_MASK);
        expect(IEEEBinary32.NEG_INFINITY)
            .toBe((IEEEBinary32.SIGN_MASK | IEEEBinary32.EXPONENT_MASK) >>> 0);

        expect(IEEEBinary64.SIGN_MASK).toBe(1n << 63n);
        expect(IEEEBinary64.NOT_SIGN_MASK).toBe((1n << 63n) - 1n);
        expect(IEEEBinary64.TRAILING_MASK).toBe((1n << 52n) - 1n);
        expect(IEEEBinary64.EXPONENT_MASK)
            .toBe(IEEEBinary64.NOT_SIGN_MASK & ~IEEEBinary64.TRAILING_MASK);
        expect(IEEEBinary64.NAN_QUIET_MASK).toBe(1n << 51n);
        expect(IEEEBinary64.NAN_PAYLOAD_MASK).toBe(IEEEBinary64.TRAILING_MASK >> 1n);
        expect(IEEEBinary64.MAX_TRAILING).toBe(IEEEBinary64.TRAILING_MASK);
        expect(IEEEBinary64.SUP_TRAILING).toBe(1n << 52n);
        expect(IEEEBinary64.MAX_SUBNORMAL).toBe(IEEEBinary64.TRAILING_MASK);
        expect(IEEEBinary64.MIN_NORMAL).toBe(IEEEBinary64.SUP_TRAILING);
        expect(IEEEBinary64.MAX_NORMAL)
            .toBe(IEEEBinary64.NOT_SIGN_MASK & ~IEEEBinary64.SUP_TRAILING);
        expect(IEEEBinary64.POS_INFINITY).toBe(IEEEBinary64.EXPONENT_MASK);
        expect(IEEEBinary64.NEG_INFINITY)
            .toBe(IEEEBinary64.SIGN_MASK | IEEEBinary64.EXPONENT_MASK);
    });

    // ---- the union: number <-> encoding -----------------------------------
    it('binary64 encoding matches a DataView of the same double', () => {
        check(anyDouble, x => {
            const b = IEEEBinary64.fromNumber(x);
            return b.encoding === bits64(x) && Object.is(b.number, x);
        });
    });

    it('binary64 decoding matches a DataView of the same bit pattern', () => {
        check(anyEncoding64, e => {
            const b = IEEEBinary64.fromEncoding(e);
            const expected = value64(e);
            return Object.is(b.number, expected)
                || (Number.isNaN(b.number) && Number.isNaN(expected));
        });
    });

    it('binary32 encoding matches a DataView of the rounded float', () => {
        check(anyDouble, x => {
            const b = IEEEBinary32.fromNumber(x);
            return b.encoding === bits32(x) && Object.is(b.number, Math.fround(x));
        });
    });

    it('binary32 decoding matches a DataView of the same bit pattern', () => {
        check(anyEncoding32, e => {
            const b = IEEEBinary32.fromEncoding(e);
            const expected = value32(e);
            return Object.is(b.number, expected)
                || (Number.isNaN(b.number) && Number.isNaN(expected));
        });
    });

    // ---- sign / biased / trailing ----------------------------------------
    it('binary64 getEncoding and setEncoding round trip', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 1 }), fc.integer({ min: 0, max: 2047 }),
            fc.bigInt({ min: 0n, max: (1n << 52n) - 1n })),
            ([sign, biased, trailing]) => {
                const b = IEEEBinary64.fromParts(sign, biased, trailing);
                const e = b.getEncoding();
                return e.sign === sign && e.biased === biased && e.trailing === trailing
                    && b.encoding === ((BigInt(sign) << 63n) | (BigInt(biased) << 52n)
                        | trailing);
            });
    });

    it('binary32 getEncoding and setEncoding round trip', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 1 }), fc.integer({ min: 0, max: 255 }),
            fc.integer({ min: 0, max: (1 << 23) - 1 })),
            ([sign, biased, trailing]) => {
                const b = IEEEBinary32.fromParts(sign, biased, trailing);
                const e = b.getEncoding();
                return e.sign === sign && e.biased === biased && e.trailing === trailing
                    && b.encoding === (((sign << 31) | (biased << 23) | trailing) >>> 0);
            });
    });

    it('the three fields reconstruct the whole encoding', () => {
        check(anyEncoding64, e => {
            const b = IEEEBinary64.fromEncoding(e);
            const p = b.getEncoding();
            return IEEEBinary64.fromParts(p.sign, p.biased, p.trailing).encoding === e;
        });
        check(anyEncoding32, e => {
            const b = IEEEBinary32.fromEncoding(e);
            const p = b.getEncoding();
            return IEEEBinary32.fromParts(p.sign, p.biased, p.trailing).encoding === (e >>> 0);
        });
    });

    // ---- classification ---------------------------------------------------
    it('binary64 classification agrees with the value it decodes to', () => {
        check(anyEncoding64, e => {
            const b = IEEEBinary64.fromEncoding(e);
            const c = b.getClassification();
            const x = b.number;
            const negative = (e & IEEEBinary64.SIGN_MASK) !== 0n;
            if (Number.isNaN(x)) {
                return c === IEEEClassification.QUIET_NAN
                    || c === IEEEClassification.SIGNALING_NAN;
            }
            if (x === Number.POSITIVE_INFINITY) {
                return c === IEEEClassification.POS_INFINITY;
            }
            if (x === Number.NEGATIVE_INFINITY) {
                return c === IEEEClassification.NEG_INFINITY;
            }
            if (x === 0) {
                return c === (negative ? IEEEClassification.NEG_ZERO
                    : IEEEClassification.POS_ZERO);
            }
            const subnormal = Math.abs(x) < 2 ** -1022;
            if (subnormal) {
                return c === (negative ? IEEEClassification.NEG_SUBNORMAL
                    : IEEEClassification.POS_SUBNORMAL);
            }
            return c === (negative ? IEEEClassification.NEG_NORMAL
                : IEEEClassification.POS_NORMAL);
        });
    });

    it('the Is* predicates partition every encoding exactly as the classification', () => {
        const expectConsistent = (b: IEEEBinary32 | IEEEBinary64) => {
            const c = b.getClassification();
            const zero = c === IEEEClassification.NEG_ZERO || c === IEEEClassification.POS_ZERO;
            const sub = c === IEEEClassification.NEG_SUBNORMAL
                || c === IEEEClassification.POS_SUBNORMAL;
            const normal = c === IEEEClassification.NEG_NORMAL
                || c === IEEEClassification.POS_NORMAL;
            const inf = c === IEEEClassification.NEG_INFINITY
                || c === IEEEClassification.POS_INFINITY;
            const nan = c === IEEEClassification.QUIET_NAN
                || c === IEEEClassification.SIGNALING_NAN;
            expect(b.isZero()).toBe(zero);
            expect(b.isSubnormal()).toBe(sub);
            expect(b.isNormal()).toBe(normal);
            expect(b.isInfinite()).toBe(inf);
            expect(b.isNaN()).toBe(nan);
            expect(b.isFinite()).toBe(zero || sub || normal);
            // Exactly one category holds.
            expect([zero, sub, normal, inf, nan].filter(Boolean).length).toBe(1);
            // A NaN is quiet or signaling, and never both; upstream treats a
            // NaN whose payload bits are all zero and whose quiet bit is zero
            // as neither (that encoding is an infinity, so it cannot occur).
            if (nan) {
                expect(b.isQuietNaN() !== b.isSignalingNaN()).toBe(true);
                expect(b.isQuietNaN()).toBe(c === IEEEClassification.QUIET_NAN);
            } else {
                expect(b.isQuietNaN()).toBe(false);
                expect(b.isSignalingNaN()).toBe(false);
            }
        };
        check(anyEncoding64, e => { expectConsistent(IEEEBinary64.fromEncoding(e)); });
        check(anyEncoding32, e => { expectConsistent(IEEEBinary32.fromEncoding(e)); });
    });

    it('isSignMinus is the sign bit for every encoding, zeros included', () => {
        check(anyEncoding64, e =>
            IEEEBinary64.fromEncoding(e).isSignMinus() === ((e & IEEEBinary64.SIGN_MASK) !== 0n));
        check(anyEncoding32, e =>
            IEEEBinary32.fromEncoding(e).isSignMinus() === (((e >>> 31) & 1) === 1));
    });

    // ---- neighbours -------------------------------------------------------
    it('binary64 next-up is the least double greater than the value', () => {
        check(anyDouble.filter(x => Number.isFinite(x)), x => {
            const b = IEEEBinary64.fromNumber(x);
            const up = value64(b.getNextUp());
            expect(up > x || (x === 0 && up === Number.MIN_VALUE)).toBe(true);
            // Nothing lies strictly between x and up: the midpoint rounds to
            // one of the two endpoints.
            const mid = x / 2 + up / 2;
            expect(mid === x || mid === up).toBe(true);
        });
    });

    it('binary64 next-down is the greatest double less than the value', () => {
        check(anyDouble.filter(x => Number.isFinite(x)), x => {
            const b = IEEEBinary64.fromNumber(x);
            const down = value64(b.getNextDown());
            expect(down < x || (x === 0 && down === -Number.MIN_VALUE)).toBe(true);
            const mid = x / 2 + down / 2;
            expect(mid === x || mid === down).toBe(true);
        });
    });

    it('next-up and next-down invert each other away from the boundaries', () => {
        check(anyEncoding64, e => {
            const b = IEEEBinary64.fromEncoding(e);
            if (!b.isFinite()) { return true; }
            const up = b.getNextUp();
            const back = IEEEBinary64.fromEncoding(up).getNextDown();
            // Round trips exactly, except that both zeros step up to
            // MIN_SUBNORMAL, whose next-down is +0.
            return back === e || (b.isZero() && back === IEEEBinary64.POS_ZERO);
        });
        check(anyEncoding32, e => {
            const b = IEEEBinary32.fromEncoding(e);
            if (!b.isFinite()) { return true; }
            const down = b.getNextDown();
            const back = IEEEBinary32.fromEncoding(down).getNextUp();
            return back === (e >>> 0)
                || (b.isZero() && back === (IEEEBinary32.SIGN_MASK >>> 0));
        });
    });

    it('the neighbour steps match the monotone bit ordering', () => {
        // For a nonnegative finite encoding the next-up adds one to the bit
        // pattern; for a negative finite encoding it subtracts one from the
        // magnitude bits. This is the classic characterisation, independent
        // of the upstream branch structure.
        check(anyEncoding64, e => {
            const b = IEEEBinary64.fromEncoding(e);
            if (!b.isFinite()) { return true; }
            const magnitude = e & IEEEBinary64.NOT_SIGN_MASK;
            const negative = (e & IEEEBinary64.SIGN_MASK) !== 0n;
            const expectedUp = negative
                ? (magnitude === 0n ? IEEEBinary64.MIN_SUBNORMAL
                    : IEEEBinary64.SIGN_MASK | (magnitude - 1n))
                : magnitude + 1n;
            const expectedDown = negative
                ? IEEEBinary64.SIGN_MASK | (magnitude + 1n)
                : (magnitude === 0n ? IEEEBinary64.SIGN_MASK | IEEEBinary64.MIN_SUBNORMAL
                    : magnitude - 1n);
            return b.getNextUp() === expectedUp && b.getNextDown() === expectedDown;
        });
    });

    it('the infinities and NaNs follow the documented neighbour rules', () => {
        const posInf64 = IEEEBinary64.fromEncoding(IEEEBinary64.POS_INFINITY);
        const negInf64 = IEEEBinary64.fromEncoding(IEEEBinary64.NEG_INFINITY);
        expect(posInf64.getNextUp()).toBe(IEEEBinary64.POS_INFINITY);
        expect(posInf64.getNextDown()).toBe(IEEEBinary64.MAX_NORMAL);
        expect(negInf64.getNextDown()).toBe(IEEEBinary64.NEG_INFINITY);
        expect(negInf64.getNextUp())
            .toBe(IEEEBinary64.SIGN_MASK | IEEEBinary64.MAX_NORMAL);

        const posInf32 = IEEEBinary32.fromEncoding(IEEEBinary32.POS_INFINITY);
        const negInf32 = IEEEBinary32.fromEncoding(IEEEBinary32.NEG_INFINITY);
        expect(posInf32.getNextUp()).toBe(IEEEBinary32.POS_INFINITY);
        expect(posInf32.getNextDown()).toBe(IEEEBinary32.MAX_NORMAL);
        expect(negInf32.getNextDown()).toBe(IEEEBinary32.NEG_INFINITY);
        expect(negInf32.getNextUp())
            .toBe((IEEEBinary32.SIGN_MASK | IEEEBinary32.MAX_NORMAL) >>> 0);

        // Every NaN encoding is its own next-up and next-down.
        check(fc.bigInt({ min: 1n, max: (1n << 52n) - 1n }), trailing => {
            for (const sign of [0, 1]) {
                const b = IEEEBinary64.fromParts(sign, 2047, trailing);
                if (b.getNextUp() !== b.encoding || b.getNextDown() !== b.encoding) {
                    return false;
                }
            }
            return true;
        });
    });

    it('stepping up from -MIN_SUBNORMAL reaches -0, and down from MIN_SUBNORMAL reaches +0',
        () => {
            const negMin = IEEEBinary64.fromEncoding(
                IEEEBinary64.SIGN_MASK | IEEEBinary64.MIN_SUBNORMAL);
            expect(negMin.getNextUp()).toBe(IEEEBinary64.NEG_ZERO);
            const posMin = IEEEBinary64.fromEncoding(IEEEBinary64.MIN_SUBNORMAL);
            expect(posMin.getNextDown()).toBe(IEEEBinary64.POS_ZERO);
            // The subnormal/normal boundary is crossed by a single step.
            expect(IEEEBinary64.fromEncoding(IEEEBinary64.MAX_SUBNORMAL).getNextUp())
                .toBe(IEEEBinary64.MIN_NORMAL);
            expect(IEEEBinary64.fromEncoding(IEEEBinary64.MIN_NORMAL).getNextDown())
                .toBe(IEEEBinary64.MAX_SUBNORMAL);
        });

    it('binary32 neighbours agree with binary64 neighbours of the same float value', () => {
        check(fc.integer({ min: 1, max: 0x7F7FFFFE }), e => {
            const b32 = IEEEBinary32.fromEncoding(e);
            const up32 = value32(b32.getNextUp());
            // Stepping in binary32 must skip over many binary64 neighbours,
            // but the ordering is the same and no float lies between.
            return up32 > b32.number
                && Math.fround((b32.number / 2 + up32 / 2)) === b32.number
                || Math.fround((b32.number / 2 + up32 / 2)) === up32;
        });
    });
});
