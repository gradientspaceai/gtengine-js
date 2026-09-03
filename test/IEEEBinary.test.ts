import { describe, it, expect } from 'vitest';
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
