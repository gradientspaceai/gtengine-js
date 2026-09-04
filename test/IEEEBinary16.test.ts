import { describe, it, expect } from 'vitest';
import { IEEEBinary16 } from '../src/IEEEBinary16.js';
import { IEEEClassification } from '../src/IEEEBinary.js';
import { check, fc, finite } from './helpers/arbitraries.js';

// An independent reference implementation of binary16, written from the IEEE
// 754 definition rather than from the upstream bit manipulations.
//
// decodeHalf builds the value of an encoding from the (sign, biased,
// trailing) fields directly. encodeHalf finds the nearest representable
// value by a binary search over the sorted table of nonnegative binary16
// values, breaking ties to even, with a virtual entry 2^16 standing for the
// value that would follow MAX_NORMAL (so that inputs at or beyond the
// midpoint between MAX_NORMAL and 2^16 round to infinity).

function decodeHalf(encoding: number): number {
    const sign = (encoding & 0x8000) !== 0 ? -1 : 1;
    const biased = (encoding & 0x7C00) >>> 10;
    const trailing = encoding & 0x03FF;
    if (biased === 0) {
        return sign * trailing * Math.pow(2, -24);
    }
    if (biased === 31) {
        return trailing === 0 ? sign * Infinity : NaN;
    }
    return sign * (1024 + trailing) * Math.pow(2, biased - 25);
}

// values[k] is the value of encoding k for 0 <= k <= 0x7BFF; the entry at
// 0x7C00 is the "one past MAX_NORMAL" value 2^16.
const refValues: number[] = [];
for (let e = 0; e <= 0x7BFF; ++e) {
    refValues.push(decodeHalf(e));
}
refValues.push(Math.pow(2, 16));

function encodeHalf(value: number): number {
    if (Number.isNaN(value)) {
        // The reference is used only for finite inputs and infinities.
        throw new Error('reference encodeHalf does not handle NaN');
    }
    const signBit = (value < 0 || Object.is(value, -0)) ? 0x8000 : 0x0000;
    const x = Math.abs(value);
    if (x === Infinity) {
        return signBit | 0x7C00;
    }

    // Binary search for the largest index with refValues[i] <= x.
    let lo = 0;
    let hi = refValues.length - 1;
    if (x >= refValues[hi]) {
        return signBit | 0x7C00;
    }
    while (hi - lo > 1) {
        const mid = (lo + hi) >> 1;
        if (refValues[mid] <= x) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    const dLo = x - refValues[lo];
    const dHi = refValues[hi] - x;
    let index: number;
    if (dLo < dHi) {
        index = lo;
    } else if (dHi < dLo) {
        index = hi;
    } else {
        // Ties to even: the encoding with an even low-order bit wins. The
        // virtual entry 0x7C00 has low-order bit 0, so a tie against
        // MAX_NORMAL (low-order bit 1) rounds to infinity.
        index = (lo & 1) === 0 ? lo : hi;
    }
    return signBit | index;
}

// Round a double to the nearest binary32 value, as the C++ conversion to
// float does before the binary16 conversion.
const f32 = Math.fround;

describe('IEEEBinary16 constants', () => {
    it('matches the binary16 parameters of the upstream base class', () => {
        expect(IEEEBinary16.NUM_ENCODING_BITS).toBe(16);
        expect(IEEEBinary16.NUM_EXPONENT_BITS).toBe(5);
        expect(IEEEBinary16.NUM_SIGNIFICAND_BITS).toBe(11);
        expect(IEEEBinary16.NUM_TRAILING_BITS).toBe(10);
        expect(IEEEBinary16.EXPONENT_BIAS).toBe(15);
        expect(IEEEBinary16.MAX_BIASED_EXPONENT).toBe(31);
        expect(IEEEBinary16.MIN_SUB_EXPONENT).toBe(-14);
        expect(IEEEBinary16.MIN_EXPONENT).toBe(-24);
        expect(IEEEBinary16.SIGN_SHIFT).toBe(15);
        expect(IEEEBinary16.SIGN_MASK).toBe(0x8000);
        expect(IEEEBinary16.EXPONENT_MASK).toBe(0x7C00);
        expect(IEEEBinary16.TRAILING_MASK).toBe(0x03FF);
        expect(IEEEBinary16.MIN_SUBNORMAL).toBe(0x0001);
        expect(IEEEBinary16.MAX_SUBNORMAL).toBe(0x03FF);
        expect(IEEEBinary16.MIN_NORMAL).toBe(0x0400);
        expect(IEEEBinary16.MAX_NORMAL).toBe(0x7BFF);
        expect(IEEEBinary16.POS_INFINITY).toBe(0x7C00);
        expect(IEEEBinary16.NEG_INFINITY).toBe(0xFC00);
    });
});

describe('IEEEBinary16 known bit patterns', () => {
    const cases: [number, number][] = [
        [0.0, 0x0000],
        [1.0, 0x3C00],
        [-1.0, 0xBC00],
        [2.0, 0x4000],
        [-2.0, 0xC000],
        [0.5, 0x3800],
        [-0.5, 0xB800],
        [3.0, 0x4200],
        [65504.0, 0x7BFF],      // MAX_NORMAL
        [-65504.0, 0xFBFF],
        [Math.pow(2, -14), 0x0400],   // MIN_NORMAL
        [Math.pow(2, -24), 0x0001],   // MIN_SUBNORMAL
        [Math.pow(2, -14) - Math.pow(2, -24), 0x03FF],  // MAX_SUBNORMAL
        [Infinity, 0x7C00],
        [-Infinity, 0xFC00]
    ];

    for (const [value, encoding] of cases) {
        it(`encodes ${value} as 0x${encoding.toString(16).toUpperCase()}`, () => {
            expect(IEEEBinary16.fromNumber(value).encoding).toBe(encoding);
            expect(IEEEBinary16.fromEncoding(encoding).number).toBe(value);
        });
    }

    it('encodes negative zero with the sign bit set', () => {
        const h = IEEEBinary16.fromNumber(-0);
        expect(h.encoding).toBe(0x8000);
        expect(Object.is(h.number, -0)).toBe(true);
        expect(h.isZero()).toBe(true);
        expect(h.isSignMinus()).toBe(true);
    });

    it('encodes 0.1 as the nearest binary16 value 0x2E66', () => {
        const h = IEEEBinary16.fromNumber(0.1);
        expect(h.encoding).toBe(0x2E66);
        // 0x2E66 = 0.0999755859375
        expect(h.number).toBeCloseTo(0.1, 3);
        expect(h.number).toBe(decodeHalf(0x2E66));
    });
});

describe('IEEEBinary16 special values', () => {
    it('classifies zeros, subnormals, normals, infinities and NaNs', () => {
        expect(IEEEBinary16.fromEncoding(0x0000).getClassification())
            .toBe(IEEEClassification.POS_ZERO);
        expect(IEEEBinary16.fromEncoding(0x8000).getClassification())
            .toBe(IEEEClassification.NEG_ZERO);
        expect(IEEEBinary16.fromEncoding(0x0001).getClassification())
            .toBe(IEEEClassification.POS_SUBNORMAL);
        expect(IEEEBinary16.fromEncoding(0x8001).getClassification())
            .toBe(IEEEClassification.NEG_SUBNORMAL);
        expect(IEEEBinary16.fromEncoding(0x3C00).getClassification())
            .toBe(IEEEClassification.POS_NORMAL);
        expect(IEEEBinary16.fromEncoding(0xBC00).getClassification())
            .toBe(IEEEClassification.NEG_NORMAL);
        expect(IEEEBinary16.fromEncoding(0x7C00).getClassification())
            .toBe(IEEEClassification.POS_INFINITY);
        expect(IEEEBinary16.fromEncoding(0xFC00).getClassification())
            .toBe(IEEEClassification.NEG_INFINITY);
        // The quiet bit is 0x0200.
        expect(IEEEBinary16.fromEncoding(0x7E00).getClassification())
            .toBe(IEEEClassification.QUIET_NAN);
        expect(IEEEBinary16.fromEncoding(0x7C01).getClassification())
            .toBe(IEEEClassification.SIGNALING_NAN);
    });

    it('answers the predicate queries consistently with the classification', () => {
        const bad: number[] = [];
        for (let e = 0; e <= 0xFFFF; ++e) {
            const h = IEEEBinary16.fromEncoding(e);
            const biased = h.getBiased();
            const trailing = h.getTrailing();
            if (h.isZero() !== (biased === 0 && trailing === 0)
                || h.isSubnormal() !== (biased === 0 && trailing > 0)
                || h.isNormal() !== (biased > 0 && biased < 31)
                || h.isFinite() !== (biased < 31)
                || h.isInfinite() !== (biased === 31 && trailing === 0)
                || h.isNaN() !== (biased === 31 && trailing !== 0)
                || h.isSignMinus() !== ((e & 0x8000) !== 0)) {
                bad.push(e);
            }
        }
        expect(bad.slice(0, 5)).toEqual([]);
    });

    it('converts NaNs to NaN values and preserves the payload round-trip', () => {
        for (const e of [0x7C01, 0x7E00, 0x7FFF, 0xFC01, 0xFE00, 0xFFFF]) {
            const h = IEEEBinary16.fromEncoding(e);
            expect(Number.isNaN(h.number)).toBe(true);
            const bits32 = IEEEBinary16.convert16To32(e);
            expect(IEEEBinary16.convert32To16(bits32)).toBe(e);
        }
    });

    it('preserves the quiet-NaN bit through the 16 -> 32 conversion', () => {
        // The 16-quiet bit is 0x0200 and it becomes the 32-quiet bit
        // 0x00400000 after the 13-bit shift.
        const quiet32 = IEEEBinary16.convert16To32(0x7E00);
        expect((quiet32 & 0x00400000) >>> 0).toBe(0x00400000);
        const signaling32 = IEEEBinary16.convert16To32(0x7C01);
        expect((signaling32 & 0x00400000) >>> 0).toBe(0);
    });
});

describe('IEEEBinary16 encode/decode fields', () => {
    it('splits and reassembles the sign, biased exponent and trailing bits', () => {
        const bad: number[] = [];
        for (let e = 0; e <= 0xFFFF; ++e) {
            const parts = IEEEBinary16.fromEncoding(e).getEncoding();
            if (parts.sign !== ((e & 0x8000) >>> 15)
                || parts.biased !== ((e & 0x7C00) >>> 10)
                || parts.trailing !== (e & 0x03FF)
                || IEEEBinary16.fromParts(parts.sign, parts.biased, parts.trailing)
                    .encoding !== e) {
                bad.push(e);
            }
        }
        expect(bad.slice(0, 5)).toEqual([]);
    });
});

describe('IEEEBinary16 conversions', () => {
    it('round-trips every 16-bit encoding through binary32', () => {
        const bad: number[] = [];
        for (let e = 0; e <= 0xFFFF; ++e) {
            const bits32 = IEEEBinary16.convert16To32(e);
            if (IEEEBinary16.convert32To16(bits32) !== e) {
                bad.push(e);
            }
        }
        expect(bad.slice(0, 5)).toEqual([]);
    });

    it('decodes every finite encoding to the reference value', () => {
        const bad: number[] = [];
        for (let e = 0; e <= 0xFFFF; ++e) {
            const value = IEEEBinary16.fromEncoding(e).number;
            const expected = decodeHalf(e);
            if (Number.isNaN(expected)) {
                if (!Number.isNaN(value)) {
                    bad.push(e);
                }
            } else if (value !== expected) {
                bad.push(e);
            }
        }
        expect(bad.slice(0, 5)).toEqual([]);
    });

    it('rounds every representable binary16 value back to itself', () => {
        const bad: number[] = [];
        for (let e = 0; e <= 0xFFFF; ++e) {
            const biased = (e & 0x7C00) >>> 10;
            if (biased === 31 && (e & 0x03FF) !== 0) {
                continue;  // NaN, tested separately
            }
            if (IEEEBinary16.fromNumber(decodeHalf(e)).encoding !== e) {
                bad.push(e);
            }
        }
        expect(bad.slice(0, 5)).toEqual([]);
    });

    it('rounds midpoints between adjacent binary16 values to even', () => {
        // Midpoints in the normal range around 1.
        // 0x3C00 = 1, 0x3C01 = 1 + 2^-10, midpoint rounds to the even 0x3C00.
        const mid0 = f32((decodeHalf(0x3C00) + decodeHalf(0x3C01)) / 2);
        expect(IEEEBinary16.fromNumber(mid0).encoding).toBe(0x3C00);
        // 0x3C01 (odd) and 0x3C02 (even): the midpoint rounds up.
        const mid1 = f32((decodeHalf(0x3C01) + decodeHalf(0x3C02)) / 2);
        expect(IEEEBinary16.fromNumber(mid1).encoding).toBe(0x3C02);
        // Subnormal midpoints: 0x0002 (even) and 0x0003 (odd).
        const mid2 = f32((decodeHalf(0x0002) + decodeHalf(0x0003)) / 2);
        expect(IEEEBinary16.fromNumber(mid2).encoding).toBe(0x0002);
        const mid3 = f32((decodeHalf(0x0003) + decodeHalf(0x0004)) / 2);
        expect(IEEEBinary16.fromNumber(mid3).encoding).toBe(0x0004);
        // Zero and MIN_SUBNORMAL: the midpoint 2^-25 rounds to zero (even).
        expect(IEEEBinary16.fromNumber(Math.pow(2, -25)).encoding).toBe(0x0000);
        // Just above 2^-25 rounds to MIN_SUBNORMAL.
        expect(IEEEBinary16.fromNumber(f32(Math.pow(2, -25) * 1.001)).encoding).toBe(0x0001);
        // MAX_SUBNORMAL and MIN_NORMAL are adjacent; their midpoint rounds
        // to the even MIN_NORMAL (trailing 0).
        const mid4 = f32((decodeHalf(0x03FF) + decodeHalf(0x0400)) / 2);
        expect(IEEEBinary16.fromNumber(mid4).encoding).toBe(0x0400);
    });

    it('overflows to infinity at and beyond the MAX_NORMAL/2^16 midpoint', () => {
        const maxNormal = 65504;
        const midpoint = f32((maxNormal + 65536) / 2);  // 65520
        expect(midpoint).toBe(65520);
        expect(IEEEBinary16.fromNumber(midpoint).encoding).toBe(0x7C00);
        expect(IEEEBinary16.fromNumber(-midpoint).encoding).toBe(0xFC00);
        // Just below the midpoint rounds to MAX_NORMAL.
        const below = f32(65519);
        expect(IEEEBinary16.fromNumber(below).encoding).toBe(0x7BFF);
        expect(IEEEBinary16.fromNumber(1e30).encoding).toBe(0x7C00);
        expect(IEEEBinary16.fromNumber(-1e30).encoding).toBe(0xFC00);
    });

    it('underflows tiny magnitudes to signed zero', () => {
        expect(IEEEBinary16.fromNumber(1e-30).encoding).toBe(0x0000);
        expect(IEEEBinary16.fromNumber(-1e-30).encoding).toBe(0x8000);
        // A binary32 subnormal converts to zero.
        expect(IEEEBinary16.fromNumber(f32(Math.pow(2, -140))).encoding).toBe(0x0000);
    });

    it('matches the reference encoder on a dense sweep of binary32 values', () => {
        // Sweep the 32-bit encodings that cover the whole binary16 range,
        // including the subnormal and overflow boundaries. The mismatches
        // are collected and asserted once so the loop stays fast.
        const scratch = new ArrayBuffer(4);
        const u = new Uint32Array(scratch);
        const v = new Float32Array(scratch);
        const mismatches: string[] = [];
        let tested = 0;
        for (let bits = 0x33000000; bits <= 0x47800000; bits += 0x00000401) {
            u[0] = bits;
            const value = v[0];
            ++tested;
            const got = IEEEBinary16.fromNumber(value).encoding;
            const want = encodeHalf(value);
            if (got !== want) {
                mismatches.push(`+${value}: got 0x${got.toString(16)}, want 0x${want.toString(16)}`);
            }
            const gotNeg = IEEEBinary16.fromNumber(-value).encoding;
            const wantNeg = encodeHalf(-value);
            if (gotNeg !== wantNeg) {
                mismatches.push(`-${value}: got 0x${gotNeg.toString(16)}, want 0x${wantNeg.toString(16)}`);
            }
        }
        expect(tested).toBeGreaterThan(50000);
        expect(mismatches.slice(0, 5)).toEqual([]);
    });

    it('matches the reference encoder on pseudo-random values', () => {
        let seed = 987654321;
        const next = () => {
            seed = (seed * 1103515245 + 12345) & 0x7FFFFFFF;
            return seed / 0x7FFFFFFF;
        };
        const mismatches: string[] = [];
        for (let i = 0; i < 5000; ++i) {
            // Cover magnitudes from far below the subnormal range to far
            // above the overflow threshold.
            const exponent = -30 + next() * 50;
            const value = f32((next() * 2 - 1) * Math.pow(2, exponent));
            const got = IEEEBinary16.fromNumber(value).encoding;
            const want = encodeHalf(value);
            if (got !== want) {
                mismatches.push(`${value}: got 0x${got.toString(16)}, want 0x${want.toString(16)}`);
            }
        }
        expect(mismatches.slice(0, 5)).toEqual([]);
    });
});

describe('IEEEBinary16 neighbors', () => {
    it('gives the next-up and next-down for the special encodings', () => {
        expect(IEEEBinary16.fromEncoding(0x0000).getNextUp()).toBe(0x0001);
        expect(IEEEBinary16.fromEncoding(0x8000).getNextUp()).toBe(0x0001);
        expect(IEEEBinary16.fromEncoding(0x8001).getNextUp()).toBe(0x8000);
        expect(IEEEBinary16.fromEncoding(0x03FF).getNextUp()).toBe(0x0400);
        expect(IEEEBinary16.fromEncoding(0x7BFF).getNextUp()).toBe(0x7C00);
        expect(IEEEBinary16.fromEncoding(0x7C00).getNextUp()).toBe(0x7C00);
        expect(IEEEBinary16.fromEncoding(0xFC00).getNextUp()).toBe(0xFBFF);

        expect(IEEEBinary16.fromEncoding(0x0000).getNextDown()).toBe(0x8001);
        expect(IEEEBinary16.fromEncoding(0x8000).getNextDown()).toBe(0x8001);
        expect(IEEEBinary16.fromEncoding(0x0001).getNextDown()).toBe(0x0000);
        expect(IEEEBinary16.fromEncoding(0x0400).getNextDown()).toBe(0x03FF);
        expect(IEEEBinary16.fromEncoding(0x7C00).getNextDown()).toBe(0x7BFF);
        expect(IEEEBinary16.fromEncoding(0xFC00).getNextDown()).toBe(0xFC00);

        // NaNs are their own neighbors.
        expect(IEEEBinary16.fromEncoding(0x7E00).getNextUp()).toBe(0x7E00);
        expect(IEEEBinary16.fromEncoding(0x7E00).getNextDown()).toBe(0x7E00);
    });

    it('makes next-up and next-down inverse on the finite encodings', () => {
        const bad: number[] = [];
        for (let e = 0; e <= 0xFFFF; ++e) {
            const h = IEEEBinary16.fromEncoding(e);
            if (!h.isFinite()) {
                continue;
            }
            const back = IEEEBinary16.fromEncoding(h.getNextUp()).getNextDown();
            // The zeros collapse: next-up of -0 and +0 are both
            // MIN_SUBNORMAL, whose next-down is +0.
            const want = h.isZero() ? 0x0000 : e;
            if (back !== want) {
                bad.push(e);
            }
        }
        expect(bad.slice(0, 5)).toEqual([]);
    });

    it('orders next-up strictly above the number', () => {
        const bad: number[] = [];
        for (let e = 0; e <= 0xFFFF; ++e) {
            const h = IEEEBinary16.fromEncoding(e);
            if (!h.isFinite()) {
                continue;
            }
            if (!(IEEEBinary16.fromEncoding(h.getNextUp()).number > h.number)) {
                bad.push(e);
            }
        }
        expect(bad.slice(0, 5)).toEqual([]);
    });
});

describe('IEEEBinary16 arithmetic and comparison', () => {
    it('negates by flipping the sign bit', () => {
        expect(IEEEBinary16.negate(IEEEBinary16.fromNumber(2)).encoding).toBe(0xC000);
        expect(IEEEBinary16.negate(IEEEBinary16.fromEncoding(0x0000)).encoding).toBe(0x8000);
        expect(IEEEBinary16.negate(IEEEBinary16.fromEncoding(0x7C00)).encoding).toBe(0xFC00);
    });

    it('computes the arithmetic operations in binary32', () => {
        const a = IEEEBinary16.fromNumber(1.5);
        const b = IEEEBinary16.fromNumber(0.25);
        expect(IEEEBinary16.add(a, b)).toBe(1.75);
        expect(IEEEBinary16.sub(a, b)).toBe(1.25);
        expect(IEEEBinary16.mul(a, b)).toBe(0.375);
        expect(IEEEBinary16.div(a, b)).toBe(6);
        // Mixed with a float operand.
        expect(IEEEBinary16.add(a, 0.5)).toBe(2);
        expect(IEEEBinary16.mul(2, b)).toBe(0.5);
    });

    it('keeps the extra precision of the binary32 result', () => {
        // MAX_NORMAL + MAX_NORMAL overflows binary16 but not binary32.
        const m = IEEEBinary16.fromNumber(65504);
        expect(IEEEBinary16.add(m, m)).toBe(131008);
    });

    it('updates in place with the assignment operators', () => {
        const x = IEEEBinary16.fromNumber(1);
        IEEEBinary16.addAssign(x, IEEEBinary16.fromNumber(0.5));
        expect(x.number).toBe(1.5);
        IEEEBinary16.mulAssign(x, 2);
        expect(x.number).toBe(3);
        IEEEBinary16.subAssign(x, 1);
        expect(x.number).toBe(2);
        IEEEBinary16.divAssign(x, 4);
        expect(x.number).toBe(0.5);
    });

    it('compares by the converted values', () => {
        const a = IEEEBinary16.fromNumber(1);
        const b = IEEEBinary16.fromNumber(2);
        const posZero = IEEEBinary16.fromEncoding(0x0000);
        const negZero = IEEEBinary16.fromEncoding(0x8000);

        expect(IEEEBinary16.lessThan(a, b)).toBe(true);
        expect(IEEEBinary16.greaterThan(b, a)).toBe(true);
        expect(IEEEBinary16.lessThanOrEqual(a, a)).toBe(true);
        expect(IEEEBinary16.greaterThanOrEqual(a, a)).toBe(true);
        expect(IEEEBinary16.equals(a, a)).toBe(true);
        expect(IEEEBinary16.notEquals(a, b)).toBe(true);
        // -0 == +0 as floating-point values even though the encodings differ.
        expect(IEEEBinary16.equals(posZero, negZero)).toBe(true);
        expect(posZero.encoding).not.toBe(negZero.encoding);
        // A NaN compares unequal to itself.
        const nan = IEEEBinary16.fromEncoding(0x7E00);
        expect(IEEEBinary16.equals(nan, nan)).toBe(false);
        expect(IEEEBinary16.notEquals(nan, nan)).toBe(true);
    });
});

describe('IEEEBinary16 math wrappers', () => {
    it('evaluates the std:: functions in binary32 and rounds to binary16', () => {
        const one = IEEEBinary16.fromNumber(1);
        expect(IEEEBinary16.sqrt(IEEEBinary16.fromNumber(4)).number).toBe(2);
        expect(IEEEBinary16.fabs(IEEEBinary16.fromNumber(-3)).number).toBe(3);
        expect(IEEEBinary16.floor(IEEEBinary16.fromNumber(2.5)).number).toBe(2);
        expect(IEEEBinary16.ceil(IEEEBinary16.fromNumber(2.5)).number).toBe(3);
        expect(IEEEBinary16.exp2(IEEEBinary16.fromNumber(3)).number).toBe(8);
        expect(IEEEBinary16.log2(IEEEBinary16.fromNumber(8)).number).toBe(3);
        expect(IEEEBinary16.acos(one).number).toBe(0);
        expect(IEEEBinary16.cos(IEEEBinary16.fromNumber(0)).number).toBe(1);
        expect(IEEEBinary16.sin(IEEEBinary16.fromNumber(0)).number).toBe(0);
        expect(IEEEBinary16.pow(IEEEBinary16.fromNumber(2),
            IEEEBinary16.fromNumber(10)).number).toBe(1024);
        expect(IEEEBinary16.fmod(IEEEBinary16.fromNumber(7),
            IEEEBinary16.fromNumber(3)).number).toBe(1);
        expect(IEEEBinary16.ldexp(one, 4).number).toBe(16);
    });

    it('evaluates the gte:: Functions.h overloads', () => {
        expect(IEEEBinary16.sqr(IEEEBinary16.fromNumber(3)).number).toBe(9);
        expect(IEEEBinary16.invsqrt(IEEEBinary16.fromNumber(4)).number).toBe(0.5);
        expect(IEEEBinary16.isign(IEEEBinary16.fromNumber(-7))).toBe(-1);
        expect(IEEEBinary16.isign(IEEEBinary16.fromNumber(0))).toBe(0);
        expect(IEEEBinary16.sign(IEEEBinary16.fromNumber(7)).number).toBe(1);
        expect(IEEEBinary16.saturate(IEEEBinary16.fromNumber(3)).number).toBe(1);
        expect(IEEEBinary16.saturate(IEEEBinary16.fromNumber(-3)).number).toBe(0);
        expect(IEEEBinary16.clamp(IEEEBinary16.fromNumber(5),
            IEEEBinary16.fromNumber(0), IEEEBinary16.fromNumber(2)).number).toBe(2);
        expect(IEEEBinary16.sinpi(IEEEBinary16.fromNumber(0.5)).number).toBe(1);
        expect(IEEEBinary16.cospi(IEEEBinary16.fromNumber(0)).number).toBe(1);
        expect(IEEEBinary16.exp10(IEEEBinary16.fromNumber(2)).number).toBe(100);
    });
});

describe('IEEEBinary16 upstream quirk: NaN with a small payload', () => {
    it('converts a binary32 NaN whose high payload bits are zero to infinity', () => {
        // Convert32To16 keeps only the high-order 9 bits of the 23-bit
        // trailing significand. A binary32 signaling NaN whose payload is
        // below 2^13 therefore maps to a zero 16-trailing significand, which
        // is the encoding of infinity rather than a NaN. The port preserves
        // this upstream behavior.
        const nanBits = 0x7F800001;  // signaling NaN, payload 1
        expect(IEEEBinary16.convert32To16(nanBits)).toBe(0x7C00);
        const negNanBits = 0xFF800001 >>> 0;
        expect(IEEEBinary16.convert32To16(negNanBits)).toBe(0xFC00);
        // A NaN with a large enough payload survives.
        expect(IEEEBinary16.convert32To16(0x7FC00000)).toBe(0x7E00);
    });

    it('converts the canonical JavaScript NaN to a binary16 quiet NaN', () => {
        const h = IEEEBinary16.fromNumber(NaN);
        expect(h.isNaN()).toBe(true);
        expect(h.isQuietNaN()).toBe(true);
    });
});

// ---------------------------------------------------------------------------
// Verification block (V16).
// ---------------------------------------------------------------------------

// Bit-pattern views for the binary32 reference computations.
const refBuffer = new ArrayBuffer(4);
const refF32 = new Float32Array(refBuffer);
const refU32 = new Uint32Array(refBuffer);

function bitsOf32(value: number): number {
    refF32[0] = value;
    return refU32[0] >>> 0;
}

function valueOf32(bits: number): number {
    refU32[0] = bits >>> 0;
    return refF32[0];
}

// Every uint32 bit pattern, drawn uniformly.
const bits32 = fc.integer({ min: -0x80000000, max: 0x7FFFFFFF })
    .map(x => x >>> 0);

// Every uint16 encoding, drawn uniformly.
const encoding16 = fc.integer({ min: 0, max: 0xFFFF });

describe('IEEEBinary16 verification', () => {
    it('decodes every one of the 65536 encodings to the IEEE 754 value', () => {
        for (let e = 0; e <= 0xFFFF; ++e) {
            const x = IEEEBinary16.fromEncoding(e);
            const expected = decodeHalf(e);
            if (Number.isNaN(expected)) {
                expect(Number.isNaN(x.number)).toBe(true);
            } else {
                // Object.is separates -0 from +0, which is what the sign bit
                // of the encoding demands.
                expect(Object.is(x.number, expected)).toBe(true);
            }
            // Every binary16 value is exactly a binary32 value.
            expect(Object.is(Math.fround(x.number), x.number)).toBe(true);
        }
    }, 30000);

    it('re-encodes every non-NaN value back to the same encoding', () => {
        for (let e = 0; e <= 0xFFFF; ++e) {
            const x = IEEEBinary16.fromEncoding(e);
            if (x.isNaN()) { continue; }
            expect(IEEEBinary16.fromNumber(x.number).encoding).toBe(e);
            // Idempotence of the encoding pipeline.
            expect(IEEEBinary16.convert32To16(
                IEEEBinary16.convert16To32(e))).toBe(e);
        }
    }, 30000);

    it('splits and rebuilds the sign, biased exponent and trailing fields', () => {
        for (let e = 0; e <= 0xFFFF; ++e) {
            const x = IEEEBinary16.fromEncoding(e);
            const { sign, biased, trailing } = x.getEncoding();
            expect(sign).toBe(e >>> 15);
            expect(biased).toBe((e >>> 10) & 0x1F);
            expect(trailing).toBe(e & 0x3FF);
            expect(IEEEBinary16.fromParts(sign, biased, trailing).encoding)
                .toBe(e);
        }
    }, 30000);

    it('classifies every encoding consistently with its fields', () => {
        for (let e = 0; e <= 0xFFFF; ++e) {
            const x = IEEEBinary16.fromEncoding(e);
            const sign = e >>> 15;
            const biased = (e >>> 10) & 0x1F;
            const trailing = e & 0x3FF;
            let expected: IEEEClassification;
            if (biased === 0) {
                expected = trailing === 0
                    ? (sign ? IEEEClassification.NEG_ZERO
                        : IEEEClassification.POS_ZERO)
                    : (sign ? IEEEClassification.NEG_SUBNORMAL
                        : IEEEClassification.POS_SUBNORMAL);
            } else if (biased < 31) {
                expected = sign ? IEEEClassification.NEG_NORMAL
                    : IEEEClassification.POS_NORMAL;
            } else if (trailing === 0) {
                expected = sign ? IEEEClassification.NEG_INFINITY
                    : IEEEClassification.POS_INFINITY;
            } else {
                expected = (trailing & 0x200) !== 0
                    ? IEEEClassification.QUIET_NAN
                    : IEEEClassification.SIGNALING_NAN;
            }
            expect(x.getClassification()).toBe(expected);
            expect(x.isZero()).toBe(biased === 0 && trailing === 0);
            expect(x.isSignMinus()).toBe(sign === 1);
            expect(x.isSubnormal()).toBe(biased === 0 && trailing > 0);
            expect(x.isNormal()).toBe(biased > 0 && biased < 31);
            expect(x.isFinite()).toBe(biased < 31);
            expect(x.isInfinite()).toBe(biased === 31 && trailing === 0);
            expect(x.isNaN()).toBe(biased === 31 && trailing !== 0);
        }
    }, 30000);

    it('orders next-up and next-down by the value ordering of the encodings', () => {
        // The nonnegative encodings 0x0000..0x7C00 are increasing in value and
        // the negative ones mirror them, which is what the branchy upstream
        // code implements.
        const order: number[] = [];
        for (let e = 0x7C00; e >= 0x0001; --e) { order.push(0x8000 | e); }
        order.push(0x8000);            // -0
        order.push(0x0000);            // +0
        for (let e = 0x0001; e <= 0x7C00; ++e) { order.push(e); }

        for (let i = 0; i < order.length; ++i) {
            const x = IEEEBinary16.fromEncoding(order[i]);
            // -0 and +0 are adjacent in the list but equal in value; upstream
            // maps both to MIN_SUBNORMAL going up and to -MIN_SUBNORMAL going
            // down, so treat the pair as a single step.
            const isZero = order[i] === 0x0000 || order[i] === 0x8000;
            if (isZero) {
                expect(x.getNextUp()).toBe(IEEEBinary16.MIN_SUBNORMAL);
                expect(x.getNextDown())
                    .toBe(IEEEBinary16.SIGN_MASK | IEEEBinary16.MIN_SUBNORMAL);
                continue;
            }
            const up = i + 1 < order.length ? order[i + 1] : IEEEBinary16.POS_INFINITY;
            const down = i > 0 ? order[i - 1] : IEEEBinary16.NEG_INFINITY;
            expect(x.getNextUp()).toBe(up);
            expect(x.getNextDown()).toBe(down);
        }

        // Every NaN is its own neighbor in both directions.
        for (let t = 1; t <= 0x3FF; ++t) {
            for (const s of [0x0000, 0x8000]) {
                const e = s | 0x7C00 | t;
                expect(IEEEBinary16.fromEncoding(e).getNextUp()).toBe(e);
                expect(IEEEBinary16.fromEncoding(e).getNextDown()).toBe(e);
            }
        }
    }, 30000);

    it('flips only the sign bit when negating, for every encoding', () => {
        for (let e = 0; e <= 0xFFFF; ++e) {
            const n = IEEEBinary16.negate(IEEEBinary16.fromEncoding(e));
            expect(n.encoding).toBe(e ^ 0x8000);
            expect(IEEEBinary16.negate(n).encoding).toBe(e);
        }
    }, 30000);

    it('rounds every binary32 bit pattern the way the reference encoder does', () => {
        check(bits32, bits => {
            const value = valueOf32(bits);
            const encoded = IEEEBinary16.convert32To16(bits);
            if (Number.isNaN(value)) {
                // Upstream keeps only the high-order 9 payload bits, so a NaN
                // whose payload lives in the low 13 bits becomes an infinity
                // (upstream bug, preserved by the port).
                const payload = (bits & 0x007FFFFF) >>> 13;
                expect(encoded).toBe(((bits >>> 16) & 0x8000) | 0x7C00 | payload);
                expect(IEEEBinary16.fromEncoding(encoded).isNaN())
                    .toBe(payload !== 0);
                return;
            }
            expect(encoded).toBe(encodeHalf(value));
        });
    });

    it('rounds every double to the nearest binary16 of its binary32 rounding', () => {
        // The C++ constructor from double first converts to float, so the
        // conversion is a double rounding by construction.
        check(finite(-70000, 70000), value => {
            const x = IEEEBinary16.fromNumber(value);
            expect(x.encoding).toBe(encodeHalf(f32(value)));
        });
    });

    it('evaluates the binary operators in binary32', () => {
        check(fc.tuple(encoding16, encoding16), ([e0, e1]) => {
            const x = IEEEBinary16.fromEncoding(e0);
            const y = IEEEBinary16.fromEncoding(e1);
            const a = x.number;
            const b = y.number;
            // For binary32 operands the binary64 result rounds to the same
            // binary32 value the C++ float arithmetic produces (2*24+1 <= 53).
            const ops: Array<[number, number]> = [
                [IEEEBinary16.add(x, y), f32(a + b)],
                [IEEEBinary16.sub(x, y), f32(a - b)],
                [IEEEBinary16.mul(x, y), f32(a * b)],
                [IEEEBinary16.div(x, y), f32(a / b)]
            ];
            for (const [got, want] of ops) {
                if (Number.isNaN(want)) {
                    expect(Number.isNaN(got)).toBe(true);
                } else {
                    expect(Object.is(got, want)).toBe(true);
                }
                // The result stays a binary32 value, as the C++ float return
                // type requires.
                expect(Object.is(f32(got), got) || Number.isNaN(got)).toBe(true);
            }
            // The mixed IEEEBinary16/number overloads agree.
            expect(Object.is(IEEEBinary16.add(x, b), IEEEBinary16.add(x, y))
                || Number.isNaN(b)).toBe(true);
            expect(Object.is(IEEEBinary16.mul(a, y), IEEEBinary16.mul(x, y))
                || Number.isNaN(a)).toBe(true);
        });
    });

    it('rounds the compound assignments back into binary16', () => {
        check(fc.tuple(encoding16, encoding16), ([e0, e1]) => {
            const y = IEEEBinary16.fromEncoding(e1);
            for (const [assign, binary] of [
                [IEEEBinary16.addAssign, IEEEBinary16.add],
                [IEEEBinary16.subAssign, IEEEBinary16.sub],
                [IEEEBinary16.mulAssign, IEEEBinary16.mul],
                [IEEEBinary16.divAssign, IEEEBinary16.div]
            ] as Array<[typeof IEEEBinary16.addAssign, typeof IEEEBinary16.add]>) {
                const x = IEEEBinary16.fromEncoding(e0);
                const expected = IEEEBinary16.fromNumber(binary(x, y)).encoding;
                const returned = assign(x, y);
                // Upstream mutates the left operand and returns a reference.
                expect(returned).toBe(x);
                expect(x.encoding).toBe(expected);
            }
        });
    });

    it('compares by the converted values, with NaN unordered', () => {
        check(fc.tuple(encoding16, encoding16), ([e0, e1]) => {
            const x = IEEEBinary16.fromEncoding(e0);
            const y = IEEEBinary16.fromEncoding(e1);
            const a = x.number;
            const b = y.number;
            expect(IEEEBinary16.equals(x, y)).toBe(a === b);
            expect(IEEEBinary16.notEquals(x, y)).toBe(a !== b);
            expect(IEEEBinary16.lessThan(x, y)).toBe(a < b);
            expect(IEEEBinary16.lessThanOrEqual(x, y)).toBe(a <= b);
            expect(IEEEBinary16.greaterThan(x, y)).toBe(a > b);
            expect(IEEEBinary16.greaterThanOrEqual(x, y)).toBe(a >= b);
            if (x.isNaN() || y.isNaN()) {
                expect(IEEEBinary16.equals(x, y)).toBe(false);
                expect(IEEEBinary16.lessThan(x, y)).toBe(false);
                expect(IEEEBinary16.greaterThan(x, y)).toBe(false);
            }
            // +0 and -0 compare equal even though the encodings differ.
            if ((e0 & 0x7FFF) === 0 && (e1 & 0x7FFF) === 0) {
                expect(IEEEBinary16.equals(x, y)).toBe(true);
            }
        });
    });

    it('decomposes every finite nonzero encoding exactly with frexp', () => {
        for (let e = 0; e <= 0xFFFF; ++e) {
            const x = IEEEBinary16.fromEncoding(e);
            const { result, exponent } = IEEEBinary16.frexp(x);
            if (x.isZero()) {
                expect(Object.is(result.number, x.number)).toBe(true);
                expect(exponent).toBe(0);
                continue;
            }
            if (!x.isFinite()) {
                expect(exponent).toBe(0);
                if (x.isNaN()) {
                    expect(result.isNaN()).toBe(true);
                } else {
                    expect(Object.is(result.number, x.number)).toBe(true);
                }
                continue;
            }
            const f = result.number;
            expect(Math.abs(f)).toBeGreaterThanOrEqual(0.5);
            expect(Math.abs(f)).toBeLessThan(1);
            // f has at most as many significand bits as x, and the scaling is
            // a power of two, so the reconstruction is exact.
            expect(f * Math.pow(2, exponent)).toBe(x.number);
            expect(exponent).toBeGreaterThanOrEqual(-23);
            expect(exponent).toBeLessThanOrEqual(16);
        }
    }, 30000);

    it('evaluates the math wrappers in binary32 and rounds to binary16', () => {
        check(encoding16, e => {
            const x = IEEEBinary16.fromEncoding(e);
            const v = x.number;
            const cases: Array<[IEEEBinary16, number]> = [
                [IEEEBinary16.fabs(x), Math.abs(v)],
                [IEEEBinary16.ceil(x), Math.ceil(v)],
                [IEEEBinary16.floor(x), Math.floor(v)],
                [IEEEBinary16.sqrt(x), Math.sqrt(v)],
                [IEEEBinary16.sin(x), Math.sin(v)],
                [IEEEBinary16.cos(x), Math.cos(v)],
                [IEEEBinary16.atan(x), Math.atan(v)],
                [IEEEBinary16.exp(x), Math.exp(v)],
                [IEEEBinary16.log(x), Math.log(v)]
            ];
            for (const [got, want] of cases) {
                expect(got.encoding).toBe(
                    IEEEBinary16.fromNumber(f32(want)).encoding);
            }
            // fabs never changes anything but the sign bit.
            expect(IEEEBinary16.fabs(x).encoding).toBe(
                x.isNaN() ? IEEEBinary16.fromNumber(Math.abs(v)).encoding
                    : (e & 0x7FFF));
        });
    });

    it('round-trips the 16 -> 32 -> 16 conversions for every encoding', () => {
        check(encoding16, e => {
            const bits = IEEEBinary16.convert16To32(e);
            expect(bits).toBe(bits >>> 0);
            expect(IEEEBinary16.convert32To16(bits)).toBe(e);
            // The binary32 encoding agrees with the platform conversion for
            // everything but NaN, whose payload the platform may canonicalize.
            const x = IEEEBinary16.fromEncoding(e);
            if (!x.isNaN()) {
                expect(bitsOf32(x.number)).toBe(bits);
            }
        });
    });
});
