import { describe, it, expect } from 'vitest';
import { IEEEBinary16 } from '../src/IEEEBinary16.js';
import { IEEEClassification } from '../src/IEEEBinary.js';

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
