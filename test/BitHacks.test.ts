import { describe, it, expect } from 'vitest';
import { BitHacks } from '../src/BitHacks.js';

// Independent reference implementations used for cross-checks.
function refLeadingBit(value: number): number {
    const v = value >>> 0;
    if (v === 0) {
        return 0;
    }
    return 31 - Math.clz32(v);
}

function refTrailingBit(value: number): number {
    const v = value >>> 0;
    if (v === 0) {
        return 0;
    }
    let i = 0;
    while (((v >>> i) & 1) === 0) {
        ++i;
    }
    return i;
}

describe('BitHacks.isPowerOfTwo', () => {
    it('matches the known values', () => {
        expect(BitHacks.isPowerOfTwo(0)).toBe(false);
        expect(BitHacks.isPowerOfTwo(1)).toBe(true);
        expect(BitHacks.isPowerOfTwo(2)).toBe(true);
        expect(BitHacks.isPowerOfTwo(3)).toBe(false);
        expect(BitHacks.isPowerOfTwo(4)).toBe(true);
        expect(BitHacks.isPowerOfTwo(6)).toBe(false);
        expect(BitHacks.isPowerOfTwo(0x80000000)).toBe(true);
        expect(BitHacks.isPowerOfTwo(0xFFFFFFFF)).toBe(false);
    });

    it('agrees with a popcount test over every power of two', () => {
        for (let e = 0; e < 32; ++e) {
            const value = (1 << e) >>> 0;
            expect(BitHacks.isPowerOfTwo(value)).toBe(true);
            if (e > 1) {
                expect(BitHacks.isPowerOfTwo(value + 1)).toBe(false);
            }
        }
    });

    it('treats a negative int32 as its two-complement bit pattern', () => {
        // -2147483648 is 0x80000000, a power of two as an unsigned value.
        expect(BitHacks.isPowerOfTwo(-2147483648)).toBe(true);
        expect(BitHacks.isPowerOfTwo(-1)).toBe(false);
    });
});

describe('BitHacks.log2OfPowerOfTwo', () => {
    it('returns the exponent of every power of two', () => {
        for (let e = 0; e < 32; ++e) {
            expect(BitHacks.log2OfPowerOfTwo((1 << e) >>> 0)).toBe(e);
        }
    });

    it('returns the documented value zero for input zero', () => {
        expect(BitHacks.log2OfPowerOfTwo(0)).toBe(0);
    });
});

describe('BitHacks.getLeadingBit / getTrailingBit', () => {
    it('matches the upstream documentation example', () => {
        // 10 is binary 1010, so the leading bit is at index 3 and the
        // trailing bit is at index 1. NOTE: the upstream comment in
        // BitHacks.h claims GetTrailingBit(10) = 2, which is wrong; the
        // code (here and upstream) returns 1.
        expect(BitHacks.getLeadingBit(10)).toBe(3);
        expect(BitHacks.getTrailingBit(10)).toBe(1);
    });

    it('is exact for the isolated bits', () => {
        for (let e = 0; e < 32; ++e) {
            const value = (1 << e) >>> 0;
            expect(BitHacks.getLeadingBit(value)).toBe(e);
            expect(BitHacks.getTrailingBit(value)).toBe(e);
        }
    });

    it('returns the invalid value 0 for input 0', () => {
        expect(BitHacks.getLeadingBit(0)).toBe(0);
        expect(BitHacks.getTrailingBit(0)).toBe(0);
    });

    it('cross-checks against clz32 and a bit scan for random inputs', () => {
        let state = 123456789;
        for (let i = 0; i < 5000; ++i) {
            // xorshift32 for reproducible pseudorandom 32-bit values.
            state ^= state << 13;
            state ^= state >>> 17;
            state ^= state << 5;
            const value = state >>> 0;
            expect(BitHacks.getLeadingBit(value)).toBe(refLeadingBit(value));
            expect(BitHacks.getTrailingBit(value)).toBe(refTrailingBit(value));
        }
    });

    it('brackets the value between the powers of the leading bit', () => {
        for (let value = 1; value < 2000; ++value) {
            const leading = BitHacks.getLeadingBit(value);
            expect(2 ** leading).toBeLessThanOrEqual(value);
            expect(value).toBeLessThan(2 ** (leading + 1));
        }
    });
});

describe('BitHacks 64-bit variants', () => {
    it('is exact for the isolated bits', () => {
        for (let e = 0; e < 64; ++e) {
            const value = 1n << BigInt(e);
            expect(BitHacks.getLeadingBit64(value)).toBe(e);
            expect(BitHacks.getTrailingBit64(value)).toBe(e);
        }
    });

    it('splits at the 32-bit boundary as upstream does', () => {
        // The low word is used for the trailing bit, the high word for the
        // leading bit.
        const value = (1n << 40n) | (1n << 5n);
        expect(BitHacks.getLeadingBit64(value)).toBe(40);
        expect(BitHacks.getTrailingBit64(value)).toBe(5);
    });

    it('returns 0 for input 0', () => {
        expect(BitHacks.getLeadingBit64(0n)).toBe(0);
        expect(BitHacks.getTrailingBit64(0n)).toBe(0);
    });

    it('reproduces the int64 -> uint64 cast for negative inputs', () => {
        // -1 is 0xFFFFFFFFFFFFFFFF.
        expect(BitHacks.getLeadingBit64(-1n)).toBe(63);
        expect(BitHacks.getTrailingBit64(-1n)).toBe(0);
        // The minimum int64 is 0x8000000000000000.
        expect(BitHacks.getLeadingBit64(-(2n ** 63n))).toBe(63);
        expect(BitHacks.getTrailingBit64(-(2n ** 63n))).toBe(63);
    });

    it('cross-checks against a bigint bit scan for random inputs', () => {
        let state = 987654321n;
        const mask = (1n << 64n) - 1n;
        for (let i = 0; i < 500; ++i) {
            state = (state * 6364136223846793005n + 1442695040888963407n) & mask;
            if (state === 0n) {
                continue;
            }
            let leading = 0;
            for (let b = 63; b >= 0; --b) {
                if ((state >> BigInt(b)) & 1n) {
                    leading = b;
                    break;
                }
            }
            let trailing = 0;
            for (let b = 0; b < 64; ++b) {
                if ((state >> BigInt(b)) & 1n) {
                    trailing = b;
                    break;
                }
            }
            expect(BitHacks.getLeadingBit64(state)).toBe(leading);
            expect(BitHacks.getTrailingBit64(state)).toBe(trailing);
        }
    });
});

describe('BitHacks rounding to powers of two', () => {
    it('matches the documented boundary behavior', () => {
        expect(BitHacks.roundUpToPowerOfTwo(0)).toBe(1);
        expect(BitHacks.roundDownToPowerOfTwo(0)).toBe(0);
        // Input larger than 2^31 rounds up to 2^32.
        expect(BitHacks.roundUpToPowerOfTwo(0xFFFFFFFF)).toBe(4294967296);
        expect(BitHacks.roundDownToPowerOfTwo(0xFFFFFFFF)).toBe(2147483648);
    });

    it('is the identity on powers of two', () => {
        for (let e = 0; e < 32; ++e) {
            const value = (1 << e) >>> 0;
            expect(BitHacks.roundUpToPowerOfTwo(value)).toBe(value);
            expect(BitHacks.roundDownToPowerOfTwo(value)).toBe(value);
        }
    });

    it('brackets non-powers of two', () => {
        for (let value = 1; value < 3000; ++value) {
            const down = BitHacks.roundDownToPowerOfTwo(value);
            const up = BitHacks.roundUpToPowerOfTwo(value);
            expect(down).toBeLessThanOrEqual(value);
            expect(value).toBeLessThanOrEqual(up);
            expect(BitHacks.isPowerOfTwo(down)).toBe(true);
            expect(BitHacks.isPowerOfTwo(up)).toBe(true);
            if (BitHacks.isPowerOfTwo(value)) {
                expect(down).toBe(up);
            } else {
                expect(up).toBe(2 * down);
            }
        }
    });

    it('matches Math.ceil/Math.floor of the base-2 logarithm', () => {
        for (const value of [3, 5, 17, 63, 64, 65, 1000, 100000, 0x7FFFFFFF]) {
            expect(BitHacks.roundDownToPowerOfTwo(value)).toBe(2 ** Math.floor(Math.log2(value)));
            expect(BitHacks.roundUpToPowerOfTwo(value)).toBe(2 ** Math.ceil(Math.log2(value)));
        }
    });
});
