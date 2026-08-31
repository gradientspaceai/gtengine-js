// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BitHacks.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The leadingBit table in getLeadingBit and the trailingBit table in
// getTrailingBit are based on De Bruijn sequences. The leadingBit table is
// taken from
// https://stackoverflow.com/questions/17027878/algorithm-to-find-the-most-significant-bit
// The trailingBit table is taken from
// https://www.dotnetperls.com/trailing-bits
//
// Port notes:
//   - Upstream provides uint32_t and int32_t overloads of each function. The
//     int32_t overloads immediately static_cast to uint32_t, which is a
//     no-op on the bit pattern, so the port has a single 32-bit function per
//     operation that applies the JavaScript ToUint32 conversion ('>>> 0').
//     Passing a negative 32-bit integer therefore behaves exactly like the
//     C++ int32_t overload (the two's-complement pattern is used).
//   - The 64-bit overloads (uint64_t and int64_t) become the *64 functions
//     that take a bigint. BigInt.asUintN(64, value) reproduces the
//     int64_t -> uint64_t cast for negative inputs.
//   - Upstream guards the int32_t/int64_t overloads with LogAssert when
//     GTE_THROW_ON_BITHACKS_ERROR is defined; that symbol is not defined by
//     default, so the port does not throw either. As upstream documents, a
//     zero input to getLeadingBit/getTrailingBit has no leading or trailing
//     bit and the returned 0 is invalid.
//   - Upstream RoundUpToPowerOfTwo returns uint64_t because the result can be
//     2^32. A JavaScript number represents 2^32 exactly, so the port returns
//     a number.

// De Bruijn table for the leading-bit lookup.
const leadingBitTable: readonly number[] = [
    0, 9, 1, 10, 13, 21, 2, 29,
    11, 14, 16, 18, 22, 25, 3, 30,
    8, 12, 20, 28, 15, 17, 24, 7,
    19, 27, 23, 6, 26, 5, 4, 31
];

// De Bruijn table for the trailing-bit lookup.
const trailingBitTable: readonly number[] = [
    0, 1, 28, 2, 29, 14, 24, 3,
    30, 22, 20, 15, 25, 17, 4, 8,
    31, 27, 13, 23, 21, 19, 16, 7,
    26, 12, 18, 6, 11, 5, 10, 9
];

const MASK32 = 0xFFFFFFFFn;

export class BitHacks {
    static isPowerOfTwo(value: number): boolean {
        const v = value >>> 0;
        return v > 0 && ((v & (v - 1)) >>> 0) === 0;
    }

    static log2OfPowerOfTwo(powerOfTwo: number): number {
        const v = powerOfTwo >>> 0;
        let log2 = (v & 0xAAAAAAAA) !== 0 ? 1 : 0;
        log2 |= ((v & 0xFFFF0000) !== 0 ? 1 : 0) << 4;
        log2 |= ((v & 0xFF00FF00) !== 0 ? 1 : 0) << 3;
        log2 |= ((v & 0xF0F0F0F0) !== 0 ? 1 : 0) << 2;
        log2 |= ((v & 0xCCCCCCCC) !== 0 ? 1 : 0) << 1;
        return log2;
    }

    // The return value of the function is the index into the 32-bit value.
    // For example, getLeadingBit(10) = 3 and getTrailingBit(10) = 2. The
    // value in binary is 0b1010. The bit locations start at 0 on the right of
    // the pattern and end at 31 on the left of the pattern. If the input
    // value is zero, there is no leading bit and no trailing bit. However,
    // the functions return 0, which is considered invalid. Try to call these
    // functions only for positive inputs.
    static getLeadingBit(value: number): number {
        let v = value >>> 0;
        v |= v >>> 1;
        v |= v >>> 2;
        v |= v >>> 4;
        v |= v >>> 8;
        v |= v >>> 16;
        // The C++ multiplication is modulo 2^32; Math.imul provides the low
        // 32 bits of the product and '>>> 27' reinterprets them as unsigned.
        const key = Math.imul(v, 0x07C4ACDD) >>> 27;
        return leadingBitTable[key];
    }

    static getLeadingBit64(value: bigint): number {
        const u = BigInt.asUintN(64, value);
        const v1 = Number((u >> 32n) & MASK32);
        if (v1 !== 0) {
            return BitHacks.getLeadingBit(v1) + 32;
        }

        const v0 = Number(u & MASK32);
        return BitHacks.getLeadingBit(v0);
    }

    static getTrailingBit(value: number): number {
        const v = value >>> 0;
        // Isolate the lowest set bit. C++ uses (value & (~value + 1u)) to
        // avoid a warning about negating an unsigned integer; the JavaScript
        // bitwise operators convert to a signed 32-bit integer first, so
        // (v & -v) is the same two's-complement idiom.
        const key = Math.imul(v & -v, 0x077CB531) >>> 27;
        return trailingBitTable[key];
    }

    static getTrailingBit64(value: bigint): number {
        const u = BigInt.asUintN(64, value);
        const v0 = Number(u & MASK32);
        if (v0 !== 0) {
            return BitHacks.getTrailingBit(v0);
        }

        const v1 = Number((u >> 32n) & MASK32);
        if (v1 !== 0) {
            return BitHacks.getTrailingBit(v1) + 32;
        }
        return 0;
    }

    // Round up to a power of two. If input is zero, the return is 1. If input
    // is larger than 2^{31}, the return is 2^{32}.
    static roundUpToPowerOfTwo(value: number): number {
        const v = value >>> 0;
        if (v > 0) {
            const leading = BitHacks.getLeadingBit(v);
            const mask = (1 << leading) >>> 0;
            if (((v & ~mask) >>> 0) === 0) {
                // value is a power of two
                return v;
            } else {
                // round up to a power of two
                return mask * 2;
            }
        } else {
            return 1;
        }
    }

    // Round down to a power of two. If input is zero, the return is 0.
    static roundDownToPowerOfTwo(value: number): number {
        const v = value >>> 0;
        if (v > 0) {
            const leading = BitHacks.getLeadingBit(v);
            return (1 << leading) >>> 0;
        } else {
            return 0;
        }
    }
}
