// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IEEEBinary.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Access to the bit-level IEEE 754 representations of floating-point
// numbers.
//
// Port notes: upstream is 'template <typename Float, typename UInt, int32_t
// NumBits, int32_t Precision> class IEEEBinary' with a union of the Float
// and UInt members, instantiated as IEEEBinary32 (float/uint32_t) and
// IEEEBinary64 (double/uint64_t). The port provides the two concrete
// classes. IEEEBinary32 stores its encoding as an unsigned 32-bit 'number'
// (maintained with '>>> 0'); IEEEBinary64 stores its encoding as a 'bigint'
// because 64-bit patterns exceed the safe integer range. The C++ union is
// ported as an 'encoding' field plus a 'number' accessor pair that converts
// through a typed-array view. The C++ constructor overloads on Float/UInt
// are ambiguous in TypeScript (both are 'number' for IEEEBinary32), so
// construction uses the static factories fromEncoding, fromNumber and
// fromParts. For IEEEBinary64, getSign/getBiased/getTrailing return 'number'
// (all three fields fit in 53 bits); whole-encoding values are 'bigint'.
// The nested enum Classification is exported as IEEEClassification, shared
// by both classes.

// The types of numbers.
export enum IEEEClassification {
    NEG_INFINITY,
    NEG_SUBNORMAL,
    NEG_NORMAL,
    NEG_ZERO,
    POS_ZERO,
    POS_SUBNORMAL,
    POS_NORMAL,
    POS_INFINITY,
    QUIET_NAN,
    SIGNALING_NAN
}

// Scratch views for converting between floating-point values and their bit
// patterns.
const scratch32Buffer = new ArrayBuffer(4);
const scratch32F = new Float32Array(scratch32Buffer);
const scratch32U = new Uint32Array(scratch32Buffer);
const scratch64Buffer = new ArrayBuffer(8);
const scratch64F = new Float64Array(scratch64Buffer);
const scratch64U = new BigUint64Array(scratch64Buffer);

function classify(sign: number, biased: number, trailing: boolean,
    maxBiased: number, quiet: boolean): IEEEClassification {
    if (biased === 0) {
        if (!trailing) {
            return sign !== 0 ? IEEEClassification.NEG_ZERO : IEEEClassification.POS_ZERO;
        }
        return sign !== 0 ? IEEEClassification.NEG_SUBNORMAL : IEEEClassification.POS_SUBNORMAL;
    }
    if (biased < maxBiased) {
        return sign !== 0 ? IEEEClassification.NEG_NORMAL : IEEEClassification.POS_NORMAL;
    }
    if (!trailing) {
        return sign !== 0 ? IEEEClassification.NEG_INFINITY : IEEEClassification.POS_INFINITY;
    }
    return quiet ? IEEEClassification.QUIET_NAN : IEEEClassification.SIGNALING_NAN;
}

// The bit-level representation of IEEE binary32 (C++ float). The encoding is
// an unsigned 32-bit integer stored in a 'number'.
export class IEEEBinary32 {
    // Special constants.
    static readonly NUM_ENCODING_BITS = 32;
    static readonly NUM_EXPONENT_BITS = 8;
    static readonly NUM_SIGNIFICAND_BITS = 24;
    static readonly NUM_TRAILING_BITS = 23;
    static readonly EXPONENT_BIAS = 127;
    static readonly MAX_BIASED_EXPONENT = 255;
    static readonly MIN_SUB_EXPONENT = -126;
    static readonly MIN_EXPONENT = -149;
    static readonly SIGN_SHIFT = 31;

    static readonly SIGN_MASK = 0x80000000;
    static readonly NOT_SIGN_MASK = 0x7FFFFFFF;
    static readonly TRAILING_MASK = 0x007FFFFF;
    static readonly EXPONENT_MASK = 0x7F800000;
    static readonly NAN_QUIET_MASK = 0x00400000;
    static readonly NAN_PAYLOAD_MASK = 0x003FFFFF;
    static readonly MAX_TRAILING = 0x007FFFFF;
    static readonly SUP_TRAILING = 0x00800000;
    static readonly POS_ZERO = 0x00000000;
    static readonly NEG_ZERO = 0x80000000;
    static readonly MIN_SUBNORMAL = 0x00000001;
    static readonly MAX_SUBNORMAL = 0x007FFFFF;
    static readonly MIN_NORMAL = 0x00800000;
    static readonly MAX_NORMAL = 0x7F7FFFFF;
    static readonly POS_INFINITY = 0x7F800000;
    static readonly NEG_INFINITY = 0xFF800000;

    // The bits of the represented number as an unsigned 32-bit integer.
    encoding: number;

    constructor(encoding: number = 0) {
        this.encoding = encoding >>> 0;
    }

    // Construction from primitive elements (ports of the C++ constructor
    // overloads on UInt, Float and (sign, biased, trailing)).
    static fromEncoding(encoding: number): IEEEBinary32 {
        return new IEEEBinary32(encoding);
    }

    static fromNumber(value: number): IEEEBinary32 {
        const result = new IEEEBinary32();
        result.number = value;
        return result;
    }

    static fromParts(sign: number, biased: number, trailing: number): IEEEBinary32 {
        const result = new IEEEBinary32();
        result.setEncoding(sign, biased, trailing);
        return result;
    }

    // The represented floating-point value (the other member of the C++
    // union). Assigning rounds the input to binary32 precision.
    get number(): number {
        scratch32U[0] = this.encoding;
        return scratch32F[0];
    }

    set number(value: number) {
        scratch32F[0] = value;
        this.encoding = scratch32U[0];
    }

    getClassification(): IEEEClassification {
        const { sign, biased, trailing } = this.getEncoding();
        return classify(sign, biased, trailing !== 0,
            IEEEBinary32.MAX_BIASED_EXPONENT,
            (trailing & IEEEBinary32.NAN_QUIET_MASK) !== 0);
    }

    isZero(): boolean {
        return this.encoding === IEEEBinary32.POS_ZERO
            || this.encoding === IEEEBinary32.NEG_ZERO;
    }

    isSignMinus(): boolean {
        return (this.encoding & IEEEBinary32.SIGN_MASK) !== 0;
    }

    isSubnormal(): boolean {
        return this.getBiased() === 0 && this.getTrailing() > 0;
    }

    isNormal(): boolean {
        const biased = this.getBiased();
        return 0 < biased && biased < IEEEBinary32.MAX_BIASED_EXPONENT;
    }

    isFinite(): boolean {
        return this.getBiased() < IEEEBinary32.MAX_BIASED_EXPONENT;
    }

    isInfinite(): boolean {
        return this.getBiased() === IEEEBinary32.MAX_BIASED_EXPONENT
            && this.getTrailing() === 0;
    }

    isNaN(): boolean {
        return this.getBiased() === IEEEBinary32.MAX_BIASED_EXPONENT
            && this.getTrailing() !== 0;
    }

    isQuietNaN(): boolean {
        const trailing = this.getTrailing();
        return this.getBiased() === IEEEBinary32.MAX_BIASED_EXPONENT
            && (trailing & IEEEBinary32.NAN_QUIET_MASK) !== 0;
    }

    isSignalingNaN(): boolean {
        const trailing = this.getTrailing();
        return this.getBiased() === IEEEBinary32.MAX_BIASED_EXPONENT
            && (trailing & IEEEBinary32.NAN_QUIET_MASK) === 0
            && (trailing & IEEEBinary32.NAN_PAYLOAD_MASK) !== 0;
    }

    // Get neighboring numbers (the encodings, as unsigned 32-bit integers).
    getNextUp(): number {
        const { sign, biased, trailing } = this.getEncoding();

        if (biased === 0) {
            if (trailing === 0) {
                // The next-up for both -0 and +0 is MIN_SUBNORMAL.
                return IEEEBinary32.MIN_SUBNORMAL;
            }
            if (sign !== 0) {
                // When trailing is 1, 'this' is -MIN_SUBNORMAL and next-up
                // is -0.
                return (IEEEBinary32.SIGN_MASK | (trailing - 1)) >>> 0;
            }
            // When trailing is MAX_TRAILING, 'this' is MAX_SUBNORMAL and
            // next-up is MIN_NORMAL.
            return trailing + 1;
        }
        if (biased < IEEEBinary32.MAX_BIASED_EXPONENT) {
            const nonnegative = (this.encoding & IEEEBinary32.NOT_SIGN_MASK);
            if (sign !== 0) {
                return (IEEEBinary32.SIGN_MASK | (nonnegative - 1)) >>> 0;
            }
            return nonnegative + 1;
        }
        if (trailing === 0) {
            if (sign !== 0) {
                // The next-up of -INFINITY is -MAX_NORMAL.
                return (IEEEBinary32.SIGN_MASK | IEEEBinary32.MAX_NORMAL) >>> 0;
            }
            // The next-up of +INFINITY is +INFINITY.
            return IEEEBinary32.POS_INFINITY;
        }
        // The number is a quiet or signaling NaN, possibly with payload.
        // Just return the number itself.
        return this.encoding;
    }

    getNextDown(): number {
        const { sign, biased, trailing } = this.getEncoding();

        if (biased === 0) {
            if (trailing === 0) {
                // The next-down for both -0 and +0 is -MIN_SUBNORMAL.
                return (IEEEBinary32.SIGN_MASK | IEEEBinary32.MIN_SUBNORMAL) >>> 0;
            }
            if (sign === 0) {
                // When trailing is 1, 'this' is MIN_SUBNORMAL and next-down
                // is +0.
                return trailing - 1;
            }
            // When trailing is MAX_TRAILING, 'this' is -MAX_SUBNORMAL and
            // next-down is -MIN_NORMAL.
            return (IEEEBinary32.SIGN_MASK | (trailing + 1)) >>> 0;
        }
        if (biased < IEEEBinary32.MAX_BIASED_EXPONENT) {
            const nonnegative = (this.encoding & IEEEBinary32.NOT_SIGN_MASK);
            if (sign === 0) {
                return nonnegative - 1;
            }
            return (IEEEBinary32.SIGN_MASK | (nonnegative + 1)) >>> 0;
        }
        if (trailing === 0) {
            if (sign === 0) {
                // The next-down of +INFINITY is +MAX_NORMAL.
                return IEEEBinary32.MAX_NORMAL;
            }
            // The next-down of -INFINITY is -INFINITY.
            return IEEEBinary32.NEG_INFINITY;
        }
        // The number is a quiet or signaling NaN, possibly with payload.
        // Just return the number itself.
        return this.encoding;
    }

    // Encode and decode the binary representation. The sign is 0 (number is
    // nonnegative) or 1 (number is negative). The biased exponent is in the
    // range [0, MAX_BIASED_EXPONENT]. The trailing significand is in the
    // range [0, MAX_TRAILING].
    getSign(): number {
        return (this.encoding & IEEEBinary32.SIGN_MASK) >>> IEEEBinary32.SIGN_SHIFT;
    }

    getBiased(): number {
        return (this.encoding & IEEEBinary32.EXPONENT_MASK) >>> IEEEBinary32.NUM_TRAILING_BITS;
    }

    getTrailing(): number {
        return this.encoding & IEEEBinary32.TRAILING_MASK;
    }

    setEncoding(sign: number, biased: number, trailing: number): void {
        this.encoding = ((sign << IEEEBinary32.SIGN_SHIFT)
            | (biased << IEEEBinary32.NUM_TRAILING_BITS) | trailing) >>> 0;
    }

    getEncoding(): { sign: number, biased: number, trailing: number } {
        return {
            sign: this.getSign(),
            biased: this.getBiased(),
            trailing: this.getTrailing()
        };
    }
}

// The bit-level representation of IEEE binary64 (C++ double, the 'number'
// type). The encoding is an unsigned 64-bit integer stored in a 'bigint'.
export class IEEEBinary64 {
    // Special constants.
    static readonly NUM_ENCODING_BITS = 64;
    static readonly NUM_EXPONENT_BITS = 11;
    static readonly NUM_SIGNIFICAND_BITS = 53;
    static readonly NUM_TRAILING_BITS = 52;
    static readonly EXPONENT_BIAS = 1023;
    static readonly MAX_BIASED_EXPONENT = 2047;
    static readonly MIN_SUB_EXPONENT = -1022;
    static readonly MIN_EXPONENT = -1074;
    static readonly SIGN_SHIFT = 63;

    static readonly SIGN_MASK = 0x8000000000000000n;
    static readonly NOT_SIGN_MASK = 0x7FFFFFFFFFFFFFFFn;
    static readonly TRAILING_MASK = 0x000FFFFFFFFFFFFFn;
    static readonly EXPONENT_MASK = 0x7FF0000000000000n;
    static readonly NAN_QUIET_MASK = 0x0008000000000000n;
    static readonly NAN_PAYLOAD_MASK = 0x0007FFFFFFFFFFFFn;
    static readonly MAX_TRAILING = 0x000FFFFFFFFFFFFFn;
    static readonly SUP_TRAILING = 0x0010000000000000n;
    static readonly POS_ZERO = 0x0000000000000000n;
    static readonly NEG_ZERO = 0x8000000000000000n;
    static readonly MIN_SUBNORMAL = 0x0000000000000001n;
    static readonly MAX_SUBNORMAL = 0x000FFFFFFFFFFFFFn;
    static readonly MIN_NORMAL = 0x0010000000000000n;
    static readonly MAX_NORMAL = 0x7FEFFFFFFFFFFFFFn;
    static readonly POS_INFINITY = 0x7FF0000000000000n;
    static readonly NEG_INFINITY = 0xFFF0000000000000n;

    // The bits of the represented number as an unsigned 64-bit integer.
    encoding: bigint;

    constructor(encoding: bigint = 0n) {
        this.encoding = BigInt.asUintN(64, encoding);
    }

    // Construction from primitive elements (ports of the C++ constructor
    // overloads on UInt, Float and (sign, biased, trailing)).
    static fromEncoding(encoding: bigint): IEEEBinary64 {
        return new IEEEBinary64(encoding);
    }

    static fromNumber(value: number): IEEEBinary64 {
        const result = new IEEEBinary64();
        result.number = value;
        return result;
    }

    static fromParts(sign: number, biased: number, trailing: bigint): IEEEBinary64 {
        const result = new IEEEBinary64();
        result.setEncoding(sign, biased, trailing);
        return result;
    }

    // The represented floating-point value (the other member of the C++
    // union).
    get number(): number {
        scratch64U[0] = this.encoding;
        return scratch64F[0];
    }

    set number(value: number) {
        scratch64F[0] = value;
        this.encoding = scratch64U[0];
    }

    getClassification(): IEEEClassification {
        const sign = this.getSign();
        const biased = this.getBiased();
        const trailing = this.getTrailing();
        return classify(sign, biased, trailing !== 0n,
            IEEEBinary64.MAX_BIASED_EXPONENT,
            (trailing & IEEEBinary64.NAN_QUIET_MASK) !== 0n);
    }

    isZero(): boolean {
        return this.encoding === IEEEBinary64.POS_ZERO
            || this.encoding === IEEEBinary64.NEG_ZERO;
    }

    isSignMinus(): boolean {
        return (this.encoding & IEEEBinary64.SIGN_MASK) !== 0n;
    }

    isSubnormal(): boolean {
        return this.getBiased() === 0 && this.getTrailing() > 0n;
    }

    isNormal(): boolean {
        const biased = this.getBiased();
        return 0 < biased && biased < IEEEBinary64.MAX_BIASED_EXPONENT;
    }

    isFinite(): boolean {
        return this.getBiased() < IEEEBinary64.MAX_BIASED_EXPONENT;
    }

    isInfinite(): boolean {
        return this.getBiased() === IEEEBinary64.MAX_BIASED_EXPONENT
            && this.getTrailing() === 0n;
    }

    isNaN(): boolean {
        return this.getBiased() === IEEEBinary64.MAX_BIASED_EXPONENT
            && this.getTrailing() !== 0n;
    }

    isQuietNaN(): boolean {
        const trailing = this.getTrailing();
        return this.getBiased() === IEEEBinary64.MAX_BIASED_EXPONENT
            && (trailing & IEEEBinary64.NAN_QUIET_MASK) !== 0n;
    }

    isSignalingNaN(): boolean {
        const trailing = this.getTrailing();
        return this.getBiased() === IEEEBinary64.MAX_BIASED_EXPONENT
            && (trailing & IEEEBinary64.NAN_QUIET_MASK) === 0n
            && (trailing & IEEEBinary64.NAN_PAYLOAD_MASK) !== 0n;
    }

    // Get neighboring numbers (the encodings, as unsigned 64-bit bigints).
    getNextUp(): bigint {
        const sign = this.getSign();
        const biased = this.getBiased();
        const trailing = this.getTrailing();

        if (biased === 0) {
            if (trailing === 0n) {
                // The next-up for both -0 and +0 is MIN_SUBNORMAL.
                return IEEEBinary64.MIN_SUBNORMAL;
            }
            if (sign !== 0) {
                // When trailing is 1, 'this' is -MIN_SUBNORMAL and next-up
                // is -0.
                return IEEEBinary64.SIGN_MASK | (trailing - 1n);
            }
            // When trailing is MAX_TRAILING, 'this' is MAX_SUBNORMAL and
            // next-up is MIN_NORMAL.
            return trailing + 1n;
        }
        if (biased < IEEEBinary64.MAX_BIASED_EXPONENT) {
            const nonnegative = this.encoding & IEEEBinary64.NOT_SIGN_MASK;
            if (sign !== 0) {
                return IEEEBinary64.SIGN_MASK | (nonnegative - 1n);
            }
            return nonnegative + 1n;
        }
        if (trailing === 0n) {
            if (sign !== 0) {
                // The next-up of -INFINITY is -MAX_NORMAL.
                return IEEEBinary64.SIGN_MASK | IEEEBinary64.MAX_NORMAL;
            }
            // The next-up of +INFINITY is +INFINITY.
            return IEEEBinary64.POS_INFINITY;
        }
        // The number is a quiet or signaling NaN, possibly with payload.
        // Just return the number itself.
        return this.encoding;
    }

    getNextDown(): bigint {
        const sign = this.getSign();
        const biased = this.getBiased();
        const trailing = this.getTrailing();

        if (biased === 0) {
            if (trailing === 0n) {
                // The next-down for both -0 and +0 is -MIN_SUBNORMAL.
                return IEEEBinary64.SIGN_MASK | IEEEBinary64.MIN_SUBNORMAL;
            }
            if (sign === 0) {
                // When trailing is 1, 'this' is MIN_SUBNORMAL and next-down
                // is +0.
                return trailing - 1n;
            }
            // When trailing is MAX_TRAILING, 'this' is -MAX_SUBNORMAL and
            // next-down is -MIN_NORMAL.
            return IEEEBinary64.SIGN_MASK | (trailing + 1n);
        }
        if (biased < IEEEBinary64.MAX_BIASED_EXPONENT) {
            const nonnegative = this.encoding & IEEEBinary64.NOT_SIGN_MASK;
            if (sign === 0) {
                return nonnegative - 1n;
            }
            return IEEEBinary64.SIGN_MASK | (nonnegative + 1n);
        }
        if (trailing === 0n) {
            if (sign === 0) {
                // The next-down of +INFINITY is +MAX_NORMAL.
                return IEEEBinary64.MAX_NORMAL;
            }
            // The next-down of -INFINITY is -INFINITY.
            return IEEEBinary64.NEG_INFINITY;
        }
        // The number is a quiet or signaling NaN, possibly with payload.
        // Just return the number itself.
        return this.encoding;
    }

    // Encode and decode the binary representation. The sign is 0 (number is
    // nonnegative) or 1 (number is negative). The biased exponent is in the
    // range [0, MAX_BIASED_EXPONENT]. The trailing significand is in the
    // range [0n, MAX_TRAILING].
    getSign(): number {
        return Number((this.encoding & IEEEBinary64.SIGN_MASK) >> 63n);
    }

    getBiased(): number {
        return Number((this.encoding & IEEEBinary64.EXPONENT_MASK) >> 52n);
    }

    getTrailing(): bigint {
        return this.encoding & IEEEBinary64.TRAILING_MASK;
    }

    setEncoding(sign: number, biased: number, trailing: bigint): void {
        this.encoding = (BigInt(sign) << 63n) | (BigInt(biased) << 52n) | trailing;
    }

    getEncoding(): { sign: number, biased: number, trailing: bigint } {
        return {
            sign: this.getSign(),
            biased: this.getBiased(),
            trailing: this.getTrailing()
        };
    }
}
