// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IEEEBinary16.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The IEEE 754 binary16 ("half") floating-point type. Upstream is
//   class IEEEBinary16 : public IEEEBinary<int16_t, uint16_t, 16, 11>
// so the class inherits the bit-field constants, the classification queries
// and the next-up/next-down neighbors from the template base and adds the
// binary16 <-> binary32 conversions.
//
// Port notes:
//   - src/IEEEBinary.ts provides the two concrete instantiations
//     IEEEBinary32 and IEEEBinary64 rather than a template base, so
//     IEEEBinary16 repeats the base-class members for the 16-bit parameters
//     (NUM_ENCODING_BITS = 16, NUM_SIGNIFICAND_BITS = 11). The shared
//     Classification enum IEEEClassification is imported.
//   - The encoding is an unsigned 16-bit integer stored in a 'number' and
//     maintained with '& 0xFFFF'.
//   - The C++ constructor overloads on float, double and uint16_t are
//     ambiguous in TypeScript, so construction uses the static factories
//     fromEncoding, fromNumber and fromParts, as IEEEBinary32 does.
//   - The C++ conversion operators 'operator float()' and 'operator double()'
//     become the 'number' accessor pair. The stored value is exactly a
//     binary32 value (every binary16 is representable in binary32), so
//     'number' is the same value C++ obtains for either conversion.
//   - The free operators (unary minus, +, -, *, /, and the compound updates)
//     and the std:: math wrappers become static methods, because the module
//     names 'negate', 'add', 'sub', 'mul', 'div', 'sign', 'sqr', ... are
//     already owned by other files under the flat library export. The
//     arithmetic operators return a binary32 result as upstream does; the
//     ports apply Math.fround so the value matches the C++ float arithmetic.
//   - Convert32To16 and Convert16To32 are private upstream; the port exposes
//     them as static methods so the bit-level behavior can be tested.

import {
    IEEEBinary32, IEEEClassification
} from './IEEEBinary.js';
import { BitHacks } from './BitHacks.js';
import {
    atandivpi, atan2divpi, clamp, cospi, exp10, invsqrt, isign, saturate,
    sign, sinpi, sqr
} from './Functions.js';

// Scratch views for converting between binary32 values and their bit
// patterns.
const scratchBuffer = new ArrayBuffer(4);
const scratchF = new Float32Array(scratchBuffer);
const scratchU = new Uint32Array(scratchBuffer);

function bitsOfFloat32(value: number): number {
    scratchF[0] = value;
    return scratchU[0] >>> 0;
}

function float32OfBits(bits: number): number {
    scratchU[0] = bits >>> 0;
    return scratchF[0];
}

export class IEEEBinary16 {
    // Members of the base class IEEEBinary<int16_t, uint16_t, 16, 11>.
    static readonly NUM_ENCODING_BITS = 16;
    static readonly NUM_EXPONENT_BITS = 5;
    static readonly NUM_SIGNIFICAND_BITS = 11;
    static readonly NUM_TRAILING_BITS = 10;
    static readonly EXPONENT_BIAS = 15;
    static readonly MAX_BIASED_EXPONENT = 31;
    static readonly MIN_SUB_EXPONENT = -14;
    static readonly MIN_EXPONENT = -24;
    static readonly SIGN_SHIFT = 15;

    static readonly SIGN_MASK = 0x8000;
    static readonly NOT_SIGN_MASK = 0x7FFF;
    static readonly TRAILING_MASK = 0x03FF;
    static readonly EXPONENT_MASK = 0x7C00;
    static readonly NAN_QUIET_MASK = 0x0200;
    static readonly NAN_PAYLOAD_MASK = 0x01FF;
    static readonly MAX_TRAILING = 0x03FF;
    static readonly SUP_TRAILING = 0x0400;
    static readonly POS_ZERO = 0x0000;
    static readonly NEG_ZERO = 0x8000;
    static readonly MIN_SUBNORMAL = 0x0001;
    static readonly MAX_SUBNORMAL = 0x03FF;
    static readonly MIN_NORMAL = 0x0400;
    static readonly MAX_NORMAL = 0x7BFF;
    static readonly POS_INFINITY = 0x7C00;
    static readonly NEG_INFINITY = 0xFC00;

    // Encodings of special numbers on the "continuous 16-bit number line"
    // as 32-bit float numbers.
    private static readonly F16_AVR_MIN_SUB_ZER = 0x33000000;  // 2^{-25}
    private static readonly F16_MIN_SUB = 0x33800000;          // 2^{-24}
    private static readonly F16_MIN_NOR = 0x38800000;          // 2^{-14}
    private static readonly F16_MAX_NOR = 0x477FE000;          // 2^{16}*(1-2^{-11})
    private static readonly F16_AVR_MAX_NOR_INF = 0x477FF000;  // 2^{16}*(1-2^{-12})

    // The amount to shift when converting between signs of 16-bit and 32-bit
    // numbers.
    private static readonly CONVERSION_SIGN_SHIFT =
        IEEEBinary32.NUM_ENCODING_BITS - IEEEBinary16.NUM_ENCODING_BITS;

    // The amount to shift when converting between trailing significands of
    // 16-bit and 32-bit numbers.
    private static readonly CONVERSION_TRAILING_SHIFT =
        IEEEBinary32.NUM_SIGNIFICAND_BITS - IEEEBinary16.NUM_SIGNIFICAND_BITS;

    // The half value for round-to-nearest-ties-to-even. The fractional part
    // in the rounding is shifted left so that the leading bit is the
    // high-order bit of a 32-bit unsigned integer.
    private static readonly FRACTION_HALF = IEEEBinary32.SIGN_MASK;

    // The bits of the represented number as an unsigned 16-bit integer.
    encoding: number;

    constructor(encoding: number = 0) {
        this.encoding = encoding & 0xFFFF;
    }

    // Construction from primitive elements (the ports of the C++ constructor
    // overloads on uint16_t, float/double and (sign, biased, trailing)).
    static fromEncoding(encoding: number): IEEEBinary16 {
        return new IEEEBinary16(encoding);
    }

    static fromNumber(value: number): IEEEBinary16 {
        const result = new IEEEBinary16();
        result.number = value;
        return result;
    }

    static fromParts(sign: number, biased: number, trailing: number): IEEEBinary16 {
        const result = new IEEEBinary16();
        result.setEncoding(sign, biased, trailing);
        return result;
    }

    clone(): IEEEBinary16 {
        return new IEEEBinary16(this.encoding);
    }

    // The ports of 'operator float()' and 'operator double()'.
    get number(): number {
        return float32OfBits(IEEEBinary16.convert16To32(this.encoding));
    }

    set number(value: number) {
        this.encoding = IEEEBinary16.convert32To16(bitsOfFloat32(value));
    }

    // Classification of the number.
    getClassification(): IEEEClassification {
        const sign = this.getSign();
        const biased = this.getBiased();
        const trailing = this.getTrailing();

        if (biased === 0) {
            if (trailing === 0) {
                return sign !== 0 ? IEEEClassification.NEG_ZERO : IEEEClassification.POS_ZERO;
            }
            return sign !== 0 ? IEEEClassification.NEG_SUBNORMAL : IEEEClassification.POS_SUBNORMAL;
        }
        if (biased < IEEEBinary16.MAX_BIASED_EXPONENT) {
            return sign !== 0 ? IEEEClassification.NEG_NORMAL : IEEEClassification.POS_NORMAL;
        }
        if (trailing === 0) {
            return sign !== 0 ? IEEEClassification.NEG_INFINITY : IEEEClassification.POS_INFINITY;
        }
        return (trailing & IEEEBinary16.NAN_QUIET_MASK) !== 0
            ? IEEEClassification.QUIET_NAN : IEEEClassification.SIGNALING_NAN;
    }

    isZero(): boolean {
        return this.encoding === IEEEBinary16.POS_ZERO
            || this.encoding === IEEEBinary16.NEG_ZERO;
    }

    isSignMinus(): boolean {
        return this.getSign() === 1;
    }

    isSubnormal(): boolean {
        return this.getBiased() === 0 && this.getTrailing() > 0;
    }

    isNormal(): boolean {
        const biased = this.getBiased();
        return biased > 0 && biased < IEEEBinary16.MAX_BIASED_EXPONENT;
    }

    isFinite(): boolean {
        return this.getBiased() < IEEEBinary16.MAX_BIASED_EXPONENT;
    }

    isInfinite(): boolean {
        return this.getBiased() === IEEEBinary16.MAX_BIASED_EXPONENT
            && this.getTrailing() === 0;
    }

    isNaN(): boolean {
        return this.getBiased() === IEEEBinary16.MAX_BIASED_EXPONENT
            && this.getTrailing() !== 0;
    }

    isQuietNaN(): boolean {
        return this.getBiased() === IEEEBinary16.MAX_BIASED_EXPONENT
            && (this.getTrailing() & IEEEBinary16.NAN_QUIET_MASK) !== 0;
    }

    isSignalingNaN(): boolean {
        const trailing = this.getTrailing();
        return this.getBiased() === IEEEBinary16.MAX_BIASED_EXPONENT
            && (trailing & IEEEBinary16.NAN_QUIET_MASK) === 0
            && (trailing & IEEEBinary16.NAN_PAYLOAD_MASK) !== 0;
    }

    // Get neighboring numbers (the encodings, as unsigned 16-bit integers).
    getNextUp(): number {
        const sign = this.getSign();
        const biased = this.getBiased();
        const trailing = this.getTrailing();

        if (biased === 0) {
            if (trailing === 0) {
                // The next-up for both -0 and +0 is MIN_SUBNORMAL.
                return IEEEBinary16.MIN_SUBNORMAL;
            }
            if (sign !== 0) {
                // When trailing is 1, 'this' is -MIN_SUBNORMAL and next-up
                // is -0.
                return (IEEEBinary16.SIGN_MASK | (trailing - 1)) & 0xFFFF;
            }
            // When trailing is MAX_TRAILING, 'this' is MAX_SUBNORMAL and
            // next-up is MIN_NORMAL.
            return trailing + 1;
        }
        if (biased < IEEEBinary16.MAX_BIASED_EXPONENT) {
            const nonnegative = this.encoding & IEEEBinary16.NOT_SIGN_MASK;
            if (sign !== 0) {
                return (IEEEBinary16.SIGN_MASK | (nonnegative - 1)) & 0xFFFF;
            }
            return nonnegative + 1;
        }
        if (trailing === 0) {
            if (sign !== 0) {
                // The next-up of -INFINITY is -MAX_NORMAL.
                return (IEEEBinary16.SIGN_MASK | IEEEBinary16.MAX_NORMAL) & 0xFFFF;
            }
            // The next-up of +INFINITY is +INFINITY.
            return IEEEBinary16.POS_INFINITY;
        }
        // The number is a quiet or signaling NaN, possibly with payload.
        // Just return the number itself.
        return this.encoding;
    }

    getNextDown(): number {
        const sign = this.getSign();
        const biased = this.getBiased();
        const trailing = this.getTrailing();

        if (biased === 0) {
            if (trailing === 0) {
                // The next-down for both -0 and +0 is -MIN_SUBNORMAL.
                return (IEEEBinary16.SIGN_MASK | IEEEBinary16.MIN_SUBNORMAL) & 0xFFFF;
            }
            if (sign === 0) {
                // When trailing is 1, 'this' is MIN_SUBNORMAL and next-down
                // is +0.
                return trailing - 1;
            }
            // When trailing is MAX_TRAILING, 'this' is -MAX_SUBNORMAL and
            // next-down is -MIN_NORMAL.
            return (IEEEBinary16.SIGN_MASK | (trailing + 1)) & 0xFFFF;
        }
        if (biased < IEEEBinary16.MAX_BIASED_EXPONENT) {
            const nonnegative = this.encoding & IEEEBinary16.NOT_SIGN_MASK;
            if (sign === 0) {
                return nonnegative - 1;
            }
            return (IEEEBinary16.SIGN_MASK | (nonnegative + 1)) & 0xFFFF;
        }
        if (trailing === 0) {
            if (sign === 0) {
                // The next-down of +INFINITY is +MAX_NORMAL.
                return IEEEBinary16.MAX_NORMAL;
            }
            // The next-down of -INFINITY is -INFINITY.
            return IEEEBinary16.NEG_INFINITY;
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
        return (this.encoding & IEEEBinary16.SIGN_MASK) >>> IEEEBinary16.SIGN_SHIFT;
    }

    getBiased(): number {
        return (this.encoding & IEEEBinary16.EXPONENT_MASK) >>> IEEEBinary16.NUM_TRAILING_BITS;
    }

    getTrailing(): number {
        return this.encoding & IEEEBinary16.TRAILING_MASK;
    }

    setEncoding(sign: number, biased: number, trailing: number): void {
        this.encoding = ((sign << IEEEBinary16.SIGN_SHIFT)
            | (biased << IEEEBinary16.NUM_TRAILING_BITS) | trailing) & 0xFFFF;
    }

    getEncoding(): { sign: number, biased: number, trailing: number } {
        return {
            sign: this.getSign(),
            biased: this.getBiased(),
            trailing: this.getTrailing()
        };
    }

    // Support for conversions between 16-bit and 32-bit numbers. The inputs
    // and outputs are bit encodings, not numeric values.
    static convert32To16(inEncoding: number): number {
        // In the comments of this function, x refers to the 32-bit
        // floating-point number corresponding to inEncoding and y refers to
        // the 16-bit floating-point number that x is converted to.

        const enc = inEncoding >>> 0;

        // Extract the channels for x.
        const sign32 = (enc & IEEEBinary32.SIGN_MASK) >>> 0;
        const biased32 = (enc & IEEEBinary32.EXPONENT_MASK) >>> IEEEBinary32.NUM_TRAILING_BITS;
        let trailing32 = (enc & IEEEBinary32.TRAILING_MASK) >>> 0;
        const nonneg32 = (enc & IEEEBinary32.NOT_SIGN_MASK) >>> 0;

        // Generate the channels for y.
        const sign16 = (sign32 >>> IEEEBinary16.CONVERSION_SIGN_SHIFT) & 0xFFFF;
        let biased16: number;
        let trailing16: number;
        let frcpart: number;

        if (biased32 === 0) {
            // x is zero or 32-subnormal, the nearest y is zero.
            return sign16;
        }

        if (biased32 < IEEEBinary32.MAX_BIASED_EXPONENT) {
            // x is 32-normal.
            if (nonneg32 <= IEEEBinary16.F16_AVR_MIN_SUB_ZER) {
                // x <= 2^{-25}, the nearest y is zero.
                return sign16;
            }

            if (nonneg32 <= IEEEBinary16.F16_MIN_SUB) {
                // 2^{-25} < x <= 2^{-24}, the nearest y is 16-min-subnormal.
                return (sign16 | IEEEBinary16.MIN_SUBNORMAL) & 0xFFFF;
            }

            if (nonneg32 < IEEEBinary16.F16_MIN_NOR) {
                // 2^{-24} < x < 2^{-14}, compute nearest 16-bit subnormal y
                // using round-to-nearest-ties-to-even.
                //
                // y = 0.s9 ... s0 * 2^{14}
                // x = 1.t22 ... t0 * 2^e, where -24 <= e <= -15
                //   = (0.1 t22 ... t0 * 2^{e+15}) * 2^{-14}
                //   = (0.1 t22 ... t0 * 2^p) * 2^{-14}
                // where p = e+15 with -9 <= p <= 0. The term
                // (0.1 t22 ... t0 * 2^p) must be rounded to
                // 0.s9 ... s0 * 2^{-14}.
                const p = biased32 - IEEEBinary32.EXPONENT_BIAS + IEEEBinary16.EXPONENT_BIAS;

                // x is 32-normal, so there is an implied 1-bit that must
                // first be appended to the 32-trailing significand to obtain
                // all the bits necessary for the 16-trailing significand for
                // the 16-subnormal y. The resulting number is
                // 000000001 t22 ... t0.
                trailing32 = (trailing32 | IEEEBinary32.SUP_TRAILING) >>> 0;

                // Get the integer part.
                const rshift = -IEEEBinary16.MIN_SUB_EXPONENT - p;
                trailing16 = (trailing32 >>> rshift) & 0xFFFF;

                // Get the fractional part.
                const lshift = IEEEBinary32.NUM_ENCODING_BITS + IEEEBinary16.MIN_SUB_EXPONENT + p;
                frcpart = (trailing32 << lshift) >>> 0;

                // Round to nearest with ties to even.
                if (frcpart > IEEEBinary16.FRACTION_HALF
                    || (frcpart === IEEEBinary16.FRACTION_HALF && (trailing16 & 1) !== 0)) {
                    // If there is a carry into the exponent, the nearest is
                    // actually 16-min-normal 1.0*2^{-14}, so the high-order
                    // bit of trailing16 makes biased16 equal to 1 and the
                    // result is correct.
                    ++trailing16;
                }
                return (sign16 | trailing16) & 0xFFFF;
            }

            if (nonneg32 <= IEEEBinary16.F16_MAX_NOR) {
                // 2^{-14} <= x <= 1.1111111111*2^{15}, compute nearest
                // 16-bit normal y using round-to-nearest-ties-to-even.

                // The exponents of x and y are the same, although the biased
                // exponents are different because of different
                // exponent-bias parameters.
                const e = biased32 - IEEEBinary32.EXPONENT_BIAS;
                biased16 = ((e + IEEEBinary16.EXPONENT_BIAS) << IEEEBinary16.NUM_TRAILING_BITS)
                    & 0xFFFF;

                // Let x = 1.t22...t0 * 2^e and y = 1.s9...s0 * 2^e. Both x
                // and y have an implied leading 1-bit (both are normal), so
                // we can ignore it. The number 0.t22...t0 must be rounded to
                // the number 0.s9...s0.

                // Get the integer part.
                trailing16 = (trailing32 >>> IEEEBinary16.CONVERSION_TRAILING_SHIFT) & 0xFFFF;

                // Get the fractional part.
                const lshift = IEEEBinary32.NUM_ENCODING_BITS + IEEEBinary16.MIN_SUB_EXPONENT + 1;
                frcpart = (trailing32 << lshift) >>> 0;

                // Round to nearest with ties to even.
                if (frcpart > IEEEBinary16.FRACTION_HALF
                    || (frcpart === IEEEBinary16.FRACTION_HALF && (trailing16 & 1) !== 0)) {
                    // If there is a carry into the exponent, the addition of
                    // trailing16 to biased16 (rather than or-ing) produces
                    // the correct result.
                    ++trailing16;
                }
                return (sign16 | (biased16 + trailing16)) & 0xFFFF;
            }

            if (nonneg32 < IEEEBinary16.F16_AVR_MAX_NOR_INF) {
                // 1.1111111111*2^{15} < x < (MAX_NORMAL+INFINITY)/2, so the
                // number is closest to 16-max-normal.
                return (sign16 | IEEEBinary16.MAX_NORMAL) & 0xFFFF;
            }

            // nonneg32 >= (MAX_NORMAL+INFINITY)/2, so convert to 16-infinite.
            return (sign16 | IEEEBinary16.POS_INFINITY) & 0xFFFF;
        }

        if (trailing32 === 0) {
            // The number is 32-infinite. Convert to 16-infinite.
            return (sign16 | IEEEBinary16.POS_INFINITY) & 0xFFFF;
        }

        // The number is 32-NaN. Convert to 16-NaN with 16-payload the
        // high-order 9 bits of the 32-payload. The 32-quiet-NaN mask bit is
        // copied in the conversion.
        const maskPayload = (trailing32 >>> IEEEBinary16.CONVERSION_TRAILING_SHIFT) & 0xFFFF;
        return (sign16 | IEEEBinary16.EXPONENT_MASK | maskPayload) & 0xFFFF;
    }

    static convert16To32(inEncoding: number): number {
        const enc = inEncoding & 0xFFFF;

        // Extract the channels for the IEEEBinary16 number.
        const sign16 = enc & IEEEBinary16.SIGN_MASK;
        const biased16 = (enc & IEEEBinary16.EXPONENT_MASK) >>> IEEEBinary16.NUM_TRAILING_BITS;
        const trailing16 = enc & IEEEBinary16.TRAILING_MASK;

        // Generate the channels for the binary32 number.
        const sign32 = (sign16 << IEEEBinary16.CONVERSION_SIGN_SHIFT) >>> 0;
        let biased32: number;
        let trailing32: number;

        if (biased16 === 0) {
            if (trailing16 === 0) {
                // The number is 16-zero. Convert to 32-zero.
                return sign32;
            }
            // The number is 16-subnormal. Convert to 32-normal.
            trailing32 = trailing16;
            const leading = BitHacks.getLeadingBit(trailing32);
            const shift = 23 - leading;
            biased32 = IEEEBinary32.EXPONENT_BIAS - 1 - shift;
            trailing32 = ((trailing32 << shift) & IEEEBinary32.TRAILING_MASK) >>> 0;
            return (sign32 | (biased32 << IEEEBinary32.NUM_TRAILING_BITS) | trailing32) >>> 0;
        }

        if (biased16 < IEEEBinary16.MAX_BIASED_EXPONENT) {
            // The number is 16-normal. Convert to 32-normal.
            biased32 = biased16 - IEEEBinary16.EXPONENT_BIAS + IEEEBinary32.EXPONENT_BIAS;
            trailing32 = (trailing16 << IEEEBinary16.CONVERSION_TRAILING_SHIFT) >>> 0;
            return (sign32 | (biased32 << IEEEBinary32.NUM_TRAILING_BITS) | trailing32) >>> 0;
        }

        if (trailing16 === 0) {
            // The number is 16-infinite. Convert to 32-infinite.
            return (sign32 | IEEEBinary32.EXPONENT_MASK) >>> 0;
        }

        // The number is 16-NaN. Convert to 32-NaN with 32-payload whose
        // high-order 9 bits are from the 16-payload. The 16-quiet-NaN mask
        // bit is copied in the conversion.
        const maskPayload = (trailing16 << IEEEBinary16.CONVERSION_TRAILING_SHIFT) >>> 0;
        return (sign32 | IEEEBinary32.EXPONENT_MASK | maskPayload) >>> 0;
    }

    // Comparison. Upstream compares the values after conversion to float.
    static equals(x: IEEEBinary16, y: IEEEBinary16): boolean {
        return x.number === y.number;
    }

    static notEquals(x: IEEEBinary16, y: IEEEBinary16): boolean {
        return x.number !== y.number;
    }

    static lessThan(x: IEEEBinary16, y: IEEEBinary16): boolean {
        return x.number < y.number;
    }

    static lessThanOrEqual(x: IEEEBinary16, y: IEEEBinary16): boolean {
        return x.number <= y.number;
    }

    static greaterThan(x: IEEEBinary16, y: IEEEBinary16): boolean {
        return x.number > y.number;
    }

    static greaterThanOrEqual(x: IEEEBinary16, y: IEEEBinary16): boolean {
        return x.number >= y.number;
    }

    // Arithmetic operations (high-precision). The unary minus flips the sign
    // bit and stays in binary16; the binary operators produce a binary32
    // result, as upstream does.
    static negate(x: IEEEBinary16): IEEEBinary16 {
        return new IEEEBinary16((x.encoding ^ IEEEBinary16.SIGN_MASK) & 0xFFFF);
    }

    static add(x: IEEEBinary16 | number, y: IEEEBinary16 | number): number {
        return Math.fround(toFloat(x) + toFloat(y));
    }

    static sub(x: IEEEBinary16 | number, y: IEEEBinary16 | number): number {
        return Math.fround(toFloat(x) - toFloat(y));
    }

    static mul(x: IEEEBinary16 | number, y: IEEEBinary16 | number): number {
        return Math.fround(toFloat(x) * toFloat(y));
    }

    static div(x: IEEEBinary16 | number, y: IEEEBinary16 | number): number {
        return Math.fround(toFloat(x) / toFloat(y));
    }

    // Arithmetic updates. Upstream mutates the left operand in place, so the
    // ports mutate 'x' and return it.
    static addAssign(x: IEEEBinary16, y: IEEEBinary16 | number): IEEEBinary16 {
        x.number = IEEEBinary16.add(x, y);
        return x;
    }

    static subAssign(x: IEEEBinary16, y: IEEEBinary16 | number): IEEEBinary16 {
        x.number = IEEEBinary16.sub(x, y);
        return x;
    }

    static mulAssign(x: IEEEBinary16, y: IEEEBinary16 | number): IEEEBinary16 {
        x.number = IEEEBinary16.mul(x, y);
        return x;
    }

    static divAssign(x: IEEEBinary16, y: IEEEBinary16 | number): IEEEBinary16 {
        x.number = IEEEBinary16.div(x, y);
        return x;
    }

    // The ports of the std:: math overloads for IEEEBinary16. Upstream
    // evaluates in binary32 and converts the result back to binary16; the
    // ports apply Math.fround to the double-precision result for the same
    // effect.
    static acos(x: IEEEBinary16): IEEEBinary16 { return half(Math.acos(x.number)); }
    static acosh(x: IEEEBinary16): IEEEBinary16 { return half(Math.acosh(x.number)); }
    static asin(x: IEEEBinary16): IEEEBinary16 { return half(Math.asin(x.number)); }
    static asinh(x: IEEEBinary16): IEEEBinary16 { return half(Math.asinh(x.number)); }
    static atan(x: IEEEBinary16): IEEEBinary16 { return half(Math.atan(x.number)); }
    static atanh(x: IEEEBinary16): IEEEBinary16 { return half(Math.atanh(x.number)); }
    static atan2(y: IEEEBinary16, x: IEEEBinary16): IEEEBinary16 {
        return half(Math.atan2(y.number, x.number));
    }
    static ceil(x: IEEEBinary16): IEEEBinary16 { return half(Math.ceil(x.number)); }
    static cos(x: IEEEBinary16): IEEEBinary16 { return half(Math.cos(x.number)); }
    static cosh(x: IEEEBinary16): IEEEBinary16 { return half(Math.cosh(x.number)); }
    static exp(x: IEEEBinary16): IEEEBinary16 { return half(Math.exp(x.number)); }
    static exp2(x: IEEEBinary16): IEEEBinary16 { return half(Math.pow(2, x.number)); }
    static fabs(x: IEEEBinary16): IEEEBinary16 { return half(Math.abs(x.number)); }
    static floor(x: IEEEBinary16): IEEEBinary16 { return half(Math.floor(x.number)); }
    static fmod(x: IEEEBinary16, y: IEEEBinary16): IEEEBinary16 {
        return half(x.number % y.number);
    }
    // The port of 'std::frexp(IEEEBinary16 x, int32_t* exponent)'. The
    // upstream output parameter becomes a field of the returned object, the
    // BSNumber.frexp/BSRational.frexp precedent. The decomposition is
    // x = result * 2^exponent with |result| in [1/2, 1); zero, infinity and
    // NaN are returned unchanged with exponent 0, as std::frexp does.
    static frexp(x: IEEEBinary16): { result: IEEEBinary16, exponent: number } {
        const value = x.number;
        if (value === 0 || !Number.isFinite(value)) {
            return { result: x.clone(), exponent: 0 };
        }
        // Every binary16 value is exactly a binary32 normal (a 16-subnormal
        // converts to a 32-normal), so the biased exponent of the binary32
        // encoding determines the decomposition exactly.
        const biased = (bitsOfFloat32(value) & IEEEBinary32.EXPONENT_MASK)
            >>> IEEEBinary32.NUM_TRAILING_BITS;
        const exponent = biased - IEEEBinary32.EXPONENT_BIAS + 1;
        return { result: half(value * Math.pow(2, -exponent)), exponent };
    }
    static ldexp(x: IEEEBinary16, exponent: number): IEEEBinary16 {
        return half(x.number * Math.pow(2, exponent));
    }
    static log(x: IEEEBinary16): IEEEBinary16 { return half(Math.log(x.number)); }
    static log2(x: IEEEBinary16): IEEEBinary16 { return half(Math.log2(x.number)); }
    static log10(x: IEEEBinary16): IEEEBinary16 { return half(Math.log10(x.number)); }
    static pow(x: IEEEBinary16, y: IEEEBinary16): IEEEBinary16 {
        return half(Math.pow(x.number, y.number));
    }
    static sin(x: IEEEBinary16): IEEEBinary16 { return half(Math.sin(x.number)); }
    static sinh(x: IEEEBinary16): IEEEBinary16 { return half(Math.sinh(x.number)); }
    static sqrt(x: IEEEBinary16): IEEEBinary16 { return half(Math.sqrt(x.number)); }
    static tan(x: IEEEBinary16): IEEEBinary16 { return half(Math.tan(x.number)); }
    static tanh(x: IEEEBinary16): IEEEBinary16 { return half(Math.tanh(x.number)); }

    // The ports of the gte:: Functions.h overloads for IEEEBinary16.
    static atandivpi(x: IEEEBinary16): IEEEBinary16 { return half(atandivpi(x.number)); }
    static atan2divpi(y: IEEEBinary16, x: IEEEBinary16): IEEEBinary16 {
        return half(atan2divpi(y.number, x.number));
    }
    static clamp(x: IEEEBinary16, xmin: IEEEBinary16, xmax: IEEEBinary16): IEEEBinary16 {
        return half(clamp(x.number, xmin.number, xmax.number));
    }
    static cospi(x: IEEEBinary16): IEEEBinary16 { return half(cospi(x.number)); }
    static exp10(x: IEEEBinary16): IEEEBinary16 { return half(exp10(x.number)); }
    static invsqrt(x: IEEEBinary16): IEEEBinary16 { return half(invsqrt(x.number)); }
    static isign(x: IEEEBinary16): number { return isign(x.number); }
    static saturate(x: IEEEBinary16): IEEEBinary16 { return half(saturate(x.number)); }
    static sign(x: IEEEBinary16): IEEEBinary16 { return half(sign(x.number)); }
    static sinpi(x: IEEEBinary16): IEEEBinary16 { return half(sinpi(x.number)); }
    static sqr(x: IEEEBinary16): IEEEBinary16 { return half(sqr(x.number)); }
}

function toFloat(x: IEEEBinary16 | number): number {
    return typeof x === 'number' ? x : x.number;
}

function half(value: number): IEEEBinary16 {
    return IEEEBinary16.fromNumber(Math.fround(value));
}
