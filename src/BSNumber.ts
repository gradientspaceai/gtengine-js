// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BSNumber.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The class BSNumber (binary scientific number) provides exact arithmetic for
// robust algorithms, typically those for which we need to know the exact sign
// of determinants. A nonzero number is represented as
//   x = sign * uinteger * 2^{biasedExponent}
// where 'uinteger' is a positive odd integer. The number zero is represented
// by sign = 0, biasedExponent = 0 and uinteger = 0.
//
// Because the unsigned integer is always odd, the representation is unique.
// Addition, subtraction and multiplication of BSNumber objects are exact; the
// number of bits of the unsigned integer grows as needed. Division is not
// supported (that is BSRational's job), so hasDivisionOperator is false.
//
// The document at the following link describes the design, implementation and
// use of BSNumber and BSRational.
//   https://www.geometrictools.com/Documentation/ArbitraryPrecision.pdf
//
// Port notes:
// - Upstream is 'template <typename UInteger> class BSNumber', where UInteger
//   is one of UIntegerAP32 (arbitrary precision) or UIntegerFP32<N> (fixed
//   precision), both built on 32-bit blocks. Per PORTING.md, the port uses the
//   JavaScript 'bigint' as the unsigned-integer storage, so the template
//   parameter disappears and BSNumber is arbitrary precision. The UInteger
//   interface maps as follows.
//     GetNumBits()            -> numBitsOf(u), the bit length of the bigint
//     operator==, operator<   -> bigint === and <
//     Add, Sub, Mul           -> bigint +, -, *
//     ShiftLeft(n, s)         -> n << BigInt(s)
//     ShiftRightToOdd(n)      -> shiftRightToOdd(n) below
//     RoundUp()               -> shiftRightToOdd(w + 1n)
//     GetPrefix(numRequested) -> getPrefix() below
//   The fixed-precision variant (and hence UInteger::GetMaxSize(), used by
//   upstream Convert to reject over-large precisions) has no analogue; bigint
//   is always arbitrary precision, so that check is dropped.
// - C++ has both BSNumber(float) and BSNumber(double); the port has
//   BSNumber.fromNumber (the binary64 conversion, which is exact for every
//   'number' including all the integer constructors) and
//   BSNumber.fromFloat32 (the binary32 conversion, for the float path).
//   Likewise operator float / operator double become toFloat32() / toNumber().
//   The int64_t/uint64_t constructors, which need more than 53 bits of
//   integer, become BSNumber.fromBigInt.
// - Operators become named methods: ==, <, ... -> equals, lessThan, ...;
//   +, -, * -> add, sub, mul; unary - -> negated() (returning a new object),
//   distinct from the upstream in-place Negate() -> negate().
// - Upstream Convert(input, precision, mode, output) becomes the exported
//   convertBSNumber(input, precision, mode) returning the result; the <cfenv>
//   rounding-mode flags become the BSNumberRoundingMode enum.
// - The 'std::' math functions (acos, sqrt, floor, ...) and the gte:: helper
//   functions (clamp, sign, sqr, ...) that upstream overloads for BSNumber
//   round-trip through double. They are ported as static methods of BSNumber
//   rather than module functions because the module-level names are already
//   taken by the 'number' versions in Functions.ts (see the global export
//   uniqueness rule in PORTING.md).
// - The GTE_VALIDATE_BSNUMBER helper IsValid() is always available as
//   isValid(); it is not called automatically (upstream calls it only when the
//   conditional define is enabled), but it is useful in tests.
// - Write/Read serialize the C++ block representation to a binary stream; the
//   bigint storage has no such layout, so they are omitted.

import { logAssert, logError } from './Logger.js';
import { IEEEBinary32, IEEEBinary64 } from './IEEEBinary.js';
import type { ArbitraryPrecisionNumber } from './TypeTraits.js';
import {
    atandivpi, atan2divpi, clamp, cospi, exp10, invsqrt, isign, saturate,
    sign, sinpi, sqr
} from './Functions.js';

// The number of bits of a nonnegative bigint; the port of
// UInteger::GetNumBits(). The bit length of zero is 0.
function numBitsOf(u: bigint): number {
    if (u <= 0n) {
        return 0;
    }
    // Convert to a hexadecimal string, which is exact and fast for large
    // numbers, then account for the leading nibble.
    const hex = u.toString(16);
    let numBits = 4 * (hex.length - 1);
    let leading = Number.parseInt(hex[0], 16);
    while (leading > 0) {
        ++numBits;
        leading >>= 1;
    }
    return numBits;
}

// The number of trailing zero bits of a positive bigint; the shift returned
// by UInteger::ShiftRightToOdd(...).
function trailingZerosOf(u: bigint): number {
    if (u === 0n) {
        return 0;
    }
    let shift = 0;
    // Strip 32 bits at a time, then the remaining bits one at a time.
    while ((u & 0xFFFFFFFFn) === 0n) {
        u >>= 32n;
        shift += 32;
    }
    let block = Number(u & 0xFFFFFFFFn);
    while ((block & 1) === 0) {
        block >>>= 1;
        ++shift;
    }
    return shift;
}

// The port of 'shift = result.ShiftRightToOdd(number)': the odd number
// obtained by removing the trailing zero bits, together with the number of
// bits removed.
function shiftRightToOdd(u: bigint): { odd: bigint, shift: number } {
    const shift = trailingZerosOf(u);
    return { odd: u >> BigInt(shift), shift };
}

// The rounding modes of upstream convertBSNumber. Upstream uses the <cfenv>
// flags FE_TONEAREST, FE_DOWNWARD, FE_TOWARDZERO and FE_UPWARD, whose values
// are implementation defined; the names are preserved.
export enum BSNumberRoundingMode {
    FE_TONEAREST = 0,
    FE_DOWNWARD = 1,
    FE_TOWARDZERO = 2,
    FE_UPWARD = 3
}

export class BSNumber implements ArbitraryPrecisionNumber {
    // See TypeTraits.ts. BSNumber is arbitrary precision and, unlike
    // BSRational, has no division operator.
    readonly isArbitraryPrecision = true as const;
    readonly hasDivisionOperator = false;

    // The number 0 is represented by mSign = 0, mBiasedExponent = 0 and
    // mUInteger = 0n. For nonzero numbers, mSign is +1 or -1 and mUInteger is
    // a positive odd integer.
    private mSign: number;
    private mBiasedExponent: number;
    private mUInteger: bigint;

    // The default constructor generates the zero BSNumber.
    constructor() {
        this.mSign = 0;
        this.mBiasedExponent = 0;
        this.mUInteger = 0n;
    }

    // Construction from the raw representation. WARNING: as with the upstream
    // block of SetSign/SetBiasedExponent/GetUInteger().CopyFrom calls, the
    // caller is responsible for supplying a valid representation (mUInteger
    // positive and odd when sign != 0, all members zero when sign == 0).
    static fromParts(sign: number, biasedExponent: number, uinteger: bigint): BSNumber {
        const result = new BSNumber();
        result.mSign = sign;
        result.mBiasedExponent = biasedExponent;
        result.mUInteger = uinteger;
        return result;
    }

    // Conversion from a 'number' using the IEEE binary64 representation. This
    // is exact: every finite double is a binary scientific number. It is the
    // port of BSNumber(double) and also of the int32_t/uint32_t constructors,
    // whose inputs are exactly representable as doubles.
    static fromNumber(value: number): BSNumber {
        const x = IEEEBinary64.fromNumber(value);
        const result = new BSNumber();
        result.convertFrom(x.getSign(), x.getBiased(), x.getTrailing(),
            IEEEBinary64.NUM_TRAILING_BITS, IEEEBinary64.EXPONENT_BIAS,
            IEEEBinary64.MAX_BIASED_EXPONENT, IEEEBinary64.MIN_SUB_EXPONENT,
            IEEEBinary64.SUP_TRAILING);
        return result;
    }

    // Conversion from a 'number' rounded to the IEEE binary32 representation;
    // the port of BSNumber(float).
    static fromFloat32(value: number): BSNumber {
        const x = IEEEBinary32.fromNumber(value);
        const result = new BSNumber();
        result.convertFrom(x.getSign(), x.getBiased(), BigInt(x.getTrailing()),
            IEEEBinary32.NUM_TRAILING_BITS, IEEEBinary32.EXPONENT_BIAS,
            IEEEBinary32.MAX_BIASED_EXPONENT, IEEEBinary32.MIN_SUB_EXPONENT,
            BigInt(IEEEBinary32.SUP_TRAILING));
        return result;
    }

    // Exact conversion from an integer of any size; the port of the
    // int64_t/uint64_t constructors, which use BitHacks::GetTrailingBit to
    // normalize the magnitude to an odd number.
    static fromBigInt(value: bigint): BSNumber {
        const result = new BSNumber();
        if (value === 0n) {
            return result;
        }
        const magnitude = (value < 0n ? -value : value);
        result.mSign = (value < 0n ? -1 : 1);
        const { odd, shift } = shiftRightToOdd(magnitude);
        result.mBiasedExponent = shift;
        result.mUInteger = odd;
        return result;
    }

    // The number must be of the form "x" or "+x" or "-x", where x is a
    // positive integer with nonzero leading digit.
    static fromString(numberString: string): BSNumber {
        logAssert(numberString.length > 0, 'A number must be specified.');

        // Get the leading '+' or '-' if it exists.
        let intNumber: string;
        let signValue: number;
        if (numberString[0] === '+') {
            intNumber = numberString.substring(1);
            signValue = +1;
            logAssert(intNumber.length > 0, 'Invalid number format.');
        } else if (numberString[0] === '-') {
            intNumber = numberString.substring(1);
            signValue = -1;
            logAssert(intNumber.length > 0, 'Invalid number format.');
        } else {
            intNumber = numberString;
            signValue = +1;
        }

        const result = BSNumber.convertToInteger(intNumber);
        if (result.mSign !== 0) {
            result.mSign = signValue;
        }
        return result;
    }

    // C++ copy assignment; TypeScript objects alias, so clone explicitly.
    clone(): BSNumber {
        return BSNumber.fromParts(this.mSign, this.mBiasedExponent, this.mUInteger);
    }

    // Conversions to floating-point. These always use the default rounding
    // mode, round-to-nearest-ties-to-even.
    toNumber(): number {
        const { s, e, t } = this.convertTo(
            IEEEBinary64.NUM_SIGNIFICAND_BITS, IEEEBinary64.EXPONENT_BIAS,
            IEEEBinary64.MAX_BIASED_EXPONENT, IEEEBinary64.MIN_EXPONENT,
            IEEEBinary64.MIN_SUB_EXPONENT, IEEEBinary64.SUP_TRAILING);
        return IEEEBinary64.fromParts(s, e, t).number;
    }

    toFloat32(): number {
        const { s, e, t } = this.convertTo(
            IEEEBinary32.NUM_SIGNIFICAND_BITS, IEEEBinary32.EXPONENT_BIAS,
            IEEEBinary32.MAX_BIASED_EXPONENT, IEEEBinary32.MIN_EXPONENT,
            IEEEBinary32.MIN_SUB_EXPONENT, BigInt(IEEEBinary32.SUP_TRAILING));
        return IEEEBinary32.fromParts(s, e, Number(t)).number;
    }

    // Member access.
    //
    // A block of calls involving setSign, setBiasedExponent and setUInteger
    // implies a deferred creation of a number; the intermediate states can be
    // invalid, which is why these setters throw no exceptions.
    setSign(signValue: number): void {
        this.mSign = signValue;
    }

    getSign(): number {
        return this.mSign;
    }

    // In-place negation (upstream Negate).
    negate(): void {
        this.mSign = -this.mSign;
    }

    setBiasedExponent(biasedExponent: number): void {
        this.mBiasedExponent = biasedExponent;
    }

    getBiasedExponent(): number {
        return this.mBiasedExponent;
    }

    setExponent(exponent: number): void {
        this.mBiasedExponent = exponent - numBitsOf(this.mUInteger) + 1;
    }

    getExponent(): number {
        return this.mBiasedExponent + numBitsOf(this.mUInteger) - 1;
    }

    getUInteger(): bigint {
        return this.mUInteger;
    }

    setUInteger(uinteger: bigint): void {
        this.mUInteger = uinteger;
    }

    // The port of UInteger::GetNumBits() for the stored unsigned integer.
    getNumBits(): number {
        return numBitsOf(this.mUInteger);
    }

    // The invariant for a nonzero BSNumber is that the unsigned integer part
    // is a positive odd number; for zero, all members are zero. Upstream
    // exposes this only when GTE_VALIDATE_BSNUMBER is defined.
    isValid(): boolean {
        if (this.mSign !== 0) {
            return this.mUInteger > 0n && (this.mUInteger & 1n) === 1n;
        }
        return this.mBiasedExponent === 0 && this.mUInteger === 0n;
    }

    // Comparisons.
    equals(other: BSNumber): boolean {
        return (this.mSign === other.mSign ? BSNumber.equalIgnoreSign(this, other) : false);
    }

    notEquals(other: BSNumber): boolean {
        return !this.equals(other);
    }

    lessThan(other: BSNumber): boolean {
        if (this.mSign > 0) {
            if (other.mSign <= 0) {
                return false;
            }
            // Both numbers are positive.
            return BSNumber.lessThanIgnoreSign(this, other);
        } else if (this.mSign < 0) {
            if (other.mSign >= 0) {
                return true;
            }
            // Both numbers are negative.
            return BSNumber.lessThanIgnoreSign(other, this);
        } else {
            return other.mSign > 0;
        }
    }

    lessThanOrEqual(other: BSNumber): boolean {
        return !other.lessThan(this);
    }

    greaterThan(other: BSNumber): boolean {
        return other.lessThan(this);
    }

    greaterThanOrEqual(other: BSNumber): boolean {
        return !this.lessThan(other);
    }

    // Unary operations. The upstream unary operator+ is the identity and is
    // omitted; unary operator- is negated().
    negated(): BSNumber {
        const result = this.clone();
        result.mSign = -result.mSign;
        return result;
    }

    // Arithmetic. The results are exact.
    add(n1: BSNumber): BSNumber {
        const n0: BSNumber = this;

        if (n0.mSign === 0) {
            return n1.clone();
        }

        if (n1.mSign === 0) {
            return n0.clone();
        }

        if (n0.mSign > 0) {
            if (n1.mSign > 0) {
                // n0 + n1 = |n0| + |n1|
                return BSNumber.addIgnoreSign(n0, n1, +1);
            } else { // n1.mSign < 0
                if (!BSNumber.equalIgnoreSign(n0, n1)) {
                    if (BSNumber.lessThanIgnoreSign(n1, n0)) {
                        // n0 + n1 = |n0| - |n1| > 0
                        return BSNumber.subIgnoreSign(n0, n1, +1);
                    } else {
                        // n0 + n1 = -(|n1| - |n0|) < 0
                        return BSNumber.subIgnoreSign(n1, n0, -1);
                    }
                }
                // else n0 + n1 = 0
            }
        } else { // n0.mSign < 0
            if (n1.mSign < 0) {
                // n0 + n1 = -(|n0| + |n1|)
                return BSNumber.addIgnoreSign(n0, n1, -1);
            } else { // n1.mSign > 0
                if (!BSNumber.equalIgnoreSign(n0, n1)) {
                    if (BSNumber.lessThanIgnoreSign(n1, n0)) {
                        // n0 + n1 = -(|n0| - |n1|) < 0
                        return BSNumber.subIgnoreSign(n0, n1, -1);
                    } else {
                        // n0 + n1 = |n1| - |n0| > 0
                        return BSNumber.subIgnoreSign(n1, n0, +1);
                    }
                }
                // else n0 + n1 = 0
            }
        }

        return new BSNumber();  // = 0
    }

    sub(n1: BSNumber): BSNumber {
        const n0: BSNumber = this;

        if (n0.mSign === 0) {
            return n1.negated();
        }

        if (n1.mSign === 0) {
            return n0.clone();
        }

        if (n0.mSign > 0) {
            if (n1.mSign < 0) {
                // n0 - n1 = |n0| + |n1|
                return BSNumber.addIgnoreSign(n0, n1, +1);
            } else { // n1.mSign > 0
                if (!BSNumber.equalIgnoreSign(n0, n1)) {
                    if (BSNumber.lessThanIgnoreSign(n1, n0)) {
                        // n0 - n1 = |n0| - |n1| > 0
                        return BSNumber.subIgnoreSign(n0, n1, +1);
                    } else {
                        // n0 - n1 = -(|n1| - |n0|) < 0
                        return BSNumber.subIgnoreSign(n1, n0, -1);
                    }
                }
                // else n0 - n1 = 0
            }
        } else { // n0.mSign < 0
            if (n1.mSign > 0) {
                // n0 - n1 = -(|n0| + |n1|)
                return BSNumber.addIgnoreSign(n0, n1, -1);
            } else { // n1.mSign < 0
                if (!BSNumber.equalIgnoreSign(n0, n1)) {
                    if (BSNumber.lessThanIgnoreSign(n1, n0)) {
                        // n0 - n1 = -(|n0| - |n1|) < 0
                        return BSNumber.subIgnoreSign(n0, n1, -1);
                    } else {
                        // n0 - n1 = |n1| - |n0| > 0
                        return BSNumber.subIgnoreSign(n1, n0, +1);
                    }
                }
                // else n0 - n1 = 0
            }
        }

        return new BSNumber();  // = 0
    }

    mul(other: BSNumber): BSNumber {
        const result = new BSNumber();  // = 0
        const resultSign = this.mSign * other.mSign;
        if (resultSign !== 0) {
            result.mSign = resultSign;
            result.mBiasedExponent = this.mBiasedExponent + other.mBiasedExponent;
            // The product of two odd numbers is odd, so no normalization is
            // needed here.
            result.mUInteger = this.mUInteger * other.mUInteger;
        }
        return result;
    }

    // Helper for converting a string to a BSNumber. The string must be valid
    // for a nonnegative integer without a leading '+' sign.
    private static convertToInteger(numberString: string): BSNumber {
        let digit = numberString.charCodeAt(numberString.length - 1) - 0x30;
        logAssert(0 <= digit && digit <= 9, 'Invalid number format.');
        let x = BSNumber.fromNumber(digit);
        if (numberString.length > 1) {
            logAssert(/^[1-9]/.test(numberString), 'Invalid number format.');
            logAssert(/^[0-9]+$/.test(numberString), 'Invalid number format.');
            const ten = BSNumber.fromNumber(10);
            let pow10 = BSNumber.fromNumber(10);
            for (let i = 1, j = numberString.length - 2; i < numberString.length; ++i, --j) {
                digit = numberString.charCodeAt(j) - 0x30;
                if (digit > 0) {
                    x = x.add(BSNumber.fromNumber(digit).mul(pow10));
                }
                pow10 = pow10.mul(ten);
            }
        }
        return x;
    }

    // Helpers for equals, lessThan, add, sub.
    private static equalIgnoreSign(n0: BSNumber, n1: BSNumber): boolean {
        return n0.mBiasedExponent === n1.mBiasedExponent && n0.mUInteger === n1.mUInteger;
    }

    private static lessThanIgnoreSign(n0: BSNumber, n1: BSNumber): boolean {
        const e0 = n0.getExponent(), e1 = n1.getExponent();
        if (e0 < e1) {
            return true;
        }
        if (e0 > e1) {
            return false;
        }
        // The port of UInteger::operator<, which is NOT the integer
        // comparison of the stored magnitudes. We get here with
        // n0 = 1.u * 2^p and n1 = 1.v * 2^p, so the unsigned integers must be
        // compared as if they are left-aligned with each other; although the
        // numbers have the same exponent, it is possible that n0 < n1 while
        // numBits(1u) > numBits(1v). Upstream aligns the leading 1-bits and
        // compares one 32-bit block at a time, which is the same as shifting
        // both magnitudes to a common bit length and comparing them.
        const b0 = numBitsOf(n0.mUInteger);
        const b1 = numBitsOf(n1.mUInteger);
        const maxBits = Math.max(b0, b1);
        return (n0.mUInteger << BigInt(maxBits - b0))
            < (n1.mUInteger << BigInt(maxBits - b1));
    }

    // Add two positive numbers.
    private static addIgnoreSign(n0: BSNumber, n1: BSNumber, resultSign: number): BSNumber {
        const result = new BSNumber();

        const diff = n0.mBiasedExponent - n1.mBiasedExponent;
        if (diff > 0) {
            result.mUInteger = (n0.mUInteger << BigInt(diff)) + n1.mUInteger;
            result.mBiasedExponent = n1.mBiasedExponent;
        } else if (diff < 0) {
            result.mUInteger = n0.mUInteger + (n1.mUInteger << BigInt(-diff));
            result.mBiasedExponent = n0.mBiasedExponent;
        } else {
            // Both integers are odd, so the sum is even; normalize it.
            const { odd, shift } = shiftRightToOdd(n0.mUInteger + n1.mUInteger);
            result.mUInteger = odd;
            result.mBiasedExponent = n0.mBiasedExponent + shift;
        }

        result.mSign = resultSign;
        return result;
    }

    // Subtract two positive numbers where n0 > n1.
    private static subIgnoreSign(n0: BSNumber, n1: BSNumber, resultSign: number): BSNumber {
        const result = new BSNumber();

        const diff = n0.mBiasedExponent - n1.mBiasedExponent;
        if (diff > 0) {
            result.mUInteger = (n0.mUInteger << BigInt(diff)) - n1.mUInteger;
            result.mBiasedExponent = n1.mBiasedExponent;
        } else if (diff < 0) {
            result.mUInteger = n0.mUInteger - (n1.mUInteger << BigInt(-diff));
            result.mBiasedExponent = n0.mBiasedExponent;
        } else {
            // Both integers are odd, so the difference is even; normalize it.
            const { odd, shift } = shiftRightToOdd(n0.mUInteger - n1.mUInteger);
            result.mUInteger = odd;
            result.mBiasedExponent = n0.mBiasedExponent + shift;
        }

        result.mSign = resultSign;
        return result;
    }

    // Support for conversions from floating-point numbers to BSNumber. The
    // IEEE-specific constants are passed in rather than selected by a C++
    // template parameter.
    private convertFrom(s: number, e: number, t: bigint, numTrailingBits: number,
        exponentBias: number, maxBiasedExponent: number, minSubExponent: number,
        supTrailing: bigint): void {
        if (e === 0) {
            if (t === 0n) {
                // x = (-1)^s * 0
                this.mSign = 0;
                this.mBiasedExponent = 0;
                this.mUInteger = 0n;
            } else {
                // Subnormal numbers: x = (-1)^s * 0.t * 2^{1-EXPONENT_BIAS}
                const last = trailingZerosOf(t);
                const diff = numTrailingBits - last;
                this.mSign = (s > 0 ? -1 : 1);
                this.mBiasedExponent = minSubExponent - diff;
                this.mUInteger = t >> BigInt(last);
            }
        } else if (e < maxBiasedExponent) {
            // Normal numbers: x = (-1)^s * 1.t * 2^{e-EXPONENT_BIAS}
            if (t > 0n) {
                const last = trailingZerosOf(t);
                const diff = numTrailingBits - last;
                this.mSign = (s > 0 ? -1 : 1);
                this.mBiasedExponent = e - exponentBias - diff;
                this.mUInteger = (t | supTrailing) >> BigInt(last);
            } else {
                this.mSign = (s > 0 ? -1 : 1);
                this.mBiasedExponent = e - exponentBias;
                this.mUInteger = 1n;
            }
        } else {
            // e == MAX_BIASED_EXPONENT, special numbers.
            if (t === 0n) {
                // Infinities. BSNumber has no representation for these, so
                // return (-1)^s * 2^{1+EXPONENT_BIAS} for a graceful exit
                // (the upstream default; GTE_THROW_ON_CONVERT_FROM_INFINITY_
                // OR_NAN would throw instead).
                this.mSign = (s > 0 ? -1 : 1);
                this.mBiasedExponent = 1 + exponentBias;
                this.mUInteger = 1n;
            } else {
                // Not-a-number (NaN). Return 0 for a graceful exit.
                this.mSign = 0;
                this.mBiasedExponent = 0;
                this.mUInteger = 0n;
            }
        }
    }

    // Support for conversions from BSNumber to floating-point numbers. The
    // returned (sign, biased exponent, trailing significand) triple is fed to
    // the corresponding IEEEBinary* encoder.
    private convertTo(numSignificandBits: number, exponentBias: number,
        maxBiasedExponent: number, minExponent: number, minSubExponent: number,
        supTrailing: bigint): { s: number, e: number, t: bigint } {
        const s = (this.mSign < 0 ? 1 : 0);
        let e: number;
        let t: bigint;

        if (this.mSign !== 0) {
            // The conversions use round-to-nearest-ties-to-even semantics.
            const exponent = this.getExponent();
            if (exponent < minExponent) {
                if (exponent < minExponent - 1 || numBitsOf(this.mUInteger) === 1) {
                    // x = 1.0*2^{MIN_EXPONENT-1}, round to zero.
                    e = 0;
                    t = 0n;
                } else {
                    // Round to min subnormal.
                    e = 0;
                    t = 1n;
                }
            } else if (exponent < minSubExponent) {
                // The second input is in {0, ..., NUM_TRAILING_BITS-1}.
                t = this.getTrailing(numSignificandBits, 0,
                    minSubExponent - exponent - 1);
                if (t & supTrailing) {
                    // Leading NUM_SIGNIFICAND_BITS bits were all 1, so round
                    // to min normal.
                    e = 1;
                    t = 0n;
                } else {
                    e = 0;
                }
            } else if (exponent <= exponentBias) {
                e = exponent + exponentBias;
                t = this.getTrailing(numSignificandBits, 1, 0);
                if (t & (supTrailing << 1n)) {
                    // Carry-out occurred, so increase the exponent by 1 and
                    // shift right to compensate.
                    ++e;
                    t >>= 1n;
                }
                // Eliminate the leading 1 (implied for normals).
                t &= ~supTrailing;
            } else {
                // Set to infinity.
                e = maxBiasedExponent;
                t = 0n;
            }
        } else {
            // The input is zero.
            e = 0;
            t = 0n;
        }

        return { s, e, t };
    }

    private getTrailing(numSignificandBits: number, normal: number, sigma: number): bigint {
        const numRequested = numSignificandBits + normal;

        // We need numRequested bits to determine the rounding direction.
        // These are stored in the high-order bits of 'prefix'.
        const prefix = this.getPrefix();

        // The first bit index after the implied binary point for rounding.
        const diff = numRequested - sigma;
        const roundBitIndex = 64 - diff;

        // Determine the value based on round-to-nearest-ties-to-even.
        const mask = 1n << BigInt(roundBitIndex);
        let round: bigint;
        if (prefix & mask) {
            // The first bit of the remainder is 1.
            if (numBitsOf(this.mUInteger) === diff) {
                // The first bit of the remainder is the lowest-order bit of
                // the unsigned integer. Apply the ties-to-even rule.
                if (prefix & (mask << 1n)) {
                    // The last bit of the trailing significand is odd, so
                    // round up.
                    round = 1n;
                } else {
                    // The last bit of the trailing significand is even, so
                    // round down.
                    round = 0n;
                }
            } else {
                // The first bit of the remainder is not the lowest-order bit
                // of the unsigned integer. The remainder as a fraction is
                // larger than 1/2, so round up.
                round = 1n;
            }
        } else {
            // The first bit of the remainder is 0, so round down.
            round = 0n;
        }

        // Get the unrounded trailing significand and apply the rounding.
        const trailing = prefix >> BigInt(roundBitIndex + 1);
        return trailing + round;
    }

    // The port of UInteger::GetPrefix(numRequested): the leading bits of the
    // nonzero unsigned integer, with the leading 1-bit placed at bit 63 and
    // the low-order bits zero-filled when fewer than 64 bits are available.
    // Upstream stops copying 32-bit blocks once 'numRequested' bits have been
    // consumed, which zeroes the remaining low-order bits of the result; the
    // callers never read those bits (they read at most the leading
    // numRequested + 1 bits), so the parameter is dropped here.
    private getPrefix(): bigint {
        const numBits = numBitsOf(this.mUInteger);
        if (numBits >= 64) {
            return this.mUInteger >> BigInt(numBits - 64);
        }
        return BigInt.asUintN(64, this.mUInteger << BigInt(64 - numBits));
    }

    // The upstream overloads of the std:: math functions and of the gte::
    // helper functions for BSNumber. Except for fabs, frexp and ldexp, they
    // convert to double, call the double version and convert back, so the
    // results are not exact. They are static methods here to avoid colliding
    // with the module-level 'number' versions in Functions.ts.
    static acos(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.acos(x.toNumber()));
    }

    static acosh(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.acosh(x.toNumber()));
    }

    static asin(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.asin(x.toNumber()));
    }

    static asinh(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.asinh(x.toNumber()));
    }

    static atan(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.atan(x.toNumber()));
    }

    static atanh(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.atanh(x.toNumber()));
    }

    static atan2(y: BSNumber, x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.atan2(y.toNumber(), x.toNumber()));
    }

    static ceil(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.ceil(x.toNumber()));
    }

    static cos(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.cos(x.toNumber()));
    }

    static cosh(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.cosh(x.toNumber()));
    }

    static exp(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.exp(x.toNumber()));
    }

    static exp2(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.pow(2, x.toNumber()));
    }

    // Exact: the absolute value only clears the sign.
    static fabs(x: BSNumber): BSNumber {
        return (x.getSign() >= 0 ? x.clone() : x.negated());
    }

    static floor(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.floor(x.toNumber()));
    }

    static fmod(x: BSNumber, y: BSNumber): BSNumber {
        return BSNumber.fromNumber(x.toNumber() % y.toNumber());
    }

    // Exact: x = result * 2^exponent with result in [1/2, 1). The C++
    // 'int32_t* exponent' output parameter becomes a field of the returned
    // object.
    static frexp(x: BSNumber): { result: BSNumber, exponent: number } {
        if (x.getSign() !== 0) {
            const result = x.clone();
            const exponent = result.getExponent() + 1;
            result.setExponent(-1);
            return { result, exponent };
        }
        return { result: new BSNumber(), exponent: 0 };
    }

    // Exact: multiplication by a power of two only shifts the exponent.
    static ldexp(x: BSNumber, exponent: number): BSNumber {
        const result = x.clone();
        result.setBiasedExponent(result.getBiasedExponent() + exponent);
        return result;
    }

    static log(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.log(x.toNumber()));
    }

    static log2(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.log2(x.toNumber()));
    }

    static log10(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.log10(x.toNumber()));
    }

    static pow(x: BSNumber, y: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.pow(x.toNumber(), y.toNumber()));
    }

    // The IEEE remainder: x - n*y where n is the integer nearest x/y, ties
    // resolved to the even n (the port of std::remainder).
    static remainder(x: BSNumber, y: BSNumber): BSNumber {
        const dx = x.toNumber();
        const dy = y.toNumber();
        const ratio = dx / dy;
        let n = Math.round(ratio);
        if (Math.abs(ratio - Math.trunc(ratio)) === 0.5 && (n % 2 !== 0)) {
            // Math.round breaks ties upward; std::remainder breaks them to
            // the even integer.
            n -= 1;
        }
        return BSNumber.fromNumber(dx - n * dy);
    }

    static sin(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.sin(x.toNumber()));
    }

    static sinh(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.sinh(x.toNumber()));
    }

    static sqrt(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.sqrt(x.toNumber()));
    }

    static tan(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.tan(x.toNumber()));
    }

    static tanh(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(Math.tanh(x.toNumber()));
    }

    static atandivpi(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(atandivpi(x.toNumber()));
    }

    static atan2divpi(y: BSNumber, x: BSNumber): BSNumber {
        return BSNumber.fromNumber(atan2divpi(y.toNumber(), x.toNumber()));
    }

    static clamp(x: BSNumber, xmin: BSNumber, xmax: BSNumber): BSNumber {
        return BSNumber.fromNumber(clamp(x.toNumber(), xmin.toNumber(), xmax.toNumber()));
    }

    static cospi(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(cospi(x.toNumber()));
    }

    static exp10(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(exp10(x.toNumber()));
    }

    static invsqrt(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(invsqrt(x.toNumber()));
    }

    static isign(x: BSNumber): number {
        return isign(x.toNumber());
    }

    static saturate(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(saturate(x.toNumber()));
    }

    static sign(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(sign(x.toNumber()));
    }

    static sinpi(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(sinpi(x.toNumber()));
    }

    static sqr(x: BSNumber): BSNumber {
        return BSNumber.fromNumber(sqr(x.toNumber()));
    }

    // Compute u * v + w exactly.
    static fma(u: BSNumber, v: BSNumber, w: BSNumber): BSNumber {
        return u.mul(v).add(w);
    }

    // Sum of products (SOP) u * v + w * z, computed exactly.
    static robustSOP(u: BSNumber, v: BSNumber, w: BSNumber, z: BSNumber): BSNumber {
        return u.mul(v).add(w.mul(z));
    }

    // Difference of products (DOP) u * v - w * z, computed exactly.
    static robustDOP(u: BSNumber, v: BSNumber, w: BSNumber, z: BSNumber): BSNumber {
        return u.mul(v).sub(w.mul(z));
    }
}

// Explicit conversion to a user-specified precision (upstream Convert). The
// output is the input rounded to 'precision' bits using the specified
// rounding mode. Upstream also rejects a precision that exceeds the storage
// capacity of a fixed-precision UInteger; the bigint storage of this port has
// no such limit, so only the positivity of 'precision' is checked.
export function convertBSNumber(input: BSNumber, precision: number,
    roundingMode: BSNumberRoundingMode): BSNumber {
    if (precision <= 0) {
        logError('Precision must be positive.');
    }

    if (input.getSign() === 0) {
        return new BSNumber();
    }

    // Let p = precision and n+1 be the number of bits of the input. Compute
    // n+1-p. If it is nonpositive, then the requested precision is already
    // satisfied by the input.
    const inW = input.getUInteger();
    const np1mp = input.getNumBits() - precision;
    if (np1mp <= 0) {
        return input.clone();
    }

    // At this point, the requested number of bits is smaller than the number
    // of bits of the input. Round the input to the smaller number of bits
    // using the specified rounding mode. The leading 'precision' bits are
    // u_n ... u_{n-p+1}; 'lastBit' is u_{n-p+1} and the remainder starts with
    // the bit u_{n-p}.
    let outW = inW >> BigInt(np1mp);
    const lastBit = Number(outW & 1n);
    const positive = Number((inW >> BigInt(np1mp - 1)) & 1n);

    const signValue = input.getSign();
    let outExponent = input.getExponent();
    const precisionM1 = precision - 1;

    // The port of 'outExponent += outW.RoundUp()', which adds one to the
    // integer and normalizes it back to an odd number.
    const roundUp = (): void => {
        const { odd, shift } = shiftRightToOdd(outW + 1n);
        outW = odd;
        outExponent += shift;
    };

    if (roundingMode === BSNumberRoundingMode.FE_TONEAREST) {
        // Determine whether u_{n-p} is positive.
        if (positive !== 0 && (np1mp > 1 || lastBit === 1)) {
            roundUp();
        }
        // else round down, equivalent to truncating the r bits
    } else if (roundingMode === BSNumberRoundingMode.FE_UPWARD) {
        // The remainder r must be positive because n-p >= 0 and u_0 = 1.
        if (signValue > 0) {
            roundUp();
        }
        // else round down, equivalent to truncating the r bits
    } else if (roundingMode === BSNumberRoundingMode.FE_DOWNWARD) {
        // The remainder r must be positive because n-p >= 0 and u_0 = 1.
        if (signValue < 0) {
            // Round down. This is the round-up operation applied to w, but
            // the final sign is negative which amounts to rounding down.
            roundUp();
        }
        // else round down, equivalent to truncating the r bits
    } else if (roundingMode !== BSNumberRoundingMode.FE_TOWARDZERO) {
        // Currently, no additional implementation-dependent modes are
        // supported for rounding.
        logError('Implementation-dependent rounding mode not supported.');
    }
    // else roundingMode == FE_TOWARDZERO. Truncate the r bits, which requires
    // no additional work.

    // Shift the bits if necessary to obtain the invariant that BSNumber
    // objects have bit patterns that are odd integers.
    if (outW > 0n && (outW & 1n) === 0n) {
        const { odd, shift } = shiftRightToOdd(outW);
        outW = odd;
        outExponent += shift;
    }

    // Do not use setExponent(outExponent) at this step. The number of
    // requested bits is 'precision' but the number of bits of outW will be
    // different when round-up occurs, and setExponent uses that count.
    return BSNumber.fromParts(signValue, outExponent - precisionM1, outW);
}
