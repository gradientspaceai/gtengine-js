// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BSRational.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The class BSRational (binary scientific rational) is the exact rational
// number type built from a pair of BSNumber objects,
//   x = numerator / denominator
// where the denominator is chosen to be positive, which allows some
// simplification of the comparisons. Addition, subtraction, multiplication
// and division are all exact; no common factors are removed, so the number of
// bits grows as needed.
//
// The document at the following link describes the design, implementation and
// use of BSNumber and BSRational.
//   https://www.geometrictools.com/Documentation/ArbitraryPrecision.pdf
//
// Port notes:
// - Upstream is 'template <typename UInteger> class BSRational'. As with the
//   BSNumber port, the UInteger template parameter disappears because the
//   unsigned integers are JavaScript 'bigint' values, so BSRational is always
//   arbitrary precision.
// - The C++ constructor overloads (float, double, int32_t, uint32_t, int64_t,
//   uint64_t, BSNumber, each with an optional matching denominator, plus
//   std::string) become the static factories fromNumber, fromFloat32,
//   fromBigInt, fromBSNumber and fromString, each taking an optional
//   denominator argument. The default constructor generates zero (0/1).
// - Operators become named methods, mirroring BSNumber.ts: ==, <, ... ->
//   equals, lessThan, ...; +, -, *, / -> add, sub, mul, div; unary - ->
//   negated() (a new object), distinct from the in-place Negate() ->
//   negate(). operator float / operator double -> toFloat32() / toNumber().
// - Upstream has four Convert overloads. They become
//     convertBSRationalToBSNumber(input, precision, mode) -> BSNumber
//     convertBSRational(input, precision, mode)           -> BSRational
//     convertBSRationalToNumber(input, mode)              -> number
//     convertBSRationalToFloat32(input, mode)             -> number
//   using the BSNumberRoundingMode enum of BSNumber.ts in place of the
//   <cfenv> flags. Upstream also rejects a precision that exceeds the storage
//   capacity of a fixed-precision UInteger; the bigint storage of this port
//   has no such limit, so only the positivity of 'precision' is checked.
// - The 'std::' math functions (acos, sqrt, floor, ...) and the gte:: helper
//   functions (clamp, sign, sqr, ...) that upstream overloads for BSRational
//   are static methods of BSRational, exactly as the BSNumber port does them,
//   because the module-level names are taken by the 'number' versions in
//   Functions.ts (see the global export uniqueness rule in PORTING.md).
// - Write/Read serialize the C++ block representation to a binary stream; the
//   bigint storage has no such layout, so they are omitted (as in BSNumber).

import { logAssert, logError } from './Logger.js';
import { BSNumber, BSNumberRoundingMode } from './BSNumber.js';
import { IEEEBinary32, IEEEBinary64 } from './IEEEBinary.js';
import type { ArbitraryPrecisionNumber } from './TypeTraits.js';
import {
    atandivpi, atan2divpi, clamp, cospi, exp10, invsqrt, isign, saturate,
    sign, sinpi, sqr
} from './Functions.js';

// The port of 'shift = result.ShiftRightToOdd(number)': the odd number
// obtained by removing the trailing zero bits of a positive bigint, together
// with the number of bits removed. BSNumber.ts keeps its copy of this helper
// module-private, so it is repeated here; it is used only by
// convertBSRationalToBSNumber, which ports UInteger::RoundUp() as
// shiftRightToOdd(w + 1n).
function shiftRightToOdd(u: bigint): { odd: bigint, shift: number } {
    if (u === 0n) {
        return { odd: 0n, shift: 0 };
    }
    let shift = 0;
    let v = u;
    // Strip 32 bits at a time, then the remaining bits one at a time.
    while ((v & 0xFFFFFFFFn) === 0n) {
        v >>= 32n;
        shift += 32;
    }
    let block = Number(v & 0xFFFFFFFFn);
    while ((block & 1) === 0) {
        block >>>= 1;
        ++shift;
    }
    return { odd: u >> BigInt(shift), shift };
}

export class BSRational implements ArbitraryPrecisionNumber {
    // See TypeTraits.ts. BSRational is arbitrary precision and, unlike
    // BSNumber, has a division operator.
    readonly isArbitraryPrecision = true as const;
    readonly hasDivisionOperator = true;

    private mNumerator: BSNumber;
    private mDenominator: BSNumber;

    // The default constructor generates the zero BSRational, 0/1.
    constructor() {
        this.mNumerator = new BSNumber();
        this.mDenominator = BSNumber.fromNumber(1);
    }

    // The port of BSRational(BSNumber) and BSRational(BSNumber, BSNumber).
    // The inputs are copied (C++ value semantics).
    static fromBSNumber(numerator: BSNumber, denominator?: BSNumber): BSRational {
        const result = new BSRational();
        result.mNumerator = numerator.clone();
        if (denominator === undefined) {
            result.mDenominator = BSNumber.fromNumber(1);
            return result;
        }

        result.mDenominator = denominator.clone();
        logAssert(result.mDenominator.getSign() !== 0, 'Division by zero.');
        if (result.mDenominator.getSign() < 0) {
            result.mNumerator.setSign(-result.mNumerator.getSign());
            result.mDenominator.setSign(1);
        }

        // Set the exponent of the denominator to zero, but you can do so only
        // by modifying the biased exponent. Adjust the numerator accordingly.
        // This prevents large growth of the exponents in both numerator and
        // denominator simultaneously.
        //
        // Upstream adjusts the numerator's biased exponent unconditionally.
        // When the numerator is zero, that produces a BSNumber with sign 0
        // and a nonzero biased exponent, which violates the BSNumber
        // invariant (the value is still zero, so no arithmetic is corrupted).
        // The port skips the adjustment for a zero numerator, matching the
        // guard the BSNumber port applies in fromString.
        if (result.mNumerator.getSign() !== 0) {
            result.mNumerator.setBiasedExponent(
                result.mNumerator.getBiasedExponent() - result.mDenominator.getExponent());
        }
        result.mDenominator.setBiasedExponent(-(result.mDenominator.getNumBits() - 1));
        return result;
    }

    // Conversion from 'number' values using the IEEE binary64 representation,
    // which is exact. This is the port of BSRational(double) and
    // BSRational(double, double) and also of the int32_t/uint32_t
    // constructors, whose inputs are exactly representable as doubles.
    static fromNumber(numerator: number, denominator?: number): BSRational {
        return BSRational.fromBSNumber(BSNumber.fromNumber(numerator),
            denominator === undefined ? undefined : BSNumber.fromNumber(denominator));
    }

    // Conversion from 'number' values rounded to the IEEE binary32
    // representation; the port of BSRational(float) and
    // BSRational(float, float).
    static fromFloat32(numerator: number, denominator?: number): BSRational {
        return BSRational.fromBSNumber(BSNumber.fromFloat32(numerator),
            denominator === undefined ? undefined : BSNumber.fromFloat32(denominator));
    }

    // Exact conversion from integers of any size; the port of the
    // int64_t/uint64_t constructors.
    static fromBigInt(numerator: bigint, denominator?: bigint): BSRational {
        return BSRational.fromBSNumber(BSNumber.fromBigInt(numerator),
            denominator === undefined ? undefined : BSNumber.fromBigInt(denominator));
    }

    // The number must be of the form "x", "x.", "x.y" or ".y", optionally
    // preceded by '+' or '-', where x and y are strings of decimal digits.
    static fromString(numberString: string): BSRational {
        logAssert(numberString.length > 0, 'A number must be specified.');

        // Get the leading '+' or '-' if it exists.
        let fpNumber: string;
        let signValue: number;
        if (numberString[0] === '+') {
            fpNumber = numberString.substring(1);
            signValue = +1;
            logAssert(fpNumber.length > 0, 'Invalid number format.');
        } else if (numberString[0] === '-') {
            fpNumber = numberString.substring(1);
            signValue = -1;
            logAssert(fpNumber.length > 0, 'Invalid number format.');
        } else {
            fpNumber = numberString;
            signValue = +1;
        }

        const result = new BSRational();
        const decimal = fpNumber.indexOf('.');
        if (decimal >= 0) {
            if (decimal > 0) {
                // The number is "x.y" (with y possibly empty, in which case
                // the fractional part is zero and this handles "x." as well;
                // upstream has a separate "x." branch that is unreachable
                // because the index of '.' is always less than the length of
                // the string).
                const intPart = BSNumber.fromString(fpNumber.substring(0, decimal));
                const frcPart = BSRational.convertToFraction(fpNumber.substring(decimal + 1));
                result.mNumerator = intPart.mul(frcPart.mDenominator).add(frcPart.mNumerator);
                result.mDenominator = frcPart.mDenominator;
            } else {
                // The number is ".y".
                const frcPart = BSRational.convertToFraction(fpNumber.substring(1));
                result.mNumerator = frcPart.mNumerator;
                result.mDenominator = frcPart.mDenominator;
            }
        } else {
            // The number is "x".
            result.mNumerator = BSNumber.fromString(fpNumber);
            result.mDenominator = BSNumber.fromNumber(1);
        }

        // Upstream calls mNumerator.SetSign(sign) unconditionally, which for
        // a zero numerator (for example, "-0.0") produces a BSNumber with
        // sign -1 and a zero unsigned integer; that violates the invariant
        // and makes the number compare unequal to zero. The port keeps the
        // zero, matching the guard the BSNumber port applies in fromString.
        if (result.mNumerator.getSign() !== 0) {
            result.mNumerator.setSign(signValue);
        }
        return result;
    }

    // C++ copy assignment; TypeScript objects alias, so clone explicitly.
    clone(): BSRational {
        const result = new BSRational();
        result.mNumerator = this.mNumerator.clone();
        result.mDenominator = this.mDenominator.clone();
        return result;
    }

    // Conversions to floating-point. These always use the default rounding
    // mode, round-to-nearest-ties-to-even.
    toNumber(): number {
        return convertBSRationalToNumber(this, BSNumberRoundingMode.FE_TONEAREST);
    }

    toFloat32(): number {
        return convertBSRationalToFloat32(this, BSNumberRoundingMode.FE_TONEAREST);
    }

    // Member access.
    setSign(signValue: number): void {
        this.mNumerator.setSign(signValue);
        this.mDenominator.setSign(1);
    }

    getSign(): number {
        return this.mNumerator.getSign() * this.mDenominator.getSign();
    }

    // In-place negation (upstream Negate).
    negate(): void {
        this.mNumerator.negate();
    }

    // The returned objects alias the internal storage, as the non-const
    // upstream accessors do; write through them only when you intend to
    // modify this number (see BSRational.frexp and BSRational.ldexp).
    getNumerator(): BSNumber {
        return this.mNumerator;
    }

    getDenominator(): BSNumber {
        return this.mDenominator;
    }

    // Comparisons. The denominators are positive, so the cross-multiplied
    // comparisons need no sign adjustments.
    equals(r: BSRational): boolean {
        // Do inexpensive sign tests first for optimum performance.
        if (this.mNumerator.getSign() !== r.mNumerator.getSign()) {
            return false;
        }
        if (this.mNumerator.getSign() === 0) {
            // The numbers are both zero.
            return true;
        }

        return this.mNumerator.mul(r.mDenominator).equals(this.mDenominator.mul(r.mNumerator));
    }

    notEquals(r: BSRational): boolean {
        return !this.equals(r);
    }

    lessThan(r: BSRational): boolean {
        // Do inexpensive sign tests first for optimum performance.
        const s0 = this.mNumerator.getSign();
        const s1 = r.mNumerator.getSign();
        if (s0 > 0) {
            if (s1 <= 0) {
                return false;
            }
        } else if (s0 === 0) {
            return s1 > 0;
        } else { // s0 < 0
            if (s1 >= 0) {
                return true;
            }
        }

        return this.mNumerator.mul(r.mDenominator).lessThan(this.mDenominator.mul(r.mNumerator));
    }

    lessThanOrEqual(r: BSRational): boolean {
        return !r.lessThan(this);
    }

    greaterThan(r: BSRational): boolean {
        return r.lessThan(this);
    }

    greaterThanOrEqual(r: BSRational): boolean {
        return !this.lessThan(r);
    }

    // Unary operations. The upstream unary operator+ is the identity and is
    // omitted; unary operator- is negated().
    negated(): BSRational {
        return BSRational.fromBSNumber(this.mNumerator.negated(), this.mDenominator);
    }

    // Arithmetic. The results are exact.
    add(r: BSRational): BSRational {
        const product0 = this.mNumerator.mul(r.mDenominator);
        const product1 = this.mDenominator.mul(r.mNumerator);
        const numerator = product0.add(product1);

        // Complex expressions can lead to 0/denom, where denom is not 1.
        if (numerator.getSign() !== 0) {
            const denominator = this.mDenominator.mul(r.mDenominator);
            return BSRational.fromBSNumber(numerator, denominator);
        }
        return new BSRational();
    }

    sub(r: BSRational): BSRational {
        const product0 = this.mNumerator.mul(r.mDenominator);
        const product1 = this.mDenominator.mul(r.mNumerator);
        const numerator = product0.sub(product1);

        // Complex expressions can lead to 0/denom, where denom is not 1.
        if (numerator.getSign() !== 0) {
            const denominator = this.mDenominator.mul(r.mDenominator);
            return BSRational.fromBSNumber(numerator, denominator);
        }
        return new BSRational();
    }

    mul(r: BSRational): BSRational {
        const numerator = this.mNumerator.mul(r.mNumerator);

        // Complex expressions can lead to 0/denom, where denom is not 1.
        if (numerator.getSign() !== 0) {
            const denominator = this.mDenominator.mul(r.mDenominator);
            return BSRational.fromBSNumber(numerator, denominator);
        }
        return new BSRational();
    }

    div(r: BSRational): BSRational {
        logAssert(r.mNumerator.getSign() !== 0, 'Division by zero in BSRational::operator/.');

        const numerator = this.mNumerator.mul(r.mDenominator);

        // Complex expressions can lead to 0/denom, where denom is not 1.
        if (numerator.getSign() !== 0) {
            const denominator = this.mDenominator.mul(r.mNumerator);
            if (denominator.getSign() < 0) {
                numerator.setSign(-numerator.getSign());
                denominator.setSign(1);
            }
            return BSRational.fromBSNumber(numerator, denominator);
        }
        return new BSRational();
    }

    // Helper for converting a string to a BSRational, where the string is the
    // fractional part "y" of the string "x.y".
    private static convertToFraction(numberString: string): BSRational {
        logAssert(/^[0-9]*$/.test(numberString), 'Invalid number format.');
        let y = new BSRational();
        const ten = BSRational.fromNumber(10);
        let pow10 = BSRational.fromNumber(10);
        for (let i = 0; i < numberString.length; ++i) {
            const digit = numberString.charCodeAt(i) - 0x30;
            if (digit > 0) {
                y = y.add(BSRational.fromNumber(digit).div(pow10));
            }
            pow10 = pow10.mul(ten);
        }
        return y;
    }

    // The upstream overloads of the std:: math functions and of the gte::
    // helper functions for BSRational. Except for fabs, frexp and ldexp, they
    // convert to double, call the double version and convert back, so the
    // results are not exact.
    static acos(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.acos(x.toNumber()));
    }

    static acosh(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.acosh(x.toNumber()));
    }

    static asin(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.asin(x.toNumber()));
    }

    static asinh(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.asinh(x.toNumber()));
    }

    static atan(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.atan(x.toNumber()));
    }

    static atanh(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.atanh(x.toNumber()));
    }

    static atan2(y: BSRational, x: BSRational): BSRational {
        return BSRational.fromNumber(Math.atan2(y.toNumber(), x.toNumber()));
    }

    static ceil(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.ceil(x.toNumber()));
    }

    static cos(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.cos(x.toNumber()));
    }

    static cosh(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.cosh(x.toNumber()));
    }

    static exp(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.exp(x.toNumber()));
    }

    static exp2(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.pow(2, x.toNumber()));
    }

    // Exact: the absolute value only drops the sign of the numerator.
    static fabs(x: BSRational): BSRational {
        return (x.getSign() >= 0 ? x.clone() : x.negated());
    }

    static floor(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.floor(x.toNumber()));
    }

    static fmod(x: BSRational, y: BSRational): BSRational {
        const dx = x.toNumber();
        const dy = y.toNumber();
        return BSRational.fromNumber(dx % dy);
    }

    // Exact: x = result * 2^exponent with |result| in [1/2, 1).
    static frexp(x: BSRational): { result: BSRational, exponent: number } {
        if (x.getSign() === 0) {
            return { result: new BSRational(), exponent: 0 };
        }

        const result = x.clone();
        const numer = result.getNumerator();
        const denom = result.getDenominator();
        let e = numer.getExponent() - denom.getExponent();
        numer.setExponent(0);
        denom.setExponent(0);
        const saveSign = numer.getSign();
        numer.setSign(1);
        if (numer.greaterThanOrEqual(denom)) {
            ++e;
            numer.setExponent(-1);
        }
        numer.setSign(saveSign);
        return { result, exponent: e };
    }

    // Exact: multiplication by a power of two.
    static ldexp(x: BSRational, exponent: number): BSRational {
        const result = x.clone();
        const numerator = result.getNumerator();
        numerator.setBiasedExponent(numerator.getBiasedExponent() + exponent);
        return result;
    }

    static log(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.log(x.toNumber()));
    }

    static log2(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.log2(x.toNumber()));
    }

    static log10(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.log10(x.toNumber()));
    }

    static pow(x: BSRational, y: BSRational): BSRational {
        return BSRational.fromNumber(Math.pow(x.toNumber(), y.toNumber()));
    }

    // The port of std::remainder: x - n*y where n is the integer nearest to
    // x/y, ties to even. Upstream converts both operands to double and calls
    // std::remainder, whose result is exact (and always representable as a
    // double). JavaScript has no built-in equivalent ('%' is std::fmod), so
    // the exact integer computation of BSNumber.remainder is reused. Forming
    // the quotient as a double instead loses exactness above 2^53: for
    // remainder(1e17, 3) the double quotient rounds down by one, and
    // 'dx - n*dy' evaluates to 0 where std::remainder returns 1.
    static remainder(x: BSRational, y: BSRational): BSRational {
        const dx = x.toNumber();
        const dy = y.toNumber();
        if (!Number.isFinite(dx) || Number.isNaN(dy) || dy === 0) {
            // std::remainder is NaN here; BSRational has no NaN
            // representation and upstream's conversion of a NaN yields zero.
            return new BSRational();
        }
        if (!Number.isFinite(dy)) {
            // std::remainder(x, +-infinity) = x.
            return BSRational.fromNumber(dx);
        }
        return BSRational.fromBSNumber(BSNumber.remainder(
            BSNumber.fromNumber(dx), BSNumber.fromNumber(dy)));
    }

    static sin(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.sin(x.toNumber()));
    }

    static sinh(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.sinh(x.toNumber()));
    }

    static sqrt(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.sqrt(x.toNumber()));
    }

    static tan(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.tan(x.toNumber()));
    }

    static tanh(x: BSRational): BSRational {
        return BSRational.fromNumber(Math.tanh(x.toNumber()));
    }

    static atandivpi(x: BSRational): BSRational {
        return BSRational.fromNumber(atandivpi(x.toNumber()));
    }

    static atan2divpi(y: BSRational, x: BSRational): BSRational {
        return BSRational.fromNumber(atan2divpi(y.toNumber(), x.toNumber()));
    }

    static clamp(x: BSRational, xmin: BSRational, xmax: BSRational): BSRational {
        return BSRational.fromNumber(clamp(x.toNumber(), xmin.toNumber(), xmax.toNumber()));
    }

    static cospi(x: BSRational): BSRational {
        return BSRational.fromNumber(cospi(x.toNumber()));
    }

    static exp10(x: BSRational): BSRational {
        return BSRational.fromNumber(exp10(x.toNumber()));
    }

    static invsqrt(x: BSRational): BSRational {
        return BSRational.fromNumber(invsqrt(x.toNumber()));
    }

    static isign(x: BSRational): number {
        return isign(x.toNumber());
    }

    static saturate(x: BSRational): BSRational {
        return BSRational.fromNumber(saturate(x.toNumber()));
    }

    static sign(x: BSRational): BSRational {
        return BSRational.fromNumber(sign(x.toNumber()));
    }

    static sinpi(x: BSRational): BSRational {
        return BSRational.fromNumber(sinpi(x.toNumber()));
    }

    static sqr(x: BSRational): BSRational {
        return BSRational.fromNumber(sqr(x.toNumber()));
    }

    // Compute u * v + w exactly.
    static fma(u: BSRational, v: BSRational, w: BSRational): BSRational {
        return u.mul(v).add(w);
    }

    // Sum of products (SOP) u*v+w*z, computed exactly.
    static robustSOP(u: BSRational, v: BSRational, w: BSRational, z: BSRational): BSRational {
        return u.mul(v).add(w.mul(z));
    }

    // Difference of products (DOP) u*v-w*z, computed exactly.
    static robustDOP(u: BSRational, v: BSRational, w: BSRational, z: BSRational): BSRational {
        return u.mul(v).sub(w.mul(z));
    }
}

// Explicit conversion of a BSRational to a BSNumber with a user-specified
// precision (the number of bits of the significand) and rounding mode. This
// is the port of the upstream Convert(input, precision, roundingMode,
// BSNumber& output).
export function convertBSRationalToBSNumber(input: BSRational, precision: number,
    roundingMode: BSNumberRoundingMode): BSNumber {
    if (precision <= 0) {
        logError('Precision must be positive.');
    }

    if (input.getSign() === 0) {
        return new BSNumber();
    }

    let n = input.getNumerator().clone();
    const d = input.getDenominator().clone();

    // The ratio is abstractly of the form n/d = (1.u*2^p)/(1.v*2^q). Convert
    // to the form
    //   (1.u/1.v)*2^{p-q}, if 1.u >= 1.v
    //   2*(1.u/1.v)*2^{p-q-1} if 1.u < 1.v
    // which are in the interval [1,2).
    const signValue = n.getSign() * d.getSign();
    n.setSign(1);
    d.setSign(1);
    let pmq = n.getExponent() - d.getExponent();
    n.setExponent(0);
    d.setExponent(0);
    if (n.lessThan(d)) {
        n.setExponent(n.getExponent() + 1);
        --pmq;
    }

    // Let p = precision. At this time, n/d = 1.c in [1,2). Define the
    // sequence of bits w = 1c = w_{p-1} w_{p-2} ... w_0 r, where w_{p-1} = 1.
    // The bits r after w_0 are used for rounding based on the user-specified
    // rounding mode.
    //
    // Compute p bits for w, the leading bit guaranteed to be 1 and occurring
    // at bit index precision-1. Upstream fills the 32-bit blocks of a
    // UInteger; here the bits are accumulated directly in a bigint, so the
    // bit written at loop index i is bit i of w.
    const two = BSNumber.fromNumber(2);
    let w = 0n;
    const precisionM1 = precision - 1;
    let lastBit = -1;
    for (let i = precisionM1; i >= 0; --i) {
        if (n.lessThan(d)) {
            n = two.mul(n);
            lastBit = 0;
        } else {
            n = two.mul(n.sub(d));
            w |= (1n << BigInt(i));
            lastBit = 1;
            if (n.getSign() === 0) {
                // The input rational has a finite number of bits in its
                // representation, so it is exactly a BSNumber.
                if (i > 0) {
                    // The number n is zero for the remainder of the loop, so
                    // the last bit of the p-bit precision pattern is a zero.
                    // There is no need to continue looping.
                    lastBit = 0;
                }
                break;
            }
        }
    }

    // The port of 'pmq += w.RoundUp()', which adds one to the integer and
    // normalizes it back to an odd number.
    const roundUp = (): void => {
        const rounded = shiftRightToOdd(w + 1n);
        w = rounded.odd;
        pmq += rounded.shift;
    };

    // At this point as a sequence of bits, r = n/d = r0 r1 ...
    if (roundingMode === BSNumberRoundingMode.FE_TONEAREST) {
        n = n.sub(d);
        if (n.getSign() > 0 || (n.getSign() === 0 && lastBit === 1)) {
            roundUp();
        }
        // else round down, equivalent to truncating the r bits
    } else if (roundingMode === BSNumberRoundingMode.FE_UPWARD) {
        if (n.getSign() > 0 && signValue > 0) {
            roundUp();
        }
        // else round down, equivalent to truncating the r bits
    } else if (roundingMode === BSNumberRoundingMode.FE_DOWNWARD) {
        if (n.getSign() > 0 && signValue < 0) {
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
    if ((w & 1n) === 0n) {
        const shifted = shiftRightToOdd(w);
        w = shifted.odd;
        pmq += shifted.shift;
    }

    // Do not use setExponent(pmq) at this step. The number of requested bits
    // is 'precision' but the number of bits of w will be different when
    // round-up occurs, and setExponent uses that count.
    return BSNumber.fromParts(signValue, pmq - precisionM1, w);
}

// The port of the Convert overload whose output is a BSRational. It is used
// to avoid having to expose BSNumber in the APConversion class as well as
// other places where BSRational is passed via a template parameter named
// Rational.
export function convertBSRational(input: BSRational, precision: number,
    roundingMode: BSNumberRoundingMode): BSRational {
    return BSRational.fromBSNumber(
        convertBSRationalToBSNumber(input, precision, roundingMode));
}

// Convert to a 'number' (the C++ double) using the specified rounding mode.
export function convertBSRationalToNumber(input: BSRational,
    roundingMode: BSNumberRoundingMode): number {
    return convertBSRationalToBSNumber(input,
        IEEEBinary64.NUM_SIGNIFICAND_BITS, roundingMode).toNumber();
}

// Convert to a 'number' rounded to the IEEE binary32 representation (the C++
// float) using the specified rounding mode.
export function convertBSRationalToFloat32(input: BSRational,
    roundingMode: BSNumberRoundingMode): number {
    return convertBSRationalToBSNumber(input,
        IEEEBinary32.NUM_SIGNIFICAND_BITS, roundingMode).toFloat32();
}
