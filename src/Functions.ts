// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Functions.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// These functions are convenient for some applications.
//
// Port notes:
//   - Upstream declares a 'float' overload and a 'double' overload of each
//     function. TypeScript numbers are IEEE binary64, so only the 'double'
//     overloads are ported (see PORTING.md).
//   - Upstream function names FMA, RobustSOP and RobustDOP become fma,
//     robustSOP and robustDOP.
//   - JavaScript has no fused-multiply-add operator, so fma() below is a
//     software implementation that returns the correctly rounded value of
//     u * v + w computed with a single rounding, matching std::fma. Upstream
//     offers GTE_DISCARD_FMA to fall back to 'u * v + w' with two roundings;
//     this port always uses the exact implementation because robustSOP and
//     robustDOP are meaningless without it.

import { GTE_C_INV_PI, GTE_C_PI, GTE_C_LN_10 } from './Constants.js';

export function atandivpi(x: number): number {
    return Math.atan(x) * GTE_C_INV_PI;
}

export function atan2divpi(y: number, x: number): number {
    return Math.atan2(y, x) * GTE_C_INV_PI;
}

export function clamp(x: number, xmin: number, xmax: number): number {
    return (x <= xmin ? xmin : (x >= xmax ? xmax : x));
}

export function cospi(x: number): number {
    return Math.cos(x * GTE_C_PI);
}

export function exp10(x: number): number {
    return Math.exp(x * GTE_C_LN_10);
}

export function invsqrt(x: number): number {
    return 1.0 / Math.sqrt(x);
}

export function isign(x: number): number {
    return (x > 0.0 ? 1 : (x < 0.0 ? -1 : 0));
}

export function saturate(x: number): number {
    return (x <= 0.0 ? 0.0 : (x >= 1.0 ? 1.0 : x));
}

export function sign(x: number): number {
    return (x > 0.0 ? 1.0 : (x < 0.0 ? -1.0 : 0.0));
}

export function sinpi(x: number): number {
    return Math.sin(x * GTE_C_PI);
}

export function sqr(x: number): number {
    return x * x;
}

// Compute u * v + w as a single operation; that is, the exact value of
// u * v + w is computed and rounded once (round to nearest, ties to even).
// This is the semantics of std::fma for IEEE binary64.
export function fma(u: number, v: number, w: number): number {
    // Non-finite inputs. When u or v is infinite, the product u * v is
    // computed without rounding error (it is an infinity, or a NaN for the
    // invalid operation 0 * infinity), so the hardware expression produces
    // the IEEE 754 fusedMultiplyAdd result. When u and v are finite and w is
    // infinite, the result is w; the hardware expression could return NaN
    // instead if the finite product overflowed to the opposite infinity.
    if (Number.isNaN(u) || Number.isNaN(v) || Number.isNaN(w)) {
        return Number.NaN;
    }
    if (!Number.isFinite(u) || !Number.isFinite(v)) {
        return u * v + w;
    }
    if (!Number.isFinite(w)) {
        return w;
    }
    if (u === 0 || v === 0) {
        // The product is an exact (signed) zero, so a single rounding of
        // (+/-0) + w is what the hardware addition already computes.
        return u * v + w;
    }
    if (w === 0) {
        // The exact sum is the exact product, so one rounding of the product
        // is the answer.
        return u * v;
    }

    const du = decompose(u);
    const dv = decompose(v);
    const dw = decompose(w);

    // The exact product is (-1)^ps * mp * 2^ep.
    const ps = du.s ^ dv.s;
    const mp = du.m * dv.m;
    const ep = du.e + dv.e;

    // Align the product and the addend to a common power of two and add
    // exactly. The result is n * 2^e with n a signed big integer.
    const e = Math.min(ep, dw.e);
    const sp = ps !== 0 ? -mp : mp;
    const sw = dw.s !== 0 ? -dw.m : dw.m;
    const n = (sp << BigInt(ep - e)) + (sw << BigInt(dw.e - e));
    return roundScaledInteger(n, e);
}

// Robust sum of products (SOP) u * v + w * z.
export function robustSOP(u: number, v: number, w: number, z: number): number {
    const productWZ = w * z;
    const roundingError = fma(w, z, -productWZ);
    const result = fma(u, v, productWZ) + roundingError;
    return result;
}

// Robust difference of products (DOP) u * v - w * z.
export function robustDOP(u: number, v: number, w: number, z: number): number {
    const productWZ = w * z;
    const roundingError = fma(-w, z, productWZ);
    const result = fma(u, v, -productWZ) + roundingError;
    return result;
}

// Support for the software fused-multiply-add. --------------------------

const fmaView = new DataView(new ArrayBuffer(8));

// Split a finite nonzero number x into (s, m, e) with
// x = (-1)^s * m * 2^e, where m is a nonnegative big integer.
function decompose(x: number): { s: number; m: bigint; e: number } {
    fmaView.setFloat64(0, x);
    const bits = fmaView.getBigUint64(0);
    const s = Number(bits >> 63n);
    const biasedExponent = Number((bits >> 52n) & 0x7FFn);
    const trailing = bits & 0xFFFFFFFFFFFFFn;
    if (biasedExponent === 0) {
        // Subnormal (or zero): the implied leading bit is 0 and the exponent
        // is that of the minimum subnormal, 2^-1074.
        return { s, m: trailing, e: -1074 };
    }
    // Normal: restore the implied leading 1 bit. The significand is scaled by
    // 2^52, hence the exponent bias of 1023 + 52 = 1075.
    return { s, m: trailing | 0x10000000000000n, e: biasedExponent - 1075 };
}

// The number of bits in the binary representation of a positive big integer.
function bitLength(a: bigint): number {
    const hex = a.toString(16);
    return (hex.length - 1) * 4 + (32 - Math.clz32(Number.parseInt(hex[0], 16)));
}

// Round n * 2^e (n a signed big integer) to the nearest binary64 value, ties
// to even. This is the single rounding required by fma.
function roundScaledInteger(n: bigint, e: number): number {
    if (n === 0n) {
        // The exact result is zero. Round to nearest returns +0 for a zero
        // sum of nonzero operands of opposite sign.
        return 0;
    }

    const negative = n < 0n;
    let a = negative ? -n : n;
    const numBits = bitLength(a);

    // The exact magnitude lies in [2^exponent, 2^(exponent+1)).
    const exponent = numBits - 1 + e;

    // The number of low-order bits of 'a' that must be discarded. Normal
    // numbers keep 53 significant bits; subnormal numbers are multiples of
    // 2^-1074.
    const shift = exponent >= -1022 ? numBits - 53 : -1074 - e;
    let scale = e;
    if (shift > 0) {
        const bigShift = BigInt(shift);
        const q = a >> bigShift;
        const remainder = a - (q << bigShift);
        const half = 1n << (bigShift - 1n);
        // Round to nearest, ties to even.
        a = (remainder > half || (remainder === half && (q & 1n) === 1n)) ? q + 1n : q;
        scale = e + shift;
        if (a === 0n) {
            return negative ? -0 : 0;
        }
    }

    // Overflow to infinity when the rounded magnitude is at least 2^1024.
    if (bitLength(a) - 1 + scale > 1023) {
        return negative ? Number.NEGATIVE_INFINITY : Number.POSITIVE_INFINITY;
    }

    const magnitude = ldexp(Number(a), scale);
    return negative ? -magnitude : magnitude;
}

// Compute m * 2^exp without intermediate overflow or underflow. The callers
// guarantee that the exact result is representable, so no double rounding
// occurs.
function ldexp(m: number, exp: number): number {
    let r = m;
    let k = exp;
    while (k > 1023) {
        r *= 2 ** 1023;
        k -= 1023;
    }
    while (k < -1022) {
        r *= 2 ** -1022;
        k += 1022;
    }
    return r * 2 ** k;
}
