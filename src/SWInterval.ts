// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SWInterval.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Software (SW) interval arithmetic: each operation computes in round-to-
// nearest floating-point and then widens the resulting bounds by one ulp in
// the appropriate directions so that the returned interval encloses the
// exact real result. Upstream widens with std::nextafter(value, -max) and
// std::nextafter(value, +max).
//
// PORT DEVIATION (directed rounding / nextafter): JavaScript has neither
// C++ fesetround directed rounding nor std::nextafter. This port implements
// the ulp-adjustment helpers nextDown(x) and nextUp(x) via IEEE-754 binary64
// bit manipulation (Float64Array/BigUint64Array views, the same technique as
// the ported IEEEBinary.ts). They reproduce std::nextafter(x, -MAX_VALUE)
// and std::nextafter(x, +MAX_VALUE) exactly, including the behavior at
// zeros, subnormals, infinities (nextDown(+inf) = MAX_VALUE, nextUp(-inf) =
// -MAX_VALUE, matching std::nextafter stepping from an infinity toward
// +/-max) and NaN. The enclosure property is therefore preserved: every
// computed bound is moved one ulp outward from the round-to-nearest result,
// so result intervals are conservative (at most one ulp wider per bound than
// upstream's would be with true directed rounding — identical to upstream
// SWInterval, which uses the same nextafter widening). The one caveat,
// inherited from upstream, is overflow: if a bound computation overflows to
// an infinity, nextafter toward +/-max pulls it back to +/-MAX_VALUE.
// FPInterval.h (a separate batch), which upstream implements with
// fesetround directed rounding, will reuse these nextUp/nextDown helpers
// and widen each bound by one ulp where upstream changes the rounding mode.
//
// Port notes: only T = number (C++ double) is ported. The C++ free
// operators become instance methods returning new objects (BSPrecision
// precedent):
//   operator- (unary) -> negate
//   operator+ -> add, operator- -> sub, operator* -> mul, operator/ -> div
//     (each accepts an SWInterval or a scalar 'number', porting the
//      operator(interval, T) overloads)
//   operator-(T, interval) -> static scalarSub
//   operator/(T, interval) -> static scalarDiv
//     (the non-commutative scalar-on-the-left overloads; s + v and s * v
//      are covered by v.add(s) and v.mul(s), which produce identical
//      endpoints)
// The compound assignment operators (+=, -=, *=, /=) are not ported; use
// u = u.add(v) and so on. operator[] becomes get(i). The static leaf-node
// operations Add/Sub/Mul/Div keep their upstream 2-argument and internal
// 4-argument overloads as the static methods add/sub/mul/div; Mul2,
// Reciprocal, ReciprocalDown, ReciprocalUp and Reals become mul2,
// reciprocal, reciprocalDown, reciprocalUp and reals.
//
// The upstream macro GTE_THROW_ON_INVALID_SWINTERVAL (disabled by default)
// is ported as the static flag SWInterval.throwOnInvalid: enable it to trap
// construction of an interval [e0, e1] with e0 > e1.

import { logAssert } from './Logger.js';

// Scratch views for converting between floating-point values and their bit
// patterns.
const scratchBuffer = new ArrayBuffer(8);
const scratchF = new Float64Array(scratchBuffer);
const scratchU = new BigUint64Array(scratchBuffer);

// The port of std::nextafter(x, -std::numeric_limits<double>::max()): the
// next representable binary64 value from x in the direction of -MAX_VALUE.
export function nextDown(x: number): number {
    if (Number.isNaN(x)) {
        return x;
    }
    if (x === -Number.MAX_VALUE) {
        // nextafter(from, to) returns 'to' when from equals to.
        return x;
    }
    if (x === 0) {
        // The next-down for both -0 and +0 is -MIN_SUBNORMAL.
        return -Number.MIN_VALUE;
    }
    scratchF[0] = x;
    if (x > 0) {
        // Includes nextDown(+infinity) = +MAX_VALUE.
        scratchU[0] -= 1n;
    } else {
        // nextDown(-infinity) steps toward -MAX_VALUE, as std::nextafter
        // does.
        if (x === Number.NEGATIVE_INFINITY) {
            return -Number.MAX_VALUE;
        }
        scratchU[0] += 1n;
    }
    return scratchF[0];
}

// The port of std::nextafter(x, +std::numeric_limits<double>::max()): the
// next representable binary64 value from x in the direction of +MAX_VALUE.
export function nextUp(x: number): number {
    if (Number.isNaN(x)) {
        return x;
    }
    if (x === Number.MAX_VALUE) {
        // nextafter(from, to) returns 'to' when from equals to.
        return x;
    }
    if (x === 0) {
        // The next-up for both -0 and +0 is +MIN_SUBNORMAL.
        return Number.MIN_VALUE;
    }
    scratchF[0] = x;
    if (x > 0) {
        // nextUp(+infinity) steps toward +MAX_VALUE, as std::nextafter does.
        if (x === Number.POSITIVE_INFINITY) {
            return Number.MAX_VALUE;
        }
        scratchU[0] += 1n;
    } else {
        // Includes nextUp(-infinity) = -MAX_VALUE.
        scratchU[0] -= 1n;
    }
    return scratchF[0];
}

// The SWInterval [e0, e1] must satisfy e0 <= e1. Set
// SWInterval.throwOnInvalid to trap invalid construction where e0 > e1.
export class SWInterval {
    // Convenient constants.
    static readonly zero: number = 0;
    static readonly one: number = 1;
    static readonly max: number = Number.MAX_VALUE;
    static readonly inf: number = Number.POSITIVE_INFINITY;

    // Port of the upstream macro GTE_THROW_ON_INVALID_SWINTERVAL.
    static throwOnInvalid = false;

    private mEndpoints: [number, number];

    // Construction. All such intervals are conceptually immutable once
    // created; the endpoints cannot be modified outside the arithmetic
    // operations. new SWInterval() creates [0, 0], new SWInterval(e)
    // creates the degenerate interval [e, e] and new SWInterval(e0, e1)
    // creates [e0, e1].
    constructor(e0: number = 0, e1?: number) {
        if (e1 === undefined) {
            this.mEndpoints = [e0, e0];
        } else {
            this.mEndpoints = [e0, e1];
            if (SWInterval.throwOnInvalid) {
                logAssert(e0 <= e1, 'Invalid SWInterval.');
            }
        }
    }

    // The port of the C++ constructor taking std::array<T, 2>.
    static fromEndpoints(endpoint: readonly [number, number]): SWInterval {
        return new SWInterval(endpoint[0], endpoint[1]);
    }

    // C++ copy construction/assignment copies by value; TS objects alias,
    // so copies are made explicit.
    clone(): SWInterval {
        return new SWInterval(this.mEndpoints[0], this.mEndpoints[1]);
    }

    // Member access (the port of operator[]). It is only possible to read
    // the endpoints.
    get(i: number): number {
        return this.mEndpoints[i];
    }

    getEndpoints(): [number, number] {
        return [this.mEndpoints[0], this.mEndpoints[1]];
    }

    // Arithmetic operations to compute intervals at the leaf nodes of an
    // expression tree (the 2-argument overloads, where u and v are the raw
    // floating-point variables of the expression). The 4-argument overloads
    // are the upstream internal-use functions consumed by the instance
    // operator methods at the interior nodes of the expression tree.
    static add(u: number, v: number): SWInterval;
    static add(u0: number, u1: number, v0: number, v1: number): SWInterval;
    static add(a0: number, a1: number, a2?: number, a3?: number): SWInterval {
        if (a2 === undefined || a3 === undefined) {
            const add = a0 + a1;
            return new SWInterval(nextDown(add), nextUp(add));
        }
        return new SWInterval(nextDown(a0 + a2), nextUp(a1 + a3));
    }

    static sub(u: number, v: number): SWInterval;
    static sub(u0: number, u1: number, v0: number, v1: number): SWInterval;
    static sub(a0: number, a1: number, a2?: number, a3?: number): SWInterval {
        if (a2 === undefined || a3 === undefined) {
            const sub = a0 - a1;
            return new SWInterval(nextDown(sub), nextUp(sub));
        }
        return new SWInterval(nextDown(a0 - a3), nextUp(a1 - a2));
    }

    static mul(u: number, v: number): SWInterval;
    static mul(u0: number, u1: number, v0: number, v1: number): SWInterval;
    static mul(a0: number, a1: number, a2?: number, a3?: number): SWInterval {
        if (a2 === undefined || a3 === undefined) {
            const mul = a0 * a1;
            return new SWInterval(nextDown(mul), nextUp(mul));
        }
        return new SWInterval(nextDown(a0 * a2), nextUp(a1 * a3));
    }

    static mul2(u0: number, u1: number, v0: number, v1: number): SWInterval {
        const u0mv1 = nextDown(u0 * v1);
        const u1mv0 = nextDown(u1 * v0);
        const u0mv0 = nextUp(u0 * v0);
        const u1mv1 = nextUp(u1 * v1);
        return new SWInterval(Math.min(u0mv1, u1mv0), Math.max(u0mv0, u1mv1));
    }

    static div(u: number, v: number): SWInterval;
    static div(u0: number, u1: number, v0: number, v1: number): SWInterval;
    static div(a0: number, a1: number, a2?: number, a3?: number): SWInterval {
        if (a2 === undefined || a3 === undefined) {
            if (a1 !== SWInterval.zero) {
                const div = a0 / a1;
                return new SWInterval(nextDown(div), nextUp(div));
            }
            // Division by zero does not lead to a determinate SWInterval.
            // Return the entire set of real numbers.
            return SWInterval.reals();
        }
        return new SWInterval(nextDown(a0 / a3), nextUp(a1 / a2));
    }

    static reciprocal(v0: number, v1: number): SWInterval {
        return new SWInterval(
            nextDown(SWInterval.one / v1), nextUp(SWInterval.one / v0));
    }

    static reciprocalDown(v: number): SWInterval {
        const recpv = nextDown(SWInterval.one / v);
        return new SWInterval(recpv, SWInterval.inf);
    }

    static reciprocalUp(v: number): SWInterval {
        const recpv = nextUp(SWInterval.one / v);
        return new SWInterval(-SWInterval.inf, recpv);
    }

    static reals(): SWInterval {
        return new SWInterval(-SWInterval.inf, SWInterval.inf);
    }

    // Unary operator-. Negation of [e0, e1] produces [-e1, -e0]. This
    // operation needs to be supported in the sense of negating a "number"
    // in an arithmetic expression.
    negate(): SWInterval {
        return new SWInterval(-this.mEndpoints[1], -this.mEndpoints[0]);
    }

    // Addition operations.
    add(v: SWInterval | number): SWInterval {
        const u = this.mEndpoints;
        if (typeof v === 'number') {
            return SWInterval.add(u[0], u[1], v, v);
        }
        return SWInterval.add(u[0], u[1], v.mEndpoints[0], v.mEndpoints[1]);
    }

    // Subtraction operations.
    sub(v: SWInterval | number): SWInterval {
        const u = this.mEndpoints;
        if (typeof v === 'number') {
            return SWInterval.sub(u[0], u[1], v, v);
        }
        return SWInterval.sub(u[0], u[1], v.mEndpoints[0], v.mEndpoints[1]);
    }

    // The port of C++ operator-(T u, SWInterval const& v).
    static scalarSub(u: number, v: SWInterval): SWInterval {
        return SWInterval.sub(u, u, v.mEndpoints[0], v.mEndpoints[1]);
    }

    // Multiplication operations.
    mul(v: SWInterval | number): SWInterval {
        const u = this.mEndpoints;
        if (typeof v === 'number') {
            if (v >= SWInterval.zero) {
                return SWInterval.mul(u[0], u[1], v, v);
            } else {
                return SWInterval.mul(u[1], u[0], v, v);
            }
        }
        const w = v.mEndpoints;
        if (u[0] >= SWInterval.zero) {
            if (w[0] >= SWInterval.zero) {
                return SWInterval.mul(u[0], u[1], w[0], w[1]);
            } else if (w[1] <= SWInterval.zero) {
                return SWInterval.mul(u[1], u[0], w[0], w[1]);
            } else { // w[0] < 0 < w[1]
                return SWInterval.mul(u[1], u[1], w[0], w[1]);
            }
        } else if (u[1] <= SWInterval.zero) {
            if (w[0] >= SWInterval.zero) {
                return SWInterval.mul(u[0], u[1], w[1], w[0]);
            } else if (w[1] <= SWInterval.zero) {
                return SWInterval.mul(u[1], u[0], w[1], w[0]);
            } else { // w[0] < 0 < w[1]
                return SWInterval.mul(u[0], u[0], w[1], w[0]);
            }
        } else { // u[0] < 0 < u[1]
            if (w[0] >= SWInterval.zero) {
                return SWInterval.mul(u[0], u[1], w[1], w[1]);
            } else if (w[1] <= SWInterval.zero) {
                return SWInterval.mul(u[1], u[0], w[0], w[0]);
            } else { // w[0] < 0 < w[1]
                return SWInterval.mul2(u[0], u[1], w[0], w[1]);
            }
        }
    }

    // Division operations. If the divisor SWInterval is [v0, v1] with
    // v0 < 0 < v1, then the returned SWInterval is (-inf, +inf) instead of
    // Union((-inf, 1/v0), (1/v1, +inf)). An application should try to avoid
    // this case by branching based on [v0, 0] and [0, v1].
    div(v: SWInterval | number): SWInterval {
        const u = this.mEndpoints;
        if (typeof v === 'number') {
            if (v > SWInterval.zero) {
                return SWInterval.div(u[0], u[1], v, v);
            } else if (v < SWInterval.zero) {
                return SWInterval.div(u[1], u[0], v, v);
            } else { // v = 0
                return SWInterval.reals();
            }
        }
        const w = v.mEndpoints;
        if (w[0] > SWInterval.zero || w[1] < SWInterval.zero) {
            return this.mul(SWInterval.reciprocal(w[0], w[1]));
        } else {
            if (w[0] === SWInterval.zero) {
                return this.mul(SWInterval.reciprocalDown(w[1]));
            } else if (w[1] === SWInterval.zero) {
                return this.mul(SWInterval.reciprocalUp(w[0]));
            } else { // w[0] < 0 < w[1]
                return SWInterval.reals();
            }
        }
    }

    // The port of C++ operator/(T u, SWInterval const& v). The scalar u is
    // multiplied by the reciprocal interval; scalar multiplication is
    // commutative, so reciprocal(...).mul(u) produces the same endpoints as
    // the upstream u * Reciprocal(...).
    static scalarDiv(u: number, v: SWInterval): SWInterval {
        const w = v.mEndpoints;
        if (w[0] > SWInterval.zero || w[1] < SWInterval.zero) {
            return SWInterval.reciprocal(w[0], w[1]).mul(u);
        } else {
            if (w[0] === SWInterval.zero) {
                return SWInterval.reciprocalDown(w[1]).mul(u);
            } else if (w[1] === SWInterval.zero) {
                return SWInterval.reciprocalUp(w[0]).mul(u);
            } else { // w[0] < 0 < w[1]
                return SWInterval.reals();
            }
        }
    }
}
