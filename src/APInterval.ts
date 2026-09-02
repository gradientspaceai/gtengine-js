// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) APInterval.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Arbitrary-precision (AP) interval arithmetic. Unlike SWInterval, whose
// endpoints are IEEE binary64 numbers that must be widened by one ulp after
// every operation, the endpoints here are exact arbitrary-precision numbers,
// so each operation on the endpoints is exact and no rounding compensation
// is needed. The interval [e0,e1] must satisfy e0 <= e1.
//
// Port notes:
//   * Upstream is the class template APInterval<APType>, where APType must
//     be an arbitrary-precision type (the static_assert on
//     is_arbitrary_precision<APType>). The division, reciprocal and Reals
//     members are enabled only when APType has a division operator, which
//     among the ported arbitrary-precision types means BSRational (BSNumber
//     declares hasDivisionOperator = false and QFNumber is ported over
//     'number' coefficients only). Division is integral to this class - the
//     entire operator/ family plus Reals() depends on it - so, following the
//     precedent that a template layer with a single meaningful instantiation
//     is dropped (B34), this port is the concrete BSRational instantiation.
//   * The C++ free operators become instance methods returning new objects
//     (SWInterval/BSPrecision precedent):
//       operator- (unary) -> negate
//       operator+ -> add, operator- -> sub, operator* -> mul, operator/ ->
//         div (each accepts an APInterval or a BSRational scalar, porting
//         the operator(interval, APType) overloads)
//       operator-(APType, interval) -> static scalarSub
//       operator/(APType, interval) -> static scalarDiv
//         (the non-commutative scalar-on-the-left overloads; s + v and s * v
//          are covered by v.add(s) and v.mul(s), which produce identical
//          endpoints)
//     The compound assignment operators (+=, -=, *=, /=) are not ported; use
//     u = u.add(v) and so on. operator[] becomes get(i). The static leaf-node
//     operations Add/Sub/Mul/Div keep their upstream 2-argument and internal
//     4-argument overloads as the static methods add/sub/mul/div; Mul2,
//     Reciprocal, ReciprocalDown, ReciprocalUp and Reals become mul2,
//     reciprocal, reciprocalDown, reciprocalUp and reals.
//   * C++ copies endpoints by value. TypeScript objects alias and BSRational
//     has mutating members (setSign, negate), so the constructor clones its
//     inputs and get()/getEndpoints() return clones. An APInterval is
//     therefore immutable in practice, as upstream documents.
//   * Sign tests against zero are performed with BSRational.getSign() rather
//     than by comparing against a constructed zero. The two are equivalent
//     for finite endpoints, and getSign() additionally orders the
//     "infinite" endpoints correctly (see the note on reals() below).
//   * The upstream macro GTE_THROW_ON_INVALID_APINTERVAL (disabled by
//     default) is ported as the static flag APInterval.throwOnInvalid:
//     enable it to trap construction of an interval [e0,e1] with e0 > e1.
//
// Infinite endpoints: an indeterminate result (division by an interval that
// straddles zero, or by the scalar zero) is the whole real line, which
// upstream represents by BSRational endpoints whose sign has been set to -2
// and +2 by SetSign. Neither BSNumber nor BSRational has a representation
// for infinity; sign +-2 is a sentinel outside the class invariant (the
// unsigned integer of such a number is zero, which the invariant reserves
// for sign 0). It behaves correctly for the sign tests this class performs
// - getSign() returns -2 or +2, and comparisons against zero take the
// zero-sign branch - but it is not ordered correctly against a finite
// nonzero endpoint of the same sign, and arithmetic on it is meaningless.
// See "Upstream bug suspects" in the port PR. Use isInfinite(e) to detect
// such an endpoint; an interval containing one carries no information other
// than "indeterminate".

import { BSRational } from './BSRational';
import { logAssert } from './Logger';

// Is the endpoint one of the -infinity/+infinity sentinels produced by
// APInterval.reals(), reciprocalDown() and reciprocalUp()? The sentinel is
// a zero-magnitude BSRational whose sign has been set to -2 or +2. A
// sentinel scaled by a nonzero finite number keeps the zero magnitude and
// multiplies the sign, so the test accepts any |sign| >= 2.
export function isInfinite(e: BSRational): boolean {
    return Math.abs(e.getSign()) >= 2;
}

// The port of 'APType posinf(0); posinf.SetSign(+2);'.
function positiveInfinity(): BSRational {
    const posinf = new BSRational();
    posinf.setSign(+2);
    return posinf;
}

// The port of 'APType neginf(0); neginf.SetSign(-2);'.
function negativeInfinity(): BSRational {
    const neginf = new BSRational();
    neginf.setSign(-2);
    return neginf;
}

// The minimum and maximum of two arbitrary-precision numbers (the port of
// std::min and std::max as used by Mul2).
function apMin(u: BSRational, v: BSRational): BSRational {
    return v.lessThan(u) ? v : u;
}

function apMax(u: BSRational, v: BSRational): BSRational {
    return u.lessThan(v) ? v : u;
}

export class APInterval {
    // Port of the upstream macro GTE_THROW_ON_INVALID_APINTERVAL.
    static throwOnInvalid = false;

    private mEndpoints: [BSRational, BSRational];

    // Construction. This is the only way to create an interval. All such
    // intervals are immutable once created. new APInterval() creates [0,0],
    // new APInterval(e) creates the degenerate interval [e,e] and
    // new APInterval(e0, e1) creates [e0,e1].
    constructor(e0?: BSRational, e1?: BSRational) {
        if (e0 === undefined) {
            this.mEndpoints = [new BSRational(), new BSRational()];
        } else if (e1 === undefined) {
            this.mEndpoints = [e0.clone(), e0.clone()];
        } else {
            this.mEndpoints = [e0.clone(), e1.clone()];
            if (APInterval.throwOnInvalid) {
                logAssert(e0.lessThanOrEqual(e1), 'Invalid interval.');
            }
        }
    }

    // The port of the C++ constructor taking std::array<APType, 2>.
    static fromEndpoints(endpoint: readonly [BSRational, BSRational]): APInterval {
        return new APInterval(endpoint[0], endpoint[1]);
    }

    // C++ copy construction/assignment copies by value; TS objects alias, so
    // copies are made explicit.
    clone(): APInterval {
        return new APInterval(this.mEndpoints[0], this.mEndpoints[1]);
    }

    // Member access (the port of operator[]). It is only possible to read
    // the endpoints. You cannot modify the endpoints outside the arithmetic
    // operations, so a clone is returned.
    get(i: number): BSRational {
        return this.mEndpoints[i].clone();
    }

    getEndpoints(): [BSRational, BSRational] {
        return [this.mEndpoints[0].clone(), this.mEndpoints[1].clone()];
    }

    // Arithmetic operations to compute intervals at the leaf nodes of an
    // expression tree (the 2-argument overloads, where u and v are the raw
    // arbitrary-precision variables of the expression). The 4-argument
    // overloads are the upstream internal-use functions consumed by the
    // instance operator methods at the interior nodes of the expression
    // tree.
    static add(u: BSRational, v: BSRational): APInterval;
    static add(u0: BSRational, u1: BSRational, v0: BSRational, v1: BSRational): APInterval;
    static add(a0: BSRational, a1: BSRational, a2?: BSRational, a3?: BSRational): APInterval {
        if (a2 === undefined || a3 === undefined) {
            return new APInterval(a0.add(a1));
        }
        return new APInterval(a0.add(a2), a1.add(a3));
    }

    static sub(u: BSRational, v: BSRational): APInterval;
    static sub(u0: BSRational, u1: BSRational, v0: BSRational, v1: BSRational): APInterval;
    static sub(a0: BSRational, a1: BSRational, a2?: BSRational, a3?: BSRational): APInterval {
        if (a2 === undefined || a3 === undefined) {
            return new APInterval(a0.sub(a1));
        }
        return new APInterval(a0.sub(a3), a1.sub(a2));
    }

    static mul(u: BSRational, v: BSRational): APInterval;
    static mul(u0: BSRational, u1: BSRational, v0: BSRational, v1: BSRational): APInterval;
    static mul(a0: BSRational, a1: BSRational, a2?: BSRational, a3?: BSRational): APInterval {
        if (a2 === undefined || a3 === undefined) {
            return new APInterval(a0.mul(a1));
        }
        return new APInterval(a0.mul(a2), a1.mul(a3));
    }

    static mul2(u0: BSRational, u1: BSRational, v0: BSRational, v1: BSRational): APInterval {
        const u0mv1 = u0.mul(v1);
        const u1mv0 = u1.mul(v0);
        const u0mv0 = u0.mul(v0);
        const u1mv1 = u1.mul(v1);
        return new APInterval(apMin(u0mv1, u1mv0), apMax(u0mv0, u1mv1));
    }

    static div(u: BSRational, v: BSRational): APInterval;
    static div(u0: BSRational, u1: BSRational, v0: BSRational, v1: BSRational): APInterval;
    static div(a0: BSRational, a1: BSRational, a2?: BSRational, a3?: BSRational): APInterval {
        if (a2 === undefined || a3 === undefined) {
            if (a1.getSign() !== 0) {
                return new APInterval(a0.div(a1));
            }
            // Division by zero does not lead to a determinate interval. Just
            // return the entire set of real numbers.
            return APInterval.reals();
        }
        return new APInterval(a0.div(a3), a1.div(a2));
    }

    static reciprocal(v0: BSRational, v1: BSRational): APInterval {
        const one = BSRational.fromNumber(1);
        return new APInterval(one.div(v1), one.div(v0));
    }

    static reciprocalDown(v: BSRational): APInterval {
        const recpv = BSRational.fromNumber(1).div(v);
        return new APInterval(recpv, positiveInfinity());
    }

    static reciprocalUp(v: BSRational): APInterval {
        const recpv = BSRational.fromNumber(1).div(v);
        return new APInterval(negativeInfinity(), recpv);
    }

    static reals(): APInterval {
        return new APInterval(negativeInfinity(), positiveInfinity());
    }

    // Unary operator-. Negation of [e0,e1] produces [-e1,-e0]. This
    // operation needs to be supported in the sense of negating a "number"
    // in an arithmetic expression.
    negate(): APInterval {
        return new APInterval(this.mEndpoints[1].negated(), this.mEndpoints[0].negated());
    }

    // Addition operations.
    add(v: APInterval | BSRational): APInterval {
        const u = this.mEndpoints;
        if (v instanceof APInterval) {
            return APInterval.add(u[0], u[1], v.mEndpoints[0], v.mEndpoints[1]);
        }
        return APInterval.add(u[0], u[1], v, v);
    }

    // Subtraction operations.
    sub(v: APInterval | BSRational): APInterval {
        const u = this.mEndpoints;
        if (v instanceof APInterval) {
            return APInterval.sub(u[0], u[1], v.mEndpoints[0], v.mEndpoints[1]);
        }
        return APInterval.sub(u[0], u[1], v, v);
    }

    // The port of C++ operator-(APType const& u, APInterval const& v).
    static scalarSub(u: BSRational, v: APInterval): APInterval {
        return APInterval.sub(u, u, v.mEndpoints[0], v.mEndpoints[1]);
    }

    // Multiplication operations.
    mul(v: APInterval | BSRational): APInterval {
        const u = this.mEndpoints;
        if (!(v instanceof APInterval)) {
            if (v.getSign() >= 0) {
                return APInterval.mul(u[0], u[1], v, v);
            } else {
                return APInterval.mul(u[1], u[0], v, v);
            }
        }
        const w = v.mEndpoints;
        if (u[0].getSign() >= 0) {
            if (w[0].getSign() >= 0) {
                return APInterval.mul(u[0], u[1], w[0], w[1]);
            } else if (w[1].getSign() <= 0) {
                return APInterval.mul(u[1], u[0], w[0], w[1]);
            } else { // w[0] < 0 < w[1]
                return APInterval.mul(u[1], u[1], w[0], w[1]);
            }
        } else if (u[1].getSign() <= 0) {
            if (w[0].getSign() >= 0) {
                return APInterval.mul(u[0], u[1], w[1], w[0]);
            } else if (w[1].getSign() <= 0) {
                return APInterval.mul(u[1], u[0], w[1], w[0]);
            } else { // w[0] < 0 < w[1]
                return APInterval.mul(u[0], u[0], w[1], w[0]);
            }
        } else { // u[0] < 0 < u[1]
            if (w[0].getSign() >= 0) {
                return APInterval.mul(u[0], u[1], w[1], w[1]);
            } else if (w[1].getSign() <= 0) {
                return APInterval.mul(u[1], u[0], w[0], w[0]);
            } else { // w[0] < 0 < w[1]
                return APInterval.mul2(u[0], u[1], w[0], w[1]);
            }
        }
    }

    // Division operations. If the divisor interval is [v0,v1] with
    // v0 < 0 < v1, then the returned interval is (-infinity,+infinity)
    // instead of Union((-infinity,1/v0),(1/v1,+infinity)). An application
    // should try to avoid this case by branching based on [v0,0] and [0,v1].
    div(v: APInterval | BSRational): APInterval {
        const u = this.mEndpoints;
        if (!(v instanceof APInterval)) {
            if (v.getSign() > 0) {
                return APInterval.div(u[0], u[1], v, v);
            } else if (v.getSign() < 0) {
                return APInterval.div(u[1], u[0], v, v);
            } else { // v = 0
                return APInterval.reals();
            }
        }
        const w = v.mEndpoints;
        if (w[0].getSign() > 0 || w[1].getSign() < 0) {
            return this.mul(APInterval.reciprocal(w[0], w[1]));
        } else {
            if (w[0].getSign() === 0) {
                return this.mul(APInterval.reciprocalDown(w[1]));
            } else if (w[1].getSign() === 0) {
                return this.mul(APInterval.reciprocalUp(w[0]));
            } else { // w[0] < 0 < w[1]
                return APInterval.reals();
            }
        }
    }

    // The port of C++ operator/(APType const& u, APInterval const& v). The
    // scalar u is treated as the degenerate interval [u,u], for which
    // 'u * W' and 'W.mul(u)' produce identical endpoints.
    static scalarDiv(u: BSRational, v: APInterval): APInterval {
        const w = v.mEndpoints;
        if (w[0].getSign() > 0 || w[1].getSign() < 0) {
            return APInterval.reciprocal(w[0], w[1]).mul(u);
        } else {
            if (w[0].getSign() === 0) {
                return APInterval.reciprocalDown(w[1]).mul(u);
            } else if (w[1].getSign() === 0) {
                return APInterval.reciprocalUp(w[0]).mul(u);
            } else { // w[0] < 0 < w[1]
                return APInterval.reals();
            }
        }
    }
}
