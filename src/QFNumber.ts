// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) QFNumber.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Class QFNumber is an implementation for quadratic fields with N >= 1
// square root terms. The theory and details are provided in
// https://www.geometrictools.com/Documentation/QuadraticFields.pdf
//
// Port notes: upstream is the recursive template QFNumber<T, N>, where the
// N = 1 specialization stores coefficients of type T and the general case
// stores coefficients of type QFNumber<T, N - 1>. This port uses the
// number-based path (T = number; upstream typically instantiates T with a
// rational type, which may replace 'number' in a later revision). The
// compile-time recursion becomes a runtime one: a coefficient is either a
// 'number' (the N = 1 case) or a nested QFNumber (the N >= 2 case), captured
// by the QFCoefficient union type. A quadratic field number is
// x[0] + x[1] * sqrt(d) where both coefficients have the same nesting depth.
//
// The C++ free operators become instance methods returning new objects
// (BSPrecision precedent):
//   operator- (unary) -> negate
//   operator+  -> add,  operator-  -> sub,  operator*  -> mul,
//   operator/  -> div  (each accepts a QFNumber or a scalar 'number',
//                       porting the operator(q, s) overloads)
//   operator-(s, q) -> static scalarSub, operator/(s, q) -> static scalarDiv
//     (the non-commutative scalar-on-the-left overloads; s + q and s * q are
//      covered by q.add(s) and q.mul(s))
//   operator== -> equals,      operator!= -> notEquals,
//   operator<  -> lessThan,    operator<= -> lessThanEqual,
//   operator>  -> greaterThan, operator>= -> greaterThanEqual
// The compound assignment operators (+=, -=, *=, /=) are not ported; use
// q = q.add(v) and so on. C++ assignment copies by value; use clone() where
// a copy is needed.
//
// Upstream marks QFNumber as arbitrary precision in TypeTraits.h, so the
// class implements the ArbitraryPrecisionNumber marker interface. Division
// is supported (hasDivisionOperator is true) for the number-based path.
//
// The upstream macro GTE_ASSERT_ON_QFNUMBER_MISMATCHED_D (disabled by
// default) is ported as the static flag QFNumber.assertOnMismatchedD: set it
// to true if you want the logging system to trap when arithmetic operations
// are performed on two quadratic field numbers that do not share the same
// value d.

import { logAssert, logError } from './Logger';
import type { ArbitraryPrecisionNumber } from './TypeTraits';

// A coefficient of a quadratic field number: 'number' at recursion depth
// N = 1, a nested QFNumber at depth N >= 2.
export type QFCoefficient = number | QFNumber;

// Coefficient arithmetic. In upstream the two operands of +, -, * and / on
// coefficients always have the same type (both T or both QFNumber<T, N-1>)
// except for the scalar operator overloads, where the scalar is a T combined
// with a coefficient at any depth.
function cNeg(a: QFCoefficient): QFCoefficient {
    return typeof a === 'number' ? -a : a.negate();
}

function cAdd(a: QFCoefficient, b: QFCoefficient | number): QFCoefficient {
    if (typeof a === 'number') {
        if (typeof b === 'number') {
            return a + b;
        }
        // s + q = q + s (addition of a scalar is commutative).
        return b.add(a);
    }
    return a.add(b);
}

function cSub(a: QFCoefficient, b: QFCoefficient | number): QFCoefficient {
    if (typeof a === 'number') {
        if (typeof b === 'number') {
            return a - b;
        }
        return QFNumber.scalarSub(a, b);
    }
    return a.sub(b);
}

function cMul(a: QFCoefficient, b: QFCoefficient | number): QFCoefficient {
    if (typeof a === 'number') {
        if (typeof b === 'number') {
            return a * b;
        }
        // s * q = q * s (multiplication by a scalar is commutative).
        return b.mul(a);
    }
    return a.mul(b);
}

function cDiv(a: QFCoefficient, b: QFCoefficient | number): QFCoefficient {
    if (typeof a === 'number') {
        if (typeof b === 'number') {
            return a / b;
        }
        return QFNumber.scalarDiv(a, b);
    }
    return a.div(b);
}

// Coefficient comparisons. Upstream compares coefficients only at equal
// depth (both T or both QFNumber<T, N-1>).
function cEqual(a: QFCoefficient, b: QFCoefficient): boolean {
    if (typeof a === 'number' && typeof b === 'number') {
        return a === b;
    }
    if (typeof a !== 'number' && typeof b !== 'number') {
        return a.equals(b);
    }
    logError('Mismatched coefficient depth.');
}

function cLessThan(a: QFCoefficient, b: QFCoefficient): boolean {
    if (typeof a === 'number' && typeof b === 'number') {
        return a < b;
    }
    if (typeof a !== 'number' && typeof b !== 'number') {
        return a.lessThan(b);
    }
    logError('Mismatched coefficient depth.');
}

function cGreaterThan(a: QFCoefficient, b: QFCoefficient): boolean {
    return cLessThan(b, a);
}

function cLessThanEqual(a: QFCoefficient, b: QFCoefficient): boolean {
    return !cLessThan(b, a);
}

function cGreaterThanEqual(a: QFCoefficient, b: QFCoefficient): boolean {
    return !cLessThan(a, b);
}

// Arithmetic for quadratic fields. The quadratic field number is
// x[0] + x[1] * sqrt(d).
export class QFNumber implements ArbitraryPrecisionNumber {
    // Port of the upstream macro GTE_ASSERT_ON_QFNUMBER_MISMATCHED_D:
    // enable to trap arithmetic operations and comparisons on two quadratic
    // field numbers that do not share the same value d.
    static assertOnMismatchedD = false;

    // Arbitrary-precision marker (upstream TypeTraits.h marks QFNumber as
    // arbitrary precision; the number-based path supports division).
    readonly isArbitraryPrecision = true as const;
    readonly hasDivisionOperator: boolean = true;

    // The quadratic field number is x[0] + x[1] * sqrt(d).
    x: [QFCoefficient, QFCoefficient];
    d: number;

    // Create z = x0 + x1 * sqrt(d). The zero-argument form creates
    // z = 0 + 0 * sqrt(0). For N >= 2 numbers, pass QFNumber coefficients
    // of equal depth.
    constructor(x0: QFCoefficient = 0, x1: QFCoefficient = 0, d: number = 0) {
        this.x = [x0, x1];
        this.d = d;
    }

    // Create z = 0 + 0 * sqrt(d) = 0 (the port of the C++ constructor
    // 'explicit QFNumber(T const& inD)', which is ambiguous with the
    // coefficient constructor in TypeScript).
    static fromD(d: number): QFNumber {
        return new QFNumber(0, 0, d);
    }

    // C++ assignment copies by value; TS objects alias, so copies are made
    // explicit. The copy is deep (nested coefficients are cloned).
    clone(): QFNumber {
        const c0 = typeof this.x[0] === 'number' ? this.x[0] : this.x[0].clone();
        const c1 = typeof this.x[1] === 'number' ? this.x[1] : this.x[1].clone();
        return new QFNumber(c0, c1, this.d);
    }

    private static checkD(q0: QFNumber, q1: QFNumber): void {
        if (QFNumber.assertOnMismatchedD) {
            logAssert(q0.d === q1.d, 'Mismatched d-value.');
        }
    }

    // Unary operator-.
    negate(): QFNumber {
        return new QFNumber(cNeg(this.x[0]), cNeg(this.x[1]), this.d);
    }

    // Arithmetic operations between elements of a quadratic field must occur
    // only when the d-values are the same. To trap mismatches, read the
    // comments about assertOnMismatchedD at the beginning of this file.
    add(other: QFNumber | number): QFNumber {
        if (typeof other === 'number') {
            return new QFNumber(cAdd(this.x[0], other), this.x[1], this.d);
        }
        QFNumber.checkD(this, other);
        return new QFNumber(
            cAdd(this.x[0], other.x[0]),
            cAdd(this.x[1], other.x[1]),
            this.d);
    }

    sub(other: QFNumber | number): QFNumber {
        if (typeof other === 'number') {
            return new QFNumber(cSub(this.x[0], other), this.x[1], this.d);
        }
        QFNumber.checkD(this, other);
        return new QFNumber(
            cSub(this.x[0], other.x[0]),
            cSub(this.x[1], other.x[1]),
            this.d);
    }

    mul(other: QFNumber | number): QFNumber {
        if (typeof other === 'number') {
            return new QFNumber(
                cMul(this.x[0], other), cMul(this.x[1], other), this.d);
        }
        QFNumber.checkD(this, other);
        return new QFNumber(
            cAdd(cMul(this.x[0], other.x[0]),
                cMul(cMul(this.x[1], other.x[1]), this.d)),
            cAdd(cMul(this.x[0], other.x[1]),
                cMul(this.x[1], other.x[0])),
            this.d);
    }

    div(other: QFNumber | number): QFNumber {
        if (typeof other === 'number') {
            return new QFNumber(
                cDiv(this.x[0], other), cDiv(this.x[1], other), this.d);
        }
        QFNumber.checkD(this, other);
        const denom = cSub(cMul(other.x[0], other.x[0]),
            cMul(cMul(other.x[1], other.x[1]), this.d));
        const numer0 = cSub(cMul(this.x[0], other.x[0]),
            cMul(cMul(this.x[1], other.x[1]), this.d));
        const numer1 = cSub(cMul(this.x[1], other.x[0]),
            cMul(this.x[0], other.x[1]));
        return new QFNumber(cDiv(numer0, denom), cDiv(numer1, denom), this.d);
    }

    // The port of C++ operator-(T const& s, QFNumber const& q).
    static scalarSub(s: number, q: QFNumber): QFNumber {
        return new QFNumber(cSub(s, q.x[0]), cNeg(q.x[1]), q.d);
    }

    // The port of C++ operator/(T const& s, QFNumber const& q).
    static scalarDiv(s: number, q: QFNumber): QFNumber {
        const denom = cSub(cMul(q.x[0], q.x[0]),
            cMul(cMul(q.x[1], q.x[1]), q.d));
        const x0 = cDiv(cMul(q.x[0], s), denom);
        const x1 = cNeg(cDiv(cMul(q.x[1], s), denom));
        return new QFNumber(x0, x1, q.d);
    }

    // Comparisons between numbers of a quadratic field must occur only when
    // the d-values are the same. To trap mismatches, read the comments about
    // assertOnMismatchedD at the beginning of this file.
    equals(other: QFNumber): boolean {
        QFNumber.checkD(this, other);
        if (this.d === 0 || cEqual(this.x[1], other.x[1])) {
            return cEqual(this.x[0], other.x[0]);
        } else if (cGreaterThan(this.x[1], other.x[1])) {
            if (cGreaterThanEqual(this.x[0], other.x[0])) {
                return false;
            } else { // this.x[0] < other.x[0]
                const diff = this.sub(other);
                return cEqual(cMul(diff.x[0], diff.x[0]),
                    cMul(cMul(diff.x[1], diff.x[1]), diff.d));
            }
        } else { // this.x[1] < other.x[1]
            if (cLessThanEqual(this.x[0], other.x[0])) {
                return false;
            } else { // this.x[0] > other.x[0]
                const diff = this.sub(other);
                return cEqual(cMul(diff.x[0], diff.x[0]),
                    cMul(cMul(diff.x[1], diff.x[1]), diff.d));
            }
        }
    }

    notEquals(other: QFNumber): boolean {
        return !this.equals(other);
    }

    lessThan(other: QFNumber): boolean {
        QFNumber.checkD(this, other);
        if (this.d === 0 || cEqual(this.x[1], other.x[1])) {
            return cLessThan(this.x[0], other.x[0]);
        } else if (cGreaterThan(this.x[1], other.x[1])) {
            if (cGreaterThanEqual(this.x[0], other.x[0])) {
                return false;
            } else { // this.x[0] < other.x[0]
                const diff = this.sub(other);
                return cGreaterThan(cMul(diff.x[0], diff.x[0]),
                    cMul(cMul(diff.x[1], diff.x[1]), diff.d));
            }
        } else { // this.x[1] < other.x[1]
            if (cLessThanEqual(this.x[0], other.x[0])) {
                return true;
            } else { // this.x[0] > other.x[0]
                const diff = this.sub(other);
                return cLessThan(cMul(diff.x[0], diff.x[0]),
                    cMul(cMul(diff.x[1], diff.x[1]), diff.d));
            }
        }
    }

    greaterThan(other: QFNumber): boolean {
        return other.lessThan(this);
    }

    lessThanEqual(other: QFNumber): boolean {
        return !other.lessThan(this);
    }

    greaterThanEqual(other: QFNumber): boolean {
        return !this.lessThan(other);
    }
}
