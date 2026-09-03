// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Polynomial1.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A polynomial p(t) = c[0] + c[1] * t + ... + c[degree] * t^degree with
// real-valued coefficients, stored in order of increasing power of t.
//
// Port notes:
// - The two C++ constructors (degree, and std::initializer_list) are
//   ambiguous in TypeScript, so the initializer-list form becomes the static
//   factory fromCoefficients(...) (constructor-factory precedent). The
//   remaining constructor takes the degree and zero-fills, as upstream.
// - operator[] -> get(i) / set(i, value) (Vector precedent).
// - Comparison operators -> equals, notEquals, lessThan, lessThanOrEqual,
//   greaterThan, greaterThanOrEqual. As with std::vector, the order is
//   lexicographic on the coefficient arrays, and a proper prefix compares
//   less than the longer array.
// - operator()(t) -> evaluate(t).
// - The free arithmetic operators become instance methods returning new
//   objects: unary operator- -> negate; operator+ -> add; operator- -> sub;
//   operator* -> mul; operator/ -> div. The polynomial-and-scalar forms are
//   overloads of add/sub/mul (the commuted scalar forms are the same
//   operation); operator-(scalar, p) -> subFrom(scalar) and
//   operator/(p, scalar) -> div(scalar). The compound assignment operators
//   (+=, -=, *=, /=) are omitted: upstream defines each as the corresponding
//   binary operation followed by assignment, so p = p.add(q) is the port.
// - Divide writes the quotient and remainder to reference parameters
//   upstream; the port returns { quotient, remainder }.
// - The free functions GreatestCommonDivisor and SquareFreeFactorization
//   become greatestCommonDivisor and squareFreeFactorization;
//   SquareFreeFactorization appends to an output vector upstream, and the
//   port returns the array of factors. Upstream's do-while loop in
//   SquareFreeFactorization can spin forever for floating-point Real; the
//   port caps the iterations and throws instead (see the comment there).

import { logAssert, logError } from './Logger.js';

export class Polynomial1 {
    // The class is designed so that mCoefficient.length >= 1.
    protected mCoefficient: number[];

    // Create a polynomial of the specified degree with all coefficients set
    // to zero (to ensure initialization). You are responsible for setting
    // the coefficients, presumably with the degree-term set to a nonzero
    // number.
    constructor(degree: number = 0) {
        logAssert(degree >= 0 && Number.isInteger(degree), 'Invalid degree.');
        this.mCoefficient = new Array<number>(degree + 1).fill(0);
    }

    // Create a polynomial from coefficients listed in order of increasing
    // power of t. The degree is the number of values minus 1, but then
    // adjusted so that coefficient[degree] is not zero (unless all the
    // values are zero). This is the port of the initializer-list
    // constructor, for which C++ guarantees values.size() > 0.
    static fromCoefficients(values: number[]): Polynomial1 {
        logAssert(values.length > 0, 'Invalid number of coefficients.');
        const result = new Polynomial1(values.length - 1);
        for (let i = 0; i < values.length; ++i) {
            result.mCoefficient[i] = values[i];
        }
        result.eliminateLeadingZeros();
        return result;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Polynomial1 {
        const result = new Polynomial1(this.getDegree());
        for (let i = 0; i < this.mCoefficient.length; ++i) {
            result.mCoefficient[i] = this.mCoefficient[i];
        }
        return result;
    }

    // Support for partial construction, where the default constructor is
    // used when the degree is not yet known. Upstream calls std::vector
    // resize, which preserves the existing coefficients and value-initializes
    // (that is, zeros) any new ones.
    setDegree(degree: number): void {
        logAssert(degree >= 0 && Number.isInteger(degree), 'Invalid degree.');
        const size = degree + 1;
        while (this.mCoefficient.length < size) {
            this.mCoefficient.push(0);
        }
        this.mCoefficient.length = size;
    }

    // Set all coefficients to the specified value.
    setCoefficients(value: number): void {
        this.mCoefficient.fill(value);
    }

    // Member access.
    getDegree(): number {
        // By design, mCoefficient.length > 0.
        return this.mCoefficient.length - 1;
    }

    // The internal coefficient array (the port of the const-reference
    // accessor). Do not modify the returned array; use set(i, value) or
    // clone().
    getCoefficients(): number[] {
        return this.mCoefficient;
    }

    get(i: number): number {
        return this.mCoefficient[i];
    }

    set(i: number, value: number): void {
        this.mCoefficient[i] = value;
    }

    // Comparisons (lexicographic on the coefficient arrays).
    equals(p: Polynomial1): boolean {
        if (this.mCoefficient.length !== p.mCoefficient.length) {
            return false;
        }
        for (let i = 0; i < this.mCoefficient.length; ++i) {
            if (this.mCoefficient[i] !== p.mCoefficient[i]) {
                return false;
            }
        }
        return true;
    }

    notEquals(p: Polynomial1): boolean {
        return !this.equals(p);
    }

    lessThan(p: Polynomial1): boolean {
        const n = Math.min(this.mCoefficient.length, p.mCoefficient.length);
        for (let i = 0; i < n; ++i) {
            if (this.mCoefficient[i] < p.mCoefficient[i]) {
                return true;
            }
            if (this.mCoefficient[i] > p.mCoefficient[i]) {
                return false;
            }
        }
        return this.mCoefficient.length < p.mCoefficient.length;
    }

    lessThanOrEqual(p: Polynomial1): boolean {
        return !p.lessThan(this);
    }

    greaterThan(p: Polynomial1): boolean {
        return p.lessThan(this);
    }

    greaterThanOrEqual(p: Polynomial1): boolean {
        return !this.lessThan(p);
    }

    // Evaluate the polynomial using Horner's method.
    evaluate(t: number): number {
        let i = this.mCoefficient.length;
        let result = this.mCoefficient[--i];
        for (--i; i >= 0; --i) {
            result *= t;
            result += this.mCoefficient[i];
        }
        return result;
    }

    // Compute the derivative of the polynomial.
    getDerivative(): Polynomial1 {
        const degree = this.getDegree();
        if (degree > 0) {
            const result = new Polynomial1(degree - 1);
            for (let i0 = 0, i1 = 1; i0 < degree; ++i0, ++i1) {
                result.mCoefficient[i0] = this.mCoefficient[i1] * i1;
            }
            return result;
        } else {
            const result = new Polynomial1(0);
            result.mCoefficient[0] = 0;
            return result;
        }
    }

    // Inversion (invpoly[i] = poly[degree - i] for 0 <= i <= degree).
    getInversion(): Polynomial1 {
        const degree = this.getDegree();
        const result = new Polynomial1(degree);
        for (let i = 0; i <= degree; ++i) {
            result.mCoefficient[i] = this.mCoefficient[degree - i];
        }
        return result;
    }

    // Translation. If 'this' is p(t), return p(t - t0).
    getTranslation(t0: number): Polynomial1 {
        const factor = Polynomial1.fromCoefficients([-t0, 1]);  // f(t) = t - t0
        const degree = this.getDegree();
        let result = Polynomial1.fromCoefficients([this.mCoefficient[degree]]);
        for (let i = 1, j = degree - 1; i <= degree; ++i, --j) {
            result = factor.mul(result).add(this.mCoefficient[j]);
        }
        return result;
    }

    // Eliminate any leading zeros in the polynomial, except in the case the
    // degree is 0 and the coefficient is 0. The elimination is necessary
    // when arithmetic operations cause a decrease in the degree of the
    // result. For example, (1 + x + x^2) + (1 + 2*x - x^2) = (2 + 3*x). The
    // inputs both have degree 2, so the result is created with degree 2.
    // After the addition we find that the degree is in fact 1 and resize the
    // array of coefficients. This function is called internally by the
    // arithmetic operations, but it is exposed in the public interface in
    // case you need it for your own purposes.
    eliminateLeadingZeros(): void {
        const size = this.mCoefficient.length;
        if (size > 1) {
            let leading: number;
            for (leading = size - 1; leading > 0; --leading) {
                if (this.mCoefficient[leading] !== 0) {
                    break;
                }
            }

            this.mCoefficient.length = ++leading;
        }
    }

    // If 'this' is P(t) and the divisor is D(t) with
    // degree(P) >= degree(D), then P(t) = Q(t)*D(t) + R(t) where Q(t) is the
    // quotient with degree(Q) = degree(P) - degree(D) and R(t) is the
    // remainder with degree(R) < degree(D). If this routine is called with
    // degree(P) < degree(D), then Q = 0 and R = P are returned.
    divide(divisor: Polynomial1): { quotient: Polynomial1; remainder: Polynomial1 } {
        const quotient = new Polynomial1(0);
        const remainder = new Polynomial1(0);
        const divisorDegree = divisor.getDegree();
        const quotientDegree = this.getDegree() - divisorDegree;
        if (quotientDegree >= 0) {
            quotient.setDegree(quotientDegree);

            // Temporary storage for the remainder.
            const tmp = this.clone();

            // Do the division using the Euclidean algorithm.
            const inv = 1 / divisor.mCoefficient[divisorDegree];
            for (let i = quotientDegree; i >= 0; --i) {
                let j = divisorDegree + i;
                quotient.mCoefficient[i] = inv * tmp.mCoefficient[j];
                for (j--; j >= i; j--) {
                    tmp.mCoefficient[j] -=
                        quotient.mCoefficient[i] * divisor.mCoefficient[j - i];
                }
            }

            // Calculate the correct degree for the remainder.
            if (divisorDegree >= 1) {
                let remainderDegree = divisorDegree - 1;
                while (remainderDegree > 0 && tmp.mCoefficient[remainderDegree] === 0) {
                    --remainderDegree;
                }

                remainder.setDegree(remainderDegree);
                for (let i = 0; i <= remainderDegree; ++i) {
                    remainder.mCoefficient[i] = tmp.mCoefficient[i];
                }
            } else {
                remainder.setDegree(0);
                remainder.mCoefficient[0] = 0;
            }
            return { quotient, remainder };
        } else {
            quotient.setDegree(0);
            quotient.mCoefficient[0] = 0;
            return { quotient, remainder: this.clone() };
        }
    }

    // Scale the polynomial so the highest-degree term has coefficient 1.
    makeMonic(): void {
        this.eliminateLeadingZeros();
        const last = this.mCoefficient.length - 1;
        if (this.mCoefficient[last] !== 1) {
            const degree = this.getDegree();
            const invLeading = 1 / this.mCoefficient[last];
            this.mCoefficient[last] = 1;
            for (let i = 0; i < degree; ++i) {
                this.mCoefficient[i] *= invLeading;
            }
        }
    }

    // Unary operator-.
    negate(): Polynomial1 {
        const degree = this.getDegree();
        const result = new Polynomial1(degree);
        for (let i = 0; i <= degree; ++i) {
            result.mCoefficient[i] = -this.mCoefficient[i];
        }
        return result;
    }

    // operator+(p0, p1) and operator+(p, scalar) [= operator+(scalar, p)].
    add(p: Polynomial1): Polynomial1;
    add(scalar: number): Polynomial1;
    add(arg: Polynomial1 | number): Polynomial1 {
        if (typeof arg === 'number') {
            const degree = this.getDegree();
            const result = new Polynomial1(degree);
            result.mCoefficient[0] = this.mCoefficient[0] + arg;
            for (let i = 1; i <= degree; ++i) {
                result.mCoefficient[i] = this.mCoefficient[i];
            }
            return result;
        }

        const p0Degree = this.getDegree(), p1Degree = arg.getDegree();
        let i: number;
        if (p0Degree >= p1Degree) {
            const result = new Polynomial1(p0Degree);
            for (i = 0; i <= p1Degree; ++i) {
                result.mCoefficient[i] = this.mCoefficient[i] + arg.mCoefficient[i];
            }
            for (/**/; i <= p0Degree; ++i) {
                result.mCoefficient[i] = this.mCoefficient[i];
            }
            result.eliminateLeadingZeros();
            return result;
        } else {
            const result = new Polynomial1(p1Degree);
            for (i = 0; i <= p0Degree; ++i) {
                result.mCoefficient[i] = this.mCoefficient[i] + arg.mCoefficient[i];
            }
            for (/**/; i <= p1Degree; ++i) {
                result.mCoefficient[i] = arg.mCoefficient[i];
            }
            result.eliminateLeadingZeros();
            return result;
        }
    }

    // operator-(p0, p1) and operator-(p, scalar).
    sub(p: Polynomial1): Polynomial1;
    sub(scalar: number): Polynomial1;
    sub(arg: Polynomial1 | number): Polynomial1 {
        if (typeof arg === 'number') {
            const degree = this.getDegree();
            const result = new Polynomial1(degree);
            result.mCoefficient[0] = this.mCoefficient[0] - arg;
            for (let i = 1; i <= degree; ++i) {
                result.mCoefficient[i] = this.mCoefficient[i];
            }
            return result;
        }

        const p0Degree = this.getDegree(), p1Degree = arg.getDegree();
        let i: number;
        if (p0Degree >= p1Degree) {
            const result = new Polynomial1(p0Degree);
            for (i = 0; i <= p1Degree; ++i) {
                result.mCoefficient[i] = this.mCoefficient[i] - arg.mCoefficient[i];
            }
            for (/**/; i <= p0Degree; ++i) {
                result.mCoefficient[i] = this.mCoefficient[i];
            }
            result.eliminateLeadingZeros();
            return result;
        } else {
            const result = new Polynomial1(p1Degree);
            for (i = 0; i <= p0Degree; ++i) {
                result.mCoefficient[i] = this.mCoefficient[i] - arg.mCoefficient[i];
            }
            for (/**/; i <= p1Degree; ++i) {
                result.mCoefficient[i] = -arg.mCoefficient[i];
            }
            result.eliminateLeadingZeros();
            return result;
        }
    }

    // operator-(scalar, p), that is, scalar - this.
    subFrom(scalar: number): Polynomial1 {
        const degree = this.getDegree();
        const result = new Polynomial1(degree);
        result.mCoefficient[0] = scalar - this.mCoefficient[0];
        for (let i = 1; i <= degree; ++i) {
            result.mCoefficient[i] = -this.mCoefficient[i];
        }
        return result;
    }

    // operator*(p0, p1) and operator*(p, scalar) [= operator*(scalar, p)].
    mul(p: Polynomial1): Polynomial1;
    mul(scalar: number): Polynomial1;
    mul(arg: Polynomial1 | number): Polynomial1 {
        if (typeof arg === 'number') {
            const degree = this.getDegree();
            const result = new Polynomial1(degree);
            for (let i = 0; i <= degree; ++i) {
                result.mCoefficient[i] = arg * this.mCoefficient[i];
            }
            return result;
        }

        const p0Degree = this.getDegree(), p1Degree = arg.getDegree();
        const result = new Polynomial1(p0Degree + p1Degree);
        result.setCoefficients(0);
        for (let i0 = 0; i0 <= p0Degree; ++i0) {
            for (let i1 = 0; i1 <= p1Degree; ++i1) {
                result.mCoefficient[i0 + i1] += this.mCoefficient[i0] * arg.mCoefficient[i1];
            }
        }
        return result;
    }

    // operator/(p, scalar).
    div(scalar: number): Polynomial1 {
        logAssert(scalar !== 0, 'Division by zero.');

        const degree = this.getDegree();
        const invScalar = 1 / scalar;
        const result = new Polynomial1(degree);
        for (let i = 0; i <= degree; ++i) {
            result.mCoefficient[i] = invScalar * this.mCoefficient[i];
        }
        return result;
    }
}

// Compute the greatest common divisor of two polynomials. The returned
// polynomial has leading coefficient 1 (except when zero-valued polynomials
// are passed to the function).
export function greatestCommonDivisor(p0: Polynomial1, p1: Polynomial1): Polynomial1 {
    // The numerator should be the polynomial of larger degree.
    let a: Polynomial1, b: Polynomial1;
    if (p0.getDegree() >= p1.getDegree()) {
        a = p0.clone();
        b = p1.clone();
    } else {
        a = p1.clone();
        b = p0.clone();
    }

    const zero = Polynomial1.fromCoefficients([0]);
    if (a.equals(zero) || b.equals(zero)) {
        return a.notEquals(zero) ? a : zero;
    }

    // Make the polynomials monic to keep the coefficients a reasonable size
    // when computing with floating-point coefficients.
    a.makeMonic();
    b.makeMonic();

    for (; ;) {
        const { remainder: r } = a.divide(b);
        if (r.notEquals(zero)) {
            // a = q * b + r, so gcd(a,b) = gcd(b,r).
            a = b;
            b = r;
            b.makeMonic();
        } else {
            b.makeMonic();
            break;
        }
    }

    return b;
}

// Factor f = factor[0]*factor[1]^2*factor[2]^3*...*factor[n-1]^n according to
// the square-free factorization algorithm
// https://en.wikipedia.org/wiki/Square-free_polynomial
//
// The algorithm is exact only for exact arithmetic. With floating-point
// coefficients, greatestCommonDivisor tests remainders against exactly zero,
// so a rounding error of a few ULP in a division can make it report the
// constant 1 for a pair of polynomials whose true GCD is nonconstant. When
// that happens the degree of b never decreases and upstream's do-while loop
// never exits. Example: f = (t-2)^2 * (t+2) * (t+3) * (t-4), for which
// gcd(f, f') rounds to t - 1.9999999999999998; every later remainder is
// O(1e-14) rather than 0 and the loop spins forever with a = 1 and d growing
// without bound. In exact arithmetic the loop runs at most degree(f) times
// (one iteration per multiplicity), so the port caps the iterations at
// degree(f) + 1 and throws rather than hanging the caller.
export function squareFreeFactorization(f: Polynomial1): Polynomial1[] {
    // In the calls to divide(...), we know that the divisor exactly divides
    // the numerator, so r = 0 after all such calls.
    const factors: Polynomial1[] = [];
    const fder = f.getDerivative();

    let a = greatestCommonDivisor(f, fder);
    let b = f.divide(a).quotient;      // b = f / a
    let c = fder.divide(a).quotient;   // c = fder / a
    let d = c.sub(b.getDerivative());

    const maxIterations = f.getDegree() + 1;
    let iteration = 0;
    do {
        if (++iteration > maxIterations) {
            logError('The square-free factorization did not converge. The ' +
                'floating-point GCD computation is numerically unstable for ' +
                'this polynomial.');
        }
        a = greatestCommonDivisor(b, d);
        factors.push(a);
        b = b.divide(a).quotient;      // b = b / a
        c = d.divide(a).quotient;      // c = d / a
        d = c.sub(b.getDerivative());
    } while (b.getDegree() > 0);

    return factors;
}
