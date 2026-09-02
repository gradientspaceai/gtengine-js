// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) PolynomialRoot.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Support for computing roots of polynomials of degrees 1, 2, 3, or 4.
//
// Port notes:
//   - The upstream comparison operators become the methods 'equals' and
//     'lessThan'; both compare only the root estimate x, as upstream (the
//     multiplicity m plays no part in the ordering).
//   - 'PolynomialRootBisect' has the output reference parameters xMin and
//     xMax; the port takes them as values and returns { xMin, xMax }.
//   - The upstream static_assert restricting T to 'float' or 'double' is the
//     port's implicit 'number'.

// The port of struct PolynomialRoot<T>.
export class PolynomialRoot {
    // x is the root estimate and m is the multiplicity of x. The object is
    // invalid when m is 0.
    x: number;
    m: number;

    constructor(x: number = 0, m: number = 0) {
        this.x = x;
        this.m = m;
    }

    equals(other: PolynomialRoot): boolean {
        return this.x === other.x;
    }

    lessThan(other: PolynomialRoot): boolean {
        return this.x < other.x;
    }
}

// Compute a tight interval [xMin,xMax] for a root to the polynomial F(x).
// The inputs signFMin and signFMax are in {-1,1} and are the theoretical
// signs of F(xMin) and F(xMax) for the initial xMin and xMax. They are
// required to have opposite signs. Bisection is performed using
// floating-point arithmetic for speed.
export function polynomialRootBisect(F: (x: number) => number,
    signFMin: number, signFMax: number, xMin: number, xMax: number):
    { xMin: number, xMax: number } {
    const zero = 0;
    const fMin = F(xMin);
    const trueSignFMin = (fMin > zero ? +1 : (fMin < zero ? -1 : 0));
    if (trueSignFMin !== signFMin) {
        // Floating-point rounding errors prevent the correct classification
        // of the multiplicity of roots.
        return { xMin: xMin, xMax: xMin };
    }

    const fMax = F(xMax);
    const trueSignFMax = (fMax > zero ? +1 : (fMax < zero ? -1 : 0));
    if (trueSignFMax !== signFMax) {
        // Floating-point rounding errors prevent the correct classification
        // of the multiplicity of roots.
        return { xMin: xMax, xMax: xMax };
    }

    // The signs are correct for bisection. The iteration algorithm terminates
    // when the function value at the midpoint is 0. Or it terminates when the
    // midpoint of the current interval equals one of the interval endpoints,
    // at which time the interval endpoints are consecutive floating-point
    // numbers. The upper bound maxBisections is sufficiently large to ensure
    // the loop terminates, but the typical number of iterations is much
    // smaller.
    const maxBisections = 4096;
    for (let iteration = 1; iteration < maxBisections; ++iteration) {
        const x = 0.5 * (xMin + xMax);
        const f = F(x);

        if (x === xMin || x === xMax) {
            // The floating-point numbers xMin and xMax are consecutive, in
            // which case subdivision cannot produce a floating-point number
            // between them. Return the bounding interval to the caller for
            // further processing.
            return { xMin: xMin, xMax: xMax };
        }

        const signF = (f > zero ? +1 : (f < zero ? -1 : 0));
        if (signF === 0) {
            // The function is exactly zero and a root is found.
            return { xMin: x, xMax: x };
        }

        // Update the correct endpoint to the midpoint.
        if (signF === signFMin) {
            xMin = x;
        } else { // signF === signFMax
            xMax = x;
        }
    }

    return { xMin: xMin, xMax: xMax };
}
