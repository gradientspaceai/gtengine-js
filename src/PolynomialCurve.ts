// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) PolynomialCurve.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A parametric curve X(t) whose N components are each a polynomial in t.
// The derivatives of the components are polynomials themselves, computed
// once when a component polynomial is set, so evaluation of the jet is just
// polynomial evaluation.
//
// Port notes (following ParametricCurve/BezierCurve/TCBSplineCurve):
// - Upstream 'template <int32_t N, typename Real>' becomes a runtime
//   dimension passed as the first constructor argument.
// - The two C++ constructors become TypeScript constructor overloads: the
//   'components' argument is optional. The default form creates a curve with
//   all components the constant zero (all degree-0 polynomials), as upstream.
//   C++ enforces components.size() == N at compile time via std::array; the
//   port asserts it at runtime.
// - The component polynomials are copied in and out is by reference to the
//   stored objects (the port of the const-reference accessors), matching the
//   Polynomial1 getCoefficients() precedent: do not modify the returned
//   polynomial; use setPolynomial or clone() instead.
//
// Deviation from upstream: the upstream constructors never set
// mConstructed = true, so 'operator bool' reports failure for a successfully
// built curve. The port sets it (the TCBSplineCurve precedent), so
// isConstructed() is meaningful. This does not affect evaluate, which (like
// TCBSplineCurve, unlike BezierCurve/BSplineCurve) does not test the flag.

import { logAssert } from './Logger';
import { ParametricCurve } from './ParametricCurve';
import { Polynomial1 } from './Polynomial1';
import { Vector } from './Vector';

export class PolynomialCurve extends ParametricCurve {
    protected mPolynomial: Polynomial1[];
    protected mDer1Polynomial: Polynomial1[];
    protected mDer2Polynomial: Polynomial1[];
    protected mDer3Polynomial: Polynomial1[];

    // Construction. The default form creates a polynomial curve with all
    // components set to the constant zero (all degree-0 polynomials). You
    // can set these to other polynomials using setPolynomial.
    constructor(dimension: number, tmin: number, tmax: number,
        components?: readonly Polynomial1[]) {
        super(dimension, tmin, tmax);

        logAssert(dimension > 0 && Number.isInteger(dimension),
            'Invalid dimension.');

        this.mPolynomial = new Array<Polynomial1>(dimension);
        this.mDer1Polynomial = new Array<Polynomial1>(dimension);
        this.mDer2Polynomial = new Array<Polynomial1>(dimension);
        this.mDer3Polynomial = new Array<Polynomial1>(dimension);
        for (let i = 0; i < dimension; ++i) {
            // The port of default construction of std::array<Polynomial1, N>:
            // each component is the degree-0 polynomial with coefficient 0.
            this.mPolynomial[i] = new Polynomial1(0);
            this.mDer1Polynomial[i] = new Polynomial1(0);
            this.mDer2Polynomial[i] = new Polynomial1(0);
            this.mDer3Polynomial[i] = new Polynomial1(0);
        }

        if (components !== undefined) {
            logAssert(components.length === dimension,
                'Invalid number of components.');

            for (let i = 0; i < dimension; ++i) {
                this.setPolynomial(i, components[i]);
            }
        }

        this.mConstructed = true;
    }

    // Member access. The polynomial is copied in and its first-, second- and
    // third-order derivatives are computed and stored.
    setPolynomial(i: number, poly: Polynomial1): void {
        logAssert(i >= 0 && i < this.mDimension && Number.isInteger(i),
            'Invalid index.');

        this.mPolynomial[i] = poly.clone();
        this.mDer1Polynomial[i] = this.mPolynomial[i].getDerivative();
        this.mDer2Polynomial[i] = this.mDer1Polynomial[i].getDerivative();
        this.mDer3Polynomial[i] = this.mDer2Polynomial[i].getDerivative();
    }

    // The returned polynomials are the stored objects (the port of the
    // upstream const-reference accessors); do not modify them.
    getPolynomial(i: number): Polynomial1 {
        return this.mPolynomial[i];
    }

    getDer1Polynomial(i: number): Polynomial1 {
        return this.mDer1Polynomial[i];
    }

    getDer2Polynomial(i: number): Polynomial1 {
        return this.mDer2Polynomial[i];
    }

    getDer3Polynomial(i: number): Polynomial1 {
        return this.mDer3Polynomial[i];
    }

    // Evaluation of the curve. The function supports derivative calculation
    // through order 3; that is, order <= 3 is required. If you want only the
    // position, pass in order of 0. If you want the position and first
    // derivative, pass in order of 1, and so on. The output array 'jet' must
    // have enough storage to support the maximum order. The values are
    // ordered as: position, first derivative, second derivative, third
    // derivative.
    override evaluate(t: number, order: number, jet: Vector[]): void {
        const n = this.mDimension;

        for (let i = 0; i < n; ++i) {
            jet[0].values[i] = this.mPolynomial[i].evaluate(t);
        }

        if (order >= 1) {
            for (let i = 0; i < n; ++i) {
                jet[1].values[i] = this.mDer1Polynomial[i].evaluate(t);
            }

            if (order >= 2) {
                for (let i = 0; i < n; ++i) {
                    jet[2].values[i] = this.mDer2Polynomial[i].evaluate(t);
                }

                // Upstream tests 'order == 3' here rather than 'order >= 3',
                // so an order larger than the documented maximum of 3 leaves
                // jet[3] untouched. The precondition is order <= 3, so the
                // quirk is preserved.
                if (order === 3) {
                    for (let i = 0; i < n; ++i) {
                        jet[3].values[i] = this.mDer3Polynomial[i].evaluate(t);
                    }
                }
            }
        }
    }
}
