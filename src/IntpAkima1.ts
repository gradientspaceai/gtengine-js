// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpAkima1.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Akima interpolation is described in
// https://en.wikipedia.org/wiki/Akima_spline
//
// Port notes: upstream nested class IntpAkima1::Polynomial is exported as
// IntpAkima1Polynomial (global export uniqueness) so derived classes in
// other files can use it. Its 'operator[]' coefficient reference accessor is
// ported as getCoeff(i)/setCoeff(i, value). The overloaded 'operator()(x)'
// and 'operator()(order, x)' of both the interpolator and the polynomial are
// ported as the overloaded method 'evaluate'. The virtual
// 'Lookup(x, index&, dx&)' is ported as the abstract method 'lookup(x)'
// returning { index, dx }. ComputeDerivative(Real* slope) receives a pointer
// into a slope array; the port takes the array and the offset of the
// four-element window.

import { logAssert } from './Logger.js';

// P(x) = c[0] + c[1]*x + c[2]*x^2 + c[3]*x^3
export class IntpAkima1Polynomial {
    private mCoeff: number[];

    constructor() {
        this.mCoeff = [0, 0, 0, 0];
    }

    getCoeff(i: number): number {
        return this.mCoeff[i];
    }

    setCoeff(i: number, value: number): void {
        this.mCoeff[i] = value;
    }

    // evaluate(x) evaluates the polynomial. evaluate(order, x) evaluates the
    // derivative of the specified order (order zero is the function itself);
    // the returned value is zero when order >= 4.
    evaluate(x: number): number;
    evaluate(order: number, x: number): number;
    evaluate(arg0: number, arg1?: number): number {
        if (arg1 === undefined) {
            const x = arg0;
            return this.mCoeff[0] + x * (this.mCoeff[1] + x * (this.mCoeff[2] + x * this.mCoeff[3]));
        }

        const order = arg0;
        const x = arg1;
        switch (order) {
            case 0:
                return this.mCoeff[0] + x * (this.mCoeff[1] + x * (this.mCoeff[2] + x * this.mCoeff[3]));
            case 1:
                return this.mCoeff[1] + x * (2 * this.mCoeff[2] + x * 3 * this.mCoeff[3]);
            case 2:
                return 2 * this.mCoeff[2] + x * 6 * this.mCoeff[3];
            case 3:
                return 6 * this.mCoeff[3];
        }

        return 0;
    }
}

export abstract class IntpAkima1 {
    protected mQuantity: number;
    protected mF: readonly number[];
    protected mPoly: IntpAkima1Polynomial[];

    // Construction (abstract base class). The samples F are aliased, not
    // copied.
    protected constructor(quantity: number, F: readonly number[]) {
        // At least three data points are needed to construct the estimates
        // of the boundary derivatives.
        logAssert(quantity >= 3, 'Invalid input to IntpAkima1 constructor.');

        this.mQuantity = quantity;
        this.mF = F;
        this.mPoly = [];
        for (let i = 0; i < quantity - 1; ++i) {
            this.mPoly.push(new IntpAkima1Polynomial());
        }
    }

    // Member access.
    getQuantity(): number {
        return this.mQuantity;
    }

    getF(): readonly number[] {
        return this.mF;
    }

    abstract getXMin(): number;

    abstract getXMax(): number;

    // Evaluate the function and its derivatives. The functions clamp the
    // inputs to xmin <= x <= xmax. evaluate(x) is for function evaluation.
    // evaluate(order, x) is for function or derivative evaluations; the
    // 'order' argument is the order of the derivative or zero for the
    // function itself.
    evaluate(x: number): number;
    evaluate(order: number, x: number): number;
    evaluate(arg0: number, arg1?: number): number {
        if (arg1 === undefined) {
            let x = arg0;
            x = Math.min(Math.max(x, this.getXMin()), this.getXMax());
            const { index, dx } = this.lookup(x);
            return this.mPoly[index].evaluate(dx);
        }

        const order = arg0;
        let x = arg1;
        x = Math.min(Math.max(x, this.getXMin()), this.getXMax());
        const { index, dx } = this.lookup(x);
        return this.mPoly[index].evaluate(order, dx);
    }

    // Estimate the derivative from the four consecutive slopes
    // slope[offset], ..., slope[offset + 3].
    protected computeDerivative(slope: readonly number[], offset: number): number {
        if (slope[offset + 1] !== slope[offset + 2]) {
            if (slope[offset] !== slope[offset + 1]) {
                if (slope[offset + 2] !== slope[offset + 3]) {
                    const ad0 = Math.abs(slope[offset + 3] - slope[offset + 2]);
                    const ad1 = Math.abs(slope[offset] - slope[offset + 1]);
                    return (ad0 * slope[offset + 1] + ad1 * slope[offset + 2]) / (ad0 + ad1);
                }
                else {
                    return slope[offset + 2];
                }
            }
            else {
                if (slope[offset + 2] !== slope[offset + 3]) {
                    return slope[offset + 1];
                }
                else {
                    return 0.5 * (slope[offset + 1] + slope[offset + 2]);
                }
            }
        }
        else {
            return slope[offset + 1];
        }
    }

    protected abstract lookup(x: number): { index: number, dx: number };
}
