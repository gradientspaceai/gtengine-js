// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpAkimaUniform2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Akima interpolation is described in
// https://en.wikipedia.org/wiki/Akima_spline
// The interpolator is for uniformly spaced (x,y)-values. The input samples
// F must be stored in row-major order to represent f(x,y); that is,
// F[c + xBound*r] corresponds to f(x,y), where c is the index corresponding
// to x and r is the index corresponding to y.
//
// Port notes: the nested private Polynomial class is module-private (not
// exported). Its coefficient reference accessor 'A(ix, iy)' is ported as
// getA(ix, iy)/setA(ix, iy, value), and its overloaded 'operator()' is
// ported as the overloaded method 'evaluate'. The interpolator's overloaded
// 'operator()(x, y)' and 'operator()(xOrder, yOrder, x, y)' are ported as
// the overloaded method 'evaluate'. XLookup/YLookup with output references
// are ported as methods returning { index, delta }. C++ 'array[iy][ix]'
// accesses on Array2 become 'array.get(ix, iy)' per the Array2 port, and
// ComputeDerivative(Real* slope) with pointer offsets receives the four
// window values directly.

import { Array2 } from './Array2.js';
import { logAssert } from './Logger.js';

// P(x,y) = (1,x,x^2,x^3)*A*(1,y,y^2,y^3). The matrix term A[ix][iy]
// corresponds to the polynomial term x^{ix} y^{iy}.
class Polynomial {
    // mCoeff[ix][iy]
    private mCoeff: number[][];

    constructor() {
        this.mCoeff = [];
        for (let i = 0; i < 4; ++i) {
            this.mCoeff.push([0, 0, 0, 0]);
        }
    }

    getA(ix: number, iy: number): number {
        return this.mCoeff[ix][iy];
    }

    setA(ix: number, iy: number, value: number): void {
        this.mCoeff[ix][iy] = value;
    }

    evaluate(x: number, y: number): number;
    evaluate(xOrder: number, yOrder: number, x: number, y: number): number;
    evaluate(arg0: number, arg1: number, arg2?: number, arg3?: number): number {
        if (arg2 === undefined || arg3 === undefined) {
            const x = arg0;
            const y = arg1;
            const B: number[] = [0, 0, 0, 0];
            for (let i = 0; i <= 3; ++i) {
                B[i] = this.mCoeff[i][0] + y * (this.mCoeff[i][1]
                    + y * (this.mCoeff[i][2] + y * this.mCoeff[i][3]));
            }

            return B[0] + x * (B[1] + x * (B[2] + x * B[3]));
        }

        const xOrder = arg0;
        const yOrder = arg1;
        const x = arg2;
        const y = arg3;

        const xPow: number[] = [0, 0, 0, 0];
        switch (xOrder) {
            case 0:
                xPow[0] = 1;
                xPow[1] = x;
                xPow[2] = x * x;
                xPow[3] = x * x * x;
                break;
            case 1:
                xPow[0] = 0;
                xPow[1] = 1;
                xPow[2] = 2 * x;
                xPow[3] = 3 * x * x;
                break;
            case 2:
                xPow[0] = 0;
                xPow[1] = 0;
                xPow[2] = 2;
                xPow[3] = 6 * x;
                break;
            case 3:
                xPow[0] = 0;
                xPow[1] = 0;
                xPow[2] = 0;
                xPow[3] = 6;
                break;
            default:
                return 0;
        }

        const yPow: number[] = [0, 0, 0, 0];
        switch (yOrder) {
            case 0:
                yPow[0] = 1;
                yPow[1] = y;
                yPow[2] = y * y;
                yPow[3] = y * y * y;
                break;
            case 1:
                yPow[0] = 0;
                yPow[1] = 1;
                yPow[2] = 2 * y;
                yPow[3] = 3 * y * y;
                break;
            case 2:
                yPow[0] = 0;
                yPow[1] = 0;
                yPow[2] = 2;
                yPow[3] = 6 * y;
                break;
            case 3:
                yPow[0] = 0;
                yPow[1] = 0;
                yPow[2] = 0;
                yPow[3] = 6;
                break;
            default:
                return 0;
        }

        let p = 0;
        for (let iy = 0; iy <= 3; ++iy) {
            for (let ix = 0; ix <= 3; ++ix) {
                p += this.mCoeff[ix][iy] * xPow[ix] * yPow[iy];
            }
        }

        return p;
    }
}

export class IntpAkimaUniform2 {
    private mXBound: number;
    private mYBound: number;
    private mQuantity: number;
    private mXMin: number;
    private mXMax: number;
    private mXSpacing: number;
    private mYMin: number;
    private mYMax: number;
    private mYSpacing: number;
    private mF: readonly number[];
    private mPoly: Array2<Polynomial>;

    // Construction. The samples F are aliased, not copied.
    constructor(xBound: number, yBound: number, xMin: number, xSpacing: number,
        yMin: number, ySpacing: number, F: readonly number[]) {
        // At least a 3x3 block of data points is needed to construct the
        // estimates of the boundary derivatives.
        logAssert(xBound >= 3 && yBound >= 3 && F.length >= xBound * yBound,
            'Invalid input.');
        logAssert(xSpacing > 0 && ySpacing > 0, 'Invalid input.');

        this.mXBound = xBound;
        this.mYBound = yBound;
        this.mQuantity = xBound * yBound;
        this.mXMin = xMin;
        this.mXSpacing = xSpacing;
        this.mYMin = yMin;
        this.mYSpacing = ySpacing;
        this.mF = F;

        this.mXMax = xMin + xSpacing * (xBound - 1);
        this.mYMax = yMin + ySpacing * (yBound - 1);

        this.mPoly = new Array2<Polynomial>(xBound - 1, yBound - 1);
        for (let iy = 0; iy < yBound - 1; ++iy) {
            for (let ix = 0; ix < xBound - 1; ++ix) {
                this.mPoly.set(ix, iy, new Polynomial());
            }
        }

        // Create a 2D wrapper for the 1D samples.
        const Fmap = new Array2<number>(xBound, yBound, F.slice(0, this.mQuantity));

        // Construct first-order derivatives.
        const FX = new Array2<number>(xBound, yBound);
        const FY = new Array2<number>(xBound, yBound);
        this.getFX(Fmap, FX);
        this.getFY(Fmap, FY);

        // Construct second-order derivatives.
        const FXY = new Array2<number>(xBound, yBound);
        this.getFXY(Fmap, FXY);

        // Construct polynomials.
        this.getPolynomials(Fmap, FX, FY, FXY);
    }

    // Member access.
    getXBound(): number {
        return this.mXBound;
    }

    getYBound(): number {
        return this.mYBound;
    }

    getQuantity(): number {
        return this.mQuantity;
    }

    getF(): readonly number[] {
        return this.mF;
    }

    getXMin(): number {
        return this.mXMin;
    }

    getXMax(): number {
        return this.mXMax;
    }

    getXSpacing(): number {
        return this.mXSpacing;
    }

    getYMin(): number {
        return this.mYMin;
    }

    getYMax(): number {
        return this.mYMax;
    }

    getYSpacing(): number {
        return this.mYSpacing;
    }

    // Evaluate the function and its derivatives. The functions clamp the
    // inputs to xmin <= x <= xmax and ymin <= y <= ymax. evaluate(x, y) is
    // for function evaluation. evaluate(xOrder, yOrder, x, y) is for
    // function or derivative evaluations; the xOrder argument is the order
    // of the x-derivative and the yOrder argument is the order of the
    // y-derivative. Both orders are zero to get the function value itself.
    evaluate(x: number, y: number): number;
    evaluate(xOrder: number, yOrder: number, x: number, y: number): number;
    evaluate(arg0: number, arg1: number, arg2?: number, arg3?: number): number {
        if (arg2 === undefined || arg3 === undefined) {
            let x = arg0;
            let y = arg1;
            x = Math.min(Math.max(x, this.mXMin), this.mXMax);
            y = Math.min(Math.max(y, this.mYMin), this.mYMax);
            const { index: ix, delta: dx } = this.xLookup(x);
            const { index: iy, delta: dy } = this.yLookup(y);
            return this.mPoly.get(ix, iy).evaluate(dx, dy);
        }

        const xOrder = arg0;
        const yOrder = arg1;
        let x = arg2;
        let y = arg3;
        x = Math.min(Math.max(x, this.mXMin), this.mXMax);
        y = Math.min(Math.max(y, this.mYMin), this.mYMax);
        const { index: ix, delta: dx } = this.xLookup(x);
        const { index: iy, delta: dy } = this.yLookup(y);
        return this.mPoly.get(ix, iy).evaluate(xOrder, yOrder, dx, dy);
    }

    // Support for construction.
    private getFX(F: Array2<number>, FX: Array2<number>): void {
        const slope = new Array2<number>(this.mXBound + 3, this.mYBound);
        const invDX = 1 / this.mXSpacing;
        let ix: number, iy: number;
        for (iy = 0; iy < this.mYBound; ++iy) {
            for (ix = 0; ix < this.mXBound - 1; ++ix) {
                slope.set(ix + 2, iy, (F.get(ix + 1, iy) - F.get(ix, iy)) * invDX);
            }

            slope.set(1, iy, 2 * slope.get(2, iy) - slope.get(3, iy));
            slope.set(0, iy, 2 * slope.get(1, iy) - slope.get(2, iy));
            slope.set(this.mXBound + 1, iy,
                2 * slope.get(this.mXBound, iy) - slope.get(this.mXBound - 1, iy));
            slope.set(this.mXBound + 2, iy,
                2 * slope.get(this.mXBound + 1, iy) - slope.get(this.mXBound, iy));
        }

        for (iy = 0; iy < this.mYBound; ++iy) {
            for (ix = 0; ix < this.mXBound; ++ix) {
                FX.set(ix, iy, this.computeDerivative(
                    slope.get(ix, iy), slope.get(ix + 1, iy),
                    slope.get(ix + 2, iy), slope.get(ix + 3, iy)));
            }
        }
    }

    private getFY(F: Array2<number>, FY: Array2<number>): void {
        const slope = new Array2<number>(this.mYBound + 3, this.mXBound);
        const invDY = 1 / this.mYSpacing;
        let ix: number, iy: number;
        for (ix = 0; ix < this.mXBound; ++ix) {
            for (iy = 0; iy < this.mYBound - 1; ++iy) {
                slope.set(iy + 2, ix, (F.get(ix, iy + 1) - F.get(ix, iy)) * invDY);
            }

            slope.set(1, ix, 2 * slope.get(2, ix) - slope.get(3, ix));
            slope.set(0, ix, 2 * slope.get(1, ix) - slope.get(2, ix));
            slope.set(this.mYBound + 1, ix,
                2 * slope.get(this.mYBound, ix) - slope.get(this.mYBound - 1, ix));
            slope.set(this.mYBound + 2, ix,
                2 * slope.get(this.mYBound + 1, ix) - slope.get(this.mYBound, ix));
        }

        for (ix = 0; ix < this.mXBound; ++ix) {
            for (iy = 0; iy < this.mYBound; ++iy) {
                FY.set(ix, iy, this.computeDerivative(
                    slope.get(iy, ix), slope.get(iy + 1, ix),
                    slope.get(iy + 2, ix), slope.get(iy + 3, ix)));
            }
        }
    }

    private getFXY(F: Array2<number>, FXY: Array2<number>): void {
        const xBoundM1 = this.mXBound - 1;
        const yBoundM1 = this.mYBound - 1;
        const ix0 = xBoundM1, ix1 = ix0 - 1, ix2 = ix1 - 1;
        const iy0 = yBoundM1, iy1 = iy0 - 1, iy2 = iy1 - 1;
        let ix: number, iy: number;

        const invDXDY = 1 / (this.mXSpacing * this.mYSpacing);

        // corners
        FXY.set(0, 0, 0.25 * invDXDY * (
            9 * F.get(0, 0)
            - 12 * F.get(1, 0)
            + 3 * F.get(2, 0)
            - 12 * F.get(0, 1)
            + 16 * F.get(1, 1)
            - 4 * F.get(2, 1)
            + 3 * F.get(0, 2)
            - 4 * F.get(1, 2)
            + F.get(2, 2)));

        FXY.set(xBoundM1, 0, 0.25 * invDXDY * (
            9 * F.get(ix0, 0)
            - 12 * F.get(ix1, 0)
            + 3 * F.get(ix2, 0)
            - 12 * F.get(ix0, 1)
            + 16 * F.get(ix1, 1)
            - 4 * F.get(ix2, 1)
            + 3 * F.get(ix0, 2)
            - 4 * F.get(ix1, 2)
            + F.get(ix2, 2)));

        FXY.set(0, yBoundM1, 0.25 * invDXDY * (
            9 * F.get(0, iy0)
            - 12 * F.get(1, iy0)
            + 3 * F.get(2, iy0)
            - 12 * F.get(0, iy1)
            + 16 * F.get(1, iy1)
            - 4 * F.get(2, iy1)
            + 3 * F.get(0, iy2)
            - 4 * F.get(1, iy2)
            + F.get(2, iy2)));

        FXY.set(xBoundM1, yBoundM1, 0.25 * invDXDY * (
            9 * F.get(ix0, iy0)
            - 12 * F.get(ix1, iy0)
            + 3 * F.get(ix2, iy0)
            - 12 * F.get(ix0, iy1)
            + 16 * F.get(ix1, iy1)
            - 4 * F.get(ix2, iy1)
            + 3 * F.get(ix0, iy2)
            - 4 * F.get(ix1, iy2)
            + F.get(ix2, iy2)));

        // x-edges
        for (ix = 1; ix < xBoundM1; ++ix) {
            FXY.set(ix, 0, 0.25 * invDXDY * (
                3 * (F.get(ix - 1, 0) - F.get(ix + 1, 0))
                - 4 * (F.get(ix - 1, 1) - F.get(ix + 1, 1))
                + (F.get(ix - 1, 2) - F.get(ix + 1, 2))));

            FXY.set(ix, yBoundM1, 0.25 * invDXDY * (
                3 * (F.get(ix - 1, iy0) - F.get(ix + 1, iy0))
                - 4 * (F.get(ix - 1, iy1) - F.get(ix + 1, iy1))
                + (F.get(ix - 1, iy2) - F.get(ix + 1, iy2))));
        }

        // y-edges
        for (iy = 1; iy < yBoundM1; ++iy) {
            FXY.set(0, iy, 0.25 * invDXDY * (
                3 * (F.get(0, iy - 1) - F.get(0, iy + 1))
                - 4 * (F.get(1, iy - 1) - F.get(1, iy + 1))
                + (F.get(2, iy - 1) - F.get(2, iy + 1))));

            FXY.set(xBoundM1, iy, 0.25 * invDXDY * (
                3 * (F.get(ix0, iy - 1) - F.get(ix0, iy + 1))
                - 4 * (F.get(ix1, iy - 1) - F.get(ix1, iy + 1))
                + (F.get(ix2, iy - 1) - F.get(ix2, iy + 1))));
        }

        // interior
        for (iy = 1; iy < yBoundM1; ++iy) {
            for (ix = 1; ix < xBoundM1; ++ix) {
                FXY.set(ix, iy, 0.25 * invDXDY * (F.get(ix - 1, iy - 1) -
                    F.get(ix + 1, iy - 1) - F.get(ix - 1, iy + 1) + F.get(ix + 1, iy + 1)));
            }
        }
    }

    private getPolynomials(F: Array2<number>, FX: Array2<number>,
        FY: Array2<number>, FXY: Array2<number>): void {
        const xBoundM1 = this.mXBound - 1;
        const yBoundM1 = this.mYBound - 1;
        for (let iy = 0; iy < yBoundM1; ++iy) {
            for (let ix = 0; ix < xBoundM1; ++ix) {
                // Note the 'transposing' of the 2x2 blocks (to match
                // notation used in the polynomial definition): G[a][b] is
                // the value at (x + a*dx, y + b*dy).
                const G = [
                    [F.get(ix, iy), F.get(ix, iy + 1)],
                    [F.get(ix + 1, iy), F.get(ix + 1, iy + 1)]
                ];

                const GX = [
                    [FX.get(ix, iy), FX.get(ix, iy + 1)],
                    [FX.get(ix + 1, iy), FX.get(ix + 1, iy + 1)]
                ];

                const GY = [
                    [FY.get(ix, iy), FY.get(ix, iy + 1)],
                    [FY.get(ix + 1, iy), FY.get(ix + 1, iy + 1)]
                ];

                const GXY = [
                    [FXY.get(ix, iy), FXY.get(ix, iy + 1)],
                    [FXY.get(ix + 1, iy), FXY.get(ix + 1, iy + 1)]
                ];

                this.construct(this.mPoly.get(ix, iy), G, GX, GY, GXY);
            }
        }
    }

    // Estimate the derivative from four consecutive slopes.
    private computeDerivative(slope0: number, slope1: number,
        slope2: number, slope3: number): number {
        if (slope1 !== slope2) {
            if (slope0 !== slope1) {
                if (slope2 !== slope3) {
                    const ad0 = Math.abs(slope3 - slope2);
                    const ad1 = Math.abs(slope0 - slope1);
                    return (ad0 * slope1 + ad1 * slope2) / (ad0 + ad1);
                }
                else {
                    return slope2;
                }
            }
            else {
                if (slope2 !== slope3) {
                    return slope1;
                }
                else {
                    return 0.5 * (slope1 + slope2);
                }
            }
        }
        else {
            return slope1;
        }
    }

    private construct(poly: Polynomial, F: number[][], FX: number[][],
        FY: number[][], FXY: number[][]): void {
        const dx = this.mXSpacing;
        const dy = this.mYSpacing;
        const invDX = 1 / dx, invDX2 = invDX * invDX;
        const invDY = 1 / dy, invDY2 = invDY * invDY;
        let b0: number, b1: number, b2: number, b3: number;

        poly.setA(0, 0, F[0][0]);
        poly.setA(1, 0, FX[0][0]);
        poly.setA(0, 1, FY[0][0]);
        poly.setA(1, 1, FXY[0][0]);

        b0 = (F[1][0] - poly.evaluate(0, 0, dx, 0)) * invDX2;
        b1 = (FX[1][0] - poly.evaluate(1, 0, dx, 0)) * invDX;
        poly.setA(2, 0, 3 * b0 - b1);
        poly.setA(3, 0, (-2 * b0 + b1) * invDX);

        b0 = (F[0][1] - poly.evaluate(0, 0, 0, dy)) * invDY2;
        b1 = (FY[0][1] - poly.evaluate(0, 1, 0, dy)) * invDY;
        poly.setA(0, 2, 3 * b0 - b1);
        poly.setA(0, 3, (-2 * b0 + b1) * invDY);

        b0 = (FY[1][0] - poly.evaluate(0, 1, dx, 0)) * invDX2;
        b1 = (FXY[1][0] - poly.evaluate(1, 1, dx, 0)) * invDX;
        poly.setA(2, 1, 3 * b0 - b1);
        poly.setA(3, 1, (-2 * b0 + b1) * invDX);

        b0 = (FX[0][1] - poly.evaluate(1, 0, 0, dy)) * invDY2;
        b1 = (FXY[0][1] - poly.evaluate(1, 1, 0, dy)) * invDY;
        poly.setA(1, 2, 3 * b0 - b1);
        poly.setA(1, 3, (-2 * b0 + b1) * invDY);

        b0 = (F[1][1] - poly.evaluate(0, 0, dx, dy)) * invDX2 * invDY2;
        b1 = (FX[1][1] - poly.evaluate(1, 0, dx, dy)) * invDX * invDY2;
        b2 = (FY[1][1] - poly.evaluate(0, 1, dx, dy)) * invDX2 * invDY;
        b3 = (FXY[1][1] - poly.evaluate(1, 1, dx, dy)) * invDX * invDY;
        poly.setA(2, 2, 9 * b0 - 3 * b1 - 3 * b2 + b3);
        poly.setA(3, 2, (-6 * b0 + 3 * b1 + 2 * b2 - b3) * invDX);
        poly.setA(2, 3, (-6 * b0 + 2 * b1 + 3 * b2 - b3) * invDY);
        poly.setA(3, 3, (4 * b0 - 2 * b1 - 2 * b2 + b3) * invDX * invDY);
    }

    // Support for evaluation.
    private xLookup(x: number): { index: number, delta: number } {
        for (let xIndex = 0, xIndexP1 = 1; xIndexP1 < this.mXBound; ++xIndex, ++xIndexP1) {
            if (x < this.mXMin + this.mXSpacing * xIndexP1) {
                return { index: xIndex, delta: x - (this.mXMin + this.mXSpacing * xIndex) };
            }
        }

        const xIndex = this.mXBound - 2;
        return { index: xIndex, delta: x - (this.mXMin + this.mXSpacing * xIndex) };
    }

    private yLookup(y: number): { index: number, delta: number } {
        for (let yIndex = 0, yIndexP1 = 1; yIndexP1 < this.mYBound; ++yIndex, ++yIndexP1) {
            if (y < this.mYMin + this.mYSpacing * yIndexP1) {
                return { index: yIndex, delta: y - (this.mYMin + this.mYSpacing * yIndex) };
            }
        }

        const yIndex = this.mYBound - 2;
        return { index: yIndex, delta: y - (this.mYMin + this.mYSpacing * yIndex) };
    }
}
