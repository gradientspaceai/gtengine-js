// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpBicubic2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The interpolator is for uniformly spaced (x,y)-values. The input samples F
// must be stored in row-major order to represent f(x,y); that is,
// F[c + xBound*r] corresponds to f(x,y), where c is the index corresponding
// to x and r is the index corresponding to y. Exact interpolation is achieved
// by setting catmullRom to 'true', giving you the Catmull-Rom blending
// matrix. If a smooth interpolation is desired, set catmullRom to 'false' to
// obtain B-spline blending.
//
// Port notes: the two upstream 'operator()' overloads become TypeScript
// overloads of evaluate(...) distinguished by arity: evaluate(x, y) for the
// function value and evaluate(xOrder, yOrder, x, y) for the function or its
// derivatives. The upstream value-only overload is the xOrder = yOrder = 0
// case of the second overload evaluated with the identical floating-point
// operations (the derivative multipliers are exactly 1), so the port
// implements the general case once and forwards to it.
//
// The samples array is aliased, not copied, matching the upstream
// 'Real const* F' member.

import { logAssert } from './Logger';

export class IntpBicubic2 {
    private mXBound: number;
    private mYBound: number;
    private mQuantity: number;
    private mXMin: number;
    private mXMax: number;
    private mXSpacing: number;
    private mInvXSpacing: number;
    private mYMin: number;
    private mYMax: number;
    private mYSpacing: number;
    private mInvYSpacing: number;
    private mF: readonly number[];
    private mBlend: number[][];

    // Construction.
    constructor(xBound: number, yBound: number, xMin: number, xSpacing: number,
        yMin: number, ySpacing: number, F: readonly number[], catmullRom: boolean) {
        this.mXBound = xBound;
        this.mYBound = yBound;
        this.mQuantity = xBound * yBound;
        this.mXMin = xMin;
        this.mXSpacing = xSpacing;
        this.mYMin = yMin;
        this.mYSpacing = ySpacing;
        this.mF = F;

        // At least a 3x3 block of data points is needed to construct the
        // estimates of the boundary derivatives.
        logAssert(this.mXBound >= 3 && this.mYBound >= 3 && F.length >= this.mQuantity,
            'Invalid input.');
        logAssert(this.mXSpacing > 0 && this.mYSpacing > 0, 'Invalid input.');

        this.mXMax = this.mXMin + this.mXSpacing * (this.mXBound - 1);
        this.mInvXSpacing = 1 / this.mXSpacing;
        this.mYMax = this.mYMin + this.mYSpacing * (this.mYBound - 1);
        this.mInvYSpacing = 1 / this.mYSpacing;

        if (catmullRom) {
            this.mBlend = [
                [0, -0.5, 1, -0.5],
                [1, 0, -2.5, 1.5],
                [0, 0.5, 2, -1.5],
                [0, 0, -0.5, 0.5]
            ];
        } else {
            this.mBlend = [
                [1 / 6, -3 / 6, 3 / 6, -1 / 6],
                [4 / 6, 0 / 6, -6 / 6, 3 / 6],
                [1 / 6, 3 / 6, 3 / 6, -3 / 6],
                [0 / 6, 0 / 6, 0 / 6, 1 / 6]
            ];
        }
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
    // inputs to xmin <= x <= xmax and ymin <= y <= ymax. The 2-argument form
    // is for function evaluation. The 4-argument form is for function or
    // derivative evaluations, where xOrder is the order of the x-derivative
    // and yOrder is the order of the y-derivative. Both orders are zero to
    // get the function value itself. Orders larger than 3 produce 0.
    evaluate(x: number, y: number): number;
    evaluate(xOrder: number, yOrder: number, x: number, y: number): number;
    evaluate(a0: number, a1: number, a2?: number, a3?: number): number {
        let xOrder: number, yOrder: number, x: number, y: number;
        if (a2 === undefined || a3 === undefined) {
            xOrder = 0;
            yOrder = 0;
            x = a0;
            y = a1;
        } else {
            xOrder = a0;
            yOrder = a1;
            x = a2;
            y = a3;
        }

        // Compute x-index and clamp to image.
        const xIndex = (x - this.mXMin) * this.mInvXSpacing;
        let ix = Math.trunc(xIndex);
        if (ix < 0) {
            ix = 0;
        } else if (ix >= this.mXBound) {
            ix = this.mXBound - 1;
        }

        // Compute y-index and clamp to image.
        const yIndex = (y - this.mYMin) * this.mInvYSpacing;
        let iy = Math.trunc(yIndex);
        if (iy < 0) {
            iy = 0;
        } else if (iy >= this.mYBound) {
            iy = this.mYBound - 1;
        }

        const U: number[] = [0, 0, 0, 0];
        let dx: number, xMult: number;
        switch (xOrder) {
            case 0:
                dx = xIndex - ix;
                U[0] = 1;
                U[1] = dx;
                U[2] = dx * U[1];
                U[3] = dx * U[2];
                xMult = 1;
                break;
            case 1:
                dx = xIndex - ix;
                U[0] = 0;
                U[1] = 1;
                U[2] = 2 * dx;
                U[3] = 3 * dx * dx;
                xMult = this.mInvXSpacing;
                break;
            case 2:
                dx = xIndex - ix;
                U[0] = 0;
                U[1] = 0;
                U[2] = 2;
                U[3] = 6 * dx;
                xMult = this.mInvXSpacing * this.mInvXSpacing;
                break;
            case 3:
                U[0] = 0;
                U[1] = 0;
                U[2] = 0;
                U[3] = 6;
                xMult = this.mInvXSpacing * this.mInvXSpacing * this.mInvXSpacing;
                break;
            default:
                return 0;
        }

        const V: number[] = [0, 0, 0, 0];
        let dy: number, yMult: number;
        switch (yOrder) {
            case 0:
                dy = yIndex - iy;
                V[0] = 1;
                V[1] = dy;
                V[2] = dy * V[1];
                V[3] = dy * V[2];
                yMult = 1;
                break;
            case 1:
                dy = yIndex - iy;
                V[0] = 0;
                V[1] = 1;
                V[2] = 2 * dy;
                V[3] = 3 * dy * dy;
                yMult = this.mInvYSpacing;
                break;
            case 2:
                dy = yIndex - iy;
                V[0] = 0;
                V[1] = 0;
                V[2] = 2;
                V[3] = 6 * dy;
                yMult = this.mInvYSpacing * this.mInvYSpacing;
                break;
            case 3:
                V[0] = 0;
                V[1] = 0;
                V[2] = 0;
                V[3] = 6;
                yMult = this.mInvYSpacing * this.mInvYSpacing * this.mInvYSpacing;
                break;
            default:
                return 0;
        }

        // Compute P = M*U and Q = M*V.
        const P: number[] = [0, 0, 0, 0];
        const Q: number[] = [0, 0, 0, 0];
        for (let row = 0; row < 4; ++row) {
            P[row] = 0;
            Q[row] = 0;
            for (let col = 0; col < 4; ++col) {
                P[row] += this.mBlend[row][col] * U[col];
                Q[row] += this.mBlend[row][col] * V[col];
            }
        }

        // Compute (M*U)^t D (M*V) where D is the 4x4 subimage containing
        // (x,y).
        --ix;
        --iy;
        let result = 0;
        for (let row = 0; row < 4; ++row) {
            let yClamp = iy + row;
            if (yClamp < 0) {
                yClamp = 0;
            } else if (yClamp > this.mYBound - 1) {
                yClamp = this.mYBound - 1;
            }

            for (let col = 0; col < 4; ++col) {
                let xClamp = ix + col;
                if (xClamp < 0) {
                    xClamp = 0;
                } else if (xClamp > this.mXBound - 1) {
                    xClamp = this.mXBound - 1;
                }

                result += P[col] * Q[row] * this.mF[xClamp + this.mXBound * yClamp];
            }
        }
        result *= xMult * yMult;

        return result;
    }
}
