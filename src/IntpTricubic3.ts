// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpTricubic3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The interpolator is for uniformly spaced (x,y,z)-values. The input samples
// must be stored in lexicographical order to represent f(x,y,z); that is,
// F[c + xBound*(r + yBound*s)] corresponds to f(x,y,z), where c is the index
// corresponding to x, r is the index corresponding to y, and s is the index
// corresponding to z. Exact interpolation is achieved by setting catmullRom
// to 'true', giving you the Catmull-Rom blending matrix. If a smooth
// interpolation is desired, set catmullRom to 'false' to obtain B-spline
// blending.
//
// Port notes:
// - The two upstream 'operator()' overloads become the overloaded method
//   'evaluate': 'evaluate(x, y, z)' for the function value and
//   'evaluate(xOrder, yOrder, zOrder, x, y, z)' for a mixed partial
//   derivative (all orders zero gives the function value).
// - The samples F are referenced, not copied, matching the upstream
//   'Real const* mF'. Mutating the caller's array changes later evaluations.
// - 'static_cast<int32_t>' truncates toward zero, so it is Math.trunc here
//   (not Math.floor: the two differ for the negative indices produced by
//   inputs below the domain minimum).
// - Upstream validates only 'F != nullptr'. Because a too-short array is an
//   out-of-bounds read in C++ but a silent 'undefined' in TypeScript, the
//   port also requires F.length >= xBound*yBound*zBound.
// - The B-spline blending entries are written as the upstream quotients
//   (4/6, -3/6, ...) rather than decimal literals so the rounding matches.
// - The upstream comment claims the evaluators clamp the inputs to
//   xMin <= x <= xMax (and likewise in y and z). They do not: only the index
//   of the containing cell is clamped, while the fractional cell coordinate
//   keeps its out-of-range value. An input outside the domain therefore
//   extrapolates the boundary cell's cubic polynomial. The port preserves
//   this behavior.

import { logAssert } from './Logger.js';

export class IntpTricubic3 {
    private readonly mXBound: number;
    private readonly mYBound: number;
    private readonly mZBound: number;
    private readonly mQuantity: number;
    private readonly mXMin: number;
    private readonly mXMax: number;
    private readonly mXSpacing: number;
    private readonly mInvXSpacing: number;
    private readonly mYMin: number;
    private readonly mYMax: number;
    private readonly mYSpacing: number;
    private readonly mInvYSpacing: number;
    private readonly mZMin: number;
    private readonly mZMax: number;
    private readonly mZSpacing: number;
    private readonly mInvZSpacing: number;
    private readonly mF: ArrayLike<number>;
    private readonly mBlend: number[][];

    // Construction. The samples F are stored by reference.
    constructor(xBound: number, yBound: number, zBound: number, xMin: number,
        xSpacing: number, yMin: number, ySpacing: number, zMin: number,
        zSpacing: number, F: ArrayLike<number>, catmullRom: boolean) {
        // At least a 4x4x4 block of data points are needed to construct the
        // tricubic interpolation.
        logAssert(xBound >= 4 && yBound >= 4 && zBound >= 4
            && F !== null && F !== undefined
            && F.length >= xBound * yBound * zBound
            && xSpacing > 0 && ySpacing > 0 && zSpacing > 0,
            'Invalid input.');

        this.mXBound = xBound;
        this.mYBound = yBound;
        this.mZBound = zBound;
        this.mQuantity = xBound * yBound * zBound;
        this.mXMin = xMin;
        this.mXSpacing = xSpacing;
        this.mYMin = yMin;
        this.mYSpacing = ySpacing;
        this.mZMin = zMin;
        this.mZSpacing = zSpacing;
        this.mF = F;

        this.mXMax = this.mXMin + this.mXSpacing * (this.mXBound - 1);
        this.mInvXSpacing = 1 / this.mXSpacing;
        this.mYMax = this.mYMin + this.mYSpacing * (this.mYBound - 1);
        this.mInvYSpacing = 1 / this.mYSpacing;
        this.mZMax = this.mZMin + this.mZSpacing * (this.mZBound - 1);
        this.mInvZSpacing = 1 / this.mZSpacing;

        if (catmullRom) {
            this.mBlend = [
                [0, -0.5, 1, -0.5],
                [1, 0, -2.5, 1.5],
                [0, 0.5, 2, -1.5],
                [0, 0, -0.5, 0.5]
            ];
        }
        else {
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

    getZBound(): number {
        return this.mZBound;
    }

    getQuantity(): number {
        return this.mQuantity;
    }

    getF(): ArrayLike<number> {
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

    getZMin(): number {
        return this.mZMin;
    }

    getZMax(): number {
        return this.mZMax;
    }

    getZSpacing(): number {
        return this.mZSpacing;
    }

    // Evaluate the function and its derivatives. evaluate(x, y, z) is for
    // function evaluation. evaluate(xOrder, yOrder, zOrder, x, y, z) is for
    // function or derivative evaluations; the xOrder argument is the order of
    // the x-derivative, the yOrder argument is the order of the y-derivative
    // and the zOrder argument is the order of the z-derivative. All orders
    // are zero to get the function value itself. An order outside [0,3]
    // produces 0, the derivative of the cubic polynomial on the interior of
    // a cell. See the class comment for the treatment of inputs outside
    // [xMin,xMax] x [yMin,yMax] x [zMin,zMax].
    evaluate(x: number, y: number, z: number): number;
    evaluate(xOrder: number, yOrder: number, zOrder: number, x: number,
        y: number, z: number): number;
    evaluate(arg0: number, arg1: number, arg2: number, arg3?: number,
        arg4?: number, arg5?: number): number {
        if (arg3 === undefined || arg4 === undefined || arg5 === undefined) {
            // Compute the indices and clamp them to the image.
            const xIndex = (arg0 - this.mXMin) * this.mInvXSpacing;
            const ix = clampIndex(Math.trunc(xIndex), this.mXBound);
            const yIndex = (arg1 - this.mYMin) * this.mInvYSpacing;
            const iy = clampIndex(Math.trunc(yIndex), this.mYBound);
            const zIndex = (arg2 - this.mZMin) * this.mInvZSpacing;
            const iz = clampIndex(Math.trunc(zIndex), this.mZBound);

            const U = powers(xIndex - ix);
            const V = powers(yIndex - iy);
            const W = powers(zIndex - iz);

            return this.tensorProduct(U, V, W, ix, iy, iz);
        }

        const xOrder = arg0;
        const yOrder = arg1;
        const zOrder = arg2;

        // Compute the indices and clamp them to the image.
        const xIndex = (arg3 - this.mXMin) * this.mInvXSpacing;
        const ix = clampIndex(Math.trunc(xIndex), this.mXBound);
        const yIndex = (arg4 - this.mYMin) * this.mInvYSpacing;
        const iy = clampIndex(Math.trunc(yIndex), this.mYBound);
        const zIndex = (arg5 - this.mZMin) * this.mInvZSpacing;
        const iz = clampIndex(Math.trunc(zIndex), this.mZBound);

        const U = derivativePowers(xOrder, xIndex - ix);
        if (U === null) {
            return 0;
        }
        const xMult = multiplier(this.mInvXSpacing, xOrder);

        const V = derivativePowers(yOrder, yIndex - iy);
        if (V === null) {
            return 0;
        }
        const yMult = multiplier(this.mInvYSpacing, yOrder);

        const W = derivativePowers(zOrder, zIndex - iz);
        if (W === null) {
            return 0;
        }
        const zMult = multiplier(this.mInvZSpacing, zOrder);

        // Upstream multiplies the accumulated result by the single factor
        // xMult*yMult*zMult; keep that association for bit-identical rounding.
        return this.tensorProduct(U, V, W, ix, iy, iz) * (xMult * yMult * zMult);
    }

    // Compute P = M*U, Q = M*V and R = M*W, then the tensor product
    // (M*U)(M*V)(M*W)*D where D is the 4x4x4 subimage containing (x,y,z).
    private tensorProduct(U: number[], V: number[], W: number[],
        ix: number, iy: number, iz: number): number {
        const P = [0, 0, 0, 0], Q = [0, 0, 0, 0], R = [0, 0, 0, 0];
        for (let row = 0; row < 4; ++row) {
            const blendRow = this.mBlend[row];
            for (let col = 0; col < 4; ++col) {
                P[row] += blendRow[col] * U[col];
                Q[row] += blendRow[col] * V[col];
                R[row] += blendRow[col] * W[col];
            }
        }

        // The 4x4x4 subimage starts one sample before the containing cell.
        const ix0 = ix - 1;
        const iy0 = iy - 1;
        const iz0 = iz - 1;

        let result = 0;
        for (let slice = 0; slice < 4; ++slice) {
            const zClamp = clampIndex(iz0 + slice, this.mZBound);

            for (let row = 0; row < 4; ++row) {
                const yClamp = clampIndex(iy0 + row, this.mYBound);

                for (let col = 0; col < 4; ++col) {
                    const xClamp = clampIndex(ix0 + col, this.mXBound);

                    result += P[col] * Q[row] * R[slice] *
                        this.mF[xClamp + this.mXBound * (yClamp + this.mYBound * zClamp)];
                }
            }
        }

        return result;
    }
}

// Clamp a sample index to [0, bound-1]. Upstream writes this inline for each
// of x, y and z in both evaluators.
function clampIndex(index: number, bound: number): number {
    if (index < 0) {
        return 0;
    }
    if (index >= bound) {
        return bound - 1;
    }
    return index;
}

// The polynomial basis (1, d, d^2, d^3) evaluated at the fractional cell
// coordinate d.
function powers(d: number): number[] {
    const d2 = d * d;
    return [1, d, d2, d * d2];
}

// The 'order'-th derivative of (1, d, d^2, d^3) with respect to d. A null
// return corresponds to the upstream 'default: return (Real)0' branch, which
// applies to all orders outside [0,3].
function derivativePowers(order: number, d: number): number[] | null {
    switch (order) {
        case 0:
            return powers(d);
        case 1:
            return [0, 1, 2 * d, 3 * d * d];
        case 2:
            return [0, 0, 2, 6 * d];
        case 3:
            return [0, 0, 0, 6];
        default:
            return null;
    }
}

// The chain-rule factor (1/spacing)^order that converts a derivative with
// respect to the fractional cell coordinate into one with respect to the
// world coordinate. The products are written out rather than using Math.pow
// so the rounding matches the upstream expressions.
function multiplier(invSpacing: number, order: number): number {
    switch (order) {
        case 0:
            return 1;
        case 1:
            return invSpacing;
        case 2:
            return invSpacing * invSpacing;
        default:
            return invSpacing * invSpacing * invSpacing;
    }
}
