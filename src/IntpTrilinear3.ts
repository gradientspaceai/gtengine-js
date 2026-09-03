// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpTrilinear3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The interpolator is for uniformly spaced (x,y,z)-values. The input samples
// must be stored in lexicographical order to represent f(x,y,z); that is,
// F[c + xBound*(r + yBound*s)] corresponds to f(x,y,z), where c is the index
// corresponding to x, r is the index corresponding to y, and s is the index
// corresponding to z.
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
// - The upstream comment claims the evaluators clamp the inputs to
//   xMin <= x <= xMax (and likewise in y and z). They do not: only the index
//   of the containing cell is clamped, while the fractional cell coordinate
//   keeps its out-of-range value. Above the maximum this is harmless, because
//   both entries of that axis' stencil clamp to the last sample and the two
//   blend weights sum to one, so the boundary value is held. Below the
//   minimum the stencil still spans two distinct samples, so the boundary
//   cell's linear polynomial is extrapolated instead of clamped. The port
//   preserves this behavior.

import { logAssert } from './Logger.js';

export class IntpTrilinear3 {
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
        zSpacing: number, F: ArrayLike<number>) {
        // At least a 2x2x2 block of data points are needed to construct the
        // trilinear interpolation.
        logAssert(xBound >= 2 && yBound >= 2 && zBound >= 2
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

        this.mBlend = [
            [1, -1],
            [0, 1]
        ];
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
    // are zero to get the function value itself. An order larger than 1
    // produces 0, the derivative of the trilinear polynomial on the interior
    // of a cell. See the class comment for the treatment of inputs outside
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

            const U = [1, xIndex - ix];
            const V = [1, yIndex - iy];
            const W = [1, zIndex - iz];

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

        let U: number[], xMult: number;
        if (xOrder === 0) {
            U = [1, xIndex - ix];
            xMult = 1;
        }
        else if (xOrder === 1) {
            U = [0, 1];
            xMult = this.mInvXSpacing;
        }
        else {
            return 0;
        }

        let V: number[], yMult: number;
        if (yOrder === 0) {
            V = [1, yIndex - iy];
            yMult = 1;
        }
        else if (yOrder === 1) {
            V = [0, 1];
            yMult = this.mInvYSpacing;
        }
        else {
            return 0;
        }

        let W: number[], zMult: number;
        if (zOrder === 0) {
            W = [1, zIndex - iz];
            zMult = 1;
        }
        else if (zOrder === 1) {
            W = [0, 1];
            zMult = this.mInvZSpacing;
        }
        else {
            return 0;
        }

        // Upstream multiplies the accumulated result by the single factor
        // xMult*yMult*zMult; keep that association for bit-identical rounding.
        return this.tensorProduct(U, V, W, ix, iy, iz) * (xMult * yMult * zMult);
    }

    // Compute P = M*U, Q = M*V and R = M*W, then the tensor product
    // (M*U)(M*V)(M*W)*D where D is the 2x2x2 subimage containing (x,y,z).
    private tensorProduct(U: number[], V: number[], W: number[],
        ix: number, iy: number, iz: number): number {
        const P = [0, 0], Q = [0, 0], R = [0, 0];
        for (let row = 0; row < 2; ++row) {
            const blendRow = this.mBlend[row];
            for (let col = 0; col < 2; ++col) {
                P[row] += blendRow[col] * U[col];
                Q[row] += blendRow[col] * V[col];
                R[row] += blendRow[col] * W[col];
            }
        }

        let result = 0;
        for (let slice = 0; slice < 2; ++slice) {
            let zClamp = iz + slice;
            if (zClamp >= this.mZBound) {
                zClamp = this.mZBound - 1;
            }

            for (let row = 0; row < 2; ++row) {
                let yClamp = iy + row;
                if (yClamp >= this.mYBound) {
                    yClamp = this.mYBound - 1;
                }

                for (let col = 0; col < 2; ++col) {
                    let xClamp = ix + col;
                    if (xClamp >= this.mXBound) {
                        xClamp = this.mXBound - 1;
                    }

                    result += P[col] * Q[row] * R[slice] *
                        this.mF[xClamp + this.mXBound * (yClamp + this.mYBound * zClamp)];
                }
            }
        }

        return result;
    }
}

// Clamp a truncated sample index to [0, bound-1]. Upstream writes this inline
// for each of x, y and z in both evaluators.
function clampIndex(index: number, bound: number): number {
    if (index < 0) {
        return 0;
    }
    if (index >= bound) {
        return bound - 1;
    }
    return index;
}
