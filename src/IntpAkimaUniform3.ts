// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpAkimaUniform3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Akima interpolation is described in
// https://en.wikipedia.org/wiki/Akima_spline
// The interpolator is for uniformly spaced (x,y,z)-values. The input samples
// must be stored in lexicographical order to represent f(x,y,z); that is,
// F[c + xBound*(r + yBound*s)] corresponds to f(x,y,z), where c is the index
// corresponding to x, r is the index corresponding to y, and s is the index
// corresponding to z.
//
// Port notes:
//   - Upstream 'F[iz][iy][ix]' on an Array3<Real> becomes 'F.get(ix, iy, iz)',
//     following the Array3 port convention (index order matches the
//     constructor's bound order).
//   - The two upstream 'operator()' overloads become TypeScript overloads of
//     evaluate(...) distinguished by arity: evaluate(x, y, z) for the function
//     value and evaluate(xOrder, yOrder, zOrder, x, y, z) for the function or
//     its derivatives. The value-only form is exactly the all-orders-zero case
//     (the same floating-point operations), so the port implements the general
//     case once and forwards to it.
//   - The upstream private nested class Polynomial is a module-private class
//     here; its 4x4x4 coefficient tensor is stored in a flat Float64Array with
//     'A[ix][iy][iz]' at flat index ix + 4 * (iy + 4 * iz).
//   - Upstream passes a pointer into a slope row to ComputeDerivative; the port
//     passes the four slope values.
//   - The 2x2x2 corner blocks that upstream builds as Real[2][2][2] locals are
//     built as flat Float64Array(8) values with 'B[i][j][k]' at index
//     4*i + 2*j + k (i is the x-offset, j the y-offset, k the z-offset), which
//     preserves upstream's 'transposed' block indexing.
//   - The samples array is aliased by the mF member (getF returns it), as
//     upstream's 'Real const* F' is. The Array3 wrapper used during
//     construction is built from a copy of the first xBound*yBound*zBound
//     samples so that a longer input array is accepted, as in the
//     IntpAkimaUniform2 port.
//   - Upstream bug preserved: GetFXY, GetFXZ, GetFYZ and GetFXYZ reuse the
//     min-boundary one-sided difference coefficients at the max boundaries
//     with reflected sample indices and without negating the mask, so each
//     reflected direction flips the sign of the estimated mixed partial
//     derivative there. This is the 3D analogue of the IntpAkimaUniform2
//     issue (gtengine-js issue #58); as in that port the quirk is preserved
//     so the TypeScript results match upstream exactly.

import { Array3 } from './Array3';
import { logAssert } from './Logger';

// Flat indices into a 2x2x2 corner block; BIJK is upstream's B[i][j][k].
const B000 = 0;
const B001 = 1;
const B010 = 2;
const B011 = 3;
const B100 = 4;
const B101 = 5;
const B110 = 6;
const B111 = 7;

// P(x,y,z) = sum_{i=0}^3 sum_{j=0}^3 sum_{k=0}^3 a_{ijk} x^i y^j z^k. The
// tensor term A(ix,iy,iz) corresponds to the polynomial term
// x^{ix} y^{iy} z^{iz}.
class AkimaPolynomial3 {
    // The coefficients are zero-initialized, as upstream's constructor does.
    private mCoeff: Float64Array = new Float64Array(64);

    getA(ix: number, iy: number, iz: number): number {
        return this.mCoeff[ix + 4 * (iy + 4 * iz)];
    }

    setA(ix: number, iy: number, iz: number, value: number): void {
        this.mCoeff[ix + 4 * (iy + 4 * iz)] = value;
    }

    // Evaluate the polynomial or one of its partial derivatives. Any order
    // larger than 3 gives 0 because the polynomial has degree 3 in each
    // variable.
    evaluate(xOrder: number, yOrder: number, zOrder: number,
        x: number, y: number, z: number): number {
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

        const zPow: number[] = [0, 0, 0, 0];
        switch (zOrder) {
            case 0:
                zPow[0] = 1;
                zPow[1] = z;
                zPow[2] = z * z;
                zPow[3] = z * z * z;
                break;
            case 1:
                zPow[0] = 0;
                zPow[1] = 1;
                zPow[2] = 2 * z;
                zPow[3] = 3 * z * z;
                break;
            case 2:
                zPow[0] = 0;
                zPow[1] = 0;
                zPow[2] = 2;
                zPow[3] = 6 * z;
                break;
            case 3:
                zPow[0] = 0;
                zPow[1] = 0;
                zPow[2] = 0;
                zPow[3] = 6;
                break;
            default:
                return 0;
        }

        let p = 0;
        for (let iz = 0; iz <= 3; ++iz) {
            for (let iy = 0; iy <= 3; ++iy) {
                for (let ix = 0; ix <= 3; ++ix) {
                    p += this.mCoeff[ix + 4 * (iy + 4 * iz)] * xPow[ix] * yPow[iy] * zPow[iz];
                }
            }
        }

        return p;
    }
}

// Extract the 2x2x2 block of A whose lowest corner is (ix, iy, iz), stored so
// that upstream's B[i][j][k] is at flat index 4*i + 2*j + k.
function extractBlock(A: Array3<number>, ix: number, iy: number, iz: number): Float64Array {
    const B = new Float64Array(8);
    for (let i = 0; i < 2; ++i) {
        for (let j = 0; j < 2; ++j) {
            for (let k = 0; k < 2; ++k) {
                B[4 * i + 2 * j + k] = A.get(ix + i, iy + j, iz + k);
            }
        }
    }
    return B;
}

function makeArray3(bound0: number, bound1: number, bound2: number): Array3<number> {
    const a = new Array3<number>(bound0, bound1, bound2);
    a.fill(0);
    return a;
}

export class IntpAkimaUniform3 {
    private mXBound: number;
    private mYBound: number;
    private mZBound: number;
    private mQuantity: number;
    private mXMin: number;
    private mXMax: number;
    private mXSpacing: number;
    private mYMin: number;
    private mYMax: number;
    private mYSpacing: number;
    private mZMin: number;
    private mZMax: number;
    private mZSpacing: number;
    private mF: readonly number[];
    private mPoly: Array3<AkimaPolynomial3>;

    // Construction. The samples F are in lexicographical order with x varying
    // fastest and z varying slowest.
    constructor(xBound: number, yBound: number, zBound: number,
        xMin: number, xSpacing: number, yMin: number, ySpacing: number,
        zMin: number, zSpacing: number, F: readonly number[]) {
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

        // At least a 3x3x3 block of data points is needed to construct the
        // estimates of the boundary derivatives.
        logAssert(this.mXBound >= 3 && this.mYBound >= 3 && this.mZBound >= 3 &&
            F.length >= this.mQuantity, 'Invalid input.');
        logAssert(this.mXSpacing > 0 && this.mYSpacing > 0 && this.mZSpacing > 0,
            'Invalid input.');

        this.mXMax = this.mXMin + this.mXSpacing * (this.mXBound - 1);
        this.mYMax = this.mYMin + this.mYSpacing * (this.mYBound - 1);
        this.mZMax = this.mZMin + this.mZSpacing * (this.mZBound - 1);

        this.mPoly = new Array3<AkimaPolynomial3>(xBound - 1, yBound - 1, zBound - 1);
        for (let i = 0; i < this.mPoly.data().length; ++i) {
            this.mPoly.data()[i] = new AkimaPolynomial3();
        }

        // Create a 3D wrapper for the 1D samples.
        const Fmap = new Array3<number>(this.mXBound, this.mYBound, this.mZBound,
            F.slice(0, this.mQuantity));

        // Construct first-order derivatives.
        const FX = makeArray3(this.mXBound, this.mYBound, this.mZBound);
        const FY = makeArray3(this.mXBound, this.mYBound, this.mZBound);
        const FZ = makeArray3(this.mXBound, this.mYBound, this.mZBound);
        this.getFX(Fmap, FX);
        this.getFY(Fmap, FY);
        this.getFZ(Fmap, FZ);

        // Construct second-order derivatives.
        const FXY = makeArray3(this.mXBound, this.mYBound, this.mZBound);
        const FXZ = makeArray3(this.mXBound, this.mYBound, this.mZBound);
        const FYZ = makeArray3(this.mXBound, this.mYBound, this.mZBound);
        this.getFXY(Fmap, FXY);
        this.getFXZ(Fmap, FXZ);
        this.getFYZ(Fmap, FYZ);

        // Construct third-order derivatives.
        const FXYZ = makeArray3(this.mXBound, this.mYBound, this.mZBound);
        this.getFXYZ(Fmap, FXYZ);

        // Construct polynomials.
        this.getPolynomials(Fmap, FX, FY, FZ, FXY, FXZ, FYZ, FXYZ);
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

    getZMin(): number {
        return this.mZMin;
    }

    getZMax(): number {
        return this.mZMax;
    }

    getZSpacing(): number {
        return this.mZSpacing;
    }

    // Evaluate the function and its derivatives. The inputs are clamped to
    // xmin <= x <= xmax, ymin <= y <= ymax and zmin <= z <= zmax. The
    // 3-argument form is for function evaluation. The 6-argument form is for
    // function or derivative evaluations, where xOrder, yOrder and zOrder are
    // the orders of the x-, y- and z-derivatives. All orders are zero to get
    // the function value itself.
    evaluate(x: number, y: number, z: number): number;
    evaluate(xOrder: number, yOrder: number, zOrder: number,
        x: number, y: number, z: number): number;
    evaluate(a0: number, a1: number, a2: number,
        a3?: number, a4?: number, a5?: number): number {
        let xOrder: number, yOrder: number, zOrder: number;
        let x: number, y: number, z: number;
        if (a3 === undefined || a4 === undefined || a5 === undefined) {
            xOrder = 0;
            yOrder = 0;
            zOrder = 0;
            x = a0;
            y = a1;
            z = a2;
        } else {
            xOrder = a0;
            yOrder = a1;
            zOrder = a2;
            x = a3;
            y = a4;
            z = a5;
        }

        x = Math.min(Math.max(x, this.mXMin), this.mXMax);
        y = Math.min(Math.max(y, this.mYMin), this.mYMax);
        z = Math.min(Math.max(z, this.mZMin), this.mZMax);
        const xl = this.xLookup(x);
        const yl = this.yLookup(y);
        const zl = this.zLookup(z);
        return this.mPoly.get(xl.index, yl.index, zl.index).evaluate(
            xOrder, yOrder, zOrder, xl.d, yl.d, zl.d);
    }

    // Support for construction.
    private getFX(F: Array3<number>, FX: Array3<number>): void {
        const slope = makeArray3(this.mXBound + 3, this.mYBound, this.mZBound);
        const invDX = 1 / this.mXSpacing;
        for (let iz = 0; iz < this.mZBound; ++iz) {
            for (let iy = 0; iy < this.mYBound; ++iy) {
                for (let ix = 0; ix < this.mXBound - 1; ++ix) {
                    slope.set(ix + 2, iy, iz,
                        (F.get(ix + 1, iy, iz) - F.get(ix, iy, iz)) * invDX);
                }

                slope.set(1, iy, iz,
                    2 * slope.get(2, iy, iz) - slope.get(3, iy, iz));
                slope.set(0, iy, iz,
                    2 * slope.get(1, iy, iz) - slope.get(2, iy, iz));
                slope.set(this.mXBound + 1, iy, iz,
                    2 * slope.get(this.mXBound, iy, iz) - slope.get(this.mXBound - 1, iy, iz));
                slope.set(this.mXBound + 2, iy, iz,
                    2 * slope.get(this.mXBound + 1, iy, iz) - slope.get(this.mXBound, iy, iz));
            }
        }

        for (let iz = 0; iz < this.mZBound; ++iz) {
            for (let iy = 0; iy < this.mYBound; ++iy) {
                for (let ix = 0; ix < this.mXBound; ++ix) {
                    FX.set(ix, iy, iz, IntpAkimaUniform3.computeDerivative(
                        slope.get(ix, iy, iz), slope.get(ix + 1, iy, iz),
                        slope.get(ix + 2, iy, iz), slope.get(ix + 3, iy, iz)));
                }
            }
        }
    }

    private getFY(F: Array3<number>, FY: Array3<number>): void {
        // Upstream stores the y-slopes as slope[iz][ix][iy], so bound0 is
        // yBound + 3, bound1 is xBound and bound2 is zBound.
        const slope = makeArray3(this.mYBound + 3, this.mXBound, this.mZBound);
        const invDY = 1 / this.mYSpacing;
        for (let iz = 0; iz < this.mZBound; ++iz) {
            for (let ix = 0; ix < this.mXBound; ++ix) {
                for (let iy = 0; iy < this.mYBound - 1; ++iy) {
                    slope.set(iy + 2, ix, iz,
                        (F.get(ix, iy + 1, iz) - F.get(ix, iy, iz)) * invDY);
                }

                slope.set(1, ix, iz,
                    2 * slope.get(2, ix, iz) - slope.get(3, ix, iz));
                slope.set(0, ix, iz,
                    2 * slope.get(1, ix, iz) - slope.get(2, ix, iz));
                slope.set(this.mYBound + 1, ix, iz,
                    2 * slope.get(this.mYBound, ix, iz) - slope.get(this.mYBound - 1, ix, iz));
                slope.set(this.mYBound + 2, ix, iz,
                    2 * slope.get(this.mYBound + 1, ix, iz) - slope.get(this.mYBound, ix, iz));
            }
        }

        for (let iz = 0; iz < this.mZBound; ++iz) {
            for (let ix = 0; ix < this.mXBound; ++ix) {
                for (let iy = 0; iy < this.mYBound; ++iy) {
                    FY.set(ix, iy, iz, IntpAkimaUniform3.computeDerivative(
                        slope.get(iy, ix, iz), slope.get(iy + 1, ix, iz),
                        slope.get(iy + 2, ix, iz), slope.get(iy + 3, ix, iz)));
                }
            }
        }
    }

    private getFZ(F: Array3<number>, FZ: Array3<number>): void {
        // Upstream stores the z-slopes as slope[iy][ix][iz], so bound0 is
        // zBound + 3, bound1 is xBound and bound2 is yBound.
        const slope = makeArray3(this.mZBound + 3, this.mXBound, this.mYBound);
        const invDZ = 1 / this.mZSpacing;
        for (let iy = 0; iy < this.mYBound; ++iy) {
            for (let ix = 0; ix < this.mXBound; ++ix) {
                for (let iz = 0; iz < this.mZBound - 1; ++iz) {
                    slope.set(iz + 2, ix, iy,
                        (F.get(ix, iy, iz + 1) - F.get(ix, iy, iz)) * invDZ);
                }

                slope.set(1, ix, iy,
                    2 * slope.get(2, ix, iy) - slope.get(3, ix, iy));
                slope.set(0, ix, iy,
                    2 * slope.get(1, ix, iy) - slope.get(2, ix, iy));
                slope.set(this.mZBound + 1, ix, iy,
                    2 * slope.get(this.mZBound, ix, iy) - slope.get(this.mZBound - 1, ix, iy));
                slope.set(this.mZBound + 2, ix, iy,
                    2 * slope.get(this.mZBound + 1, ix, iy) - slope.get(this.mZBound, ix, iy));
            }
        }

        for (let iy = 0; iy < this.mYBound; ++iy) {
            for (let ix = 0; ix < this.mXBound; ++ix) {
                for (let iz = 0; iz < this.mZBound; ++iz) {
                    FZ.set(ix, iy, iz, IntpAkimaUniform3.computeDerivative(
                        slope.get(iz, ix, iy), slope.get(iz + 1, ix, iy),
                        slope.get(iz + 2, ix, iy), slope.get(iz + 3, ix, iy)));
                }
            }
        }
    }

    private getFXY(F: Array3<number>, FXY: Array3<number>): void {
        const xBoundM1 = this.mXBound - 1;
        const yBoundM1 = this.mYBound - 1;
        const ix0 = xBoundM1, ix1 = ix0 - 1, ix2 = ix1 - 1;
        const iy0 = yBoundM1, iy1 = iy0 - 1, iy2 = iy1 - 1;

        const invDXDY = 1 / (this.mXSpacing * this.mYSpacing);
        for (let iz = 0; iz < this.mZBound; ++iz) {
            // corners of z-slice
            FXY.set(0, 0, iz, 0.25 * invDXDY * (
                9 * F.get(0, 0, iz)
                - 12 * F.get(1, 0, iz)
                + 3 * F.get(2, 0, iz)
                - 12 * F.get(0, 1, iz)
                + 16 * F.get(1, 1, iz)
                - 4 * F.get(2, 1, iz)
                + 3 * F.get(0, 2, iz)
                - 4 * F.get(1, 2, iz)
                + F.get(2, 2, iz)));

            FXY.set(xBoundM1, 0, iz, 0.25 * invDXDY * (
                9 * F.get(ix0, 0, iz)
                - 12 * F.get(ix1, 0, iz)
                + 3 * F.get(ix2, 0, iz)
                - 12 * F.get(ix0, 1, iz)
                + 16 * F.get(ix1, 1, iz)
                - 4 * F.get(ix2, 1, iz)
                + 3 * F.get(ix0, 2, iz)
                - 4 * F.get(ix1, 2, iz)
                + F.get(ix2, 2, iz)));

            FXY.set(0, yBoundM1, iz, 0.25 * invDXDY * (
                9 * F.get(0, iy0, iz)
                - 12 * F.get(1, iy0, iz)
                + 3 * F.get(2, iy0, iz)
                - 12 * F.get(0, iy1, iz)
                + 16 * F.get(1, iy1, iz)
                - 4 * F.get(2, iy1, iz)
                + 3 * F.get(0, iy2, iz)
                - 4 * F.get(1, iy2, iz)
                + F.get(2, iy2, iz)));

            FXY.set(xBoundM1, yBoundM1, iz, 0.25 * invDXDY * (
                9 * F.get(ix0, iy0, iz)
                - 12 * F.get(ix1, iy0, iz)
                + 3 * F.get(ix2, iy0, iz)
                - 12 * F.get(ix0, iy1, iz)
                + 16 * F.get(ix1, iy1, iz)
                - 4 * F.get(ix2, iy1, iz)
                + 3 * F.get(ix0, iy2, iz)
                - 4 * F.get(ix1, iy2, iz)
                + F.get(ix2, iy2, iz)));

            // x-edges of z-slice
            for (let ix = 1; ix < xBoundM1; ++ix) {
                FXY.set(ix, 0, iz, 0.25 * invDXDY * (
                    3 * (F.get(ix - 1, 0, iz) - F.get(ix + 1, 0, iz)) -
                    4 * (F.get(ix - 1, 1, iz) - F.get(ix + 1, 1, iz)) +
                    (F.get(ix - 1, 2, iz) - F.get(ix + 1, 2, iz))));

                FXY.set(ix, yBoundM1, iz, 0.25 * invDXDY * (
                    3 * (F.get(ix - 1, iy0, iz) - F.get(ix + 1, iy0, iz))
                    - 4 * (F.get(ix - 1, iy1, iz) - F.get(ix + 1, iy1, iz)) +
                    (F.get(ix - 1, iy2, iz) - F.get(ix + 1, iy2, iz))));
            }

            // y-edges of z-slice
            for (let iy = 1; iy < yBoundM1; ++iy) {
                FXY.set(0, iy, iz, 0.25 * invDXDY * (
                    3 * (F.get(0, iy - 1, iz) - F.get(0, iy + 1, iz)) -
                    4 * (F.get(1, iy - 1, iz) - F.get(1, iy + 1, iz)) +
                    (F.get(2, iy - 1, iz) - F.get(2, iy + 1, iz))));

                FXY.set(xBoundM1, iy, iz, 0.25 * invDXDY * (
                    3 * (F.get(ix0, iy - 1, iz) - F.get(ix0, iy + 1, iz))
                    - 4 * (F.get(ix1, iy - 1, iz) - F.get(ix1, iy + 1, iz)) +
                    (F.get(ix2, iy - 1, iz) - F.get(ix2, iy + 1, iz))));
            }

            // interior of z-slice
            for (let iy = 1; iy < yBoundM1; ++iy) {
                for (let ix = 1; ix < xBoundM1; ++ix) {
                    FXY.set(ix, iy, iz, 0.25 * invDXDY * (
                        F.get(ix - 1, iy - 1, iz) - F.get(ix + 1, iy - 1, iz) -
                        F.get(ix - 1, iy + 1, iz) + F.get(ix + 1, iy + 1, iz)));
                }
            }
        }
    }

    private getFXZ(F: Array3<number>, FXZ: Array3<number>): void {
        const xBoundM1 = this.mXBound - 1;
        const zBoundM1 = this.mZBound - 1;
        const ix0 = xBoundM1, ix1 = ix0 - 1, ix2 = ix1 - 1;
        const iz0 = zBoundM1, iz1 = iz0 - 1, iz2 = iz1 - 1;

        const invDXDZ = 1 / (this.mXSpacing * this.mZSpacing);
        for (let iy = 0; iy < this.mYBound; ++iy) {
            // corners of y-slice
            FXZ.set(0, iy, 0, 0.25 * invDXDZ * (
                9 * F.get(0, iy, 0)
                - 12 * F.get(1, iy, 0)
                + 3 * F.get(2, iy, 0)
                - 12 * F.get(0, iy, 1)
                + 16 * F.get(1, iy, 1)
                - 4 * F.get(2, iy, 1)
                + 3 * F.get(0, iy, 2)
                - 4 * F.get(1, iy, 2)
                + F.get(2, iy, 2)));

            FXZ.set(xBoundM1, iy, 0, 0.25 * invDXDZ * (
                9 * F.get(ix0, iy, 0)
                - 12 * F.get(ix1, iy, 0)
                + 3 * F.get(ix2, iy, 0)
                - 12 * F.get(ix0, iy, 1)
                + 16 * F.get(ix1, iy, 1)
                - 4 * F.get(ix2, iy, 1)
                + 3 * F.get(ix0, iy, 2)
                - 4 * F.get(ix1, iy, 2)
                + F.get(ix2, iy, 2)));

            FXZ.set(0, iy, zBoundM1, 0.25 * invDXDZ * (
                9 * F.get(0, iy, iz0)
                - 12 * F.get(1, iy, iz0)
                + 3 * F.get(2, iy, iz0)
                - 12 * F.get(0, iy, iz1)
                + 16 * F.get(1, iy, iz1)
                - 4 * F.get(2, iy, iz1)
                + 3 * F.get(0, iy, iz2)
                - 4 * F.get(1, iy, iz2)
                + F.get(2, iy, iz2)));

            FXZ.set(xBoundM1, iy, zBoundM1, 0.25 * invDXDZ * (
                9 * F.get(ix0, iy, iz0)
                - 12 * F.get(ix1, iy, iz0)
                + 3 * F.get(ix2, iy, iz0)
                - 12 * F.get(ix0, iy, iz1)
                + 16 * F.get(ix1, iy, iz1)
                - 4 * F.get(ix2, iy, iz1)
                + 3 * F.get(ix0, iy, iz2)
                - 4 * F.get(ix1, iy, iz2)
                + F.get(ix2, iy, iz2)));

            // x-edges of y-slice
            for (let ix = 1; ix < xBoundM1; ++ix) {
                FXZ.set(ix, iy, 0, 0.25 * invDXDZ * (
                    3 * (F.get(ix - 1, iy, 0) - F.get(ix + 1, iy, 0)) -
                    4 * (F.get(ix - 1, iy, 1) - F.get(ix + 1, iy, 1)) +
                    (F.get(ix - 1, iy, 2) - F.get(ix + 1, iy, 2))));

                FXZ.set(ix, iy, zBoundM1, 0.25 * invDXDZ * (
                    3 * (F.get(ix - 1, iy, iz0) - F.get(ix + 1, iy, iz0))
                    - 4 * (F.get(ix - 1, iy, iz1) - F.get(ix + 1, iy, iz1)) +
                    (F.get(ix - 1, iy, iz2) - F.get(ix + 1, iy, iz2))));
            }

            // z-edges of y-slice
            for (let iz = 1; iz < zBoundM1; ++iz) {
                FXZ.set(0, iy, iz, 0.25 * invDXDZ * (
                    3 * (F.get(0, iy, iz - 1) - F.get(0, iy, iz + 1)) -
                    4 * (F.get(1, iy, iz - 1) - F.get(1, iy, iz + 1)) +
                    (F.get(2, iy, iz - 1) - F.get(2, iy, iz + 1))));

                FXZ.set(xBoundM1, iy, iz, 0.25 * invDXDZ * (
                    3 * (F.get(ix0, iy, iz - 1) - F.get(ix0, iy, iz + 1))
                    - 4 * (F.get(ix1, iy, iz - 1) - F.get(ix1, iy, iz + 1)) +
                    (F.get(ix2, iy, iz - 1) - F.get(ix2, iy, iz + 1))));
            }

            // interior of y-slice
            for (let iz = 1; iz < zBoundM1; ++iz) {
                for (let ix = 1; ix < xBoundM1; ++ix) {
                    FXZ.set(ix, iy, iz, 0.25 * invDXDZ * (
                        F.get(ix - 1, iy, iz - 1) - F.get(ix + 1, iy, iz - 1) -
                        F.get(ix - 1, iy, iz + 1) + F.get(ix + 1, iy, iz + 1)));
                }
            }
        }
    }

    private getFYZ(F: Array3<number>, FYZ: Array3<number>): void {
        const yBoundM1 = this.mYBound - 1;
        const zBoundM1 = this.mZBound - 1;
        const iy0 = yBoundM1, iy1 = iy0 - 1, iy2 = iy1 - 1;
        const iz0 = zBoundM1, iz1 = iz0 - 1, iz2 = iz1 - 1;

        const invDYDZ = 1 / (this.mYSpacing * this.mZSpacing);
        for (let ix = 0; ix < this.mXBound; ++ix) {
            // corners of x-slice
            FYZ.set(ix, 0, 0, 0.25 * invDYDZ * (
                9 * F.get(ix, 0, 0)
                - 12 * F.get(ix, 1, 0)
                + 3 * F.get(ix, 2, 0)
                - 12 * F.get(ix, 0, 1)
                + 16 * F.get(ix, 1, 1)
                - 4 * F.get(ix, 2, 1)
                + 3 * F.get(ix, 0, 2)
                - 4 * F.get(ix, 1, 2)
                + F.get(ix, 2, 2)));

            FYZ.set(ix, yBoundM1, 0, 0.25 * invDYDZ * (
                9 * F.get(ix, iy0, 0)
                - 12 * F.get(ix, iy1, 0)
                + 3 * F.get(ix, iy2, 0)
                - 12 * F.get(ix, iy0, 1)
                + 16 * F.get(ix, iy1, 1)
                - 4 * F.get(ix, iy2, 1)
                + 3 * F.get(ix, iy0, 2)
                - 4 * F.get(ix, iy1, 2)
                + F.get(ix, iy2, 2)));

            FYZ.set(ix, 0, zBoundM1, 0.25 * invDYDZ * (
                9 * F.get(ix, 0, iz0)
                - 12 * F.get(ix, 1, iz0)
                + 3 * F.get(ix, 2, iz0)
                - 12 * F.get(ix, 0, iz1)
                + 16 * F.get(ix, 1, iz1)
                - 4 * F.get(ix, 2, iz1)
                + 3 * F.get(ix, 0, iz2)
                - 4 * F.get(ix, 1, iz2)
                + F.get(ix, 2, iz2)));

            FYZ.set(ix, yBoundM1, zBoundM1, 0.25 * invDYDZ * (
                9 * F.get(ix, iy0, iz0)
                - 12 * F.get(ix, iy1, iz0)
                + 3 * F.get(ix, iy2, iz0)
                - 12 * F.get(ix, iy0, iz1)
                + 16 * F.get(ix, iy1, iz1)
                - 4 * F.get(ix, iy2, iz1)
                + 3 * F.get(ix, iy0, iz2)
                - 4 * F.get(ix, iy1, iz2)
                + F.get(ix, iy2, iz2)));

            // y-edges of x-slice
            for (let iy = 1; iy < yBoundM1; ++iy) {
                FYZ.set(ix, iy, 0, 0.25 * invDYDZ * (
                    3 * (F.get(ix, iy - 1, 0) - F.get(ix, iy + 1, 0)) -
                    4 * (F.get(ix, iy - 1, 1) - F.get(ix, iy + 1, 1)) +
                    (F.get(ix, iy - 1, 2) - F.get(ix, iy + 1, 2))));

                FYZ.set(ix, iy, zBoundM1, 0.25 * invDYDZ * (
                    3 * (F.get(ix, iy - 1, iz0) - F.get(ix, iy + 1, iz0))
                    - 4 * (F.get(ix, iy - 1, iz1) - F.get(ix, iy + 1, iz1)) +
                    (F.get(ix, iy - 1, iz2) - F.get(ix, iy + 1, iz2))));
            }

            // z-edges of x-slice
            for (let iz = 1; iz < zBoundM1; ++iz) {
                FYZ.set(ix, 0, iz, 0.25 * invDYDZ * (
                    3 * (F.get(ix, 0, iz - 1) - F.get(ix, 0, iz + 1)) -
                    4 * (F.get(ix, 1, iz - 1) - F.get(ix, 1, iz + 1)) +
                    (F.get(ix, 2, iz - 1) - F.get(ix, 2, iz + 1))));

                FYZ.set(ix, yBoundM1, iz, 0.25 * invDYDZ * (
                    3 * (F.get(ix, iy0, iz - 1) - F.get(ix, iy0, iz + 1))
                    - 4 * (F.get(ix, iy1, iz - 1) - F.get(ix, iy1, iz + 1)) +
                    (F.get(ix, iy2, iz - 1) - F.get(ix, iy2, iz + 1))));
            }

            // interior of x-slice
            for (let iz = 1; iz < zBoundM1; ++iz) {
                for (let iy = 1; iy < yBoundM1; ++iy) {
                    FYZ.set(ix, iy, iz, 0.25 * invDYDZ * (
                        F.get(ix, iy - 1, iz - 1) - F.get(ix, iy + 1, iz - 1) -
                        F.get(ix, iy - 1, iz + 1) + F.get(ix, iy + 1, iz + 1)));
                }
            }
        }
    }

    private getFXYZ(F: Array3<number>, FXYZ: Array3<number>): void {
        const xBoundM1 = this.mXBound - 1;
        const yBoundM1 = this.mYBound - 1;
        const zBoundM1 = this.mZBound - 1;

        const invDXDYDZ = 1 / (this.mXSpacing * this.mYSpacing * this.mZSpacing);

        // convolution masks
        //   centered difference, O(h^2)
        const CDer = [-0.5, 0, 0.5];
        //   one-sided difference, O(h^2)
        const ODer = [-1.5, 2, -0.5];
        let mask: number;

        // corners
        FXYZ.set(0, 0, 0, 0);
        FXYZ.set(xBoundM1, 0, 0, 0);
        FXYZ.set(0, yBoundM1, 0, 0);
        FXYZ.set(xBoundM1, yBoundM1, 0, 0);
        FXYZ.set(0, 0, zBoundM1, 0);
        FXYZ.set(xBoundM1, 0, zBoundM1, 0);
        FXYZ.set(0, yBoundM1, zBoundM1, 0);
        FXYZ.set(xBoundM1, yBoundM1, zBoundM1, 0);
        for (let iz = 0; iz <= 2; ++iz) {
            for (let iy = 0; iy <= 2; ++iy) {
                for (let ix = 0; ix <= 2; ++ix) {
                    mask = invDXDYDZ * ODer[ix] * ODer[iy] * ODer[iz];
                    FXYZ.set(0, 0, 0,
                        FXYZ.get(0, 0, 0) + mask * F.get(ix, iy, iz));
                    FXYZ.set(xBoundM1, 0, 0,
                        FXYZ.get(xBoundM1, 0, 0) + mask * F.get(xBoundM1 - ix, iy, iz));
                    FXYZ.set(0, yBoundM1, 0,
                        FXYZ.get(0, yBoundM1, 0) + mask * F.get(ix, yBoundM1 - iy, iz));
                    FXYZ.set(xBoundM1, yBoundM1, 0,
                        FXYZ.get(xBoundM1, yBoundM1, 0) + mask * F.get(xBoundM1 - ix, yBoundM1 - iy, iz));
                    FXYZ.set(0, 0, zBoundM1,
                        FXYZ.get(0, 0, zBoundM1) + mask * F.get(ix, iy, zBoundM1 - iz));
                    FXYZ.set(xBoundM1, 0, zBoundM1,
                        FXYZ.get(xBoundM1, 0, zBoundM1) + mask * F.get(xBoundM1 - ix, iy, zBoundM1 - iz));
                    FXYZ.set(0, yBoundM1, zBoundM1,
                        FXYZ.get(0, yBoundM1, zBoundM1) + mask * F.get(ix, yBoundM1 - iy, zBoundM1 - iz));
                    FXYZ.set(xBoundM1, yBoundM1, zBoundM1,
                        FXYZ.get(xBoundM1, yBoundM1, zBoundM1) + mask * F.get(xBoundM1 - ix, yBoundM1 - iy, zBoundM1 - iz));
                }
            }
        }

        // x-edges
        for (let ix0 = 1; ix0 < xBoundM1; ++ix0) {
            FXYZ.set(ix0, 0, 0, 0);
            FXYZ.set(ix0, yBoundM1, 0, 0);
            FXYZ.set(ix0, 0, zBoundM1, 0);
            FXYZ.set(ix0, yBoundM1, zBoundM1, 0);
            for (let iz = 0; iz <= 2; ++iz) {
                for (let iy = 0; iy <= 2; ++iy) {
                    for (let ix = 0; ix <= 2; ++ix) {
                        mask = invDXDYDZ * CDer[ix] * ODer[iy] * ODer[iz];
                        FXYZ.set(ix0, 0, 0,
                            FXYZ.get(ix0, 0, 0) + mask * F.get(ix0 + ix - 1, iy, iz));
                        FXYZ.set(ix0, yBoundM1, 0,
                            FXYZ.get(ix0, yBoundM1, 0) + mask * F.get(ix0 + ix - 1, yBoundM1 - iy, iz));
                        FXYZ.set(ix0, 0, zBoundM1,
                            FXYZ.get(ix0, 0, zBoundM1) + mask * F.get(ix0 + ix - 1, iy, zBoundM1 - iz));
                        FXYZ.set(ix0, yBoundM1, zBoundM1,
                            FXYZ.get(ix0, yBoundM1, zBoundM1) + mask * F.get(ix0 + ix - 1, yBoundM1 - iy, zBoundM1 - iz));
                    }
                }
            }
        }

        // y-edges
        for (let iy0 = 1; iy0 < yBoundM1; ++iy0) {
            FXYZ.set(0, iy0, 0, 0);
            FXYZ.set(xBoundM1, iy0, 0, 0);
            FXYZ.set(0, iy0, zBoundM1, 0);
            FXYZ.set(xBoundM1, iy0, zBoundM1, 0);
            for (let iz = 0; iz <= 2; ++iz) {
                for (let iy = 0; iy <= 2; ++iy) {
                    for (let ix = 0; ix <= 2; ++ix) {
                        mask = invDXDYDZ * ODer[ix] * CDer[iy] * ODer[iz];
                        FXYZ.set(0, iy0, 0,
                            FXYZ.get(0, iy0, 0) + mask * F.get(ix, iy0 + iy - 1, iz));
                        FXYZ.set(xBoundM1, iy0, 0,
                            FXYZ.get(xBoundM1, iy0, 0) + mask * F.get(xBoundM1 - ix, iy0 + iy - 1, iz));
                        FXYZ.set(0, iy0, zBoundM1,
                            FXYZ.get(0, iy0, zBoundM1) + mask * F.get(ix, iy0 + iy - 1, zBoundM1 - iz));
                        FXYZ.set(xBoundM1, iy0, zBoundM1,
                            FXYZ.get(xBoundM1, iy0, zBoundM1) + mask * F.get(xBoundM1 - ix, iy0 + iy - 1, zBoundM1 - iz));
                    }
                }
            }
        }

        // z-edges
        for (let iz0 = 1; iz0 < zBoundM1; ++iz0) {
            FXYZ.set(0, 0, iz0, 0);
            FXYZ.set(xBoundM1, 0, iz0, 0);
            FXYZ.set(0, yBoundM1, iz0, 0);
            FXYZ.set(xBoundM1, yBoundM1, iz0, 0);
            for (let iz = 0; iz <= 2; ++iz) {
                for (let iy = 0; iy <= 2; ++iy) {
                    for (let ix = 0; ix <= 2; ++ix) {
                        mask = invDXDYDZ * ODer[ix] * ODer[iy] * CDer[iz];
                        FXYZ.set(0, 0, iz0,
                            FXYZ.get(0, 0, iz0) + mask * F.get(ix, iy, iz0 + iz - 1));
                        FXYZ.set(xBoundM1, 0, iz0,
                            FXYZ.get(xBoundM1, 0, iz0) + mask * F.get(xBoundM1 - ix, iy, iz0 + iz - 1));
                        FXYZ.set(0, yBoundM1, iz0,
                            FXYZ.get(0, yBoundM1, iz0) + mask * F.get(ix, yBoundM1 - iy, iz0 + iz - 1));
                        FXYZ.set(xBoundM1, yBoundM1, iz0,
                            FXYZ.get(xBoundM1, yBoundM1, iz0) + mask * F.get(xBoundM1 - ix, yBoundM1 - iy, iz0 + iz - 1));
                    }
                }
            }
        }

        // xy-faces
        for (let iy0 = 1; iy0 < yBoundM1; ++iy0) {
            for (let ix0 = 1; ix0 < xBoundM1; ++ix0) {
                FXYZ.set(ix0, iy0, 0, 0);
                FXYZ.set(ix0, iy0, zBoundM1, 0);
                for (let iz = 0; iz <= 2; ++iz) {
                    for (let iy = 0; iy <= 2; ++iy) {
                        for (let ix = 0; ix <= 2; ++ix) {
                            mask = invDXDYDZ * CDer[ix] * CDer[iy] * ODer[iz];
                            FXYZ.set(ix0, iy0, 0,
                                FXYZ.get(ix0, iy0, 0) + mask * F.get(ix0 + ix - 1, iy0 + iy - 1, iz));
                            FXYZ.set(ix0, iy0, zBoundM1,
                                FXYZ.get(ix0, iy0, zBoundM1) + mask * F.get(ix0 + ix - 1, iy0 + iy - 1, zBoundM1 - iz));
                        }
                    }
                }
            }
        }

        // xz-faces
        for (let iz0 = 1; iz0 < zBoundM1; ++iz0) {
            for (let ix0 = 1; ix0 < xBoundM1; ++ix0) {
                FXYZ.set(ix0, 0, iz0, 0);
                FXYZ.set(ix0, yBoundM1, iz0, 0);
                for (let iz = 0; iz <= 2; ++iz) {
                    for (let iy = 0; iy <= 2; ++iy) {
                        for (let ix = 0; ix <= 2; ++ix) {
                            mask = invDXDYDZ * CDer[ix] * ODer[iy] * CDer[iz];
                            FXYZ.set(ix0, 0, iz0,
                                FXYZ.get(ix0, 0, iz0) + mask * F.get(ix0 + ix - 1, iy, iz0 + iz - 1));
                            FXYZ.set(ix0, yBoundM1, iz0,
                                FXYZ.get(ix0, yBoundM1, iz0) + mask * F.get(ix0 + ix - 1, yBoundM1 - iy, iz0 + iz - 1));
                        }
                    }
                }
            }
        }

        // yz-faces
        for (let iz0 = 1; iz0 < zBoundM1; ++iz0) {
            for (let iy0 = 1; iy0 < yBoundM1; ++iy0) {
                FXYZ.set(0, iy0, iz0, 0);
                FXYZ.set(xBoundM1, iy0, iz0, 0);
                for (let iz = 0; iz <= 2; ++iz) {
                    for (let iy = 0; iy <= 2; ++iy) {
                        for (let ix = 0; ix <= 2; ++ix) {
                            mask = invDXDYDZ * ODer[ix] * CDer[iy] * CDer[iz];
                            FXYZ.set(0, iy0, iz0,
                                FXYZ.get(0, iy0, iz0) + mask * F.get(ix, iy0 + iy - 1, iz0 + iz - 1));
                            FXYZ.set(xBoundM1, iy0, iz0,
                                FXYZ.get(xBoundM1, iy0, iz0) + mask * F.get(xBoundM1 - ix, iy0 + iy - 1, iz0 + iz - 1));
                        }
                    }
                }
            }
        }

        // interiors
        for (let iz0 = 1; iz0 < zBoundM1; ++iz0) {
            for (let iy0 = 1; iy0 < yBoundM1; ++iy0) {
                for (let ix0 = 1; ix0 < xBoundM1; ++ix0) {
                    FXYZ.set(ix0, iy0, iz0, 0);

                    for (let iz = 0; iz <= 2; ++iz) {
                        for (let iy = 0; iy <= 2; ++iy) {
                            for (let ix = 0; ix <= 2; ++ix) {
                                mask = invDXDYDZ * CDer[ix] * CDer[iy] * CDer[iz];
                                FXYZ.set(ix0, iy0, iz0,
                                    FXYZ.get(ix0, iy0, iz0) +
                                    mask * F.get(ix0 + ix - 1, iy0 + iy - 1, iz0 + iz - 1));
                            }
                        }
                    }
                }
            }
        }
    }

    private getPolynomials(F: Array3<number>, FX: Array3<number>,
        FY: Array3<number>, FZ: Array3<number>, FXY: Array3<number>,
        FXZ: Array3<number>, FYZ: Array3<number>, FXYZ: Array3<number>): void {
        const xBoundM1 = this.mXBound - 1;
        const yBoundM1 = this.mYBound - 1;
        const zBoundM1 = this.mZBound - 1;
        for (let iz = 0; iz < zBoundM1; ++iz) {
            for (let iy = 0; iy < yBoundM1; ++iy) {
                for (let ix = 0; ix < xBoundM1; ++ix) {
                    // Note the 'transposing' of the 2x2x2 blocks (to match
                    // the notation used in the polynomial definition).
                    const G = extractBlock(F, ix, iy, iz);
                    const GX = extractBlock(FX, ix, iy, iz);
                    const GY = extractBlock(FY, ix, iy, iz);
                    const GZ = extractBlock(FZ, ix, iy, iz);
                    const GXY = extractBlock(FXY, ix, iy, iz);
                    const GXZ = extractBlock(FXZ, ix, iy, iz);
                    const GYZ = extractBlock(FYZ, ix, iy, iz);
                    const GXYZ = extractBlock(FXYZ, ix, iy, iz);

                    this.construct(this.mPoly.get(ix, iy, iz),
                        G, GX, GY, GZ, GXY, GXZ, GYZ, GXYZ);
                }
            }
        }
    }

    private static computeDerivative(s0: number, s1: number, s2: number, s3: number): number {
        if (s1 !== s2) {
            if (s0 !== s1) {
                if (s2 !== s3) {
                    const ad0 = Math.abs(s3 - s2);
                    const ad1 = Math.abs(s0 - s1);
                    return (ad0 * s1 + ad1 * s2) / (ad0 + ad1);
                } else {
                    return s2;
                }
            } else {
                if (s2 !== s3) {
                    return s1;
                } else {
                    return 0.5 * (s1 + s2);
                }
            }
        } else {
            return s1;
        }
    }

    private construct(poly: AkimaPolynomial3,
        F: Float64Array, FX: Float64Array, FY: Float64Array, FZ: Float64Array,
        FXY: Float64Array, FXZ: Float64Array, FYZ: Float64Array,
        FXYZ: Float64Array): void {
        const dx = this.mXSpacing, dy = this.mYSpacing, dz = this.mZSpacing;
        const invDX = 1 / dx, invDX2 = invDX * invDX;
        const invDY = 1 / dy, invDY2 = invDY * invDY;
        const invDZ = 1 / dz, invDZ2 = invDZ * invDZ;
        let b0: number, b1: number, b2: number, b3: number;
        let b4: number, b5: number, b6: number, b7: number;

        poly.setA(0, 0, 0, F[B000]);
        poly.setA(1, 0, 0, FX[B000]);
        poly.setA(0, 1, 0, FY[B000]);
        poly.setA(0, 0, 1, FZ[B000]);
        poly.setA(1, 1, 0, FXY[B000]);
        poly.setA(1, 0, 1, FXZ[B000]);
        poly.setA(0, 1, 1, FYZ[B000]);
        poly.setA(1, 1, 1, FXYZ[B000]);

        // solve for Aij0
        b0 = (F[B100] - poly.evaluate(0, 0, 0, dx, 0, 0)) * invDX2;
        b1 = (FX[B100] - poly.evaluate(1, 0, 0, dx, 0, 0)) * invDX;
        poly.setA(2, 0, 0, 3 * b0 - b1);
        poly.setA(3, 0, 0, (-2 * b0 + b1) * invDX);

        b0 = (F[B010] - poly.evaluate(0, 0, 0, 0, dy, 0)) * invDY2;
        b1 = (FY[B010] - poly.evaluate(0, 1, 0, 0, dy, 0)) * invDY;
        poly.setA(0, 2, 0, 3 * b0 - b1);
        poly.setA(0, 3, 0, (-2 * b0 + b1) * invDY);

        b0 = (FY[B100] - poly.evaluate(0, 1, 0, dx, 0, 0)) * invDX2;
        b1 = (FXY[B100] - poly.evaluate(1, 1, 0, dx, 0, 0)) * invDX;
        poly.setA(2, 1, 0, 3 * b0 - b1);
        poly.setA(3, 1, 0, (-2 * b0 + b1) * invDX);

        b0 = (FX[B010] - poly.evaluate(1, 0, 0, 0, dy, 0)) * invDY2;
        b1 = (FXY[B010] - poly.evaluate(1, 1, 0, 0, dy, 0)) * invDY;
        poly.setA(1, 2, 0, 3 * b0 - b1);
        poly.setA(1, 3, 0, (-2 * b0 + b1) * invDY);

        b0 = (F[B110] - poly.evaluate(0, 0, 0, dx, dy, 0)) * invDX2 * invDY2;
        b1 = (FX[B110] - poly.evaluate(1, 0, 0, dx, dy, 0)) * invDX * invDY2;
        b2 = (FY[B110] - poly.evaluate(0, 1, 0, dx, dy, 0)) * invDX2 * invDY;
        b3 = (FXY[B110] - poly.evaluate(1, 1, 0, dx, dy, 0)) * invDX * invDY;
        poly.setA(2, 2, 0, 9 * b0 - 3 * b1 - 3 * b2 + b3);
        poly.setA(3, 2, 0, (-6 * b0 + 3 * b1 + 2 * b2 - b3) * invDX);
        poly.setA(2, 3, 0, (-6 * b0 + 2 * b1 + 3 * b2 - b3) * invDY);
        poly.setA(3, 3, 0, (4 * b0 - 2 * b1 - 2 * b2 + b3) * invDX * invDY);

        // solve for Ai0k
        b0 = (F[B001] - poly.evaluate(0, 0, 0, 0, 0, dz)) * invDZ2;
        b1 = (FZ[B001] - poly.evaluate(0, 0, 1, 0, 0, dz)) * invDZ;
        poly.setA(0, 0, 2, 3 * b0 - b1);
        poly.setA(0, 0, 3, (-2 * b0 + b1) * invDZ);

        b0 = (FZ[B100] - poly.evaluate(0, 0, 1, dx, 0, 0)) * invDX2;
        b1 = (FXZ[B100] - poly.evaluate(1, 0, 1, dx, 0, 0)) * invDX;
        poly.setA(2, 0, 1, 3 * b0 - b1);
        poly.setA(3, 0, 1, (-2 * b0 + b1) * invDX);

        b0 = (FX[B001] - poly.evaluate(1, 0, 0, 0, 0, dz)) * invDZ2;
        b1 = (FXZ[B001] - poly.evaluate(1, 0, 1, 0, 0, dz)) * invDZ;
        poly.setA(1, 0, 2, 3 * b0 - b1);
        poly.setA(1, 0, 3, (-2 * b0 + b1) * invDZ);

        b0 = (F[B101] - poly.evaluate(0, 0, 0, dx, 0, dz)) * invDX2 * invDZ2;
        b1 = (FX[B101] - poly.evaluate(1, 0, 0, dx, 0, dz)) * invDX * invDZ2;
        b2 = (FZ[B101] - poly.evaluate(0, 0, 1, dx, 0, dz)) * invDX2 * invDZ;
        b3 = (FXZ[B101] - poly.evaluate(1, 0, 1, dx, 0, dz)) * invDX * invDZ;
        poly.setA(2, 0, 2, 9 * b0 - 3 * b1 - 3 * b2 + b3);
        poly.setA(3, 0, 2, (-6 * b0 + 3 * b1 + 2 * b2 - b3) * invDX);
        poly.setA(2, 0, 3, (-6 * b0 + 2 * b1 + 3 * b2 - b3) * invDZ);
        poly.setA(3, 0, 3, (4 * b0 - 2 * b1 - 2 * b2 + b3) * invDX * invDZ);

        // solve for A0jk
        b0 = (FZ[B010] - poly.evaluate(0, 0, 1, 0, dy, 0)) * invDY2;
        b1 = (FYZ[B010] - poly.evaluate(0, 1, 1, 0, dy, 0)) * invDY;
        poly.setA(0, 2, 1, 3 * b0 - b1);
        poly.setA(0, 3, 1, (-2 * b0 + b1) * invDY);

        b0 = (FY[B001] - poly.evaluate(0, 1, 0, 0, 0, dz)) * invDZ2;
        b1 = (FYZ[B001] - poly.evaluate(0, 1, 1, 0, 0, dz)) * invDZ;
        poly.setA(0, 1, 2, 3 * b0 - b1);
        poly.setA(0, 1, 3, (-2 * b0 + b1) * invDZ);

        b0 = (F[B011] - poly.evaluate(0, 0, 0, 0, dy, dz)) * invDY2 * invDZ2;
        b1 = (FY[B011] - poly.evaluate(0, 1, 0, 0, dy, dz)) * invDY * invDZ2;
        b2 = (FZ[B011] - poly.evaluate(0, 0, 1, 0, dy, dz)) * invDY2 * invDZ;
        b3 = (FYZ[B011] - poly.evaluate(0, 1, 1, 0, dy, dz)) * invDY * invDZ;
        poly.setA(0, 2, 2, 9 * b0 - 3 * b1 - 3 * b2 + b3);
        poly.setA(0, 3, 2, (-6 * b0 + 3 * b1 + 2 * b2 - b3) * invDY);
        poly.setA(0, 2, 3, (-6 * b0 + 2 * b1 + 3 * b2 - b3) * invDZ);
        poly.setA(0, 3, 3, (4 * b0 - 2 * b1 - 2 * b2 + b3) * invDY * invDZ);

        // solve for Aij1
        b0 = (FYZ[B100] - poly.evaluate(0, 1, 1, dx, 0, 0)) * invDX2;
        b1 = (FXYZ[B100] - poly.evaluate(1, 1, 1, dx, 0, 0)) * invDX;
        poly.setA(2, 1, 1, 3 * b0 - b1);
        poly.setA(3, 1, 1, (-2 * b0 + b1) * invDX);

        b0 = (FXZ[B010] - poly.evaluate(1, 0, 1, 0, dy, 0)) * invDY2;
        b1 = (FXYZ[B010] - poly.evaluate(1, 1, 1, 0, dy, 0)) * invDY;
        poly.setA(1, 2, 1, 3 * b0 - b1);
        poly.setA(1, 3, 1, (-2 * b0 + b1) * invDY);

        b0 = (FZ[B110] - poly.evaluate(0, 0, 1, dx, dy, 0)) * invDX2 * invDY2;
        b1 = (FXZ[B110] - poly.evaluate(1, 0, 1, dx, dy, 0)) * invDX * invDY2;
        b2 = (FYZ[B110] - poly.evaluate(0, 1, 1, dx, dy, 0)) * invDX2 * invDY;
        b3 = (FXYZ[B110] - poly.evaluate(1, 1, 1, dx, dy, 0)) * invDX * invDY;
        poly.setA(2, 2, 1, 9 * b0 - 3 * b1 - 3 * b2 + b3);
        poly.setA(3, 2, 1, (-6 * b0 + 3 * b1 + 2 * b2 - b3) * invDX);
        poly.setA(2, 3, 1, (-6 * b0 + 2 * b1 + 3 * b2 - b3) * invDY);
        poly.setA(3, 3, 1, (4 * b0 - 2 * b1 - 2 * b2 + b3) * invDX * invDY);

        // solve for Ai1k
        b0 = (FXY[B001] - poly.evaluate(1, 1, 0, 0, 0, dz)) * invDZ2;
        b1 = (FXYZ[B001] - poly.evaluate(1, 1, 1, 0, 0, dz)) * invDZ;
        poly.setA(1, 1, 2, 3 * b0 - b1);
        poly.setA(1, 1, 3, (-2 * b0 + b1) * invDZ);

        b0 = (FY[B101] - poly.evaluate(0, 1, 0, dx, 0, dz)) * invDX2 * invDZ2;
        b1 = (FXY[B101] - poly.evaluate(1, 1, 0, dx, 0, dz)) * invDX * invDZ2;
        b2 = (FYZ[B101] - poly.evaluate(0, 1, 1, dx, 0, dz)) * invDX2 * invDZ;
        b3 = (FXYZ[B101] - poly.evaluate(1, 1, 1, dx, 0, dz)) * invDX * invDZ;
        poly.setA(2, 1, 2, 9 * b0 - 3 * b1 - 3 * b2 + b3);
        poly.setA(3, 1, 2, (-6 * b0 + 3 * b1 + 2 * b2 - b3) * invDX);
        poly.setA(2, 1, 3, (-6 * b0 + 2 * b1 + 3 * b2 - b3) * invDZ);
        poly.setA(3, 1, 3, (4 * b0 - 2 * b1 - 2 * b2 + b3) * invDX * invDZ);

        // solve for A1jk
        b0 = (FX[B011] - poly.evaluate(1, 0, 0, 0, dy, dz)) * invDY2 * invDZ2;
        b1 = (FXY[B011] - poly.evaluate(1, 1, 0, 0, dy, dz)) * invDY * invDZ2;
        b2 = (FXZ[B011] - poly.evaluate(1, 0, 1, 0, dy, dz)) * invDY2 * invDZ;
        b3 = (FXYZ[B011] - poly.evaluate(1, 1, 1, 0, dy, dz)) * invDY * invDZ;
        poly.setA(1, 2, 2, 9 * b0 - 3 * b1 - 3 * b2 + b3);
        poly.setA(1, 3, 2, (-6 * b0 + 3 * b1 + 2 * b2 - b3) * invDY);
        poly.setA(1, 2, 3, (-6 * b0 + 2 * b1 + 3 * b2 - b3) * invDZ);
        poly.setA(1, 3, 3, (4 * b0 - 2 * b1 - 2 * b2 + b3) * invDY * invDZ);

        // solve for remaining Aijk with i >= 2, j >= 2, k >= 2
        b0 = (F[B111] - poly.evaluate(0, 0, 0, dx, dy, dz)) * invDX2 * invDY2 * invDZ2;
        b1 = (FX[B111] - poly.evaluate(1, 0, 0, dx, dy, dz)) * invDX * invDY2 * invDZ2;
        b2 = (FY[B111] - poly.evaluate(0, 1, 0, dx, dy, dz)) * invDX2 * invDY * invDZ2;
        b3 = (FZ[B111] - poly.evaluate(0, 0, 1, dx, dy, dz)) * invDX2 * invDY2 * invDZ;
        b4 = (FXY[B111] - poly.evaluate(1, 1, 0, dx, dy, dz)) * invDX * invDY * invDZ2;
        b5 = (FXZ[B111] - poly.evaluate(1, 0, 1, dx, dy, dz)) * invDX * invDY2 * invDZ;
        b6 = (FYZ[B111] - poly.evaluate(0, 1, 1, dx, dy, dz)) * invDX2 * invDY * invDZ;
        b7 = (FXYZ[B111] - poly.evaluate(1, 1, 1, dx, dy, dz)) * invDX * invDY * invDZ;
        poly.setA(2, 2, 2, 27 * b0 - 9 * b1 - 9 * b2 -
            9 * b3 + 3 * b4 + 3 * b5 + 3 * b6 - b7);
        poly.setA(3, 2, 2, (-18 * b0 + 9 * b1 + 6 * b2 +
            6 * b3 - 3 * b4 - 3 * b5 - 2 * b6 + b7) * invDX);
        poly.setA(2, 3, 2, (-18 * b0 + 6 * b1 + 9 * b2 +
            6 * b3 - 3 * b4 - 2 * b5 - 3 * b6 + b7) * invDY);
        poly.setA(2, 2, 3, (-18 * b0 + 6 * b1 + 6 * b2 +
            9 * b3 - 2 * b4 - 3 * b5 - 3 * b6 + b7) * invDZ);
        poly.setA(3, 3, 2, (12 * b0 - 6 * b1 - 6 * b2 -
            4 * b3 + 3 * b4 + 2 * b5 + 2 * b6 - b7) * invDX * invDY);
        poly.setA(3, 2, 3, (12 * b0 - 6 * b1 - 4 * b2 -
            6 * b3 + 2 * b4 + 3 * b5 + 2 * b6 - b7) * invDX * invDZ);
        poly.setA(2, 3, 3, (12 * b0 - 4 * b1 - 6 * b2 -
            6 * b3 + 2 * b4 + 2 * b5 + 3 * b6 - b7) * invDY * invDZ);
        poly.setA(3, 3, 3, (-8 * b0 + 4 * b1 + 4 * b2 +
            4 * b3 - 2 * b4 - 2 * b5 - 2 * b6 + b7) * invDX * invDY * invDZ);
    }

    private xLookup(x: number): { index: number; d: number } {
        let xIndex = 0;
        for (let xIndexP1 = 1; xIndexP1 < this.mXBound; ++xIndex, ++xIndexP1) {
            if (x < this.mXMin + this.mXSpacing * xIndexP1) {
                return { index: xIndex, d: x - (this.mXMin + this.mXSpacing * xIndex) };
            }
        }

        --xIndex;
        return { index: xIndex, d: x - (this.mXMin + this.mXSpacing * xIndex) };
    }

    private yLookup(y: number): { index: number; d: number } {
        let yIndex = 0;
        for (let yIndexP1 = 1; yIndexP1 < this.mYBound; ++yIndex, ++yIndexP1) {
            if (y < this.mYMin + this.mYSpacing * yIndexP1) {
                return { index: yIndex, d: y - (this.mYMin + this.mYSpacing * yIndex) };
            }
        }

        --yIndex;
        return { index: yIndex, d: y - (this.mYMin + this.mYSpacing * yIndex) };
    }

    private zLookup(z: number): { index: number; d: number } {
        let zIndex = 0;
        for (let zIndexP1 = 1; zIndexP1 < this.mZBound; ++zIndex, ++zIndexP1) {
            if (z < this.mZMin + this.mZSpacing * zIndexP1) {
                return { index: zIndex, d: z - (this.mZMin + this.mZSpacing * zIndex) };
            }
        }

        --zIndex;
        return { index: zIndex, d: z - (this.mZMin + this.mZSpacing * zIndex) };
    }
}
