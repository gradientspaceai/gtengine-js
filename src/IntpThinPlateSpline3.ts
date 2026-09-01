// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpThinPlateSpline3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// WARNING. The implementation allows you to transform the inputs (x,y,z) to
// the unit cube and perform the interpolation in that space. The idea is to
// keep the floating-point numbers to order 1 for numerical stability of the
// algorithm. The classical thin-plate spline algorithm does not include this
// transformation. The interpolation is invariant to translations and
// rotations of (x,y,z) but not to scaling. The following document is about
// thin plate splines.
//   https://www.geometrictools.com/Documentation/ThinPlateSplines.pdf
//
// Port notes:
// - 'operator()(Real,Real,Real)' becomes 'evaluate(x, y, z)', matching the
//   Intp* precedent (IntpTricubic3). 'IsInitialized' becomes
//   'isInitialized()'. When the spline was not initialized, evaluate returns
//   Number.MAX_VALUE (the port of std::numeric_limits<Real>::max()).
// - The upstream WARNING comment describes the invariance of the classical
//   (untransformed) spline. With transformToUnitCube = true the spline is
//   invariant to translation and to uniform scaling but not to rotation,
//   i.e. the comment does not describe the transformed configuration. The
//   port preserves the behavior and only documents the mismatch.

import { inverse, multiplyAB, multiplyATB } from './Matrix';
import { GMatrix } from './GMatrix';
import { logAssert } from './Logger';

// Kernel(t) = -|t|
export function intpThinPlateSpline3Kernel(t: number): number {
    return -Math.abs(t);
}

export class IntpThinPlateSpline3 {
    // Input data.
    private readonly mNumPoints: number;
    private readonly mX: number[];
    private readonly mY: number[];
    private readonly mZ: number[];
    private readonly mSmooth: number;

    // Thin plate spline coefficients. The A[] coefficients are associated
    // with the Green's functions G(x,y,z,*) and the B[] coefficients are
    // associated with the affine term B[0] + B[1]*x + B[2]*y + B[3]*z.
    private readonly mA: number[];
    private readonly mB: number[];

    // Extent of input data.
    private mXMin: number;
    private mXMax: number;
    private mXInvRange: number;
    private mYMin: number;
    private mYMax: number;
    private mYInvRange: number;
    private mZMin: number;
    private mZMax: number;
    private mZInvRange: number;

    private mInitialized: boolean;

    // Construction. Data points are (x,y,z,f(x,y,z)). The smoothing
    // parameter must be nonnegative.
    constructor(numPoints: number, X: ArrayLike<number>, Y: ArrayLike<number>,
        Z: ArrayLike<number>, F: ArrayLike<number>, smooth: number,
        transformToUnitCube: boolean) {
        logAssert(numPoints >= 4 && X.length >= numPoints &&
            Y.length >= numPoints && Z.length >= numPoints &&
            F.length >= numPoints && smooth >= 0, 'Invalid input.');

        this.mNumPoints = numPoints;
        this.mX = new Array<number>(numPoints).fill(0);
        this.mY = new Array<number>(numPoints).fill(0);
        this.mZ = new Array<number>(numPoints).fill(0);
        this.mSmooth = smooth;
        this.mA = new Array<number>(numPoints).fill(0);
        this.mB = [0, 0, 0, 0];
        this.mInitialized = false;

        let i: number, row: number, col: number;

        if (transformToUnitCube) {
            // Map input (x,y,z) to the unit cube. This is not part of the
            // classical thin-plate spline algorithm, because the
            // interpolation is not invariant to scalings.
            this.mXMin = X[0];
            this.mXMax = X[0];
            this.mYMin = Y[0];
            this.mYMax = Y[0];
            this.mZMin = Z[0];
            this.mZMax = Z[0];
            for (i = 1; i < numPoints; ++i) {
                if (X[i] < this.mXMin) { this.mXMin = X[i]; }
                if (X[i] > this.mXMax) { this.mXMax = X[i]; }
                if (Y[i] < this.mYMin) { this.mYMin = Y[i]; }
                if (Y[i] > this.mYMax) { this.mYMax = Y[i]; }
                if (Z[i] < this.mZMin) { this.mZMin = Z[i]; }
                if (Z[i] > this.mZMax) { this.mZMax = Z[i]; }
            }
            this.mXInvRange = 1 / (this.mXMax - this.mXMin);
            this.mYInvRange = 1 / (this.mYMax - this.mYMin);
            this.mZInvRange = 1 / (this.mZMax - this.mZMin);
            for (i = 0; i < numPoints; ++i) {
                this.mX[i] = (X[i] - this.mXMin) * this.mXInvRange;
                this.mY[i] = (Y[i] - this.mYMin) * this.mYInvRange;
                this.mZ[i] = (Z[i] - this.mZMin) * this.mZInvRange;
            }
        } else {
            // The classical thin-plate spline uses the data as is. The
            // values mXMax, mYMax, and mZMax are not used, but they are
            // initialized anyway (to irrelevant numbers).
            this.mXMin = 0;
            this.mXMax = 1;
            this.mXInvRange = 1;
            this.mYMin = 0;
            this.mYMax = 1;
            this.mYInvRange = 1;
            this.mZMin = 0;
            this.mZMax = 1;
            this.mZInvRange = 1;
            for (i = 0; i < numPoints; ++i) {
                this.mX[i] = X[i];
                this.mY[i] = Y[i];
                this.mZ[i] = Z[i];
            }
        }

        // Compute matrix A = M + lambda*I [NxN matrix].
        const AMat = new GMatrix(numPoints, numPoints);
        for (row = 0; row < numPoints; ++row) {
            for (col = 0; col < numPoints; ++col) {
                if (row === col) {
                    AMat.set(row, col, this.mSmooth);
                } else {
                    const dx = this.mX[row] - this.mX[col];
                    const dy = this.mY[row] - this.mY[col];
                    const dz = this.mZ[row] - this.mZ[col];
                    const t = Math.sqrt(dx * dx + dy * dy + dz * dz);
                    AMat.set(row, col, intpThinPlateSpline3Kernel(t));
                }
            }
        }

        // Compute matrix B [Nx4 matrix].
        const BMat = new GMatrix(numPoints, 4);
        for (row = 0; row < numPoints; ++row) {
            BMat.set(row, 0, 1);
            BMat.set(row, 1, this.mX[row]);
            BMat.set(row, 2, this.mY[row]);
            BMat.set(row, 3, this.mZ[row]);
        }

        // Compute A^{-1}.
        const invA = inverse(AMat);
        if (!invA.invertible) {
            return;
        }
        const invAMat = invA.inverse;

        // Compute P = B^t A^{-1} [4xN matrix].
        const PMat = multiplyATB(BMat, invAMat);

        // Compute Q = P B = B^t A^{-1} B [4x4 matrix].
        const QMat = multiplyAB(PMat, BMat);

        // Compute Q^{-1}.
        const invQ = inverse(QMat);
        if (!invQ.invertible) {
            return;
        }
        const invQMat = invQ.inverse;

        // Compute P*w.
        const prod = [0, 0, 0, 0];
        for (row = 0; row < 4; ++row) {
            prod[row] = 0;
            for (i = 0; i < numPoints; ++i) {
                prod[row] += PMat.get(row, i) * F[i];
            }
        }

        // Compute 'b' vector for smooth thin plate spline.
        for (row = 0; row < 4; ++row) {
            this.mB[row] = 0;
            for (i = 0; i < 4; ++i) {
                this.mB[row] += invQMat.get(row, i) * prod[i];
            }
        }

        // Compute w - B*b.
        const tmp = new Array<number>(numPoints).fill(0);
        for (row = 0; row < numPoints; ++row) {
            tmp[row] = F[row];
            for (i = 0; i < 4; ++i) {
                tmp[row] -= BMat.get(row, i) * this.mB[i];
            }
        }

        // Compute 'a' vector for smooth thin plate spline.
        for (row = 0; row < numPoints; ++row) {
            this.mA[row] = 0;
            for (i = 0; i < numPoints; ++i) {
                this.mA[row] += invAMat.get(row, i) * tmp[i];
            }
        }

        this.mInitialized = true;
    }

    // Check this after the constructor call to see whether the thin plate
    // spline coefficients were successfully computed. If so, then calls to
    // evaluate(x, y, z) will work properly.
    isInitialized(): boolean {
        return this.mInitialized;
    }

    // Evaluate the interpolator. If isInitialized() returns false, this
    // returns Number.MAX_VALUE.
    evaluate(x: number, y: number, z: number): number {
        if (this.mInitialized) {
            // Map (x,y,z) to the unit cube.
            x = (x - this.mXMin) * this.mXInvRange;
            y = (y - this.mYMin) * this.mYInvRange;
            z = (z - this.mZMin) * this.mZInvRange;

            let result = this.mB[0] + this.mB[1] * x + this.mB[2] * y +
                this.mB[3] * z;
            for (let i = 0; i < this.mNumPoints; ++i) {
                const dx = x - this.mX[i];
                const dy = y - this.mY[i];
                const dz = z - this.mZ[i];
                const t = Math.sqrt(dx * dx + dy * dy + dz * dz);
                result += this.mA[i] * intpThinPlateSpline3Kernel(t);
            }
            return result;
        }

        return Number.MAX_VALUE;
    }

    // Compute the functional value a^T*M*a when lambda is zero or
    // lambda*w^T*(M+lambda*I)*w when lambda is positive. See the thin plate
    // splines PDF for a description of these quantities.
    computeFunctional(): number {
        let functional = 0;
        for (let row = 0; row < this.mNumPoints; ++row) {
            for (let col = 0; col < this.mNumPoints; ++col) {
                if (row === col) {
                    functional += this.mSmooth * this.mA[row] * this.mA[col];
                } else {
                    const dx = this.mX[row] - this.mX[col];
                    const dy = this.mY[row] - this.mY[col];
                    const dz = this.mZ[row] - this.mZ[col];
                    const t = Math.sqrt(dx * dx + dy * dy + dz * dz);
                    functional += intpThinPlateSpline3Kernel(t) *
                        this.mA[row] * this.mA[col];
                }
            }
        }

        if (this.mSmooth > 0) {
            functional *= this.mSmooth;
        }

        return functional;
    }
}
