// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpThinPlateSpline2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// WARNING. The implementation allows you to transform the inputs (x,y) to
// the unit square and perform the interpolation in that space. The idea is
// to keep the floating-point numbers to order 1 for numerical stability of
// the algorithm. The classical thin-plate spline algorithm does not include
// this transformation. The interpolation is invariant to translations and
// rotations of (x,y) but not to scaling. The following document is about
// thin plate splines.
//   https://www.geometrictools.com/Documentation/ThinPlateSplines.pdf
//
// Port notes:
// - 'operator()(Real,Real)' becomes 'evaluate(x, y)', matching the Intp*
//   precedent (IntpTricubic3). 'IsInitialized' becomes 'isInitialized()'.
//   When the spline was not initialized, evaluate returns Number.MAX_VALUE
//   (the port of std::numeric_limits<Real>::max()).
// - The upstream WARNING comment describes the invariance of the classical
//   (untransformed) spline. With transformToUnitSquare = true the spline is
//   invariant to translation and to uniform scaling but not to rotation,
//   i.e. the comment does not describe the transformed configuration. The
//   port preserves the behavior and only documents the mismatch.

import { inverse, multiplyAB, multiplyATB } from './Matrix.js';
import { GMatrix } from './GMatrix.js';
import { logAssert } from './Logger.js';

// Kernel(t) = t^2 * log(t^2)
export function intpThinPlateSpline2Kernel(t: number): number {
    if (t > 0) {
        const t2 = t * t;
        return t2 * Math.log(t2);
    }
    return 0;
}

export class IntpThinPlateSpline2 {
    // Input data.
    private readonly mNumPoints: number;
    private readonly mX: number[];
    private readonly mY: number[];
    private readonly mSmooth: number;

    // Thin plate spline coefficients. The A[] coefficients are associated
    // with the Green's functions G(x,y,*) and the B[] coefficients are
    // associated with the affine term B[0] + B[1]*x + B[2]*y.
    private readonly mA: number[];
    private readonly mB: number[];

    // Extent of input data.
    private mXMin: number;
    private mXMax: number;
    private mXInvRange: number;
    private mYMin: number;
    private mYMax: number;
    private mYInvRange: number;

    private mInitialized: boolean;

    // Construction. Data points are (x,y,f(x,y)). The smoothing parameter
    // must be nonnegative.
    constructor(numPoints: number, X: ArrayLike<number>, Y: ArrayLike<number>,
        F: ArrayLike<number>, smooth: number, transformToUnitSquare: boolean) {
        logAssert(numPoints >= 3 && X.length >= numPoints &&
            Y.length >= numPoints && F.length >= numPoints && smooth >= 0,
            'Invalid input.');

        this.mNumPoints = numPoints;
        this.mX = new Array<number>(numPoints).fill(0);
        this.mY = new Array<number>(numPoints).fill(0);
        this.mSmooth = smooth;
        this.mA = new Array<number>(numPoints).fill(0);
        this.mB = [0, 0, 0];
        this.mInitialized = false;

        let i: number, row: number, col: number;

        if (transformToUnitSquare) {
            // Map input (x,y) to the unit square. This is not part of the
            // classical thin-plate spline algorithm because the
            // interpolation is not invariant to scalings.
            this.mXMin = X[0];
            this.mXMax = X[0];
            this.mYMin = Y[0];
            this.mYMax = Y[0];
            for (i = 1; i < numPoints; ++i) {
                if (X[i] < this.mXMin) { this.mXMin = X[i]; }
                if (X[i] > this.mXMax) { this.mXMax = X[i]; }
                if (Y[i] < this.mYMin) { this.mYMin = Y[i]; }
                if (Y[i] > this.mYMax) { this.mYMax = Y[i]; }
            }
            this.mXInvRange = 1 / (this.mXMax - this.mXMin);
            this.mYInvRange = 1 / (this.mYMax - this.mYMin);
            for (i = 0; i < numPoints; ++i) {
                this.mX[i] = (X[i] - this.mXMin) * this.mXInvRange;
                this.mY[i] = (Y[i] - this.mYMin) * this.mYInvRange;
            }
        } else {
            // The classical thin-plate spline uses the data as is. The
            // values mXMax and mYMax are not used, but they are initialized
            // anyway (to irrelevant numbers).
            this.mXMin = 0;
            this.mXMax = 1;
            this.mXInvRange = 1;
            this.mYMin = 0;
            this.mYMax = 1;
            this.mYInvRange = 1;
            for (i = 0; i < numPoints; ++i) {
                this.mX[i] = X[i];
                this.mY[i] = Y[i];
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
                    const t = Math.sqrt(dx * dx + dy * dy);
                    AMat.set(row, col, intpThinPlateSpline2Kernel(t));
                }
            }
        }

        // Compute matrix B [Nx3 matrix].
        const BMat = new GMatrix(numPoints, 3);
        for (row = 0; row < numPoints; ++row) {
            BMat.set(row, 0, 1);
            BMat.set(row, 1, this.mX[row]);
            BMat.set(row, 2, this.mY[row]);
        }

        // Compute A^{-1}.
        const invA = inverse(AMat);
        if (!invA.invertible) {
            return;
        }
        const invAMat = invA.inverse;

        // Compute P = B^T A^{-1} [3xN matrix].
        const PMat = multiplyATB(BMat, invAMat);

        // Compute Q = P B = B^T A^{-1} B [3x3 matrix].
        const QMat = multiplyAB(PMat, BMat);

        // Compute Q^{-1}.
        const invQ = inverse(QMat);
        if (!invQ.invertible) {
            return;
        }
        const invQMat = invQ.inverse;

        // Compute P*z.
        const prod = [0, 0, 0];
        for (row = 0; row < 3; ++row) {
            prod[row] = 0;
            for (i = 0; i < numPoints; ++i) {
                prod[row] += PMat.get(row, i) * F[i];
            }
        }

        // Compute 'b' vector for smooth thin plate spline.
        for (row = 0; row < 3; ++row) {
            this.mB[row] = 0;
            for (i = 0; i < 3; ++i) {
                this.mB[row] += invQMat.get(row, i) * prod[i];
            }
        }

        // Compute z - B*b.
        const tmp = new Array<number>(numPoints).fill(0);
        for (row = 0; row < numPoints; ++row) {
            tmp[row] = F[row];
            for (i = 0; i < 3; ++i) {
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
    // evaluate(x, y) will work properly.
    isInitialized(): boolean {
        return this.mInitialized;
    }

    // Evaluate the interpolator. If isInitialized() returns false, this
    // returns Number.MAX_VALUE.
    evaluate(x: number, y: number): number {
        if (this.mInitialized) {
            // Map (x,y) to the unit square.
            x = (x - this.mXMin) * this.mXInvRange;
            y = (y - this.mYMin) * this.mYInvRange;

            let result = this.mB[0] + this.mB[1] * x + this.mB[2] * y;
            for (let i = 0; i < this.mNumPoints; ++i) {
                const dx = x - this.mX[i];
                const dy = y - this.mY[i];
                const t = Math.sqrt(dx * dx + dy * dy);
                result += this.mA[i] * intpThinPlateSpline2Kernel(t);
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
                    const t = Math.sqrt(dx * dx + dy * dy);
                    functional += intpThinPlateSpline2Kernel(t) *
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
