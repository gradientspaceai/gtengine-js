// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BSplineSurfaceFit.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Least-squares fitting of a rectangular lattice of sample points by an open
// uniform B-spline surface. The algorithm implemented here is based on the
// document
//   https://www.geometrictools.com/Documentation/BSplineSurfaceLeastSquaresFit.pdf
//
// Port notes: upstream is specialized to Vector3<Real> samples and control
// points; the port uses Vector objects of dimension 3 (the port has no
// separate Vector3 class). The C++ 'Vector3<Real> const* sampleData' pointer
// becomes a readonly array in the same row-major order,
// sample[i0 + numSamples0*i1]; the control points use the same layout.
// Upstream's Array2<Real> ATMat[d] is row major (numControls[d] rows by
// numSamples[d] columns) and its flat buffer is passed to
// BandedMatrix::SolveSystem<true>; the port uses a flat number[] with the
// same row-major indexing.

import { logAssert } from './Logger';
import { BasisFunction, BasisFunctionInput, UniqueKnot } from './BasisFunction';
import { BandedMatrix } from './BandedMatrix';
import { Vector } from './Vector';

export class BSplineSurfaceFit {
    // Input sample information.
    private mNumSamples: number[];
    private mSampleData: readonly Vector[];

    // The fitted B-spline surface, open and with uniform knots.
    private mDegree: number[];
    private mNumControls: number[];
    private mControlData: Vector[];
    private mBasis: BasisFunction[];

    // Construction. The preconditions for calling the constructor are
    //   1 <= degree0 && degree0 + 1 < numControls0 <= numSamples0
    //   1 <= degree1 && degree1 + 1 < numControls1 <= numSamples1
    // The sample data must be in row-major order. The control data is also
    // stored in row-major order.
    constructor(degree0: number, numControls0: number, numSamples0: number,
        degree1: number, numControls1: number, numSamples1: number,
        sampleData: readonly Vector[]) {
        logAssert(1 <= degree0 && degree0 + 1 < numControls0, 'Invalid degree.');
        logAssert(numControls0 <= numSamples0, 'Invalid number of controls.');
        logAssert(1 <= degree1 && degree1 + 1 < numControls1, 'Invalid degree.');
        logAssert(numControls1 <= numSamples1, 'Invalid number of controls.');
        logAssert(sampleData.length >= numSamples0 * numSamples1, 'Invalid sample data.');

        this.mSampleData = sampleData;
        this.mDegree = [degree0, degree1];
        this.mNumSamples = [numSamples0, numSamples1];
        this.mNumControls = [numControls0, numControls1];
        this.mControlData = new Array<Vector>(numControls0 * numControls1);
        for (let i = 0; i < this.mControlData.length; ++i) {
            this.mControlData[i] = new Vector(3);
        }
        this.mBasis = [new BasisFunction(), new BasisFunction()];

        const tMultiplier = [0, 0];
        let dim: number;
        for (dim = 0; dim < 2; ++dim) {
            const input = new BasisFunctionInput();
            input.numControls = this.mNumControls[dim];
            input.degree = this.mDegree[dim];
            input.uniform = true;
            input.periodic = false;
            input.numUniqueKnots = this.mNumControls[dim] - this.mDegree[dim] + 1;
            input.uniqueKnots = new Array<UniqueKnot>(input.numUniqueKnots);
            for (let i = 0; i < input.numUniqueKnots; ++i) {
                input.uniqueKnots[i] = new UniqueKnot();
            }
            input.uniqueKnots[0].t = 0;
            input.uniqueKnots[0].multiplicity = this.mDegree[dim] + 1;
            const last = input.numUniqueKnots - 1;
            const factor = 1 / last;
            for (let i = 1; i < last; ++i) {
                input.uniqueKnots[i].t = factor * i;
                input.uniqueKnots[i].multiplicity = 1;
            }
            input.uniqueKnots[last].t = 1;
            input.uniqueKnots[last].multiplicity = this.mDegree[dim] + 1;
            this.mBasis[dim].create(input);

            tMultiplier[dim] = 1 / (this.mNumSamples[dim] - 1);
        }

        // Fit the data points with a B-spline surface using a least-squares
        // error metric. The problem is of the form A0^T*A0*Q*A1^T*A1 =
        // A0^T*P*A1, where A0^T*A0 and A1^T*A1 are banded matrices, P
        // contains the sample data, and Q is the unknown matrix of control
        // points.
        let t: number;
        let i0: number, i1: number, i2: number;

        // Construct the matrices A0^T*A0 and A1^T*A1.
        const ataMat = [
            new BandedMatrix(this.mNumControls[0], this.mDegree[0] + 1, this.mDegree[0] + 1),
            new BandedMatrix(this.mNumControls[1], this.mDegree[1] + 1, this.mDegree[1] + 1)
        ];

        for (dim = 0; dim < 2; ++dim) {
            for (i0 = 0; i0 < this.mNumControls[dim]; ++i0) {
                for (i1 = 0; i1 < i0; ++i1) {
                    ataMat[dim].set(i0, i1, ataMat[dim].get(i1, i0));
                }

                let i1Max = i0 + this.mDegree[dim];
                if (i1Max >= this.mNumControls[dim]) {
                    i1Max = this.mNumControls[dim] - 1;
                }

                for (i1 = i0; i1 <= i1Max; ++i1) {
                    let value = 0;
                    for (i2 = 0; i2 < this.mNumSamples[dim]; ++i2) {
                        t = tMultiplier[dim] * i2;
                        const { minIndex: imin, maxIndex: imax } = this.mBasis[dim].evaluate(t, 0);
                        if (imin <= i0 && i0 <= imax && imin <= i1 && i1 <= imax) {
                            const b0 = this.mBasis[dim].getValue(0, i0);
                            const b1 = this.mBasis[dim].getValue(0, i1);
                            value += b0 * b1;
                        }
                    }
                    ataMat[dim].set(i0, i1, value);
                }
            }
        }

        // Construct the matrices A0^T and A1^T. A[d]^T has mNumControls[d]
        // rows and mNumSamples[d] columns, stored in row-major order.
        const atMat: number[][] = [[], []];
        for (dim = 0; dim < 2; ++dim) {
            atMat[dim] = new Array<number>(this.mNumControls[dim] * this.mNumSamples[dim]).fill(0);
            for (i0 = 0; i0 < this.mNumControls[dim]; ++i0) {
                for (i1 = 0; i1 < this.mNumSamples[dim]; ++i1) {
                    t = tMultiplier[dim] * i1;
                    const { minIndex: imin, maxIndex: imax } = this.mBasis[dim].evaluate(t, 0);
                    if (imin <= i0 && i0 <= imax) {
                        atMat[dim][i0 * this.mNumSamples[dim] + i1] = this.mBasis[dim].getValue(0, i0);
                    }
                }
            }
        }

        // Compute X0 = (A0^T*A0)^{-1}*A0^T and X1 = (A1^T*A1)^{-1}*A1^T by
        // solving the linear systems A0^T*A0*X0 = A0^T and A1^T*A1*X1 = A1^T.
        for (dim = 0; dim < 2; ++dim) {
            const solved = ataMat[dim].solveSystemMatrix(atMat[dim], this.mNumSamples[dim], true);
            logAssert(solved,
                'Failed to solve linear system in BSplineSurfaceFit constructor.');
        }

        // The control points for the fitted surface are stored in the matrix
        // Q = X0*P*X1^T, where P is the matrix of sample data.
        for (i1 = 0; i1 < this.mNumControls[1]; ++i1) {
            for (i0 = 0; i0 < this.mNumControls[0]; ++i0) {
                const sum = new Vector(3);
                for (let j1 = 0; j1 < this.mNumSamples[1]; ++j1) {
                    const x1Value = atMat[1][i1 * this.mNumSamples[1] + j1];
                    for (let j0 = 0; j0 < this.mNumSamples[0]; ++j0) {
                        const x0Value = atMat[0][i0 * this.mNumSamples[0] + j0];
                        const sample = this.mSampleData[j0 + this.mNumSamples[0] * j1];
                        const scalar = x0Value * x1Value;
                        for (let k = 0; k < 3; ++k) {
                            sum.values[k] += scalar * sample.values[k];
                        }
                    }
                }
                this.mControlData[i0 + this.mNumControls[0] * i1] = sum;
            }
        }
    }

    // Access to input sample information.
    getNumSamples(dimension: number): number {
        return this.mNumSamples[dimension];
    }

    getSampleData(): readonly Vector[] {
        return this.mSampleData;
    }

    // Access to output control point and surface information.
    getDegree(dimension: number): number {
        return this.mDegree[dimension];
    }

    getNumControls(dimension: number): number {
        return this.mNumControls[dimension];
    }

    // The returned array aliases internal storage (upstream returns a
    // pointer).
    getControlData(): Vector[] {
        return this.mControlData;
    }

    getBasis(dimension: number): BasisFunction {
        return this.mBasis[dimension];
    }

    // Evaluation of the B-spline surface. It is defined for 0 <= u <= 1 and
    // 0 <= v <= 1. If a parameter value is outside [0,1], it is clamped to
    // [0,1].
    getPosition(u: number, v: number): Vector {
        const ru = this.mBasis[0].evaluate(u, 0);
        const rv = this.mBasis[1].evaluate(v, 0);

        const position = new Vector(3);
        for (let iv = rv.minIndex; iv <= rv.maxIndex; ++iv) {
            const value1 = this.mBasis[1].getValue(0, iv);
            for (let iu = ru.minIndex; iu <= ru.maxIndex; ++iu) {
                const value0 = this.mBasis[0].getValue(0, iu);
                const control = this.mControlData[iu + this.mNumControls[0] * iv];
                const scalar = value0 * value1;
                for (let k = 0; k < 3; ++k) {
                    position.values[k] += scalar * control.values[k];
                }
            }
        }
        return position;
    }
}
