// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BSplineCurveFit.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Least-squares fitting of a sequence of sample points by an open uniform
// B-spline curve. The algorithm implemented here is based on the document
//   https://www.geometrictools.com/Documentation/BSplineCurveLeastSquaresFit.pdf
//
// Port notes: the C++ 'Real const* sampleData' pointer (a contiguous block of
// numSamples * dimension values) becomes a 'readonly number[]' with the same
// flat layout; the fitted control points are likewise a flat array of
// numControls * dimension values. Upstream's Array2<Real> ATMat is stored in
// row-major order (numControls rows by numSamples columns), and its flat
// buffer is passed to BandedMatrix::SolveSystem<true>; the port uses a flat
// number[] with the same row-major indexing, passed to
// BandedMatrix.solveSystemMatrix(..., true). The Evaluate output pointer
// becomes a returned array.

import { logAssert } from './Logger.js';
import { BasisFunction, BasisFunctionInput, UniqueKnot } from './BasisFunction.js';
import { BandedMatrix } from './BandedMatrix.js';

export class BSplineCurveFit {
    // Input sample information.
    private mDimension: number;
    private mNumSamples: number;
    private mSampleData: readonly number[];

    // The fitted B-spline curve, open and with uniform knots.
    private mDegree: number;
    private mNumControls: number;
    private mControlData: number[];
    private mBasis: BasisFunction;

    // Construction. The preconditions for calling the constructor are
    // 1 <= degree && degree < numControls <= numSamples - degree - 1.
    // The sample points are contiguous blocks of 'dimension' values stored in
    // sampleData.
    constructor(dimension: number, numSamples: number, sampleData: readonly number[],
        degree: number, numControls: number) {
        logAssert(dimension >= 1, 'Invalid dimension.');
        logAssert(1 <= degree && degree < numControls, 'Invalid degree.');
        logAssert(sampleData.length >= dimension * numSamples, 'Invalid sample data.');
        logAssert(numControls <= numSamples - degree - 1, 'Invalid number of controls.');

        this.mDimension = dimension;
        this.mNumSamples = numSamples;
        this.mSampleData = sampleData;
        this.mDegree = degree;
        this.mNumControls = numControls;
        this.mControlData = new Array<number>(dimension * numControls).fill(0);

        const input = new BasisFunctionInput();
        input.numControls = numControls;
        input.degree = degree;
        input.uniform = true;
        input.periodic = false;
        input.numUniqueKnots = numControls - degree + 1;
        input.uniqueKnots = new Array<UniqueKnot>(input.numUniqueKnots);
        for (let i = 0; i < input.numUniqueKnots; ++i) {
            input.uniqueKnots[i] = new UniqueKnot();
        }
        input.uniqueKnots[0].t = 0;
        input.uniqueKnots[0].multiplicity = degree + 1;
        const last = input.numUniqueKnots - 1;
        const factor = 1 / last;
        for (let i = 1; i < last; ++i) {
            input.uniqueKnots[i].t = factor * i;
            input.uniqueKnots[i].multiplicity = 1;
        }
        input.uniqueKnots[last].t = 1;
        input.uniqueKnots[last].multiplicity = degree + 1;
        this.mBasis = new BasisFunction();
        this.mBasis.create(input);

        // Fit the data points with a B-spline curve using a least-squares
        // error metric. The problem is of the form A^T*A*Q = A^T*P, where
        // A^T*A is a banded matrix, P contains the sample data, and Q is the
        // unknown vector of control points.
        const tMultiplier = 1 / (this.mNumSamples - 1);
        let t: number;
        let i0: number, i1: number, i2: number, j: number;

        // Construct the matrix A^T*A.
        const degp1 = this.mDegree + 1;
        const numBands = (this.mNumControls > degp1 ? degp1 : this.mDegree);
        const ataMat = new BandedMatrix(this.mNumControls, numBands, numBands);
        for (i0 = 0; i0 < this.mNumControls; ++i0) {
            for (i1 = 0; i1 < i0; ++i1) {
                ataMat.set(i0, i1, ataMat.get(i1, i0));
            }

            let i1Max = i0 + this.mDegree;
            if (i1Max >= this.mNumControls) {
                i1Max = this.mNumControls - 1;
            }

            for (i1 = i0; i1 <= i1Max; ++i1) {
                let value = 0;
                for (i2 = 0; i2 < this.mNumSamples; ++i2) {
                    t = tMultiplier * i2;
                    const { minIndex: imin, maxIndex: imax } = this.mBasis.evaluate(t, 0);
                    if (imin <= i0 && i0 <= imax && imin <= i1 && i1 <= imax) {
                        const b0 = this.mBasis.getValue(0, i0);
                        const b1 = this.mBasis.getValue(0, i1);
                        value += b0 * b1;
                    }
                }
                ataMat.set(i0, i1, value);
            }
        }

        // Construct the matrix A^T. It has mNumControls rows and mNumSamples
        // columns, stored in row-major order.
        const atMat = new Array<number>(this.mNumControls * this.mNumSamples).fill(0);
        for (i0 = 0; i0 < this.mNumControls; ++i0) {
            for (i1 = 0; i1 < this.mNumSamples; ++i1) {
                t = tMultiplier * i1;
                const { minIndex: imin, maxIndex: imax } = this.mBasis.evaluate(t, 0);
                if (imin <= i0 && i0 <= imax) {
                    atMat[i0 * this.mNumSamples + i1] = this.mBasis.getValue(0, i0);
                }
            }
        }

        // Compute X0 = (A^T*A)^{-1}*A^T by solving the linear system
        // A^T*A*X = A^T.
        const solved = ataMat.solveSystemMatrix(atMat, this.mNumSamples, true);
        logAssert(solved, 'Failed to solve linear system.');

        // The control points for the fitted curve are stored in the vector
        // Q = X0*P, where P is the vector of sample data.
        for (i0 = 0; i0 < this.mNumControls; ++i0) {
            const qBase = i0 * this.mDimension;
            for (i1 = 0; i1 < this.mNumSamples; ++i1) {
                const pBase = i1 * this.mDimension;
                const xValue = atMat[i0 * this.mNumSamples + i1];
                for (j = 0; j < this.mDimension; ++j) {
                    this.mControlData[qBase + j] += xValue * this.mSampleData[pBase + j];
                }
            }
        }

        // Set the first and last output control points to match the first and
        // last input samples. This supports the application of fitting
        // keyframe data with B-spline curves. The user expects that the curve
        // passes through the first and last positions in order to support
        // matching two consecutive keyframe sequences.
        const cEnd1 = this.mDimension * (this.mNumControls - 1);
        const sEnd1 = this.mDimension * (this.mNumSamples - 1);
        for (j = 0; j < this.mDimension; ++j) {
            this.mControlData[j] = this.mSampleData[j];
            this.mControlData[cEnd1 + j] = this.mSampleData[sEnd1 + j];
        }
    }

    // Access to input sample information.
    getDimension(): number {
        return this.mDimension;
    }

    getNumSamples(): number {
        return this.mNumSamples;
    }

    getSampleData(): readonly number[] {
        return this.mSampleData;
    }

    // Access to output control point and curve information.
    getDegree(): number {
        return this.mDegree;
    }

    getNumControls(): number {
        return this.mNumControls;
    }

    // The returned array aliases internal storage (upstream returns a
    // pointer).
    getControlData(): number[] {
        return this.mControlData;
    }

    getBasis(): BasisFunction {
        return this.mBasis;
    }

    // Evaluation of the B-spline curve. It is defined for 0 <= t <= 1. If a
    // t-value is outside [0,1], an open spline clamps it to [0,1]. The
    // returned array has 'dimension' elements.
    evaluate(t: number, order: number): number[] {
        const { minIndex: imin, maxIndex: imax } = this.mBasis.evaluate(t, order);

        const value = new Array<number>(this.mDimension);
        let source = this.mDimension * imin;
        let basisValue = this.mBasis.getValue(order, imin);
        for (let j = 0; j < this.mDimension; ++j) {
            value[j] = basisValue * this.mControlData[source++];
        }

        for (let i = imin + 1; i <= imax; ++i) {
            basisValue = this.mBasis.getValue(order, i);
            for (let j = 0; j < this.mDimension; ++j) {
                value[j] += basisValue * this.mControlData[source++];
            }
        }

        return value;
    }

    getPosition(t: number): number[] {
        return this.evaluate(t, 0);
    }
}
