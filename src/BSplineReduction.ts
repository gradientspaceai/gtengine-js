// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BSplineReduction.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The BSplineReduction class is an implementation of the algorithm in
//   https://www.geometrictools.com/Documentation/BSplineReduction.pdf
// for least-squares fitting of points in the continuous sense by an L2
// integral norm. The least-squares fitting implemented in the file
// BSplineCurveFit.h is in the discrete sense by an L2 summation. The intended
// use for this class is to take an open B-spline curve, defined by its
// control points and degree, and reduce the number of control points
// dramatically to obtain another curve that is close to the original one.
//
// The input control points must number 2 or more. The input degree must
// satisfy 1 <= degree < inControls.length. The degree of the output curve is
// the same as that of the input curve. The input fraction must be in [0,1].
// If the fraction is 1, the output curve is identical to the input curve. If
// the fraction is too small to produce a valid number of control points, the
// number of output control points is degree+1.
//
// Port notes: upstream 'template <int32_t N, typename Real>' becomes a
// runtime dimension carried by the control-point Vector objects. The
// 'operator()' with the 'outControls' output reference becomes compute(...)
// returning the array of control points, per the PORTING.md operator mapping.
// The C++ member state (mBasis, mIndex, ...) exists only to feed the
// std::function integrand; the port keeps the same members so that the
// integrand is the same function of the same state.
//
// Accuracy note (upstream behavior, preserved): the entries of A and B are
// integrals of products of B-spline basis functions, which are piecewise
// polynomials whose derivatives jump at the knots. Upstream evaluates them
// with Integration::Romberg(8, ...) over the whole overlap of the supports
// rather than knot span by knot span, so the entries carry a quadrature
// error on the order of 1e-3. In consequence the fit does not reproduce a
// curve that already lies in the output spline space exactly, and the rows
// of A^{-1}B sum to one only to that accuracy (so the result is only
// approximately invariant under translation of the control points).

import { BandedMatrix } from './BandedMatrix.js';
import { GMatrix } from './GMatrix.js';
import { Integration } from './Integration.js';
import { IntrIntervalsFI } from './IntrIntervals.js';
import { logAssert } from './Logger.js';
import { multiplyAB } from './Matrix.js';
import { Vector } from './Vector.js';

export class BSplineReduction {
    private mDegree: number;
    private mQuantity: [number, number];
    private mNumKnots: [number, number];  // N+D+2
    private mKnot: [number[], number[]];

    // For the integration-based least-squares fitting.
    private mBasis: [number, number];
    private mIndex: [number, number];

    constructor() {
        this.mDegree = 0;
        this.mQuantity = [0, 0];
        this.mNumKnots = [0, 0];
        this.mKnot = [[], []];
        this.mBasis = [0, 0];
        this.mIndex = [0, 0];
    }

    // Reduce the control points of a B-spline curve of the given degree. The
    // returned array holds the control points of the reduced curve; the knots
    // of both curves are those of an open uniform B-spline on [0,1].
    compute(inControls: readonly Vector[], degree: number,
        fraction: number): Vector[] {
        const numInControls = inControls.length;
        logAssert(numInControls >= 2 && 1 <= degree && degree < numInControls,
            'Invalid input.');

        const dimension = inControls[0].values.length;

        // Clamp the number of control points to [degree+1, quantity-1].
        let numOutControls = Math.trunc(fraction * numInControls);
        if (numOutControls >= numInControls) {
            return inControls.map(control => control.clone());
        }
        if (numOutControls < degree + 1) {
            numOutControls = degree + 1;
        }

        // Allocate output control points.
        const outControls = new Array<Vector>(numOutControls);

        // Set up basis function parameters. Function 0 corresponds to the
        // output curve. Function 1 corresponds to the input curve.
        this.mDegree = degree;
        this.mQuantity[0] = numOutControls;
        this.mQuantity[1] = numInControls;

        for (let j = 0; j <= 1; ++j) {
            this.mNumKnots[j] = this.mQuantity[j] + this.mDegree + 1;
            this.mKnot[j] = new Array<number>(this.mNumKnots[j]).fill(0);

            let i: number;
            for (i = 0; i <= this.mDegree; ++i) {
                this.mKnot[j][i] = 0;
            }

            const denom = this.mQuantity[j] - this.mDegree;
            const factor = 1 / denom;
            for (/**/; i < this.mQuantity[j]; ++i) {
                this.mKnot[j][i] = (i - this.mDegree) * factor;
            }

            for (/**/; i < this.mNumKnots[j]; ++i) {
                this.mKnot[j][i] = 1;
            }
        }

        // Construct matrix A (depends only on the output basis function).
        let value: number, tmin: number, tmax: number;
        let i0: number, i1: number;

        this.mBasis[0] = 0;
        this.mBasis[1] = 0;

        const integrand = (t: number): number => {
            const value0 = this.f(this.mBasis[0], this.mIndex[0], this.mDegree, t);
            const value1 = this.f(this.mBasis[1], this.mIndex[1], this.mDegree, t);
            return value0 * value1;
        };

        const A = new BandedMatrix(this.mQuantity[0], this.mDegree, this.mDegree);
        for (i0 = 0; i0 < this.mQuantity[0]; ++i0) {
            this.mIndex[0] = i0;
            tmax = this.maxSupport(0, i0);

            for (i1 = i0; i1 <= i0 + this.mDegree && i1 < this.mQuantity[0]; ++i1) {
                this.mIndex[1] = i1;
                tmin = this.minSupport(0, i1);

                value = Integration.romberg(8, tmin, tmax, integrand);
                A.set(i0, i1, value);
                A.set(i1, i0, value);
            }
        }

        // Construct A^{-1}. TODO (upstream): this is inefficient; use an
        // iterative scheme to invert A?
        const inverseElements = new Array<number>(this.mQuantity[0] * this.mQuantity[0]).fill(0);
        const invertible = A.computeInverse(inverseElements, true);
        logAssert(invertible, 'Failed to invert matrix.');
        const invA = new GMatrix(this.mQuantity[0], this.mQuantity[0]);
        for (i0 = 0; i0 < this.mQuantity[0]; ++i0) {
            for (i1 = 0; i1 < this.mQuantity[0]; ++i1) {
                invA.set(i0, i1, inverseElements[i0 * this.mQuantity[0] + i1]);
            }
        }

        // Construct B (depends on both input and output basis functions).
        this.mBasis[1] = 1;
        const B = new GMatrix(this.mQuantity[0], this.mQuantity[1]);
        const query = new IntrIntervalsFI();
        for (i0 = 0; i0 < this.mQuantity[0]; ++i0) {
            this.mIndex[0] = i0;
            const tmin0 = this.minSupport(0, i0);
            const tmax0 = this.maxSupport(0, i0);

            for (i1 = 0; i1 < this.mQuantity[1]; ++i1) {
                this.mIndex[1] = i1;
                const tmin1 = this.minSupport(1, i1);
                const tmax1 = this.maxSupport(1, i1);

                const result = query.find([tmin0, tmax0], [tmin1, tmax1]);
                if (result.numIntersections === 2) {
                    value = Integration.romberg(8, result.overlap[0],
                        result.overlap[1], integrand);

                    B.set(i0, i1, value);
                } else {
                    B.set(i0, i1, 0);
                }
            }
        }

        // Construct A^{-1}*B.
        const prod = multiplyAB(invA, B);

        // Construct the control points for the least-squares curve.
        for (i0 = 0; i0 < this.mQuantity[0]; ++i0) {
            outControls[i0] = new Vector(dimension);
        }
        for (i0 = 0; i0 < this.mQuantity[0]; ++i0) {
            for (i1 = 0; i1 < this.mQuantity[1]; ++i1) {
                const scalar = prod.get(i0, i1);
                const control = inControls[i1];
                for (let k = 0; k < dimension; ++k) {
                    outControls[i0].values[k] += control.values[k] * scalar;
                }
            }
        }

        return outControls;
    }

    private minSupport(basis: number, i: number): number {
        return this.mKnot[basis][i];
    }

    private maxSupport(basis: number, i: number): number {
        return this.mKnot[basis][i + 1 + this.mDegree];
    }

    // The Cox-de Boor recursion for the B-spline basis function of index i
    // and degree j for the knot vector selected by 'basis'.
    private f(basis: number, i: number, j: number, t: number): number {
        if (j > 0) {
            let result = 0;

            let denom = this.mKnot[basis][i + j] - this.mKnot[basis][i];
            if (denom > 0) {
                result += (t - this.mKnot[basis][i])
                    * this.f(basis, i, j - 1, t) / denom;
            }

            denom = this.mKnot[basis][i + j + 1] - this.mKnot[basis][i + 1];
            if (denom > 0) {
                result += (this.mKnot[basis][i + j + 1] - t)
                    * this.f(basis, i + 1, j - 1, t) / denom;
            }

            return result;
        }

        if (this.mKnot[basis][i] <= t && t < this.mKnot[basis][i + 1]) {
            return 1;
        }
        return 0;
    }
}
