// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) TCBSplineCurve.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the tension-continuity-bias (TCB) spline for a set of key frames.
// The algorithm was invented by Kochanek and Bartels and is described in
// https://www.geometrictools.com/Documentation/KBSplines.pdf
//
// Port notes (following BSplineCurve): upstream
// 'template <int32_t N, typename T>' becomes a runtime dimension passed as
// the first constructor argument and carried by the key-frame points (Vector
// objects). The C++ null pointers 'firstOutTangent' and 'lastInTangent'
// become optional arguments. All input arrays are copied (C++ value
// semantics).
//
// Deviation from upstream: the upstream constructor never sets
// mConstructed = true, so 'operator bool' reports failure for a successfully
// built curve. The port sets it, so isConstructed() is meaningful. This does
// not affect Evaluate, which (unlike the other curve classes) does not test
// mConstructed.

import { logAssert } from './Logger.js';
import { ParametricCurve } from './ParametricCurve.js';
import { Vector, div, length as vectorLength } from './Vector.js';

export class TCBSplineCurve extends ParametricCurve {
    // The constructor inputs.
    protected mPoint: Vector[];
    protected mTension: number[];
    protected mContinuity: number[];
    protected mBias: number[];
    protected mLambda: number[];

    // Tangent vectors derived from the constructor inputs.
    protected mInTangent: Vector[];
    protected mOutTangent: Vector[];

    // Polynomial coefficients. The mA[] are the degree 0 coefficients, the
    // mB[] are the degree 1 coefficients, the mC[] are the degree 2
    // coefficients and the mD[] are the degree 3 coefficients.
    protected mA: Vector[];
    protected mB: Vector[];
    protected mC: Vector[];
    protected mD: Vector[];

    // The inputs point[], time[], tension[], continuity[] and bias[] must
    // have the same number of elements n >= 2. If you want the speed to be
    // continuous for the entire spline, the input lambda[] must have n
    // elements that are all positive; otherwise lambda[] should have 0
    // elements. If you want to specify the outgoing tangent at time[0] and
    // the incoming tangent at time[n-1], pass those vectors; otherwise, the
    // boundary tangents are computed by internally duplicating the boundary
    // points, which effectively means point[-1] = point[0] and
    // point[n] = point[n-1].
    constructor(dimension: number, point: readonly Vector[],
        time: readonly number[], tension: readonly number[],
        continuity: readonly number[], bias: readonly number[],
        lambda: readonly number[], firstOutTangent?: Vector,
        lastInTangent?: Vector) {
        const numSegments = (point.length >= 2 ? point.length - 1 : 0);
        super(dimension, numSegments,
            // Upstream indexes time.data() unconditionally; the port
            // substitutes zeros when 'time' is too short so that the
            // LogAssert below (not an out-of-range read) reports the error.
            time.length >= numSegments + 1
                ? time : new Array<number>(numSegments + 1).fill(0));

        logAssert(
            point.length >= 2 &&
            time.length === point.length &&
            tension.length === point.length &&
            continuity.length === point.length &&
            bias.length === point.length &&
            (lambda.length === 0 || lambda.length === point.length),
            'Invalid size in TCBSpline constructor.');

        this.mPoint = point.map(p => p.clone());
        this.mTension = tension.slice();
        this.mContinuity = continuity.slice();
        this.mBias = bias.slice();
        this.mLambda = lambda.slice();

        const n = point.length;
        this.mInTangent = new Array<Vector>(n);
        this.mOutTangent = new Array<Vector>(n);
        for (let i = 0; i < n; ++i) {
            this.mInTangent[i] = new Vector(dimension);
            this.mOutTangent[i] = new Vector(dimension);
        }

        const numSeg = this.getNumSegments();
        this.mA = new Array<Vector>(numSeg);
        this.mB = new Array<Vector>(numSeg);
        this.mC = new Array<Vector>(numSeg);
        this.mD = new Array<Vector>(numSeg);
        for (let i = 0; i < numSeg; ++i) {
            this.mA[i] = new Vector(dimension);
            this.mB[i] = new Vector(dimension);
            this.mC[i] = new Vector(dimension);
            this.mD[i] = new Vector(dimension);
        }

        this.computeFirstTangents(firstOutTangent);
        this.computeInteriorTangents();
        this.computeLastTangents(lastInTangent);
        this.computeCoefficients();

        this.mConstructed = true;
    }

    // Member access.
    getNumKeyFrames(): number {
        return this.mPoint.length;
    }

    getPoints(): readonly Vector[] {
        return this.mPoint;
    }

    getTensions(): readonly number[] {
        return this.mTension;
    }

    getContinuities(): readonly number[] {
        return this.mContinuity;
    }

    getBiases(): readonly number[] {
        return this.mBias;
    }

    getLambdas(): readonly number[] {
        return this.mLambda;
    }

    getInTangents(): readonly Vector[] {
        return this.mInTangent;
    }

    getOutTangents(): readonly Vector[] {
        return this.mOutTangent;
    }

    // Evaluation of the curve. It is required that order <= 3, which allows
    // computing derivatives through order 3. If you want only the position,
    // pass in order of 0. If you want the position and first derivative, pass
    // in order of 1, and so on. The output array 'jet' must have enough
    // storage to support the specified order (use createJet()). The values
    // are ordered as: position, first derivative, second derivative, and so
    // on.
    override evaluate(t: number, order: number, jet: Vector[]): void {
        const { key, u } = this.getKeyInfo(t);
        const n = this.mDimension;
        const A = this.mA[key].values;
        const B = this.mB[key].values;
        const C = this.mC[key].values;
        const D = this.mD[key].values;

        // Compute the position.
        const jet0 = new Vector(n);
        for (let k = 0; k < n; ++k) {
            jet0.values[k] = A[k] + u * (B[k] + u * (C[k] + u * D[k]));
        }
        jet[0] = jet0;

        if (order >= 1) {
            // Compute the first-order derivative. Upstream divides the vector
            // with 'operator/=', which multiplies by the reciprocal of the
            // scalar and produces the ZERO vector when the scalar is zero (see
            // Vector.h); the port's div() has the same semantics. A plain
            // componentwise division would differ in the last ulp and would
            // produce infinities for a zero-length segment.
            const delta = this.mTime[key + 1] - this.mTime[key];
            const jet1 = new Vector(n);
            for (let k = 0; k < n; ++k) {
                jet1.values[k] = B[k] + u * (2 * C[k] + (3 * u) * D[k]);
            }
            jet[1] = div(jet1, delta);

            if (order >= 2) {
                // Compute the second-order derivative.
                const deltaSqr = delta * delta;
                const jet2 = new Vector(n);
                for (let k = 0; k < n; ++k) {
                    jet2.values[k] = 2 * C[k] + (6 * u) * D[k];
                }
                jet[2] = div(jet2, deltaSqr);

                if (order === 3) {
                    const deltaCub = deltaSqr * delta;
                    const jet3 = new Vector(n);
                    for (let k = 0; k < n; ++k) {
                        jet3.values[k] = 6 * D[k];
                    }
                    jet[3] = div(jet3, deltaCub);
                }
            }
        }
    }

    // Support for construction.
    protected computeFirstTangents(firstOutTangent?: Vector): void {
        const n = this.mDimension;
        const out = this.mOutTangent[0].values;
        if (firstOutTangent !== undefined) {
            for (let k = 0; k < n; ++k) {
                out[k] = firstOutTangent.values[k];
            }
        }
        else {
            const omT = 1 - this.mTension[0];
            const omC = 1 - this.mContinuity[0];
            const omB = 1 - this.mBias[0];
            const twoDelta = 2 * (this.mTime[1] - this.mTime[0]);
            const coeff = omT * omC * omB / twoDelta;
            for (let k = 0; k < n; ++k) {
                out[k] = coeff
                    * (this.mPoint[1].values[k] - this.mPoint[0].values[k]);
            }
        }

        if (this.mLambda.length > 0) {
            for (let k = 0; k < n; ++k) {
                out[k] *= this.mLambda[0];
            }
        }

        this.mInTangent[0] = this.mOutTangent[0].clone();
    }

    protected computeLastTangents(lastInTangent?: Vector): void {
        const n = this.mDimension;
        const nm1 = this.mPoint.length - 1;
        const inTan = this.mInTangent[nm1].values;
        if (lastInTangent !== undefined) {
            for (let k = 0; k < n; ++k) {
                inTan[k] = lastInTangent.values[k];
            }
        }
        else {
            const nm2 = nm1 - 1;
            const omT = 1 - this.mTension[nm1];
            const omC = 1 - this.mContinuity[nm1];
            const opB = 1 + this.mBias[nm1];
            const twoDelta = 2 * (this.mTime[nm1] - this.mTime[nm2]);
            const coeff = omT * omC * opB / twoDelta;
            for (let k = 0; k < n; ++k) {
                inTan[k] = coeff
                    * (this.mPoint[nm1].values[k] - this.mPoint[nm2].values[k]);
            }
        }

        if (this.mLambda.length > 0) {
            for (let k = 0; k < n; ++k) {
                inTan[k] *= this.mLambda[nm1];
            }
        }

        this.mOutTangent[nm1] = this.mInTangent[nm1].clone();
    }

    protected computeInteriorTangents(): void {
        const dim = this.mDimension;
        const n = this.mPoint.length;
        for (let km1 = 0, k = 1, kp1 = 2; kp1 < n; km1 = k, k = kp1, ++kp1) {
            const P0 = this.mPoint[km1].values;
            const P1 = this.mPoint[k].values;
            const P2 = this.mPoint[kp1].values;
            const omT = 1 - this.mTension[k];
            const omC = 1 - this.mContinuity[k];
            const opC = 1 + this.mContinuity[k];
            const omB = 1 - this.mBias[k];
            const opB = 1 + this.mBias[k];
            const twoDelta0 = 2 * (this.mTime[k] - this.mTime[km1]);
            const twoDelta1 = 2 * (this.mTime[kp1] - this.mTime[k]);
            const inCoeff0 = omT * omC * opB / twoDelta0;
            const inCoeff1 = omT * opC * omB / twoDelta1;
            const outCoeff0 = omT * opC * opB / twoDelta0;
            const outCoeff1 = omT * omC * omB / twoDelta1;
            const inTan = this.mInTangent[k].values;
            const outTan = this.mOutTangent[k].values;
            for (let j = 0; j < dim; ++j) {
                const P1mP0 = P1[j] - P0[j];
                const P2mP1 = P2[j] - P1[j];
                inTan[j] = inCoeff0 * P1mP0 + inCoeff1 * P2mP1;
                outTan[j] = outCoeff0 * P1mP0 + outCoeff1 * P2mP1;
            }
        }

        if (this.mLambda.length > 0) {
            for (let k = 1, kp1 = 2; kp1 < n; k = kp1, ++kp1) {
                const inLength = vectorLength(this.mInTangent[k]);
                const outLength = vectorLength(this.mOutTangent[k]);
                const common = 2 * this.mLambda[k] / (inLength + outLength);
                const inCoeff = outLength * common;
                const outCoeff = inLength * common;
                const inTan = this.mInTangent[k].values;
                const outTan = this.mOutTangent[k].values;
                for (let j = 0; j < dim; ++j) {
                    inTan[j] *= inCoeff;
                    outTan[j] *= outCoeff;
                }
            }
        }
    }

    protected computeCoefficients(): void {
        const dim = this.mDimension;
        for (let k = 0, kp1 = 1; kp1 < this.mPoint.length; k = kp1, ++kp1) {
            const P0 = this.mPoint[k].values;
            const P1 = this.mPoint[kp1].values;
            const TOut0 = this.mOutTangent[k].values;
            const TIn1 = this.mInTangent[kp1].values;
            const delta = this.mTime[kp1] - this.mTime[k];
            const A = this.mA[k].values;
            const B = this.mB[k].values;
            const C = this.mC[k].values;
            const D = this.mD[k].values;
            for (let j = 0; j < dim; ++j) {
                const P1mP0 = P1[j] - P0[j];
                A[j] = P0[j];
                B[j] = delta * TOut0[j];
                C[j] = 3 * P1mP0 - delta * (2 * TOut0[j] + TIn1[j]);
                D[j] = -2 * P1mP0 + delta * (TOut0[j] + TIn1[j]);
            }
        }
    }

    // Determine the index i for which time[i] <= t < time[i+1]. The returned
    // value u is in [0,1].
    protected getKeyInfo(t: number): { key: number, u: number } {
        const time = this.mTime;
        if (t <= time[0]) {
            return { key: 0, u: 0 };
        }

        const numSegments = this.mA.length;
        if (t < time[numSegments]) {
            for (let i = 0; i < numSegments; ++i) {
                if (t < time[i + 1]) {
                    return {
                        key: i,
                        u: (t - time[i]) / (time[i + 1] - time[i])
                    };
                }
            }
        }

        return { key: numSegments - 1, u: 1 };
    }
}
