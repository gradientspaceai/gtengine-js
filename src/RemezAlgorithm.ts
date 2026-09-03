// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RemezAlgorithm.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Remez exchange algorithm computes the minimax polynomial approximation
// P(x) of degree 'degree' to a function F(x) on the interval [xMin,xMax].
// The algorithm maintains degree+2 x-nodes. On each pass it fits the
// interpolating conditions F(x[i]) - P(x[i]) = (-1)^i * E, verifies that the
// error function E(x) = F(x) - P(x) oscillates in sign at the nodes,
// partitions [xMin,xMax] at the roots of E(x) and replaces the interior nodes
// by the local extrema of E(x) on those subintervals. This is the machinery
// used to generate the coefficient tables of the estimate-family functions.

import { logAssert } from './Logger.js';
import { GTE_C_HALF_PI } from './Constants.js';
import { Polynomial1 } from './Polynomial1.js';

// The port of 'std::function<T(T const&)>'; follows the Minimize1Function
// precedent of Minimize1.ts.
export type RemezFunction = (x: number) => number;

// Upstream execute(...) returns std::numeric_limits<size_t>::max() when the
// errors at the x-nodes stop oscillating, which terminates the algorithm. The
// port uses Number.MAX_SAFE_INTEGER as that sentinel.
export const REMEZ_FAILURE = Number.MAX_SAFE_INTEGER;

export class RemezAlgorithm {
    // Inputs to execute(...).
    private mF: RemezFunction;
    private mFDer: RemezFunction;
    private mXMin: number;
    private mXMax: number;
    private mDegree: number;
    private mMaxRemezIterations: number;
    private mMaxBisectionIterations: number;
    private mMaxBracketIterations: number;

    // Outputs from execute(...).
    private mPCoefficients: number[];
    private mEstimatedMaxError: number;
    private mXNodes: number[];
    private mErrors: number[];

    // Members used in the intermediate computations.
    private mFValues: number[];
    private mUCoefficients: number[];
    private mVCoefficients: number[];
    private mPartition: number[];

    constructor() {
        this.mF = (x: number) => x;
        this.mFDer = (x: number) => x;
        this.mXMin = 0;
        this.mXMax = 0;
        this.mDegree = 0;
        this.mMaxRemezIterations = 0;
        this.mMaxBisectionIterations = 0;
        this.mMaxBracketIterations = 0;
        this.mPCoefficients = [];
        this.mEstimatedMaxError = 0;
        this.mXNodes = [];
        this.mErrors = [];
        this.mFValues = [];
        this.mUCoefficients = [];
        this.mVCoefficients = [];
        this.mPartition = [];
    }

    // Compute the minimax polynomial of the specified degree. The return
    // value is the number of Remez iterations executed, or REMEZ_FAILURE
    // when the errors at the x-nodes stopped oscillating.
    execute(F: RemezFunction, FDer: RemezFunction, xMin: number, xMax: number,
        degree: number, maxRemezIterations: number, maxBisectionIterations: number,
        maxBracketIterations: number): number {
        logAssert(
            xMin < xMax &&
            degree > 0 &&
            maxRemezIterations > 0 &&
            maxBisectionIterations > 0 &&
            maxBracketIterations > 0,
            'Invalid input.');

        this.mF = F;
        this.mFDer = FDer;
        this.mXMin = xMin;
        this.mXMax = xMax;
        this.mDegree = degree;
        this.mMaxRemezIterations = maxRemezIterations;
        this.mMaxBisectionIterations = maxBisectionIterations;

        // NOTE (upstream quirk, preserved): mMaxBracketIterations is stored
        // by Execute(...) but never read anywhere in RemezAlgorithm.h. It is
        // kept here so the port's signature and validation match upstream.
        this.mMaxBracketIterations = maxBracketIterations;

        this.mPCoefficients = new Array<number>(this.mDegree + 1).fill(0);
        this.mEstimatedMaxError = 0;
        this.mXNodes = new Array<number>(this.mDegree + 2).fill(0);
        this.mErrors = new Array<number>(this.mDegree + 2).fill(0);

        this.mFValues = new Array<number>(this.mDegree + 2).fill(0);
        this.mUCoefficients = new Array<number>(this.mDegree + 1).fill(0);
        this.mVCoefficients = new Array<number>(this.mDegree + 1).fill(0);
        this.mEstimatedMaxError = 0;
        this.mPartition = new Array<number>(this.mDegree + 3).fill(0);

        this.computeInitialXNodes();
        let iteration = 0;
        for (iteration = 0; iteration < this.mMaxRemezIterations; ++iteration) {
            this.computeFAtXNodes();
            this.computeUCoefficients();
            this.computeVCoefficients();
            this.computeEstimatedError();
            this.computePCoefficients();
            if (this.isOscillatory()) {
                this.computePartition();
                this.computeXExtremes();
            }
            else {
                iteration = REMEZ_FAILURE;
                break;
            }
        }
        return iteration;
    }

    // The output of the algorithm. The coefficients are those of
    // P(x) = sum_{i=0}^{degree} c[i] * x^i.
    getCoefficients(): number[] {
        return this.mPCoefficients;
    }

    getEstimatedMaxError(): number {
        return this.mEstimatedMaxError;
    }

    getXNodes(): number[] {
        return this.mXNodes;
    }

    getErrors(): number[] {
        return this.mErrors;
    }

    private computeInitialXNodes(): void {
        // Get the Chebyshev nodes for the interval [-1,1].
        const numNodes = this.mXNodes.length;
        const halfPiDivDegree = GTE_C_HALF_PI / this.mDegree;
        const cosAngles = new Array<number>(numNodes).fill(0);
        cosAngles[0] = -1;
        for (let i = 1, j = 2 * this.mDegree - 1; i <= this.mDegree; ++i, j -= 2) {
            const angle = j * halfPiDivDegree;
            cosAngles[i] = Math.cos(angle);
        }
        cosAngles[numNodes - 1] = 1;
        if ((numNodes & 1) !== 0) {
            // Avoid the rounding errors when the angle is pi/2, where
            // cos(pi/2) is theoretically zero.
            cosAngles[Math.floor(numNodes / 2)] = 0;
        }

        // Transform the nodes to the interval [xMin, xMax].
        const half = 0.5;
        const center = half * (this.mXMax + this.mXMin);
        const radius = half * (this.mXMax - this.mXMin);
        this.mXNodes[0] = this.mXMin;
        for (let i = 1; i <= this.mDegree; ++i) {
            this.mXNodes[i] = center + radius * cosAngles[i];
        }
        this.mXNodes[this.mXNodes.length - 1] = this.mXMax;
    }

    private computeFAtXNodes(): void {
        for (let i = 0; i < this.mXNodes.length; ++i) {
            this.mFValues[i] = this.mF(this.mXNodes[i]);
        }
    }

    // Compute polynomial u(x) for which u(x[i]) = F(x[i]).
    private computeUCoefficients(): void {
        for (let i = 0; i < this.mUCoefficients.length; ++i) {
            this.mUCoefficients[i] = this.mFValues[i];
            for (let j = 0; j < i; ++j) {
                this.mUCoefficients[i] -= this.mUCoefficients[j];
                this.mUCoefficients[i] /= this.mXNodes[i] - this.mXNodes[j];
            }
        }
    }

    // Compute polynomial v(x) for which v(x[i]) = (-1)^i.
    private computeVCoefficients(): void {
        let sign = 1;
        for (let i = 0; i < this.mVCoefficients.length; ++i) {
            this.mVCoefficients[i] = sign;
            for (let j = 0; j < i; ++j) {
                this.mVCoefficients[i] -= this.mVCoefficients[j];
                this.mVCoefficients[i] /= this.mXNodes[i] - this.mXNodes[j];
            }
            sign = -sign;
        }
    }

    private computeEstimatedError(): void {
        const powNegOne = ((this.mDegree & 1) !== 0 ? -1 : 1);
        const xBack = this.mXNodes[this.mXNodes.length - 1];
        const fBack = this.mFValues[this.mFValues.length - 1];
        const uBack = this.evaluateU(xBack);
        const vBack = this.evaluateV(xBack);
        this.mEstimatedMaxError = (uBack - fBack) / (vBack + powNegOne);
    }

    private computePCoefficients(): void {
        // Compute the P-polynomial symbolically as a Newton polynomial in
        // order to obtain the coefficients from the t-powers.
        const constant = new Array<number>(this.mUCoefficients.length).fill(0);
        for (let i = 0; i < this.mUCoefficients.length; ++i) {
            constant[i] = this.mUCoefficients[i] -
                this.mEstimatedMaxError * this.mVCoefficients[i];
        }

        let index = this.mUCoefficients.length - 1;
        let poly = Polynomial1.fromCoefficients([constant[index]]);
        --index;
        for (let i = 1; i < this.mUCoefficients.length; ++i, --index) {
            const linear = Polynomial1.fromCoefficients([-this.mXNodes[index], 1]);
            poly = linear.mul(poly).add(constant[index]);
        }

        for (let i = 0; i < this.mPCoefficients.length; ++i) {
            // Polynomial1 eliminates leading zero coefficients, so a degree
            // drop leaves the high-order coefficients absent; they are zero.
            const c = poly.get(i);
            this.mPCoefficients[i] = (c !== undefined ? c : 0);
        }
    }

    private isOscillatory(): boolean {
        // Compute the errors F(x)-P(x) for the current nodes and verify they
        // are oscillatory.
        for (let i = 0; i < this.mXNodes.length; ++i) {
            this.mErrors[i] = this.mF(this.mXNodes[i]) - this.evaluateP(this.mXNodes[i]);
        }

        for (let i0 = 0, i1 = 1; i1 < this.mXNodes.length; i0 = i1++) {
            if ((this.mErrors[i0] > 0 && this.mErrors[i1] > 0) ||
                (this.mErrors[i0] < 0 && this.mErrors[i1] < 0)) {
                // The process terminates when the errors are not oscillatory.
                return false;
            }
        }
        return true;
    }

    private computePartition(): void {
        // Define E(x) = F(x) - P(x). Use bisection to compute the roots of
        // E(x). The algorithm partitions [xMin, xMax] into degree+2
        // subintervals, each subinterval with E(x) positive or with E(x)
        // negative. Later, the local extrema on the subintervals are computed
        // using a line-search algorithm. The extreme locations become the next
        // set of x-nodes.
        const half = 0.5;
        this.mPartition[0] = this.mXMin;
        this.mPartition[this.mPartition.length - 1] = this.mXMax;
        for (let i0 = 0, i1 = 1; i1 < this.mXNodes.length; i0 = i1++) {
            let x0 = this.mXNodes[i0];
            let x1 = this.mXNodes[i1];
            let xMid = 0;
            let eMid = 0;
            const sign0 = (this.mErrors[i0] > 0 ? 1 : -1);
            const sign1 = (this.mErrors[i1] > 0 ? 1 : -1);
            let signMid = 0;

            for (let iteration = 0; iteration < this.mMaxBisectionIterations; ++iteration) {
                xMid = half * (x0 + x1);
                if (xMid === x0 || xMid === x1) {
                    // We are at the limit of floating-point precision for the
                    // average of endpoints.
                    break;
                }

                // Update the correct endpoint to the midpoint.
                eMid = this.mF(xMid) - this.evaluateP(xMid);
                signMid = (eMid > 0 ? 1 : (eMid < 0 ? -1 : 0));
                if (signMid === sign0) {
                    x0 = xMid;
                }
                else if (signMid === sign1) {
                    x1 = xMid;
                }
                else {
                    // Found a root (numerically rounded to zero).
                    break;
                }
            }

            // It is possible that the maximum number of bisections was applied
            // without convergence. Use the last computed xMid as the root.
            this.mPartition[i1] = xMid;
        }
    }

    // Find the local extrema of E(x) on the subintervals of the partition.
    private computeXExtremes(): void {
        const nextXNodes = new Array<number>(this.mXNodes.length).fill(0);
        nextXNodes[0] = this.mXMin;
        for (let i0 = 1, i1 = 2; i0 < this.mDegree + 1; i0 = i1++) {
            nextXNodes[i0] = this.getXExtreme(this.mPartition[i0], this.mPartition[i1]);
        }
        nextXNodes[this.mDegree + 1] = this.mXMax;
        this.mXNodes = nextXNodes;
    }

    // Bisect [x0,x1] for the root of E'(x) = F'(x) - P'(x), which is the
    // location of the local extreme of E(x) on that subinterval.
    private getXExtreme(x0: number, x1: number): number {
        const half = 0.5;
        const eder0 = this.mFDer(x0) - this.evaluatePDer(x0);
        const eder1 = this.mFDer(x1) - this.evaluatePDer(x1);
        const signEDer0 = (eder0 > 0 ? 1 : (eder0 < 0 ? -1 : 0));
        const signEDer1 = (eder1 > 0 ? 1 : (eder1 < 0 ? -1 : 0));
        logAssert(
            signEDer0 * signEDer1 === -1,
            'The interval [x0,x1] does not bound a root.');

        let xmid = 0;
        let ederMid = 0;
        let signEMid = 0;
        for (let i = 0; i < this.mMaxBisectionIterations; ++i) {
            xmid = half * (x0 + x1);
            if (xmid === x0 || xmid === x1) {
                return xmid;
            }

            ederMid = this.mFDer(xmid) - this.evaluatePDer(xmid);
            signEMid = (ederMid > 0 ? 1 : (ederMid < 0 ? -1 : 0));
            if (signEMid === signEDer0) {
                x0 = xmid;
            }
            else if (signEMid === signEDer1) {
                x1 = xmid;
            }
            else {
                break;
            }
        }
        return xmid;
    }

    // Evaluate u(x) =
    //   u[0]+(x-xn[0])*(u[1]+(x-xn[1])*(u[2]+...+(x-xn[n-1])*u[n-1]))
    private evaluateU(x: number): number {
        let index = this.mUCoefficients.length - 1;
        let result = this.mUCoefficients[index];
        --index;
        for (let i = 1; i < this.mUCoefficients.length; ++i, --index) {
            result = this.mUCoefficients[index] + (x - this.mXNodes[index]) * result;
        }
        return result;
    }

    // Evaluate v(x) =
    //   v[0]+(x-xn[0])*(v[1]+(x-xn[1])*(v[2]+...+(x-xn[n-1])*v[n-1]))
    private evaluateV(x: number): number {
        let index = this.mVCoefficients.length - 1;
        let result = this.mVCoefficients[index];
        --index;
        for (let i = 1; i < this.mVCoefficients.length; ++i, --index) {
            result = this.mVCoefficients[index] + (x - this.mXNodes[index]) * result;
        }
        return result;
    }

    // Evaluate p(x) = sum_{i=0}^{n} p[i] * x^i.
    private evaluateP(x: number): number {
        let index = this.mPCoefficients.length - 1;
        let result = this.mPCoefficients[index];
        --index;
        for (let i = 1; i < this.mPCoefficients.length; ++i, --index) {
            result = this.mPCoefficients[index] + x * result;
        }
        return result;
    }

    private evaluatePDer(x: number): number {
        let index = this.mPCoefficients.length - 1;
        let result = index * this.mPCoefficients[index];
        --index;
        for (let i = 2; i < this.mPCoefficients.length; ++i, --index) {
            result = index * this.mPCoefficients[index] + x * result;
        }
        return result;
    }
}
