// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpBSplineUniform.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// IntpBSplineUniform is the class for B-spline interpolation of uniformly
// spaced N-dimensional data. The algorithm is described in
//   https://www.geometrictools.com/Documentation/BSplineInterpolation.pdf
//
// The Controls adapter allows access to your control points without regard
// to how you organize your data. You can even defer the computation of a
// control point until it is needed via the get(...) call that Controls must
// provide, and you can cache the points according to your own needs.
//
// Port notes:
// - The C++ template parameter 'Controls' is duck-typed, so the port
//   declares the adapter interface explicitly as
//   IntpBSplineUniformControls<T>. The control-point type T is arbitrary;
//   C++ requires it to support assignment, addition and scalar
//   multiplication, and since TypeScript has no operator overloading the
//   adapter also supplies those two operations as add(c0, c1) and
//   mul(c, s). The zero element is still passed to the constructor as
//   'ctZero', matching upstream.
// - The two C++ control-point accessors, 'operator()(int32_t const* tuple)'
//   for the general-dimension code and 'operator()(i0, i1, ...)' for the
//   dimension specializations, are unified as get(indices), which takes an
//   index array of the appropriate length.
// - The C++ class template IntpBSplineUniform<Real, Controls, N> has a
//   general implementation (N known at compile time, and the N = 0
//   specialization for run-time N) plus hand-optimized specializations for
//   N = 1, 2 and 3. TypeScript has no compile-time dimension, so the
//   general implementation is the single class IntpBSplineUniform<T> (the
//   run-time-dimension form, including its guard on the sizes of the
//   'order' and 't' inputs), and the specializations are the separate
//   classes IntpBSplineUniform1/2/3<T>. Upstream the specializations derive
//   from IntpBSplineUniformShared only to reuse its static helpers and
//   therefore construct the base with a do-nothing default constructor; the
//   port makes them standalone classes that call the static helpers
//   directly, so IntpBSplineUniformShared<T> has only the real constructor.
// - The anonymous enum of cache modes becomes the exported enum
//   IntpBSplineUniformCacheMode.
// - ComputeBlendingMatrix, ComputeDCoefficients, ComputePowers and GetKey
//   write to reference parameters upstream; the port returns the values.
// - Upstream's dead local 'Polynomial1<Real> sm1' in ComputeBlendingMatrix
//   is dropped.

import { logError } from './Logger.js';
import { Polynomial1 } from './Polynomial1.js';

// Support for caching the intermediate tensor product of control points
// with the blending matrices. A precached container has all elements
// precomputed before any evaluate(...) calls; the flags are all set to
// true. A cached container fills the elements on demand; the flags are
// initially false, indicating the element has not yet been computed. After
// it is computed and stored, the flag is set to true.
export enum IntpBSplineUniformCacheMode {
    NO_CACHING = 0,
    PRE_CACHING = 1,
    ON_DEMAND_CACHING = 2
}

// The adapter that provides the control points and the arithmetic on the
// control-point type T.
export interface IntpBSplineUniformControls<T> {
    // The number of elements in the specified dimension.
    getSize(dimension: number): number;

    // Get a control point based on an n-tuple lookup. The interpolator does
    // not need to know your organization; all it needs is the desired
    // control point. The 'indices' input has one element per dimension.
    get(indices: readonly number[]): T;

    // C2 = C0 + C1.
    add(c0: T, c1: T): T;

    // C1 = C0 * s.
    mul(c0: T, s: number): T;
}

// The coefficient of t^i in p, where the coefficients of the terms of
// degree larger than the degree of p are zero. Upstream indexes the
// coefficient array directly, which reads out of bounds when an arithmetic
// operation has reduced the degree of the polynomial.
function coefficient(p: Polynomial1, i: number): number {
    return i <= p.getDegree() ? p.get(i) : 0;
}

export abstract class IntpBSplineUniformShared<T> {
    // Compute the blending matrix that combines the control points and the
    // polynomial vector. The matrix A is stored in row-major order.
    static computeBlendingMatrix(degree: number): number[] {
        const degreeP1 = degree + 1;
        const A = new Array<number>(degreeP1 * degreeP1).fill(0);

        if (degree === 0) {
            A[0] = 1;
            return A;
        }

        // P_{0,0}(s)
        const P: Polynomial1[] = [];
        for (let k = 0; k < degreeP1; ++k) {
            P.push(new Polynomial1(0));
        }
        P[0].set(0, 1);

        // L0 = s/j
        const L0 = new Polynomial1(1);
        L0.set(0, 0);

        // L1(s) = (j + 1 - s)/j
        const L1 = new Polynomial1(1);

        // Compute
        //   P_{j,k}(s) = L0(s)*P_{j-1,k}(s) + L1(s)*P_{j-1,k-1}(s-1)
        // for 0 <= k <= j where 1 <= j <= degree. When k = 0,
        // P_{j-1,-1}(s) = 0, so P_{j,0}(s) = L0(s)*P_{j-1,0}(s). When k = j,
        // P_{j-1,j}(s) = 0, so P_{j,j}(s) = L1(s)*P_{j-1,j-1}(s). The
        // polynomials at level j-1 are currently stored in P[0] through
        // P[j-1]. The polynomials at level j are computed and stored in P[0]
        // through P[j]; that is, they are computed in place to reduce memory
        // usage and copying. This requires computing P[k] (level j) from
        // P[k] (level j-1) and P[k-1] (level j-1), which means we have to
        // process k = j down to k = 0.
        for (let j = 1; j <= degree; ++j) {
            const invJ = 1 / j;
            L0.set(1, invJ);
            L1.set(0, 1 + invJ);
            L1.set(1, -invJ);

            for (let k = j; k >= 0; --k) {
                let result = Polynomial1.fromCoefficients([0]);

                if (k > 0) {
                    result = result.add(L1.mul(P[k - 1].getTranslation(1)));
                }

                if (k < j) {
                    result = result.add(L0.mul(P[k]));
                }

                P[k] = result;
            }
        }

        // Compute Q_{d,k}(s) = P_{d,k}(s + k).
        const Q: Polynomial1[] = [];
        for (let k = 0; k <= degree; ++k) {
            Q.push(P[k].getTranslation(-k));
        }

        // Extract the matrix A from the Q-polynomials. Row r of A contains
        // the coefficients of Q_{d,d-r}(s).
        for (let k = 0, row = degree; k <= degree; ++k, --row) {
            for (let col = 0; col <= degree; ++col) {
                A[col + degreeP1 * row] = coefficient(Q[k], col);
            }
        }

        return A;
    }

    // Compute the coefficients for the derivative polynomial terms.
    static computeDCoefficients(degree: number):
        { dCoefficients: number[], ellMax: number[] } {
        const numDCoefficients = (degree + 1) * (degree + 2) / 2;
        const dCoefficients = new Array<number>(numDCoefficients).fill(1);

        for (let order = 1, col0 = 0, col1 = degree + 1; order <= degree; ++order) {
            ++col0;
            for (let c = order, m = 1; c <= degree; ++c, ++m, ++col0, ++col1) {
                dCoefficients[col1] = dCoefficients[col0] * m;
            }
        }

        const ellMax = new Array<number>(degree + 1).fill(0);
        ellMax[0] = degree;
        for (let i0 = 0, i1 = 1; i1 <= degree; i0 = i1++) {
            ellMax[i1] = ellMax[i0] + degree - i0;
        }

        return { dCoefficients, ellMax };
    }

    // Compute powers of ds/dt.
    static computePowers(degree: number, numControls: number, tmin: number,
        tmax: number): number[] {
        const dsdt = (numControls - degree) / (tmax - tmin);
        // Upstream sizes the array to degree+1 and then writes index 1
        // unconditionally, which overruns the buffer when degree is 0. The
        // JavaScript array simply grows, and the resulting element is the
        // ds/dt value that GetKey expects, so degree 0 is well behaved here.
        const powerDSDT = new Array<number>(degree + 1).fill(0);
        powerDSDT[0] = 1;
        powerDSDT[1] = dsdt;
        for (let i = 2, im1 = 1; i <= degree; ++i, ++im1) {
            powerDSDT[i] = powerDSDT[im1] * dsdt;
        }
        return powerDSDT;
    }

    // Determine the interval [index,index+1) corresponding to the specified
    // value of t and compute u in that interval.
    static getKey(t: number, tmin: number, tmax: number, dsdt: number,
        numControls: number, degree: number): { index: number, u: number } {
        // Compute s - d = ((c + 1 - d)/(c + 1))(t + 1/2), the index for
        // which d + index <= s < d + index + 1. Let u = s - d - index so
        // that 0 <= u < 1.
        if (t > tmin) {
            if (t < tmax) {
                const smd = dsdt * (t - tmin);
                const index = Math.floor(smd);
                return { index, u: smd - index };
            }
            else {
                // In the evaluation, s = c + 1 - d and i = c - d. This
                // causes s-d-i to be 1 in G_c(c+1-d). Effectively, the
                // selection of i extends the s-domain [d,c+1) to its
                // support [d,c+1].
                return { index: numControls - 1 - degree, u: 1 };
            }
        }
        else {
            return { index: 0, u: 0 };
        }
    }

    // Constructor inputs.
    protected mNumDimensions: number;  // N
    protected mDegree: number[];  // degree[N]
    protected mControls: IntpBSplineUniformControls<T>;
    protected mCTZero: T;
    protected mCacheMode: IntpBSplineUniformCacheMode;

    // Parameters for B-spline evaluation. All arrays have N elements.
    protected mNumLocalControls: number;  // product of (degree[]+1)
    protected mDegreeP1: number[];
    protected mNumControls: number[];
    protected mTMin: number[];
    protected mTMax: number[];
    protected mBlender: number[][];
    protected mDCoefficient: number[][];
    protected mLMax: number[][];
    protected mPowerDSDT: number[][];
    protected mITuple: number[];
    protected mJTuple: number[];
    protected mKTuple: number[];
    protected mLTuple: number[];
    protected mSumIJTuple: number[];
    protected mUTuple: number[];
    protected mPTuple: number[];

    // Support for non-cached B-spline evaluation. The outer array has N
    // elements.
    protected mPhi: number[][];

    // Support for cached B-spline evaluation.
    protected mTBound: number[];  // tbound[2*N]
    protected mComputeJTuple: number[];  // computejtuple[N]
    protected mComputeSumIJTuple: number[];  // computesumijtuple[N]
    protected mDegreeMinusOrder: number[];  // degreeminusorder[N]
    protected mTerm: T[];  // mTerm[N]
    protected mTensor: T[];  // depends on numcontrols
    protected mCached: boolean[];  // same size as mTensor

    // Abstract base class construction. The controls adapter is aliased,
    // not copied; the caller is responsible for ensuring that it persists as
    // long as the interpolator.
    protected constructor(numDimensions: number, degrees: readonly number[],
        controls: IntpBSplineUniformControls<T>, ctZero: T,
        cacheMode: IntpBSplineUniformCacheMode) {
        this.mNumDimensions = numDimensions;
        this.mControls = controls;
        this.mCTZero = ctZero;
        this.mCacheMode = cacheMode;
        this.mNumLocalControls = 0;

        const zeros = (n: number) => new Array<number>(n).fill(0);
        this.mDegree = zeros(numDimensions);
        this.mDegreeP1 = zeros(numDimensions);
        this.mNumControls = zeros(numDimensions);
        this.mTMin = zeros(numDimensions);
        this.mTMax = zeros(numDimensions);
        this.mBlender = [];
        this.mDCoefficient = [];
        this.mLMax = [];
        this.mPowerDSDT = [];
        this.mITuple = zeros(numDimensions);
        this.mJTuple = zeros(numDimensions);
        this.mKTuple = zeros(numDimensions);
        this.mLTuple = zeros(numDimensions);
        this.mSumIJTuple = zeros(numDimensions);
        this.mUTuple = zeros(numDimensions);
        this.mPTuple = zeros(numDimensions);
        this.mPhi = [];
        this.mTBound = [];
        this.mComputeJTuple = [];
        this.mComputeSumIJTuple = [];
        this.mDegreeMinusOrder = [];
        this.mTerm = [];
        this.mTensor = [];
        this.mCached = [];

        // The condition c+1 > d+1 is required so that when s = c+1-d, its
        // maximum value, we have at least two s-knots (d and d + 1).
        for (let dim = 0; dim < numDimensions; ++dim) {
            if (controls.getSize(dim) <= degrees[dim] + 1) {
                logError('Incompatible degree and number of controls.');
            }
        }

        this.mNumLocalControls = 1;
        for (let dim = 0; dim < numDimensions; ++dim) {
            this.mDegree[dim] = degrees[dim];
            this.mDegreeP1[dim] = degrees[dim] + 1;
            this.mNumLocalControls *= this.mDegreeP1[dim];
            this.mNumControls[dim] = controls.getSize(dim);
            this.mTMin[dim] = -0.5;
            this.mTMax[dim] = this.mNumControls[dim] - 0.5;
            this.mBlender.push(IntpBSplineUniformShared.computeBlendingMatrix(
                this.mDegree[dim]));
            const dc = IntpBSplineUniformShared.computeDCoefficients(this.mDegree[dim]);
            this.mDCoefficient.push(dc.dCoefficients);
            this.mLMax.push(dc.ellMax);
            this.mPowerDSDT.push(IntpBSplineUniformShared.computePowers(
                this.mDegree[dim], this.mNumControls[dim], this.mTMin[dim],
                this.mTMax[dim]));
        }

        if (this.mCacheMode === IntpBSplineUniformCacheMode.NO_CACHING) {
            for (let dim = 0; dim < numDimensions; ++dim) {
                this.mPhi.push(new Array<number>(this.mDegreeP1[dim]).fill(0));
            }
        }
        else {
            this.initializeTensors();
        }
    }

    // Member access.
    getDegree(dim: number): number {
        return this.mDegree[dim];
    }

    getNumControls(dim: number): number {
        return this.mNumControls[dim];
    }

    getTMin(dim: number): number {
        return this.mTMin[dim];
    }

    getTMax(dim: number): number {
        return this.mTMax[dim];
    }

    getCacheMode(): IntpBSplineUniformCacheMode {
        return this.mCacheMode;
    }

    // For the multidimensional tensor Phi(iTuple, kTuple), compute the
    // portion of the 1-dimensional index that corresponds to iTuple.
    protected getRowIndex(i: readonly number[]): number {
        let rowIndex = i[this.mNumDimensions - 1];
        let j1 = 2 * this.mNumDimensions - 2;
        for (let j0 = this.mNumDimensions - 2; j0 >= 0; --j0, --j1) {
            rowIndex = this.mTBound[j1] * rowIndex + i[j0];
        }
        rowIndex = this.mTBound[j1] * rowIndex;
        return rowIndex;
    }

    // For the multidimensional tensor Phi(iTuple, kTuple), combine the
    // getRowIndex(...) output with kTuple to produce the full 1-dimensional
    // index.
    protected getIndex(rowIndex: number, k: readonly number[]): number {
        let index = rowIndex + k[this.mNumDimensions - 1];
        for (let j = this.mNumDimensions - 2; j >= 0; --j) {
            index = this.mTBound[j] * index + k[j];
        }
        return index;
    }

    // Compute Phi(iTuple, kTuple). The 'index' value is an already computed
    // 1-dimensional index for the tensor.
    protected computeTensor(i: readonly number[], k: readonly number[],
        index: number): void {
        const controls = this.mControls;
        let element = this.mCTZero;
        for (let dim = 0; dim < this.mNumDimensions; ++dim) {
            this.mComputeJTuple[dim] = 0;
        }
        for (let iterate = 0; iterate < this.mNumLocalControls; ++iterate) {
            let blend = 1;
            for (let dim = 0; dim < this.mNumDimensions; ++dim) {
                blend *= this.mBlender[dim][k[dim]
                    + this.mDegreeP1[dim] * this.mComputeJTuple[dim]];
                this.mComputeSumIJTuple[dim] = i[dim] + this.mComputeJTuple[dim];
            }
            element = controls.add(element,
                controls.mul(controls.get(this.mComputeSumIJTuple), blend));

            for (let dim = 0; dim < this.mNumDimensions; ++dim) {
                if (++this.mComputeJTuple[dim] < this.mDegreeP1[dim]) {
                    break;
                }
                this.mComputeJTuple[dim] = 0;
            }
        }
        this.mTensor[index] = element;
    }

    // Allocate the containers used for caching and fill in the tensor for
    // precaching when that mode is selected.
    protected initializeTensors(): void {
        this.mTBound = new Array<number>(2 * this.mNumDimensions).fill(0);
        this.mComputeJTuple = new Array<number>(this.mNumDimensions).fill(0);
        this.mComputeSumIJTuple = new Array<number>(this.mNumDimensions).fill(0);
        this.mDegreeMinusOrder = new Array<number>(this.mNumDimensions).fill(0);
        this.mTerm = new Array<T>(this.mNumDimensions).fill(this.mCTZero);

        let current = 0;
        let numCached = 1;
        for (let dim = 0; dim < this.mNumDimensions; ++dim, ++current) {
            this.mTBound[current] = this.mDegreeP1[dim];
            numCached *= this.mTBound[current];
        }
        for (let dim = 0; dim < this.mNumDimensions; ++dim, ++current) {
            this.mTBound[current] = this.mNumControls[dim] - this.mDegree[dim];
            numCached *= this.mTBound[current];
        }
        this.mTensor = new Array<T>(numCached).fill(this.mCTZero);
        this.mCached = new Array<boolean>(numCached).fill(false);
        if (this.mCacheMode === IntpBSplineUniformCacheMode.PRE_CACHING) {
            const tuple = new Array<number>(2 * this.mNumDimensions).fill(0);
            const i = tuple.slice(this.mNumDimensions);
            const k = tuple.slice(0, this.mNumDimensions);
            for (let index = 0; index < numCached; ++index) {
                for (let dim = 0; dim < this.mNumDimensions; ++dim) {
                    k[dim] = tuple[dim];
                    i[dim] = tuple[this.mNumDimensions + dim];
                }
                this.computeTensor(i, k, index);
                for (let j = 0; j < 2 * this.mNumDimensions; ++j) {
                    if (++tuple[j] < this.mTBound[j]) {
                        break;
                    }
                    tuple[j] = 0;
                }
            }
            this.mCached.fill(true);
        }
        else {
            this.mCached.fill(false);
        }
    }

    // Evaluate the interpolator. Each element of 'order' indicates the
    // order of the derivative you want to compute. For the function value
    // itself, pass in 'order' that has all 0 elements.
    protected evaluateNoCaching(order: readonly number[], t: readonly number[]): T {
        const controls = this.mControls;
        let result = this.mCTZero;
        for (let dim = 0; dim < this.mNumDimensions; ++dim) {
            if (order[dim] < 0 || order[dim] > this.mDegree[dim]) {
                return result;
            }
        }

        for (let dim = 0; dim < this.mNumDimensions; ++dim) {
            const key = IntpBSplineUniformShared.getKey(t[dim], this.mTMin[dim],
                this.mTMax[dim], this.mPowerDSDT[dim][1], this.mNumControls[dim],
                this.mDegree[dim]);
            this.mITuple[dim] = key.index;
            this.mUTuple[dim] = key.u;
        }

        for (let dim = 0; dim < this.mNumDimensions; ++dim) {
            let jIndex = 0;
            for (let j = 0; j <= this.mDegree[dim]; ++j) {
                let kjIndex = this.mDegree[dim] + jIndex;
                let ell = this.mLMax[dim][order[dim]];
                this.mPhi[dim][j] = 0;
                for (let k = this.mDegree[dim]; k >= order[dim]; --k) {
                    this.mPhi[dim][j] = this.mPhi[dim][j] * this.mUTuple[dim]
                        + this.mBlender[dim][kjIndex--] * this.mDCoefficient[dim][ell--];
                }
                jIndex += this.mDegreeP1[dim];
            }
        }

        for (let dim = 0; dim < this.mNumDimensions; ++dim) {
            this.mJTuple[dim] = 0;
            this.mSumIJTuple[dim] = this.mITuple[dim];
            this.mPTuple[dim] = this.mPhi[dim][0];
        }
        for (let iterate = 0; iterate < this.mNumLocalControls; ++iterate) {
            let product = 1;
            for (let dim = 0; dim < this.mNumDimensions; ++dim) {
                product *= this.mPTuple[dim];
            }

            result = controls.add(result,
                controls.mul(controls.get(this.mSumIJTuple), product));

            for (let dim = 0; dim < this.mNumDimensions; ++dim) {
                if (++this.mJTuple[dim] <= this.mDegree[dim]) {
                    this.mSumIJTuple[dim] = this.mITuple[dim] + this.mJTuple[dim];
                    this.mPTuple[dim] = this.mPhi[dim][this.mJTuple[dim]];
                    break;
                }
                this.mJTuple[dim] = 0;
                this.mSumIJTuple[dim] = this.mITuple[dim];
                this.mPTuple[dim] = this.mPhi[dim][0];
            }
        }

        let adjust = 1;
        for (let dim = 0; dim < this.mNumDimensions; ++dim) {
            adjust *= this.mPowerDSDT[dim][order[dim]];
        }
        result = controls.mul(result, adjust);
        return result;
    }

    protected evaluateCaching(order: readonly number[], t: readonly number[]): T {
        const controls = this.mControls;
        let numIterates = 1;
        for (let dim = 0; dim < this.mNumDimensions; ++dim) {
            this.mDegreeMinusOrder[dim] = this.mDegree[dim] - order[dim];
            if (this.mDegreeMinusOrder[dim] < 0
                || this.mDegreeMinusOrder[dim] > this.mDegree[dim]) {
                return this.mCTZero;
            }
            numIterates *= this.mDegreeMinusOrder[dim] + 1;
        }

        for (let dim = 0; dim < this.mNumDimensions; ++dim) {
            const key = IntpBSplineUniformShared.getKey(t[dim], this.mTMin[dim],
                this.mTMax[dim], this.mPowerDSDT[dim][1], this.mNumControls[dim],
                this.mDegree[dim]);
            this.mITuple[dim] = key.index;
            this.mUTuple[dim] = key.u;
        }

        const rowIndex = this.getRowIndex(this.mITuple);
        for (let dim = 0; dim < this.mNumDimensions; ++dim) {
            this.mJTuple[dim] = 0;
            this.mKTuple[dim] = this.mDegree[dim];
            this.mLTuple[dim] = this.mLMax[dim][order[dim]];
            this.mTerm[dim] = this.mCTZero;
        }
        for (let iterate = 0; iterate < numIterates; ++iterate) {
            const index = this.getIndex(rowIndex, this.mKTuple);
            if (this.mCacheMode === IntpBSplineUniformCacheMode.ON_DEMAND_CACHING
                && !this.mCached[index]) {
                this.computeTensor(this.mITuple, this.mKTuple, index);
                this.mCached[index] = true;
            }
            this.mTerm[0] = controls.add(controls.mul(this.mTerm[0], this.mUTuple[0]),
                controls.mul(this.mTensor[index], this.mDCoefficient[0][this.mLTuple[0]]));
            for (let dim = 0; dim < this.mNumDimensions; ++dim) {
                if (++this.mJTuple[dim] <= this.mDegreeMinusOrder[dim]) {
                    --this.mKTuple[dim];
                    --this.mLTuple[dim];
                    break;
                }
                const dimp1 = dim + 1;
                if (dimp1 < this.mNumDimensions) {
                    this.mTerm[dimp1] = controls.add(
                        controls.mul(this.mTerm[dimp1], this.mUTuple[dimp1]),
                        controls.mul(this.mTerm[dim],
                            this.mDCoefficient[dimp1][this.mLTuple[dimp1]]));
                    this.mTerm[dim] = this.mCTZero;
                    this.mJTuple[dim] = 0;
                    this.mKTuple[dim] = this.mDegree[dim];
                    this.mLTuple[dim] = this.mLMax[dim][order[dim]];
                }
            }
        }
        let result = this.mTerm[this.mNumDimensions - 1];

        let adjust = 1;
        for (let dim = 0; dim < this.mNumDimensions; ++dim) {
            adjust *= this.mPowerDSDT[dim][order[dim]];
        }
        result = controls.mul(result, adjust);
        return result;
    }
}

// Implementation for B-spline interpolation whose dimension is known only
// at run time.
export class IntpBSplineUniform<T> extends IntpBSplineUniformShared<T> {
    // The caller is responsible for ensuring that the IntpBSplineUniform
    // object persists as long as the input 'controls' exists.
    constructor(degrees: readonly number[], controls: IntpBSplineUniformControls<T>,
        ctZero: T, cacheMode: IntpBSplineUniformCacheMode) {
        super(degrees.length, degrees, controls, ctZero, cacheMode);
    }

    // Evaluate the interpolator. Each element of 'order' indicates the
    // order of the derivative you want to compute. For the function value
    // itself, pass in 'order' that has all 0 elements.
    evaluate(order: readonly number[], t: readonly number[]): T {
        if (order.length >= this.mNumDimensions && t.length >= this.mNumDimensions) {
            if (this.mCacheMode === IntpBSplineUniformCacheMode.NO_CACHING) {
                return this.evaluateNoCaching(order, t);
            }
            else {
                return this.evaluateCaching(order, t);
            }
        }
        else {
            return this.mCTZero;
        }
    }
}

// Specialization for 1-dimensional data.
export class IntpBSplineUniform1<T> {
    // Constructor inputs.
    protected mDegree: number;
    protected mControls: IntpBSplineUniformControls<T>;
    protected mCTZero: T;
    protected mCacheMode: IntpBSplineUniformCacheMode;

    // Parameters for B-spline evaluation.
    protected mDegreeP1: number;
    protected mNumControls: number;
    protected mTMin: number;
    protected mTMax: number;
    protected mBlender: number[];
    protected mDCoefficient: number[];
    protected mLMax: number[];
    protected mPowerDSDT: number[];

    // Support for non-cached B-spline evaluation.
    protected mPhi: number[];

    // Support for cached B-spline evaluation.
    protected mNumTRows: number;
    protected mNumTCols: number;
    protected mTensor: T[];
    protected mCached: boolean[];

    // The caller is responsible for ensuring that the IntpBSplineUniform1
    // object persists as long as the input 'controls' exists.
    constructor(degree: number, controls: IntpBSplineUniformControls<T>, ctZero: T,
        cacheMode: IntpBSplineUniformCacheMode) {
        this.mDegree = degree;
        this.mControls = controls;
        this.mCTZero = ctZero;
        this.mCacheMode = cacheMode;

        // The condition c+1 > d+1 is required so that when s = c+1-d, its
        // maximum value, we have at least two s-knots (d and d + 1).
        if (controls.getSize(0) <= degree + 1) {
            logError('Incompatible degree and number of controls.');
        }

        this.mDegreeP1 = degree + 1;
        this.mNumControls = controls.getSize(0);
        this.mTMin = -0.5;
        this.mTMax = this.mNumControls - 0.5;
        this.mNumTRows = 0;
        this.mNumTCols = 0;
        this.mBlender = IntpBSplineUniformShared.computeBlendingMatrix(degree);
        const dc = IntpBSplineUniformShared.computeDCoefficients(degree);
        this.mDCoefficient = dc.dCoefficients;
        this.mLMax = dc.ellMax;
        this.mPowerDSDT = IntpBSplineUniformShared.computePowers(degree,
            this.mNumControls, this.mTMin, this.mTMax);
        this.mPhi = [];
        this.mTensor = [];
        this.mCached = [];
        if (cacheMode === IntpBSplineUniformCacheMode.NO_CACHING) {
            this.mPhi = new Array<number>(this.mDegreeP1).fill(0);
        }
        else {
            this.initializeTensors();
        }
    }

    // Member access. The input dim is ignored; it exists so that the
    // interface matches the general-dimension class.
    getDegree(_dim: number = 0): number {
        return this.mDegree;
    }

    getNumControls(_dim: number = 0): number {
        return this.mNumControls;
    }

    getTMin(_dim: number = 0): number {
        return this.mTMin;
    }

    getTMax(_dim: number = 0): number {
        return this.mTMax;
    }

    getCacheMode(): IntpBSplineUniformCacheMode {
        return this.mCacheMode;
    }

    // Evaluate the interpolator. The order is 0 when you want the B-spline
    // function value itself. The order is 1 for the first derivative of the
    // function, and so on.
    evaluate(order: readonly number[], t: readonly number[]): T {
        const controls = this.mControls;
        let result = this.mCTZero;
        if (0 <= order[0] && order[0] <= this.mDegree) {
            const key = IntpBSplineUniformShared.getKey(t[0], this.mTMin, this.mTMax,
                this.mPowerDSDT[1], this.mNumControls, this.mDegree);
            const i = key.index;
            const u = key.u;

            if (this.mCacheMode === IntpBSplineUniformCacheMode.NO_CACHING) {
                let jIndex = 0;
                for (let j = 0; j <= this.mDegree; ++j) {
                    let kjIndex = this.mDegree + jIndex;
                    let ell = this.mLMax[order[0]];
                    this.mPhi[j] = 0;
                    for (let k = this.mDegree; k >= order[0]; --k) {
                        this.mPhi[j] = this.mPhi[j] * u
                            + this.mBlender[kjIndex--] * this.mDCoefficient[ell--];
                    }
                    jIndex += this.mDegreeP1;
                }

                for (let j = 0; j <= this.mDegree; ++j) {
                    result = controls.add(result,
                        controls.mul(controls.get([i + j]), this.mPhi[j]));
                }
            }
            else {
                const iIndex = this.mNumTCols * i;
                let kiIndex = this.mDegree + iIndex;
                let ell = this.mLMax[order[0]];
                for (let k = this.mDegree; k >= order[0]; --k) {
                    if (this.mCacheMode === IntpBSplineUniformCacheMode.ON_DEMAND_CACHING
                        && !this.mCached[kiIndex]) {
                        this.computeTensor(i, k, kiIndex);
                        this.mCached[kiIndex] = true;
                    }

                    result = controls.add(controls.mul(result, u),
                        controls.mul(this.mTensor[kiIndex--], this.mDCoefficient[ell--]));
                }
            }

            result = controls.mul(result, this.mPowerDSDT[order[0]]);
        }
        return result;
    }

    protected computeTensor(r: number, c: number, index: number): void {
        const controls = this.mControls;
        let element = this.mCTZero;
        for (let j = 0; j <= this.mDegree; ++j) {
            element = controls.add(element, controls.mul(controls.get([r + j]),
                this.mBlender[c + this.mDegreeP1 * j]));
        }
        this.mTensor[index] = element;
    }

    protected initializeTensors(): void {
        this.mNumTRows = this.mNumControls - this.mDegree;
        this.mNumTCols = this.mDegreeP1;
        const numCached = this.mNumTRows * this.mNumTCols;
        this.mTensor = new Array<T>(numCached).fill(this.mCTZero);
        this.mCached = new Array<boolean>(numCached).fill(false);
        if (this.mCacheMode === IntpBSplineUniformCacheMode.PRE_CACHING) {
            for (let r = 0, index = 0; r < this.mNumTRows; ++r) {
                for (let c = 0; c < this.mNumTCols; ++c, ++index) {
                    this.computeTensor(r, c, index);
                }
            }
            this.mCached.fill(true);
        }
        else {
            this.mCached.fill(false);
        }
    }
}

// Specialization for 2-dimensional data.
export class IntpBSplineUniform2<T> {
    // Constructor inputs.
    protected mDegree: number[];
    protected mControls: IntpBSplineUniformControls<T>;
    protected mCTZero: T;
    protected mCacheMode: IntpBSplineUniformCacheMode;

    // Parameters for B-spline evaluation.
    protected mDegreeP1: number[];
    protected mNumControls: number[];
    protected mTMin: number[];
    protected mTMax: number[];
    protected mBlender: number[][];
    protected mDCoefficient: number[][];
    protected mLMax: number[][];
    protected mPowerDSDT: number[][];

    // Support for non-cached B-spline evaluation.
    protected mPhi: number[][];

    // Support for cached B-spline evaluation.
    protected mNumTRows: number[];
    protected mNumTCols: number[];
    protected mTensor: T[];
    protected mCached: boolean[];

    // The caller is responsible for ensuring that the IntpBSplineUniform2
    // object persists as long as the input 'controls' exists.
    constructor(degrees: readonly number[], controls: IntpBSplineUniformControls<T>,
        ctZero: T, cacheMode: IntpBSplineUniformCacheMode) {
        this.mDegree = [degrees[0], degrees[1]];
        this.mControls = controls;
        this.mCTZero = ctZero;
        this.mCacheMode = cacheMode;

        // The condition c+1 > d+1 is required so that when s = c+1-d, its
        // maximum value, we have at least two s-knots (d and d + 1).
        for (let dim = 0; dim < 2; ++dim) {
            if (controls.getSize(dim) <= this.mDegree[dim] + 1) {
                logError('Incompatible degree and number of controls.');
            }
        }

        this.mDegreeP1 = [0, 0];
        this.mNumControls = [0, 0];
        this.mTMin = [0, 0];
        this.mTMax = [0, 0];
        this.mNumTRows = [0, 0];
        this.mNumTCols = [0, 0];
        this.mBlender = [];
        this.mDCoefficient = [];
        this.mLMax = [];
        this.mPowerDSDT = [];
        this.mPhi = [];
        this.mTensor = [];
        this.mCached = [];

        for (let dim = 0; dim < 2; ++dim) {
            this.mDegreeP1[dim] = this.mDegree[dim] + 1;
            this.mNumControls[dim] = controls.getSize(dim);
            this.mTMin[dim] = -0.5;
            this.mTMax[dim] = this.mNumControls[dim] - 0.5;
            this.mBlender.push(IntpBSplineUniformShared.computeBlendingMatrix(
                this.mDegree[dim]));
            const dc = IntpBSplineUniformShared.computeDCoefficients(this.mDegree[dim]);
            this.mDCoefficient.push(dc.dCoefficients);
            this.mLMax.push(dc.ellMax);
            this.mPowerDSDT.push(IntpBSplineUniformShared.computePowers(
                this.mDegree[dim], this.mNumControls[dim], this.mTMin[dim],
                this.mTMax[dim]));
        }

        if (cacheMode === IntpBSplineUniformCacheMode.NO_CACHING) {
            for (let dim = 0; dim < 2; ++dim) {
                this.mPhi.push(new Array<number>(this.mDegreeP1[dim]).fill(0));
            }
        }
        else {
            this.initializeTensors();
        }
    }

    // Member access.
    getDegree(dim: number): number {
        return this.mDegree[dim];
    }

    getNumControls(dim: number): number {
        return this.mNumControls[dim];
    }

    getTMin(dim: number): number {
        return this.mTMin[dim];
    }

    getTMax(dim: number): number {
        return this.mTMax[dim];
    }

    getCacheMode(): IntpBSplineUniformCacheMode {
        return this.mCacheMode;
    }

    // Evaluate the interpolator. The order is (0,0) when you want the
    // B-spline function value itself. The order0 is 1 for the first
    // derivative with respect to t0 and the order1 is 1 for the first
    // derivative with respect to t1. Higher-order derivatives in other
    // t-inputs are computed similarly.
    evaluate(order: readonly number[], t: readonly number[]): T {
        const controls = this.mControls;
        let result = this.mCTZero;
        if (0 <= order[0] && order[0] <= this.mDegree[0]
            && 0 <= order[1] && order[1] <= this.mDegree[1]) {
            const i = [0, 0];
            const u = [0, 0];
            for (let dim = 0; dim < 2; ++dim) {
                const key = IntpBSplineUniformShared.getKey(t[dim], this.mTMin[dim],
                    this.mTMax[dim], this.mPowerDSDT[dim][1], this.mNumControls[dim],
                    this.mDegree[dim]);
                i[dim] = key.index;
                u[dim] = key.u;
            }

            if (this.mCacheMode === IntpBSplineUniformCacheMode.NO_CACHING) {
                for (let dim = 0; dim < 2; ++dim) {
                    let jIndex = 0;
                    for (let j = 0; j <= this.mDegree[dim]; ++j) {
                        let kjIndex = this.mDegree[dim] + jIndex;
                        let ell = this.mLMax[dim][order[dim]];
                        this.mPhi[dim][j] = 0;
                        for (let k = this.mDegree[dim]; k >= order[dim]; --k) {
                            this.mPhi[dim][j] = this.mPhi[dim][j] * u[dim]
                                + this.mBlender[dim][kjIndex--]
                                * this.mDCoefficient[dim][ell--];
                        }
                        jIndex += this.mDegreeP1[dim];
                    }
                }

                for (let j1 = 0; j1 <= this.mDegree[1]; ++j1) {
                    const phi1 = this.mPhi[1][j1];
                    for (let j0 = 0; j0 <= this.mDegree[0]; ++j0) {
                        const phi0 = this.mPhi[0][j0];
                        const phi01 = phi0 * phi1;
                        result = controls.add(result, controls.mul(
                            controls.get([i[0] + j0, i[1] + j1]), phi01));
                    }
                }
            }
            else {
                const i0i1Index = this.mNumTCols[1] * (i[0] + this.mNumTRows[0] * i[1]);
                let k1i0i1Index = this.mDegree[1] + i0i1Index;
                let ell1 = this.mLMax[1][order[1]];
                for (let k1 = this.mDegree[1]; k1 >= order[1]; --k1) {
                    let k0k1i0i1Index = this.mDegree[0] + this.mNumTCols[0] * k1i0i1Index;
                    let ell0 = this.mLMax[0][order[0]];
                    let term = this.mCTZero;
                    for (let k0 = this.mDegree[0]; k0 >= order[0]; --k0) {
                        if (this.mCacheMode === IntpBSplineUniformCacheMode.ON_DEMAND_CACHING
                            && !this.mCached[k0k1i0i1Index]) {
                            this.computeTensor(i[0], i[1], k0, k1, k0k1i0i1Index);
                            this.mCached[k0k1i0i1Index] = true;
                        }
                        term = controls.add(controls.mul(term, u[0]), controls.mul(
                            this.mTensor[k0k1i0i1Index--], this.mDCoefficient[0][ell0--]));
                    }
                    result = controls.add(controls.mul(result, u[1]),
                        controls.mul(term, this.mDCoefficient[1][ell1--]));
                    --k1i0i1Index;
                }
            }

            let adjust = 1;
            for (let dim = 0; dim < 2; ++dim) {
                adjust *= this.mPowerDSDT[dim][order[dim]];
            }
            result = controls.mul(result, adjust);
        }
        return result;
    }

    computeTensor(r0: number, r1: number, c0: number, c1: number, index: number): void {
        const controls = this.mControls;
        let element = this.mCTZero;
        for (let j1 = 0; j1 <= this.mDegree[1]; ++j1) {
            const blend1 = this.mBlender[1][c1 + this.mDegreeP1[1] * j1];
            for (let j0 = 0; j0 <= this.mDegree[0]; ++j0) {
                const blend0 = this.mBlender[0][c0 + this.mDegreeP1[0] * j0];
                const blend01 = blend0 * blend1;
                element = controls.add(element, controls.mul(
                    controls.get([r0 + j0, r1 + j1]), blend01));
            }
        }
        this.mTensor[index] = element;
    }

    protected initializeTensors(): void {
        let numCached = 1;
        for (let dim = 0; dim < 2; ++dim) {
            this.mNumTRows[dim] = this.mNumControls[dim] - this.mDegree[dim];
            this.mNumTCols[dim] = this.mDegreeP1[dim];
            numCached *= this.mNumTRows[dim] * this.mNumTCols[dim];
        }
        this.mTensor = new Array<T>(numCached).fill(this.mCTZero);
        this.mCached = new Array<boolean>(numCached).fill(false);
        if (this.mCacheMode === IntpBSplineUniformCacheMode.PRE_CACHING) {
            for (let r1 = 0, index = 0; r1 < this.mNumTRows[1]; ++r1) {
                for (let r0 = 0; r0 < this.mNumTRows[0]; ++r0) {
                    for (let c1 = 0; c1 < this.mNumTCols[1]; ++c1) {
                        for (let c0 = 0; c0 < this.mNumTCols[0]; ++c0, ++index) {
                            this.computeTensor(r0, r1, c0, c1, index);
                        }
                    }
                }
            }
            this.mCached.fill(true);
        }
        else {
            this.mCached.fill(false);
        }
    }
}

// Specialization for 3-dimensional data.
export class IntpBSplineUniform3<T> {
    // Constructor inputs.
    protected mDegree: number[];
    protected mControls: IntpBSplineUniformControls<T>;
    protected mCTZero: T;
    protected mCacheMode: IntpBSplineUniformCacheMode;

    // Parameters for B-spline evaluation.
    protected mDegreeP1: number[];
    protected mNumControls: number[];
    protected mTMin: number[];
    protected mTMax: number[];
    protected mBlender: number[][];
    protected mDCoefficient: number[][];
    protected mLMax: number[][];
    protected mPowerDSDT: number[][];

    // Support for non-cached B-spline evaluation.
    protected mPhi: number[][];

    // Support for cached B-spline evaluation.
    protected mNumTRows: number[];
    protected mNumTCols: number[];
    protected mTensor: T[];
    protected mCached: boolean[];

    // The caller is responsible for ensuring that the IntpBSplineUniform3
    // object persists as long as the input 'controls' exists.
    constructor(degrees: readonly number[], controls: IntpBSplineUniformControls<T>,
        ctZero: T, cacheMode: IntpBSplineUniformCacheMode) {
        this.mDegree = [degrees[0], degrees[1], degrees[2]];
        this.mControls = controls;
        this.mCTZero = ctZero;
        this.mCacheMode = cacheMode;

        // The condition c+1 > d+1 is required so that when s = c+1-d, its
        // maximum value, we have at least two s-knots (d and d + 1).
        for (let dim = 0; dim < 3; ++dim) {
            if (controls.getSize(dim) <= this.mDegree[dim] + 1) {
                logError('Incompatible degree and number of controls.');
            }
        }

        this.mDegreeP1 = [0, 0, 0];
        this.mNumControls = [0, 0, 0];
        this.mTMin = [0, 0, 0];
        this.mTMax = [0, 0, 0];
        this.mNumTRows = [0, 0, 0];
        this.mNumTCols = [0, 0, 0];
        this.mBlender = [];
        this.mDCoefficient = [];
        this.mLMax = [];
        this.mPowerDSDT = [];
        this.mPhi = [];
        this.mTensor = [];
        this.mCached = [];

        for (let dim = 0; dim < 3; ++dim) {
            this.mDegreeP1[dim] = this.mDegree[dim] + 1;
            this.mNumControls[dim] = controls.getSize(dim);
            this.mTMin[dim] = -0.5;
            this.mTMax[dim] = this.mNumControls[dim] - 0.5;
            this.mBlender.push(IntpBSplineUniformShared.computeBlendingMatrix(
                this.mDegree[dim]));
            const dc = IntpBSplineUniformShared.computeDCoefficients(this.mDegree[dim]);
            this.mDCoefficient.push(dc.dCoefficients);
            this.mLMax.push(dc.ellMax);
            this.mPowerDSDT.push(IntpBSplineUniformShared.computePowers(
                this.mDegree[dim], this.mNumControls[dim], this.mTMin[dim],
                this.mTMax[dim]));
        }

        if (cacheMode === IntpBSplineUniformCacheMode.NO_CACHING) {
            for (let dim = 0; dim < 3; ++dim) {
                this.mPhi.push(new Array<number>(this.mDegreeP1[dim]).fill(0));
            }
        }
        else {
            this.initializeTensors();
        }
    }

    // Member access. The input dim specifies the dimension (0, 1, 2).
    getDegree(dim: number): number {
        return this.mDegree[dim];
    }

    getNumControls(dim: number): number {
        return this.mNumControls[dim];
    }

    getTMin(dim: number): number {
        return this.mTMin[dim];
    }

    getTMax(dim: number): number {
        return this.mTMax[dim];
    }

    getCacheMode(): IntpBSplineUniformCacheMode {
        return this.mCacheMode;
    }

    // Evaluate the interpolator. The order is (0,0,0) when you want the
    // B-spline function value itself. The order0 is 1 for the first
    // derivative with respect to t0, the order1 is 1 for the first
    // derivative with respect to t1 or the order2 is 1 for the first
    // derivative with respect to t2. Higher-order derivatives in other
    // t-inputs are computed similarly.
    evaluate(order: readonly number[], t: readonly number[]): T {
        const controls = this.mControls;
        let result = this.mCTZero;
        if (0 <= order[0] && order[0] <= this.mDegree[0]
            && 0 <= order[1] && order[1] <= this.mDegree[1]
            && 0 <= order[2] && order[2] <= this.mDegree[2]) {
            const i = [0, 0, 0];
            const u = [0, 0, 0];
            for (let dim = 0; dim < 3; ++dim) {
                const key = IntpBSplineUniformShared.getKey(t[dim], this.mTMin[dim],
                    this.mTMax[dim], this.mPowerDSDT[dim][1], this.mNumControls[dim],
                    this.mDegree[dim]);
                i[dim] = key.index;
                u[dim] = key.u;
            }

            if (this.mCacheMode === IntpBSplineUniformCacheMode.NO_CACHING) {
                for (let dim = 0; dim < 3; ++dim) {
                    let jIndex = 0;
                    for (let j = 0; j <= this.mDegree[dim]; ++j) {
                        let kjIndex = this.mDegree[dim] + jIndex;
                        let ell = this.mLMax[dim][order[dim]];
                        this.mPhi[dim][j] = 0;
                        for (let k = this.mDegree[dim]; k >= order[dim]; --k) {
                            this.mPhi[dim][j] = this.mPhi[dim][j] * u[dim]
                                + this.mBlender[dim][kjIndex--]
                                * this.mDCoefficient[dim][ell--];
                        }
                        jIndex += this.mDegreeP1[dim];
                    }
                }

                for (let j2 = 0; j2 <= this.mDegree[2]; ++j2) {
                    const phi2 = this.mPhi[2][j2];
                    for (let j1 = 0; j1 <= this.mDegree[1]; ++j1) {
                        const phi1 = this.mPhi[1][j1];
                        const phi12 = phi1 * phi2;
                        for (let j0 = 0; j0 <= this.mDegree[0]; ++j0) {
                            const phi0 = this.mPhi[0][j0];
                            const phi012 = phi0 * phi12;
                            result = controls.add(result, controls.mul(
                                controls.get([i[0] + j0, i[1] + j1, i[2] + j2]), phi012));
                        }
                    }
                }
            }
            else {
                const i0i1i2Index = this.mNumTCols[2] * (i[0] + this.mNumTRows[0]
                    * (i[1] + this.mNumTRows[1] * i[2]));
                let k2i0i1i2Index = this.mDegree[2] + i0i1i2Index;
                let ell2 = this.mLMax[2][order[2]];
                for (let k2 = this.mDegree[2]; k2 >= order[2]; --k2) {
                    let k1k2i0i1i2Index = this.mDegree[1] + this.mNumTCols[1] * k2i0i1i2Index;
                    let ell1 = this.mLMax[1][order[1]];
                    let term1 = this.mCTZero;
                    for (let k1 = this.mDegree[1]; k1 >= order[1]; --k1) {
                        let k0k1k2i0i1i2Index = this.mDegree[0]
                            + this.mNumTCols[0] * k1k2i0i1i2Index;
                        let ell0 = this.mLMax[0][order[0]];
                        let term0 = this.mCTZero;
                        for (let k0 = this.mDegree[0]; k0 >= order[0]; --k0) {
                            if (this.mCacheMode === IntpBSplineUniformCacheMode.ON_DEMAND_CACHING
                                && !this.mCached[k0k1k2i0i1i2Index]) {
                                this.computeTensor(i[0], i[1], i[2], k0, k1, k2,
                                    k0k1k2i0i1i2Index);
                                this.mCached[k0k1k2i0i1i2Index] = true;
                            }

                            term0 = controls.add(controls.mul(term0, u[0]),
                                controls.mul(this.mTensor[k0k1k2i0i1i2Index--],
                                    this.mDCoefficient[0][ell0--]));
                        }
                        term1 = controls.add(controls.mul(term1, u[1]),
                            controls.mul(term0, this.mDCoefficient[1][ell1--]));
                        --k1k2i0i1i2Index;
                    }
                    result = controls.add(controls.mul(result, u[2]),
                        controls.mul(term1, this.mDCoefficient[2][ell2--]));
                    --k2i0i1i2Index;
                }
            }

            let adjust = 1;
            for (let dim = 0; dim < 3; ++dim) {
                adjust *= this.mPowerDSDT[dim][order[dim]];
            }
            result = controls.mul(result, adjust);
        }
        return result;
    }

    protected computeTensor(r0: number, r1: number, r2: number, c0: number,
        c1: number, c2: number, index: number): void {
        const controls = this.mControls;
        let element = this.mCTZero;
        for (let j2 = 0; j2 <= this.mDegree[2]; ++j2) {
            const blend2 = this.mBlender[2][c2 + this.mDegreeP1[2] * j2];
            for (let j1 = 0; j1 <= this.mDegree[1]; ++j1) {
                const blend1 = this.mBlender[1][c1 + this.mDegreeP1[1] * j1];
                const blend12 = blend1 * blend2;
                for (let j0 = 0; j0 <= this.mDegree[0]; ++j0) {
                    const blend0 = this.mBlender[0][c0 + this.mDegreeP1[0] * j0];
                    const blend012 = blend0 * blend12;
                    element = controls.add(element, controls.mul(
                        controls.get([r0 + j0, r1 + j1, r2 + j2]), blend012));
                }
            }
        }
        this.mTensor[index] = element;
    }

    protected initializeTensors(): void {
        let numCached = 1;
        for (let dim = 0; dim < 3; ++dim) {
            this.mNumTRows[dim] = this.mNumControls[dim] - this.mDegree[dim];
            this.mNumTCols[dim] = this.mDegreeP1[dim];
            numCached *= this.mNumTRows[dim] * this.mNumTCols[dim];
        }
        this.mTensor = new Array<T>(numCached).fill(this.mCTZero);
        this.mCached = new Array<boolean>(numCached).fill(false);
        if (this.mCacheMode === IntpBSplineUniformCacheMode.PRE_CACHING) {
            for (let r2 = 0, index = 0; r2 < this.mNumTRows[2]; ++r2) {
                for (let r1 = 0; r1 < this.mNumTRows[1]; ++r1) {
                    for (let r0 = 0; r0 < this.mNumTRows[0]; ++r0) {
                        for (let c2 = 0; c2 < this.mNumTCols[2]; ++c2) {
                            for (let c1 = 0; c1 < this.mNumTCols[1]; ++c1) {
                                for (let c0 = 0; c0 < this.mNumTCols[0]; ++c0, ++index) {
                                    this.computeTensor(r0, r1, r2, c0, c1, c2, index);
                                }
                            }
                        }
                    }
                }
            }
            this.mCached.fill(true);
        }
        else {
            this.mCached.fill(false);
        }
    }
}
