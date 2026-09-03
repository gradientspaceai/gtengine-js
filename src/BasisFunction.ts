// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BasisFunction.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

import { logAssert, logError } from './Logger.js';
import { Array2 } from './Array2.js';

// A unique knot value together with its multiplicity in the knot vector.
export class UniqueKnot {
    public t: number;
    public multiplicity: number;

    constructor(t: number = 0, multiplicity: number = 0) {
        this.t = t;
        this.multiplicity = multiplicity;
    }
}

export class BasisFunctionInput {
    public numControls: number;
    public degree: number;
    public uniform: boolean;
    public periodic: boolean;
    public numUniqueKnots: number;
    public uniqueKnots: UniqueKnot[];

    // With no arguments, the members are zero-initialized (upstream leaves
    // them "uninitialized" but its default constructor value-initializes).
    // With both arguments, construct an open uniform curve with t in [0,1].
    constructor(inNumControls?: number, inDegree?: number) {
        if (inNumControls === undefined || inDegree === undefined) {
            this.numControls = 0;
            this.degree = 0;
            this.uniform = false;
            this.periodic = false;
            this.numUniqueKnots = 0;
            this.uniqueKnots = [];
            return;
        }

        this.numControls = inNumControls;
        this.degree = inDegree;
        this.uniform = true;
        this.periodic = false;
        this.numUniqueKnots = this.numControls - this.degree + 1;
        this.uniqueKnots = new Array<UniqueKnot>(this.numUniqueKnots);
        for (let i = 0; i < this.numUniqueKnots; ++i) {
            this.uniqueKnots[i] = new UniqueKnot();
        }
        this.uniqueKnots[0].t = 0;
        this.uniqueKnots[0].multiplicity = this.degree + 1;
        for (let i = 1; i <= this.numUniqueKnots - 2; ++i) {
            this.uniqueKnots[i].t = i / (this.numUniqueKnots - 1);
            this.uniqueKnots[i].multiplicity = 1;
        }
        this.uniqueKnots[this.numUniqueKnots - 1].t = 1;
        this.uniqueKnots[this.numUniqueKnots - 1].multiplicity = this.degree + 1;
    }
}

export class BasisFunction {
    // Let n be the number of control points. Let d be the degree, where
    // 1 <= d <= n-1.  The number of knots is k = n + d + 1.  The knots
    // are t[i] for 0 <= i < k and must be nondecreasing, t[i] <= t[i+1],
    // but a knot value can be repeated.  Let s be the number of distinct
    // knots.  Let the distinct knots be u[j] for 0 <= j < s, so u[j] <
    // u[j+1] for all j.  The set of u[j] is called a 'breakpoint
    // sequence'.  Let m[j] >= 1 be the multiplicity; that is, if t[i] is
    // the first occurrence of u[j], then t[i+r] = t[i] for 1 <= r < m[j].
    // The multiplicities have the constraints m[0] <= d+1, m[s-1] <= d+1,
    // and m[j] <= d for 1 <= j <= s-2.  Also, k = sum_{j=0}^{s-1} m[j],
    // which says the multiplicities account for all k knots.
    //
    // Given a knot vector (t[0],...,t[n+d]), the domain of the
    // corresponding B-spline curve is the interval [t[d],t[n]].
    //
    // The corresponding B-spline or NURBS curve is characterized as
    // follows.  See "Geometric Modeling with Splines: An Introduction" by
    // Elaine Cohen, Richard F. Riesenfeld and Gershon Elber, AK Peters,
    // 2001, Natick MA.  The curve is 'open' when m[0] = m[s-1] = d+1;
    // otherwise, it is 'floating'.  An open curve is uniform when the
    // knots t[d] through t[n] are equally spaced; that is, t[i+1] - t[i]
    // are a common value for d <= i <= n-1.  By implication, s = n-d+1
    // and m[j] = 1 for 1 <= j <= s-2.  An open curve that does not
    // satisfy these conditions is said to be nonuniform.  A floating
    // curve is uniform when m[j] = 1 for 0 <= j <= s-1 and t[i+1] - t[i]
    // are a common value for 0 <= i <= k-2; otherwise, the floating curve
    // is nonuniform.
    //
    // A special case of a floating curve is a periodic curve.  The intent
    // is that the curve is closed, so the first and last control points
    // should be the same, which ensures C^{0} continuity.  Higher-order
    // continuity is obtained by repeating more control points.  If the
    // control points are P[0] through P[n-1], append the points P[0]
    // through P[d-1] to ensure C^{d-1} continuity.  Additionally, the
    // knots must be chosen properly.  You may choose t[d] through t[n] as
    // you wish.  The other knots are defined by
    //   t[i] - t[i-1] = t[n-d+i] - t[n-d+i-1]
    //   t[n+i] - t[n+i-1] = t[d+i] - t[d+i-1]
    // for 1 <= i <= d.

    // Constructor inputs and values derived from them.
    private mNumControls: number;
    private mDegree: number;
    private mTMin: number;
    private mTMax: number;
    private mTLength: number;
    private mOpen: boolean;
    private mUniform: boolean;
    private mPeriodic: boolean;
    private mUniqueKnots: UniqueKnot[];
    private mKnots: number[];

    // Lookup information for the getIndex() function.  The 't' member is a
    // unique knot value.  The 'lastIndex' member is the index in mKnots[]
    // for the last occurrence of that knot value.
    private mKeys: { t: number; lastIndex: number }[];

    // Storage for the basis functions and their first three derivatives;
    // mJet[i] is array[d+1][n+d].
    private mJet: Array2<number>[];

    // Construction.  The determination that the curve is open or floating
    // is based on the multiplicities.  The 'uniform' input is used to avoid
    // misclassifications due to floating-point rounding errors.
    // Specifically, the breakpoints might be equally spaced (uniform) as
    // real numbers, but the floating-point representations can have
    // rounding errors that cause the knot differences not to be exactly the
    // same constant.  A periodic curve can have uniform or nonuniform
    // knots.  This object makes copies of the input arrays.
    constructor(input?: BasisFunctionInput) {
        this.mNumControls = 0;
        this.mDegree = 0;
        this.mTMin = 0;
        this.mTMax = 0;
        this.mTLength = 0;
        this.mOpen = false;
        this.mUniform = false;
        this.mPeriodic = false;
        this.mUniqueKnots = [];
        this.mKnots = [];
        this.mKeys = [];
        this.mJet = [new Array2<number>(), new Array2<number>(),
            new Array2<number>(), new Array2<number>()];

        if (input !== undefined) {
            this.create(input);
        }
    }

    // Support for explicit creation in classes that have array members
    // involving BasisFunction.  This is a call-once function.
    create(input: BasisFunctionInput): void {
        logAssert(this.mNumControls === 0 && this.mDegree === 0, 'Object already created.');
        logAssert(input.numControls >= 2, 'Invalid number of control points.');
        logAssert(1 <= input.degree && input.degree < input.numControls, 'Invalid degree.');
        logAssert(input.numUniqueKnots >= 2, 'Invalid number of unique knots.');

        this.mNumControls = (input.periodic ? input.numControls + input.degree : input.numControls);
        this.mDegree = input.degree;
        this.mTMin = 0;
        this.mTMax = 0;
        this.mTLength = 0;
        this.mOpen = false;
        this.mUniform = input.uniform;
        this.mPeriodic = input.periodic;
        for (let i = 0; i < 4; ++i) {
            this.mJet[i] = new Array2<number>();
        }

        this.mUniqueKnots = new Array<UniqueKnot>(input.numUniqueKnots);
        for (let i = 0; i < input.numUniqueKnots; ++i) {
            // C++ copies by value; clone explicitly.
            this.mUniqueKnots[i] = new UniqueKnot(
                input.uniqueKnots[i].t, input.uniqueKnots[i].multiplicity);
        }

        let u = this.mUniqueKnots[0].t;
        for (let i = 1; i < input.numUniqueKnots - 1; ++i) {
            const uNext = this.mUniqueKnots[i].t;
            logAssert(u < uNext, 'Unique knots are not strictly increasing.');
            u = uNext;
        }

        const mult0 = this.mUniqueKnots[0].multiplicity;
        logAssert(mult0 >= 1 && mult0 <= this.mDegree + 1, 'Invalid first multiplicity.');

        const mult1 = this.mUniqueKnots[this.mUniqueKnots.length - 1].multiplicity;
        logAssert(mult1 >= 1 && mult1 <= this.mDegree + 1, 'Invalid last multiplicity.');

        for (let i = 1; i <= input.numUniqueKnots - 2; ++i) {
            const mult = this.mUniqueKnots[i].multiplicity;
            logAssert(mult >= 1 && mult <= this.mDegree + 1, 'Invalid interior multiplicity.');
        }

        this.mOpen = (mult0 === mult1 && mult0 === this.mDegree + 1);

        this.mKnots = new Array<number>(this.mNumControls + this.mDegree + 1);
        this.mKeys = new Array<{ t: number; lastIndex: number }>(input.numUniqueKnots);
        let sum = 0;
        for (let i = 0, j = 0; i < input.numUniqueKnots; ++i) {
            const tCommon = this.mUniqueKnots[i].t;
            const mult = this.mUniqueKnots[i].multiplicity;
            for (let k = 0; k < mult; ++k, ++j) {
                this.mKnots[j] = tCommon;
            }

            this.mKeys[i] = { t: tCommon, lastIndex: sum - 1 };
            sum += mult;
        }

        this.mTMin = this.mKnots[this.mDegree];
        this.mTMax = this.mKnots[this.mNumControls];
        this.mTLength = this.mTMax - this.mTMin;

        const numRows = this.mDegree + 1;
        const numCols = this.mNumControls + this.mDegree;
        for (let i = 0; i < 4; ++i) {
            this.mJet[i] = new Array2<number>(numCols, numRows);
            this.mJet[i].fill(0);
        }
    }

    // Member access.
    getNumControls(): number {
        return this.mNumControls;
    }

    getDegree(): number {
        return this.mDegree;
    }

    getNumUniqueKnots(): number {
        return this.mUniqueKnots.length;
    }

    getNumKnots(): number {
        return this.mKnots.length;
    }

    getMinDomain(): number {
        return this.mTMin;
    }

    getMaxDomain(): number {
        return this.mTMax;
    }

    isOpen(): boolean {
        return this.mOpen;
    }

    isUniform(): boolean {
        return this.mUniform;
    }

    isPeriodic(): boolean {
        return this.mPeriodic;
    }

    // The returned arrays alias internal storage (upstream returns const
    // pointers); do not modify them.
    getUniqueKnots(): readonly UniqueKnot[] {
        return this.mUniqueKnots;
    }

    getKnots(): readonly number[] {
        return this.mKnots;
    }

    // Evaluation of the basis function and its derivatives through
    // order 3.  For the function value only, pass order 0.  For the
    // function and first derivative, pass order 1, and so on.  The C++
    // int32_t& outputs minIndex and maxIndex are returned as an object.
    evaluate(t: number, order: number): { minIndex: number; maxIndex: number } {
        logAssert(0 <= order && order <= 3, 'Invalid order.');

        const { index: i, t: tClamped } = this.getIndex(t);
        t = tClamped;
        this.mJet[0].set(i, 0, 1);

        if (order >= 1) {
            this.mJet[1].set(i, 0, 0);
            if (order >= 2) {
                this.mJet[2].set(i, 0, 0);
                if (order >= 3) {
                    this.mJet[3].set(i, 0, 0);
                }
            }
        }

        let n0 = t - this.mKnots[i];
        let n1 = this.mKnots[i + 1] - t;
        let e0: number, e1: number, d0: number, d1: number, invD0: number, invD1: number;
        let j: number;
        for (j = 1; j <= this.mDegree; j++) {
            d0 = this.mKnots[i + j] - this.mKnots[i];
            d1 = this.mKnots[i + 1] - this.mKnots[i - j + 1];
            invD0 = (d0 > 0 ? 1 / d0 : 0);
            invD1 = (d1 > 0 ? 1 / d1 : 0);

            e0 = n0 * this.mJet[0].get(i, j - 1);
            this.mJet[0].set(i, j, e0 * invD0);
            e1 = n1 * this.mJet[0].get(i - j + 1, j - 1);
            this.mJet[0].set(i - j, j, e1 * invD1);

            if (order >= 1) {
                e0 = n0 * this.mJet[1].get(i, j - 1) + this.mJet[0].get(i, j - 1);
                this.mJet[1].set(i, j, e0 * invD0);
                e1 = n1 * this.mJet[1].get(i - j + 1, j - 1) - this.mJet[0].get(i - j + 1, j - 1);
                this.mJet[1].set(i - j, j, e1 * invD1);

                if (order >= 2) {
                    e0 = n0 * this.mJet[2].get(i, j - 1) + 2 * this.mJet[1].get(i, j - 1);
                    this.mJet[2].set(i, j, e0 * invD0);
                    e1 = n1 * this.mJet[2].get(i - j + 1, j - 1) - 2 * this.mJet[1].get(i - j + 1, j - 1);
                    this.mJet[2].set(i - j, j, e1 * invD1);

                    if (order >= 3) {
                        e0 = n0 * this.mJet[3].get(i, j - 1) + 3 * this.mJet[2].get(i, j - 1);
                        this.mJet[3].set(i, j, e0 * invD0);
                        e1 = n1 * this.mJet[3].get(i - j + 1, j - 1) - 3 * this.mJet[2].get(i - j + 1, j - 1);
                        this.mJet[3].set(i - j, j, e1 * invD1);
                    }
                }
            }
        }

        for (j = 2; j <= this.mDegree; ++j) {
            for (let k = i - j + 1; k < i; ++k) {
                n0 = t - this.mKnots[k];
                n1 = this.mKnots[k + j + 1] - t;
                d0 = this.mKnots[k + j] - this.mKnots[k];
                d1 = this.mKnots[k + j + 1] - this.mKnots[k + 1];
                invD0 = (d0 > 0 ? 1 / d0 : 0);
                invD1 = (d1 > 0 ? 1 / d1 : 0);

                e0 = n0 * this.mJet[0].get(k, j - 1);
                e1 = n1 * this.mJet[0].get(k + 1, j - 1);
                this.mJet[0].set(k, j, e0 * invD0 + e1 * invD1);

                if (order >= 1) {
                    e0 = n0 * this.mJet[1].get(k, j - 1) + this.mJet[0].get(k, j - 1);
                    e1 = n1 * this.mJet[1].get(k + 1, j - 1) - this.mJet[0].get(k + 1, j - 1);
                    this.mJet[1].set(k, j, e0 * invD0 + e1 * invD1);

                    if (order >= 2) {
                        e0 = n0 * this.mJet[2].get(k, j - 1) + 2 * this.mJet[1].get(k, j - 1);
                        e1 = n1 * this.mJet[2].get(k + 1, j - 1) - 2 * this.mJet[1].get(k + 1, j - 1);
                        this.mJet[2].set(k, j, e0 * invD0 + e1 * invD1);

                        if (order >= 3) {
                            e0 = n0 * this.mJet[3].get(k, j - 1) + 3 * this.mJet[2].get(k, j - 1);
                            e1 = n1 * this.mJet[3].get(k + 1, j - 1) - 3 * this.mJet[2].get(k + 1, j - 1);
                            this.mJet[3].set(k, j, e0 * invD0 + e1 * invD1);
                        }
                    }
                }
            }
        }

        return { minIndex: i - this.mDegree, maxIndex: i };
    }

    // Access the results of the call to evaluate(...).  The index i must
    // satisfy minIndex <= i <= maxIndex.  If it is not, the function
    // throws.  The separation of evaluation and access is based on
    // local control of the basis function; that is, only the accessible
    // values are (potentially) not zero.
    getValue(order: number, i: number): number {
        logAssert(0 <= order && order < 4, 'Invalid order.');
        logAssert(0 <= i && i < this.mNumControls + this.mDegree, 'Invalid index.');
        return this.mJet[order].get(i, this.mDegree);
    }

    // Determine the index i for which knot[i] <= t < knot[i+1].  The
    // t-value is modified (wrapped for periodic splines, clamped for
    // nonperiodic splines); the possibly modified value is returned along
    // with the index (the C++ passes t by reference).
    private getIndex(t: number): { index: number; t: number } {
        // Find the index i for which knot[i] <= t < knot[i+1].
        if (this.mPeriodic) {
            // Wrap to [tmin,tmax].
            let r = (t - this.mTMin) % this.mTLength;
            if (r < 0) {
                r += this.mTLength;
            }
            t = this.mTMin + r;
        }

        // Clamp to [tmin,tmax].  For the periodic case, this handles
        // small numerical rounding errors near the domain endpoints.
        if (t <= this.mTMin) {
            return { index: this.mDegree, t: this.mTMin };
        }
        if (t >= this.mTMax) {
            return { index: this.mNumControls - 1, t: this.mTMax };
        }

        // At this point, tmin < t < tmax.
        for (const key of this.mKeys) {
            if (t < key.t) {
                return { index: key.lastIndex, t };
            }
        }

        // We should not reach this code.
        logError('Unexpected condition.');
    }
}
