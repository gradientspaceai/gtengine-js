// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MinimizeN.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Cartesian-product domain provided to getMinimum(*) has minimum values
// stored in t0[0..d-1] and maximum values stored in t1[0..d-1], where d is
// 'dimensions'. The domain is searched along lines through the current
// estimate of the minimum location. Each such line is searched for a minimum
// using a Minimize1 object. This is called "Powell's Direction Set Method".
// The parameters 'maxLevel' and 'maxBracket' are used by Minimize1, so read
// the documentation for that class (in its source file) to understand what
// these mean. The input 'maxIterations' is the number of iterations for the
// direction-set method.

import { GVector } from './GVector.js';
import { Minimize1 } from './Minimize1.js';
import { Vector, add, div, length as vectorLength, mul, sub } from './Vector.js';

// The port of 'std::function<Real(Real const*)>'. The domain point is passed
// as an array of 'dimensions' components. Follows the Minimize1Function
// precedent of Minimize1.ts.
export type MinimizeNFunction = (t: readonly number[]) => number;

// Upstream returns the location and value of the minimum through the output
// parameters tMin and fMin; per PORTING.md these become the named fields of
// a returned object.
export interface MinimizeNResult {
    // The location of the estimated minimum, tMin[0..d-1].
    tMin: number[];

    // The value of the function at the estimated minimum.
    fMin: number;
}

export class MinimizeN {
    private mDimensions: number;
    private mFunction: MinimizeNFunction;
    private mMaxIterations: number;
    private mEpsilon: number;
    private mDirections: GVector[];
    private mDConjIndex: number;
    private mDCurrIndex: number;
    private mTCurr: GVector;
    private mTSave: GVector;
    private mFCurr: number;
    private mMinimizer: Minimize1;

    // Construction.
    constructor(dimensions: number, F: MinimizeNFunction, maxLevel: number,
        maxBracket: number, maxIterations: number, epsilon: number = 1e-06) {
        this.mDimensions = dimensions;
        this.mFunction = F;
        this.mMaxIterations = maxIterations;
        this.mEpsilon = 0;
        this.mDConjIndex = dimensions;
        this.mDCurrIndex = 0;
        this.mTCurr = new GVector(dimensions);
        this.mTSave = new GVector(dimensions);
        this.mFCurr = 0;

        // The direction set has 'dimensions'+1 entries; the last one stores
        // the conjugate direction estimated by an iteration.
        this.mDirections = new Array<GVector>(dimensions + 1);
        for (let i = 0; i <= dimensions; ++i) {
            this.mDirections[i] = new GVector(dimensions);
        }

        // The 1-dimensional minimizer searches F along the line
        // mTCurr + t * mDirections[mDCurrIndex].
        this.mMinimizer = new Minimize1((t: number): number => {
            const point = add(this.mTCurr, mul(this.mDirections[this.mDCurrIndex], t));
            return this.mFunction(point.values);
        }, maxLevel, maxBracket);

        this.setEpsilon(epsilon);
    }

    // Member access.
    setEpsilon(epsilon: number): void {
        this.mEpsilon = (epsilon > 0 ? epsilon : 0);
    }

    getEpsilon(): number {
        return this.mEpsilon;
    }

    // Find the minimum on the Cartesian-product domain whose minimum values
    // are stored in t0[0..d-1] and whose maximum values are stored in
    // t1[0..d-1], where d is 'dimensions'. An initial guess is specified in
    // tInitial[0..d-1]. The location of the minimum is the returned tMin
    // and the value of the minimum is the returned fMin.
    getMinimum(t0: readonly number[], t1: readonly number[],
        tInitial: readonly number[]): MinimizeNResult {
        // The initial guess.
        this.mFCurr = this.mFunction(tInitial);
        for (let i = 0; i < this.mDimensions; ++i) {
            this.mTSave.values[i] = tInitial[i];
            this.mTCurr.values[i] = tInitial[i];
        }

        // Initialize the direction set to the standard Euclidean basis.
        for (let i = 0; i < this.mDimensions; ++i) {
            this.mDirections[i].makeUnit(i);
        }

        const domain = { ell0: 0, ell1: 0 };
        for (let iter = 0; iter < this.mMaxIterations; ++iter) {
            // Find minimum in each direction and update current location.
            for (let i = 0; i < this.mDimensions; ++i) {
                this.mDCurrIndex = i;
                this.computeDomain(t0, t1, domain);
                const result = this.mMinimizer.getMinimum(domain.ell0, domain.ell1, 0);
                this.mFCurr = result.fMin;
                this.mTCurr = toGVector(
                    add(this.mTCurr, mul(this.mDirections[i], result.tMin)));
            }

            // Estimate a unit-length conjugate direction.
            this.mDirections[this.mDConjIndex] =
                toGVector(sub(this.mTCurr, this.mTSave));
            const length = vectorLength(this.mDirections[this.mDConjIndex]);
            if (length <= this.mEpsilon) {
                // New position did not change significantly from old one.
                // Should there be a better convergence criterion here?
                break;
            }

            this.mDirections[this.mDConjIndex] =
                toGVector(div(this.mDirections[this.mDConjIndex], length));

            // Minimize in conjugate direction.
            this.mDCurrIndex = this.mDConjIndex;
            this.computeDomain(t0, t1, domain);
            const result = this.mMinimizer.getMinimum(domain.ell0, domain.ell1, 0);
            this.mFCurr = result.fMin;
            this.mTCurr = toGVector(
                add(this.mTCurr, mul(this.mDirections[this.mDCurrIndex], result.tMin)));

            // Cycle the directions and add conjugate direction to set. The
            // shift moves the conjugate direction from slot mDimensions into
            // slot mDimensions-1 and drops the oldest direction, so slot
            // mDimensions is the scratch slot for the next pass.
            //
            // UPSTREAM BUG (fixed here): MinimizeN.h executes
            // 'mDConjIndex = 0;' at this point. From the second iteration
            // onward that writes the conjugate direction into mDirections[0],
            // where the shift below immediately overwrites it with
            // mDirections[1], while mDirections[mDimensions] is never updated
            // again and keeps the first iteration's conjugate direction. The
            // direction set then degenerates (mDirections[mDimensions-1] and
            // [mDimensions-2] become equal) and Powell's method stalls. For
            // Rosenbrock's function from (-1.2,1) on [-2,2]^2 the upstream
            // code stops at f = 3.0548 at (-0.7436, 0.5650) no matter how many
            // iterations are allowed; without the stray assignment the same
            // run converges to f = 8e-12 at (1.0000, 1.0000). mDConjIndex is
            // initialized to mDimensions by the constructor, which is the
            // value the algorithm requires, so the port simply leaves it
            // alone.
            for (let i = 0, ip1 = 1; i < this.mDimensions; ++i, ++ip1) {
                // C++ vector assignment copies; clone so the entries do not
                // alias one another.
                this.mDirections[i] = cloneGVector(this.mDirections[ip1]);
            }

            // Set parameters for next pass.
            this.mTSave = cloneGVector(this.mTCurr);
        }

        return { tMin: this.mTCurr.values.slice(0, this.mDimensions), fMin: this.mFCurr };
    }

    // The current estimate of the minimum location is mTCurr[0..d-1]. The
    // direction of the current line to search is mDirections[mDCurrIndex].
    // This line must be clipped against the Cartesian-product domain, a
    // process implemented in this function. If the line is mTCurr+s*mDCurr,
    // the clip result is the s-interval [ell0,ell1].
    private computeDomain(t0: readonly number[], t1: readonly number[],
        domain: { ell0: number; ell1: number }): void {
        domain.ell0 = -Number.MAX_VALUE;
        domain.ell1 = +Number.MAX_VALUE;

        const direction = this.mDirections[this.mDCurrIndex];
        for (let i = 0; i < this.mDimensions; ++i) {
            const value = direction.values[i];
            if (value !== 0) {
                let b0 = t0[i] - this.mTCurr.values[i];
                let b1 = t1[i] - this.mTCurr.values[i];
                const inv = 1 / value;
                if (value > 0) {
                    // The valid t-interval is [b0,b1].
                    b0 *= inv;
                    if (b0 > domain.ell0) {
                        domain.ell0 = b0;
                    }
                    b1 *= inv;
                    if (b1 < domain.ell1) {
                        domain.ell1 = b1;
                    }
                }
                else {
                    // The valid t-interval is [b1,b0].
                    b0 *= inv;
                    if (b0 < domain.ell1) {
                        domain.ell1 = b0;
                    }
                    b1 *= inv;
                    if (b1 > domain.ell0) {
                        domain.ell0 = b1;
                    }
                }
            }
        }

        // Correction if numerical errors lead to values nearly zero.
        if (domain.ell0 > 0) {
            domain.ell0 = 0;
        }
        if (domain.ell1 < 0) {
            domain.ell1 = 0;
        }
    }
}

// The Vector free functions return Vector objects; MinimizeN stores GVector
// objects (matching upstream), so wrap the results.
function toGVector(v: Vector): GVector {
    return GVector.fromArray(v.values);
}

function cloneGVector(v: GVector): GVector {
    return GVector.fromArray(v.values);
}
