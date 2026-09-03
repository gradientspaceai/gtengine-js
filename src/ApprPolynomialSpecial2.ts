// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprPolynomialSpecial2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Fit the data with a polynomial of the form
//     w = sum_{i=0}^{n-1} c[i]*x^{p[i]}
// where p[i] are distinct nonnegative powers provided by the caller. A
// least-squares fitting algorithm is used, but the input data is first
// mapped to (x,w) in [-1,1]^2 for numerical robustness.
//
// Port notes: the observation type 'std::array<Real,2>' becomes
// 'readonly number[]' (length 2), per the PORTING.md type mapping. The
// private member functions Transform and DoLeastSquares stay private
// methods because both mutate the transform state used by evaluate().

import { ApprQuery } from './ApprQuery.js';
import { GMatrix } from './GMatrix.js';
import { logAssert } from './Logger.js';
import { inverse, mulMatrix } from './Matrix.js';
import { Vector } from './Vector.js';

export class ApprPolynomialSpecial2 extends ApprQuery<readonly number[]> {
    private mDegrees: number[];
    private mParameters: number[];

    // Support for evaluation. The coefficients were generated for the
    // samples mapped to [-1,1]^2. The evaluate() function must transform x
    // to x' in [-1,1], compute w' in [-1,1], then transform w' to w.
    private mXDomain: [number, number];
    private mWDomain: [number, number];
    private mScale: [number, number];
    private mInvTwoWScale: number;

    // This array is used by evaluate() to avoid reallocation per call.
    private mXPowers: number[];

    // Initialize the model parameters to zero. The degrees must be
    // nonnegative and strictly increasing.
    constructor(degrees: readonly number[]) {
        super();
        this.mDegrees = degrees.slice();
        this.mParameters = new Array<number>(this.mDegrees.length).fill(0);

        logAssert(this.mDegrees.length > 0, 'The input array must have elements.');
        let lastDegree = -1;
        for (const degree of this.mDegrees) {
            logAssert(degree > lastDegree, 'Degrees must be increasing.');
            lastDegree = degree;
        }

        this.mXDomain = [Number.MAX_VALUE, -Number.MAX_VALUE];
        this.mWDomain = [Number.MAX_VALUE, -Number.MAX_VALUE];
        this.mScale = [0, 0];
        this.mInvTwoWScale = 0;

        // Powers of x are computed up to twice the powers when constructing
        // the fitted polynomial. Powers of x are computed up to the powers
        // for the evaluation of the fitted polynomial.
        this.mXPowers = new Array<number>(
            2 * this.mDegrees[this.mDegrees.length - 1] + 1).fill(0);
        this.mXPowers[0] = 1;
    }

    // Basic fitting algorithm. See ApprQuery.ts for the various fit(...)
    // functions that you can call.
    fitIndexed(observations: readonly (readonly number[])[],
        indices: readonly number[]): boolean {
        if (this.validIndices(observations, indices)) {
            // Transform the observations to [-1,1]^2 for numerical
            // robustness.
            const transformed = this.transform(observations, indices);

            // Fit the transformed data using a least-squares algorithm.
            return this.doLeastSquares(transformed);
        }

        this.mParameters.fill(0);
        return false;
    }

    // Get the parameters for the best fit.
    getParameters(): readonly number[] {
        return this.mParameters;
    }

    getMinimumRequired(): number {
        return this.mParameters.length;
    }

    // Compute the model error for the specified observation for the current
    // model parameters. The returned value for observation (x0,w0) is
    // |w(x0) - w0|, where w(x) is the fitted polynomial.
    error(observation: readonly number[]): number {
        const w = this.evaluate(observation[0]);
        return Math.abs(w - observation[1]);
    }

    // The port of the C++ 'dynamic_cast' guarded assignment. C++ value
    // assignment copies, so the port deep-copies the parameters.
    copyParameters(input: ApprQuery<readonly number[]>): void {
        if (input instanceof ApprPolynomialSpecial2) {
            this.mDegrees = input.mDegrees.slice();
            this.mParameters = input.mParameters.slice();
            this.mXDomain = [input.mXDomain[0], input.mXDomain[1]];
            this.mWDomain = [input.mWDomain[0], input.mWDomain[1]];
            this.mScale = [input.mScale[0], input.mScale[1]];
            this.mInvTwoWScale = input.mInvTwoWScale;
            this.mXPowers = input.mXPowers.slice();
        }
    }

    // Evaluate the polynomial. The domain interval is provided so you can
    // interpolate (x in domain) or extrapolate (x not in domain).
    getXDomain(): readonly [number, number] {
        return this.mXDomain;
    }

    evaluate(x: number): number {
        // Transform x to x' in [-1,1].
        x = -1 + 2 * this.mScale[0] * (x - this.mXDomain[0]);

        // Compute relevant powers of x.
        const jmax = this.mDegrees[this.mDegrees.length - 1];
        for (let j = 1, jm1 = 0; j <= jmax; ++j, ++jm1) {
            this.mXPowers[j] = this.mXPowers[jm1] * x;
        }

        let w = 0;
        const isup = this.mDegrees.length;
        for (let i = 0; i < isup; ++i) {
            const xp = this.mXPowers[this.mDegrees[i]];
            w += this.mParameters[i] * xp;
        }

        // Transform w from [-1,1] back to the original space.
        w = (w + 1) * this.mInvTwoWScale + this.mWDomain[0];
        return w;
    }

    // Transform the (x,w) values to (x',w') in [-1,1]^2.
    private transform(observations: readonly (readonly number[])[],
        indices: readonly number[]): number[][] {
        const numSamples = indices.length;
        const transformed: number[][] = new Array<number[]>(numSamples);

        const omin = observations[indices[0]].slice();
        const omax = omin.slice();
        for (let s = 1; s < numSamples; ++s) {
            const obs = observations[indices[s]];
            for (let i = 0; i < 2; ++i) {
                if (obs[i] < omin[i]) {
                    omin[i] = obs[i];
                }
                else if (obs[i] > omax[i]) {
                    omax[i] = obs[i];
                }
            }
        }

        this.mXDomain[0] = omin[0];
        this.mXDomain[1] = omax[0];
        this.mWDomain[0] = omin[1];
        this.mWDomain[1] = omax[1];
        for (let i = 0; i < 2; ++i) {
            this.mScale[i] = 1 / (omax[i] - omin[i]);
        }

        for (let s = 0; s < numSamples; ++s) {
            const obs = observations[indices[s]];
            const t: number[] = [0, 0];
            for (let i = 0; i < 2; ++i) {
                t[i] = -1 + 2 * this.mScale[i] * (obs[i] - omin[i]);
            }
            transformed[s] = t;
        }
        this.mInvTwoWScale = 0.5 / this.mScale[1];
        return transformed;
    }

    // The least-squares fitting algorithm for the transformed data.
    private doLeastSquares(transformed: readonly (readonly number[])[]): boolean {
        // Set up a linear system A*X = B, where X are the polynomial
        // coefficients.
        const size = this.mDegrees.length;
        const A = new GMatrix(size, size);
        const B = new Vector(size);

        const numSamples = transformed.length;
        const twoMaxXDegree = 2 * this.mDegrees[this.mDegrees.length - 1];
        for (let i = 0; i < numSamples; ++i) {
            // Compute relevant powers of x.
            const x = transformed[i][0];
            const w = transformed[i][1];
            for (let j = 1, jm1 = 0; j <= twoMaxXDegree; ++j, ++jm1) {
                this.mXPowers[j] = this.mXPowers[jm1] * x;
            }

            for (let row = 0; row < size; ++row) {
                // Update the upper-triangular portion of the symmetric
                // matrix.
                for (let col = row; col < size; ++col) {
                    A.set(row, col, A.get(row, col) +
                        this.mXPowers[this.mDegrees[row] + this.mDegrees[col]]);
                }

                // Update the right-hand side of the system.
                B.values[row] += this.mXPowers[this.mDegrees[row]] * w;
            }
        }

        // Copy the upper-triangular portion of the symmetric matrix to the
        // lower-triangular portion.
        for (let row = 0; row < size; ++row) {
            for (let col = 0; col < row; ++col) {
                A.set(row, col, A.get(col, row));
            }
        }

        // Precondition by normalizing the sums.
        const invNumSamples = 1 / numSamples;
        for (let i = 0; i < A.numElements; ++i) {
            A.values[i] *= invNumSamples;
        }
        for (let i = 0; i < size; ++i) {
            B.values[i] *= invNumSamples;
        }

        // Solve for the polynomial coefficients.
        const coefficients = mulMatrix(inverse(A).inverse, B) as Vector;
        let hasNonzero = false;
        for (let i = 0; i < size; ++i) {
            this.mParameters[i] = coefficients.values[i];
            if (coefficients.values[i] !== 0) {
                hasNonzero = true;
            }
        }
        return hasNonzero;
    }
}
