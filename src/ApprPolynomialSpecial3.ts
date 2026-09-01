// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprPolynomialSpecial3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Fit the data with a polynomial of the form
//     w = sum_{i=0}^{n-1} c[i]*x^{p[i]}*y^{q[i]}
// where <p[i],q[i]> are distinct pairs of nonnegative powers provided by the
// caller. A least-squares fitting algorithm is used, but the input data is
// first mapped to (x,y,w) in [-1,1]^3 for numerical robustness.
//
// Port notes: see ApprPolynomialSpecial2.ts for the shared conventions.
// Upstream validates that xDegrees and yDegrees are each strictly
// increasing. That is a stronger requirement than "the <p[i],q[i]> pairs are
// distinct" promised by the header comment -- a natural term set such as
// {(0,0),(1,0),(0,1)} is rejected. The port preserves the upstream
// assertions verbatim rather than relaxing them.

import { ApprQuery } from './ApprQuery';
import { GMatrix } from './GMatrix';
import { logAssert } from './Logger';
import { inverse, mulMatrix } from './Matrix';
import { Vector } from './Vector';

export class ApprPolynomialSpecial3 extends ApprQuery<readonly number[]> {
    private mXDegrees: number[];
    private mYDegrees: number[];
    private mParameters: number[];

    // Support for evaluation. The coefficients were generated for the
    // samples mapped to [-1,1]^3. The evaluate() function must transform
    // (x,y) to (x',y') in [-1,1]^2, compute w' in [-1,1], then transform w'
    // to w.
    private mXDomain: [number, number];
    private mYDomain: [number, number];
    private mWDomain: [number, number];
    private mScale: [number, number, number];
    private mInvTwoWScale: number;

    // These arrays are used by evaluate() to avoid reallocation per call.
    private mXPowers: number[];
    private mYPowers: number[];

    // Initialize the model parameters to zero. The degrees must be
    // nonnegative and strictly increasing.
    constructor(xDegrees: readonly number[], yDegrees: readonly number[]) {
        super();
        this.mXDegrees = xDegrees.slice();
        this.mYDegrees = yDegrees.slice();
        this.mParameters = new Array<number>(this.mXDegrees.length).fill(0);

        logAssert(this.mXDegrees.length === this.mYDegrees.length,
            'The input arrays must have the same size.');

        logAssert(this.mXDegrees.length > 0, 'The input array must have elements.');
        let lastDegree = -1;
        for (const degree of this.mXDegrees) {
            logAssert(degree > lastDegree, 'Degrees must be increasing.');
            lastDegree = degree;
        }

        logAssert(this.mYDegrees.length > 0, 'The input array must have elements.');
        lastDegree = -1;
        for (const degree of this.mYDegrees) {
            logAssert(degree > lastDegree, 'Degrees must be increasing.');
            lastDegree = degree;
        }

        this.mXDomain = [Number.MAX_VALUE, -Number.MAX_VALUE];
        this.mYDomain = [Number.MAX_VALUE, -Number.MAX_VALUE];
        this.mWDomain = [Number.MAX_VALUE, -Number.MAX_VALUE];
        this.mScale = [0, 0, 0];
        this.mInvTwoWScale = 0;

        // Powers of x and y are computed up to twice the powers when
        // constructing the fitted polynomial. Powers of x and y are computed
        // up to the powers for the evaluation of the fitted polynomial.
        this.mXPowers = new Array<number>(
            2 * this.mXDegrees[this.mXDegrees.length - 1] + 1).fill(0);
        this.mXPowers[0] = 1;
        this.mYPowers = new Array<number>(
            2 * this.mYDegrees[this.mYDegrees.length - 1] + 1).fill(0);
        this.mYPowers[0] = 1;
    }

    // Basic fitting algorithm. See ApprQuery.ts for the various fit(...)
    // functions that you can call.
    fitIndexed(observations: readonly (readonly number[])[],
        indices: readonly number[]): boolean {
        if (this.validIndices(observations, indices)) {
            // Transform the observations to [-1,1]^3 for numerical
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
    // model parameters. The returned value for observation (x0,y0,w0) is
    // |w(x0,y0) - w0|, where w(x,y) is the fitted polynomial.
    error(observation: readonly number[]): number {
        const w = this.evaluate(observation[0], observation[1]);
        return Math.abs(w - observation[2]);
    }

    // The port of the C++ 'dynamic_cast' guarded assignment. C++ value
    // assignment copies, so the port deep-copies the parameters.
    copyParameters(input: ApprQuery<readonly number[]>): void {
        if (input instanceof ApprPolynomialSpecial3) {
            this.mXDegrees = input.mXDegrees.slice();
            this.mYDegrees = input.mYDegrees.slice();
            this.mParameters = input.mParameters.slice();
            this.mXDomain = [input.mXDomain[0], input.mXDomain[1]];
            this.mYDomain = [input.mYDomain[0], input.mYDomain[1]];
            this.mWDomain = [input.mWDomain[0], input.mWDomain[1]];
            this.mScale = [input.mScale[0], input.mScale[1], input.mScale[2]];
            this.mInvTwoWScale = input.mInvTwoWScale;
            this.mXPowers = input.mXPowers.slice();
            this.mYPowers = input.mYPowers.slice();
        }
    }

    // Evaluate the polynomial. The domain intervals are provided so you can
    // interpolate ((x,y) in domain) or extrapolate ((x,y) not in domain).
    getXDomain(): readonly [number, number] {
        return this.mXDomain;
    }

    getYDomain(): readonly [number, number] {
        return this.mYDomain;
    }

    evaluate(x: number, y: number): number {
        // Transform (x,y) to (x',y') in [-1,1]^2.
        x = -1 + 2 * this.mScale[0] * (x - this.mXDomain[0]);
        y = -1 + 2 * this.mScale[1] * (y - this.mYDomain[0]);

        // Compute relevant powers of x and y.
        let jmax = this.mXDegrees[this.mXDegrees.length - 1];
        for (let j = 1, jm1 = 0; j <= jmax; ++j, ++jm1) {
            this.mXPowers[j] = this.mXPowers[jm1] * x;
        }

        jmax = this.mYDegrees[this.mYDegrees.length - 1];
        for (let j = 1, jm1 = 0; j <= jmax; ++j, ++jm1) {
            this.mYPowers[j] = this.mYPowers[jm1] * y;
        }

        let w = 0;
        const isup = this.mXDegrees.length;
        for (let i = 0; i < isup; ++i) {
            const xp = this.mXPowers[this.mXDegrees[i]];
            const yp = this.mYPowers[this.mYDegrees[i]];
            w += this.mParameters[i] * xp * yp;
        }

        // Transform w from [-1,1] back to the original space.
        w = (w + 1) * this.mInvTwoWScale + this.mWDomain[0];
        return w;
    }

    // Transform the (x,y,w) values to (x',y',w') in [-1,1]^3.
    private transform(observations: readonly (readonly number[])[],
        indices: readonly number[]): number[][] {
        const numSamples = indices.length;
        const transformed: number[][] = new Array<number[]>(numSamples);

        const omin = observations[indices[0]].slice();
        const omax = omin.slice();
        for (let s = 1; s < numSamples; ++s) {
            const obs = observations[indices[s]];
            for (let i = 0; i < 3; ++i) {
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
        this.mYDomain[0] = omin[1];
        this.mYDomain[1] = omax[1];
        this.mWDomain[0] = omin[2];
        this.mWDomain[1] = omax[2];
        for (let i = 0; i < 3; ++i) {
            this.mScale[i] = 1 / (omax[i] - omin[i]);
        }

        for (let s = 0; s < numSamples; ++s) {
            const obs = observations[indices[s]];
            const t: number[] = [0, 0, 0];
            for (let i = 0; i < 3; ++i) {
                t[i] = -1 + 2 * this.mScale[i] * (obs[i] - omin[i]);
            }
            transformed[s] = t;
        }
        this.mInvTwoWScale = 0.5 / this.mScale[2];
        return transformed;
    }

    // The least-squares fitting algorithm for the transformed data.
    private doLeastSquares(transformed: readonly (readonly number[])[]): boolean {
        // Set up a linear system A*X = B, where X are the polynomial
        // coefficients.
        const size = this.mXDegrees.length;
        const A = new GMatrix(size, size);
        const B = new Vector(size);

        const numSamples = transformed.length;
        const twoMaxXDegree = 2 * this.mXDegrees[this.mXDegrees.length - 1];
        const twoMaxYDegree = 2 * this.mYDegrees[this.mYDegrees.length - 1];
        for (let i = 0; i < numSamples; ++i) {
            // Compute relevant powers of x and y.
            const x = transformed[i][0];
            const y = transformed[i][1];
            const w = transformed[i][2];
            for (let j = 1, jm1 = 0; j <= twoMaxXDegree; ++j, ++jm1) {
                this.mXPowers[j] = this.mXPowers[jm1] * x;
            }
            for (let j = 1, jm1 = 0; j <= twoMaxYDegree; ++j, ++jm1) {
                this.mYPowers[j] = this.mYPowers[jm1] * y;
            }

            for (let row = 0; row < size; ++row) {
                // Update the upper-triangular portion of the symmetric
                // matrix.
                let xp: number, yp: number;
                for (let col = row; col < size; ++col) {
                    xp = this.mXPowers[this.mXDegrees[row] + this.mXDegrees[col]];
                    yp = this.mYPowers[this.mYDegrees[row] + this.mYDegrees[col]];
                    A.set(row, col, A.get(row, col) + xp * yp);
                }

                // Update the right-hand side of the system.
                xp = this.mXPowers[this.mXDegrees[row]];
                yp = this.mYPowers[this.mYDegrees[row]];
                B.values[row] += xp * yp * w;
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
