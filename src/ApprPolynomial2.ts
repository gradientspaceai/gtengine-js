// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprPolynomial2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The samples are (x[i],w[i]) for 0 <= i < S. Think of w as a function of
// x, say w = f(x). The function fits the samples with a polynomial of
// degree d, say w = sum_{i=0}^d c[i]*x^i. The method is a least-squares
// fitting algorithm. getParameters() returns the coefficients c[i] for
// 0 <= i <= d. The observation type is a number[] of length 2, which
// represents a pair (x,w).
//
// WARNING. The fitting algorithm for polynomial terms
//   (1,x,x^2,...,x^d)
// is known to be nonrobust for large degrees and for large magnitude data.
// One alternative is to use orthogonal polynomials
//   (f[0](x),...,f[d](x))
// and apply the least-squares algorithm to these. Another alternative is to
// transform
//   (x',w') = ((x-xcen)/rng, w/rng)
// where xmin = min(x[i]), xmax = max(x[i]), xcen = (xmin+xmax)/2, and
// rng = xmax-xmin. Fit the (x',w') points,
//   w' = sum_{i=0}^d c'[i]*(x')^i.
// The original polynomial is evaluated as
//   w = rng*sum_{i=0}^d c'[i]*((x-xcen)/rng)^i
//
// Port notes: the observation type 'std::array<Real,2>' becomes
// 'readonly number[]' (length 2), per the PORTING.md type mapping. The
// upstream 'Array2<Real> xPower' scratch table becomes a flat number[] with
// the same index arithmetic. The x-domain is an accumulating [min,max] pair
// that upstream never resets between fits; the port preserves that.

import { ApprQuery } from './ApprQuery';
import { GMatrix } from './GMatrix';
import { inverse, mulMatrix } from './Matrix';
import { Vector } from './Vector';

export class ApprPolynomial2 extends ApprQuery<readonly number[]> {
    private mDegree: number;
    private mSize: number;
    private mXDomain: [number, number];
    private mParameters: number[];

    // Initialize the model parameters to zero.
    constructor(degree: number) {
        super();
        this.mDegree = degree;
        this.mSize = degree + 1;
        this.mParameters = new Array<number>(this.mSize).fill(0);
        this.mXDomain = [Number.MAX_VALUE, -Number.MAX_VALUE];
    }

    // Basic fitting algorithm. See ApprQuery.ts for the various fit(...)
    // functions that you can call.
    fitIndexed(observations: readonly (readonly number[])[],
        indices: readonly number[]): boolean {
        if (this.validIndices(observations, indices)) {
            // Compute the powers of x.
            const numSamples = indices.length;
            const twoDegree = 2 * this.mDegree;
            const stride = twoDegree + 1;
            const xPower = new Array<number>(stride * numSamples).fill(0);
            for (let s = 0; s < numSamples; ++s) {
                const x = observations[indices[s]][0];
                this.mXDomain[0] = Math.min(x, this.mXDomain[0]);
                this.mXDomain[1] = Math.max(x, this.mXDomain[1]);

                xPower[s * stride] = 1;
                for (let i0 = 1; i0 <= twoDegree; ++i0) {
                    xPower[s * stride + i0] = x * xPower[s * stride + i0 - 1];
                }
            }

            // Matrix A is the Vandermonde matrix and vector B is the
            // right-hand side of the linear system A*X = B.
            const A = new GMatrix(this.mSize, this.mSize);
            const B = new Vector(this.mSize);
            for (let i0 = 0; i0 <= this.mDegree; ++i0) {
                let sum = 0;
                for (let s = 0; s < numSamples; ++s) {
                    const w = observations[indices[s]][1];
                    sum += w * xPower[s * stride + i0];
                }

                B.values[i0] = sum;

                for (let i1 = 0; i1 <= this.mDegree; ++i1) {
                    sum = 0;
                    for (let s = 0; s < numSamples; ++s) {
                        sum += xPower[s * stride + i0 + i1];
                    }

                    A.set(i0, i1, sum);
                }
            }

            // Solve for the polynomial coefficients.
            const coefficients = mulMatrix(inverse(A).inverse, B) as Vector;
            let hasNonzero = false;
            for (let i = 0; i < this.mSize; ++i) {
                this.mParameters[i] = coefficients.values[i];
                if (coefficients.values[i] !== 0) {
                    hasNonzero = true;
                }
            }
            return hasNonzero;
        }

        this.mParameters.fill(0);
        return false;
    }

    // Get the parameters for the best fit.
    getParameters(): readonly number[] {
        return this.mParameters;
    }

    getMinimumRequired(): number {
        return this.mSize;
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
        if (input instanceof ApprPolynomial2) {
            this.mDegree = input.mDegree;
            this.mSize = input.mSize;
            this.mXDomain = [input.mXDomain[0], input.mXDomain[1]];
            this.mParameters = input.mParameters.slice();
        }
    }

    // Evaluate the polynomial. The domain interval is provided so you can
    // interpolate (x in domain) or extrapolate (x not in domain).
    getXDomain(): readonly [number, number] {
        return this.mXDomain;
    }

    evaluate(x: number): number {
        let i = this.mDegree;
        let w = this.mParameters[i];
        while (--i >= 0) {
            w = this.mParameters[i] + w * x;
        }
        return w;
    }
}
