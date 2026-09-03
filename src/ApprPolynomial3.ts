// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprPolynomial3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The samples are (x[i],y[i],w[i]) for 0 <= i < S. Think of w as a function
// of x and y, say w = f(x,y). The function fits the samples with a
// polynomial of degree d0 in x and degree d1 in y, say
//   w = sum_{i=0}^{d0} sum_{j=0}^{d1} c[i][j]*x^i*y^j
// The method is a least-squares fitting algorithm. getParameters() returns
// c[i][j] = parameters[i+(d0+1)*j] for a total of (d0+1)*(d1+1)
// coefficients. The observation type is a number[] of length 3, which
// represents a triple (x,y,w).
//
// WARNING. The fitting algorithm for polynomial terms
//   (1,x,x^2,...,x^d0), (1,y,y^2,...,y^d1)
// is known to be nonrobust for large degrees and for large magnitude data.
// One alternative is to use orthogonal polynomials
//   (f[0](x),...,f[d0](x)), (g[0](y),...,g[d1](y))
// and apply the least-squares algorithm to these. Another alternative is to
// transform
//   (x',y',w') = ((x-xcen)/rng, (y-ycen)/rng, w/rng)
// where xmin = min(x[i]), xmax = max(x[i]), xcen = (xmin+xmax)/2,
// ymin = min(y[i]), ymax = max(y[i]), ycen = (ymin+ymax)/2, and
// rng = max(xmax-xmin,ymax-ymin). Fit the (x',y',w') points,
//   w' = sum_{i=0}^{d0} sum_{j=0}^{d1} c'[i][j]*(x')^i*(y')^j
// The original polynomial is evaluated as
//   w = rng * sum_{i=0}^{d0} sum_{j=0}^{d1} c'[i][j] *
//       ((x-xcen)/rng)^i * ((y-ycen)/rng)^j
//
// Port notes: see ApprPolynomial2.ts for the shared conventions (the
// observation type, the flattened power tables, the accumulating domains).

import { ApprQuery } from './ApprQuery.js';
import { GMatrix } from './GMatrix.js';
import { inverse, mulMatrix } from './Matrix.js';
import { Vector } from './Vector.js';

export class ApprPolynomial3 extends ApprQuery<readonly number[]> {
    private mXDegree: number;
    private mYDegree: number;
    private mXDegreeP1: number;
    private mYDegreeP1: number;
    private mSize: number;
    private mXDomain: [number, number];
    private mYDomain: [number, number];
    private mParameters: number[];

    // This array is used by evaluate() to avoid reallocation for each call.
    private mYCoefficient: number[];

    // Initialize the model parameters to zero.
    constructor(xDegree: number, yDegree: number) {
        super();
        this.mXDegree = xDegree;
        this.mYDegree = yDegree;
        this.mXDegreeP1 = xDegree + 1;
        this.mYDegreeP1 = yDegree + 1;
        this.mSize = this.mXDegreeP1 * this.mYDegreeP1;
        this.mParameters = new Array<number>(this.mSize).fill(0);
        this.mYCoefficient = new Array<number>(this.mYDegreeP1).fill(0);
        this.mXDomain = [Number.MAX_VALUE, -Number.MAX_VALUE];
        this.mYDomain = [Number.MAX_VALUE, -Number.MAX_VALUE];
    }

    // Basic fitting algorithm. See ApprQuery.ts for the various fit(...)
    // functions that you can call.
    fitIndexed(observations: readonly (readonly number[])[],
        indices: readonly number[]): boolean {
        if (this.validIndices(observations, indices)) {
            // Compute the powers of x and y.
            const numSamples = indices.length;
            const twoXDegree = 2 * this.mXDegree;
            const twoYDegree = 2 * this.mYDegree;
            const xStride = twoXDegree + 1;
            const yStride = twoYDegree + 1;
            const xPower = new Array<number>(xStride * numSamples).fill(0);
            const yPower = new Array<number>(yStride * numSamples).fill(0);
            for (let s = 0; s < numSamples; ++s) {
                const x = observations[indices[s]][0];
                const y = observations[indices[s]][1];
                this.mXDomain[0] = Math.min(x, this.mXDomain[0]);
                this.mXDomain[1] = Math.max(x, this.mXDomain[1]);
                this.mYDomain[0] = Math.min(y, this.mYDomain[0]);
                this.mYDomain[1] = Math.max(y, this.mYDomain[1]);

                xPower[s * xStride] = 1;
                for (let i0 = 1; i0 <= twoXDegree; ++i0) {
                    xPower[s * xStride + i0] = x * xPower[s * xStride + i0 - 1];
                }

                yPower[s * yStride] = 1;
                for (let j0 = 1; j0 <= twoYDegree; ++j0) {
                    yPower[s * yStride + j0] = y * yPower[s * yStride + j0 - 1];
                }
            }

            // Matrix A is the Vandermonde matrix and vector B is the
            // right-hand side of the linear system A*X = B.
            const A = new GMatrix(this.mSize, this.mSize);
            const B = new Vector(this.mSize);
            for (let j0 = 0; j0 <= this.mYDegree; ++j0) {
                for (let i0 = 0; i0 <= this.mXDegree; ++i0) {
                    let sum = 0;
                    const k0 = i0 + this.mXDegreeP1 * j0;
                    for (let s = 0; s < numSamples; ++s) {
                        const w = observations[indices[s]][2];
                        sum += w * xPower[s * xStride + i0] * yPower[s * yStride + j0];
                    }

                    B.values[k0] = sum;

                    for (let j1 = 0; j1 <= this.mYDegree; ++j1) {
                        for (let i1 = 0; i1 <= this.mXDegree; ++i1) {
                            sum = 0;
                            const k1 = i1 + this.mXDegreeP1 * j1;
                            for (let s = 0; s < numSamples; ++s) {
                                sum += xPower[s * xStride + i0 + i1] *
                                    yPower[s * yStride + j0 + j1];
                            }

                            A.set(k0, k1, sum);
                        }
                    }
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
    // model parameters. The returned value for observation (x0,y0,w0) is
    // |w(x0,y0) - w0|, where w(x,y) is the fitted polynomial.
    error(observation: readonly number[]): number {
        const w = this.evaluate(observation[0], observation[1]);
        return Math.abs(w - observation[2]);
    }

    // The port of the C++ 'dynamic_cast' guarded assignment. C++ value
    // assignment copies, so the port deep-copies the parameters.
    copyParameters(input: ApprQuery<readonly number[]>): void {
        if (input instanceof ApprPolynomial3) {
            this.mXDegree = input.mXDegree;
            this.mYDegree = input.mYDegree;
            this.mXDegreeP1 = input.mXDegreeP1;
            this.mYDegreeP1 = input.mYDegreeP1;
            this.mSize = input.mSize;
            this.mXDomain = [input.mXDomain[0], input.mXDomain[1]];
            this.mYDomain = [input.mYDomain[0], input.mYDomain[1]];
            this.mParameters = input.mParameters.slice();
            this.mYCoefficient = input.mYCoefficient.slice();
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
        let i0: number, i1: number;
        let w: number;

        for (i1 = 0; i1 <= this.mYDegree; ++i1) {
            i0 = this.mXDegree;
            w = this.mParameters[i0 + this.mXDegreeP1 * i1];
            while (--i0 >= 0) {
                w = this.mParameters[i0 + this.mXDegreeP1 * i1] + w * x;
            }
            this.mYCoefficient[i1] = w;
        }

        i1 = this.mYDegree;
        w = this.mYCoefficient[i1];
        while (--i1 >= 0) {
            w = this.mYCoefficient[i1] + w * y;
        }

        return w;
    }
}
