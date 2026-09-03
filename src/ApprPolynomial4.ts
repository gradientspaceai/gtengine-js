// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprPolynomial4.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The samples are (x[i],y[i],z[i],w[i]) for 0 <= i < S. Think of w as a
// function of x, y, and z, say w = f(x,y,z). The function fits the samples
// with a polynomial of degree d0 in x, degree d1 in y, and degree d2 in z,
// say
//   w = sum_{i=0}^{d0} sum_{j=0}^{d1} sum_{k=0}^{d2} c[i][j][k]*x^i*y^j*z^k
// The method is a least-squares fitting algorithm. getParameters() returns
// c[i][j][k] = parameters[i+(d0+1)*(j+(d1+1)*k)] for a total of
// (d0+1)*(d1+1)*(d2+1) coefficients. The observation type is a number[] of
// length 4, which represents a tuple (x,y,z,w).
//
// WARNING. The fitting algorithm for polynomial terms
//   (1,x,x^2,...,x^d0), (1,y,y^2,...,y^d1), (1,z,z^2,...,z^d2)
// is known to be nonrobust for large degrees and for large magnitude data.
// One alternative is to use orthogonal polynomials
//   (f[0](x),...,f[d0](x)), (g[0](y),...,g[d1](y)), (h[0](z),...,h[d2](z))
// and apply the least-squares algorithm to these. Another alternative is to
// transform
//   (x',y',z',w') = ((x-xcen)/rng, (y-ycen)/rng, (z-zcen)/rng, w/rng)
// where xcen, ycen and zcen are the midpoints of the data ranges and
// rng = max(xmax-xmin,ymax-ymin,zmax-zmin). Fit the transformed points and
// scale the result back to the original space.
//
// Port notes: see ApprPolynomial2.ts for the shared conventions (the
// observation type, the flattened power tables, the accumulating domains).

import { ApprQuery } from './ApprQuery.js';
import { GMatrix } from './GMatrix.js';
import { inverse, mulMatrix } from './Matrix.js';
import { Vector } from './Vector.js';

export class ApprPolynomial4 extends ApprQuery<readonly number[]> {
    private mXDegree: number;
    private mYDegree: number;
    private mZDegree: number;
    private mXDegreeP1: number;
    private mYDegreeP1: number;
    private mZDegreeP1: number;
    private mSize: number;
    private mXDomain: [number, number];
    private mYDomain: [number, number];
    private mZDomain: [number, number];
    private mParameters: number[];

    // These arrays are used by evaluate() to avoid reallocation per call.
    private mYZCoefficient: number[];
    private mZCoefficient: number[];

    // Initialize the model parameters to zero.
    constructor(xDegree: number, yDegree: number, zDegree: number) {
        super();
        this.mXDegree = xDegree;
        this.mYDegree = yDegree;
        this.mZDegree = zDegree;
        this.mXDegreeP1 = xDegree + 1;
        this.mYDegreeP1 = yDegree + 1;
        this.mZDegreeP1 = zDegree + 1;
        this.mSize = this.mXDegreeP1 * this.mYDegreeP1 * this.mZDegreeP1;
        this.mParameters = new Array<number>(this.mSize).fill(0);
        this.mYZCoefficient =
            new Array<number>(this.mYDegreeP1 * this.mZDegreeP1).fill(0);
        this.mZCoefficient = new Array<number>(this.mZDegreeP1).fill(0);
        this.mXDomain = [Number.MAX_VALUE, -Number.MAX_VALUE];
        this.mYDomain = [Number.MAX_VALUE, -Number.MAX_VALUE];
        this.mZDomain = [Number.MAX_VALUE, -Number.MAX_VALUE];
    }

    // Basic fitting algorithm. See ApprQuery.ts for the various fit(...)
    // functions that you can call.
    fitIndexed(observations: readonly (readonly number[])[],
        indices: readonly number[]): boolean {
        if (this.validIndices(observations, indices)) {
            // Compute the powers of x, y, and z.
            const numSamples = indices.length;
            const twoXDegree = 2 * this.mXDegree;
            const twoYDegree = 2 * this.mYDegree;
            const twoZDegree = 2 * this.mZDegree;
            const xStride = twoXDegree + 1;
            const yStride = twoYDegree + 1;
            const zStride = twoZDegree + 1;
            const xPower = new Array<number>(xStride * numSamples).fill(0);
            const yPower = new Array<number>(yStride * numSamples).fill(0);
            const zPower = new Array<number>(zStride * numSamples).fill(0);
            for (let s = 0; s < numSamples; ++s) {
                const x = observations[indices[s]][0];
                const y = observations[indices[s]][1];
                const z = observations[indices[s]][2];
                this.mXDomain[0] = Math.min(x, this.mXDomain[0]);
                this.mXDomain[1] = Math.max(x, this.mXDomain[1]);
                this.mYDomain[0] = Math.min(y, this.mYDomain[0]);
                this.mYDomain[1] = Math.max(y, this.mYDomain[1]);
                this.mZDomain[0] = Math.min(z, this.mZDomain[0]);
                this.mZDomain[1] = Math.max(z, this.mZDomain[1]);

                xPower[s * xStride] = 1;
                for (let i0 = 1; i0 <= twoXDegree; ++i0) {
                    xPower[s * xStride + i0] = x * xPower[s * xStride + i0 - 1];
                }

                yPower[s * yStride] = 1;
                for (let j0 = 1; j0 <= twoYDegree; ++j0) {
                    yPower[s * yStride + j0] = y * yPower[s * yStride + j0 - 1];
                }

                zPower[s * zStride] = 1;
                for (let k0 = 1; k0 <= twoZDegree; ++k0) {
                    zPower[s * zStride + k0] = z * zPower[s * zStride + k0 - 1];
                }
            }

            // Matrix A is the Vandermonde matrix and vector B is the
            // right-hand side of the linear system A*X = B.
            const A = new GMatrix(this.mSize, this.mSize);
            const B = new Vector(this.mSize);
            for (let k0 = 0; k0 <= this.mZDegree; ++k0) {
                for (let j0 = 0; j0 <= this.mYDegree; ++j0) {
                    for (let i0 = 0; i0 <= this.mXDegree; ++i0) {
                        let sum = 0;
                        const n0 = i0 + this.mXDegreeP1 *
                            (j0 + this.mYDegreeP1 * k0);
                        for (let s = 0; s < numSamples; ++s) {
                            const w = observations[indices[s]][3];
                            sum += w * xPower[s * xStride + i0] *
                                yPower[s * yStride + j0] *
                                zPower[s * zStride + k0];
                        }

                        B.values[n0] = sum;

                        for (let k1 = 0; k1 <= this.mZDegree; ++k1) {
                            for (let j1 = 0; j1 <= this.mYDegree; ++j1) {
                                for (let i1 = 0; i1 <= this.mXDegree; ++i1) {
                                    sum = 0;
                                    const n1 = i1 + this.mXDegreeP1 *
                                        (j1 + this.mYDegreeP1 * k1);
                                    for (let s = 0; s < numSamples; ++s) {
                                        sum += xPower[s * xStride + i0 + i1] *
                                            yPower[s * yStride + j0 + j1] *
                                            zPower[s * zStride + k0 + k1];
                                    }

                                    A.set(n0, n1, sum);
                                }
                            }
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
    // model parameters. The returned value for observation (x0,y0,z0,w0) is
    // |w(x0,y0,z0) - w0|, where w(x,y,z) is the fitted polynomial.
    error(observation: readonly number[]): number {
        const w = this.evaluate(observation[0], observation[1], observation[2]);
        return Math.abs(w - observation[3]);
    }

    // The port of the C++ 'dynamic_cast' guarded assignment. C++ value
    // assignment copies, so the port deep-copies the parameters.
    copyParameters(input: ApprQuery<readonly number[]>): void {
        if (input instanceof ApprPolynomial4) {
            this.mXDegree = input.mXDegree;
            this.mYDegree = input.mYDegree;
            this.mZDegree = input.mZDegree;
            this.mXDegreeP1 = input.mXDegreeP1;
            this.mYDegreeP1 = input.mYDegreeP1;
            this.mZDegreeP1 = input.mZDegreeP1;
            this.mSize = input.mSize;
            this.mXDomain = [input.mXDomain[0], input.mXDomain[1]];
            this.mYDomain = [input.mYDomain[0], input.mYDomain[1]];
            this.mZDomain = [input.mZDomain[0], input.mZDomain[1]];
            this.mParameters = input.mParameters.slice();
            this.mYZCoefficient = input.mYZCoefficient.slice();
            this.mZCoefficient = input.mZCoefficient.slice();
        }
    }

    // Evaluate the polynomial. The domain intervals are provided so you can
    // interpolate ((x,y,z) in domain) or extrapolate ((x,y,z) not in
    // domain).
    getXDomain(): readonly [number, number] {
        return this.mXDomain;
    }

    getYDomain(): readonly [number, number] {
        return this.mYDomain;
    }

    getZDomain(): readonly [number, number] {
        return this.mZDomain;
    }

    evaluate(x: number, y: number, z: number): number {
        let i0: number, i1: number, i2: number;
        let w: number;

        for (i2 = 0; i2 <= this.mZDegree; ++i2) {
            for (i1 = 0; i1 <= this.mYDegree; ++i1) {
                i0 = this.mXDegree;
                w = this.mParameters[i0 + this.mXDegreeP1 *
                    (i1 + this.mYDegreeP1 * i2)];
                while (--i0 >= 0) {
                    w = this.mParameters[i0 + this.mXDegreeP1 *
                        (i1 + this.mYDegreeP1 * i2)] + w * x;
                }
                this.mYZCoefficient[i1 + this.mYDegreeP1 * i2] = w;
            }
        }

        for (i2 = 0; i2 <= this.mZDegree; ++i2) {
            i1 = this.mYDegree;
            w = this.mYZCoefficient[i1 + this.mYDegreeP1 * i2];
            while (--i1 >= 0) {
                w = this.mYZCoefficient[i1 + this.mYDegreeP1 * i2] + w * y;
            }
            this.mZCoefficient[i2] = w;
        }

        i2 = this.mZDegree;
        w = this.mZCoefficient[i2];
        while (--i2 >= 0) {
            w = this.mZCoefficient[i2] + w * z;
        }

        return w;
    }
}
