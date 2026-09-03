// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprHeightLine2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Least-squares fit of a line to height data (x,f(x)). The line is of the
// form (y - yAvr) = a*(x - xAvr), where (xAvr,yAvr) is the average of the
// sample points. The return value of fit is 'true' if and only if the fit is
// successful (the input points are not degenerate to a single point). The
// parameters are ((xAvr,yAvr),(a,-1)) on success and ((0,0),(0,0)) on
// failure. The error for (x0,y0) is [a*(x0-xAvr)-(y0-yAvr)]^2.
//
// Port notes: upstream stores the parameters in a
// std::pair<Vector2,Vector2>; per PORTING.md the pair becomes an object with
// named fields (average, coefficients).

import { ApprQuery } from './ApprQuery.js';
import { Vector, dot, sub } from './Vector.js';

export interface ApprHeightLine2Parameters {
    // (xAvr,yAvr), the average of the sample points.
    average: Vector;

    // (a,-1), so that the fitted line is Dot((x,y) - average, (a,-1)) = 0.
    coefficients: Vector;
}

export class ApprHeightLine2 extends ApprQuery<Vector> {
    private mParameters: ApprHeightLine2Parameters;

    // Initialize the model parameters to zero.
    constructor() {
        super();
        this.mParameters = { average: new Vector(2), coefficients: new Vector(2) };
    }

    // Basic fitting algorithm. See ApprQuery.ts for the various fit(...)
    // functions that you can call.
    fitIndexed(points: readonly Vector[], indices: readonly number[]): boolean {
        if (this.validIndices(points, indices)) {
            const numIndices = indices.length;

            // Compute the mean of the points.
            const mean = new Vector(2);
            for (let i = 0; i < numIndices; ++i) {
                const p = points[indices[i]];
                mean.values[0] += p.values[0];
                mean.values[1] += p.values[1];
            }
            mean.values[0] /= numIndices;
            mean.values[1] /= numIndices;

            // Compute the covariance matrix of the points.
            let covar00 = 0, covar01 = 0;
            for (let i = 0; i < numIndices; ++i) {
                const p = points[indices[i]];
                const d0 = p.values[0] - mean.values[0];
                const d1 = p.values[1] - mean.values[1];
                covar00 += d0 * d0;
                covar01 += d0 * d1;
            }

            // Decompose the covariance matrix.
            if (covar00 > 0) {
                const coefficients = new Vector(2);
                coefficients.values[0] = covar01 / covar00;
                coefficients.values[1] = -1;
                this.mParameters = { average: mean, coefficients: coefficients };
                return true;
            }
        }

        this.mParameters = { average: new Vector(2), coefficients: new Vector(2) };
        return false;
    }

    // Get the parameters for the best fit.
    getParameters(): ApprHeightLine2Parameters {
        return this.mParameters;
    }

    getMinimumRequired(): number {
        return 2;
    }

    error(point: Vector): number {
        const d = dot(sub(point, this.mParameters.average),
            this.mParameters.coefficients);
        return d * d;
    }

    // The port of the C++ 'dynamic_cast' guarded assignment. C++ value
    // assignment copies, so the port deep-copies the parameters.
    copyParameters(input: ApprQuery<Vector>): void {
        if (input instanceof ApprHeightLine2) {
            this.mParameters = {
                average: input.mParameters.average.clone(),
                coefficients: input.mParameters.coefficients.clone()
            };
        }
    }
}
