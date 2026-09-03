// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprHeightPlane3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Least-squares fit of a plane to height data (x,y,f(x,y)). The plane is of
// the form (z - zAvr) = a*(x - xAvr) + b*(y - yAvr), where (xAvr,yAvr,zAvr)
// is the average of the sample points. The return value is 'true' if and
// only if the fit is successful (the input points are noncollinear). The
// parameters are ((xAvr,yAvr,zAvr),(a,b,-1)) on success and ((0,0,0),(0,0,0))
// on failure. The error for (x0,y0,z0) is
// [a*(x0-xAvr)+b*(y0-yAvr)-(z0-zAvr)]^2.
//
// Port notes: upstream stores the parameters in a
// std::pair<Vector3,Vector3>; per PORTING.md the pair becomes an object with
// named fields (average, coefficients).

import { ApprQuery } from './ApprQuery.js';
import { Vector, dot, sub } from './Vector.js';

export interface ApprHeightPlane3Parameters {
    // (xAvr,yAvr,zAvr), the average of the sample points.
    average: Vector;

    // (a,b,-1), so that the fitted plane is
    // Dot((x,y,z) - average, (a,b,-1)) = 0.
    coefficients: Vector;
}

export class ApprHeightPlane3 extends ApprQuery<Vector> {
    private mParameters: ApprHeightPlane3Parameters;

    // Initialize the model parameters to zero.
    constructor() {
        super();
        this.mParameters = { average: new Vector(3), coefficients: new Vector(3) };
    }

    // Basic fitting algorithm. See ApprQuery.ts for the various fit(...)
    // functions that you can call.
    fitIndexed(points: readonly Vector[], indices: readonly number[]): boolean {
        if (this.validIndices(points, indices)) {
            const numIndices = indices.length;

            // Compute the mean of the points.
            const mean = new Vector(3);
            for (let i = 0; i < numIndices; ++i) {
                const p = points[indices[i]];
                for (let d = 0; d < 3; ++d) {
                    mean.values[d] += p.values[d];
                }
            }
            for (let d = 0; d < 3; ++d) {
                mean.values[d] /= numIndices;
            }

            // Compute the covariance matrix of the points.
            let covar00 = 0, covar01 = 0, covar02 = 0;
            let covar11 = 0, covar12 = 0;
            for (let i = 0; i < numIndices; ++i) {
                const p = points[indices[i]];
                const d0 = p.values[0] - mean.values[0];
                const d1 = p.values[1] - mean.values[1];
                const d2 = p.values[2] - mean.values[2];
                covar00 += d0 * d0;
                covar01 += d0 * d1;
                covar02 += d0 * d2;
                covar11 += d1 * d1;
                covar12 += d1 * d2;
            }

            // Decompose the covariance matrix.
            const det = covar00 * covar11 - covar01 * covar01;
            if (det !== 0) {
                const invDet = 1 / det;
                const coefficients = new Vector(3);
                coefficients.values[0] = (covar11 * covar02 - covar01 * covar12) * invDet;
                coefficients.values[1] = (covar00 * covar12 - covar01 * covar02) * invDet;
                coefficients.values[2] = -1;
                this.mParameters = { average: mean, coefficients: coefficients };
                return true;
            }
        }

        this.mParameters = { average: new Vector(3), coefficients: new Vector(3) };
        return false;
    }

    // Get the parameters for the best fit.
    getParameters(): ApprHeightPlane3Parameters {
        return this.mParameters;
    }

    getMinimumRequired(): number {
        return 3;
    }

    error(point: Vector): number {
        const d = dot(sub(point, this.mParameters.average),
            this.mParameters.coefficients);
        return d * d;
    }

    // The port of the C++ 'dynamic_cast' guarded assignment. C++ value
    // assignment copies, so the port deep-copies the parameters.
    copyParameters(input: ApprQuery<Vector>): void {
        if (input instanceof ApprHeightPlane3) {
            this.mParameters = {
                average: input.mParameters.average.clone(),
                coefficients: input.mParameters.coefficients.clone()
            };
        }
    }
}
