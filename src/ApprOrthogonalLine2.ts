// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprOrthogonalLine2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Least-squares fit of a line to (x,y) data by using distance measurements
// orthogonal to the proposed line. The return value is 'true' if and only if
// the fit is unique (the fit always succeeds; the return is 'true' when the
// maximum eigenvalue is unique). The parameters value is a line with (P,D) =
// (origin,direction). The error for S = (x0,y0) is
// (S-P)^T*(I - D*D^T)*(S-P).

import { ApprQuery } from './ApprQuery.js';
import { Line } from './Line.js';
import { SymmetricEigensolver2x2 } from './SymmetricEigensolver2x2.js';
import { Vector, dot } from './Vector.js';

export class ApprOrthogonalLine2 extends ApprQuery<Vector> {
    private mParameters: Line;

    // Initialize the model parameters to zero.
    constructor() {
        super();
        this.mParameters = Line.fromOriginDirection(new Vector(2), new Vector(2));
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
            let covar00 = 0, covar01 = 0, covar11 = 0;
            for (let i = 0; i < numIndices; ++i) {
                const p = points[indices[i]];
                const d0 = p.values[0] - mean.values[0];
                const d1 = p.values[1] - mean.values[1];
                covar00 += d0 * d0;
                covar01 += d0 * d1;
                covar11 += d1 * d1;
            }

            // Solve the eigensystem.
            const es = new SymmetricEigensolver2x2();
            const result = es.solve(covar00, covar01, covar11, +1);

            // The line direction is the eigenvector in the direction of
            // largest variance of the points.
            this.mParameters.origin = mean;
            this.mParameters.direction = Vector.fromArray(result.evecs[1]);

            // The fitted line is unique when the maximum eigenvalue has
            // multiplicity 1.
            return result.evals[0] < result.evals[1];
        }

        this.mParameters = Line.fromOriginDirection(new Vector(2), new Vector(2));
        return false;
    }

    // Get the parameters for the best fit.
    getParameters(): Line {
        return this.mParameters;
    }

    getMinimumRequired(): number {
        return 2;
    }

    error(point: Vector): number {
        const diff = new Vector(2);
        diff.values[0] = point.values[0] - this.mParameters.origin.values[0];
        diff.values[1] = point.values[1] - this.mParameters.origin.values[1];
        const sqrlen = dot(diff, diff);
        const d = dot(diff, this.mParameters.direction);
        return Math.abs(sqrlen - d * d);
    }

    // The port of the C++ 'dynamic_cast' guarded assignment. C++ value
    // assignment copies, so the port deep-copies the line.
    copyParameters(input: ApprQuery<Vector>): void {
        if (input instanceof ApprOrthogonalLine2) {
            this.mParameters = input.mParameters.clone();
        }
    }
}
