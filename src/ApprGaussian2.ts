// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprGaussian2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Fit points with a Gaussian distribution. The center is the mean of the
// points, the axes are the eigenvectors of the covariance matrix and the
// extents are the eigenvalues of the covariance matrix and are returned in
// increasing order. An oriented box is used to store the mean, axes and
// extents.
//
// Port notes: the (numPoints, Vector2 const*, numIndices, int32_t const*)
// pointer pairs collapse to arrays (see ApprQuery.ts). The observation type
// is a 2D Vector.

import { ApprQuery } from './ApprQuery';
import { OrientedBox } from './OrientedBox';
import { SymmetricEigensolver2x2 } from './SymmetricEigensolver2x2';
import { Vector, dot } from './Vector';

export class ApprGaussian2 extends ApprQuery<Vector> {
    private mParameters: OrientedBox;

    // Initialize the model parameters to zero.
    constructor() {
        super();
        this.mParameters = new OrientedBox(2);
        this.setZeroParameters();
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
            const invSize = 1 / numIndices;
            mean.values[0] *= invSize;
            mean.values[1] *= invSize;

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
            covar00 *= invSize;
            covar01 *= invSize;
            covar11 *= invSize;

            // Solve the eigensystem.
            const es = new SymmetricEigensolver2x2();
            const result = es.solve(covar00, covar01, covar11, +1);
            this.mParameters.center = mean;
            this.mParameters.axis[0] = Vector.fromArray(result.evecs[0]);
            this.mParameters.axis[1] = Vector.fromArray(result.evecs[1]);
            this.mParameters.extent = Vector.fromArray(result.evals);
            return true;
        }

        this.setZeroParameters();
        return false;
    }

    // Get the parameters for the best fit.
    getParameters(): OrientedBox {
        return this.mParameters;
    }

    getMinimumRequired(): number {
        return 2;
    }

    error(point: Vector): number {
        const diff = new Vector(2);
        diff.values[0] = point.values[0] - this.mParameters.center.values[0];
        diff.values[1] = point.values[1] - this.mParameters.center.values[1];
        let error = 0;
        for (let i = 0; i < 2; ++i) {
            if (this.mParameters.extent.values[i] > 0) {
                const ratio = dot(diff, this.mParameters.axis[i])
                    / this.mParameters.extent.values[i];
                error += ratio * ratio;
            }
        }
        return error;
    }

    // The port of the C++ 'dynamic_cast' guarded assignment. C++ value
    // assignment copies, so the port deep-copies the oriented box.
    copyParameters(input: ApprQuery<Vector>): void {
        if (input instanceof ApprGaussian2) {
            this.mParameters = input.mParameters.clone();
        }
    }

    private setZeroParameters(): void {
        this.mParameters.center = new Vector(2);
        this.mParameters.axis[0] = new Vector(2);
        this.mParameters.axis[1] = new Vector(2);
        this.mParameters.extent = new Vector(2);
    }
}
