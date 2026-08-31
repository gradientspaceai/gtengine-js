// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprGaussian3.h
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
// Port notes: the pointer/count pairs collapse to arrays (see ApprQuery.ts).
// The observation type is a 3D Vector.

import { ApprQuery } from './ApprQuery';
import { OrientedBox } from './OrientedBox';
import { SymmetricEigensolver3x3 } from './SymmetricEigensolver3x3';
import { Vector, dot } from './Vector';

export class ApprGaussian3 extends ApprQuery<Vector> {
    private mParameters: OrientedBox;

    // Initialize the model parameters to zero.
    constructor() {
        super();
        this.mParameters = new OrientedBox(3);
        this.setZeroParameters();
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
            const invSize = 1 / numIndices;
            for (let d = 0; d < 3; ++d) {
                mean.values[d] *= invSize;
            }

            // Compute the covariance matrix of the points.
            let covar00 = 0, covar01 = 0, covar02 = 0;
            let covar11 = 0, covar12 = 0, covar22 = 0;
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
                covar22 += d2 * d2;
            }
            covar00 *= invSize;
            covar01 *= invSize;
            covar02 *= invSize;
            covar11 *= invSize;
            covar12 *= invSize;
            covar22 *= invSize;

            // Solve the eigensystem.
            const es = new SymmetricEigensolver3x3();
            const result = es.solve(covar00, covar01, covar02, covar11,
                covar12, covar22, false, +1);
            this.mParameters.center = mean;
            this.mParameters.axis[0] = Vector.fromArray(result.evecs[0]);
            this.mParameters.axis[1] = Vector.fromArray(result.evecs[1]);
            this.mParameters.axis[2] = Vector.fromArray(result.evecs[2]);
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
        const diff = new Vector(3);
        for (let d = 0; d < 3; ++d) {
            diff.values[d] = point.values[d] - this.mParameters.center.values[d];
        }
        let error = 0;
        for (let i = 0; i < 3; ++i) {
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
        if (input instanceof ApprGaussian3) {
            this.mParameters = input.mParameters.clone();
        }
    }

    private setZeroParameters(): void {
        this.mParameters.center = new Vector(3);
        this.mParameters.axis[0] = new Vector(3);
        this.mParameters.axis[1] = new Vector(3);
        this.mParameters.axis[2] = new Vector(3);
        this.mParameters.extent = new Vector(3);
    }
}
