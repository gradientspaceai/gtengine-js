// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprOrthogonalPlane3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Least-squares fit of a plane to (x,y,z) data by using distance
// measurements orthogonal to the proposed plane. The return value is 'true'
// if and only if the fit is unique (the fit always succeeds; the return is
// 'true' when the minimum eigenvalue is unique). The parameters value is
// (P,N) = (origin,normal). The error for S = (x0,y0,z0) is |Dot(N,S-P)|.
//
// Port notes: upstream stores the parameters in a
// std::pair<Vector3,Vector3>; per PORTING.md the pair becomes an object with
// named fields (origin, normal).

import { ApprQuery } from './ApprQuery';
import { SymmetricEigensolver3x3 } from './SymmetricEigensolver3x3';
import { Vector, dot, sub } from './Vector';

export interface ApprOrthogonalPlane3Parameters {
    // P, the average of the sample points.
    origin: Vector;

    // N, the unit-length plane normal.
    normal: Vector;
}

export class ApprOrthogonalPlane3 extends ApprQuery<Vector> {
    private mParameters: ApprOrthogonalPlane3Parameters;

    // Initialize the model parameters to zero.
    constructor() {
        super();
        this.mParameters = { origin: new Vector(3), normal: new Vector(3) };
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

            // The plane normal is the eigenvector in the direction of
            // smallest variance of the points.
            this.mParameters = {
                origin: mean,
                normal: Vector.fromArray(result.evecs[0])
            };

            // The fitted plane is unique when the minimum eigenvalue has
            // multiplicity 1.
            return result.evals[0] < result.evals[1];
        }

        this.mParameters = { origin: new Vector(3), normal: new Vector(3) };
        return false;
    }

    // Get the parameters for the best fit.
    getParameters(): ApprOrthogonalPlane3Parameters {
        return this.mParameters;
    }

    getMinimumRequired(): number {
        return 3;
    }

    error(point: Vector): number {
        return Math.abs(dot(this.mParameters.normal,
            sub(point, this.mParameters.origin)));
    }

    // The port of the C++ 'dynamic_cast' guarded assignment. C++ value
    // assignment copies, so the port deep-copies the parameters.
    copyParameters(input: ApprQuery<Vector>): void {
        if (input instanceof ApprOrthogonalPlane3) {
            this.mParameters = {
                origin: input.mParameters.origin.clone(),
                normal: input.mParameters.normal.clone()
            };
        }
    }
}
