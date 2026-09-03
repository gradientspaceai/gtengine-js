// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprGreatCircle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Least-squares fit of a great circle to unit-length vectors (x,y,z) by
// using distance measurements orthogonal (and measured along great circles)
// to the proposed great circle. The input points are unit length. The
// returned value is unit length, call it N. The fitted great circle is
// defined by Dot(N,X) = 0, where X is a unit-length vector on the great
// circle.
//
// Port notes: the (numPoints, Vector3 const*) pointer pair collapses to a
// single array of 3D Vectors, and the output reference parameters become the
// return values (a Vector for ApprGreatCircle3, an object literal for
// ApprGreatArc3). Following the Appr* precedent, 'operator()' becomes
// 'compute'.

import { GTE_C_TWO_PI } from './Constants.js';
import { logAssert } from './Logger.js';
import { SymmetricEigensolver3x3 } from './SymmetricEigensolver3x3.js';
import { Vector, add, dot, mul, normalize } from './Vector.js';
import { computeOrthogonalComplement3 } from './Vector3.js';

export class ApprGreatCircle3 {
    // The returned unit-length vector N defines the fitted great circle by
    // Dot(N,X) = 0.
    compute(points: readonly Vector[]): Vector {
        const numPoints = points.length;
        logAssert(numPoints > 0, 'ApprGreatCircle3: no points.');

        // Compute the covariance matrix of the vectors.
        let covar00 = 0, covar01 = 0, covar02 = 0;
        let covar11 = 0, covar12 = 0, covar22 = 0;
        for (let i = 0; i < numPoints; ++i) {
            const diff = points[i];
            logAssert(diff.size === 3, 'ApprGreatCircle3: points must be 3D.');
            covar00 += diff.values[0] * diff.values[0];
            covar01 += diff.values[0] * diff.values[1];
            covar02 += diff.values[0] * diff.values[2];
            covar11 += diff.values[1] * diff.values[1];
            covar12 += diff.values[1] * diff.values[2];
            covar22 += diff.values[2] * diff.values[2];
        }

        const invNumPoints = 1 / numPoints;
        covar00 *= invNumPoints;
        covar01 *= invNumPoints;
        covar02 *= invNumPoints;
        covar11 *= invNumPoints;
        covar12 *= invNumPoints;
        covar22 *= invNumPoints;

        // Solve the eigensystem. The eigenvector associated with the
        // smallest eigenvalue is the great-circle normal.
        const es = new SymmetricEigensolver3x3();
        const result = es.solve(covar00, covar01, covar02, covar11, covar12,
            covar22, false, +1);
        return Vector.fromArray(result.evecs[0]);
    }
}

// The result of ApprGreatArc3.compute. The endpoints arcEnd0 and arcEnd1 are
// perpendicular to normal. When the arc is viewed by looking at the plane of
// the great circle with a view direction of -normal, the arc is traversed
// counterclockwise starting at arcEnd0 and ending at arcEnd1.
export interface ApprGreatArc3Result {
    normal: Vector;
    arcEnd0: Vector;
    arcEnd1: Vector;
}

// In addition to the least-squares fit of a great circle, the input vectors
// are projected onto that circle. The sector of smallest angle (possibly
// obtuse) that contains the points is computed. The endpoints of the arc of
// the sector are returned.
export class ApprGreatArc3 {
    compute(points: readonly Vector[]): ApprGreatArc3Result {
        const numPoints = points.length;
        logAssert(numPoints > 0, 'ApprGreatArc3: no points.');

        // Get the least-squares great circle for the vectors. The circle is
        // on the plane Dot(N,X) = 0. Generate a basis from N.
        const basis: Vector[] = [new Vector(3), new Vector(3), new Vector(3)];
        basis[0] = new ApprGreatCircle3().compute(points);
        computeOrthogonalComplement3(1, basis);

        // The vectors are X[i] = u[i]*U + v[i]*V + w[i]*N. The projections
        // are
        //   P[i] = (u[i]*U + v[i]*V)/sqrt(u[i]*u[i] + v[i]*v[i])
        // The great circle is parameterized by
        //   C(t) = cos(t)*U + sin(t)*V
        // Compute the angles t in [-pi,pi] for the projections onto the
        // great circle. It is not necessary to normalize (u[i],v[i]),
        // instead computing t = atan2(v[i],u[i]). Each item is (u, v, angle).
        const items: Array<[number, number, number]> = new Array(numPoints);
        for (let i = 0; i < numPoints; ++i) {
            const u = dot(basis[1], points[i]);
            const v = dot(basis[2], points[i]);
            items[i] = [u, v, Math.atan2(v, u)];
        }
        items.sort((item0, item1) => item0[2] - item1[2]);

        // Locate the pair of consecutive angles whose difference is a
        // maximum. Effectively, we are constructing a cone of minimum angle
        // that contains the unit-length vectors.
        const numPointsM1 = numPoints - 1;
        let maxDiff = GTE_C_TWO_PI + items[0][2] - items[numPointsM1][2];
        let end0 = 0, end1 = numPointsM1;
        for (let i0 = 0, i1 = 1; i0 < numPointsM1; i0 = i1++) {
            const diff = items[i1][2] - items[i0][2];
            if (diff > maxDiff) {
                maxDiff = diff;
                end0 = i1;
                end1 = i0;
            }
        }

        const normal = basis[0];
        const arcEnd0 = add(mul(items[end0][0], basis[1]),
            mul(items[end0][1], basis[2]));
        const arcEnd1 = add(mul(items[end1][0], basis[1]),
            mul(items[end1][1], basis[2]));
        normalize(arcEnd0);
        normalize(arcEnd1);
        return { normal, arcEnd0, arcEnd1 };
    }
}
