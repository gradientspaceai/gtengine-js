// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprSphere3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Least-squares fit of a sphere to a set of points. The algorithms are
// described in Section 5 of
//   https://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf
// fitUsingLengths uses the algorithm of Section 5.1.
// fitUsingSquaredLengths uses the algorithm of Section 5.2.
//
// Port notes: see ApprCircle2.ts; the (numPoints, Vector3 const*) pointer
// pairs collapse to arrays of 3D Vectors and the 'Sphere3& sphere' output
// reference stays an in/out parameter object because fitUsingLengths reads
// the incoming center as the initial guess for its iteration.

import { Hypersphere } from './Hypersphere';
import { Vector, dot, length } from './Vector';

export class ApprSphere3 {
    // The return value is 'true' when the linear system of the algorithm is
    // solvable, 'false' otherwise. If 'false' is returned, the sphere center
    // and radius are set to zero values.
    fitUsingSquaredLengths(points: readonly Vector[], sphere: Hypersphere): boolean {
        const numPoints = points.length;

        // Compute the average of the data points.
        let a0 = 0, a1 = 0, a2 = 0;
        for (let i = 0; i < numPoints; ++i) {
            a0 += points[i].values[0];
            a1 += points[i].values[1];
            a2 += points[i].values[2];
        }
        const invNumPoints = 1 / numPoints;
        a0 *= invNumPoints;
        a1 *= invNumPoints;
        a2 *= invNumPoints;

        // Compute the covariance matrix M of the Y[i] = X[i]-A and the
        // right-hand side R of the linear system M*(C-A) = R.
        let m00 = 0, m01 = 0, m02 = 0, m11 = 0, m12 = 0, m22 = 0;
        let r0 = 0, r1 = 0, r2 = 0;
        for (let i = 0; i < numPoints; ++i) {
            const y0 = points[i].values[0] - a0;
            const y1 = points[i].values[1] - a1;
            const y2 = points[i].values[2] - a2;
            const y0y0 = y0 * y0;
            const y0y1 = y0 * y1;
            const y0y2 = y0 * y2;
            const y1y1 = y1 * y1;
            const y1y2 = y1 * y2;
            const y2y2 = y2 * y2;
            m00 += y0y0;
            m01 += y0y1;
            m02 += y0y2;
            m11 += y1y1;
            m12 += y1y2;
            m22 += y2y2;
            const sumSqr = y0y0 + y1y1 + y2y2;
            r0 += sumSqr * y0;
            r1 += sumSqr * y1;
            r2 += sumSqr * y2;
        }
        r0 *= 0.5;
        r1 *= 0.5;
        r2 *= 0.5;

        // Solve the linear system M*(C-A) = R for the center C.
        const cof00 = m11 * m22 - m12 * m12;
        const cof01 = m02 * m12 - m01 * m22;
        const cof02 = m01 * m12 - m02 * m11;
        const det = m00 * cof00 + m01 * cof01 + m02 * cof02;
        if (det !== 0) {
            const cof11 = m00 * m22 - m02 * m02;
            const cof12 = m01 * m02 - m00 * m12;
            const cof22 = m00 * m11 - m01 * m01;
            sphere.center.values[0] = a0 + (cof00 * r0 + cof01 * r1 + cof02 * r2) / det;
            sphere.center.values[1] = a1 + (cof01 * r0 + cof11 * r1 + cof12 * r2) / det;
            sphere.center.values[2] = a2 + (cof02 * r0 + cof12 * r1 + cof22 * r2) / det;
            let rsqr = 0;
            for (let i = 0; i < numPoints; ++i) {
                const d0 = points[i].values[0] - sphere.center.values[0];
                const d1 = points[i].values[1] - sphere.center.values[1];
                const d2 = points[i].values[2] - sphere.center.values[2];
                rsqr += d0 * d0 + d1 * d1 + d2 * d2;
            }
            rsqr *= invNumPoints;
            sphere.radius = Math.sqrt(rsqr);
            return true;
        }
        else {
            sphere.center.values[0] = 0;
            sphere.center.values[1] = 0;
            sphere.center.values[2] = 0;
            sphere.radius = 0;
            return false;
        }
    }

    // Fit the points using lengths to drive the least-squares algorithm. If
    // initialCenterIsAverage is 'true', the initial guess for the sphere
    // center is the average of the data points. If it is 'false', the
    // incoming sphere center is used as-is to start the iterative algorithm.
    // (The upstream documentation states the opposite for the two cases; the
    // code, which the port preserves, keys the average off 'true'.) Starting
    // from a supplied center tends to converge more rapidly than starting
    // from the average of the points, which is slow to converge when the
    // points are clustered along a small solid angle. Either way, this
    // algorithm can be much slower than fitUsingSquaredLengths.
    //
    // The value epsilon may be chosen as a positive number for the
    // comparison of consecutive estimated sphere centers, terminating the
    // iterations when the center difference has length less than or equal to
    // epsilon.
    //
    // The return value is the number of iterations used. If it is the input
    // maxIterations, you can either accept the result or polish the result
    // by calling the function again with initialCenterIsAverage set to
    // 'false', so that the current center is the initial guess.
    fitUsingLengths(points: readonly Vector[], maxIterations: number,
        initialCenterIsAverage: boolean, sphere: Hypersphere,
        epsilon: number = 0): number {
        const numPoints = points.length;

        // Compute the average of the data points.
        const average = points[0].clone();
        for (let i = 1; i < numPoints; ++i) {
            average.values[0] += points[i].values[0];
            average.values[1] += points[i].values[1];
            average.values[2] += points[i].values[2];
        }
        const invNumPoints = 1 / numPoints;
        average.values[0] *= invNumPoints;
        average.values[1] *= invNumPoints;
        average.values[2] *= invNumPoints;

        // The initial guess for the center.
        if (initialCenterIsAverage) {
            sphere.center = average.clone();
        }

        const epsilonSqr = epsilon * epsilon;
        let iteration = 0;
        for (iteration = 0; iteration < maxIterations; ++iteration) {
            // Update the iterates.
            const current = sphere.center.clone();

            // Compute average L, dL/da, dL/db, dL/dc.
            let lenAverage = 0;
            const derLenAverage = new Vector(3);
            for (let i = 0; i < numPoints; ++i) {
                const diff = new Vector(3);
                diff.values[0] = points[i].values[0] - sphere.center.values[0];
                diff.values[1] = points[i].values[1] - sphere.center.values[1];
                diff.values[2] = points[i].values[2] - sphere.center.values[2];
                const len = length(diff);
                if (len > 0) {
                    lenAverage += len;
                    const invLength = 1 / len;
                    derLenAverage.values[0] -= invLength * diff.values[0];
                    derLenAverage.values[1] -= invLength * diff.values[1];
                    derLenAverage.values[2] -= invLength * diff.values[2];
                }
            }
            lenAverage *= invNumPoints;
            derLenAverage.values[0] *= invNumPoints;
            derLenAverage.values[1] *= invNumPoints;
            derLenAverage.values[2] *= invNumPoints;

            for (let d = 0; d < 3; ++d) {
                sphere.center.values[d] = average.values[d]
                    + lenAverage * derLenAverage.values[d];
            }
            sphere.radius = lenAverage;

            const diff = new Vector(3);
            for (let d = 0; d < 3; ++d) {
                diff.values[d] = sphere.center.values[d] - current.values[d];
            }
            const diffSqrLen = dot(diff, diff);
            if (diffSqrLen <= epsilonSqr) {
                break;
            }
        }

        return iteration;
    }
}
