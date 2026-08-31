// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprCircle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Least-squares fit of a circle to a set of points. The algorithms are
// described in Section 5 of
//   https://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf
// fitUsingLengths uses the algorithm of Section 5.1.
// fitUsingSquaredLengths uses the algorithm of Section 5.2.
//
// Port notes: the C++ functions take (numPoints, Vector2 const*) pointer
// pairs; TypeScript arrays carry their length, so each pair collapses to a
// single array of 2D Vectors. The 'Circle2& circle' output reference stays an
// in/out parameter object rather than becoming a returned literal, because
// fitUsingLengths reads the incoming center as the initial guess for its
// iteration.

import { Hypersphere } from './Hypersphere';
import { Vector, dot, length } from './Vector';

export class ApprCircle2 {
    // The return value is 'true' when the linear system of the algorithm is
    // solvable, 'false' otherwise. If 'false' is returned, the circle center
    // and radius are set to zero values.
    fitUsingSquaredLengths(points: readonly Vector[], circle: Hypersphere): boolean {
        const numPoints = points.length;

        // Compute the average of the data points.
        let a0 = 0, a1 = 0;
        for (let i = 0; i < numPoints; ++i) {
            a0 += points[i].values[0];
            a1 += points[i].values[1];
        }
        const invNumPoints = 1 / numPoints;
        a0 *= invNumPoints;
        a1 *= invNumPoints;

        // Compute the covariance matrix M of the Y[i] = X[i]-A and the
        // right-hand side R of the linear system M*(C-A) = R.
        let m00 = 0, m01 = 0, m11 = 0;
        let r0 = 0, r1 = 0;
        for (let i = 0; i < numPoints; ++i) {
            const y0 = points[i].values[0] - a0;
            const y1 = points[i].values[1] - a1;
            const y0y0 = y0 * y0;
            const y0y1 = y0 * y1;
            const y1y1 = y1 * y1;
            m00 += y0y0;
            m01 += y0y1;
            m11 += y1y1;
            r0 += (y0y0 + y1y1) * y0;
            r1 += (y0y0 + y1y1) * y1;
        }
        r0 *= 0.5;
        r1 *= 0.5;

        // Solve the linear system M*(C-A) = R for the center C.
        const det = m00 * m11 - m01 * m01;
        if (det !== 0) {
            circle.center.values[0] = a0 + (m11 * r0 - m01 * r1) / det;
            circle.center.values[1] = a1 + (m00 * r1 - m01 * r0) / det;
            let rsqr = 0;
            for (let i = 0; i < numPoints; ++i) {
                const d0 = points[i].values[0] - circle.center.values[0];
                const d1 = points[i].values[1] - circle.center.values[1];
                rsqr += d0 * d0 + d1 * d1;
            }
            rsqr *= invNumPoints;
            circle.radius = Math.sqrt(rsqr);
            return true;
        }
        else {
            circle.center.values[0] = 0;
            circle.center.values[1] = 0;
            circle.radius = 0;
            return false;
        }
    }

    // Fit the points using lengths to drive the least-squares algorithm. If
    // initialCenterIsAverage is 'true', the initial guess for the circle
    // center is the average of the data points. If it is 'false', the
    // incoming circle center is used as-is to start the iterative algorithm.
    // (The upstream documentation states the opposite for the two cases; the
    // code, which the port preserves, keys the average off 'true'.) Starting
    // from a supplied center tends to converge more rapidly than starting
    // from the average of the points, which is slow to converge when the
    // points are clustered along a small arc. Either way, this algorithm can
    // be much slower than fitUsingSquaredLengths.
    //
    // The value epsilon may be chosen as a positive number for the
    // comparison of consecutive estimated circle centers, terminating the
    // iterations when the center difference has length less than or equal to
    // epsilon.
    //
    // The return value is the number of iterations used. If it is the input
    // maxIterations, you can either accept the result or polish the result
    // by calling the function again with initialCenterIsAverage set to
    // 'false', so that the current center is the initial guess.
    fitUsingLengths(points: readonly Vector[], maxIterations: number,
        initialCenterIsAverage: boolean, circle: Hypersphere,
        epsilon: number = 0): number {
        const numPoints = points.length;

        // Compute the average of the data points.
        const average = points[0].clone();
        for (let i = 1; i < numPoints; ++i) {
            average.values[0] += points[i].values[0];
            average.values[1] += points[i].values[1];
        }
        const invNumPoints = 1 / numPoints;
        average.values[0] *= invNumPoints;
        average.values[1] *= invNumPoints;

        // The initial guess for the center.
        if (initialCenterIsAverage) {
            circle.center = average.clone();
        }

        const epsilonSqr = epsilon * epsilon;
        let iteration = 0;
        for (iteration = 0; iteration < maxIterations; ++iteration) {
            // Update the iterates.
            const current = circle.center.clone();

            // Compute average L, dL/da, dL/db.
            let lenAverage = 0;
            const derLenAverage = new Vector(2);
            for (let i = 0; i < numPoints; ++i) {
                const diff = new Vector(2);
                diff.values[0] = points[i].values[0] - circle.center.values[0];
                diff.values[1] = points[i].values[1] - circle.center.values[1];
                const len = length(diff);
                if (len > 0) {
                    lenAverage += len;
                    const invLength = 1 / len;
                    derLenAverage.values[0] -= invLength * diff.values[0];
                    derLenAverage.values[1] -= invLength * diff.values[1];
                }
            }
            lenAverage *= invNumPoints;
            derLenAverage.values[0] *= invNumPoints;
            derLenAverage.values[1] *= invNumPoints;

            circle.center.values[0] = average.values[0] + lenAverage * derLenAverage.values[0];
            circle.center.values[1] = average.values[1] + lenAverage * derLenAverage.values[1];
            circle.radius = lenAverage;

            const diff = new Vector(2);
            diff.values[0] = circle.center.values[0] - current.values[0];
            diff.values[1] = circle.center.values[1] - current.values[1];
            const diffSqrLen = dot(diff, diff);
            if (diffSqrLen <= epsilonSqr) {
                break;
            }
        }

        return iteration;
    }
}
