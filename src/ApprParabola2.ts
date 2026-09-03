// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprParabola2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The algorithm for least-squares fitting of a point set by a parabola is
// described in
//   https://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf
//
// Port notes:
// * Upstream is templated on T in {float, double, BSRational<*>}; the port
//   is the double instantiation (see PORTING.md).
// * The (numPoints, Vector2 const*) overloads collapse into the array
//   overloads; the 'std::array<T,3>& u' output reference and the
//   'T* meanSquareError' optional output become members of the returned
//   result object. The error is always computed (upstream computes it only
//   when the pointer is non-null, an optimization that matters only for the
//   rational instantiation).

import { logAssert } from './Logger.js';
import { LinearSystem } from './LinearSystem.js';
import { Matrix } from './Matrix.js';
import { Vector } from './Vector.js';

// The result of ApprParabola2.fit. The parabola is
//   y = u[0] * x^2 + u[1] * x + u[2]
// The 'success' member is the upstream boolean return value: the linear
// system was solvable. On failure the coefficients are all zero.
export interface ApprParabola2FitResult {
    success: boolean;
    u: number[];
    meanSquareError: number;
}

// The result of ApprParabola2.fitRobust. The parabola is
//   y - b = v[0] * (x - a)^2 + v[1] * (x - a) + v[2]
// where (a,b) is the average of the samples.
export interface ApprParabola2FitRobustResult {
    success: boolean;
    average: Vector;
    v: number[];
    meanSquareError: number;
}

export class ApprParabola2 {
    // Fit with y = u0*x^2 + u1*x + u2. The code uses a specialized 3x3
    // linear system solver. This is faster than Gaussian elimination, an
    // LDL^T decomposition, or a 3x3 eigensystem solver.
    static fit(points: readonly Vector[]): ApprParabola2FitResult {
        const numPoints = points.length;
        logAssert(numPoints >= 3, 'Insufficient points to fit with a parabola.');

        const A = new Matrix(3, 3);  // The constructor creates the zero matrix.
        const B = new Vector(3);

        for (let i = 0; i < numPoints; ++i) {
            const point = points[i].values;
            const x2 = point[0] * point[0];
            const x3 = point[0] * x2;
            const x4 = x2 * x2;
            const x2y = x2 * point[1];
            const xy = point[0] * point[1];

            A.set(0, 0, A.get(0, 0) + x4);
            A.set(0, 1, A.get(0, 1) + x3);
            A.set(0, 2, A.get(0, 2) + x2);
            A.set(1, 2, A.get(1, 2) + point[0]);

            B.values[0] += x2y;
            B.values[1] += xy;
            B.values[2] += point[1];
        }

        finishSystem(A, B, numPoints);

        const { X, invertible } = LinearSystem.solve3x3(A, B);
        const u = [X.values[0], X.values[1], X.values[2]];
        let meanSquareError = 0;
        if (invertible) {
            let totalSqrError = 0;
            for (let i = 0; i < numPoints; ++i) {
                const point = points[i].values;
                const error =
                    u[0] * point[0] * point[0] +
                    u[1] * point[0] +
                    u[2] - point[1];
                totalSqrError += error * error;
            }
            // Upstream computes sqrt(totalSqrError)/numPoints, which is
            // neither the mean square error nor the root-mean-square error.
            // The quirk is preserved (it is monotone in totalSqrError, so it
            // still orders fits) and reported as an upstream bug suspect.
            meanSquareError = Math.sqrt(totalSqrError) / numPoints;
        }
        return { success: invertible, u, meanSquareError };
    }

    // Fit with y-b = v0*(x-a)^2 + v1*(x-a) + v2, where the average of the n
    // samples is (a,b) = [sum_{i=0}{n-1} (x_i,y_i)]/n. To convert back to the
    // u-polynomial output by fit(...): u0 = v0, u1 = v1 - 2*v0*a, and
    // u2 = v0*a^2 - v1*a + v2 + b. fitRobust is more expensive to compute
    // than fit, but the effect of rounding errors is mitigated.
    static fitRobust(points: readonly Vector[]): ApprParabola2FitRobustResult {
        const numPoints = points.length;
        logAssert(numPoints >= 3, 'Insufficient points to fit with a parabola.');

        const A = new Matrix(3, 3);  // The constructor creates the zero matrix.
        const B = new Vector(3);

        // Compute the mean of the points.
        const average = new Vector(2);
        for (let i = 0; i < numPoints; ++i) {
            average.values[0] += points[i].values[0];
            average.values[1] += points[i].values[1];
        }
        average.values[0] /= numPoints;
        average.values[1] /= numPoints;

        for (let i = 0; i < numPoints; ++i) {
            const d0 = points[i].values[0] - average.values[0];
            const d1 = points[i].values[1] - average.values[1];
            const x2 = d0 * d0;
            const x3 = d0 * x2;
            const x4 = x2 * x2;
            const x2y = x2 * d1;
            const xy = d0 * d1;

            A.set(0, 0, A.get(0, 0) + x4);
            A.set(0, 1, A.get(0, 1) + x3);
            A.set(0, 2, A.get(0, 2) + x2);
            A.set(1, 2, A.get(1, 2) + d0);

            B.values[0] += x2y;
            B.values[1] += xy;
            B.values[2] += d1;
        }

        finishSystem(A, B, numPoints);

        const { X, invertible } = LinearSystem.solve3x3(A, B);
        const v = [X.values[0], X.values[1], X.values[2]];
        let meanSquareError = 0;
        if (invertible) {
            let totalSqrError = 0;
            for (let i = 0; i < numPoints; ++i) {
                const d0 = points[i].values[0] - average.values[0];
                const d1 = points[i].values[1] - average.values[1];
                const error = v[0] * d0 * d0 + v[1] * d0 + v[2] - d1;
                totalSqrError += error * error;
            }
            // See the comment in fit(...) about this formula.
            meanSquareError = Math.sqrt(totalSqrError) / numPoints;
        }
        return { success: invertible, average, v, meanSquareError };
    }
}

// Fill in the symmetric entries of A and scale A and B by dividing by the
// number of points. The scaling reduces the magnitude of the numbers to help
// with numerical conditioning. The value A(2,2) is already scaled to 1 (the
// flat index into A(2,2) is 8).
function finishSystem(A: Matrix, B: Vector, numPoints: number): void {
    A.set(1, 0, A.get(0, 1));
    A.set(1, 1, A.get(0, 2));
    A.set(2, 0, A.get(0, 2));
    A.set(2, 1, A.get(1, 2));
    A.set(2, 2, 1);

    for (let i = 0; i < 8; ++i) {
        A.values[i] /= numPoints;
    }

    for (let i = 0; i < 3; ++i) {
        B.values[i] /= numPoints;
    }
}
