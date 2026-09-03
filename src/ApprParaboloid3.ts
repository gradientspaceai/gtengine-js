// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprParaboloid3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The algorithm for least-squares fitting of a point set by a paraboloid is
// described in
//   https://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf
//
// Port notes:
// * Upstream is templated on T in {float, double, BSRational<*>}; the port
//   is the double instantiation (see PORTING.md). Upstream dispatches the
//   6x6 solve on is_arbitrary_precision<T>: the floating-point branch uses
//   LDLTDecomposition and the rational branch uses LinearSystem (Gaussian
//   elimination). Only the floating-point branch is ported. The
//   LinearSystem-based solve is kept as 'solveGaussian' so an exact path can
//   be added later.
// * The (numPoints, Vector3 const*) overloads collapse into the array
//   overloads; the 'std::array<T,6>& u' output reference and the
//   'T* meanSquareError' optional output become members of the returned
//   result object. The error is always computed (upstream computes it only
//   when the pointer is non-null, an optimization that matters only for the
//   rational instantiation).

import { LDLTDecomposition } from './LDLTDecomposition.js';
import { LinearSystem } from './LinearSystem.js';
import { logAssert } from './Logger.js';
import { Matrix } from './Matrix.js';
import { Vector } from './Vector.js';

// The result of ApprParaboloid3.fit. The paraboloid is
//   z = u[0]*x^2 + u[1]*x*y + u[2]*y^2 + u[3]*x + u[4]*y + u[5]
// The 'success' member is the upstream boolean return value: the linear
// system was solvable. On failure the coefficients are all zero.
export interface ApprParaboloid3FitResult {
    success: boolean;
    u: number[];
    meanSquareError: number;
}

// The result of ApprParaboloid3.fitRobust. The paraboloid is
//   z-c = v0*(x-a)^2 + v1*(x-a)*(y-b) + v2*(y-b)^2 + v3*(x-a) + v4*(y-b) + v5
// where (a,b,c) is the average of the samples.
export interface ApprParaboloid3FitRobustResult {
    success: boolean;
    average: Vector;
    v: number[];
    meanSquareError: number;
}

export class ApprParaboloid3 {
    // Fit with z = u0*x^2 + u1*x*y + u2*y^2 + u3*x + u4*y + u5.
    static fit(points: readonly Vector[]): ApprParaboloid3FitResult {
        const numPoints = points.length;
        logAssert(numPoints >= 6, 'Insufficient points to fit with a paraboloid.');

        const A = new Matrix(6, 6);  // The constructor creates the zero matrix.
        const B = new Vector(6);

        for (let i = 0; i < numPoints; ++i) {
            accumulate(A, B, points[i].values[0], points[i].values[1],
                points[i].values[2]);
        }

        finishSystem(A, B, numPoints);

        const { success, X } = solve(A, B);
        const u = X.values.slice(0, 6);
        let meanSquareError = 0;
        if (success) {
            let totalSqrError = 0;
            for (let i = 0; i < numPoints; ++i) {
                const p = points[i].values;
                const error =
                    u[0] * p[0] * p[0] +
                    u[1] * p[0] * p[1] +
                    u[2] * p[1] * p[1] +
                    u[3] * p[0] +
                    u[4] * p[1] +
                    u[5] - p[2];
                totalSqrError += error * error;
            }
            // Upstream computes sqrt(totalSqrError)/numPoints, which is
            // neither the mean square error nor the root-mean-square error.
            // The quirk is preserved (it is monotone in totalSqrError, so it
            // still orders fits) and reported as an upstream bug suspect.
            meanSquareError = Math.sqrt(totalSqrError) / numPoints;
        }
        return { success, u, meanSquareError };
    }

    // Fit with z-c = v0*(x-a)^2 + v1*(x-a)*(y-b) + v2*(y-b)^2 + v3*(x-a)
    // + v4*(y-b) + v5, where the average of the n samples is (a,b,c) =
    // [sum_{i=0}^{n-1} (x_i,y_i,z_i)]/n. To convert back to the u-polynomial
    // output by fit(...): u0 = v0, u1 = v1, u2 = v2, u3 = v3 - v0*2*a - v1*b,
    // u4 = v4 - v1*a - v2*2*b, and
    // u5 = v0*a^2 + v1*a*b + v2*b^2 - v3*a - v4*b + v5 + c.
    // fitRobust is more expensive to compute than fit, but the effect of
    // rounding errors is mitigated.
    static fitRobust(points: readonly Vector[]): ApprParaboloid3FitRobustResult {
        const numPoints = points.length;
        logAssert(numPoints >= 6, 'Insufficient points to fit with a paraboloid.');

        const A = new Matrix(6, 6);  // The constructor creates the zero matrix.
        const B = new Vector(6);

        // Compute the mean of the points.
        const average = new Vector(3);
        for (let i = 0; i < numPoints; ++i) {
            average.values[0] += points[i].values[0];
            average.values[1] += points[i].values[1];
            average.values[2] += points[i].values[2];
        }
        average.values[0] /= numPoints;
        average.values[1] /= numPoints;
        average.values[2] /= numPoints;

        for (let i = 0; i < numPoints; ++i) {
            accumulate(A, B,
                points[i].values[0] - average.values[0],
                points[i].values[1] - average.values[1],
                points[i].values[2] - average.values[2]);
        }

        finishSystem(A, B, numPoints);

        const { success, X } = solve(A, B);
        const v = X.values.slice(0, 6);
        let meanSquareError = 0;
        if (success) {
            let totalSqrError = 0;
            for (let i = 0; i < numPoints; ++i) {
                const d0 = points[i].values[0] - average.values[0];
                const d1 = points[i].values[1] - average.values[1];
                const d2 = points[i].values[2] - average.values[2];
                const error =
                    v[0] * d0 * d0 +
                    v[1] * d0 * d1 +
                    v[2] * d1 * d1 +
                    v[3] * d0 +
                    v[4] * d1 +
                    v[5] - d2;
                totalSqrError += error * error;
            }
            // See the comment in fit(...) about this formula.
            meanSquareError = Math.sqrt(totalSqrError) / numPoints;
        }
        return { success, average, v, meanSquareError };
    }
}

// Accumulate the moments of one sample (x,y,z) into the linear system.
function accumulate(A: Matrix, B: Vector, x: number, y: number,
    z: number): void {
    const x2 = x * x;
    const xy = x * y;
    const y2 = y * y;
    const zx = z * x;
    const zy = z * y;
    const x3 = x * x2;
    const x2y = x2 * y;
    const xy2 = x * y2;
    const y3 = y * y2;
    const zx2 = z * x2;
    const zxy = z * xy;
    const zy2 = z * y2;
    const x4 = x2 * x2;
    const x3y = x3 * y;
    const x2y2 = x2 * y2;
    const xy3 = x * y3;
    const y4 = y2 * y2;

    A.set(0, 0, A.get(0, 0) + x4);
    A.set(0, 1, A.get(0, 1) + x3y);
    A.set(0, 2, A.get(0, 2) + x2y2);
    A.set(0, 3, A.get(0, 3) + x3);
    A.set(0, 4, A.get(0, 4) + x2y);
    A.set(0, 5, A.get(0, 5) + x2);
    A.set(1, 2, A.get(1, 2) + xy3);
    A.set(1, 4, A.get(1, 4) + xy2);
    A.set(1, 5, A.get(1, 5) + xy);
    A.set(2, 2, A.get(2, 2) + y4);
    A.set(2, 4, A.get(2, 4) + y3);
    A.set(2, 5, A.get(2, 5) + y2);
    A.set(3, 5, A.get(3, 5) + x);
    A.set(4, 5, A.get(4, 5) + y);

    B.values[0] += zx2;
    B.values[1] += zxy;
    B.values[2] += zy2;
    B.values[3] += zx;
    B.values[4] += zy;
    B.values[5] += z;
}

// Fill in the entries of A determined by the symmetry of the moments, then
// scale A and B by dividing by the number of points. The scaling reduces the
// magnitude of the numbers to help with numerical conditioning. The value
// A(5,5) is already scaled to 1 (the flat index into A(5,5) is 35).
function finishSystem(A: Matrix, B: Vector, numPoints: number): void {
    A.set(1, 0, A.get(0, 1));
    A.set(1, 1, A.get(0, 2));
    A.set(1, 3, A.get(0, 4));
    A.set(2, 0, A.get(0, 2));
    A.set(2, 1, A.get(1, 2));
    A.set(2, 3, A.get(1, 4));
    A.set(3, 0, A.get(0, 3));
    A.set(3, 1, A.get(1, 3));
    A.set(3, 2, A.get(2, 3));
    A.set(3, 3, A.get(0, 5));
    A.set(3, 4, A.get(1, 5));
    A.set(4, 0, A.get(0, 4));
    A.set(4, 1, A.get(1, 4));
    A.set(4, 2, A.get(2, 4));
    A.set(4, 3, A.get(3, 4));
    A.set(4, 4, A.get(2, 5));
    A.set(5, 0, A.get(0, 5));
    A.set(5, 1, A.get(1, 5));
    A.set(5, 2, A.get(2, 5));
    A.set(5, 3, A.get(3, 5));
    A.set(5, 4, A.get(4, 5));
    A.set(5, 5, 1);

    for (let i = 0; i < 35; ++i) {
        A.values[i] /= numPoints;
    }

    for (let i = 0; i < 6; ++i) {
        B.values[i] /= numPoints;
    }
}

// The floating-point solve. The LDLTDecomposition class avoids the dynamic
// memory allocation of the GaussianElimination-based LinearSystem solver.
function solve(A: Matrix, B: Vector): { success: boolean, X: Vector } {
    return new LDLTDecomposition(6).solve(A, B);
}

// The solve used by upstream for the arbitrary-precision instantiation,
// where the LDLT approach has too many arithmetic operations to complete in
// a reasonable amount of time. It is exposed for callers who prefer Gaussian
// elimination (for example, when A is positive semidefinite and the LDLT
// factoring encounters a zero pivot).
export function apprParaboloid3SolveGaussian(A: Matrix, B: Vector):
    { success: boolean, X: Vector } {
    const { X, invertible } = LinearSystem.solve(6, A.values, B.values);
    return {
        success: invertible,
        X: invertible ? Vector.fromArray(X) : new Vector(6)
    };
}
