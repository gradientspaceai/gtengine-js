// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprEllipse2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// An ellipse is defined implicitly by (X-C)^T * M * (X-C) = 1, where C is
// the center, M is a positive definite matrix and X is any point on the
// ellipse. The code implements a nonlinear least-squares fitting algorithm
// for the error function
//   F(C,M) = sum_{i=0}^{n-1} ((X[i] - C)^T * M * (X[i] - C) - 1)^2
// for n data points X[0] through X[n-1]. An Ellipse2 object has member
// 'center' that corresponds to C. It also has axes with unit-length
// directions 'axis[]' and corresponding axis half-lengths 'extent[]'. The
// matrix is M = sum_{i=0}^1 axis[i] * axis[i]^T / extent[i]^2, where axis[i]
// is a 2x1 vector and axis[i]^T is a 1x2 vector.
//
// The minimizer uses a 2-step gradient descent algorithm.
//
// Given the current (C,M), locate a minimum of
//   G(t) = F(C - t * dF(C,M)/dC, M)
// for t > 0. The function G(t) >= 0 is a polynomial of degree 4 with
// derivative G'(t) that is a polynomial of degree 3. G'(t) must have a
// positive root because G(0) > 0 and G'(0) < 0 and the G-coefficient of t^4
// is positive. The positive root T that produces the smallest G-value is used
// to update the center C' = C - T * dF/dC(C,M).
//
// Given the current (C,M), locate a minimum of
//   H(t) = F(C, M - t * dF(C,M)/dM)
// for t > 0. The function H(t) >= 0 is a polynomial of degree 2 with
// derivative H'(t) that is a polynomial of degree 1. H'(t) must have a
// positive root because H(0) > 0 and H'(0) < 0 and the H-coefficient of t^2
// is positive. The positive root T that produces the smallest G-value is used
// to update the matrix M' = M - T * dF/dC(C,M) as long as M' is positive
// definite. If M' is not positive definite, the root is halved for a finite
// number of steps until M' is positive definite.
//
// Port notes: 'operator()' becomes 'compute' (the Appr* precedent). The
// 'Ellipse2& ellipse' output reference stays an in/out parameter object
// because it is also the initial guess when
// 'useEllipseForInitialGuess' is true. Upstream ignores the boolean returned
// by GetContainer for the initial oriented box; the port's
// getContainerOrientedBox2 returns null on failure and the port raises an
// error rather than reading an uninitialized box as upstream does.

import { getContainerOrientedBox2 } from './ContOrientedBox2.js';
import type { Ellipse2 } from './Hyperellipsoid.js';
import { logAssert } from './Logger.js';
import { Matrix, addMatrix, divMatrix, mulMatrix, outerProduct } from './Matrix.js';
import { determinant2x2 } from './Matrix2x2.js';
import { RootsPolynomial } from './RootsPolynomial.js';
import { SymmetricEigensolver2x2 } from './SymmetricEigensolver2x2.js';
import { Vector, add, dot, mul, normalize, sub } from './Vector.js';

export class ApprEllipse2 {
    // If you want this function to compute the initial guess for the
    // ellipse, set 'useEllipseForInitialGuess' to false. An oriented
    // bounding box containing the points is used to start the minimizer. Set
    // 'useEllipseForInitialGuess' to true if you want the initial guess to
    // be the input ellipse. This is useful if you want to repeat the query.
    // The returned value is the error function value for the output
    // 'ellipse'.
    compute(points: readonly Vector[], numIterations: number,
        useEllipseForInitialGuess: boolean, ellipse: Ellipse2): number {
        logAssert(points.length > 0, 'ApprEllipse2: no points.');
        for (const point of points) {
            logAssert(point.size === 2, 'ApprEllipse2: points must be 2D.');
        }

        let C: Vector;
        let M = new Matrix(2, 2);  // the zero matrix
        if (useEllipseForInitialGuess) {
            C = ellipse.center.clone();
            for (let i = 0; i < 2; ++i) {
                const product = outerProduct(ellipse.axis[i], ellipse.axis[i]);
                M = addMatrix(M, divMatrix(product,
                    ellipse.extent.values[i] * ellipse.extent.values[i]));
            }
        }
        else {
            const box = getContainerOrientedBox2(points);
            logAssert(box !== null,
                'ApprEllipse2: failed to compute the initial oriented box.');
            C = box.center.clone();
            for (let i = 0; i < 2; ++i) {
                const product = outerProduct(box.axis[i], box.axis[i]);
                M = addMatrix(M, divMatrix(product,
                    box.extent.values[i] * box.extent.values[i]));
            }
        }

        let error = ApprEllipse2.errorFunction(points, C, M);
        for (let i = 0; i < numIterations; ++i) {
            const updateM = ApprEllipse2.updateMatrix(points, C, M);
            error = updateM.error;
            M = updateM.M;
            const updateC = ApprEllipse2.updateCenter(points, M, C);
            error = updateC.error;
            C = updateC.C;
        }

        // Extract the ellipse axes and extents.
        const solver = new SymmetricEigensolver2x2();
        const { evals, evecs } = solver.solve(M.get(0, 0), M.get(0, 1),
            M.get(1, 1), +1);

        ellipse.center = C;
        for (let i = 0; i < 2; ++i) {
            ellipse.axis[i] = Vector.fromArray(evecs[i]);
            ellipse.extent.values[i] = 1 / Math.sqrt(evals[i]);
        }

        return error;
    }

    // The port of UpdateCenter. The C++ in/out reference 'C' becomes an
    // input plus the 'C' field of the returned object.
    private static updateCenter(points: readonly Vector[], M: Matrix,
        C: Vector): { error: number, C: Vector } {
        const epsilon = 1e-06;

        const MDelta = new Array<Vector>(points.length);
        const a = new Array<number>(points.length);
        const invQuantity = 1 / points.length;
        let negDFDC = new Vector(2);
        let aMean = 0, aaMean = 0;
        for (let i = 0; i < points.length; ++i) {
            const Delta = sub(points[i], C);
            MDelta[i] = mulMatrix(M, Delta);
            a[i] = dot(Delta, MDelta[i]) - 1;
            aMean += a[i];
            aaMean += a[i] * a[i];
            negDFDC = add(negDFDC, mul(a[i], MDelta[i]));
        }
        aMean *= invQuantity;
        aaMean *= invQuantity;
        if (normalize(negDFDC) < epsilon) {
            return { error: aaMean, C };
        }

        let bMean = 0, abMean = 0, bbMean = 0;
        const c = dot(negDFDC, mulMatrix(M, negDFDC));
        for (let i = 0; i < points.length; ++i) {
            const b = dot(negDFDC, MDelta[i]);
            bMean += b;
            abMean += a[i] * b;
            bbMean += b * b;
        }
        bMean *= invQuantity;
        abMean *= invQuantity;
        bbMean *= invQuantity;

        // Compute the coefficients of the quartic polynomial q(t) that
        // represents the error function on the given line in the gradient
        // descent minimization.
        const q = [
            aaMean,
            -4 * abMean,
            4 * bbMean + 2 * c * aMean,
            -4 * c * bMean,
            c * c
        ];

        // Compute the coefficients of q'(t).
        const dq = [q[1], 2 * q[2], 3 * q[3], 4 * q[4]];

        // Compute the roots of q'(t).
        const rmMap = RootsPolynomial.solveCubic(dq[0], dq[1], dq[2], dq[3]);

        // Choose the root that leads to the minimum along the gradient
        // descent line and update the center to that point.
        let minError = aaMean;
        let minRoot = 0;
        for (const rm of rmMap) {
            const root = rm.root;
            if (root > 0) {
                const error = q[0] + root * (q[1] + root * (q[2] +
                    root * (q[3] + root * q[4])));
                if (error < minError) {
                    minError = error;
                    minRoot = root;
                }
            }
        }

        if (minRoot > 0) {
            return { error: minError, C: add(C, mul(minRoot, negDFDC)) };
        }
        return { error: aaMean, C };
    }

    // The port of UpdateMatrix. The C++ in/out reference 'M' becomes an
    // input plus the 'M' field of the returned object.
    private static updateMatrix(points: readonly Vector[], C: Vector,
        M: Matrix): { error: number, M: Matrix } {
        const epsilon = 1e-06;

        const Delta = new Array<Vector>(points.length);
        const a = new Array<number>(points.length);
        const invQuantity = 1 / points.length;
        const negDFDM = new Matrix(2, 2);  // zero matrix, symmetric

        let aaMean = 0;
        for (let i = 0; i < points.length; ++i) {
            Delta[i] = sub(points[i], C);
            a[i] = dot(Delta[i], mulMatrix(M, Delta[i])) - 1;
            const twoA = 2 * a[i];
            const d0 = Delta[i].values[0];
            const d1 = Delta[i].values[1];
            negDFDM.set(0, 0, negDFDM.get(0, 0) - a[i] * d0 * d0);
            negDFDM.set(0, 1, negDFDM.get(0, 1) - twoA * d0 * d1);
            negDFDM.set(1, 1, negDFDM.get(1, 1) - a[i] * d1 * d1);
            aaMean += a[i] * a[i];
        }
        aaMean *= invQuantity;

        // Normalize the matrix as if it were a vector of numbers.
        const length = Math.sqrt(
            negDFDM.get(0, 0) * negDFDM.get(0, 0) +
            negDFDM.get(0, 1) * negDFDM.get(0, 1) +
            negDFDM.get(1, 1) * negDFDM.get(1, 1));
        if (length < epsilon) {
            return { error: aaMean, M };
        }
        const invLength = 1 / length;
        negDFDM.set(0, 0, negDFDM.get(0, 0) * invLength);
        negDFDM.set(0, 1, negDFDM.get(0, 1) * invLength);
        negDFDM.set(1, 1, negDFDM.get(1, 1) * invLength);

        // Fill in the lower triangular portion because negDFDM is a
        // symmetric matrix.
        negDFDM.set(1, 0, negDFDM.get(0, 1));

        let abMean = 0, bbMean = 0;
        for (let i = 0; i < points.length; ++i) {
            const b = dot(Delta[i], mulMatrix(negDFDM, Delta[i]));
            abMean += a[i] * b;
            bbMean += b * b;
        }
        abMean *= invQuantity;
        bbMean *= invQuantity;

        // Compute the coefficients of the quadratic polynomial q(t) that
        // represents the error function on the given line in the gradient
        // descent minimization.
        const q = [aaMean, 2 * abMean, bbMean];

        // Compute the coefficients of q'(t).
        const dq = [q[1], 2 * q[2]];

        // Compute the root as long as it is positive and M + root * negDFDM
        // is a positive definite matrix.
        let root = -dq[0] / dq[1];
        if (root > 0) {
            // Use Sylvester's criterion for testing positive definiteness.
            // A for(;;) loop terminates for floating-point arithmetic but
            // not for rational arithmetic. Limit the number of iterations so
            // that the loop terminates for rational arithmetic but 'return'
            // occurs for floating-point arithmetic.
            for (let k = 0; k < 2048; ++k) {
                const nextM = addMatrix(M, mulMatrix(negDFDM, root));
                if (nextM.get(0, 0) > 0) {
                    const det = determinant2x2(nextM);
                    if (det > 0) {
                        const minError = q[0] + root * (q[1] + root * q[2]);
                        return { error: minError, M: nextM };
                    }
                }
                root *= 0.5;
            }
        }
        return { error: aaMean, M };
    }

    private static errorFunction(points: readonly Vector[], C: Vector,
        M: Matrix): number {
        let error = 0;
        for (const P of points) {
            const Delta = sub(P, C);
            const a = dot(Delta, mulMatrix(M, Delta)) - 1;
            error += a * a;
        }
        error /= points.length;
        return error;
    }
}
