// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprQuadratic2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: the C++ classes have a single 'operator()', so per PORTING.md
// they become 'compute(...)'. The '(numPoints, Vector2 const*)' pointer pair
// collapses to a single array of 2D Vectors. The 'std::array<Real,6>&'
// output reference of ApprQuadratic2 becomes a field of the returned object
// (the return value being the fit measure); ApprQuadraticCircle2 keeps the
// 'Circle2&' as an output parameter object it writes into, matching
// ApprCircle2.ts.
//
// Note that the header comment claims M = (sum_i V[i])(sum_i V[i])^T, but
// the code accumulates M = sum_i V[i]*V[i]^T (the correct normal-equation
// matrix). The port implements the code, not the comment.

import { Hypersphere } from './Hypersphere.js';
import { Matrix } from './Matrix.js';
import { SymmetricEigensolver } from './SymmetricEigensolver.js';
import { Vector, dot } from './Vector.js';

export interface ApprQuadratic2Result {
    // The coefficients C[0..5] of the quadratic fit (a unit-length
    // eigenvector of M).
    coefficients: number[];

    // A nonnegative measure of the fit (the minimum eigenvalue of M);
    // 0 means an exact fit, positive otherwise.
    minEigenvalue: number;
}

// The quadratic fit is
//   0 = C[0] + C[1]*x + C[2]*y + C[3]*x^2 + C[4]*x*y + C[5]*y^2
// which has one degree of freedom in the coefficients. Eliminate the degree
// of freedom by minimizing the quadratic form E(C) = C^T M C subject to
// Length(C) = 1 with M = sum_i V[i]*V[i]^T where
//   V = (1, x, y, x^2, x*y, y^2)
// The minimum value is the smallest eigenvalue of M and C is a
// corresponding unit length eigenvector.
//
// Canonical forms. The quadratic equation can be factored into
// P^T A P + B^T P + K = 0 where P = (x,y), K = C[0], B = (C[1],C[2]) and A
// is a 2x2 symmetric matrix with A00 = C[3], A01 = C[4]/2 and A11 = C[5].
// Using an eigendecomposition, matrix A = R^T D R where R is orthogonal and
// D is diagonal. Define V = R*P = (v0,v1), E = R*B = (e0,e1),
// D = diag(d0,d1) and f = K to obtain
//   d0 v0^2 + d1 v1^2 + e0 v0 + e1 v1 + f = 0
// The classification depends on the signs of the d_i.
export class ApprQuadratic2 {
    compute(points: readonly Vector[]): ApprQuadratic2Result {
        const numPoints = points.length;
        const M = new Matrix(6, 6);  // the constructor sets M to zero
        for (let i = 0; i < numPoints; ++i) {
            const x = points[i].values[0];
            const y = points[i].values[1];
            const x2 = x * x;
            const y2 = y * y;
            const xy = x * y;
            const x3 = x * x2;
            const xy2 = x * y2;
            const x2y = x * xy;
            const y3 = y * y2;
            const x4 = x * x3;
            const x2y2 = x * xy2;
            const x3y = x * x2y;
            const y4 = y * y3;
            const xy3 = x * y3;

            // M(0, 0) += 1
            M.set(0, 1, M.get(0, 1) + x);
            M.set(0, 2, M.get(0, 2) + y);
            M.set(0, 3, M.get(0, 3) + x2);
            M.set(0, 4, M.get(0, 4) + xy);
            M.set(0, 5, M.get(0, 5) + y2);

            // M(1, 1) += x2    [M(0,3)]
            // M(1, 2) += xy    [M(0,4)]
            M.set(1, 3, M.get(1, 3) + x3);
            M.set(1, 4, M.get(1, 4) + x2y);
            M.set(1, 5, M.get(1, 5) + xy2);

            // M(2, 2) += y2    [M(0,5)]
            // M(2, 3) += x2y   [M(1,4)]
            // M(2, 4) += xy2   [M(1,5)]
            M.set(2, 5, M.get(2, 5) + y3);

            M.set(3, 3, M.get(3, 3) + x4);
            M.set(3, 4, M.get(3, 4) + x3y);
            M.set(3, 5, M.get(3, 5) + x2y2);

            // M(4, 4) += x2y2  [M(3,5)]
            M.set(4, 5, M.get(4, 5) + xy3);

            M.set(5, 5, M.get(5, 5) + y4);
        }

        const rNumPoints = numPoints;
        M.set(0, 0, rNumPoints);
        M.set(1, 1, M.get(0, 3));  // x2
        M.set(1, 2, M.get(0, 4));  // xy
        M.set(2, 2, M.get(0, 5));  // y2
        M.set(2, 3, M.get(1, 4));  // x2y
        M.set(2, 4, M.get(1, 5));  // xy2
        M.set(4, 4, M.get(3, 5));  // x2y2

        for (let row = 0; row < 6; ++row) {
            for (let col = 0; col < row; ++col) {
                M.set(row, col, M.get(col, row));
            }
        }

        for (let row = 0; row < 6; ++row) {
            for (let col = 0; col < 6; ++col) {
                M.set(row, col, M.get(row, col) / rNumPoints);
            }
        }

        M.set(0, 0, 1);

        const es = new SymmetricEigensolver(6, 1024);
        es.solve(M.values, +1);
        const coefficients = es.getEigenvector(0);

        // For an exact fit, numeric round-off errors might make the minimum
        // eigenvalue just slightly negative. Return the clamped value
        // because the application might rely on the return value being
        // nonnegative.
        return {
            coefficients: coefficients,
            minEigenvalue: Math.max(es.getEigenvalue(0), 0)
        };
    }
}

// If you believe your points are nearly circular, use this class. The circle
// is of the form
//   C'[0] + C'[1]*x + C'[2]*y + C'[3]*(x^2 + y^2) = 0
// where Length(C') = 1. The coefficients used are
//   C = (C'[0] / C'[3], C'[1] / C'[3], C'[2] / C'[3])
//     = (C[0], C[1], C[2])
// so the fitted circle is
//   C[0] + C[1]*x + C[2]*y + x^2 + y^2 = 0
// The center is (xc,yc) = -(C[1],C[2])/2 and the radius is
// r = sqrt(xc * xc + yc * yc - C[0]).
export class ApprQuadraticCircle2 {
    // The fitted circle is written into 'circle'. The return value is the
    // (clamped nonnegative) minimum eigenvalue, a measure of the fit.
    compute(points: readonly Vector[], circle: Hypersphere): number {
        const numPoints = points.length;
        const M = new Matrix(4, 4);  // the constructor sets M to zero
        for (let i = 0; i < numPoints; ++i) {
            const x = points[i].values[0];
            const y = points[i].values[1];
            const x2 = x * x;
            const y2 = y * y;
            const xy = x * y;
            const r2 = x2 + y2;
            const xr2 = x * r2;
            const yr2 = y * r2;
            const r4 = r2 * r2;

            // M(0, 0) += 1
            M.set(0, 1, M.get(0, 1) + x);
            M.set(0, 2, M.get(0, 2) + y);
            M.set(0, 3, M.get(0, 3) + r2);

            M.set(1, 1, M.get(1, 1) + x2);
            M.set(1, 2, M.get(1, 2) + xy);
            M.set(1, 3, M.get(1, 3) + xr2);

            M.set(2, 2, M.get(2, 2) + y2);
            M.set(2, 3, M.get(2, 3) + yr2);

            M.set(3, 3, M.get(3, 3) + r4);
        }

        const rNumPoints = numPoints;
        M.set(0, 0, rNumPoints);

        for (let row = 0; row < 4; ++row) {
            for (let col = 0; col < row; ++col) {
                M.set(row, col, M.get(col, row));
            }
        }

        for (let row = 0; row < 4; ++row) {
            for (let col = 0; col < 4; ++col) {
                M.set(row, col, M.get(row, col) / rNumPoints);
            }
        }

        M.set(0, 0, 1);

        const es = new SymmetricEigensolver(4, 1024);
        es.solve(M.values, +1);
        const evector = es.getEigenvector(0);

        const coefficients: number[] = [0, 0, 0];
        for (let row = 0; row < 3; ++row) {
            coefficients[row] = evector[row] / evector[3];
        }

        // Clamp the radius to nonnegative values in case rounding errors
        // cause sqrRadius to be slightly negative.
        const negHalf = -0.5;
        circle.center.values[0] = negHalf * coefficients[1];
        circle.center.values[1] = negHalf * coefficients[2];
        const sqrRadius = dot(circle.center, circle.center) - coefficients[0];
        circle.radius = Math.sqrt(Math.max(sqrRadius, 0));

        // For an exact fit, numeric round-off errors might make the minimum
        // eigenvalue just slightly negative. Return the clamped value
        // because the application might rely on the return value being
        // nonnegative.
        return Math.max(es.getEigenvalue(0), 0);
    }
}
