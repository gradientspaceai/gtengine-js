// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprQuadratic3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: see ApprQuadratic2.ts for the shared conventions ('operator()'
// -> 'compute', pointer-plus-count -> array, the 'std::array' output
// reference becoming a field of the returned object). As in ApprQuadratic2,
// the header comment's M = (sum_i V_i)(sum_i V_i)^T is a typo for the
// M = sum_i V_i*V_i^T the code actually accumulates.
//
// Unlike ApprQuadratic2/ApprQuadraticCircle2/ApprQuadraticSphere3, upstream
// ApprQuadratic3 does not reset M(0,0) to 1 after the division by the number
// of points. That is harmless: M(0,0) was set to numPoints beforehand, so
// the division already leaves exactly 1 there. The port preserves the
// upstream code as written.

import { Hypersphere } from './Hypersphere.js';
import { Matrix } from './Matrix.js';
import { SymmetricEigensolver } from './SymmetricEigensolver.js';
import { Vector, dot } from './Vector.js';

export interface ApprQuadratic3Result {
    // The coefficients C[0..9] of the quadratic fit (a unit-length
    // eigenvector of M).
    coefficients: number[];

    // A nonnegative measure of the fit (the minimum eigenvalue of M);
    // 0 means an exact fit, positive otherwise.
    minEigenvalue: number;
}

// The quadratic fit is
//   0 = C[0] + C[1]*x + C[2]*y + C[3]*z + C[4]*x^2 + C[5]*x*y
//       + C[6]*x*z + C[7]*y^2 + C[8]*y*z + C[9]*z^2
// which has one degree of freedom in the coefficients. Eliminate the degree
// of freedom by minimizing the quadratic form E(C) = C^T M C subject to
// Length(C) = 1 with M = sum_i V_i*V_i^T where
//   V = (1, x, y, z, x^2, x*y, x*z, y^2, y*z, z^2)
// The minimum value is the smallest eigenvalue of M and C is a
// corresponding unit length eigenvector.
//
// Canonical forms. The quadratic equation can be factored into
// P^T A P + B^T P + K = 0 where P = (x,y,z), K = C[0], B = (C[1],C[2],C[3])
// and A is a 3x3 symmetric matrix with A00 = C[4], A01 = C[5]/2,
// A02 = C[6]/2, A11 = C[7], A12 = C[8]/2 and A22 = C[9]. Using an
// eigendecomposition, matrix A = R^T D R where R is orthogonal and D is
// diagonal. Define V = R*P = (v0,v1,v2), E = R*B = (e0,e1,e2),
// D = diag(d0,d1,d2) and f = K to obtain
//   d0 v0^2 + d1 v1^2 + d2 v2^2 + e0 v0 + e1 v1 + e2 v2 + f = 0
// The classification depends on the signs of the d_i. See QuadricSurface.ts
// for determining the type of quadric surface.
export class ApprQuadratic3 {
    compute(points: readonly Vector[]): ApprQuadratic3Result {
        const numPoints = points.length;
        const M = new Matrix(10, 10);  // the constructor sets M to zero
        for (let i = 0; i < numPoints; ++i) {
            const x = points[i].values[0];
            const y = points[i].values[1];
            const z = points[i].values[2];
            const x2 = x * x;
            const y2 = y * y;
            const z2 = z * z;
            const xy = x * y;
            const xz = x * z;
            const yz = y * z;
            const x3 = x * x2;
            const xy2 = x * y2;
            const xz2 = x * z2;
            const x2y = x * xy;
            const x2z = x * xz;
            const xyz = x * yz;
            const y3 = y * y2;
            const yz2 = y * z2;
            const y2z = y * yz;
            const z3 = z * z2;
            const x4 = x * x3;
            const x2y2 = x * xy2;
            const x2z2 = x * xz2;
            const x3y = x * x2y;
            const x3z = x * x2z;
            const x2yz = x * xyz;
            const y4 = y * y3;
            const y2z2 = y * yz2;
            const xy3 = x * y3;
            const xy2z = x * y2z;
            const y3z = y * y2z;
            const z4 = z * z3;
            const xyz2 = x * yz2;
            const xz3 = x * z3;
            const yz3 = y * z3;

            // M(0, 0) += 1
            M.set(0, 1, M.get(0, 1) + x);
            M.set(0, 2, M.get(0, 2) + y);
            M.set(0, 3, M.get(0, 3) + z);
            M.set(0, 4, M.get(0, 4) + x2);
            M.set(0, 5, M.get(0, 5) + xy);
            M.set(0, 6, M.get(0, 6) + xz);
            M.set(0, 7, M.get(0, 7) + y2);
            M.set(0, 8, M.get(0, 8) + yz);
            M.set(0, 9, M.get(0, 9) + z2);

            // M(1, 1) += x2    [M(0,4)]
            // M(1, 2) += xy    [M(0,5)]
            // M(1, 3) += xz    [M(0,6)]
            M.set(1, 4, M.get(1, 4) + x3);
            M.set(1, 5, M.get(1, 5) + x2y);
            M.set(1, 6, M.get(1, 6) + x2z);
            M.set(1, 7, M.get(1, 7) + xy2);
            M.set(1, 8, M.get(1, 8) + xyz);
            M.set(1, 9, M.get(1, 9) + xz2);

            // M(2, 2) += y2    [M(0,7)]
            // M(2, 3) += yz    [M(0,8)]
            // M(2, 4) += x2y   [M(1,5)]
            M.set(2, 5, M.get(2, 5) + xy2);
            // M(2, 6) += xyz   [M(1,8)]
            M.set(2, 7, M.get(2, 7) + y3);
            M.set(2, 8, M.get(2, 8) + y2z);
            M.set(2, 9, M.get(2, 9) + yz2);

            // M(3, 3) += z2    [M(0,9)]
            // M(3, 4) += x2z   [M(1,6)]
            // M(3, 5) += xyz   [M(1,8)]
            // M(3, 6) += xz2   [M(1,9)]
            // M(3, 7) += y2z   [M(2,8)]
            // M(3, 8) += yz2   [M(2,9)]
            M.set(3, 9, M.get(3, 9) + z3);

            M.set(4, 4, M.get(4, 4) + x4);
            M.set(4, 5, M.get(4, 5) + x3y);
            M.set(4, 6, M.get(4, 6) + x3z);
            M.set(4, 7, M.get(4, 7) + x2y2);
            M.set(4, 8, M.get(4, 8) + x2yz);
            M.set(4, 9, M.get(4, 9) + x2z2);

            // M(5, 5) += x2y2  [M(4,7)]
            // M(5, 6) += x2yz  [M(4,8)]
            M.set(5, 7, M.get(5, 7) + xy3);
            M.set(5, 8, M.get(5, 8) + xy2z);
            M.set(5, 9, M.get(5, 9) + xyz2);

            // M(6, 6) += x2z2  [M(4,9)]
            // M(6, 7) += xy2z  [M(5,8)]
            // M(6, 8) += xyz2  [M(5,9)]
            M.set(6, 9, M.get(6, 9) + xz3);

            M.set(7, 7, M.get(7, 7) + y4);
            M.set(7, 8, M.get(7, 8) + y3z);
            M.set(7, 9, M.get(7, 9) + y2z2);

            // M(8, 8) += y2z2  [M(7,9)]
            M.set(8, 9, M.get(8, 9) + yz3);

            M.set(9, 9, M.get(9, 9) + z4);
        }

        const rNumPoints = numPoints;
        M.set(0, 0, rNumPoints);
        M.set(1, 1, M.get(0, 4));  // x2
        M.set(1, 2, M.get(0, 5));  // xy
        M.set(1, 3, M.get(0, 6));  // xz
        M.set(2, 2, M.get(0, 7));  // y2
        M.set(2, 3, M.get(0, 8));  // yz
        M.set(2, 4, M.get(1, 5));  // x2y
        M.set(2, 6, M.get(1, 8));  // xyz
        M.set(3, 3, M.get(0, 9));  // z2
        M.set(3, 4, M.get(1, 6));  // x2z
        M.set(3, 5, M.get(1, 8));  // xyz
        M.set(3, 6, M.get(1, 9));  // xz2
        M.set(3, 7, M.get(2, 8));  // y2z
        M.set(3, 8, M.get(2, 9));  // yz2
        M.set(5, 5, M.get(4, 7));  // x2y2
        M.set(5, 6, M.get(4, 8));  // x2yz
        M.set(6, 6, M.get(4, 9));  // x2z2
        M.set(6, 7, M.get(5, 8));  // xy2z
        M.set(6, 8, M.get(5, 9));  // xyz2
        M.set(8, 8, M.get(7, 9));  // y2z2

        for (let row = 0; row < 10; ++row) {
            for (let col = 0; col < row; ++col) {
                M.set(row, col, M.get(col, row));
            }
        }

        for (let row = 0; row < 10; ++row) {
            for (let col = 0; col < 10; ++col) {
                M.set(row, col, M.get(row, col) / rNumPoints);
            }
        }

        const es = new SymmetricEigensolver(10, 1024);
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

// If you believe your points are nearly spherical, use this class. The
// sphere is of the form
//   C'[0] + C'[1]*x + C'[2]*y + C'[3]*z + C'[4]*(x^2 + y^2 + z^2) = 0
// where Length(C') = 1. The coefficients used are
//   C = (C'[0] / C'[4], C'[1] / C'[4], C'[2] / C'[4], C'[3] / C'[4])
//     = (C[0], C[1], C[2], C[3])
// so the fitted sphere is
//   C[0] + C[1]*x + C[2]*y + C[3]*z + x^2 + y^2 + z^2 = 0
// The center is (xc,yc,zc) = -0.5*(C[1],C[2],C[3]) and the radius is
// r = sqrt(xc * xc + yc * yc + zc * zc - C[0]).
export class ApprQuadraticSphere3 {
    // The fitted sphere is written into 'sphere'. The return value is the
    // (clamped nonnegative) minimum eigenvalue, a measure of the fit.
    compute(points: readonly Vector[], sphere: Hypersphere): number {
        const numPoints = points.length;
        const M = new Matrix(5, 5);  // the constructor sets M to zero
        for (let i = 0; i < numPoints; ++i) {
            const x = points[i].values[0];
            const y = points[i].values[1];
            const z = points[i].values[2];
            const x2 = x * x;
            const y2 = y * y;
            const z2 = z * z;
            const xy = x * y;
            const xz = x * z;
            const yz = y * z;
            const r2 = x2 + y2 + z2;
            const xr2 = x * r2;
            const yr2 = y * r2;
            const zr2 = z * r2;
            const r4 = r2 * r2;

            // M(0, 0) += 1
            M.set(0, 1, M.get(0, 1) + x);
            M.set(0, 2, M.get(0, 2) + y);
            M.set(0, 3, M.get(0, 3) + z);
            M.set(0, 4, M.get(0, 4) + r2);

            M.set(1, 1, M.get(1, 1) + x2);
            M.set(1, 2, M.get(1, 2) + xy);
            M.set(1, 3, M.get(1, 3) + xz);
            M.set(1, 4, M.get(1, 4) + xr2);

            M.set(2, 2, M.get(2, 2) + y2);
            M.set(2, 3, M.get(2, 3) + yz);
            M.set(2, 4, M.get(2, 4) + yr2);

            M.set(3, 3, M.get(3, 3) + z2);
            M.set(3, 4, M.get(3, 4) + zr2);

            M.set(4, 4, M.get(4, 4) + r4);
        }

        const rNumPoints = numPoints;
        M.set(0, 0, rNumPoints);

        for (let row = 0; row < 5; ++row) {
            for (let col = 0; col < row; ++col) {
                M.set(row, col, M.get(col, row));
            }
        }

        for (let row = 0; row < 5; ++row) {
            for (let col = 0; col < 5; ++col) {
                M.set(row, col, M.get(row, col) / rNumPoints);
            }
        }

        M.set(0, 0, 1);

        const es = new SymmetricEigensolver(5, 1024);
        es.solve(M.values, +1);
        const evector = es.getEigenvector(0);

        const coefficients: number[] = [0, 0, 0, 0];
        for (let row = 0; row < 4; ++row) {
            coefficients[row] = evector[row] / evector[4];
        }

        // Clamp the radius to nonnegative values in case rounding errors
        // cause sqrRadius to be slightly negative.
        const negHalf = -0.5;
        sphere.center.values[0] = negHalf * coefficients[1];
        sphere.center.values[1] = negHalf * coefficients[2];
        sphere.center.values[2] = negHalf * coefficients[3];
        const sqrRadius = dot(sphere.center, sphere.center) - coefficients[0];
        sphere.radius = Math.sqrt(Math.max(sqrRadius, 0));

        // For an exact fit, numeric round-off errors might make the minimum
        // eigenvalue just slightly negative. Return the clamped value
        // because the application might rely on the return value being
        // nonnegative.
        return Math.max(es.getEigenvalue(0), 0);
    }
}
