// gtengine-js: TypeScript port of Geometric Tools Engine (GTE)
// LieGroupsAlgebras.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Lie groups and Lie algebras are useful for representing special classes of
// matrices found in applications. Implementations are provided for rotations
// in 2D and 3D and for rigid motions (rotations and translations) in 2D and
// 3D. The mathematical details are found in
//   https://www.geometrictools.com/Documentation/LieGroupsAlgebras.pdf
//
// Port notes:
// - Upstream's 'using AlgebraType/AdjointType/GroupType' member aliases
//   become exported type aliases named '<Class><Alias>'; they are doc-only
//   because the port uses runtime-dimension Vector/Matrix.
// - The upstream classes have only static members, so the port keeps them as
//   classes with static methods. No free functions are exported, which keeps
//   the library-wide flat export free of collisions.
// - The private F0/F1/F2 helpers are duplicated verbatim in each upstream
//   class; the port has a single module-private copy of each. The
//   small-angle thresholds (|t| > 1/16) and the degree-16 RotC*Estimate
//   polynomials are transcribed unchanged.
// - The two 'GeodesicPath' overloads are distinguished at runtime by the
//   type of the third argument (the group element M1 versus the precomputed
//   Lie algebra element log(M1*Inverse(M0))).
// - Only the GTE_USE_MAT_VEC branches are ported, matching the rest of the
//   library.

import {
    Matrix, addMatrix, mulMatrix, multiplyAB, multiplyABT, hliftMatrix,
    hprojectMatrix
} from './Matrix.js';
import { inverse2x2 } from './Matrix2x2.js';
import { inverse3x3, trace3x3 } from './Matrix3x3.js';
import { inverse4x4 } from './Matrix4x4.js';
import { GTE_C_PI } from './Constants.js';
import { logAssert } from './Logger.js';
import {
    rotC0Estimate, rotC1Estimate, rotC4Estimate
} from './RotationEstimate.js';
import { Vector, dot, mul, normalize } from './Vector.js';

// Compute sin(t)/t. For small |t| the closed form loses significant digits,
// so a degree-16 polynomial estimate of sin(t)/t is used instead.
function f0(t: number): number {
    if (Math.abs(t) > 0.0625) {
        return Math.sin(t) / t;
    }
    return rotC0Estimate(t, 16);
}

// Compute (1 - cos(t))/t^2, with the small-|t| estimate as in f0.
function f1(t: number): number {
    if (Math.abs(t) > 0.0625) {
        return (1 - Math.cos(t)) / t / t;
    }
    return rotC1Estimate(t, 16);
}

// Compute (t - sin(t))/t^3, with the small-|t| estimate as in f0.
function f2(t: number): number {
    if (Math.abs(t) > 0.0625) {
        return (t - Math.sin(t)) / t / t / t;
    }
    return rotC4Estimate(t, 16);
}

function assertSquare(M: Matrix, n: number, name: string): void {
    logAssert(M.numRows === n && M.numCols === n,
        `${name}: expecting a ${n}x${n} matrix.`);
}

function assertSize(v: Vector, n: number, name: string): void {
    logAssert(v.size === n, `${name}: expecting a ${n}-tuple.`);
}

// SO(2) is the Lie group for rotations in 2D. so(2) is the corresponding
// Lie algebra for SO(2) and is a 1D quantity x = (angle). The 2x2 rotation
// matrix M is generated from x by constructing a 2x2 generator G = x*G0,
// where
//   G0 = {{ 0, -1 },{ 1, 0 }}
// and then computing the power series M = exp(L(x)). For the sake of
// notation, exp(x) is used to denote exp(L(x)). The 2x2 rotation matrix is
//   M = {{ cos(x), -sin(x) }, { sin(x), cos(x) }}
// The adjoint matrix is the 1x1 identity matrix
//   A(M) = 1

// n = 2, k = 1, x = (angle)
export type LieSO2AlgebraType = number;  // kx1
export type LieSO2AdjointType = number;  // kxk
export type LieSO2GroupType = Matrix;    // nxn (2x2)

export class LieSO2 {
    // Compute the Lie group element X from the Lie algebra element x using
    // X = L(x).
    static toGroup(x: LieSO2AlgebraType): Matrix {
        const X = new Matrix(2, 2);
        X.set(0, 0, 0);
        X.set(0, 1, -x);
        X.set(1, 0, x);
        X.set(1, 1, 0);
        return X;
    }

    // Compute the Lie algebra element x from the Lie group element X using
    // x = L^{-1}(X).
    static toAlgebra(X: Matrix): LieSO2AlgebraType {
        assertSquare(X, 2, 'LieSO2');
        return X.get(1, 0);
    }

    // Compute the exponential map of the Lie algebra element x to produce
    // the Lie group element Y = exp(X) = exp(L(x)).
    static exp(x: LieSO2AlgebraType): Matrix {
        const Y = new Matrix(2, 2);
        const sn = Math.sin(x);
        const cs = Math.cos(x);
        Y.set(0, 0, cs);
        Y.set(0, 1, -sn);
        Y.set(1, 0, sn);
        Y.set(1, 1, cs);
        return Y;
    }

    // Compute the logarithm map of the Lie group element Y to produce the
    // Lie algebra element x corresponding to the Lie group element X.
    static log(Y: Matrix): LieSO2AlgebraType {
        assertSquare(Y, 2, 'LieSO2');
        return Math.atan2(Y.get(1, 0), Y.get(0, 0));
    }

    // Compute the adjoint matrix A(M) from the Lie group element M.
    static adjoint(_M: Matrix): LieSO2AdjointType {
        return 1;
    }

    // Compute log(M1*Inverse(M0)) to reduce computation time when you want
    // to evaluate geodesicPath for multiple values of t. For a rotation
    // matrix M0, Inverse(M0) is equal to Transpose(M0), which avoids a
    // general inversion of M0.
    static logM1M0Inv(M0: Matrix, M1: Matrix): LieSO2AlgebraType {
        assertSquare(M0, 2, 'LieSO2');
        assertSquare(M1, 2, 'LieSO2');
        return LieSO2.log(multiplyABT(M1, M0));
    }

    // Compute a point on the geodesic path from M0 to M1. When the third
    // argument is the group element M1, log(M1*Inverse(M0)) is computed for
    // each call; use that form for a single value of t. When the third
    // argument is the precomputed Lie algebra element log(M1*Inverse(M0)),
    // no logarithm is computed; use that form for multiple values of t.
    static geodesicPath(t: number, M0: Matrix, M1: Matrix): Matrix;
    static geodesicPath(t: number, M0: Matrix,
        logM1M0Inv: LieSO2AlgebraType): Matrix;
    static geodesicPath(t: number, M0: Matrix,
        arg: Matrix | LieSO2AlgebraType): Matrix {
        const logValue = (typeof arg === 'number'
            ? arg : LieSO2.logM1M0Inv(M0, arg));
        return multiplyAB(LieSO2.exp(t * logValue), M0);
    }
}

// SE(2) is the Lie group for rigid motions in 2D. se(2) is the
// corresponding Lie algebra for SE(2) and is a 3D quantity
// x = (angle; u0, u1), where (angle) is for the rotation matrix and
// (u0, u1) is for the translation vector. The 3x3 rigid motion M is
// generated from x by constructing a 3x3 generator
// G = x0 * G0 + x1 * G1 + x2 * G2, where
//   G0 = {{ 0, -1, 0 }, { 1, 0, 0 }, { 0, 0, 0 }}
//   G1 = {{ 0, 0, 1 }, { 0, 0, 0 }, { 0, 0, 0 }}
//   G2 = {{ 0, 0, 0 }, { 0, 0, 1 }, { 0, 0, 0 }}
// and then computing the power series M = exp(L(x)). The rigid motion
// matrix is
//   M = {{ R, T }, { 0, 1 }}
// where R is the 2x2 rotation matrix, T is the 2x1 translation vector, 0 is
// the 1x2 zero vector and 1 is a scalar. The adjoint matrix is
//   A(M) = {{ 1, 0 }, { Perp(T), R }}
// where T = (t0,t1) and Perp(T) = (t1,-t0), both 2x1 vectors but written as
// 2-tuples.

// n = 3, k = 3, x = (angle; u0, u1)
export type LieSE2AlgebraType = Vector;   // kx1 (3-tuple)
export type LieSE2AdjointType = Matrix;   // kxk (3x3)
export type LieSE2GroupType = Matrix;     // nxn (3x3)

export class LieSE2 {
    // Compute the Lie group element X from the Lie algebra element x using
    // X = L(x).
    static toGroup(x: Vector): Matrix {
        assertSize(x, 3, 'LieSE2');
        const X = new Matrix(3, 3);
        X.set(0, 0, 0);
        X.set(0, 1, -x.values[0]);
        X.set(0, 2, x.values[1]);
        X.set(1, 0, x.values[0]);
        X.set(1, 1, 0);
        X.set(1, 2, x.values[2]);
        X.set(2, 0, 0);
        X.set(2, 1, 0);
        X.set(2, 2, 0);
        return X;
    }

    // Compute the Lie algebra element x from the Lie group element X using
    // x = L^{-1}(X).
    static toAlgebra(X: Matrix): Vector {
        assertSquare(X, 3, 'LieSE2');
        return Vector.fromArray([X.get(1, 0), X.get(0, 2), X.get(1, 2)]);
    }

    // Compute the exponential map of the Lie algebra element x to produce
    // the Lie group element Y = exp(X) = exp(L(x)).
    static exp(x: Vector): Matrix {
        assertSize(x, 3, 'LieSE2');
        const sn = Math.sin(x.values[0]);
        const cs = Math.cos(x.values[0]);

        // Compute sin(t)/t.
        const a0 = f0(x.values[0]);

        // Compute (1 - cos(t))/t = t * (1-cos(t))/t^2.
        const a1 = x.values[0] * f1(x.values[0]);

        const trn0 = a0 * x.values[1] - a1 * x.values[2];
        const trn1 = a1 * x.values[1] + a0 * x.values[2];

        const Y = new Matrix(3, 3);
        Y.set(0, 0, cs);
        Y.set(0, 1, -sn);
        Y.set(0, 2, trn0);
        Y.set(1, 0, sn);
        Y.set(1, 1, cs);
        Y.set(1, 2, trn1);
        Y.set(2, 0, 0);
        Y.set(2, 1, 0);
        Y.set(2, 2, 1);
        return Y;
    }

    // Compute the logarithm map of the Lie group element Y to produce the
    // Lie algebra element x corresponding to the Lie group element X.
    static log(Y: Matrix): Vector {
        assertSquare(Y, 3, 'LieSE2');
        const x = new Vector(3);
        x.values[0] = Math.atan2(Y.get(1, 0), Y.get(0, 0));

        // Compute sin(t)/t.
        const a0 = f0(x.values[0]);

        // Compute (1 - cos(t))/t = t * (1-cos(t))/t^2.
        const a1 = x.values[0] * f1(x.values[0]);

        const V = new Matrix(2, 2);
        V.set(0, 0, a0);
        V.set(1, 0, a1);
        V.set(0, 1, -a1);
        V.set(1, 1, a0);

        const inverseV = inverse2x2(V).inverse;
        x.values[1] = inverseV.get(0, 0) * Y.get(0, 2)
            + inverseV.get(0, 1) * Y.get(1, 2);
        x.values[2] = inverseV.get(1, 0) * Y.get(0, 2)
            + inverseV.get(1, 1) * Y.get(1, 2);
        return x;
    }

    // Compute the adjoint matrix A(M) from the Lie group element M.
    static adjoint(M: Matrix): Matrix {
        assertSquare(M, 3, 'LieSE2');
        const A = new Matrix(3, 3);
        A.set(0, 0, 1);
        A.set(0, 1, 0);
        A.set(0, 2, 0);
        A.set(1, 0, M.get(1, 2));
        A.set(1, 1, M.get(0, 0));
        A.set(1, 2, M.get(0, 1));
        A.set(2, 0, -M.get(0, 2));
        A.set(2, 1, M.get(1, 0));
        A.set(2, 2, M.get(1, 1));
        return A;
    }

    // Compute log(M1*Inverse(M0)) to reduce computation time when you want
    // to evaluate geodesicPath for multiple values of t.
    static logM1M0Inv(M0: Matrix, M1: Matrix): Vector {
        assertSquare(M0, 3, 'LieSE2');
        assertSquare(M1, 3, 'LieSE2');
        return LieSE2.log(multiplyAB(M1, inverse3x3(M0).inverse));
    }

    // Compute a point on the geodesic path from M0 to M1. See the comment
    // for LieSO2.geodesicPath about the two forms of the third argument.
    static geodesicPath(t: number, M0: Matrix, M1: Matrix): Matrix;
    static geodesicPath(t: number, M0: Matrix, logM1M0Inv: Vector): Matrix;
    static geodesicPath(t: number, M0: Matrix, arg: Matrix | Vector): Matrix {
        const logValue = (arg instanceof Matrix
            ? LieSE2.logM1M0Inv(M0, arg) : arg);
        return multiplyAB(LieSE2.exp(mul(logValue, t)), M0);
    }
}

// SO(3) is the Lie group for rotations in 3D. so(3) is the corresponding
// Lie algebra for SO(3) and is a 3D quantity x = (x0,x1,x2). The 3x3
// rotation matrix M is generated from x by constructing a 3x3 generator
// G = x0*G0 + x1*G1 + x2*G2, where
//   G0 = {{ 0, 0, 0 }, { 0, 0, -1 }, { 0, 1, 0 }}
//   G1 = {{ 0, 0, 1 }, { 0, 0, 0 }, { -1, 0, 0 }}
//   G2 = {{ 0, -1, 0 }, { 1, 0, 0 }, { 0, 0, 0 }}
// and then computing the power series M = exp(L(x)). The rotation matrix is
//   M = I + (sin(angle)/angle) * S + ((1 - cos(angle))/angle^2) * S^2
// where angle is the length of x. The adjoint matrix is
//   A(M) = M

// n = 3, k = 3, x = (s0,s1,s2)
export type LieSO3AlgebraType = Vector;   // kx1 (3-tuple)
export type LieSO3AdjointType = Matrix;   // kxk (3x3)
export type LieSO3GroupType = Matrix;     // nxn (3x3)

export class LieSO3 {
    // Compute the Lie group element X from the Lie algebra element x using
    // X = L(x).
    static toGroup(x: Vector): Matrix {
        assertSize(x, 3, 'LieSO3');
        const X = new Matrix(3, 3);
        X.set(0, 0, 0);
        X.set(0, 1, -x.values[2]);
        X.set(0, 2, x.values[1]);
        X.set(1, 0, x.values[2]);
        X.set(1, 1, 0);
        X.set(1, 2, -x.values[0]);
        X.set(2, 0, -x.values[1]);
        X.set(2, 1, x.values[0]);
        X.set(2, 2, 0);
        return X;
    }

    // Compute the Lie algebra element x from the Lie group element X using
    // x = L^{-1}(X).
    static toAlgebra(X: Matrix): Vector {
        assertSquare(X, 3, 'LieSO3');
        return Vector.fromArray([X.get(2, 1), X.get(0, 2), X.get(1, 0)]);
    }

    // Compute the exponential map of the Lie algebra element x to produce
    // the Lie group element Y = exp(X) = exp(L(x)).
    static exp(x: Vector): Matrix {
        assertSize(x, 3, 'LieSO3');
        const sqrAngle = dot(x, x);
        const angle = Math.sqrt(sqrAngle);
        if (angle > 0) {
            const Y = LieSO3.toGroup(x);
            const Ysqr = multiplyAB(Y, Y);

            // Compute sin(t)/t.
            const a0 = f0(angle);

            // Compute (1 - cos(t))/t^2.
            const a1 = f1(angle);

            return addMatrix(addMatrix(Matrix.identity(3, 3),
                mulMatrix(Y, a0)), mulMatrix(Ysqr, a1));
        }
        return Matrix.identity(3, 3);
    }

    // Compute the logarithm map of the Lie group element Y to produce the
    // Lie algebra element x corresponding to the Lie group element X.
    static log(Y: Matrix): Vector {
        assertSquare(Y, 3, 'LieSO3');
        const x = new Vector(3);

        const arg = 0.5 * (trace3x3(Y) - 1);  // in [-1,1]
        if (arg > -1) {
            if (arg < 1) {
                // 0 < angle < pi
                const angle = Math.acos(arg);
                // G = (angle / (2*sin(angle)) * (Y - Y^T)
                const multiplier = 0.5 / f0(angle);
                x.values[0] = multiplier * (Y.get(2, 1) - Y.get(1, 2));
                x.values[1] = multiplier * (Y.get(0, 2) - Y.get(2, 0));
                x.values[2] = multiplier * (Y.get(1, 0) - Y.get(0, 1));
            } else {
                // arg = 1, angle = 0, Y is the identity, G is zero
                x.makeZero();
            }
        } else {
            // arg = -1, angle = pi. Observe that Y = I + (2/pi^2) * G^2.
            // Consider x as a 3x1 vector; then x * x^T = (pi^2/2)*(Y + I).
            // The right-hand side is a symmetric matrix with positive
            // diagonal entries and rank 1. Choose the row of Y + I that has
            // the largest diagonal term and normalize that row. Multiply it
            // by pi/sqrt(2) to obtain x from which G = toGroup(x). The
            // vector -x is also a candidate but irrelevant here because x
            // and -x produce the same rotation matrix. Knowing Y+I is
            // symmetric and wanting to avoid bias, use (Y(i,j)+Y(j,i))/2 for
            // the off-diagonal entries rather than Y(i,j).
            if (Y.get(0, 0) >= Y.get(1, 1)) {
                if (Y.get(0, 0) >= Y.get(2, 2)) {
                    // r00 is maximum diagonal term
                    x.values[0] = Y.get(0, 0) + 1;
                    x.values[1] = 0.5 * (Y.get(0, 1) + Y.get(1, 0));
                    x.values[2] = 0.5 * (Y.get(0, 2) + Y.get(2, 0));
                } else {
                    // r22 is maximum diagonal term
                    x.values[0] = 0.5 * (Y.get(2, 0) + Y.get(0, 2));
                    x.values[1] = 0.5 * (Y.get(2, 1) + Y.get(1, 2));
                    x.values[2] = Y.get(2, 2) + 1;
                }
            } else {
                if (Y.get(1, 1) >= Y.get(2, 2)) {
                    // r11 is maximum diagonal term
                    x.values[0] = 0.5 * (Y.get(1, 0) + Y.get(0, 1));
                    x.values[1] = Y.get(1, 1) + 1;
                    x.values[2] = 0.5 * (Y.get(1, 2) + Y.get(2, 1));
                } else {
                    // r22 is maximum diagonal term
                    x.values[0] = 0.5 * (Y.get(2, 0) + Y.get(0, 2));
                    x.values[1] = 0.5 * (Y.get(2, 1) + Y.get(1, 2));
                    x.values[2] = Y.get(2, 2) + 1;
                }
            }

            if (normalize(x) > 0) {
                // Upstream bug (fixed here): upstream scales the normalized
                // row by GTE_C_PI * GTE_C_INV_SQRT_2. Row i of Y + I equals
                // (2/pi^2) * x[i] * x, so normalizing it gives the unit
                // vector x/|x| and the correct scale factor is |x| = pi.
                // Upstream's extra 1/sqrt(2) returns a Lie algebra element
                // of length pi/sqrt(2) = 2.2214..., so exp(log(Y)) != Y for
                // every rotation by pi.
                const angle = GTE_C_PI;
                for (let i = 0; i < 3; ++i) {
                    x.values[i] *= angle;
                }
            } else {
                x.makeZero();
            }
        }

        return x;
    }

    // Compute the adjoint matrix A(M) from the Lie group element M.
    static adjoint(M: Matrix): Matrix {
        assertSquare(M, 3, 'LieSO3');
        return M.clone();
    }

    // Compute log(M1*Inverse(M0)) to reduce computation time when you want
    // to evaluate geodesicPath for multiple values of t. For a rotation
    // matrix M0, Inverse(M0) is equal to Transpose(M0), which avoids a
    // general inversion of M0.
    static logM1M0Inv(M0: Matrix, M1: Matrix): Vector {
        assertSquare(M0, 3, 'LieSO3');
        assertSquare(M1, 3, 'LieSO3');
        return LieSO3.log(multiplyABT(M1, M0));
    }

    // Compute a point on the geodesic path from M0 to M1. See the comment
    // for LieSO2.geodesicPath about the two forms of the third argument.
    static geodesicPath(t: number, M0: Matrix, M1: Matrix): Matrix;
    static geodesicPath(t: number, M0: Matrix, logM1M0Inv: Vector): Matrix;
    static geodesicPath(t: number, M0: Matrix, arg: Matrix | Vector): Matrix {
        const logValue = (arg instanceof Matrix
            ? LieSO3.logM1M0Inv(M0, arg) : arg);
        return multiplyAB(LieSO3.exp(mul(logValue, t)), M0);
    }
}

// SE(3) is the Lie group for rigid motions in 3D. se(3) is the
// corresponding Lie algebra for SE(3) and is a 6D quantity
// x = (s0,s1,s2;u0,u1,u2), where (s0,s1,s2) corresponds to the rotation
// matrix and (u0,u1,u2) corresponds to the translation vector. The 4x4
// rigid motion is generated from x by constructing a 4x4 generator
// G = x0 * G0 + x1 * G1 + x2 * G2 + x3 * G3 + x4 * G4 + x5 * G5, where
//   G0 = {{0,0,0,0},{0,0,-1,0},{0,1,0,0},{0,0,0,0}}
//   G1 = {{0,0,1,0},{0,0,0,0},{-1,0,0,0},{0,0,0,0}}
//   G2 = {{0,-1,0,0},{1,0,0,0},{0,0,0,0},{0,0,0,0}}
//   G3 = {{0,0,0,1},{0,0,0,0},{0,0,0,0},{0,0,0,0}}
//   G4 = {{0,0,0,0},{0,0,0,1},{0,0,0,0},{0,0,0,0}}
//   G5 = {{0,0,0,0},{0,0,0,0},{0,0,0,1},{0,0,0,0}}
// and then computing the power series M = exp(L(x)). The rigid motion
// matrix is
//   M = {{ R, T }, { 0, 1 }}
// where R is the 3x3 rotation matrix, T is the 3x1 translation vector, 0 is
// the 1x3 zero vector and 1 is a scalar. The adjoint matrix is
//   A(M) = {{ R, 0 }, { Skew(T)*R, R }}
// where Skew(T) = {{ 0, -T2, T1 }, { T2, 0, -T0 }, { -T1, T0, 0 }} and 0 is
// the 3x3 zero matrix.

// n = 4, k = 6, x = (s0,s1,s2,u0,u1,u2)
export type LieSE3AlgebraType = Vector;   // kx1 (6-tuple)
export type LieSE3AdjointType = Matrix;   // kxk (6x6)
export type LieSE3GroupType = Matrix;     // nxn (4x4)

export class LieSE3 {
    // Compute the Lie group element X from the Lie algebra element x using
    // X = L(x).
    static toGroup(x: Vector): Matrix {
        assertSize(x, 6, 'LieSE3');
        const X = new Matrix(4, 4);
        X.set(0, 0, 0);
        X.set(0, 1, -x.values[2]);
        X.set(0, 2, x.values[1]);
        X.set(0, 3, x.values[3]);
        X.set(1, 0, x.values[2]);
        X.set(1, 1, 0);
        X.set(1, 2, -x.values[0]);
        X.set(1, 3, x.values[4]);
        X.set(2, 0, -x.values[1]);
        X.set(2, 1, x.values[0]);
        X.set(2, 2, 0);
        X.set(2, 3, x.values[5]);
        X.set(3, 0, 0);
        X.set(3, 1, 0);
        X.set(3, 2, 0);
        X.set(3, 3, 0);
        return X;
    }

    // Compute the Lie algebra element x from the Lie group element X using
    // x = L^{-1}(X).
    static toAlgebra(X: Matrix): Vector {
        assertSquare(X, 4, 'LieSE3');
        return Vector.fromArray([
            X.get(2, 1), X.get(0, 2), X.get(1, 0),
            X.get(0, 3), X.get(1, 3), X.get(2, 3)
        ]);
    }

    // Compute the exponential map of the Lie algebra element x to produce
    // the Lie group element Y = exp(X) = exp(L(x)).
    static exp(x: Vector): Matrix {
        assertSize(x, 6, 'LieSE3');
        const s = Vector.fromArray([x.values[0], x.values[1], x.values[2]]);
        const u = Vector.fromArray([x.values[3], x.values[4], x.values[5]]);
        const S = LieSO3.toGroup(s);
        const Ssqr = multiplyAB(S, S);
        const sqrAngle = dot(s, s);
        const angle = Math.sqrt(sqrAngle);
        const a0 = f0(angle);
        const a1 = f1(angle);
        const a2 = f2(angle);
        const I3 = Matrix.identity(3, 3);
        const R = addMatrix(addMatrix(I3, mulMatrix(S, a0)),
            mulMatrix(Ssqr, a1));
        const V = addMatrix(addMatrix(I3, mulMatrix(S, a1)),
            mulMatrix(Ssqr, a2));
        const trn = mulMatrix(V, u);
        const Y = hliftMatrix(R);
        Y.set(0, 3, trn.values[0]);
        Y.set(1, 3, trn.values[1]);
        Y.set(2, 3, trn.values[2]);
        return Y;
    }

    // Compute the logarithm map of the Lie group element Y to produce the
    // Lie algebra element x corresponding to the Lie group element X.
    static log(Y: Matrix): Vector {
        assertSquare(Y, 4, 'LieSE3');
        const R = hprojectMatrix(Y);
        const s = LieSO3.log(R);
        const S = LieSO3.toGroup(s);
        const Ssqr = multiplyAB(S, S);
        const sqrAngle = dot(s, s);
        const angle = Math.sqrt(sqrAngle);
        const a1 = f1(angle);
        const a2 = f2(angle);
        const I3 = Matrix.identity(3, 3);
        const V = addMatrix(addMatrix(I3, mulMatrix(S, a1)),
            mulMatrix(Ssqr, a2));
        const inverseV = inverse3x3(V).inverse;
        const trn = Vector.fromArray([Y.get(0, 3), Y.get(1, 3), Y.get(2, 3)]);
        const u = mulMatrix(inverseV, trn);
        return Vector.fromArray([
            s.values[0], s.values[1], s.values[2],
            u.values[0], u.values[1], u.values[2]
        ]);
    }

    // Compute the adjoint matrix A(M) from the Lie group element M.
    static adjoint(M: Matrix): Matrix {
        assertSquare(M, 4, 'LieSE3');
        const R = hprojectMatrix(M);
        const skewT = LieSO3.toGroup(Vector.fromArray(
            [M.get(0, 3), M.get(1, 3), M.get(2, 3)]));
        const skewTR = multiplyAB(skewT, R);
        const A = new Matrix(6, 6);
        for (let row = 0, rowp3 = 3; row < 3; ++row, ++rowp3) {
            for (let col = 0, colp3 = 3; col < 3; ++col, ++colp3) {
                A.set(row, col, M.get(row, col));
                A.set(row, colp3, 0);
                A.set(rowp3, col, skewTR.get(row, col));
                A.set(rowp3, colp3, M.get(row, col));
            }
        }
        return A;
    }

    // Compute log(M1*Inverse(M0)) to reduce computation time when you want
    // to evaluate geodesicPath for multiple values of t.
    static logM1M0Inv(M0: Matrix, M1: Matrix): Vector {
        assertSquare(M0, 4, 'LieSE3');
        assertSquare(M1, 4, 'LieSE3');
        return LieSE3.log(multiplyAB(M1, inverse4x4(M0).inverse));
    }

    // Compute a point on the geodesic path from M0 to M1. See the comment
    // for LieSO2.geodesicPath about the two forms of the third argument.
    static geodesicPath(t: number, M0: Matrix, M1: Matrix): Matrix;
    static geodesicPath(t: number, M0: Matrix, logM1M0Inv: Vector): Matrix;
    static geodesicPath(t: number, M0: Matrix, arg: Matrix | Vector): Matrix {
        const logValue = (arg instanceof Matrix
            ? LieSE3.logM1M0Inv(M0, arg) : arg);
        return multiplyAB(LieSE3.exp(mul(logValue, t)), M0);
    }
}
