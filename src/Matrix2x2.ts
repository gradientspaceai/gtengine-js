// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Matrix2x2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: upstream's 'template <typename Real> using Matrix2x2 =
// Matrix<2, 2, Real>' is a compile-time alias. The port's Matrix has runtime
// dimensions, so Matrix2x2 is a plain type alias for Matrix and the free
// functions of this file assert that their inputs are 2-by-2.
//
// Fixed-size naming scheme (set by this batch, B56): the free functions of
// Matrix2x2.h, Matrix3x3.h and Matrix4x4.h have identical upstream names
// (Inverse, Adjoint, Determinant, Trace, DoTransform, SetBasis, GetBasis)
// that would collide with each other and with Matrix.ts's inverse and
// determinant under the library-wide flat export. Every fixed-size free
// function is therefore suffixed with its matrix size:
//   inverse2x2, adjoint2x2, determinant2x2, trace2x2, doTransform2x2,
//   setBasis2x2, getBasis2x2, makeRotation2x2, getRotationAngle2x2
// and likewise ...3x3 and ...4x4. The suffix is applied uniformly, including
// to names that happen not to collide today (makeReflection4x4, ...), so
// that later fixed-size headers can follow the same rule mechanically.
//
// Multiplication convention: upstream selects matrix-on-the-left
// (GTE_USE_MAT_VEC, the GTE default) or vector-on-the-left
// (GTE_USE_VEC_MAT) at compile time. The port always uses GTE_USE_MAT_VEC,
// so doTransform2x2(M, V) is M*V, doTransform2x2(A, B) is A*B, the basis
// vectors are the columns, and the rotation matrix is R(t) = {{c,-s},{s,c}}.
//
// Upstream's optional 'bool* reportInvertibility' output parameter becomes a
// field of the object returned by inverse2x2, matching Matrix.ts's inverse.

import { logAssert } from './Logger.js';
import { Matrix, multiplyAB, mulMatrix } from './Matrix.js';
import { Vector } from './Vector.js';

// The port of upstream's template alias. The port's Matrix carries its
// dimensions at run time, so this is documentation rather than a distinct
// type; use Matrix.fromArray(2, 2, [...]) to build one from row-major values.
export type Matrix2x2 = Matrix;

function assert2x2(M: Matrix): void {
    logAssert(M.numRows === 2 && M.numCols === 2, 'Matrix must be 2-by-2.');
}

// Create a rotation matrix from an angle (in radians). The matrix is
//   R(t) = {{c,-s},{s,c}}
// where c = cos(t), s = sin(t), and the inner-brace pairs are rows of the
// matrix. The rotation is written into the 2-by-2 matrix 'rotation', as
// upstream's output reference parameter.
export function makeRotation2x2(angle: number, rotation: Matrix): void {
    assert2x2(rotation);
    const cs = Math.cos(angle);
    const sn = Math.sin(angle);
    rotation.set(0, 0, cs);
    rotation.set(0, 1, -sn);
    rotation.set(1, 0, sn);
    rotation.set(1, 1, cs);
}

// Get the angle (radians) from a rotation matrix. The caller is responsible
// for ensuring the matrix is a rotation.
export function getRotationAngle2x2(rotation: Matrix): number {
    assert2x2(rotation);
    return Math.atan2(rotation.get(1, 0), rotation.get(0, 0));
}

// Geometric operations. When M is not invertible, the returned matrix is
// zero and 'invertible' is false.
export function inverse2x2(M: Matrix): { inverse: Matrix, invertible: boolean } {
    assert2x2(M);
    const det = M.get(0, 0) * M.get(1, 1) - M.get(0, 1) * M.get(1, 0);
    if (det !== 0) {
        const invDet = 1 / det;
        const inv = Matrix.fromArray(2, 2, [
            M.get(1, 1) * invDet, -M.get(0, 1) * invDet,
            -M.get(1, 0) * invDet, M.get(0, 0) * invDet
        ]);
        return { inverse: inv, invertible: true };
    }
    return { inverse: new Matrix(2, 2), invertible: false };
}

export function adjoint2x2(M: Matrix): Matrix {
    assert2x2(M);
    return Matrix.fromArray(2, 2, [
        M.get(1, 1), -M.get(0, 1),
        -M.get(1, 0), M.get(0, 0)
    ]);
}

export function determinant2x2(M: Matrix): number {
    assert2x2(M);
    return M.get(0, 0) * M.get(1, 1) - M.get(0, 1) * M.get(1, 0);
}

export function trace2x2(M: Matrix): number {
    assert2x2(M);
    return M.get(0, 0) + M.get(1, 1);
}

// Multiply M and V (or A and B) according to the multiplication convention.
// The port uses GTE_USE_MAT_VEC, so the results are M*V and A*B.
export function doTransform2x2(M: Matrix, V: Vector): Vector;
export function doTransform2x2(A: Matrix, B: Matrix): Matrix;
export function doTransform2x2(M: Matrix, arg: Vector | Matrix): Vector | Matrix {
    assert2x2(M);
    if (arg instanceof Matrix) {
        assert2x2(arg);
        return multiplyAB(M, arg);
    }
    logAssert(arg.size === 2, 'Vector must have 2 components.');
    return mulMatrix(M, arg);
}

// The columns of an invertible matrix form a basis for the range of the
// matrix (GTE_USE_MAT_VEC). These functions allow you to access the basis
// vectors. The caller is responsible for ensuring that the matrix is
// invertible (although the inverse is not calculated by these functions).
export function setBasis2x2(M: Matrix, i: number, V: Vector): void {
    assert2x2(M);
    logAssert(V.size === 2, 'Vector must have 2 components.');
    M.setCol(i, V);
}

export function getBasis2x2(M: Matrix, i: number): Vector {
    assert2x2(M);
    return M.getCol(i);
}
