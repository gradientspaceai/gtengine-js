// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Matrix3x3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: see Matrix2x2.ts for the fixed-size naming scheme (every free
// function of the fixed-size matrix headers is suffixed with its matrix
// size) and for the multiplication convention (the port always uses
// GTE_USE_MAT_VEC, so doTransform3x3 is M*V and A*B, and the basis vectors
// are the columns).

import { logAssert } from './Logger';
import { Matrix, multiplyAB, mulMatrix } from './Matrix';
import { Vector } from './Vector';

// The port of upstream's template alias 'using Matrix3x3 = Matrix<3,3,Real>'.
export type Matrix3x3 = Matrix;

function assert3x3(M: Matrix): void {
    logAssert(M.numRows === 3 && M.numCols === 3, 'Matrix must be 3-by-3.');
}

// Geometric operations. The inverse is the adjoint divided by the
// determinant, with the determinant expanded along the first row using the
// cofactors c00, c10 and c20. When M is not invertible, the returned matrix
// is zero and 'invertible' is false.
export function inverse3x3(M: Matrix): { inverse: Matrix, invertible: boolean } {
    assert3x3(M);
    const c00 = M.get(1, 1) * M.get(2, 2) - M.get(1, 2) * M.get(2, 1);
    const c10 = M.get(1, 2) * M.get(2, 0) - M.get(1, 0) * M.get(2, 2);
    const c20 = M.get(1, 0) * M.get(2, 1) - M.get(1, 1) * M.get(2, 0);
    const det = M.get(0, 0) * c00 + M.get(0, 1) * c10 + M.get(0, 2) * c20;
    if (det !== 0) {
        const invDet = 1 / det;
        const inv = Matrix.fromArray(3, 3, [
            c00 * invDet,
            (M.get(0, 2) * M.get(2, 1) - M.get(0, 1) * M.get(2, 2)) * invDet,
            (M.get(0, 1) * M.get(1, 2) - M.get(0, 2) * M.get(1, 1)) * invDet,
            c10 * invDet,
            (M.get(0, 0) * M.get(2, 2) - M.get(0, 2) * M.get(2, 0)) * invDet,
            (M.get(0, 2) * M.get(1, 0) - M.get(0, 0) * M.get(1, 2)) * invDet,
            c20 * invDet,
            (M.get(0, 1) * M.get(2, 0) - M.get(0, 0) * M.get(2, 1)) * invDet,
            (M.get(0, 0) * M.get(1, 1) - M.get(0, 1) * M.get(1, 0)) * invDet
        ]);
        return { inverse: inv, invertible: true };
    }
    return { inverse: new Matrix(3, 3), invertible: false };
}

export function adjoint3x3(M: Matrix): Matrix {
    assert3x3(M);
    return Matrix.fromArray(3, 3, [
        M.get(1, 1) * M.get(2, 2) - M.get(1, 2) * M.get(2, 1),
        M.get(0, 2) * M.get(2, 1) - M.get(0, 1) * M.get(2, 2),
        M.get(0, 1) * M.get(1, 2) - M.get(0, 2) * M.get(1, 1),
        M.get(1, 2) * M.get(2, 0) - M.get(1, 0) * M.get(2, 2),
        M.get(0, 0) * M.get(2, 2) - M.get(0, 2) * M.get(2, 0),
        M.get(0, 2) * M.get(1, 0) - M.get(0, 0) * M.get(1, 2),
        M.get(1, 0) * M.get(2, 1) - M.get(1, 1) * M.get(2, 0),
        M.get(0, 1) * M.get(2, 0) - M.get(0, 0) * M.get(2, 1),
        M.get(0, 0) * M.get(1, 1) - M.get(0, 1) * M.get(1, 0)
    ]);
}

export function determinant3x3(M: Matrix): number {
    assert3x3(M);
    const c00 = M.get(1, 1) * M.get(2, 2) - M.get(1, 2) * M.get(2, 1);
    const c10 = M.get(1, 2) * M.get(2, 0) - M.get(1, 0) * M.get(2, 2);
    const c20 = M.get(1, 0) * M.get(2, 1) - M.get(1, 1) * M.get(2, 0);
    return M.get(0, 0) * c00 + M.get(0, 1) * c10 + M.get(0, 2) * c20;
}

export function trace3x3(M: Matrix): number {
    assert3x3(M);
    return M.get(0, 0) + M.get(1, 1) + M.get(2, 2);
}

// Multiply M and V (or A and B) according to the multiplication convention.
// The port uses GTE_USE_MAT_VEC, so the results are M*V and A*B.
export function doTransform3x3(M: Matrix, V: Vector): Vector;
export function doTransform3x3(A: Matrix, B: Matrix): Matrix;
export function doTransform3x3(M: Matrix, arg: Vector | Matrix): Vector | Matrix {
    assert3x3(M);
    if (arg instanceof Matrix) {
        assert3x3(arg);
        return multiplyAB(M, arg);
    }
    logAssert(arg.size === 3, 'Vector must have 3 components.');
    return mulMatrix(M, arg);
}

// The columns of an invertible matrix form a basis for the range of the
// matrix (GTE_USE_MAT_VEC). The caller is responsible for ensuring that the
// matrix is invertible (although the inverse is not calculated here).
export function setBasis3x3(M: Matrix, i: number, V: Vector): void {
    assert3x3(M);
    logAssert(V.size === 3, 'Vector must have 3 components.');
    M.setCol(i, V);
}

export function getBasis3x3(M: Matrix, i: number): Vector {
    assert3x3(M);
    return M.getCol(i);
}
