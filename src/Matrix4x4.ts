// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Matrix4x4.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: see Matrix2x2.ts for the fixed-size naming scheme (every free
// function of the fixed-size matrix headers is suffixed with its matrix
// size) and for the multiplication convention (the port always uses
// GTE_USE_MAT_VEC, so doTransform4x4 is M*V and A*B, the basis vectors are
// the columns, and the special matrices below are the ones shown in the
// upstream comments).

import { logAssert } from './Logger';
import { Matrix, multiplyAB, mulMatrix } from './Matrix';
import { Vector, dot, sub } from './Vector';

// The port of upstream's template alias 'using Matrix4x4 = Matrix<4,4,Real>'.
export type Matrix4x4 = Matrix;

function assert4x4(M: Matrix): void {
    logAssert(M.numRows === 4 && M.numCols === 4, 'Matrix must be 4-by-4.');
}

function assertVector4(V: Vector): void {
    logAssert(V.size === 4, 'Vector must have 4 components.');
}

// The 2-by-2 sub-determinants of the upper two rows (a0..a5) and of the
// lower two rows (b0..b5), shared by inverse4x4, adjoint4x4 and
// determinant4x4.
function subDeterminants(M: Matrix): { a: number[], b: number[] } {
    const a = [
        M.get(0, 0) * M.get(1, 1) - M.get(0, 1) * M.get(1, 0),
        M.get(0, 0) * M.get(1, 2) - M.get(0, 2) * M.get(1, 0),
        M.get(0, 0) * M.get(1, 3) - M.get(0, 3) * M.get(1, 0),
        M.get(0, 1) * M.get(1, 2) - M.get(0, 2) * M.get(1, 1),
        M.get(0, 1) * M.get(1, 3) - M.get(0, 3) * M.get(1, 1),
        M.get(0, 2) * M.get(1, 3) - M.get(0, 3) * M.get(1, 2)
    ];
    const b = [
        M.get(2, 0) * M.get(3, 1) - M.get(2, 1) * M.get(3, 0),
        M.get(2, 0) * M.get(3, 2) - M.get(2, 2) * M.get(3, 0),
        M.get(2, 0) * M.get(3, 3) - M.get(2, 3) * M.get(3, 0),
        M.get(2, 1) * M.get(3, 2) - M.get(2, 2) * M.get(3, 1),
        M.get(2, 1) * M.get(3, 3) - M.get(2, 3) * M.get(3, 1),
        M.get(2, 2) * M.get(3, 3) - M.get(2, 3) * M.get(3, 2)
    ];
    return { a, b };
}

// The row-major elements of the adjugate (classical adjoint) of M.
function adjointElements(M: Matrix, a: number[], b: number[]): number[] {
    return [
        +M.get(1, 1) * b[5] - M.get(1, 2) * b[4] + M.get(1, 3) * b[3],
        -M.get(0, 1) * b[5] + M.get(0, 2) * b[4] - M.get(0, 3) * b[3],
        +M.get(3, 1) * a[5] - M.get(3, 2) * a[4] + M.get(3, 3) * a[3],
        -M.get(2, 1) * a[5] + M.get(2, 2) * a[4] - M.get(2, 3) * a[3],
        -M.get(1, 0) * b[5] + M.get(1, 2) * b[2] - M.get(1, 3) * b[1],
        +M.get(0, 0) * b[5] - M.get(0, 2) * b[2] + M.get(0, 3) * b[1],
        -M.get(3, 0) * a[5] + M.get(3, 2) * a[2] - M.get(3, 3) * a[1],
        +M.get(2, 0) * a[5] - M.get(2, 2) * a[2] + M.get(2, 3) * a[1],
        +M.get(1, 0) * b[4] - M.get(1, 1) * b[2] + M.get(1, 3) * b[0],
        -M.get(0, 0) * b[4] + M.get(0, 1) * b[2] - M.get(0, 3) * b[0],
        +M.get(3, 0) * a[4] - M.get(3, 1) * a[2] + M.get(3, 3) * a[0],
        -M.get(2, 0) * a[4] + M.get(2, 1) * a[2] - M.get(2, 3) * a[0],
        -M.get(1, 0) * b[3] + M.get(1, 1) * b[1] - M.get(1, 2) * b[0],
        +M.get(0, 0) * b[3] - M.get(0, 1) * b[1] + M.get(0, 2) * b[0],
        -M.get(3, 0) * a[3] + M.get(3, 1) * a[1] - M.get(3, 2) * a[0],
        +M.get(2, 0) * a[3] - M.get(2, 1) * a[1] + M.get(2, 2) * a[0]
    ];
}

// Geometric operations. When M is not invertible, the returned matrix is
// zero and 'invertible' is false.
export function inverse4x4(M: Matrix): { inverse: Matrix, invertible: boolean } {
    assert4x4(M);
    const { a, b } = subDeterminants(M);
    const det = a[0] * b[5] - a[1] * b[4] + a[2] * b[3] + a[3] * b[2]
        - a[4] * b[1] + a[5] * b[0];
    if (det !== 0) {
        const invDet = 1 / det;
        const elements = adjointElements(M, a, b);
        for (let i = 0; i < 16; ++i) {
            elements[i] *= invDet;
        }
        return { inverse: Matrix.fromArray(4, 4, elements), invertible: true };
    }
    return { inverse: new Matrix(4, 4), invertible: false };
}

export function adjoint4x4(M: Matrix): Matrix {
    assert4x4(M);
    const { a, b } = subDeterminants(M);
    return Matrix.fromArray(4, 4, adjointElements(M, a, b));
}

export function determinant4x4(M: Matrix): number {
    assert4x4(M);
    const { a, b } = subDeterminants(M);
    return a[0] * b[5] - a[1] * b[4] + a[2] * b[3] + a[3] * b[2]
        - a[4] * b[1] + a[5] * b[0];
}

export function trace4x4(M: Matrix): number {
    assert4x4(M);
    return M.get(0, 0) + M.get(1, 1) + M.get(2, 2) + M.get(3, 3);
}

// Multiply M and V (or A and B) according to the multiplication convention.
// The port uses GTE_USE_MAT_VEC, so the results are M*V and A*B.
export function doTransform4x4(M: Matrix, V: Vector): Vector;
export function doTransform4x4(A: Matrix, B: Matrix): Matrix;
export function doTransform4x4(M: Matrix, arg: Vector | Matrix): Vector | Matrix {
    assert4x4(M);
    if (arg instanceof Matrix) {
        assert4x4(arg);
        return multiplyAB(M, arg);
    }
    assertVector4(arg);
    return mulMatrix(M, arg);
}

// The columns of an invertible matrix form a basis for the range of the
// matrix (GTE_USE_MAT_VEC). The caller is responsible for ensuring that the
// matrix is invertible (although the inverse is not calculated here).
export function setBasis4x4(M: Matrix, i: number, V: Vector): void {
    assert4x4(M);
    assertVector4(V);
    M.setCol(i, V);
}

export function getBasis4x4(M: Matrix, i: number): Vector {
    assert4x4(M);
    return M.getCol(i);
}

// Special matrices, shown here using the GTE_USE_MAT_VEC multiplication
// convention.

// The projection plane is Dot(N,X-P) = 0 where N is a 3-by-1 unit-length
// normal vector and P is a 3-by-1 point on the plane. The projection is
// oblique to the plane, in the direction of the 3-by-1 vector D.
// Necessarily Dot(N,D) is not zero for this projection to make sense. Given
// a 3-by-1 point U, compute the intersection of the line U+t*D with the
// plane to obtain t = -Dot(N,U-P)/Dot(N,D); then
//
//   projection(U) = P + [I - D*N^T/Dot(N,D)]*(U-P)
//
// A 4-by-4 homogeneous transformation representing the projection is
//
//       +-                               -+
//   M = | D*N^T - Dot(N,D)*I   -Dot(N,P)D |
//       |          0^T          -Dot(N,D) |
//       +-                               -+
//
// where M applies to [U^T 1]^T by M*[U^T 1]^T. The matrix is chosen so that
// M(3,3) > 0 whenever Dot(N,D) < 0; the projection is onto the "positive
// side" of the plane.
export function makeObliqueProjection4x4(origin: Vector, normal: Vector,
    direction: Vector): Matrix {
    assertVector4(origin);
    assertVector4(normal);
    assertVector4(direction);

    const M = new Matrix(4, 4);
    const dotND = dot(normal, direction);
    const dotNO = dot(origin, normal);

    M.set(0, 0, direction.values[0] * normal.values[0] - dotND);
    M.set(0, 1, direction.values[0] * normal.values[1]);
    M.set(0, 2, direction.values[0] * normal.values[2]);
    M.set(0, 3, -dotNO * direction.values[0]);
    M.set(1, 0, direction.values[1] * normal.values[0]);
    M.set(1, 1, direction.values[1] * normal.values[1] - dotND);
    M.set(1, 2, direction.values[1] * normal.values[2]);
    M.set(1, 3, -dotNO * direction.values[1]);
    M.set(2, 0, direction.values[2] * normal.values[0]);
    M.set(2, 1, direction.values[2] * normal.values[1]);
    M.set(2, 2, direction.values[2] * normal.values[2] - dotND);
    M.set(2, 3, -dotNO * direction.values[2]);
    M.set(3, 0, 0);
    M.set(3, 1, 0);
    M.set(3, 2, 0);
    M.set(3, 3, -dotND);

    return M;
}

// The perspective projection of a point onto a plane is
//
//     +-                                                 -+
// M = | Dot(N,E-P)*I - E*N^T    -(Dot(N,E-P)*I - E*N^T)*E |
//     |        -N^T                      Dot(N,E)         |
//     +-                                                 -+
//
// where E is the eye point, P is a point on the plane, and N is a
// unit-length plane normal.
export function makePerspectiveProjection4x4(origin: Vector, normal: Vector,
    eye: Vector): Matrix {
    assertVector4(origin);
    assertVector4(normal);
    assertVector4(eye);

    const M = new Matrix(4, 4);
    const dotND = dot(normal, sub(eye, origin));

    M.set(0, 0, dotND - eye.values[0] * normal.values[0]);
    M.set(0, 1, -eye.values[0] * normal.values[1]);
    M.set(0, 2, -eye.values[0] * normal.values[2]);
    M.set(0, 3, -(M.get(0, 0) * eye.values[0] + M.get(0, 1) * eye.values[1]
        + M.get(0, 2) * eye.values[2]));
    M.set(1, 0, -eye.values[1] * normal.values[0]);
    M.set(1, 1, dotND - eye.values[1] * normal.values[1]);
    M.set(1, 2, -eye.values[1] * normal.values[2]);
    M.set(1, 3, -(M.get(1, 0) * eye.values[0] + M.get(1, 1) * eye.values[1]
        + M.get(1, 2) * eye.values[2]));
    M.set(2, 0, -eye.values[2] * normal.values[0]);
    M.set(2, 1, -eye.values[2] * normal.values[1]);
    M.set(2, 2, dotND - eye.values[2] * normal.values[2]);
    M.set(2, 3, -(M.get(2, 0) * eye.values[0] + M.get(2, 1) * eye.values[1]
        + M.get(2, 2) * eye.values[2]));
    M.set(3, 0, -normal.values[0]);
    M.set(3, 1, -normal.values[1]);
    M.set(3, 2, -normal.values[2]);
    M.set(3, 3, dot(eye, normal));

    return M;
}

// The reflection of a point through a plane is
//     +-                         -+
// M = | I-2*N*N^T    2*Dot(N,P)*N |
//     |     0^T            1      |
//     +-                         -+
//
// where P is a point on the plane and N is a unit-length plane normal.
export function makeReflection4x4(origin: Vector, normal: Vector): Matrix {
    assertVector4(origin);
    assertVector4(normal);

    const M = new Matrix(4, 4);
    const twoDotNO = 2 * dot(origin, normal);

    M.set(0, 0, 1 - 2 * normal.values[0] * normal.values[0]);
    M.set(0, 1, -2 * normal.values[0] * normal.values[1]);
    M.set(0, 2, -2 * normal.values[0] * normal.values[2]);
    M.set(0, 3, twoDotNO * normal.values[0]);
    M.set(1, 0, M.get(0, 1));
    M.set(1, 1, 1 - 2 * normal.values[1] * normal.values[1]);
    M.set(1, 2, -2 * normal.values[1] * normal.values[2]);
    M.set(1, 3, twoDotNO * normal.values[1]);
    M.set(2, 0, M.get(0, 2));
    M.set(2, 1, M.get(1, 2));
    M.set(2, 2, 1 - 2 * normal.values[2] * normal.values[2]);
    M.set(2, 3, twoDotNO * normal.values[2]);
    M.set(3, 0, 0);
    M.set(3, 1, 0);
    M.set(3, 2, 0);
    M.set(3, 3, 1);

    return M;
}
