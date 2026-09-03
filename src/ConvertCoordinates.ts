// gtengine-js: TypeScript port of Geometric Tools Engine (GTE)
// ConvertCoordinates.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Convert points and transformations between two coordinate systems. The
// mathematics involves a change of basis. See
//   https://www.geometrictools.com/Documentation/ConvertingBetweenCoordinateSystems.pdf
// for the details. Typical usage for 3D conversion is shown next.
//
// // Linear change of basis. The columns of U are the basis vectors for the
// // source coordinate system. A vector X = { x0, x1, x2 } in the source
// // coordinate system is represented by
// //   X = x0*(1,0,0) + x1*(0,1,0) + x2*(0,0,1)
// // The Cartesian coordinates for the point are the combination of these
// // terms,
// //   X = (x0, x1, x2)
// // The columns of V are the basis vectors for the target coordinate
// // system. A vector Y = { y0, y1, y2 } in the target coordinate system is
// // represented by
// //   Y = y0*(1,0,0) + y1*(0,0,1) + y2*(0,1,0)
// // The Cartesian coordinates for the vector are the combination of these
// // terms,
// //   Y = (y0, y2, y1)
// // The call Y = convert.uToV(X) computes y0, y1 and y2 so that the
// // Cartesian coordinates for X and for Y are the same. For example,
// //   X = { 1.0, 2.0, 3.0 }
// //     = 1.0*(1,0,0) + 2.0*(0,1,0) + 3.0*(0,0,1)
// //     = (1, 2, 3)
// //   Y = { 1.0, 3.0, 2.0 }
// //     = 1.0*(1,0,0) + 3.0*(0,0,1) + 2.0*(0,1,0)
// //     = (1, 2, 3)
// // X and Y represent the same vector (equal Cartesian coordinates) but
// // have different representations in the source and target coordinates.
//
// const convert = new ConvertCoordinates(3);
// const U = new Matrix(3, 3), V = new Matrix(3, 3);
// setBasis3x3(U, 0, Vector.fromArray([1, 0, 0]));
// setBasis3x3(U, 1, Vector.fromArray([0, 1, 0]));
// setBasis3x3(U, 2, Vector.fromArray([0, 0, 1]));
// setBasis3x3(V, 0, Vector.fromArray([1, 0, 0]));
// setBasis3x3(V, 1, Vector.fromArray([0, 0, 1]));
// setBasis3x3(V, 2, Vector.fromArray([0, 1, 0]));
// convert.compute(U, true, V, true);
// convert.isRightHandedU();  // true
// convert.isRightHandedV();  // false
// let X = Vector.fromArray([1, 2, 3]);
// let Y = convert.uToV(X);   // { 1, 3, 2 }
// // mulMatrix(U, X) equals mulMatrix(V, Y)
// Y = Vector.fromArray([0, 1, 2]);
// X = convert.vToU(Y);       // { 0, 2, 1 }
// // mulMatrix(U, X) equals mulMatrix(V, Y)
// const c = 0.6, s = 0.8;    // c*c + s*s = 1
// const A = new Matrix(3, 3);
// setBasis3x3(A, 0, Vector.fromArray([c, s, 0]));
// setBasis3x3(A, 1, Vector.fromArray([-s, c, 0]));
// setBasis3x3(A, 2, Vector.fromArray([0, 0, 1]));
// const B = convert.uToV(A);
//   // getBasis3x3(B, 0) = { c, 0, s}
//   // getBasis3x3(B, 1) = { 0, 1, 0}
//   // getBasis3x3(B, 2) = {-s, 0, c}
// X = mulMatrix(A, X);  // U is vector-on-right
// Y = mulMatrix(B, Y);  // V is vector-on-right
// // mulMatrix(U, X) equals mulMatrix(V, Y)
//
// // Affine change of basis. The first three columns of U are the basis
// // vectors for the source coordinate system and must have last components
// // set to 0. The last column is the origin for that system and must have
// // last component set to 1. A point X = { x0, x1, x2, 1 } in the source
// // coordinate system is represented by
// //   X = x0*(-1,0,0,0) + x1*(0,0,1,0) + x2*(0,-1,0,0) + 1*(1,2,3,1)
// // The Cartesian coordinates for the point are the combination of these
// // terms,
// //   X = (-x0 + 1, -x2 + 2, x1 + 3, 1)
// // The first three columns of V are the basis vectors for the target
// // coordinate system and must have last components set to 0. The last
// // column is the origin for that system and must have last component set
// // to 1. A point Y = { y0, y1, y2, 1 } in the target coordinate system is
// // represented by
// //   Y = y0*(0,1,0,0) + y1*(-1,0,0,0) + y2*(0,0,1,0) + 1*(4,5,6,1)
// // The Cartesian coordinates for the point are the combination of these
// // terms,
// //   Y = (-y1 + 4, y0 + 5, y2 + 6, 1)
// // The call Y = convert.uToV(X) computes y0, y1 and y2 so that the
// // Cartesian coordinates for X and for Y are the same. For example,
// //   X = { -1.0, 4.0, -3.0, 1.0 }
// //     = -1.0*(-1,0,0,0) + 4.0*(0,0,1,0) - 3.0*(0,-1,0,0) + 1.0*(1,2,3,1)
// //     = (2, 5, 7, 1)
// //   Y = { 0.0, 2.0, 1.0, 1.0 }
// //     = 0.0*(0,1,0,0) + 2.0*(-1,0,0,0) + 1.0*(0,0,1,0) + 1.0*(4,5,6,1)
// //     = (2, 5, 7, 1)
// // X and Y represent the same point (equal Cartesian coordinates) but have
// // different representations in the source and target affine coordinates.
//
// The corresponding 4D example, with convert.compute(U, true, V, false),
// gives isRightHandedU() == false, isRightHandedV() == true, and the
// transformation B = convert.uToV(A) is applied as Y = mulMatrix(Y, B)
// because V is vector-on-left. Both examples are exercised in the tests.
//
// Upstream documentation bug (preserved behavior, corrected comment): the
// header's affine example lists the four tuples { 1, 0, 0, 0 },
// { 0, c, s, 0 }, { 0, -s, c, 0 } and { 2.0, -0.9, -2.6, 1 } as
// "B.GetCol(0..3)". Because vectorOnRightV is false in that example, uToV
// transposes the product, so those tuples are the ROWS of B, not its
// columns. The code is correct; only the comment is wrong.
//
// Port notes:
// - Upstream's compile-time dimension 'template <int32_t N, typename Real>'
//   becomes a constructor argument, matching the port's runtime-dimensioned
//   Matrix and Vector.
// - operator() becomes compute(), returning upstream's bool.
// - GetC/GetInverseC return the live internal matrices (upstream returns
//   const references); do not mutate them.

import { GaussianElimination } from './GaussianElimination.js';
import { logAssert } from './Logger.js';
import { Matrix, multiplyAB, multiplyATB, mulMatrix, transpose } from './Matrix.js';
import { Vector } from './Vector.js';

export class ConvertCoordinates {
    private mN: number;
    // C = U^{-1}*V, C^{-1} = V^{-1}*U
    private mC: Matrix;
    private mInverseC: Matrix;
    private mIsVectorOnRightU: boolean;
    private mIsVectorOnRightV: boolean;
    private mIsRightHandedU: boolean;
    private mIsRightHandedV: boolean;

    // Construction of the change of basis matrix. The implementation
    // supports both linear change of basis (n = 3 for 3D) and affine change
    // of basis (n = 4 for 3D).
    constructor(n: number) {
        logAssert(n > 0, 'Invalid dimension.');
        this.mN = n;
        this.mC = Matrix.identity(n, n);
        this.mInverseC = Matrix.identity(n, n);
        this.mIsVectorOnRightU = true;
        this.mIsVectorOnRightV = true;
        this.mIsRightHandedU = true;
        this.mIsRightHandedV = true;
    }

    // Compute a change of basis between two coordinate systems. The return
    // value is true iff U and V are invertible. The matrix-vector
    // multiplication conventions affect the conversion of matrix
    // transformations. The Boolean inputs indicate how you want the matrices
    // to be interpreted when applied as transformations of a vector.
    compute(U: Matrix, vectorOnRightU: boolean, V: Matrix,
        vectorOnRightV: boolean): boolean {
        const n = this.mN;
        logAssert(U.numRows === n && U.numCols === n
            && V.numRows === n && V.numCols === n, 'Mismatched sizes.');

        // Initialize in case of early exit.
        this.mC = Matrix.identity(n, n);
        this.mInverseC = Matrix.identity(n, n);
        this.mIsVectorOnRightU = true;
        this.mIsVectorOnRightV = true;
        this.mIsRightHandedU = true;
        this.mIsRightHandedV = true;

        const resultU = new GaussianElimination().compute(n, U.values,
            { wantInverse: true });
        if (!resultU.invertible) {
            return false;
        }

        const resultV = new GaussianElimination().compute(n, V.values,
            { wantInverse: true });
        if (!resultV.invertible) {
            return false;
        }

        const inverseU = Matrix.fromArray(n, n, resultU.inverseM as number[]);
        const inverseV = Matrix.fromArray(n, n, resultV.inverseM as number[]);

        this.mC = multiplyAB(inverseU, V);
        this.mInverseC = multiplyAB(inverseV, U);
        this.mIsVectorOnRightU = vectorOnRightU;
        this.mIsVectorOnRightV = vectorOnRightV;
        this.mIsRightHandedU = (resultU.determinant > 0);
        this.mIsRightHandedV = (resultV.determinant > 0);
        return true;
    }

    // Member access.
    getDimension(): number {
        return this.mN;
    }

    getC(): Matrix {
        return this.mC;
    }

    getInverseC(): Matrix {
        return this.mInverseC;
    }

    isVectorOnRightU(): boolean {
        return this.mIsVectorOnRightU;
    }

    isVectorOnRightV(): boolean {
        return this.mIsVectorOnRightV;
    }

    isRightHandedU(): boolean {
        return this.mIsRightHandedU;
    }

    isRightHandedV(): boolean {
        return this.mIsRightHandedV;
    }

    // Convert points and transformations between coordinate systems. The
    // names of the systems are U and V to make it clear which inputs of
    // compute() they are associated with. The X vector stores coordinates
    // for the U-system and the Y vector stores coordinates for the V-system.
    //
    // For a vector, Y = C^{-1}*X.
    //
    // For a transformation, the output is computed from the tables below,
    // where the superscript T denotes the transpose operator.
    // vectorOnRightU = true:  transformation is X' = A*X
    // vectorOnRightU = false: transformation is (X')^T = X^T*A
    // vectorOnRightV = true:  transformation is Y' = B*Y
    // vectorOnRightV = false: transformation is (Y')^T = Y^T*B
    //
    // vectorOnRightU  | vectorOnRightV  | output
    // ----------------+-----------------+---------------------
    // true            | true            | C^{-1} * A * C
    // true            | false           | (C^{-1} * A * C)^T
    // false           | true            | C^{-1} * A^T * C
    // false           | false           | (C^{-1} * A^T * C)^T
    uToV(X: Vector): Vector;
    uToV(A: Matrix): Matrix;
    uToV(arg: Vector | Matrix): Vector | Matrix {
        if (arg instanceof Matrix) {
            const A = arg;
            logAssert(A.numRows === this.mN && A.numCols === this.mN,
                'Mismatched sizes.');
            const product = (this.mIsVectorOnRightU
                ? multiplyAB(multiplyAB(this.mInverseC, A), this.mC)
                : multiplyAB(this.mInverseC, multiplyATB(A, this.mC)));
            return (this.mIsVectorOnRightV ? product : transpose(product));
        }

        logAssert(arg.size === this.mN, 'Mismatched sizes.');
        return mulMatrix(this.mInverseC, arg);
    }

    // For a vector, X = C*Y.
    //
    // For a transformation:
    // vectorOnRightU  | vectorOnRightV  | output
    // ----------------+-----------------+---------------------
    // true            | true            | C * B * C^{-1}
    // true            | false           | C * B^T * C^{-1}
    // false           | true            | (C * B * C^{-1})^T
    // false           | false           | (C * B^T * C^{-1})^T
    vToU(Y: Vector): Vector;
    vToU(B: Matrix): Matrix;
    vToU(arg: Vector | Matrix): Vector | Matrix {
        if (arg instanceof Matrix) {
            const B = arg;
            logAssert(B.numRows === this.mN && B.numCols === this.mN,
                'Mismatched sizes.');
            const product = (this.mIsVectorOnRightV
                ? multiplyAB(multiplyAB(this.mC, B), this.mInverseC)
                : multiplyAB(this.mC, multiplyATB(B, this.mInverseC)));
            return (this.mIsVectorOnRightU ? product : transpose(product));
        }

        logAssert(arg.size === this.mN, 'Mismatched sizes.');
        return mulMatrix(this.mC, arg);
    }
}
