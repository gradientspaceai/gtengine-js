// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Matrix.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: upstream 'template <int32_t NumRows, int32_t NumCols, typename
// Real> class Matrix' becomes a class with runtime dimensions backed by a
// flat number[], exactly as Vector.h's compile-time N became a runtime
// dimension. GMatrix.ts subclasses this class, adding only the resizing and
// the bounds-checked/std::vector comparison behavior of upstream GMatrix.
//
// Canonical translations:
// - 'Matrix<R, C, Real> M' -> 'new Matrix(r, c)' (zero-filled, as upstream's
//   default constructor) or 'Matrix.fromArray(r, c, [...])' (row-major
//   values, the port of the std::array/initializer-list constructors).
// - 'M(r, c)' -> 'M.get(r, c)' / 'M.set(r, c, x)'.
// - 'M[i]' (storage-order-independent 1D access) -> 'M.getFlat(i)' /
//   'M.setFlat(i, x)'; the backing array is exposed as 'M.values' for hot
//   inner loops (do not change its length except through GMatrix.setSize).
// - Comparison operators -> equals, notEquals, lessThan, lessThanOrEqual,
//   greaterThan, greaterThanOrEqual (lexicographic over the row-major
//   elements, matching std::array's operators).
// - Value semantics: C++ assignment copies; use 'M.clone()' where upstream
//   copies a matrix.
//
// Storage: upstream selects row-major or column-major at compile time
// (GTE_USE_ROW_MAJOR / GTE_USE_COL_MAJOR). The port always uses row major,
// the GTE default, which is also the default of the GaussianElimination
// port used by inverse() and determinant().
//
// Free-function naming: the operators would collide with Vector.h's exported
// negate/add/sub/mul/div and hlift/hproject under the library-wide flat
// export, so the matrix versions are suffixed with the type context:
// negateMatrix, addMatrix, subMatrix, mulMatrix, divMatrix, hliftMatrix and
// hprojectMatrix. The remaining upstream free functions keep their upstream
// names in camelCase (l1Norm, l2Norm, lInfinityNorm, inverse, determinant,
// transpose, multiplyAB, multiplyABT, multiplyATB, multiplyATBT, multiplyMD,
// multiplyDM, outerProduct, makeDiagonal). Compound assignments (+=, -=, ...)
// have no in-place ports; write 'M0 = addMatrix(M0, M1)'. Unary 'operator+'
// is the identity and has no port.
//
// Dimension mismatches (compile errors upstream for Matrix, LogError calls
// for GMatrix) throw via logAssert/logError in both cases.

import { GaussianElimination } from './GaussianElimination';
import { logAssert, logError } from './Logger';
import { Vector } from './Vector';

export class Matrix {
    // The matrix elements in row-major order. The contents may be read and
    // written directly, but the length must never be changed (except by
    // GMatrix.setSize).
    readonly values: number[];

    protected mNumRows: number;
    protected mNumCols: number;

    // Create a numRows-by-numCols matrix with all elements 0. (Upstream's
    // default constructor calls MakeZero, so this matches.)
    constructor(numRows: number, numCols: number) {
        logAssert(numRows >= 0 && numCols >= 0, 'Invalid size.');
        this.mNumRows = numRows;
        this.mNumCols = numCols;
        this.values = new Array<number>(numRows * numCols).fill(0);
    }

    // Create a matrix whose elements are the input values, specified in
    // row-major order regardless of the storage scheme. At most
    // numRows*numCols values are copied; any remaining elements are zero
    // (the port of both the std::array constructor and the
    // std::initializer_list constructor).
    static fromArray(numRows: number, numCols: number,
        values: readonly number[]): Matrix {
        const result = new Matrix(numRows, numCols);
        const numValues = Math.min(values.length, numRows * numCols);
        for (let i = 0; i < numValues; ++i) {
            result.values[i] = values[i];
        }
        return result;
    }

    // Special matrices. All elements are 0.
    static zero(numRows: number, numCols: number): Matrix {
        return new Matrix(numRows, numCols);
    }

    // Element (r,c) is 1, all others are 0. If either of r or c is invalid,
    // the zero matrix is created. This is a convenience for creating the
    // standard Euclidean basis matrices.
    static unit(numRows: number, numCols: number, r: number,
        c: number): Matrix {
        const result = new Matrix(numRows, numCols);
        result.makeUnit(r, c);
        return result;
    }

    // Diagonal entries 1, others 0, even when nonsquare.
    static identity(numRows: number, numCols: number): Matrix {
        const result = new Matrix(numRows, numCols);
        result.makeIdentity();
        return result;
    }

    // Member access for which the storage representation is transparent.
    get numRows(): number {
        return this.mNumRows;
    }

    get numCols(): number {
        return this.mNumCols;
    }

    get numElements(): number {
        return this.values.length;
    }

    getNumRows(): number {
        return this.mNumRows;
    }

    getNumCols(): number {
        return this.mNumCols;
    }

    getNumElements(): number {
        return this.values.length;
    }

    // The port of GMatrix's GetSize(numRows, numCols) output parameters.
    getSize(): { numRows: number, numCols: number } {
        return { numRows: this.mNumRows, numCols: this.mNumCols };
    }

    // The matrix entry in row r and column c. As upstream's fixed-size
    // Matrix, the indices are not range checked (GMatrix overrides these to
    // throw on invalid indices).
    get(r: number, c: number): number {
        return this.values[c + this.mNumCols * r];
    }

    set(r: number, c: number, value: number): void {
        this.values[c + this.mNumCols * r] = value;
    }

    // Member access by 1-dimensional index. NOTE: These accessors are useful
    // for the manipulation of matrix entries when it does not matter whether
    // storage is row-major or column-major. Do not use constructs such as
    // M.values[c + numCols*r] that expose the storage convention.
    getFlat(i: number): number {
        return this.values[i];
    }

    setFlat(i: number, value: number): void {
        this.values[i] = value;
    }

    // Member access by rows or by columns. The input vectors must have the
    // correct number of elements for the matrix size.
    setRow(r: number, vec: Vector): void {
        for (let c = 0; c < this.mNumCols; ++c) {
            this.set(r, c, vec.values[c]);
        }
    }

    setCol(c: number, vec: Vector): void {
        for (let r = 0; r < this.mNumRows; ++r) {
            this.set(r, c, vec.values[r]);
        }
    }

    getRow(r: number): Vector {
        const vec = new Vector(this.mNumCols);
        for (let c = 0; c < this.mNumCols; ++c) {
            vec.values[c] = this.get(r, c);
        }
        return vec;
    }

    getCol(c: number): Vector {
        const vec = new Vector(this.mNumRows);
        for (let r = 0; r < this.mNumRows; ++r) {
            vec.values[r] = this.get(r, c);
        }
        return vec;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Matrix {
        return Matrix.fromArray(this.mNumRows, this.mNumCols, this.values);
    }

    // All elements are 0.
    makeZero(): void {
        this.values.fill(0);
    }

    // Element (r,c) is 1, all others zero. As upstream's fixed-size Matrix,
    // an invalid (r,c) produces the zero matrix (GMatrix overrides this to
    // throw).
    makeUnit(r: number, c: number): void {
        this.makeZero();
        if (0 <= r && r < this.mNumRows && 0 <= c && c < this.mNumCols) {
            this.set(r, c, 1);
        }
    }

    // Diagonal entries 1, others 0, even when nonsquare.
    makeIdentity(): void {
        this.makeZero();
        const numDiagonal =
            (this.mNumRows <= this.mNumCols ? this.mNumRows : this.mNumCols);
        for (let i = 0; i < numDiagonal; ++i) {
            this.set(i, i, 1);
        }
    }

    // Comparisons for sorted containers and geometric ordering. Upstream
    // compares the storage arrays lexicographically; the sizes are equal by
    // construction, so comparing matrices of different dimensions (a compile
    // error upstream) throws here. GMatrix overrides these with upstream
    // GMatrix's dimension-aware semantics.
    protected compare(mat: Matrix): number {
        logAssert(this.mNumRows === mat.mNumRows
            && this.mNumCols === mat.mNumCols, 'Matrix: mismatched sizes.');
        for (let i = 0; i < this.values.length; ++i) {
            if (this.values[i] < mat.values[i]) {
                return -1;
            }
            if (this.values[i] > mat.values[i]) {
                return +1;
            }
        }
        return 0;
    }

    equals(mat: Matrix): boolean {
        return this.compare(mat) === 0;
    }

    notEquals(mat: Matrix): boolean {
        return this.compare(mat) !== 0;
    }

    lessThan(mat: Matrix): boolean {
        return this.compare(mat) < 0;
    }

    lessThanOrEqual(mat: Matrix): boolean {
        return this.compare(mat) <= 0;
    }

    greaterThan(mat: Matrix): boolean {
        return this.compare(mat) > 0;
    }

    greaterThanOrEqual(mat: Matrix): boolean {
        return this.compare(mat) >= 0;
    }
}

function assertSameSize(M0: Matrix, M1: Matrix): void {
    logAssert(M0.numRows === M1.numRows && M0.numCols === M1.numCols,
        'Mismatched sizes');
}

// Unary operations. The port of unary 'operator-'.
export function negateMatrix(M: Matrix): Matrix {
    const result = new Matrix(M.numRows, M.numCols);
    for (let i = 0; i < M.numElements; ++i) {
        result.values[i] = -M.values[i];
    }
    return result;
}

// Linear-algebraic operations.
export function addMatrix(M0: Matrix, M1: Matrix): Matrix {
    assertSameSize(M0, M1);
    const result = new Matrix(M0.numRows, M0.numCols);
    for (let i = 0; i < M0.numElements; ++i) {
        result.values[i] = M0.values[i] + M1.values[i];
    }
    return result;
}

export function subMatrix(M0: Matrix, M1: Matrix): Matrix {
    assertSameSize(M0, M1);
    const result = new Matrix(M0.numRows, M0.numCols);
    for (let i = 0; i < M0.numElements; ++i) {
        result.values[i] = M0.values[i] - M1.values[i];
    }
    return result;
}

// The port of upstream's 'operator*' overloads: 'M * scalar', 'scalar * M',
// 'M * V' (the vector M*V), 'V * M' (the vector V^T*M) and 'A * B' (the
// matrix product, which forwards to multiplyAB).
export function mulMatrix(M: Matrix, scalar: number): Matrix;
export function mulMatrix(scalar: number, M: Matrix): Matrix;
export function mulMatrix(M: Matrix, V: Vector): Vector;
export function mulMatrix(V: Vector, M: Matrix): Vector;
export function mulMatrix(A: Matrix, B: Matrix): Matrix;
export function mulMatrix(arg0: Matrix | Vector | number,
    arg1: Matrix | Vector | number): Matrix | Vector {
    if (typeof arg0 === 'number' || typeof arg1 === 'number') {
        const M = (typeof arg0 === 'number' ? arg1 : arg0) as Matrix;
        const scalar = (typeof arg0 === 'number' ? arg0 : arg1) as number;
        const result = new Matrix(M.numRows, M.numCols);
        for (let i = 0; i < M.numElements; ++i) {
            result.values[i] = M.values[i] * scalar;
        }
        return result;
    }

    if (arg0 instanceof Matrix) {
        if (arg1 instanceof Matrix) {
            return multiplyAB(arg0, arg1);
        }

        // M*V
        const M = arg0;
        const V = arg1;
        logAssert(V.size === M.numCols, 'Mismatched sizes.');
        const result = new Vector(M.numRows);
        for (let r = 0; r < M.numRows; ++r) {
            result.values[r] = 0;
            for (let c = 0; c < M.numCols; ++c) {
                result.values[r] += M.get(r, c) * V.values[c];
            }
        }
        return result;
    }

    // V^T*M
    const V = arg0;
    const M = arg1 as Matrix;
    logAssert(V.size === M.numRows, 'Mismatched sizes.');
    const result = new Vector(M.numCols);
    for (let c = 0; c < M.numCols; ++c) {
        result.values[c] = 0;
        for (let r = 0; r < M.numRows; ++r) {
            result.values[c] += V.values[r] * M.get(r, c);
        }
    }
    return result;
}

// The port of 'M / scalar'. As upstream's fixed-size Matrix, division by
// zero produces the zero matrix, and the division is performed as
// multiplication by 1/scalar. (Upstream's GMatrix operator/= instead throws
// "Division by zero."; see the port notes in GMatrix.ts.)
export function divMatrix(M: Matrix, scalar: number): Matrix {
    const result = new Matrix(M.numRows, M.numCols);
    if (scalar !== 0) {
        const invScalar = 1 / scalar;
        for (let i = 0; i < M.numElements; ++i) {
            result.values[i] = M.values[i] * invScalar;
        }
    }
    return result;
}

// Geometric operations.
export function l1Norm(M: Matrix): number {
    let sum = 0;
    for (let i = 0; i < M.numElements; ++i) {
        sum += Math.abs(M.values[i]);
    }
    return sum;
}

export function l2Norm(M: Matrix): number {
    let sum = 0;
    for (let i = 0; i < M.numElements; ++i) {
        sum += M.values[i] * M.values[i];
    }
    return Math.sqrt(sum);
}

export function lInfinityNorm(M: Matrix): number {
    let maxAbsElement = 0;
    for (let i = 0; i < M.numElements; ++i) {
        const absElement = Math.abs(M.values[i]);
        if (absElement > maxAbsElement) {
            maxAbsElement = absElement;
        }
    }
    return maxAbsElement;
}

// The inverse of the square matrix M, computed by Gaussian elimination with
// full pivoting. Upstream's optional 'reportInvertibility' output parameter
// becomes a field of the returned object. When M is not invertible, the
// returned matrix is zero (the GaussianElimination behavior).
export function inverse(M: Matrix):
    { inverse: Matrix, invertible: boolean } {
    if (M.numRows === M.numCols) {
        const result = new GaussianElimination().compute(M.numRows, M.values,
            { wantInverse: true });
        const invM = new Matrix(M.numRows, M.numCols);
        const elements = result.inverseM as number[];
        for (let i = 0; i < invM.numElements; ++i) {
            invM.values[i] = elements[i];
        }
        return { inverse: invM, invertible: result.invertible };
    }
    return logError('Matrix must be square.');
}

export function determinant(M: Matrix): number {
    if (M.numRows === M.numCols) {
        return new GaussianElimination().compute(M.numRows, M.values).determinant;
    }
    return logError('Matrix must be square.');
}

// M^T
export function transpose(M: Matrix): Matrix {
    const result = new Matrix(M.numCols, M.numRows);
    for (let r = 0; r < M.numRows; ++r) {
        for (let c = 0; c < M.numCols; ++c) {
            result.set(c, r, M.get(r, c));
        }
    }
    return result;
}

// A*B
export function multiplyAB(A: Matrix, B: Matrix): Matrix {
    logAssert(A.numCols === B.numRows, 'Mismatched sizes.');
    const result = new Matrix(A.numRows, B.numCols);
    const numCommon = A.numCols;
    for (let r = 0; r < result.numRows; ++r) {
        for (let c = 0; c < result.numCols; ++c) {
            let value = 0;
            for (let i = 0; i < numCommon; ++i) {
                value += A.get(r, i) * B.get(i, c);
            }
            result.set(r, c, value);
        }
    }
    return result;
}

// A*B^T
export function multiplyABT(A: Matrix, B: Matrix): Matrix {
    logAssert(A.numCols === B.numCols, 'Mismatched sizes.');
    const result = new Matrix(A.numRows, B.numRows);
    const numCommon = A.numCols;
    for (let r = 0; r < result.numRows; ++r) {
        for (let c = 0; c < result.numCols; ++c) {
            let value = 0;
            for (let i = 0; i < numCommon; ++i) {
                value += A.get(r, i) * B.get(c, i);
            }
            result.set(r, c, value);
        }
    }
    return result;
}

// A^T*B
export function multiplyATB(A: Matrix, B: Matrix): Matrix {
    logAssert(A.numRows === B.numRows, 'Mismatched sizes.');
    const result = new Matrix(A.numCols, B.numCols);
    const numCommon = A.numRows;
    for (let r = 0; r < result.numRows; ++r) {
        for (let c = 0; c < result.numCols; ++c) {
            let value = 0;
            for (let i = 0; i < numCommon; ++i) {
                value += A.get(i, r) * B.get(i, c);
            }
            result.set(r, c, value);
        }
    }
    return result;
}

// A^T*B^T
export function multiplyATBT(A: Matrix, B: Matrix): Matrix {
    logAssert(A.numRows === B.numCols, 'Mismatched sizes.');
    const result = new Matrix(A.numCols, B.numRows);
    const numCommon = A.numRows;
    for (let r = 0; r < result.numRows; ++r) {
        for (let c = 0; c < result.numCols; ++c) {
            let value = 0;
            for (let i = 0; i < numCommon; ++i) {
                value += A.get(i, r) * B.get(c, i);
            }
            result.set(r, c, value);
        }
    }
    return result;
}

// M*D, where D is the diagonal matrix whose diagonal entries are the
// components of the vector D.
export function multiplyMD(M: Matrix, D: Vector): Matrix {
    logAssert(D.size === M.numCols, 'Mismatched sizes.');
    const result = new Matrix(M.numRows, M.numCols);
    for (let r = 0; r < result.numRows; ++r) {
        for (let c = 0; c < result.numCols; ++c) {
            result.set(r, c, M.get(r, c) * D.values[c]);
        }
    }
    return result;
}

// D*M, where D is the diagonal matrix whose diagonal entries are the
// components of the vector D.
export function multiplyDM(D: Vector, M: Matrix): Matrix {
    logAssert(D.size === M.numRows, 'Mismatched sizes.');
    const result = new Matrix(M.numRows, M.numCols);
    for (let r = 0; r < result.numRows; ++r) {
        for (let c = 0; c < result.numCols; ++c) {
            result.set(r, c, D.values[r] * M.get(r, c));
        }
    }
    return result;
}

// U*V^T, U is N-by-1, V is M-by-1, result is N-by-M.
export function outerProduct(U: Vector, V: Vector): Matrix {
    const result = new Matrix(U.size, V.size);
    for (let r = 0; r < result.numRows; ++r) {
        for (let c = 0; c < result.numCols; ++c) {
            result.set(r, c, U.values[r] * V.values[c]);
        }
    }
    return result;
}

// Initialization of M to a diagonal matrix whose diagonal entries are the
// components of D. The port uses upstream GMatrix.h's version, which also
// handles nonsquare M (the number of diagonal entries is the smaller of the
// number of rows and columns); it agrees with Matrix.h's square-only version
// whenever that one applies. M is modified in place, as upstream.
export function makeDiagonal(D: Vector, M: Matrix): void {
    const numDiagonal = (M.numRows <= M.numCols ? M.numRows : M.numCols);
    M.makeZero();
    for (let i = 0; i < numDiagonal; ++i) {
        M.set(i, i, D.values[i]);
    }
}

// Create an (N+1)-by-(N+1) matrix H by setting the upper N-by-N block to the
// input N-by-N matrix and all other entries to 0 except for the last row and
// last column entry which is set to 1. (Named hliftMatrix because Vector.h
// owns the name 'hlift'.)
export function hliftMatrix(M: Matrix): Matrix {
    logAssert(M.numRows === M.numCols, 'Matrix must be square.');
    const n = M.numRows;
    const result = new Matrix(n + 1, n + 1);
    result.makeIdentity();
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            result.set(r, c, M.get(r, c));
        }
    }
    return result;
}

// Extract the upper (N-1)-by-(N-1) block of the input N-by-N matrix. (Named
// hprojectMatrix because Vector.h owns the name 'hproject'.)
export function hprojectMatrix(M: Matrix): Matrix {
    logAssert(M.numRows === M.numCols, 'Matrix must be square.');
    const n = M.numRows;
    logAssert(n >= 2, 'Invalid matrix dimension.');
    const result = new Matrix(n - 1, n - 1);
    for (let r = 0; r < n - 1; ++r) {
        for (let c = 0; c < n - 1; ++c) {
            result.set(r, c, M.get(r, c));
        }
    }
    return result;
}
