// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) GMatrix.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: upstream GMatrix<Real> duplicates the entire
// Matrix<NumRows, NumCols, Real> API with a std::vector-backed, resizable
// table. The port's Matrix class already has runtime dimensions, so GMatrix
// is a subclass that adds only the resizing behavior (setSize), the
// range-checked element access and the dimension-aware comparison semantics
// of upstream GMatrix.
//
// All of Matrix's module free functions (negateMatrix, addMatrix, subMatrix,
// mulMatrix, divMatrix, l1Norm, l2Norm, lInfinityNorm, inverse, determinant,
// transpose, multiplyAB, multiplyABT, multiplyATB, multiplyATBT, multiplyMD,
// multiplyDM, outerProduct, makeDiagonal) apply to GMatrix unchanged and are
// NOT duplicated here; import them from './Matrix.js' (upstream GMatrix.h
// re-implements them with identical algorithms). They accept GMatrix inputs
// and return plain Matrix results; call GMatrix.fromMatrix(result) when a
// resizable matrix is needed.
//
// Behavioral differences from Matrix, matching upstream GMatrix.h:
// - The default constructor produces the 0-by-0 matrix (upstream GMatrix()).
// - setSize(numRows, numCols) resizes the table; a nonpositive numRows or
//   numCols clears it to 0-by-0. New elements are zero (the port of
//   std::vector::resize, which value-initializes).
// - get/set throw "Invalid index." for an out-of-range (r,c); the fixed-size
//   Matrix accessors are unchecked, as upstream.
// - makeUnit throws "Invalid index." for an invalid (r,c); the fixed-size
//   Matrix makeUnit silently produces the zero matrix, as upstream.
// - setRow/setCol/getRow/getCol validate the index ("Invalid index.") and the
//   vector size ("Mismatched sizes."); the fixed-size versions do not, since
//   upstream checks those at compile time.
// - Comparisons require equal dimensions: matrices of different dimensions
//   are unequal and are ordered by neither < nor >, whereas the fixed-size
//   Matrix comparisons throw on mismatched dimensions.
//
// Known deviation (documented rather than duplicated, since the free
// function names are owned by Matrix.ts): 'M / scalar' with scalar 0 throws
// upstream ("Division by zero.") but the shared divMatrix() from Matrix.ts
// returns the zero matrix (Matrix.h behavior).

import { logError } from './Logger.js';
import { Matrix } from './Matrix.js';
import { Vector } from './Vector.js';

export class GMatrix extends Matrix {
    // Create a numRows-by-numCols matrix (default 0-by-0, matching the
    // upstream default constructor) with all elements 0. When r and c are
    // provided, element (r,c) is set to 1 (the port of
    // GMatrix(numRows, numCols, r, c)).
    constructor(numRows: number = 0, numCols: number = 0, r?: number,
        c?: number) {
        const valid = numRows > 0 && numCols > 0;
        super(valid ? numRows : 0, valid ? numCols : 0);
        if (r !== undefined && c !== undefined) {
            this.makeUnit(r, c);
        }
    }

    // Create a GMatrix whose elements are the input values, specified in
    // row-major order. At most numRows*numCols values are copied; any
    // remaining elements are zero.
    static override fromArray(numRows: number, numCols: number,
        values: readonly number[]): GMatrix {
        const result = new GMatrix(numRows, numCols);
        const numValues = Math.min(values.length, result.numElements);
        for (let i = 0; i < numValues; ++i) {
            result.values[i] = values[i];
        }
        return result;
    }

    // Create a GMatrix that copies the elements of any Matrix. This is a
    // port convenience: the shared free functions in Matrix.ts return plain
    // Matrix objects.
    static fromMatrix(mat: Matrix): GMatrix {
        return GMatrix.fromArray(mat.numRows, mat.numCols, mat.values);
    }

    // All elements are 0.
    static override zero(numRows: number, numCols: number): GMatrix {
        return new GMatrix(numRows, numCols);
    }

    // Element (r,c) is 1, all others are 0.
    static override unit(numRows: number, numCols: number, r: number,
        c: number): GMatrix {
        const result = new GMatrix(numRows, numCols);
        result.makeUnit(r, c);
        return result;
    }

    // Diagonal entries 1, others 0, even when nonsquare.
    static override identity(numRows: number, numCols: number): GMatrix {
        const result = new GMatrix(numRows, numCols);
        result.makeIdentity();
        return result;
    }

    // Resize the table. When numRows and numCols are both positive, the
    // table has numRows*numCols elements; otherwise the matrix becomes
    // 0-by-0. Upstream uses std::vector::resize, which preserves the leading
    // elements of the old table (whose row/column meaning changes when the
    // number of columns changes) and value-initializes the new ones; the
    // port does the same, zeroing the new elements.
    setSize(numRows: number, numCols: number): void {
        if (numRows > 0 && numCols > 0) {
            this.mNumRows = numRows;
            this.mNumCols = numCols;
            const newSize = numRows * numCols;
            const oldSize = this.values.length;
            this.values.length = newSize;
            for (let i = oldSize; i < newSize; ++i) {
                this.values[i] = 0;
            }
        } else {
            this.mNumRows = 0;
            this.mNumCols = 0;
            this.values.length = 0;
        }
    }

    // A deep copy (the port of C++ copy construction/assignment).
    override clone(): GMatrix {
        return GMatrix.fromArray(this.mNumRows, this.mNumCols, this.values);
    }

    // Range-checked element access (upstream GMatrix operator()).
    override get(r: number, c: number): number {
        if (0 <= r && r < this.mNumRows && 0 <= c && c < this.mNumCols) {
            return this.values[c + this.mNumCols * r];
        }
        return logError('Invalid index.');
    }

    override set(r: number, c: number, value: number): void {
        if (0 <= r && r < this.mNumRows && 0 <= c && c < this.mNumCols) {
            this.values[c + this.mNumCols * r] = value;
            return;
        }
        logError('Invalid index.');
    }

    // Member access by rows or by columns. The input vectors must have the
    // correct number of elements for the matrix size.
    override setRow(r: number, vec: Vector): void {
        if (0 <= r && r < this.mNumRows) {
            if (vec.size === this.mNumCols) {
                for (let c = 0; c < this.mNumCols; ++c) {
                    this.set(r, c, vec.values[c]);
                }
                return;
            }
            logError('Mismatched sizes.');
        }
        logError('Invalid index.');
    }

    override setCol(c: number, vec: Vector): void {
        if (0 <= c && c < this.mNumCols) {
            if (vec.size === this.mNumRows) {
                for (let r = 0; r < this.mNumRows; ++r) {
                    this.set(r, c, vec.values[r]);
                }
                return;
            }
            logError('Mismatched sizes.');
        }
        logError('Invalid index.');
    }

    override getRow(r: number): Vector {
        if (0 <= r && r < this.mNumRows) {
            const vec = new Vector(this.mNumCols);
            for (let c = 0; c < this.mNumCols; ++c) {
                vec.values[c] = this.get(r, c);
            }
            return vec;
        }
        return logError('Invalid index.');
    }

    override getCol(c: number): Vector {
        if (0 <= c && c < this.mNumCols) {
            const vec = new Vector(this.mNumRows);
            for (let r = 0; r < this.mNumRows; ++r) {
                vec.values[r] = this.get(r, c);
            }
            return vec;
        }
        return logError('Invalid index.');
    }

    // Element (r,c) is 1, all others zero. Unlike the fixed-size Matrix
    // makeUnit, an invalid (r,c) throws.
    override makeUnit(r: number, c: number): void {
        if (0 <= r && r < this.mNumRows && 0 <= c && c < this.mNumCols) {
            this.makeZero();
            this.values[c + this.mNumCols * r] = 1;
            return;
        }
        logError('Invalid index.');
    }

    // Comparisons (for use by sorted containers). Unlike the fixed-size
    // Matrix comparisons, which throw on mismatched dimensions, these follow
    // upstream GMatrix: the dimensions must agree for any of the comparisons
    // to hold, so matrices of different dimensions are unequal and are
    // ordered by neither < nor > (an upstream quirk that is not a strict weak
    // ordering; preserved here).
    private sameSize(mat: Matrix): boolean {
        return this.mNumRows === mat.numRows && this.mNumCols === mat.numCols;
    }

    private elementCompare(mat: Matrix): number {
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

    override equals(mat: Matrix): boolean {
        return this.sameSize(mat) && this.elementCompare(mat) === 0;
    }

    override notEquals(mat: Matrix): boolean {
        return !this.equals(mat);
    }

    override lessThan(mat: Matrix): boolean {
        return this.sameSize(mat) && this.elementCompare(mat) < 0;
    }

    override lessThanOrEqual(mat: Matrix): boolean {
        return this.sameSize(mat) && this.elementCompare(mat) <= 0;
    }

    override greaterThan(mat: Matrix): boolean {
        return this.sameSize(mat) && this.elementCompare(mat) > 0;
    }

    override greaterThanOrEqual(mat: Matrix): boolean {
        return this.sameSize(mat) && this.elementCompare(mat) >= 0;
    }
}
