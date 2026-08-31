// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BandedMatrix.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

import { LexicoArray2 } from './LexicoArray2';

// A square matrix stored as a diagonal band plus 'numLBands' subdiagonal
// bands and 'numUBands' superdiagonal bands. Entries outside the bands are
// zero and are not stored.
//
// Port notes:
// - Upstream's 'Real& operator()(int32_t r, int32_t c)' returns a reference
//   that is used for reading and for writing. TypeScript has no references,
//   so it is ported as the pair get(r, c) / set(r, c, value). Upstream
//   returns a reference to the scratch member mZero for indices that are in
//   range but outside the bands (and for out-of-range indices); writes
//   through that reference are discarded. The port reproduces that: get
//   returns 0 and set is a silent no-op for those indices.
// - The compile-time 'template <bool RowMajor>' parameter of the matrix-RHS
//   solvers and of ComputeInverse becomes a runtime 'rowMajor' argument,
//   matching the LexicoArray2 port. It defaults to true (row major), which
//   is the GTE default (GTE_USE_ROW_MAJOR).
// - Upstream's 'BandedMatrix<Real> tmpA = *this;' copy in ComputeInverse is
//   ported as an explicit clone(), because TypeScript objects alias.
export class BandedMatrix {
    private mSize: number;
    private mDBand: number[];
    private mLBands: number[][];
    private mUBands: number[][];

    // The constructor requires size > 0, 0 <= numLBands < size and
    // 0 <= numUBands < size. As upstream, invalid arguments produce a matrix
    // of size 0 rather than an exception.
    constructor(size: number, numLBands: number, numUBands: number) {
        this.mSize = size;
        this.mDBand = [];
        this.mLBands = [];
        this.mUBands = [];

        if (size > 0
            && 0 <= numLBands && numLBands < size
            && 0 <= numUBands && numUBands < size) {
            this.mDBand = new Array<number>(size).fill(0);

            if (numLBands > 0) {
                let numElements = size - 1;
                for (let i = 0; i < numLBands; ++i) {
                    this.mLBands.push(new Array<number>(numElements--).fill(0));
                }
            }

            if (numUBands > 0) {
                let numElements = size - 1;
                for (let i = 0; i < numUBands; ++i) {
                    this.mUBands.push(new Array<number>(numElements--).fill(0));
                }
            }
        } else {
            // Invalid argument to the BandedMatrix constructor.
            this.mSize = 0;
        }
    }

    // Member access.
    getSize(): number {
        return this.mSize;
    }

    getDBand(): number[] {
        return this.mDBand;
    }

    getLBands(): number[][] {
        return this.mLBands;
    }

    getUBands(): number[][] {
        return this.mUBands;
    }

    // A deep copy of the matrix.
    clone(): BandedMatrix {
        const copy = new BandedMatrix(this.mSize, this.mLBands.length,
            this.mUBands.length);
        copy.mDBand = this.mDBand.slice();
        copy.mLBands = this.mLBands.map(band => band.slice());
        copy.mUBands = this.mUBands.map(band => band.slice());
        return copy;
    }

    // Read the (r, c) entry. In-range indices outside the stored bands and
    // out-of-range indices produce 0 (upstream returns a reference to its
    // mZero scratch member).
    get(r: number, c: number): number {
        if (0 <= r && r < this.mSize && 0 <= c && c < this.mSize) {
            let band = c - r;
            if (band > 0) {
                const numUBands = this.mUBands.length;
                if (--band < numUBands && r < this.mSize - 1 - band) {
                    return this.mUBands[band][r];
                }
            } else if (band < 0) {
                band = -band;
                const numLBands = this.mLBands.length;
                if (--band < numLBands && c < this.mSize - 1 - band) {
                    return this.mLBands[band][c];
                }
            } else {
                return this.mDBand[r];
            }
        }
        // The entry is zero and is not stored.
        return 0;
    }

    // Write the (r, c) entry. In-range indices outside the stored bands and
    // out-of-range indices discard the value (upstream writes through a
    // reference to its mZero scratch member).
    set(r: number, c: number, value: number): void {
        if (0 <= r && r < this.mSize && 0 <= c && c < this.mSize) {
            let band = c - r;
            if (band > 0) {
                const numUBands = this.mUBands.length;
                if (--band < numUBands && r < this.mSize - 1 - band) {
                    this.mUBands[band][r] = value;
                }
            } else if (band < 0) {
                band = -band;
                const numLBands = this.mLBands.length;
                if (--band < numLBands && c < this.mSize - 1 - band) {
                    this.mLBands[band][c] = value;
                }
            } else {
                this.mDBand[r] = value;
            }
        }
    }

    // Factor the square banded matrix A into A = L*L^T, where L is a
    // lower-triangular matrix (L^T is an upper-triangular matrix). This is an
    // LU decomposition that allows for stable inversion of A to solve A*X = B.
    // The return value is true iff the factoring is successful (L is
    // invertible). If successful, the matrix contains the Cholesky
    // factorization: L in the lower-triangular part and L^T in the
    // upper-triangular part.
    choleskyFactor(): boolean {
        if (this.mDBand.length === 0 || this.mLBands.length !== this.mUBands.length) {
            // Invalid number of bands.
            return false;
        }

        const sizeM1 = this.mSize - 1;
        const numBands = this.mLBands.length;

        let k: number, kMax: number;
        for (let i = 0; i < this.mSize; ++i) {
            let jMin = i - numBands;
            if (jMin < 0) {
                jMin = 0;
            }

            let j: number;
            for (j = jMin; j < i; ++j) {
                kMax = j + numBands;
                if (kMax > sizeM1) {
                    kMax = sizeM1;
                }

                for (k = i; k <= kMax; ++k) {
                    this.set(k, i, this.get(k, i) - this.get(i, j) * this.get(k, j));
                }
            }

            kMax = j + numBands;
            if (kMax > sizeM1) {
                kMax = sizeM1;
            }

            for (k = 0; k < i; ++k) {
                this.set(k, i, this.get(i, k));
            }

            const diagonal = this.get(i, i);
            if (diagonal <= 0) {
                return false;
            }
            const invSqrt = 1 / Math.sqrt(diagonal);
            for (k = i; k <= kMax; ++k) {
                this.set(k, i, this.get(k, i) * invSqrt);
            }
        }

        return true;
    }

    // Solve the linear system A*X = B, where A is an NxN banded matrix and B
    // is an Nx1 vector. The unknown X is also Nx1. The input to this function
    // is B. The output X is computed and stored in B. The return value is
    // true iff the system has a solution. The matrix A and the vector B are
    // both modified by this function. If successful, A contains the Cholesky
    // factorization.
    solveSystem(bVector: number[]): boolean {
        return this.choleskyFactor()
            && this.solveLowerVector(bVector)
            && this.solveUpperVector(bVector);
    }

    // Solve the linear system A*X = B, where A is an NxN banded matrix and B
    // is an NxM matrix stored flat with the ordering given by 'rowMajor'. The
    // unknown X is also NxM. The input to this function is B. The output X is
    // computed and stored in B. The return value is true iff the system has a
    // solution. The matrix A and the matrix B are both modified.
    solveSystemMatrix(bMatrix: number[], numBColumns: number,
        rowMajor: boolean = true): boolean {
        return this.choleskyFactor()
            && this.solveLowerMatrix(bMatrix, numBColumns, rowMajor)
            && this.solveUpperMatrix(bMatrix, numBColumns, rowMajor);
    }

    // Compute the inverse of the banded matrix. The return value is true when
    // the matrix is invertible, in which case the 'inverse' output is valid.
    // The return value is false when the matrix is not invertible, in which
    // case 'inverse' is invalid and should not be used. The 'inverse' array
    // must have mSize*mSize elements and uses the storage order given by
    // 'rowMajor'.
    computeInverse(inverse: number[], rowMajor: boolean = true): boolean {
        const invA = new LexicoArray2(rowMajor, this.mSize, this.mSize, inverse);

        const tmpA = this.clone();
        for (let row = 0; row < this.mSize; ++row) {
            for (let col = 0; col < this.mSize; ++col) {
                if (row !== col) {
                    invA.set(row, col, 0);
                } else {
                    invA.set(row, row, 1);
                }
            }
        }

        // Forward elimination.
        for (let row = 0; row < this.mSize; ++row) {
            // The pivot must be nonzero in order to proceed.
            const diag = tmpA.get(row, row);
            if (diag === 0) {
                return false;
            }

            const invDiag = 1 / diag;
            tmpA.set(row, row, 1);

            // Multiply the row to be consistent with a diagonal term of 1.
            const colMin = row + 1;
            let colMax = colMin + this.mUBands.length;
            if (colMax > this.mSize) {
                colMax = this.mSize;
            }

            let c: number;
            for (c = colMin; c < colMax; ++c) {
                tmpA.set(row, c, tmpA.get(row, c) * invDiag);
            }
            for (c = 0; c <= row; ++c) {
                invA.set(row, c, invA.get(row, c) * invDiag);
            }

            // Reduce the remaining rows.
            const rowMin = row + 1;
            let rowMax = rowMin + this.mLBands.length;
            if (rowMax > this.mSize) {
                rowMax = this.mSize;
            }

            for (let r = rowMin; r < rowMax; ++r) {
                const mult = tmpA.get(r, row);
                tmpA.set(r, row, 0);
                for (c = colMin; c < colMax; ++c) {
                    tmpA.set(r, c, tmpA.get(r, c) - mult * tmpA.get(row, c));
                }
                for (c = 0; c <= row; ++c) {
                    invA.set(r, c, invA.get(r, c) - mult * invA.get(row, c));
                }
            }
        }

        // Backward elimination.
        for (let row = this.mSize - 1; row >= 1; --row) {
            const rowMax = row - 1;
            let rowMin = row - this.mUBands.length;
            if (rowMin < 0) {
                rowMin = 0;
            }

            for (let r = rowMax; r >= rowMin; --r) {
                const mult = tmpA.get(r, row);
                tmpA.set(r, row, 0);
                for (let c = 0; c < this.mSize; ++c) {
                    invA.set(r, c, invA.get(r, c) - mult * invA.get(row, c));
                }
            }
        }

        return true;
    }

    // The linear system is L*U*X = B, where A = L*U and U = L^T. Reduce this
    // to U*X = L^{-1}*B. The return value is true iff the operation is
    // successful.
    private solveLowerVector(dataVector: number[]): boolean {
        const size = this.mDBand.length;
        for (let r = 0; r < size; ++r) {
            const lowerRR = this.get(r, r);
            if (lowerRR > 0) {
                for (let c = 0; c < r; ++c) {
                    const lowerRC = this.get(r, c);
                    dataVector[r] -= lowerRC * dataVector[c];
                }
                dataVector[r] /= lowerRR;
            } else {
                return false;
            }
        }
        return true;
    }

    // The linear system is U*X = L^{-1}*B. Reduce this to
    // X = U^{-1}*L^{-1}*B. The return value is true iff the operation is
    // successful.
    private solveUpperVector(dataVector: number[]): boolean {
        const size = this.mDBand.length;
        for (let r = size - 1; r >= 0; --r) {
            const upperRR = this.get(r, r);
            if (upperRR > 0) {
                for (let c = r + 1; c < size; ++c) {
                    const upperRC = this.get(r, c);
                    dataVector[r] -= upperRC * dataVector[c];
                }
                dataVector[r] /= upperRR;
            } else {
                return false;
            }
        }
        return true;
    }

    // The matrix-right-hand-side counterpart of solveLowerVector.
    private solveLowerMatrix(dataMatrix: number[], numColumns: number,
        rowMajor: boolean): boolean {
        const data = new LexicoArray2(rowMajor, this.mSize, numColumns, dataMatrix);

        for (let r = 0; r < this.mSize; ++r) {
            const lowerRR = this.get(r, r);
            if (lowerRR > 0) {
                for (let c = 0; c < r; ++c) {
                    const lowerRC = this.get(r, c);
                    for (let bCol = 0; bCol < numColumns; ++bCol) {
                        data.set(r, bCol, data.get(r, bCol) - lowerRC * data.get(c, bCol));
                    }
                }

                const inverse = 1 / lowerRR;
                for (let bCol = 0; bCol < numColumns; ++bCol) {
                    data.set(r, bCol, data.get(r, bCol) * inverse);
                }
            } else {
                return false;
            }
        }
        return true;
    }

    // The matrix-right-hand-side counterpart of solveUpperVector.
    private solveUpperMatrix(dataMatrix: number[], numColumns: number,
        rowMajor: boolean): boolean {
        const data = new LexicoArray2(rowMajor, this.mSize, numColumns, dataMatrix);

        for (let r = this.mSize - 1; r >= 0; --r) {
            const upperRR = this.get(r, r);
            if (upperRR > 0) {
                for (let c = r + 1; c < this.mSize; ++c) {
                    const upperRC = this.get(r, c);
                    for (let bCol = 0; bCol < numColumns; ++bCol) {
                        data.set(r, bCol, data.get(r, bCol) - upperRC * data.get(c, bCol));
                    }
                }

                const inverse = 1 / upperRR;
                for (let bCol = 0; bCol < numColumns; ++bCol) {
                    data.set(r, bCol, data.get(r, bCol) * inverse);
                }
            } else {
                return false;
            }
        }
        return true;
    }
}
