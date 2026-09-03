// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) LDLTDecomposition.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Factor a positive symmetric matrix A = L * D * L^T, where L is a lower
// triangular matrix with diagonal entries all 1 (L is lower unit triangular)
// and where D is a diagonal matrix with diagonal entries all positive.
//
// Port notes:
// - Upstream has two specializations of each class, one with the sizes known
//   at compile time (Matrix<N,N,T>, std::array-of-blocks) and one with the
//   sizes known only at run time (GMatrix<T>, std::vector-of-blocks). The
//   port's Matrix already carries runtime dimensions (see Matrix.ts), so the
//   two specializations collapse into a single class whose size is a
//   constructor argument. The run-time specialization's LogAssert size
//   validation is kept, since it is the stricter of the two.
// - Upstream's output reference parameters become fields of a returned
//   object: Factor(A, L, D) -> factor(A) returning { success, L, D }.
// - The two Solve overloads become distinct methods, per the port's naming
//   rules: Solve(L, D, B, X) -> solveFactored(L, D, B) returning X, and
//   Solve(A, B, X) -> solve(A, B) returning { success, X }.
// - BlockMatrix is the run-time layout: a flat array of numBlocks*numBlocks
//   blocks in row-major order, block (row,col) being element
//   'col + row * numBlocks'. Each block is a blockSize-by-blockSize matrix.
//   BlockVector is an array of numBlocks vectors, each of blockSize
//   components.
// - Upstream's 'verifySize' flags on the run-time block methods are ported as
//   optional boolean arguments that default to true, as upstream.
// - Although the class documentation says A must be positive definite, the
//   scalar factorization fails only on an exactly zero pivot; a symmetric
//   indefinite matrix with nonzero leading principal minors factors fine and
//   produces a D with negative diagonal entries. The port preserves that.
// - Upstream bug (fixed here): the run-time
//   BlockLDLTDecomposition<T>::Convert(BlockVector, GVector) verifies each
//   block vector with 'current.GetSize() == NumBlocks'; the block vectors
//   have BlockSize components, so the assertion rejects valid input whenever
//   BlockSize != NumBlocks. The port checks blockSize.

import { logAssert } from './Logger.js';
import {
    Matrix, inverse, mulMatrix, multiplyABT, subMatrix
} from './Matrix.js';
import { Vector, sub } from './Vector.js';

// An array of numBlocks*numBlocks blockSize-by-blockSize matrices, in
// row-major order over the blocks.
export type LDLTBlockMatrix = Matrix[];

// An array of numBlocks vectors, each with blockSize components.
export type LDLTBlockVector = Vector[];

export class LDLTDecomposition {
    // The size of the square matrices this object factors.
    readonly N: number;

    constructor(n: number) {
        logAssert(n > 0, 'Invalid size.');
        this.N = n;
    }

    // The matrix A must be positive definite. The implementation uses only
    // the lower-triangular portion of A. On output, L is lower unit
    // triangular and D is diagonal. When a pivot D(j,j) is zero, success is
    // false and L and D contain the partial results computed so far
    // (upstream behavior).
    factor(A: Matrix): { success: boolean, L: Matrix, D: Matrix } {
        const n = this.N;
        logAssert(A.numRows === n && A.numCols === n, 'Invalid size.');

        const L = new Matrix(n, n);
        const D = new Matrix(n, n);

        for (let j = 0; j < n; ++j) {
            let Djj = A.get(j, j);
            for (let k = 0; k < j; ++k) {
                const Ljk = L.get(j, k);
                const Dkk = D.get(k, k);
                Djj -= Ljk * Ljk * Dkk;
            }
            D.set(j, j, Djj);
            if (Djj === 0) {
                return { success: false, L, D };
            }

            L.set(j, j, 1);
            for (let i = j + 1; i < n; ++i) {
                let Lij = A.get(i, j);
                for (let k = 0; k < j; ++k) {
                    const Lik = L.get(i, k);
                    const Ljk = L.get(j, k);
                    const Dkk = D.get(k, k);
                    Lij -= Lik * Ljk * Dkk;
                }

                Lij /= Djj;
                L.set(i, j, Lij);
            }
        }
        return { success: true, L, D };
    }

    // Solve A*X = B for positive definite A = L * D * L^T with factoring
    // before the call.
    solveFactored(L: Matrix, D: Matrix, B: Vector): Vector {
        const n = this.N;
        logAssert(L.numRows === n && L.numCols === n
            && D.numRows === n && D.numCols === n && B.size === n,
            'Invalid size.');

        const X = new Vector(n);

        // Solve L * Z = L * (D * L^T * X) = B for Z.
        for (let r = 0; r < n; ++r) {
            let value = B.get(r);
            for (let c = 0; c < r; ++c) {
                value -= L.get(r, c) * X.get(c);
            }
            X.set(r, value);
        }

        // Solve D * Y = D * (L^T * X) = Z for Y.
        for (let r = 0; r < n; ++r) {
            X.set(r, X.get(r) / D.get(r, r));
        }

        // Solve L^T * Y = Z for X.
        for (let r = n - 1; r >= 0; --r) {
            let value = X.get(r);
            for (let c = r + 1; c < n; ++c) {
                value -= L.get(c, r) * X.get(c);
            }
            X.set(r, value);
        }

        return X;
    }

    // Solve A*X = B for positive semidefinite A = L * D * L^T with factoring
    // during the call. If A has a zero eigenvalue, the factoring fails, so
    // the returned 'success' indicates whether the solver succeeded. On
    // failure X is the zero vector.
    solve(A: Matrix, B: Vector): { success: boolean, X: Vector } {
        const n = this.N;
        logAssert(A.numRows === n && A.numCols === n && B.size === n,
            'Invalid size.');

        const { success, L, D } = this.factor(A);
        if (success) {
            return { success, X: this.solveFactored(L, D, B) };
        }
        return { success, X: new Vector(n) };
    }
}

// Let B represent the block size and N represent the number of blocks. The
// matrix A is (N*B)-by-(N*B) but partitioned into an N-by-N matrix of blocks,
// each block of size B-by-B and stored in row-major order. The value N*B is
// numDimensions.
export class BlockLDLTDecomposition {
    readonly blockSize: number;
    readonly numBlocks: number;
    readonly numDimensions: number;

    constructor(blockSize: number, numBlocks: number) {
        logAssert(blockSize > 0 && numBlocks > 0, 'Invalid size.');
        this.blockSize = blockSize;
        this.numBlocks = numBlocks;
        this.numDimensions = blockSize * numBlocks;
    }

    // Compute the 1-dimensional index of the block matrix in a 2-dimensional
    // BlockMatrix object.
    private getIndex(row: number, col: number): number {
        return col + row * this.numBlocks;
    }

    // Treating the matrix as a 2D table of scalars with numDimensions rows
    // and numDimensions columns, look up the correct block that stores the
    // requested element and return its value. NOTE: You are responsible for
    // ensuring that M has numBlocks-by-numBlocks elements, each of them
    // blockSize-by-blockSize.
    get(M: LDLTBlockMatrix, row: number, col: number,
        verifySize: boolean = true): number {
        const numBlocks = this.numBlocks;
        if (verifySize) {
            logAssert(M.length === numBlocks * numBlocks, 'Invalid size.');
        }

        const b = this.blockSize;
        const b0 = Math.floor(col / b);
        const b1 = Math.floor(row / b);
        const i0 = col - b * b0;
        const i1 = row - b * b1;
        const MBlock = M[this.getIndex(b1, b0)];

        if (verifySize) {
            logAssert(MBlock.numRows === b && MBlock.numCols === b,
                'Invalid size.');
        }

        return MBlock.get(i1, i0);
    }

    set(M: LDLTBlockMatrix, row: number, col: number, value: number,
        verifySize: boolean = true): void {
        const numBlocks = this.numBlocks;
        if (verifySize) {
            logAssert(M.length === numBlocks * numBlocks, 'Invalid size.');
        }

        const b = this.blockSize;
        const b0 = Math.floor(col / b);
        const b1 = Math.floor(row / b);
        const i0 = col - b * b0;
        const i1 = row - b * b1;
        const MBlock = M[this.getIndex(b1, b0)];

        if (verifySize) {
            logAssert(MBlock.numRows === b && MBlock.numCols === b,
                'Invalid size.');
        }

        MBlock.set(i1, i0, value);
    }

    // Convert from a matrix to a block matrix.
    convertMatrixToBlock(M: Matrix, verifySize: boolean = true):
        LDLTBlockMatrix {
        const numDim = this.numDimensions;
        if (verifySize) {
            logAssert(M.numRows === numDim && M.numCols === numDim,
                'Invalid size.');
        }

        const numBlocks = this.numBlocks, blockSize = this.blockSize;
        const MBlock: LDLTBlockMatrix = new Array<Matrix>(
            numBlocks * numBlocks);
        for (let r = 0, rb = 0, index = 0; r < numBlocks; ++r, rb += blockSize) {
            for (let c = 0, cb = 0; c < numBlocks; ++c, cb += blockSize, ++index) {
                const current = new Matrix(blockSize, blockSize);
                for (let j = 0; j < blockSize; ++j) {
                    for (let i = 0; i < blockSize; ++i) {
                        current.set(j, i, M.get(rb + j, cb + i));
                    }
                }
                MBlock[index] = current;
            }
        }
        return MBlock;
    }

    // Convert from a vector to a block vector.
    convertVectorToBlock(V: Vector, verifySize: boolean = true):
        LDLTBlockVector {
        if (verifySize) {
            logAssert(V.size === this.numDimensions, 'Invalid size.');
        }

        const numBlocks = this.numBlocks, blockSize = this.blockSize;
        const VBlock: LDLTBlockVector = new Array<Vector>(numBlocks);
        for (let r = 0, rb = 0; r < numBlocks; ++r, rb += blockSize) {
            const current = new Vector(blockSize);
            for (let j = 0; j < blockSize; ++j) {
                current.set(j, V.get(rb + j));
            }
            VBlock[r] = current;
        }
        return VBlock;
    }

    // Convert from a block matrix to a matrix.
    convertBlockToMatrix(MBlock: LDLTBlockMatrix,
        verifySize: boolean = true): Matrix {
        const numBlocks = this.numBlocks, blockSize = this.blockSize;
        if (verifySize) {
            logAssert(MBlock.length === numBlocks * numBlocks, 'Invalid size.');
            for (const current of MBlock) {
                logAssert(current.numRows === blockSize
                    && current.numCols === blockSize, 'Invalid size.');
            }
        }

        const M = new Matrix(this.numDimensions, this.numDimensions);
        for (let r = 0, rb = 0, index = 0; r < numBlocks; ++r, rb += blockSize) {
            for (let c = 0, cb = 0; c < numBlocks; ++c, cb += blockSize, ++index) {
                const current = MBlock[index];
                for (let j = 0; j < blockSize; ++j) {
                    for (let i = 0; i < blockSize; ++i) {
                        M.set(rb + j, cb + i, current.get(j, i));
                    }
                }
            }
        }
        return M;
    }

    // Convert from a block vector to a vector.
    convertBlockToVector(VBlock: LDLTBlockVector,
        verifySize: boolean = true): Vector {
        const numBlocks = this.numBlocks, blockSize = this.blockSize;
        if (verifySize) {
            logAssert(VBlock.length === numBlocks, 'Invalid size.');
            for (const current of VBlock) {
                // Upstream compares with NumBlocks here; see the port notes.
                logAssert(current.size === blockSize, 'Invalid size.');
            }
        }

        const V = new Vector(this.numDimensions);
        for (let r = 0, rb = 0; r < numBlocks; ++r, rb += blockSize) {
            const current = VBlock[r];
            for (let j = 0; j < blockSize; ++j) {
                V.set(rb + j, current.get(j));
            }
        }
        return V;
    }

    // The block matrix A must be positive definite. The implementation uses
    // only the lower-triangular blocks of A. On output, the block matrix L is
    // lower unit triangular (diagonal blocks are BxB identity matrices) and
    // the block matrix D is diagonal (diagonal blocks are BxB diagonal
    // matrices; the off-diagonal blocks of D are zero).
    factor(A: LDLTBlockMatrix, verifySize: boolean = true):
        { success: boolean, L: LDLTBlockMatrix, D: LDLTBlockMatrix } {
        const numBlocks = this.numBlocks, blockSize = this.blockSize;
        if (verifySize) {
            logAssert(A.length === numBlocks * numBlocks, 'Invalid size.');
            for (let i = 0; i < A.length; ++i) {
                logAssert(A[i].numRows === blockSize
                    && A[i].numCols === blockSize, 'Invalid size.');
            }
        }

        const L: LDLTBlockMatrix = new Array<Matrix>(A.length);
        const D: LDLTBlockMatrix = new Array<Matrix>(A.length);
        for (let i = 0; i < L.length; ++i) {
            L[i] = new Matrix(blockSize, blockSize);
            D[i] = new Matrix(blockSize, blockSize);
        }

        for (let j = 0; j < numBlocks; ++j) {
            let Djj = A[this.getIndex(j, j)].clone();
            for (let k = 0; k < j; ++k) {
                const Ljk = L[this.getIndex(j, k)];
                const Dkk = D[this.getIndex(k, k)];
                Djj = subMatrix(Djj, multiplyABT(mulMatrix(Ljk, Dkk), Ljk));
            }
            D[this.getIndex(j, j)] = Djj;
            const { inverse: invDjj, invertible } = inverse(Djj);
            if (!invertible) {
                return { success: false, L, D };
            }

            L[this.getIndex(j, j)].makeIdentity();
            for (let i = j + 1; i < numBlocks; ++i) {
                let Lij = A[this.getIndex(i, j)].clone();
                for (let k = 0; k < j; ++k) {
                    const Lik = L[this.getIndex(i, k)];
                    const Ljk = L[this.getIndex(j, k)];
                    const Dkk = D[this.getIndex(k, k)];
                    Lij = subMatrix(Lij, multiplyABT(mulMatrix(Lik, Dkk), Ljk));
                }
                Lij = mulMatrix(Lij, invDjj);
                L[this.getIndex(i, j)] = Lij;
            }
        }
        return { success: true, L, D };
    }

    // Solve A*X = B for positive definite A = L * D * L^T with factoring
    // before the call.
    solveFactored(L: LDLTBlockMatrix, D: LDLTBlockMatrix, B: LDLTBlockVector,
        verifySize: boolean = true): LDLTBlockVector {
        const numBlocks = this.numBlocks, blockSize = this.blockSize;
        if (verifySize) {
            const LDsize = numBlocks * numBlocks;
            logAssert(L.length === LDsize && D.length === LDsize
                && B.length === numBlocks, 'Invalid size.');

            for (let i = 0; i < L.length; ++i) {
                logAssert(L[i].numRows === blockSize
                    && L[i].numCols === blockSize
                    && D[i].numRows === blockSize
                    && D[i].numCols === blockSize, 'Invalid size.');
            }

            for (let i = 0; i < B.length; ++i) {
                logAssert(B[i].size === blockSize, 'Invalid size.');
            }
        }

        // Solve L * Z = L * (D * L^T * X) = B for Z.
        const X: LDLTBlockVector = new Array<Vector>(numBlocks);
        for (let r = 0; r < numBlocks; ++r) {
            X[r] = B[r].clone();
            for (let c = 0; c < r; ++c) {
                X[r] = sub(X[r], mulMatrix(L[this.getIndex(r, c)], X[c]));
            }
        }

        // Solve D * Y = D * (L^T * X) = Z for Y.
        for (let r = 0; r < numBlocks; ++r) {
            X[r] = mulMatrix(inverse(D[this.getIndex(r, r)]).inverse, X[r]);
        }

        // Solve L^T * Y = Z for X.
        for (let r = numBlocks - 1; r >= 0; --r) {
            for (let c = r + 1; c < numBlocks; ++c) {
                X[r] = sub(X[r], mulMatrix(X[c], L[this.getIndex(c, r)]));
            }
        }

        return X;
    }

    // Solve A*X = B for positive semidefinite A = L * D * L^T with factoring
    // during the call. If A has a zero eigenvalue, the factoring fails, so
    // the returned 'success' indicates whether the solver succeeded. On
    // failure X is an empty block vector.
    solve(A: LDLTBlockMatrix, B: LDLTBlockVector, verifySize: boolean = true):
        { success: boolean, X: LDLTBlockVector } {
        const numBlocks = this.numBlocks, blockSize = this.blockSize;
        if (verifySize) {
            logAssert(A.length === numBlocks * numBlocks
                && B.length === numBlocks, 'Invalid size.');

            for (let i = 0; i < A.length; ++i) {
                logAssert(A[i].numRows === blockSize
                    && A[i].numCols === blockSize, 'Invalid size.');
            }

            for (let i = 0; i < B.length; ++i) {
                logAssert(B[i].size === blockSize, 'Invalid size.');
            }
        }

        const { success, L, D } = this.factor(A, false);
        if (success) {
            return { success, X: this.solveFactored(L, D, B, false) };
        }
        return { success, X: [] };
    }
}
