// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) CholeskyDecomposition.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Factor a positive definite symmetric matrix A = L * L^T, where L is lower
// triangular.
//
// Port notes:
// - Upstream has two specializations of each class, one with the sizes known
//   at compile time (Matrix<N,N,Real>, std::array-of-blocks) and one with the
//   sizes known only at run time (GMatrix<Real>, std::vector-of-blocks). The
//   port's Matrix already carries runtime dimensions (see Matrix.ts), so the
//   two specializations collapse into a single class whose size is a
//   constructor argument, exactly as the fixed-size Vector/Matrix templates
//   collapsed in the foundation port. The run-time specialization's size
//   validation (LogError/LogAssert) is kept, since it is the stricter of the
//   two.
// - BlockMatrix is the run-time layout: a flat array of numBlocks*numBlocks
//   blocks in row-major order, block (row,col) being element
//   'col + row * numBlocks'. Each block is a blockSize-by-blockSize matrix.
//   BlockVector is an array of numBlocks vectors, each of blockSize
//   components.
// - Factor and the solvers work in place, as upstream: Factor overwrites the
//   lower-triangular portion of A with L, and solveLower/solveUpper overwrite
//   their vector argument with the solution.
// - Upstream bug (fixed here): the run-time specialization
//   BlockCholeskyDecomposition<Real,0,0> indexes the scalars *inside* a block
//   with its block-level helper GetIndex(row,col) = col + row*NumBlocks (see
//   SolveLower, SolveUpper, LowerTriangularSolver and SubtractiveUpdate). The
//   linear index of element (row,col) of a blockSize-by-blockSize block is
//   col + row*BlockSize, so the run-time results are wrong whenever
//   BlockSize != NumBlocks. The compile-time specialization uses the correct
//   'block(row,col)' access; the port follows it and uses get/set on the
//   block throughout.

import { logAssert, logError } from './Logger';
import { Matrix } from './Matrix';
import { Vector } from './Vector';

// An array of numBlocks*numBlocks blockSize-by-blockSize matrices, in
// row-major order over the blocks.
export type CholeskyBlockMatrix = Matrix[];

// An array of numBlocks vectors, each with blockSize components.
export type CholeskyBlockVector = Vector[];

export class CholeskyDecomposition {
    // The size of the square matrices this object factors.
    readonly N: number;

    constructor(n: number) {
        logAssert(n > 0, 'Invalid size.');
        this.N = n;
    }

    // On input, A is symmetric. Only the lower-triangular portion is
    // modified. On output, the lower-triangular portion is L where
    // A = L * L^T. The return value is false when a diagonal entry of the
    // partially factored matrix is not positive, in which case A is left
    // partially modified (upstream behavior).
    factor(A: Matrix): boolean {
        const n = this.N;
        if (A.numRows !== n || A.numCols !== n) {
            return logError('Matrix must be square.');
        }

        for (let c = 0; c < n; ++c) {
            if (A.get(c, c) <= 0) {
                return false;
            }
            A.set(c, c, Math.sqrt(A.get(c, c)));

            for (let r = c + 1; r < n; ++r) {
                A.set(r, c, A.get(r, c) / A.get(c, c));
            }

            for (let k = c + 1; k < n; ++k) {
                for (let r = k; r < n; ++r) {
                    A.set(r, k, A.get(r, k) - A.get(r, c) * A.get(k, c));
                }
            }
        }
        return true;
    }

    // Solve L*Y = B, where L is lower triangular and invertible. The input
    // value of Y is B. On output, Y is the solution.
    solveLower(L: Matrix, Y: Vector): void {
        const n = this.N;
        if (L.numRows !== n || L.numCols !== n || Y.size !== n) {
            return logError('Invalid size.');
        }

        for (let r = 0; r < n; ++r) {
            let value = Y.get(r);
            for (let c = 0; c < r; ++c) {
                value -= L.get(r, c) * Y.get(c);
            }
            Y.set(r, value / L.get(r, r));
        }
    }

    // Solve L^T*X = Y, where L is lower triangular (L^T is upper triangular)
    // and invertible. The input value of X is Y. On output, X is the
    // solution.
    solveUpper(L: Matrix, X: Vector): void {
        const n = this.N;
        if (L.numRows !== n || L.numCols !== n || X.size !== n) {
            return logError('Invalid size.');
        }

        for (let r = n - 1; r >= 0; --r) {
            let value = X.get(r);
            for (let c = r + 1; c < n; ++c) {
                value -= L.get(c, r) * X.get(c);
            }
            X.set(r, value / L.get(r, r));
        }
    }
}

// Let B represent the block size and N represent the number of blocks. The
// matrix A is (N*B)-by-(N*B) but partitioned into an N-by-N matrix of blocks,
// each block of size B-by-B. The value N*B is numDimensions.
export class BlockCholeskyDecomposition {
    readonly blockSize: number;
    readonly numBlocks: number;
    readonly numDimensions: number;

    // The decomposer has size blockSize.
    private readonly mDecomposer: CholeskyDecomposition;

    constructor(blockSize: number, numBlocks: number) {
        logAssert(blockSize > 0 && numBlocks > 0, 'Invalid input.');
        this.blockSize = blockSize;
        this.numBlocks = numBlocks;
        this.numDimensions = blockSize * numBlocks;
        this.mDecomposer = new CholeskyDecomposition(blockSize);
    }

    // Compute the 1-dimensional index of the block matrix in a 2-dimensional
    // BlockMatrix object.
    private getIndex(row: number, col: number): number {
        return col + row * this.numBlocks;
    }

    // Treating the matrix as a 2D table of scalars with numDimensions rows
    // and numDimensions columns, look up the correct block that stores the
    // requested element and return its value.
    get(M: CholeskyBlockMatrix, row: number, col: number): number {
        const b = this.blockSize;
        const b0 = Math.floor(col / b), b1 = Math.floor(row / b);
        const i0 = col % b, i1 = row % b;
        const block = M[this.getIndex(b1, b0)];
        return block.get(i1, i0);
    }

    set(M: CholeskyBlockMatrix, row: number, col: number, value: number): void {
        const b = this.blockSize;
        const b0 = Math.floor(col / b), b1 = Math.floor(row / b);
        const i0 = col % b, i1 = row % b;
        const block = M[this.getIndex(b1, b0)];
        block.set(i1, i0, value);
    }

    // Factor the block matrix in place. On output, the lower-triangular
    // blocks store L, where A = L * L^T. The diagonal blocks store the
    // lower-triangular Cholesky factor of the corresponding Schur complement.
    factor(A: CholeskyBlockMatrix): boolean {
        const numBlocks = this.numBlocks;
        for (let c = 0; c < numBlocks; ++c) {
            if (!this.mDecomposer.factor(A[this.getIndex(c, c)])) {
                return false;
            }

            for (let r = c + 1; r < numBlocks; ++r) {
                this.lowerTriangularSolver(r, c, A);
            }

            for (let k = c + 1; k < numBlocks; ++k) {
                for (let r = k; r < numBlocks; ++r) {
                    this.subtractiveUpdate(r, k, c, A);
                }
            }
        }
        return true;
    }

    // Solve L*Y = B, where L is an invertible lower-triangular block matrix
    // whose diagonal blocks are lower-triangular matrices. The input B is a
    // block vector of commensurate size. The input value of Y is B. On
    // output, Y is the solution.
    solveLower(L: CholeskyBlockMatrix, Y: CholeskyBlockVector): void {
        const numBlocks = this.numBlocks, blockSize = this.blockSize;
        for (let r = 0; r < numBlocks; ++r) {
            const Yr = Y[r];
            for (let c = 0; c < r; ++c) {
                const Lrc = L[this.getIndex(r, c)];
                const Yc = Y[c];
                for (let i = 0; i < blockSize; ++i) {
                    let value = Yr.get(i);
                    for (let j = 0; j < blockSize; ++j) {
                        value -= Lrc.get(i, j) * Yc.get(j);
                    }
                    Yr.set(i, value);
                }
            }
            this.mDecomposer.solveLower(L[this.getIndex(r, r)], Yr);
        }
    }

    // Solve L^T*X = Y, where L is an invertible lower-triangular block matrix
    // (L^T is an upper-triangular block matrix) whose diagonal blocks are
    // lower-triangular matrices. The input value of X is Y. On output, X is
    // the solution.
    solveUpper(L: CholeskyBlockMatrix, X: CholeskyBlockVector): void {
        const numBlocks = this.numBlocks, blockSize = this.blockSize;
        for (let r = numBlocks - 1; r >= 0; --r) {
            const Xr = X[r];
            for (let c = r + 1; c < numBlocks; ++c) {
                const Lcr = L[this.getIndex(c, r)];
                const Xc = X[c];
                for (let i = 0; i < blockSize; ++i) {
                    let value = Xr.get(i);
                    for (let j = 0; j < blockSize; ++j) {
                        value -= Lcr.get(j, i) * Xc.get(j);
                    }
                    Xr.set(i, value);
                }
            }
            this.mDecomposer.solveUpper(L[this.getIndex(r, r)], Xr);
        }
    }

    // Solve G(c,c)*G(r,c)^T = A(r,c)^T for G(r,c). The matrices G(c,c) and
    // A(r,c) are known quantities, and G(c,c) occupies the lower triangular
    // portion of A(c,c). The solver stores its results in-place, so A(r,c)
    // stores the G(r,c) result.
    private lowerTriangularSolver(r: number, c: number,
        A: CholeskyBlockMatrix): void {
        const blockSize = this.blockSize;
        const Acc = A[this.getIndex(c, c)];
        const Arc = A[this.getIndex(r, c)];
        for (let j = 0; j < blockSize; ++j) {
            for (let i = 0; i < j; ++i) {
                const Lji = Acc.get(j, i);
                for (let k = 0; k < blockSize; ++k) {
                    Arc.set(k, j, Arc.get(k, j) - Lji * Arc.get(k, i));
                }
            }

            const Ljj = Acc.get(j, j);
            for (let k = 0; k < blockSize; ++k) {
                Arc.set(k, j, Arc.get(k, j) / Ljj);
            }
        }
    }

    private subtractiveUpdate(r: number, k: number, c: number,
        A: CholeskyBlockMatrix): void {
        const blockSize = this.blockSize;
        const Arc = A[this.getIndex(r, c)];
        const Akc = A[this.getIndex(k, c)];
        const Ark = A[this.getIndex(r, k)];
        for (let j = 0; j < blockSize; ++j) {
            for (let i = 0; i < blockSize; ++i) {
                let value = Ark.get(j, i);
                for (let m = 0; m < blockSize; ++m) {
                    value -= Arc.get(j, m) * Akc.get(i, m);
                }
                Ark.set(j, i, value);
            }
        }
    }
}
