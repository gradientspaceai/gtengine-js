// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) LinearSystem.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Solve linear systems of equations where the matrix A is NxN. The
// 'invertible' member of the returned object is true when A is invertible;
// in that case the solution X is valid. When it is false, A is not
// invertible and X is invalid, so do not use it. When a matrix is passed as
// a flat array of numbers, the storage order is selected by the 'rowMajor'
// argument (upstream selects it at compile time with GTE_USE_ROW_MAJOR or
// GTE_USE_COL_MAJOR; the port defaults to row major, the GTE default, as
// GaussianElimination.ts does).
//
// The linear solvers that use the conjugate gradient algorithm are based on
// the discussion in "Matrix Computations, 2nd edition" by G. H. Golub and
// Charles F. Van Loan, The Johns Hopkins Press, Baltimore MD, Fourth
// Printing 1993.
//
// Port notes:
// - Upstream's overloaded static Solve members become distinctly named
//   statics, following the fixed-size naming precedent: solve2x2, solve3x3
//   and solve4x4 for the fixed-size matrix overloads, solve for the general
//   Nx1 right-hand side and solveMultiple for the NxM right-hand side.
//   Matrix2x2/3x3/4x4 are the port's doc-only aliases of Matrix, so runtime
//   size assertions replace the C++ compile-time sizes.
// - Output-parameter results become returned objects: { X, invertible } for
//   the invertible-matrix solvers, { X, solved } for the tridiagonal solvers
//   (whose false return means a zero pivot was encountered, not that a
//   determinant was computed), and { X, iterations } for the conjugate
//   gradient solvers.
// - The SparseMatrix std::map keyed by (i,j) becomes an array of
//   { row, col, value } entries. The port sorts a copy of the entries by
//   (row, col) once per solve so the accumulation order - and therefore the
//   floating-point result - matches upstream's std::map iteration order.
// - Upstream's SolveTridiagonal/SolveConstantTridiagonal index B[0] and
//   allocate a std::vector of size N-1 without validating N; N = 0 makes the
//   size computation underflow (size_t) and reads B[0] out of bounds. The
//   port asserts N >= 1 instead.
// - In SolveSymmetricCG upstream recomputes the loop invariant Dot(N, B, B)
//   on every iteration; the port hoists it out of the loop, which is
//   value-preserving.

import { GaussianElimination } from './GaussianElimination';
import { LexicoArray2 } from './LexicoArray2';
import { logAssert } from './Logger';
import { Matrix, mulMatrix } from './Matrix';
import { inverse2x2 } from './Matrix2x2';
import { inverse3x3 } from './Matrix3x3';
import { inverse4x4 } from './Matrix4x4';
import { Vector } from './Vector';

// The solution of A*X = B for a fixed-size system.
export interface LinearSystemVectorResult {
    X: Vector;
    invertible: boolean;
}

// The solution of A*X = B, where X is stored as a flat array.
export interface LinearSystemArrayResult {
    X: number[];
    invertible: boolean;
}

// The solution of a tridiagonal system. The 'solved' member is false when a
// zero pivot was encountered, in which case X is invalid.
export interface LinearSystemTridiagonalResult {
    X: number[];
    solved: boolean;
}

// The result of a conjugate gradient solve. The 'iterations' member is
// upstream's return value: the number of the iteration at which the
// residual test succeeded, or maxIterations + 1 when it never did.
export interface LinearSystemCGResult {
    X: number[];
    iterations: number;
}

// One nonzero entry of a sparse symmetric matrix. Only one of (i,j) and
// (j,i) should be stored, since the matrix is symmetric.
export interface LinearSystemSparseEntry {
    row: number;
    col: number;
    value: number;
}

export type LinearSystemSparseMatrix = readonly LinearSystemSparseEntry[];

export class LinearSystem {
    // Solve 2x2, 3x3 and 4x4 systems by inverting the matrix directly. This
    // avoids the overhead of Gaussian elimination in small dimensions.
    static solve2x2(A: Matrix, B: Vector): LinearSystemVectorResult {
        logAssert(B.size === 2, 'Invalid size.');
        const { inverse, invertible } = inverse2x2(A);
        return {
            X: invertible ? mulMatrix(inverse, B) : Vector.zero(2),
            invertible
        };
    }

    static solve3x3(A: Matrix, B: Vector): LinearSystemVectorResult {
        logAssert(B.size === 3, 'Invalid size.');
        const { inverse, invertible } = inverse3x3(A);
        return {
            X: invertible ? mulMatrix(inverse, B) : Vector.zero(3),
            invertible
        };
    }

    static solve4x4(A: Matrix, B: Vector): LinearSystemVectorResult {
        logAssert(B.size === 4, 'Invalid size.');
        const { inverse, invertible } = inverse4x4(A);
        return {
            X: invertible ? mulMatrix(inverse, B) : Vector.zero(4),
            invertible
        };
    }

    // Solve A*X = B, where A is NxN, B is Nx1 and the solution X is Nx1.
    static solve(n: number, A: readonly number[], B: readonly number[],
        rowMajor: boolean = true): LinearSystemArrayResult {
        const result = new GaussianElimination().compute(n, A, { B, rowMajor });
        return {
            X: result.X as number[],
            invertible: result.invertible
        };
    }

    // Solve A*X = B, where A is NxN, B is NxM and the solution X is NxM.
    static solveMultiple(n: number, m: number, A: readonly number[],
        B: readonly number[], rowMajor: boolean = true):
        LinearSystemArrayResult {
        const result = new GaussianElimination().compute(n, A,
            { C: B, numCols: m, rowMajor });
        return {
            X: result.Y as number[],
            invertible: result.invertible
        };
    }

    // Solve A*X = B, where A is tridiagonal. The function expects the
    // subdiagonal, diagonal and superdiagonal of A. The diagonal input must
    // have N elements. The subdiagonal and superdiagonal inputs must have
    // N-1 elements.
    static solveTridiagonal(n: number, subdiagonal: readonly number[],
        diagonal: readonly number[], superdiagonal: readonly number[],
        B: readonly number[]): LinearSystemTridiagonalResult {
        logAssert(n >= 1, 'Invalid size.');
        const X = new Array<number>(n).fill(0);
        if (diagonal[0] === 0) {
            return { X, solved: false };
        }

        const tmp = new Array<number>(n - 1).fill(0);
        let expr = diagonal[0];
        let invExpr = 1 / expr;
        X[0] = B[0] * invExpr;

        let i0: number, i1: number;
        for (i0 = 0, i1 = 1; i1 < n; ++i0, ++i1) {
            tmp[i0] = superdiagonal[i0] * invExpr;
            expr = diagonal[i1] - subdiagonal[i0] * tmp[i0];
            if (expr === 0) {
                return { X, solved: false };
            }
            invExpr = 1 / expr;
            X[i1] = (B[i1] - subdiagonal[i0] * X[i0]) * invExpr;
        }

        for (i0 = n - 1, i1 = n - 2; i1 >= 0; --i0, --i1) {
            X[i1] -= tmp[i1] * X[i0];
        }
        return { X, solved: true };
    }

    // Solve A*X = B, where A is tridiagonal, the subdiagonal elements are a
    // constant, the diagonal elements are a constant and the superdiagonal
    // elements are a constant.
    static solveConstantTridiagonal(n: number, subdiagonal: number,
        diagonal: number, superdiagonal: number, B: readonly number[]):
        LinearSystemTridiagonalResult {
        logAssert(n >= 1, 'Invalid size.');
        const X = new Array<number>(n).fill(0);
        if (diagonal === 0) {
            return { X, solved: false };
        }

        const tmp = new Array<number>(n - 1).fill(0);
        let expr = diagonal;
        let invExpr = 1 / expr;
        X[0] = B[0] * invExpr;

        let i0: number, i1: number;
        for (i0 = 0, i1 = 1; i1 < n; ++i0, ++i1) {
            tmp[i0] = superdiagonal * invExpr;
            expr = diagonal - subdiagonal * tmp[i0];
            if (expr === 0) {
                return { X, solved: false };
            }
            invExpr = 1 / expr;
            X[i1] = (B[i1] - subdiagonal * X[i0]) * invExpr;
        }

        for (i0 = n - 1, i1 = n - 2; i1 >= 0; --i0, --i1) {
            X[i1] -= tmp[i1] * X[i0];
        }
        return { X, solved: true };
    }

    // Solve A*X = B using the conjugate gradient method, where A is NxN and
    // symmetric. You must specify the maximum number of iterations and a
    // tolerance for terminating the iterations. A reasonable choice for the
    // tolerance is 1e-08 for double-precision arithmetic.
    static solveSymmetricCG(n: number, A: readonly number[],
        B: readonly number[], maxIterations: number, tolerance: number,
        rowMajor: boolean = true): LinearSystemCGResult {
        const matA = new LexicoArray2(rowMajor, n, n, A.slice());
        return solveCG(n, (X, P) => mulDense(n, matA, X, P), B, maxIterations,
            tolerance);
    }

    // Solve A*X = B using the conjugate gradient method, where A is sparse
    // and symmetric. Only one of the entries (i,j) and (j,i) should be
    // provided, since A is symmetric.
    static solveSymmetricCGSparse(n: number, A: LinearSystemSparseMatrix,
        B: readonly number[], maxIterations: number, tolerance: number):
        LinearSystemCGResult {
        // Replicate the std::map iteration order of upstream, so the
        // accumulation order of the products is the same.
        const entries = A.slice().sort(
            (e0, e1) => (e0.row !== e1.row ? e0.row - e1.row
                : e0.col - e1.col));
        return solveCG(n, (X, P) => mulSparse(n, entries, X, P), B,
            maxIterations, tolerance);
    }
}

// Support for the conjugate gradient method.
type MulFunction = (X: readonly number[], P: number[]) => void;

function dotArray(n: number, U: readonly number[],
    V: readonly number[]): number {
    let sum = 0;
    for (let i = 0; i < n; ++i) {
        sum += U[i] * V[i];
    }
    return sum;
}

function mulDense(n: number, matA: LexicoArray2, X: readonly number[],
    P: number[]): void {
    P.fill(0);
    for (let row = 0; row < n; ++row) {
        for (let col = 0; col < n; ++col) {
            P[row] += matA.get(row, col) * X[col];
        }
    }
}

function mulSparse(n: number, A: LinearSystemSparseMatrix,
    X: readonly number[], P: number[]): void {
    P.fill(0);
    for (const element of A) {
        const i = element.row, j = element.col, value = element.value;
        P[i] += value * X[j];
        if (i !== j) {
            P[j] += value * X[i];
        }
    }
}

function updateX(n: number, X: number[], alpha: number,
    P: readonly number[]): void {
    for (let i = 0; i < n; ++i) {
        X[i] += alpha * P[i];
    }
}

function updateR(n: number, R: number[], alpha: number,
    W: readonly number[]): void {
    for (let i = 0; i < n; ++i) {
        R[i] -= alpha * W[i];
    }
}

function updateP(n: number, P: number[], beta: number,
    R: readonly number[]): void {
    for (let i = 0; i < n; ++i) {
        P[i] = R[i] + beta * P[i];
    }
}

// The body shared by the dense and sparse conjugate gradient solvers; the
// only difference upstream is which Mul overload is called.
function solveCG(n: number, mul: MulFunction, B: readonly number[],
    maxIterations: number, tolerance: number): LinearSystemCGResult {
    // The first iteration.
    const X = new Array<number>(n).fill(0);
    const R = new Array<number>(n).fill(0);
    const P = new Array<number>(n).fill(0);
    const W = new Array<number>(n).fill(0);
    for (let i = 0; i < n; ++i) {
        R[i] = B[i];
    }
    let rho0 = dotArray(n, R, R);
    for (let i = 0; i < n; ++i) {
        P[i] = R[i];
    }
    mul(P, W);
    let alpha = rho0 / dotArray(n, P, W);
    updateX(n, X, alpha, P);
    updateR(n, R, alpha, W);
    let rho1 = dotArray(n, R, R);

    // The right-hand side does not change, so its length is a loop
    // invariant (upstream recomputes it on every iteration).
    const root1 = Math.sqrt(dotArray(n, B, B));

    // The remaining iterations.
    let iteration: number;
    for (iteration = 1; iteration <= maxIterations; ++iteration) {
        const root0 = Math.sqrt(rho1);
        if (root0 <= tolerance * root1) {
            break;
        }

        const beta = rho1 / rho0;
        updateP(n, P, beta, R);
        mul(P, W);
        alpha = rho1 / dotArray(n, P, W);
        updateX(n, X, alpha, P);
        updateR(n, R, alpha, W);
        rho0 = rho1;
        rho1 = dotArray(n, R, R);
    }
    return { X, iterations: iteration };
}
