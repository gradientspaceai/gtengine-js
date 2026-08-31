// gtengine-js: TypeScript port of Geometric Tools Engine (GTE)
// GaussianElimination.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

import { LexicoArray2 } from './LexicoArray2';
import { logError } from './Logger';

// Gaussian elimination with full pivoting. The input matrix M must be NxN.
// The storage convention for element lookup is chosen by the 'rowMajor'
// option (upstream selects it at compile time with GTE_USE_ROW_MAJOR or
// GTE_USE_COL_MAJOR; the port defaults to row major, the GTE default).
//
// If you want the inverse of M, set 'wantInverse'. If you want to solve
// M*X = B for X, where X and B are Nx1, pass B. If you want to solve M*Y = C
// for Y, where Y and C are NxK, pass C and pass K as 'numCols'. In all cases
// pass N as 'numRows'.
//
// Port notes:
// - Upstream signals "compute this output" by passing a non-null pointer and
//   writes the results through those pointers. The port takes an options
//   object and returns the outputs in a result object; 'inverseM', 'X' and
//   'Y' are null when the corresponding input was not requested.
// - Upstream's is_arbitrary_precision dispatch of the private Set() helper
//   exists only to choose memcpy/memset versus element-wise assignment. Both
//   branches have identical semantics, so the port keeps one element-wise
//   implementation.

export interface GaussianEliminationOptions {
    // Compute the inverse of M.
    wantInverse?: boolean;
    // Right-hand side of M*X = B, an Nx1 vector.
    B?: readonly number[] | null;
    // Right-hand side of M*Y = C, an NxK matrix in the 'rowMajor' order.
    C?: readonly number[] | null;
    // The number K of columns of C and Y. Required when C is provided.
    numCols?: number;
    // Storage order of M, inverseM, C and Y. Defaults to row major.
    rowMajor?: boolean;
}

export interface GaussianEliminationResult {
    // True iff M is invertible. When false, all outputs are zero-filled and
    // the determinant is zero.
    invertible: boolean;
    determinant: number;
    // The NxN inverse of M, or null when 'wantInverse' was not set.
    inverseM: number[] | null;
    // The Nx1 solution of M*X = B, or null when B was not provided.
    X: number[] | null;
    // The NxK solution of M*Y = C, or null when C was not provided.
    Y: number[] | null;
}

export class GaussianElimination {
    // Port of upstream's operator(). Returns the requested outputs; the
    // 'invertible' field is upstream's bool return value.
    compute(numRows: number, M: readonly number[],
        options: GaussianEliminationOptions = {}): GaussianEliminationResult {
        const wantInverse = options.wantInverse === true;
        const B = options.B ?? null;
        const C = options.C ?? null;
        const numCols = options.numCols ?? 0;
        const rowMajor = options.rowMajor ?? true;

        const numElements = numRows * numRows;
        if (numRows <= 0 || M.length < numElements || (C !== null && numCols < 1)) {
            logError('Invalid input.');
        }

        // Upstream allocates a scratch inverse when the caller does not want
        // one, because the elimination is performed in place on that buffer.
        const inverseM = new Array<number>(numElements).fill(0);
        for (let i = 0; i < numElements; ++i) {
            inverseM[i] = M[i];
        }

        let X: number[] | null = null;
        if (B !== null) {
            X = new Array<number>(numRows).fill(0);
            for (let i = 0; i < numRows; ++i) {
                X[i] = B[i];
            }
        }

        let Y: number[] | null = null;
        if (C !== null) {
            const numYElements = numRows * numCols;
            Y = new Array<number>(numYElements).fill(0);
            for (let i = 0; i < numYElements; ++i) {
                Y[i] = C[i];
            }
        }

        const matInvM = new LexicoArray2(rowMajor, numRows, numRows, inverseM);
        const matY = new LexicoArray2(rowMajor, numRows, numCols, Y ?? []);

        const colIndex = new Array<number>(numRows).fill(0);
        const rowIndex = new Array<number>(numRows).fill(0);
        const pivoted = new Array<boolean>(numRows).fill(false);

        let odd = false;
        let determinant = 1;

        // Elimination by full pivoting.
        let i1: number, i2: number, row = 0, col = 0;
        for (let i0 = 0; i0 < numRows; ++i0) {
            // Search the matrix (excluding pivoted rows) for the maximum
            // absolute entry.
            let maxValue = 0;
            for (i1 = 0; i1 < numRows; ++i1) {
                if (!pivoted[i1]) {
                    for (i2 = 0; i2 < numRows; ++i2) {
                        if (!pivoted[i2]) {
                            const value = matInvM.get(i1, i2);
                            const absValue = (value >= 0 ? value : -value);
                            if (absValue > maxValue) {
                                maxValue = absValue;
                                row = i1;
                                col = i2;
                            }
                        }
                    }
                }
            }

            if (maxValue === 0) {
                // The matrix is not invertible.
                return {
                    invertible: false,
                    determinant: 0,
                    inverseM: wantInverse ? new Array<number>(numElements).fill(0) : null,
                    X: B !== null ? new Array<number>(numRows).fill(0) : null,
                    Y: C !== null ? new Array<number>(numRows * numCols).fill(0) : null
                };
            }

            pivoted[col] = true;

            // Swap rows so that the pivot entry is in row 'col'.
            if (row !== col) {
                odd = !odd;
                for (let i = 0; i < numRows; ++i) {
                    const save = matInvM.get(row, i);
                    matInvM.set(row, i, matInvM.get(col, i));
                    matInvM.set(col, i, save);
                }

                if (X !== null) {
                    const save = X[row];
                    X[row] = X[col];
                    X[col] = save;
                }

                if (Y !== null) {
                    for (let i = 0; i < numCols; ++i) {
                        const save = matY.get(row, i);
                        matY.set(row, i, matY.get(col, i));
                        matY.set(col, i, save);
                    }
                }
            }

            // Keep track of the permutations of the rows.
            rowIndex[i0] = row;
            colIndex[i0] = col;

            // Scale the row so that the pivot entry is 1.
            const diagonal = matInvM.get(col, col);
            determinant *= diagonal;
            const inv = 1 / diagonal;
            matInvM.set(col, col, 1);
            for (i2 = 0; i2 < numRows; ++i2) {
                matInvM.set(col, i2, matInvM.get(col, i2) * inv);
            }

            if (X !== null) {
                X[col] *= inv;
            }

            if (Y !== null) {
                for (i2 = 0; i2 < numCols; ++i2) {
                    matY.set(col, i2, matY.get(col, i2) * inv);
                }
            }

            // Zero out the pivot column locations in the other rows.
            for (i1 = 0; i1 < numRows; ++i1) {
                if (i1 !== col) {
                    const save = matInvM.get(i1, col);
                    matInvM.set(i1, col, 0);
                    for (i2 = 0; i2 < numRows; ++i2) {
                        matInvM.set(i1, i2, matInvM.get(i1, i2) - matInvM.get(col, i2) * save);
                    }

                    if (X !== null) {
                        X[i1] -= X[col] * save;
                    }

                    if (Y !== null) {
                        for (i2 = 0; i2 < numCols; ++i2) {
                            matY.set(i1, i2, matY.get(i1, i2) - matY.get(col, i2) * save);
                        }
                    }
                }
            }
        }

        if (wantInverse) {
            // Reorder rows to undo any permutations in Gaussian elimination.
            for (i1 = numRows - 1; i1 >= 0; --i1) {
                if (rowIndex[i1] !== colIndex[i1]) {
                    for (i2 = 0; i2 < numRows; ++i2) {
                        const save = matInvM.get(i2, rowIndex[i1]);
                        matInvM.set(i2, rowIndex[i1], matInvM.get(i2, colIndex[i1]));
                        matInvM.set(i2, colIndex[i1], save);
                    }
                }
            }
        }

        if (odd) {
            determinant = -determinant;
        }

        return {
            invertible: true,
            determinant,
            inverseM: wantInverse ? inverseM : null,
            X,
            Y
        };
    }
}
