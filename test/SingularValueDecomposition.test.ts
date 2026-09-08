import { describe, it, expect } from 'vitest';
import { SingularValueDecomposition } from '../src/SingularValueDecomposition.js';
import { SymmetricEigensolver } from '../src/SymmetricEigensolver.js';
import {
    check, expectClose, fc, scaled
} from './helpers/arbitraries.js';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Multiply a (rows x inner) by an (inner x cols), both row-major.
function multiply(rows: number, inner: number, cols: number,
    A: readonly number[], B: readonly number[]): number[] {
    const C = new Array<number>(rows * cols).fill(0);
    for (let r = 0; r < rows; ++r) {
        for (let c = 0; c < cols; ++c) {
            let sum = 0;
            for (let k = 0; k < inner; ++k) {
                sum += A[k + inner * r] * B[c + cols * k];
            }
            C[c + cols * r] = sum;
        }
    }
    return C;
}

// The transpose of a (rows x cols) row-major matrix.
function transpose(rows: number, cols: number, A: readonly number[]): number[] {
    const T = new Array<number>(rows * cols).fill(0);
    for (let r = 0; r < rows; ++r) {
        for (let c = 0; c < cols; ++c) {
            T[r + rows * c] = A[c + cols * r];
        }
    }
    return T;
}

function maxDiff(a: readonly number[], b: readonly number[]): number {
    let m = 0;
    for (let i = 0; i < a.length; ++i) {
        m = Math.max(m, Math.abs(a[i] - b[i]));
    }
    return m;
}

// The deviation of Q^T*Q from the identity for an NxN row-major Q.
function orthogonalityError(n: number, q: readonly number[]): number {
    const product = multiply(n, n, n, transpose(n, n, q), q);
    let m = 0;
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            m = Math.max(m, Math.abs(product[c + n * r] - (r === c ? 1 : 0)));
        }
    }
    return m;
}

// A = U * S * V^T, which is the identity U^T*A*V = S rearranged.
function reconstruct(numRows: number, numCols: number, u: readonly number[],
    s: readonly number[], v: readonly number[]): number[] {
    const us = multiply(numRows, numRows, numCols, u, s);
    return multiply(numRows, numCols, numCols, us, transpose(numCols, numCols, v));
}

function randomMatrix(numRows: number, numCols: number,
    rand: () => number): number[] {
    const a = new Array<number>(numRows * numCols).fill(0);
    for (let i = 0; i < a.length; ++i) {
        a[i] = 2 * rand() - 1;
    }
    return a;
}

// The singular values of A computed independently as the square roots of the
// eigenvalues of A^T*A.
function referenceSingularValues(numRows: number, numCols: number,
    a: readonly number[]): number[] {
    const ata = multiply(numCols, numRows, numCols, transpose(numRows, numCols, a), a);
    const solver = new SymmetricEigensolver(numCols, 4096);
    solver.solve(ata, -1);
    return solver.getEigenvalues().map((value) => Math.sqrt(Math.max(value, 0)));
}

describe('SingularValueDecomposition', () => {
    it('rejects invalid constructor arguments', () => {
        expect(() => new SingularValueDecomposition(3, 1, 32)).toThrow('Invalid input.');
        expect(() => new SingularValueDecomposition(2, 3, 32)).toThrow('Invalid input.');
        expect(() => new SingularValueDecomposition(3, 3, 0)).toThrow('Invalid input.');
        expect(() => new SingularValueDecomposition(3, 3, 32)).not.toThrow();
    });

    it('rejects invalid inputs to solve', () => {
        const svd = new SingularValueDecomposition(2, 2, 32);
        expect(() => svd.solve([1, 0, 0, 1], 0)).toThrow('Invalid input to Solve.');
        expect(() => svd.solve([1, 0, 0], 8)).toThrow('Invalid input to Solve.');
    });

    it('decomposes a diagonal matrix into sorted nonnegative values', () => {
        // The singular values of diag(2,-5,3) are 5, 3 and 2.
        const n = 3;
        const a = [
            2, 0, 0,
            0, -5, 0,
            0, 0, 3
        ];
        const svd = new SingularValueDecomposition(n, n, 64);
        svd.solve(a);
        const values = svd.getSingularValues();
        expect(values[0]).toBeCloseTo(5, 13);
        expect(values[1]).toBeCloseTo(3, 13);
        expect(values[2]).toBeCloseTo(2, 13);

        const u = svd.getU();
        const v = svd.getV();
        const s = svd.getS();
        expect(orthogonalityError(n, u)).toBeLessThan(1e-13);
        expect(orthogonalityError(n, v)).toBeLessThan(1e-13);
        expect(maxDiff(reconstruct(n, n, u, s, v), a)).toBeLessThan(1e-12);
    });

    it('computes the known singular values of a 2x2 shear', () => {
        // A = {{1,1},{0,1}} has singular values (sqrt(5) + 1)/2 and
        // (sqrt(5) - 1)/2, whose product is |det(A)| = 1.
        const n = 2;
        const a = [1, 1, 0, 1];
        const svd = new SingularValueDecomposition(n, n, 64);
        const iterations = svd.solve(a);
        expect(iterations).toBeLessThan(64);

        const sqrt5 = Math.sqrt(5);
        expect(svd.getSingularValue(0)).toBeCloseTo(0.5 * (sqrt5 + 1), 13);
        expect(svd.getSingularValue(1)).toBeCloseTo(0.5 * (sqrt5 - 1), 13);
        expect(svd.getSingularValue(0) * svd.getSingularValue(1)).toBeCloseTo(1, 13);

        const u = svd.getU();
        const v = svd.getV();
        expect(maxDiff(reconstruct(n, n, u, svd.getS(), v), a)).toBeLessThan(1e-13);
    });

    it('computes the known decomposition of a 3x2 matrix', () => {
        // The last row is zero, so the decomposition is driven entirely by
        // the leading 2x2 block. The singular values are cross-checked
        // against the square roots of the eigenvalues of A^T*A.
        const numRows = 3;
        const numCols = 2;
        const a = [
            1, 3,
            2, 4,
            0, 0
        ];
        const svd = new SingularValueDecomposition(numRows, numCols, 64);
        svd.solve(a);
        const values = svd.getSingularValues();
        const reference = referenceSingularValues(numRows, numCols, a);
        expect(values[0]).toBeCloseTo(reference[0], 12);
        expect(values[1]).toBeCloseTo(reference[1], 12);
        expect(values[0]).toBeGreaterThanOrEqual(values[1]);

        const u = svd.getU();
        const v = svd.getV();
        const s = svd.getS();
        // S is MxN with the singular values on the diagonal and zeros
        // elsewhere, including the trailing M-N rows.
        for (let r = 0; r < numRows; ++r) {
            for (let c = 0; c < numCols; ++c) {
                if (r !== c) {
                    expect(s[c + numCols * r]).toBe(0);
                }
            }
        }
        expect(orthogonalityError(numRows, u)).toBeLessThan(1e-13);
        expect(orthogonalityError(numCols, v)).toBeLessThan(1e-13);
        expect(maxDiff(reconstruct(numRows, numCols, u, s, v), a)).toBeLessThan(1e-12);
    });

    it('handles the zero matrix', () => {
        const numRows = 4;
        const numCols = 3;
        const a = new Array<number>(numRows * numCols).fill(0);
        const svd = new SingularValueDecomposition(numRows, numCols, 64);
        expect(svd.solve(a)).toBe(0);
        expect(svd.getSingularValues()).toEqual([0, 0, 0]);
        expect(orthogonalityError(numRows, svd.getU())).toBeLessThan(1e-13);
        expect(orthogonalityError(numCols, svd.getV())).toBeLessThan(1e-13);
    });

    it('handles a rank-deficient matrix', () => {
        // Rows 2 and 3 are multiples of row 1, so the rank is 1 and the only
        // nonzero singular value is the Frobenius norm of A.
        const numRows = 3;
        const numCols = 3;
        const a = [
            1, 2, 3,
            2, 4, 6,
            -1, -2, -3
        ];
        const svd = new SingularValueDecomposition(numRows, numCols, 128);
        svd.solve(a);
        const values = svd.getSingularValues();
        let frobenius = 0;
        for (const element of a) {
            frobenius += element * element;
        }
        frobenius = Math.sqrt(frobenius);
        expect(values[0]).toBeCloseTo(frobenius, 12);
        expect(values[1]).toBeLessThan(1e-13);
        expect(values[2]).toBeLessThan(1e-13);

        const u = svd.getU();
        const v = svd.getV();
        expect(orthogonalityError(numRows, u)).toBeLessThan(1e-13);
        expect(orthogonalityError(numCols, v)).toBeLessThan(1e-13);
        expect(maxDiff(reconstruct(numRows, numCols, u, svd.getS(), v), a))
            .toBeLessThan(1e-12);
    });

    it('handles a matrix with an interior zero singular value', () => {
        // The bidiagonalization of this matrix has a zero diagonal entry in
        // the middle, which forces the decoupling path that chases the
        // neighboring superdiagonal entry out of that row.
        const n = 3;
        const a = [
            4, 1, 0,
            0, 0, 0,
            0, 0, 2
        ];
        const svd = new SingularValueDecomposition(n, n, 128);
        svd.solve(a);
        const values = svd.getSingularValues();
        const reference = referenceSingularValues(n, n, a);
        for (let i = 0; i < n; ++i) {
            expect(values[i]).toBeCloseTo(reference[i], 11);
        }
        expect(values[2]).toBeLessThan(1e-13);

        const u = svd.getU();
        const v = svd.getV();
        expect(orthogonalityError(n, u)).toBeLessThan(1e-13);
        expect(orthogonalityError(n, v)).toBeLessThan(1e-13);
        expect(maxDiff(reconstruct(n, n, u, svd.getS(), v), a)).toBeLessThan(1e-12);
    });

    it('handles input that is already bidiagonal (degenerate Householder steps)', () => {
        // Every Householder step operates on a subcolumn or subrow that is
        // already zero except for its leading entry, so the reflections are
        // degenerate. The U- and V-matrices must still be orthogonal and
        // reproduce A.
        const n = 5;
        const a = [
            2, 3, 0, 0, 0,
            0, -1, 4, 0, 0,
            0, 0, 5, 1, 0,
            0, 0, 0, 2, -6,
            0, 0, 0, 0, 3
        ];
        const svd = new SingularValueDecomposition(n, n, 512);
        svd.solve(a);

        const u = svd.getU();
        const v = svd.getV();
        const s = svd.getS();
        expect(orthogonalityError(n, u)).toBeLessThan(1e-12);
        expect(orthogonalityError(n, v)).toBeLessThan(1e-12);
        expect(maxDiff(reconstruct(n, n, u, s, v), a)).toBeLessThan(1e-11);

        const reference = referenceSingularValues(n, n, a);
        const values = svd.getSingularValues();
        for (let i = 0; i < n; ++i) {
            expect(values[i]).toBeCloseTo(reference[i], 10);
        }
    });

    it('handles repeated singular values', () => {
        // 2*I has the singular value 2 with multiplicity 4.
        const n = 4;
        const a = new Array<number>(n * n).fill(0);
        for (let d = 0; d < n; ++d) {
            a[d + n * d] = 2;
        }
        const svd = new SingularValueDecomposition(n, n, 64);
        svd.solve(a);
        for (const value of svd.getSingularValues()) {
            expect(value).toBeCloseTo(2, 13);
        }
        expect(maxDiff(reconstruct(n, n, svd.getU(), svd.getS(), svd.getV()), a))
            .toBeLessThan(1e-12);
    });

    it('exposes columns of U and V consistent with the full matrices', () => {
        const numRows = 5;
        const numCols = 3;
        const a = randomMatrix(numRows, numCols, makeRandom(4242));
        const svd = new SingularValueDecomposition(numRows, numCols, 256);
        svd.solve(a);

        const u = svd.getU();
        for (let index = 0; index < numRows; ++index) {
            const uColumn = svd.getUColumn(index);
            expect(uColumn.length).toBe(numRows);
            for (let row = 0; row < numRows; ++row) {
                expect(uColumn[row]).toBe(u[index + numRows * row]);
            }
        }

        const v = svd.getV();
        for (let index = 0; index < numCols; ++index) {
            const vColumn = svd.getVColumn(index);
            expect(vColumn.length).toBe(numCols);
            for (let row = 0; row < numCols; ++row) {
                expect(vColumn[row]).toBe(v[index + numCols * row]);
            }
        }

        const values = svd.getSingularValues();
        for (let index = 0; index < numCols; ++index) {
            expect(svd.getSingularValue(index)).toBe(values[index]);
        }

        // A*v_i = sigma_i * u_i for each i < N.
        for (let i = 0; i < numCols; ++i) {
            const vi = svd.getVColumn(i);
            const ui = svd.getUColumn(i);
            for (let row = 0; row < numRows; ++row) {
                let sum = 0;
                for (let c = 0; c < numCols; ++c) {
                    sum += a[c + numCols * row] * vi[c];
                }
                expect(sum).toBeCloseTo(values[i] * ui[row], 11);
            }
        }
    });

    it('returns copies of the internal matrices', () => {
        const svd = new SingularValueDecomposition(2, 2, 32);
        svd.solve([1, 2, 3, 4]);
        const u0 = svd.getU();
        u0[0] = 12345;
        expect(svd.getU()[0]).not.toBe(12345);
    });

    it('reconstructs random matrices of several shapes', () => {
        const rand = makeRandom(90210);
        const shapes: Array<[number, number]> = [
            [2, 2], [3, 2], [3, 3], [4, 2], [4, 3], [4, 4],
            [6, 3], [6, 5], [7, 7], [9, 4], [10, 8]
        ];
        for (const [numRows, numCols] of shapes) {
            for (let trial = 0; trial < 3; ++trial) {
                const a = randomMatrix(numRows, numCols, rand);
                const svd = new SingularValueDecomposition(numRows, numCols, 4096);
                const iterations = svd.solve(a);
                expect(iterations).not.toBe(SingularValueDecomposition.invalid);

                const u = svd.getU();
                const v = svd.getV();
                const s = svd.getS();
                expect(orthogonalityError(numRows, u)).toBeLessThan(1e-11);
                expect(orthogonalityError(numCols, v)).toBeLessThan(1e-11);
                expect(maxDiff(reconstruct(numRows, numCols, u, s, v), a))
                    .toBeLessThan(1e-11);

                const values = svd.getSingularValues();
                for (let i = 0; i < numCols; ++i) {
                    expect(values[i]).toBeGreaterThanOrEqual(0);
                    if (i > 0) {
                        expect(values[i - 1]).toBeGreaterThanOrEqual(values[i]);
                    }
                }

                // Cross-check against the eigenvalues of A^T*A.
                const reference = referenceSingularValues(numRows, numCols, a);
                for (let i = 0; i < numCols; ++i) {
                    expect(values[i]).toBeCloseTo(reference[i], 9);
                }
            }
        }
    });

    it('reproduces the sum of squares and the determinant identities', () => {
        // The squared Frobenius norm is the sum of squared singular values,
        // and for a square matrix |det(A)| is the product of the values.
        const n = 5;
        const a = randomMatrix(n, n, makeRandom(31337));
        const svd = new SingularValueDecomposition(n, n, 1024);
        svd.solve(a);
        const values = svd.getSingularValues();

        let frobeniusSqr = 0;
        for (const element of a) {
            frobeniusSqr += element * element;
        }
        const sumSqr = values.reduce((sum, value) => sum + value * value, 0);
        expect(sumSqr).toBeCloseTo(frobeniusSqr, 11);

        // A^T*A has eigenvalues equal to the squared singular values.
        const reference = referenceSingularValues(n, n, a);
        for (let i = 0; i < n; ++i) {
            expect(values[i]).toBeCloseTo(reference[i], 9);
        }
    });

    it('accepts a custom multiplier for the cutoffs', () => {
        const n = 4;
        const a = randomMatrix(n, n, makeRandom(2024));
        const tight = new SingularValueDecomposition(n, n, 1024);
        tight.solve(a, 1);
        const loose = new SingularValueDecomposition(n, n, 1024);
        loose.solve(a, 1024);
        // Both cutoffs are tiny multiples of the unit roundoff, so the
        // singular values agree to nearly full precision.
        expect(maxDiff(tight.getSingularValues(), loose.getSingularValues()))
            .toBeLessThan(1e-11);
    });

    it('reports failure to converge within the iteration budget', () => {
        const numRows = 8;
        const numCols = 6;
        const a = randomMatrix(numRows, numCols, makeRandom(555));
        const svd = new SingularValueDecomposition(numRows, numCols, 1);
        expect(svd.solve(a)).toBe(SingularValueDecomposition.invalid);
    });

    it('can be reused for several matrices', () => {
        const n = 3;
        const svd = new SingularValueDecomposition(n, n, 256);

        svd.solve([2, 0, 0, 0, 3, 0, 0, 0, 6]);
        let values = svd.getSingularValues();
        expect(values[0]).toBeCloseTo(6, 13);
        expect(values[2]).toBeCloseTo(2, 13);

        const a = randomMatrix(n, n, makeRandom(8080));
        svd.solve(a);
        expect(maxDiff(reconstruct(n, n, svd.getU(), svd.getS(), svd.getV()), a))
            .toBeLessThan(1e-12);

        svd.solve([2, 0, 0, 0, 3, 0, 0, 0, 6]);
        values = svd.getSingularValues();
        expect(values[0]).toBeCloseTo(6, 13);
        expect(values[1]).toBeCloseTo(3, 13);
        expect(values[2]).toBeCloseTo(2, 13);
    });
});

// ---------------------------------------------------------------------------
// Verification (V39): independent review against upstream
// SingularValueDecomposition.h.
// ---------------------------------------------------------------------------

describe('SingularValueDecomposition verification', () => {
    // A shape (M >= N >= 2) with a row-major M-by-N matrix of that shape.
    const shaped = fc.tuple(fc.integer({ min: 2, max: 5 }),
        fc.integer({ min: 0, max: 3 })).chain(([numCols, extra]) => {
            const numRows = numCols + extra;
            return fc.tuple(fc.constant(numRows), fc.constant(numCols),
                fc.array(scaled(-5, 5),
                    { minLength: numRows * numCols,
                        maxLength: numRows * numCols }));
        });

    function multiply(numRows: number, numCommon: number, numCols: number,
        A: readonly number[], B: readonly number[]): number[] {
        const C = new Array<number>(numRows * numCols).fill(0);
        for (let r = 0; r < numRows; ++r) {
            for (let c = 0; c < numCols; ++c) {
                let sum = 0;
                for (let k = 0; k < numCommon; ++k) {
                    sum += A[k + numCommon * r] * B[c + numCols * k];
                }
                C[c + numCols * r] = sum;
            }
        }
        return C;
    }

    function transpose(numRows: number, numCols: number,
        A: readonly number[]): number[] {
        const T = new Array<number>(numRows * numCols).fill(0);
        for (let r = 0; r < numRows; ++r) {
            for (let c = 0; c < numCols; ++c) {
                T[r + numRows * c] = A[c + numCols * r];
            }
        }
        return T;
    }

    function lInfinity(A: readonly number[]): number {
        let worst = 0;
        for (const value of A) { worst = Math.max(worst, Math.abs(value)); }
        return worst;
    }

    // ||M^T*M - I||_inf for an n-by-n matrix M.
    function orthonormalityError(n: number, M: readonly number[]): number {
        let worst = 0;
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                let sum = 0;
                for (let k = 0; k < n; ++k) {
                    sum += M[r + n * k] * M[c + n * k];
                }
                worst = Math.max(worst, Math.abs(sum - (r === c ? 1 : 0)));
            }
        }
        return worst;
    }

    it('reconstructs A = U*S*V^T with orthogonal U and V', () => {
        check(shaped, ([numRows, numCols, A]) => {
            const svd = new SingularValueDecomposition(numRows, numCols, 64);
            const iterations = svd.solve(A);
            expect(iterations).not
                .toBe(SingularValueDecomposition.invalid);

            const U = svd.getU(), S = svd.getS(), V = svd.getV();
            const US = multiply(numRows, numRows, numCols, U, S);
            const VT = transpose(numCols, numCols, V);
            const USVT = multiply(numRows, numCols, numCols, US, VT);
            const scale = Math.max(lInfinity(A), 1);
            for (let i = 0; i < numRows * numCols; ++i) {
                // Golub and Van Loan expect |U^T*A*V - S| to be about
                // unitRoundoff*|A|; the bound below is that with room for the
                // O(n) accumulation of the reconstruction products.
                expectClose(USVT[i], A[i], 1e-11 * scale, 0);
            }

            // Both factors are products of Householder reflections and
            // Givens rotations, so they are orthogonal to the same order.
            expect(orthonormalityError(numRows, U)).toBeLessThan(1e-12);
            expect(orthonormalityError(numCols, V)).toBeLessThan(1e-12);
        }, 100);
    });

    it('reports nonnegative singular values in descending order, off the '
        + 'diagonal zero', () => {
        check(shaped, ([numRows, numCols, A]) => {
            const svd = new SingularValueDecomposition(numRows, numCols, 64);
            expect(svd.solve(A)).not.toBe(SingularValueDecomposition.invalid);

            const S = svd.getS();
            const values = svd.getSingularValues();
            expect(values.length).toBe(numCols);
            for (let i = 0; i < numCols; ++i) {
                expect(values[i]).toBeGreaterThanOrEqual(0);
                if (i > 0) {
                    expect(values[i]).toBeLessThanOrEqual(values[i - 1]);
                }
                expect(svd.getSingularValue(i)).toBe(values[i]);
            }
            for (let r = 0; r < numRows; ++r) {
                for (let c = 0; c < numCols; ++c) {
                    if (r !== c) {
                        expect(S[c + numCols * r]).toBe(0);
                    } else {
                        expect(S[c + numCols * r]).toBe(values[r]);
                    }
                }
            }
        }, 100);
    });

    it('agrees with the eigenvalues of A^T*A', () => {
        check(shaped, ([numRows, numCols, A]) => {
            // The singular values are the square roots of the eigenvalues of
            // A^T*A. SymmetricEigensolver is a different algorithm
            // (tridiagonalization plus QR iteration), so this is an
            // independent computation.
            const ATA = new Array<number>(numCols * numCols).fill(0);
            for (let r = 0; r < numCols; ++r) {
                for (let c = 0; c < numCols; ++c) {
                    let sum = 0;
                    for (let k = 0; k < numRows; ++k) {
                        sum += A[r + numCols * k] * A[c + numCols * k];
                    }
                    ATA[c + numCols * r] = sum;
                }
            }
            const eigen = new SymmetricEigensolver(numCols, 4096);
            expect(eigen.solve(ATA, -1)).toBeGreaterThanOrEqual(0);
            const eigenvalues = eigen.getEigenvalues();

            const svd = new SingularValueDecomposition(numRows, numCols, 64);
            expect(svd.solve(A)).not.toBe(SingularValueDecomposition.invalid);
            const values = svd.getSingularValues();

            const scale = Math.max(lInfinity(A), 1);
            for (let i = 0; i < numCols; ++i) {
                // The eigenvalues come back in decreasing order (sortType
                // -1), matching the singular-value order. Forming A^T*A
                // squares the condition number, so a small singular value is
                // resolved only to about sqrt(eps)*|A|.
                const expected = Math.sqrt(Math.max(eigenvalues[i], 0));
                expectClose(values[i], expected, 1e-7 * scale, 1e-7);
            }
        }, 100);
    });

    it('places the exact singular values of a scaled orthogonal-like input',
        () => {
            check(fc.tuple(fc.integer({ min: 2, max: 5 }),
                fc.integer({ min: 0, max: 3 })).chain(([numCols, extra]) => {
                    const numRows = numCols + extra;
                    return fc.tuple(fc.constant(numRows),
                        fc.constant(numCols),
                        fc.array(scaled(0.25, 8),
                            { minLength: numCols, maxLength: numCols }));
                }), ([numRows, numCols, d]) => {
                // A is M-by-N with the diagonal entries d and zeros
                // elsewhere, so its singular values are the |d| sorted in
                // descending order.
                const A = new Array<number>(numRows * numCols).fill(0);
                for (let i = 0; i < numCols; ++i) {
                    A[i + numCols * i] = d[i];
                }
                const svd = new SingularValueDecomposition(numRows, numCols,
                    64);
                expect(svd.solve(A)).not
                    .toBe(SingularValueDecomposition.invalid);
                const expected = d.map(Math.abs).sort((p, q) => q - p);
                const values = svd.getSingularValues();
                for (let i = 0; i < numCols; ++i) {
                    expectClose(values[i], expected[i], 1e-12, 1e-12);
                }
            }, 100);
        });

    it('scales its singular values with the matrix', () => {
        check(fc.tuple(shaped, scaled(0.05, 20), fc.boolean()),
            ([[numRows, numCols, A], factor, negate]) => {
                const scale = negate ? -factor : factor;
                const scaled2 = A.map(value => scale * value);
                const base = new SingularValueDecomposition(numRows, numCols,
                    64);
                const other = new SingularValueDecomposition(numRows, numCols,
                    64);
                expect(base.solve(A)).not
                    .toBe(SingularValueDecomposition.invalid);
                expect(other.solve(scaled2)).not
                    .toBe(SingularValueDecomposition.invalid);
                const v0 = base.getSingularValues();
                const v1 = other.getSingularValues();
                for (let i = 0; i < numCols; ++i) {
                    // The algorithm is not scale invariant bit for bit (the
                    // cutoffs are relative to |B|, but the arithmetic
                    // differs), so compare with a relative tolerance.
                    expectClose(v1[i], Math.abs(scale) * v0[i], 1e-9, 1e-9);
                }
            }, 100);
    });

    it('reports a rank deficiency as a zero singular value', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 4 }),
            fc.integer({ min: 0, max: 2 })).chain(([numCols, extra]) => {
                const numRows = numCols + extra;
                return fc.tuple(fc.constant(numRows), fc.constant(numCols),
                    fc.array(scaled(-5, 5),
                        { minLength: numRows * numCols,
                            maxLength: numRows * numCols }),
                    scaled(-3, 3));
            }), ([numRows, numCols, A, factor]) => {
            // Replace the last column by 'factor' times the first, which
            // makes the matrix rank deficient.
            const M = A.slice();
            for (let r = 0; r < numRows; ++r) {
                M[(numCols - 1) + numCols * r] = factor * M[numCols * r];
            }
            const svd = new SingularValueDecomposition(numRows, numCols, 64);
            expect(svd.solve(M)).not
                .toBe(SingularValueDecomposition.invalid);
            const values = svd.getSingularValues();
            const scale = Math.max(lInfinity(M), 1);
            // The smallest singular value is zero up to the round-off of the
            // decomposition.
            expect(values[numCols - 1]).toBeLessThan(1e-11 * scale);
        }, 100);
    });

    it('exposes columns of U and V that match the full matrices', () => {
        check(shaped, ([numRows, numCols, A]) => {
            const svd = new SingularValueDecomposition(numRows, numCols, 64);
            expect(svd.solve(A)).not.toBe(SingularValueDecomposition.invalid);
            const U = svd.getU(), V = svd.getV();
            for (let index = 0; index < numRows; ++index) {
                const column = svd.getUColumn(index);
                for (let row = 0; row < numRows; ++row) {
                    expect(column[row]).toBe(U[index + numRows * row]);
                }
            }
            for (let index = 0; index < numCols; ++index) {
                const column = svd.getVColumn(index);
                for (let row = 0; row < numCols; ++row) {
                    expect(column[row]).toBe(V[index + numCols * row]);
                }
            }
        }, 60);
    });

    it('rejects the shapes upstream excludes and out-of-range indices',
        () => {
            // Upstream requires N >= 2 and M >= N, so the 1xN and Nx1 cases
            // are not supported at all.
            expect(() => new SingularValueDecomposition(1, 1, 16))
                .toThrow('Invalid input.');
            expect(() => new SingularValueDecomposition(4, 1, 16))
                .toThrow('Invalid input.');
            expect(() => new SingularValueDecomposition(2, 4, 16))
                .toThrow('Invalid input.');
            expect(() => new SingularValueDecomposition(4, 2, 0))
                .toThrow('Invalid input.');

            const svd = new SingularValueDecomposition(3, 2, 32);
            svd.solve([1, 2, 3, 4, 5, 6]);
            expect(() => svd.getSingularValue(2)).toThrow('Invalid index');
            // Upstream's index is a size_t, so a negative value wraps and
            // trips the same assertion.
            expect(() => svd.getSingularValue(-1)).toThrow('Invalid index');
            expect(() => svd.getUColumn(-1)).toThrow('Invalid index');
            expect(() => svd.getUColumn(3)).toThrow('Invalid index');
            expect(() => svd.getVColumn(-1)).toThrow('Invalid index');
            expect(() => svd.getVColumn(2)).toThrow('Invalid index');
        });
});
