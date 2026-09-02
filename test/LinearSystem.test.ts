import { describe, it, expect } from 'vitest';
import { LinearSystem } from '../src/LinearSystem';
import type {
    LinearSystemSparseEntry, LinearSystemSparseMatrix
} from '../src/LinearSystem';
import { Matrix } from '../src/Matrix';
import { Vector } from '../src/Vector';

// A deterministic pseudorandom generator so the randomized cross-checks are
// reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state + 0x6D2B79F5) >>> 0;
        let t = state;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

// Compute A*X for a row-major NxN matrix A and an Nx1 vector X.
function mulRowMajor(n: number, A: readonly number[],
    X: readonly number[]): number[] {
    const P = new Array<number>(n).fill(0);
    for (let row = 0; row < n; ++row) {
        for (let col = 0; col < n; ++col) {
            P[row] += A[col + n * row] * X[col];
        }
    }
    return P;
}

// Build the full row-major matrix of a tridiagonal system.
function tridiagonalToDense(n: number, sub: readonly number[],
    diag: readonly number[], sup: readonly number[]): number[] {
    const A = new Array<number>(n * n).fill(0);
    for (let i = 0; i < n; ++i) {
        A[i + n * i] = diag[i];
    }
    for (let i = 0; i + 1 < n; ++i) {
        A[i + n * (i + 1)] = sub[i];
        A[(i + 1) + n * i] = sup[i];
    }
    return A;
}

// A symmetric positive definite matrix, row major.
function randomSPD(n: number, rng: () => number): number[] {
    const M = new Array<number>(n * n).fill(0);
    for (let i = 0; i < n * n; ++i) {
        M[i] = 2 * rng() - 1;
    }
    const A = new Array<number>(n * n).fill(0);
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            let sum = 0;
            for (let k = 0; k < n; ++k) {
                sum += M[k + n * r] * M[k + n * c];
            }
            A[c + n * r] = sum + (r === c ? n : 0);
        }
    }
    return A;
}

describe('LinearSystem fixed-size solvers', () => {
    it('solves a 2x2 system with a known solution', () => {
        // [2 1][x] = [5], solution (1, 3).
        // [1 3][y]   [10]
        const A = Matrix.fromArray(2, 2, [2, 1, 1, 3]);
        const B = Vector.fromArray([5, 10]);
        const { X, invertible } = LinearSystem.solve2x2(A, B);
        expect(invertible).toBe(true);
        expect(X.get(0)).toBeCloseTo(1, 12);
        expect(X.get(1)).toBeCloseTo(3, 12);
    });

    it('solves a 3x3 system with a known solution', () => {
        // The identity-like system with solution (1, 2, 3).
        const A = Matrix.fromArray(3, 3, [
            2, 0, 1,
            0, 3, -1,
            1, 1, 1
        ]);
        const expected = [1, 2, 3];
        const B = Vector.fromArray([
            2 * 1 + 0 * 2 + 1 * 3,
            0 * 1 + 3 * 2 - 1 * 3,
            1 * 1 + 1 * 2 + 1 * 3
        ]);
        const { X, invertible } = LinearSystem.solve3x3(A, B);
        expect(invertible).toBe(true);
        for (let i = 0; i < 3; ++i) {
            expect(X.get(i)).toBeCloseTo(expected[i], 12);
        }
    });

    it('solves a 4x4 system and reproduces the right-hand side', () => {
        const A = Matrix.fromArray(4, 4, [
            4, 1, 0, 2,
            1, 5, 2, 0,
            0, 2, 6, 1,
            2, 0, 1, 3
        ]);
        const expected = [1, -2, 3, 0.5];
        const B = new Vector(4);
        for (let r = 0; r < 4; ++r) {
            let sum = 0;
            for (let c = 0; c < 4; ++c) {
                sum += A.get(r, c) * expected[c];
            }
            B.set(r, sum);
        }
        const { X, invertible } = LinearSystem.solve4x4(A, B);
        expect(invertible).toBe(true);
        for (let i = 0; i < 4; ++i) {
            expect(X.get(i)).toBeCloseTo(expected[i], 10);
        }
    });

    it('reports singular fixed-size matrices and returns the zero vector',
        () => {
            const A2 = Matrix.fromArray(2, 2, [1, 2, 2, 4]);
            const r2 = LinearSystem.solve2x2(A2, Vector.fromArray([1, 1]));
            expect(r2.invertible).toBe(false);
            expect(r2.X.values).toEqual([0, 0]);

            const A3 = Matrix.fromArray(3, 3, [1, 2, 3, 2, 4, 6, 0, 1, 0]);
            const r3 = LinearSystem.solve3x3(A3, Vector.fromArray([1, 2, 3]));
            expect(r3.invertible).toBe(false);
            expect(r3.X.values).toEqual([0, 0, 0]);

            const A4 = Matrix.fromArray(4, 4, [
                1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 1, 0
            ]);
            const r4 = LinearSystem.solve4x4(A4,
                Vector.fromArray([1, 2, 3, 4]));
            expect(r4.invertible).toBe(false);
            expect(r4.X.values).toEqual([0, 0, 0, 0]);
        });

    it('agrees with the general solver on randomized fixed-size systems',
        () => {
            const rng = makeRng(20260901);
            for (let trial = 0; trial < 60; ++trial) {
                for (const n of [2, 3, 4]) {
                    const A = new Array<number>(n * n);
                    for (let i = 0; i < n * n; ++i) {
                        A[i] = 2 * rng() - 1;
                    }
                    // Bias the diagonal so the matrix is well conditioned.
                    for (let i = 0; i < n; ++i) {
                        A[i + n * i] += 3;
                    }
                    const B = new Array<number>(n);
                    for (let i = 0; i < n; ++i) {
                        B[i] = 2 * rng() - 1;
                    }

                    const general = LinearSystem.solve(n, A, B);
                    expect(general.invertible).toBe(true);

                    const M = Matrix.fromArray(n, n, A);
                    const V = Vector.fromArray(B);
                    const fixed = n === 2 ? LinearSystem.solve2x2(M, V)
                        : n === 3 ? LinearSystem.solve3x3(M, V)
                            : LinearSystem.solve4x4(M, V);
                    expect(fixed.invertible).toBe(true);
                    for (let i = 0; i < n; ++i) {
                        expect(fixed.X.get(i)).toBeCloseTo(general.X[i], 9);
                    }
                }
            }
        });
});

describe('LinearSystem general solvers', () => {
    it('solves A*X = B for a known 3x3 system', () => {
        const A = [
            2, 1, 1,
            1, 3, 2,
            1, 0, 0
        ];
        const B = [2 * 2 + 1 * 3, 1 * 2 + 3 * 3, 2];
        const { X, invertible } = LinearSystem.solve(3, A, B);
        expect(invertible).toBe(true);
        expect(X[0]).toBeCloseTo(2, 10);
        expect(X[1]).toBeCloseTo(3, 10);
        expect(X[2]).toBeCloseTo(0, 10);
    });

    it('reports a singular matrix', () => {
        const A = [1, 1, 2, 2];
        const { invertible } = LinearSystem.solve(2, A, [1, 2]);
        expect(invertible).toBe(false);
    });

    it('honors the column-major storage option', () => {
        // Row-major [1 2; 3 4] is column-major [1 3; 2 4].
        const A = [1, 2, 3, 4];
        const B = [5, 11];
        const rowMajor = LinearSystem.solve(2, A, B, true);
        const colMajor = LinearSystem.solve(2, A, B, false);
        // Row major: [1 2; 3 4] X = (5, 11) => X = (1, 2).
        expect(rowMajor.X[0]).toBeCloseTo(1, 10);
        expect(rowMajor.X[1]).toBeCloseTo(2, 10);
        // Column major: [1 3; 2 4] X = (5, 11) => X = (6.5, -0.5).
        expect(colMajor.X[0]).toBeCloseTo(6.5, 10);
        expect(colMajor.X[1]).toBeCloseTo(-0.5, 10);
    });

    it('solves systems with multiple right-hand sides', () => {
        const rng = makeRng(777);
        for (let trial = 0; trial < 20; ++trial) {
            const n = 5, m = 3;
            const A = new Array<number>(n * n);
            for (let i = 0; i < n * n; ++i) {
                A[i] = 2 * rng() - 1;
            }
            for (let i = 0; i < n; ++i) {
                A[i + n * i] += 4;
            }
            const B = new Array<number>(n * m);
            for (let i = 0; i < n * m; ++i) {
                B[i] = 2 * rng() - 1;
            }

            const { X, invertible } = LinearSystem.solveMultiple(n, m, A, B);
            expect(invertible).toBe(true);
            // Verify A*X = B column by column.
            for (let col = 0; col < m; ++col) {
                const xCol = new Array<number>(n);
                for (let r = 0; r < n; ++r) {
                    xCol[r] = X[col + m * r];
                }
                const p = mulRowMajor(n, A, xCol);
                for (let r = 0; r < n; ++r) {
                    expect(p[r]).toBeCloseTo(B[col + m * r], 9);
                }
            }
        }
    });
});

describe('LinearSystem tridiagonal solvers', () => {
    it('solves a known tridiagonal system', () => {
        // The 1D Laplacian [-1 2 -1] of size 4 applied to (1, 2, 3, 4).
        const n = 4;
        const sub = [-1, -1, -1];
        const diag = [2, 2, 2, 2];
        const sup = [-1, -1, -1];
        const expected = [1, 2, 3, 4];
        const dense = tridiagonalToDense(n, sub, diag, sup);
        const B = mulRowMajor(n, dense, expected);

        const { X, solved } = LinearSystem.solveTridiagonal(n, sub, diag, sup,
            B);
        expect(solved).toBe(true);
        for (let i = 0; i < n; ++i) {
            expect(X[i]).toBeCloseTo(expected[i], 10);
        }
    });

    it('agrees with the general solver on randomized tridiagonal systems',
        () => {
            const rng = makeRng(4242);
            for (let trial = 0; trial < 40; ++trial) {
                const n = 6;
                const sub = new Array<number>(n - 1);
                const sup = new Array<number>(n - 1);
                const diag = new Array<number>(n);
                for (let i = 0; i + 1 < n; ++i) {
                    sub[i] = 2 * rng() - 1;
                    sup[i] = 2 * rng() - 1;
                }
                for (let i = 0; i < n; ++i) {
                    // Diagonally dominant, so no zero pivots occur.
                    diag[i] = 4 + rng();
                }
                const B = new Array<number>(n);
                for (let i = 0; i < n; ++i) {
                    B[i] = 2 * rng() - 1;
                }

                const tri = LinearSystem.solveTridiagonal(n, sub, diag, sup,
                    B);
                expect(tri.solved).toBe(true);
                const dense = tridiagonalToDense(n, sub, diag, sup);
                const general = LinearSystem.solve(n, dense, B);
                expect(general.invertible).toBe(true);
                for (let i = 0; i < n; ++i) {
                    expect(tri.X[i]).toBeCloseTo(general.X[i], 10);
                }
            }
        });

    it('matches the constant-coefficient solver', () => {
        const n = 7;
        const a = -1, b = 3, c = -0.5;
        const sub = new Array<number>(n - 1).fill(a);
        const diag = new Array<number>(n).fill(b);
        const sup = new Array<number>(n - 1).fill(c);
        const B = [1, -2, 3, -4, 5, -6, 7];

        const varying = LinearSystem.solveTridiagonal(n, sub, diag, sup, B);
        const constant = LinearSystem.solveConstantTridiagonal(n, a, b, c, B);
        expect(varying.solved).toBe(true);
        expect(constant.solved).toBe(true);
        for (let i = 0; i < n; ++i) {
            expect(constant.X[i]).toBeCloseTo(varying.X[i], 12);
        }

        // Verify the solution against the dense system.
        const dense = tridiagonalToDense(n, sub, diag, sup);
        const p = mulRowMajor(n, dense, constant.X);
        for (let i = 0; i < n; ++i) {
            expect(p[i]).toBeCloseTo(B[i], 10);
        }
    });

    it('handles the degenerate 1x1 case and zero pivots', () => {
        const single = LinearSystem.solveTridiagonal(1, [], [4], [], [8]);
        expect(single.solved).toBe(true);
        expect(single.X).toEqual([2]);

        const singleConst = LinearSystem.solveConstantTridiagonal(1, 0, 4, 0,
            [8]);
        expect(singleConst.solved).toBe(true);
        expect(singleConst.X).toEqual([2]);

        // diagonal[0] = 0 fails immediately.
        const zeroFirst = LinearSystem.solveTridiagonal(3, [1, 1], [0, 2, 2],
            [1, 1], [1, 1, 1]);
        expect(zeroFirst.solved).toBe(false);
        expect(LinearSystem.solveConstantTridiagonal(3, 1, 0, 1,
            [1, 1, 1]).solved).toBe(false);

        // A zero pivot arises during the forward sweep: diagonal[1] equals
        // subdiagonal[0]*superdiagonal[0]/diagonal[0].
        const zeroPivot = LinearSystem.solveTridiagonal(3, [2, 1], [1, 4, 5],
            [2, 1], [1, 1, 1]);
        expect(zeroPivot.solved).toBe(false);

        expect(() => LinearSystem.solveTridiagonal(0, [], [], [], []))
            .toThrow();
        expect(() => LinearSystem.solveConstantTridiagonal(0, 1, 1, 1, []))
            .toThrow();
    });
});

describe('LinearSystem conjugate gradient solvers', () => {
    it('solves a symmetric positive definite system', () => {
        // [4 1; 1 3] X = (1, 2) has solution (1/11, 7/11).
        const A = [4, 1, 1, 3];
        const B = [1, 2];
        const { X, iterations } = LinearSystem.solveSymmetricCG(2, A, B, 32,
            1e-12);
        expect(X[0]).toBeCloseTo(1 / 11, 10);
        expect(X[1]).toBeCloseTo(7 / 11, 10);
        // Conjugate gradient is exact in at most N steps for an NxN system.
        expect(iterations).toBeLessThanOrEqual(3);
    });

    it('agrees with the general solver on randomized SPD systems', () => {
        const rng = makeRng(31337);
        for (let trial = 0; trial < 20; ++trial) {
            const n = 6;
            const A = randomSPD(n, rng);
            const B = new Array<number>(n);
            for (let i = 0; i < n; ++i) {
                B[i] = 2 * rng() - 1;
            }

            const cg = LinearSystem.solveSymmetricCG(n, A, B, 200, 1e-12);
            const general = LinearSystem.solve(n, A, B);
            expect(general.invertible).toBe(true);
            for (let i = 0; i < n; ++i) {
                expect(cg.X[i]).toBeCloseTo(general.X[i], 8);
            }

            // The residual A*X - B is small.
            const p = mulRowMajor(n, A, cg.X);
            for (let i = 0; i < n; ++i) {
                expect(p[i]).toBeCloseTo(B[i], 8);
            }
        }
    });

    it('gives the same answers for the sparse and dense representations',
        () => {
            const rng = makeRng(90210);
            for (let trial = 0; trial < 20; ++trial) {
                const n = 5;
                const A = randomSPD(n, rng);
                const B = new Array<number>(n);
                for (let i = 0; i < n; ++i) {
                    B[i] = 2 * rng() - 1;
                }

                // Store only the upper triangle (including the diagonal),
                // shuffled, to exercise the deterministic sort.
                const entries: LinearSystemSparseEntry[] = [];
                for (let r = n - 1; r >= 0; --r) {
                    for (let c = n - 1; c >= r; --c) {
                        entries.push({ row: r, col: c, value: A[c + n * r] });
                    }
                }

                const dense = LinearSystem.solveSymmetricCG(n, A, B, 200,
                    1e-12);
                const sparse: LinearSystemSparseMatrix = entries;
                const sparseResult = LinearSystem.solveSymmetricCGSparse(n,
                    sparse, B, 200, 1e-12);
                for (let i = 0; i < n; ++i) {
                    expect(sparseResult.X[i]).toBeCloseTo(dense.X[i], 9);
                }

                const p = mulRowMajor(n, A, sparseResult.X);
                for (let i = 0; i < n; ++i) {
                    expect(p[i]).toBeCloseTo(B[i], 8);
                }
            }
        });

    it('stops after maxIterations when the tolerance is unreachable', () => {
        const rng = makeRng(5);
        const n = 6;
        const A = randomSPD(n, rng);
        const B = new Array<number>(n);
        for (let i = 0; i < n; ++i) {
            B[i] = 2 * rng() - 1;
        }
        // A tolerance of 0 cannot be met before conjugate gradient has taken
        // its N exact steps, so with maxIterations = 2 the loop runs to
        // exhaustion and the reported iteration count is maxIterations + 1.
        const maxIterations = 2;
        const { iterations } = LinearSystem.solveSymmetricCG(n, A, B,
            maxIterations, 0);
        expect(iterations).toBe(maxIterations + 1);
    });
});
