import { describe, it, expect } from 'vitest';
import { GaussianElimination } from '../src/GaussianElimination.js';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Row-major matrix multiply of an NxN by an NxK.
function multiply(n: number, k: number, A: readonly number[],
    B: readonly number[]): number[] {
    const C = new Array<number>(n * k).fill(0);
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < k; ++c) {
            let sum = 0;
            for (let i = 0; i < n; ++i) {
                sum += A[i + n * r] * B[c + k * i];
            }
            C[c + k * r] = sum;
        }
    }
    return C;
}

// Determinant by cofactor expansion, an independent reference computation.
function cofactorDeterminant(n: number, A: readonly number[]): number {
    if (n === 1) {
        return A[0];
    }
    let det = 0;
    let sign = 1;
    for (let c = 0; c < n; ++c) {
        const minor: number[] = [];
        for (let r = 1; r < n; ++r) {
            for (let cc = 0; cc < n; ++cc) {
                if (cc !== c) {
                    minor.push(A[cc + n * r]);
                }
            }
        }
        det += sign * A[c] * cofactorDeterminant(n - 1, minor);
        sign = -sign;
    }
    return det;
}

function transpose(n: number, A: readonly number[]): number[] {
    const T = new Array<number>(n * n).fill(0);
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            T[r + n * c] = A[c + n * r];
        }
    }
    return T;
}

describe('GaussianElimination', () => {
    const ge = new GaussianElimination();

    it('inverts a 1x1 matrix', () => {
        const result = ge.compute(1, [4], { wantInverse: true });
        expect(result.invertible).toBe(true);
        expect(result.determinant).toBe(4);
        expect(result.inverseM).toEqual([0.25]);
    });

    it('inverts a 2x2 matrix and computes its determinant', () => {
        // M = [[4, 7], [2, 6]], det = 10, inverse = [[0.6, -0.7], [-0.2, 0.4]].
        const result = ge.compute(2, [4, 7, 2, 6], { wantInverse: true });
        expect(result.invertible).toBe(true);
        expect(result.determinant).toBeCloseTo(10, 12);
        const inv = result.inverseM as number[];
        expect(inv[0]).toBeCloseTo(0.6, 12);
        expect(inv[1]).toBeCloseTo(-0.7, 12);
        expect(inv[2]).toBeCloseTo(-0.2, 12);
        expect(inv[3]).toBeCloseTo(0.4, 12);
    });

    it('computes the determinant of a permutation matrix (odd number of swaps)', () => {
        // The 2x2 exchange matrix has determinant -1.
        const swap2 = ge.compute(2, [0, 1, 1, 0], {});
        expect(swap2.invertible).toBe(true);
        expect(swap2.determinant).toBeCloseTo(-1, 12);

        // The 3x3 cyclic permutation matrix has determinant +1.
        const cyc3 = ge.compute(3, [0, 1, 0, 0, 0, 1, 1, 0, 0], {});
        expect(cyc3.invertible).toBe(true);
        expect(cyc3.determinant).toBeCloseTo(1, 12);
    });

    it('solves M*X = B', () => {
        // 2x + y = 5, x + 3y = 10 -> x = 1, y = 3.
        const result = ge.compute(2, [2, 1, 1, 3], { B: [5, 10] });
        expect(result.invertible).toBe(true);
        expect(result.inverseM).toBeNull();
        const X = result.X as number[];
        expect(X[0]).toBeCloseTo(1, 12);
        expect(X[1]).toBeCloseTo(3, 12);
    });

    it('solves M*Y = C for a matrix right-hand side', () => {
        const M = [2, 1, 1, 3];
        // C is 2x2 row major: [[5, 2], [10, 1]].
        const C = [5, 2, 10, 1];
        const result = ge.compute(2, M, { C, numCols: 2 });
        expect(result.invertible).toBe(true);
        const Y = result.Y as number[];
        const MY = multiply(2, 2, M, Y);
        for (let i = 0; i < 4; ++i) {
            expect(MY[i]).toBeCloseTo(C[i], 12);
        }
    });

    it('computes the inverse, the solution and the determinant together', () => {
        const M = [2, 0, 1, 1, 3, 2, 0, 1, 4];
        const B = [3, 6, 5];
        const C = [1, 0, 0, 1, 1, 1];
        const result = ge.compute(3, M,
            { wantInverse: true, B, C, numCols: 2 });
        expect(result.invertible).toBe(true);
        expect(result.determinant).toBeCloseTo(cofactorDeterminant(3, M), 10);

        const inv = result.inverseM as number[];
        const identity = multiply(3, 3, M, inv);
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 3; ++c) {
                expect(identity[c + 3 * r]).toBeCloseTo(r === c ? 1 : 0, 10);
            }
        }

        const MX = multiply(3, 1, M, result.X as number[]);
        for (let i = 0; i < 3; ++i) {
            expect(MX[i]).toBeCloseTo(B[i], 10);
        }

        const MY = multiply(3, 2, M, result.Y as number[]);
        for (let i = 0; i < 6; ++i) {
            expect(MY[i]).toBeCloseTo(C[i], 10);
        }
    });

    it('reports a singular matrix and zero-fills the outputs', () => {
        // The second row is twice the first.
        const result = ge.compute(2, [1, 2, 2, 4],
            { wantInverse: true, B: [1, 2], C: [1, 0, 0, 1], numCols: 2 });
        expect(result.invertible).toBe(false);
        expect(result.determinant).toBe(0);
        expect(result.inverseM).toEqual([0, 0, 0, 0]);
        expect(result.X).toEqual([0, 0]);
        expect(result.Y).toEqual([0, 0, 0, 0]);
    });

    it('reports the zero matrix as singular', () => {
        const result = ge.compute(3, new Array<number>(9).fill(0),
            { wantInverse: true });
        expect(result.invertible).toBe(false);
        expect(result.determinant).toBe(0);
        expect(result.inverseM).toEqual(new Array<number>(9).fill(0));
    });

    it('does not modify the input matrix', () => {
        const M = [4, 7, 2, 6];
        const copy = M.slice();
        ge.compute(2, M, { wantInverse: true });
        expect(M).toEqual(copy);
    });

    it('honors column-major storage', () => {
        const rowMajorM = [2, 0, 1, 1, 3, 2, 0, 1, 4];
        const colMajorM = transpose(3, rowMajorM);

        const asRow = ge.compute(3, rowMajorM, { wantInverse: true });
        const asCol = ge.compute(3, colMajorM,
            { wantInverse: true, rowMajor: false });

        expect(asCol.determinant).toBeCloseTo(asRow.determinant, 10);
        // The column-major inverse is the transpose of the row-major one.
        const expected = transpose(3, asRow.inverseM as number[]);
        const actual = asCol.inverseM as number[];
        for (let i = 0; i < 9; ++i) {
            expect(actual[i]).toBeCloseTo(expected[i], 10);
        }
    });

    it('rejects invalid input', () => {
        expect(() => ge.compute(0, [])).toThrow(/Invalid input/);
        expect(() => ge.compute(-1, [])).toThrow(/Invalid input/);
        expect(() => ge.compute(2, [1, 2, 3])).toThrow(/Invalid input/);
        expect(() => ge.compute(2, [1, 0, 0, 1], { C: [1, 1], numCols: 0 }))
            .toThrow(/Invalid input/);
    });

    it('inverts random matrices (randomized cross-check)', () => {
        const rand = makeRandom(2468);
        for (let n = 1; n <= 6; ++n) {
            for (let trial = 0; trial < 20; ++trial) {
                const M: number[] = [];
                for (let i = 0; i < n * n; ++i) {
                    M.push(2 * rand() - 1);
                }
                // Make the matrix diagonally dominant so it is well
                // conditioned and certainly invertible.
                for (let r = 0; r < n; ++r) {
                    M[r + n * r] += n;
                }

                const B: number[] = [];
                for (let i = 0; i < n; ++i) {
                    B.push(2 * rand() - 1);
                }

                const result = ge.compute(n, M, { wantInverse: true, B });
                expect(result.invertible).toBe(true);

                const identity = multiply(n, n, M, result.inverseM as number[]);
                for (let r = 0; r < n; ++r) {
                    for (let c = 0; c < n; ++c) {
                        expect(identity[c + n * r])
                            .toBeCloseTo(r === c ? 1 : 0, 9);
                    }
                }

                const MX = multiply(n, 1, M, result.X as number[]);
                for (let i = 0; i < n; ++i) {
                    expect(MX[i]).toBeCloseTo(B[i], 9);
                }

                const reference = cofactorDeterminant(n, M);
                expect(result.determinant / reference).toBeCloseTo(1, 8);
            }
        }
    });
});
