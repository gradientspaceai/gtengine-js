import { describe, it, expect } from 'vitest';
import { GaussianElimination } from '../src/GaussianElimination.js';
import { Matrix } from '../src/Matrix.js';
import { check, expectClose, fc, invertibleMatrix, scaled } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (group V38).
// ---------------------------------------------------------------------------
describe('GaussianElimination verification', () => {
    const ge = new GaussianElimination();

    // Row-major flat array of an n x n Matrix.
    function flat(m: Matrix): number[] {
        const a: number[] = [];
        for (let r = 0; r < m.numRows; ++r) {
            for (let c = 0; c < m.numCols; ++c) { a.push(m.get(r, c)); }
        }
        return a;
    }

    function mulFlat(n: number, a: readonly number[], b: readonly number[],
        bCols: number): number[] {
        const out = new Array<number>(n * bCols).fill(0);
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < bCols; ++c) {
                let sum = 0;
                for (let k = 0; k < n; ++k) {
                    sum += a[r * n + k] * b[k * bCols + c];
                }
                out[r * bCols + c] = sum;
            }
        }
        return out;
    }

    // ||A||_inf, used to scale residual tolerances.
    function normInf(n: number, a: readonly number[]): number {
        let best = 0;
        for (let r = 0; r < n; ++r) {
            let sum = 0;
            for (let c = 0; c < n; ++c) { sum += Math.abs(a[r * n + c]); }
            best = Math.max(best, sum);
        }
        return best;
    }

    const sizedInvertible = (n: number) => invertibleMatrix(n, 1e-2)
        .map(m => ({ n, a: flat(m) }));

    const anyInvertible = fc.oneof(sizedInvertible(2), sizedInvertible(3),
        sizedInvertible(4));

    it('the inverse is a two-sided inverse with a conditioned residual', () => {
        check(anyInvertible, ({ n, a }) => {
            const r = ge.compute(n, a, { wantInverse: true });
            expect(r.invertible).toBe(true);
            const inv = r.inverseM as number[];
            // ||A^-1 A - I|| <= c * cond(A) * eps. Estimate cond with the
            // infinity norms of A and of the computed inverse.
            const cond = normInf(n, a) * normInf(n, inv);
            const tol = 1e-12 * Math.max(cond, 1);
            const left = mulFlat(n, a, inv, n);
            const right = mulFlat(n, inv, a, n);
            for (let i = 0; i < n * n; ++i) {
                const target = (i % n === Math.floor(i / n)) ? 1 : 0;
                expect(Math.abs(left[i] - target)).toBeLessThanOrEqual(tol);
                expect(Math.abs(right[i] - target)).toBeLessThanOrEqual(tol);
            }
        });
    });

    it('the determinant matches a cofactor expansion', () => {
        check(fc.oneof(sizedInvertible(2), sizedInvertible(3)), ({ n, a }) => {
            const r = ge.compute(n, a, {});
            const expected = n === 2
                ? a[0] * a[3] - a[1] * a[2]
                : a[0] * (a[4] * a[8] - a[5] * a[7])
                - a[1] * (a[3] * a[8] - a[5] * a[6])
                + a[2] * (a[3] * a[7] - a[4] * a[6]);
            const scale = Math.pow(normInf(n, a), n);
            expectClose(r.determinant, expected, 1e-11 * scale, 1e-11);
        });
    });

    it('the determinant is multiplicative and flips sign on a row swap', () => {
        check(fc.tuple(sizedInvertible(3), sizedInvertible(3)),
            ([{ a }, { a: b }]) => {
                const n = 3;
                const detA = ge.compute(n, a, {}).determinant;
                const detB = ge.compute(n, b, {}).determinant;
                const detAB = ge.compute(n, mulFlat(n, a, b, n), {}).determinant;
                const scale = Math.pow(normInf(n, a) * normInf(n, b), n);
                expectClose(detAB, detA * detB, 1e-10 * scale, 1e-10);

                // Swapping two rows negates the determinant exactly (full
                // pivoting picks the same pivots, only the parity changes).
                const swapped = a.slice();
                for (let c = 0; c < n; ++c) {
                    swapped[0 * n + c] = a[1 * n + c];
                    swapped[1 * n + c] = a[0 * n + c];
                }
                expectClose(ge.compute(n, swapped, {}).determinant, -detA,
                    1e-11 * Math.pow(normInf(n, a), n), 1e-11);
            }, 100);
    });

    it('solves M*X = B and M*Y = C consistently with the inverse', () => {
        check(fc.tuple(anyInvertible,
            fc.array(scaled(-6, 6), { minLength: 12, maxLength: 12 })),
            ([{ n, a }, pool]) => {
                const b = pool.slice(0, n);
                const numCols = 2;
                // C in row-major order.
                const c = pool.slice(0, n * numCols);

                const r = ge.compute(n, a, {
                    wantInverse: true, B: b, C: c, numCols
                });
                expect(r.invertible).toBe(true);
                const x = r.X as number[];
                const y = r.Y as number[];
                const inv = r.inverseM as number[];
                const cond = normInf(n, a) * normInf(n, inv);
                const tol = 1e-11 * Math.max(cond, 1)
                    * Math.max(...b.map(Math.abs), ...c.map(Math.abs), 1);

                // M*X = B and M*Y = C.
                const mx = mulFlat(n, a, x, 1);
                for (let i = 0; i < n; ++i) {
                    expect(Math.abs(mx[i] - b[i])).toBeLessThanOrEqual(tol);
                }
                const my = mulFlat(n, a, y, numCols);
                for (let i = 0; i < n * numCols; ++i) {
                    expect(Math.abs(my[i] - c[i])).toBeLessThanOrEqual(tol);
                }
                // X = M^-1 * B.
                const ib = mulFlat(n, inv, b, 1);
                for (let i = 0; i < n; ++i) {
                    expect(Math.abs(x[i] - ib[i])).toBeLessThanOrEqual(tol);
                }
                // Requesting only some outputs leaves the others null.
                const only = ge.compute(n, a, { B: b });
                expect(only.inverseM).toBeNull();
                expect(only.Y).toBeNull();
                expect(only.X).not.toBeNull();
            });
    });

    it('column-major storage transposes the problem', () => {
        check(sizedInvertible(3), ({ n, a }) => {
            // Reading a row-major buffer as column major is reading A^T, so
            // the column-major inverse of the same buffer is the row-major
            // inverse of A^T, that is (A^-1)^T laid out column major, which
            // is (A^-1) laid out row major.
            const rowMajor = ge.compute(n, a, { wantInverse: true });
            const colMajor = ge.compute(n, a,
                { wantInverse: true, rowMajor: false });
            expect(rowMajor.invertible).toBe(true);
            expect(colMajor.invertible).toBe(true);
            const scale = 1e-9 * Math.max(
                ...(rowMajor.inverseM as number[]).map(Math.abs), 1);
            for (let i = 0; i < n * n; ++i) {
                expectClose((rowMajor.inverseM as number[])[i],
                    (colMajor.inverseM as number[])[i], scale, 1e-9);
            }
            expectClose(rowMajor.determinant, colMajor.determinant,
                1e-9 * Math.abs(rowMajor.determinant) + 1e-12, 1e-9);
        });
    });

    it('a singular matrix zero-fills every requested output', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 4 }),
            fc.array(scaled(-5, 5), { minLength: 16, maxLength: 16 }),
            fc.array(scaled(-5, 5), { minLength: 8, maxLength: 8 })),
            ([n, pool, rhs]) => {
                // A matrix with a zero row is singular and is detected
                // exactly: the elimination subtracts multiples of the pivot
                // row scaled by that row's own (zero) entry, so the zero row
                // stays zero and the pivot search eventually finds only
                // zeros. (A matrix with two equal rows is mathematically
                // singular but not detected: 1/pivot*pivot is not exactly 1,
                // so the eliminated row keeps round-off residuals and the
                // pivot search still finds a nonzero maximum.)
                const a = new Array<number>(n * n).fill(0);
                for (let r = 0; r < n; ++r) {
                    if (r === 1) { continue; }
                    for (let c = 0; c < n; ++c) {
                        a[r * n + c] = pool[r * n + c];
                    }
                }
                const b = rhs.slice(0, n);
                const r = ge.compute(n, a, {
                    wantInverse: true, B: b, C: b, numCols: 1
                });
                expect(r.invertible).toBe(false);
                expect(r.determinant).toBe(0);
                expect(r.inverseM).toEqual(new Array<number>(n * n).fill(0));
                expect(r.X).toEqual(new Array<number>(n).fill(0));
                expect(r.Y).toEqual(new Array<number>(n).fill(0));
            });
    });

    it('never modifies the input arrays', () => {
        check(fc.tuple(anyInvertible,
            fc.array(scaled(-6, 6), { minLength: 8, maxLength: 8 })),
            ([{ n, a }, pool]) => {
                const aCopy = a.slice();
                const b = pool.slice(0, n);
                const bCopy = b.slice();
                ge.compute(n, a, { wantInverse: true, B: b, C: b, numCols: 1 });
                expect(a).toEqual(aCopy);
                expect(b).toEqual(bCopy);
            });
    });

    it('rejects invalid input', () => {
        check(fc.integer({ min: -3, max: 0 }), n => {
            expect(() => ge.compute(n, [1])).toThrow('Invalid input.');
        });
        expect(() => ge.compute(3, [1, 2, 3])).toThrow('Invalid input.');
        expect(() => ge.compute(2, [1, 0, 0, 1], { C: [1, 1], numCols: 0 }))
            .toThrow('Invalid input.');
    });

    it('upstream defect (preserved, #375): a subnormal pivot reports '
        + 'invertible with NaN entries', () => {
            // The pivot search only rejects an exactly zero maximum entry, so
            // a subnormal pivot passes, 1/pivot overflows and the row scaling
            // produces infinity * 0 = NaN.
            const denormal = ge.compute(2,
                [5e-324, 1e-323, 3e-324, 7e-324], { wantInverse: true });
            expect(denormal.invertible).toBe(true);
            // The determinant is not a number the caller can use either: the
            // scaled pivots overflow.
            expect(Number.isFinite(denormal.determinant)).toBe(false);
            expect((denormal.inverseM as number[]).some(x => !Number.isFinite(x)))
                .toBe(true);

            // Why the port does not "fix" this by reporting non-invertible:
            // full pivoting takes the largest entry first, so a subnormal
            // pivot does not imply an underflowing determinant. Here the
            // determinant is representable and correct while the inverse
            // (about 1e320) is not.
            const mixed = ge.compute(2, [1e300, 0, 0, 1e-320],
                { wantInverse: true });
            expect(mixed.invertible).toBe(true);
            expectClose(mixed.determinant, 1e300 * 1e-320, 0, 1e-10);
            expect((mixed.inverseM as number[]).some(x => !Number.isFinite(x)))
                .toBe(true);
        });
});
