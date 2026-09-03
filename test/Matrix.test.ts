import { describe, it, expect } from 'vitest';
import {
    check, finite, nonzero, vector, matrix, wellScaledMatrix, expectClose, fc
} from './helpers/arbitraries.js';
import {
    Matrix, negateMatrix, addMatrix, subMatrix, mulMatrix, divMatrix,
    l1Norm, l2Norm, lInfinityNorm, inverse, determinant, transpose,
    multiplyAB, multiplyABT, multiplyATB, multiplyATBT, multiplyMD,
    multiplyDM, outerProduct, makeDiagonal, hliftMatrix, hprojectMatrix
} from '../src/Matrix.js';
import { Vector, dot, length } from '../src/Vector.js';

function expectMatrixClose(actual: Matrix, expected: Matrix,
    tolerance: number = 1e-12): void {
    expect(actual.numRows).toBe(expected.numRows);
    expect(actual.numCols).toBe(expected.numCols);
    for (let i = 0; i < expected.numElements; ++i) {
        expect(Math.abs(actual.values[i] - expected.values[i]))
            .toBeLessThanOrEqual(tolerance);
    }
}

// A deterministic pseudorandom generator so the randomized cross-checks are
// reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296 - 0.5;
    };
}

function randomMatrix(numRows: number, numCols: number,
    rand: () => number): Matrix {
    const M = new Matrix(numRows, numCols);
    for (let i = 0; i < M.numElements; ++i) {
        M.values[i] = 4 * rand();
    }
    return M;
}

// An independent determinant by cofactor expansion, used to cross-check the
// Gaussian-elimination-based determinant().
function cofactorDeterminant(M: Matrix): number {
    const n = M.numRows;
    if (n === 1) {
        return M.get(0, 0);
    }

    let result = 0;
    let sign = 1;
    for (let c = 0; c < n; ++c) {
        const minor = new Matrix(n - 1, n - 1);
        for (let r1 = 1; r1 < n; ++r1) {
            for (let c1 = 0, c2 = 0; c1 < n; ++c1) {
                if (c1 !== c) {
                    minor.set(r1 - 1, c2, M.get(r1, c1));
                    ++c2;
                }
            }
        }
        result += sign * M.get(0, c) * cofactorDeterminant(minor);
        sign = -sign;
    }
    return result;
}

describe('Matrix construction and special matrices', () => {
    it('the constructor produces a zero-filled table', () => {
        const M = new Matrix(2, 3);
        expect(M.numRows).toBe(2);
        expect(M.numCols).toBe(3);
        expect(M.numElements).toBe(6);
        expect(M.values).toEqual([0, 0, 0, 0, 0, 0]);
        expect(M.getSize()).toEqual({ numRows: 2, numCols: 3 });
        expect(M.getNumRows()).toBe(2);
        expect(M.getNumCols()).toBe(3);
        expect(M.getNumElements()).toBe(6);
    });

    it('fromArray copies row-major values', () => {
        const M = Matrix.fromArray(2, 3, [1, 2, 3, 4, 5, 6]);
        expect(M.get(0, 0)).toBe(1);
        expect(M.get(0, 2)).toBe(3);
        expect(M.get(1, 0)).toBe(4);
        expect(M.get(1, 2)).toBe(6);
    });

    it('fromArray zero-fills the remaining elements, as the upstream '
        + 'initializer-list constructor does', () => {
        const M = Matrix.fromArray(2, 3, [1, 2, 3, 4]);
        expect(M.values).toEqual([1, 2, 3, 4, 0, 0]);
    });

    it('fromArray ignores extra values', () => {
        const M = Matrix.fromArray(2, 2, [1, 2, 3, 4, 5, 6]);
        expect(M.values).toEqual([1, 2, 3, 4]);
    });

    it('zero, unit and identity', () => {
        expect(Matrix.zero(2, 2).values).toEqual([0, 0, 0, 0]);
        expect(Matrix.unit(2, 3, 1, 2).values).toEqual([0, 0, 0, 0, 0, 1]);
        // Nonsquare identity has ones on the diagonal only.
        expect(Matrix.identity(2, 3).values).toEqual([1, 0, 0, 0, 1, 0]);
        expect(Matrix.identity(3, 2).values).toEqual([1, 0, 0, 1, 0, 0]);
    });

    it('makeUnit with an invalid index produces the zero matrix', () => {
        const M = Matrix.fromArray(2, 2, [1, 2, 3, 4]);
        M.makeUnit(2, 0);
        expect(M.values).toEqual([0, 0, 0, 0]);
        M.makeUnit(0, -1);
        expect(M.values).toEqual([0, 0, 0, 0]);
    });

    it('makeZero and makeIdentity', () => {
        const M = Matrix.fromArray(2, 2, [1, 2, 3, 4]);
        M.makeIdentity();
        expect(M.values).toEqual([1, 0, 0, 1]);
        M.makeZero();
        expect(M.values).toEqual([0, 0, 0, 0]);
    });

    it('a negative dimension throws', () => {
        expect(() => new Matrix(-1, 2)).toThrow('Invalid size.');
    });
});

describe('Matrix element and row/column access', () => {
    it('get/set and getFlat/setFlat agree with row-major storage', () => {
        const M = new Matrix(2, 3);
        M.set(1, 2, 7);
        expect(M.get(1, 2)).toBe(7);
        expect(M.getFlat(5)).toBe(7);
        M.setFlat(0, 9);
        expect(M.get(0, 0)).toBe(9);
    });

    it('getRow and getCol', () => {
        const M = Matrix.fromArray(2, 3, [1, 2, 3, 4, 5, 6]);
        expect(M.getRow(1).values).toEqual([4, 5, 6]);
        expect(M.getCol(2).values).toEqual([3, 6]);
    });

    it('setRow and setCol', () => {
        const M = new Matrix(2, 3);
        M.setRow(0, Vector.fromArray([1, 2, 3]));
        M.setCol(1, Vector.fromArray([8, 9]));
        expect(M.values).toEqual([1, 8, 3, 0, 9, 0]);
    });

    it('clone is a deep copy', () => {
        const M = Matrix.fromArray(2, 2, [1, 2, 3, 4]);
        const N = M.clone();
        N.set(0, 0, 100);
        expect(M.get(0, 0)).toBe(1);
        expect(N.get(0, 0)).toBe(100);
    });
});

describe('Matrix comparisons', () => {
    it('are lexicographic over the row-major elements', () => {
        const A = Matrix.fromArray(2, 2, [1, 2, 3, 4]);
        const B = Matrix.fromArray(2, 2, [1, 2, 3, 5]);
        expect(A.equals(A.clone())).toBe(true);
        expect(A.notEquals(B)).toBe(true);
        expect(A.lessThan(B)).toBe(true);
        expect(A.lessThanOrEqual(B)).toBe(true);
        expect(B.greaterThan(A)).toBe(true);
        expect(B.greaterThanOrEqual(A)).toBe(true);
        expect(A.greaterThan(B)).toBe(false);
        expect(A.lessThanOrEqual(A.clone())).toBe(true);
        expect(A.greaterThanOrEqual(A.clone())).toBe(true);
    });

    it('throw on mismatched dimensions (a compile error upstream)', () => {
        const A = new Matrix(2, 2);
        const B = new Matrix(2, 3);
        expect(() => A.equals(B)).toThrow('Matrix: mismatched sizes.');
    });
});

describe('Matrix arithmetic', () => {
    const A = Matrix.fromArray(2, 2, [1, 2, 3, 4]);
    const B = Matrix.fromArray(2, 2, [5, 6, 7, 8]);

    it('negate', () => {
        expect(negateMatrix(A).values).toEqual([-1, -2, -3, -4]);
    });

    it('add and sub', () => {
        expect(addMatrix(A, B).values).toEqual([6, 8, 10, 12]);
        expect(subMatrix(B, A).values).toEqual([4, 4, 4, 4]);
        expect(subMatrix(A, A).values).toEqual([0, 0, 0, 0]);
    });

    it('add and sub throw on mismatched sizes', () => {
        const C = new Matrix(3, 2);
        expect(() => addMatrix(A, C)).toThrow('Mismatched sizes');
        expect(() => subMatrix(A, C)).toThrow('Mismatched sizes');
    });

    it('scalar multiplication accepts either argument order', () => {
        expect(mulMatrix(A, 2).values).toEqual([2, 4, 6, 8]);
        expect(mulMatrix(2, A).values).toEqual([2, 4, 6, 8]);
    });

    it('scalar division; division by zero produces the zero matrix', () => {
        expect(divMatrix(A, 2).values).toEqual([0.5, 1, 1.5, 2]);
        expect(divMatrix(A, 0).values).toEqual([0, 0, 0, 0]);
    });

    it('identities: A + (-A) = 0 and (A + B) - B = A', () => {
        expect(addMatrix(A, negateMatrix(A)).values).toEqual([0, 0, 0, 0]);
        expectMatrixClose(subMatrix(addMatrix(A, B), B), A);
    });
});

describe('Matrix products', () => {
    it('M*V and V^T*M with known values', () => {
        const M = Matrix.fromArray(2, 3, [1, 2, 3, 4, 5, 6]);
        const v3 = Vector.fromArray([1, 0, -1]);
        const v2 = Vector.fromArray([1, 2]);
        expect(mulMatrix(M, v3).values).toEqual([-2, -2]);
        expect(mulMatrix(v2, M).values).toEqual([9, 12, 15]);
    });

    it('M*V and V^T*M throw on mismatched sizes', () => {
        const M = Matrix.fromArray(2, 3, [1, 2, 3, 4, 5, 6]);
        expect(() => mulMatrix(M, Vector.fromArray([1, 2])))
            .toThrow('Mismatched sizes.');
        expect(() => mulMatrix(Vector.fromArray([1, 2, 3]), M))
            .toThrow('Mismatched sizes.');
    });

    it('the identity acts as the multiplicative identity', () => {
        const M = Matrix.fromArray(2, 3, [1, 2, 3, 4, 5, 6]);
        expectMatrixClose(multiplyAB(Matrix.identity(2, 2), M), M);
        expectMatrixClose(multiplyAB(M, Matrix.identity(3, 3)), M);
        const v = Vector.fromArray([1, 2, 3]);
        expect(mulMatrix(Matrix.identity(3, 3), v).values).toEqual(v.values);
    });

    it('A*B with known values; mulMatrix forwards to multiplyAB', () => {
        const A = Matrix.fromArray(2, 3, [1, 2, 3, 4, 5, 6]);
        const B = Matrix.fromArray(3, 2, [7, 8, 9, 10, 11, 12]);
        const AB = multiplyAB(A, B);
        expect(AB.numRows).toBe(2);
        expect(AB.numCols).toBe(2);
        expect(AB.values).toEqual([58, 64, 139, 154]);
        expect(mulMatrix(A, B).values).toEqual(AB.values);
    });

    it('A*B throws on mismatched sizes', () => {
        const A = new Matrix(2, 3);
        expect(() => multiplyAB(A, new Matrix(2, 3)))
            .toThrow('Mismatched sizes.');
    });

    it('matrix multiplication is associative and distributive', () => {
        const rand = makeRandom(12345);
        const A = randomMatrix(3, 4, rand);
        const B = randomMatrix(4, 2, rand);
        const C = randomMatrix(2, 3, rand);
        const D = randomMatrix(4, 2, rand);
        expectMatrixClose(multiplyAB(multiplyAB(A, B), C),
            multiplyAB(A, multiplyAB(B, C)), 1e-12);
        expectMatrixClose(multiplyAB(A, addMatrix(B, D)),
            addMatrix(multiplyAB(A, B), multiplyAB(A, D)), 1e-12);
    });

    it('the transposed products agree with transpose + multiplyAB', () => {
        const rand = makeRandom(999);
        const A = randomMatrix(3, 4, rand);
        const B = randomMatrix(2, 4, rand);
        const C = randomMatrix(3, 2, rand);
        expectMatrixClose(multiplyABT(A, B), multiplyAB(A, transpose(B)));
        expectMatrixClose(multiplyATB(A, C), multiplyAB(transpose(A), C));
        const D = randomMatrix(2, 3, rand);
        expectMatrixClose(multiplyATBT(A, D),
            multiplyAB(transpose(A), transpose(D)));
    });

    it('the transposed products throw on mismatched sizes', () => {
        const A = new Matrix(3, 4);
        expect(() => multiplyABT(A, new Matrix(2, 3)))
            .toThrow('Mismatched sizes.');
        expect(() => multiplyATB(A, new Matrix(2, 3)))
            .toThrow('Mismatched sizes.');
        expect(() => multiplyATBT(A, new Matrix(3, 2)))
            .toThrow('Mismatched sizes.');
    });

    it('multiplyMD and multiplyDM agree with multiplication by a diagonal '
        + 'matrix', () => {
        const M = Matrix.fromArray(2, 3, [1, 2, 3, 4, 5, 6]);
        const dCols = Vector.fromArray([2, 3, 4]);
        const dRows = Vector.fromArray([5, 6]);
        const Dc = new Matrix(3, 3);
        makeDiagonal(dCols, Dc);
        const Dr = new Matrix(2, 2);
        makeDiagonal(dRows, Dr);
        expectMatrixClose(multiplyMD(M, dCols), multiplyAB(M, Dc));
        expectMatrixClose(multiplyDM(dRows, M), multiplyAB(Dr, M));
        expect(multiplyMD(M, dCols).values).toEqual([2, 6, 12, 8, 15, 24]);
        expect(multiplyDM(dRows, M).values).toEqual([5, 10, 15, 24, 30, 36]);
    });

    it('multiplyMD and multiplyDM throw on mismatched sizes', () => {
        const M = new Matrix(2, 3);
        expect(() => multiplyMD(M, Vector.fromArray([1, 2])))
            .toThrow('Mismatched sizes.');
        expect(() => multiplyDM(Vector.fromArray([1, 2, 3]), M))
            .toThrow('Mismatched sizes.');
    });

    it('outerProduct is U*V^T', () => {
        const U = Vector.fromArray([1, 2, 3]);
        const V = Vector.fromArray([4, 5]);
        const M = outerProduct(U, V);
        expect(M.numRows).toBe(3);
        expect(M.numCols).toBe(2);
        expect(M.values).toEqual([4, 5, 8, 10, 12, 15]);
        // The outer product of nonzero vectors has rank 1, so a 2x2 minor
        // has zero determinant.
        expect(determinant(Matrix.fromArray(2, 2,
            [M.get(0, 0), M.get(0, 1), M.get(1, 0), M.get(1, 1)])))
            .toBeCloseTo(0, 12);
    });
});

describe('Matrix norms', () => {
    const M = Matrix.fromArray(2, 2, [3, -4, 0, 12]);

    it('L1, L2 and L-infinity norms have their known values', () => {
        expect(l1Norm(M)).toBe(19);
        expect(l2Norm(M)).toBeCloseTo(Math.sqrt(9 + 16 + 144), 12);
        expect(lInfinityNorm(M)).toBe(12);
    });

    it('the norms of the zero matrix are zero', () => {
        const Z = new Matrix(3, 3);
        expect(l1Norm(Z)).toBe(0);
        expect(l2Norm(Z)).toBe(0);
        expect(lInfinityNorm(Z)).toBe(0);
    });

    it('the L2 norm of the NxN identity is sqrt(N)', () => {
        expect(l2Norm(Matrix.identity(4, 4))).toBeCloseTo(2, 12);
    });
});

describe('Matrix transpose', () => {
    it('swaps the dimensions and the elements', () => {
        const M = Matrix.fromArray(2, 3, [1, 2, 3, 4, 5, 6]);
        const T = transpose(M);
        expect(T.numRows).toBe(3);
        expect(T.numCols).toBe(2);
        expect(T.values).toEqual([1, 4, 2, 5, 3, 6]);
    });

    it('is an involution', () => {
        const rand = makeRandom(4242);
        const M = randomMatrix(3, 5, rand);
        expectMatrixClose(transpose(transpose(M)), M, 0);
    });

    it('(A*B)^T = B^T * A^T', () => {
        const rand = makeRandom(777);
        const A = randomMatrix(3, 4, rand);
        const B = randomMatrix(4, 2, rand);
        expectMatrixClose(transpose(multiplyAB(A, B)),
            multiplyAB(transpose(B), transpose(A)));
    });
});

describe('Matrix determinant', () => {
    it('has its known values', () => {
        expect(determinant(Matrix.fromArray(1, 1, [5]))).toBeCloseTo(5, 12);
        expect(determinant(Matrix.fromArray(2, 2, [1, 2, 3, 4])))
            .toBeCloseTo(-2, 12);
        expect(determinant(Matrix.fromArray(3, 3,
            [6, 1, 1, 4, -2, 5, 2, 8, 7]))).toBeCloseTo(-306, 12);
        expect(determinant(Matrix.identity(4, 4))).toBeCloseTo(1, 12);
    });

    it('is zero for a singular matrix', () => {
        expect(determinant(Matrix.fromArray(2, 2, [1, 2, 2, 4]))).toBe(0);
        expect(determinant(new Matrix(3, 3))).toBe(0);
    });

    it('agrees with cofactor expansion on random matrices', () => {
        const rand = makeRandom(20260830);
        for (let trial = 0; trial < 20; ++trial) {
            const n = 2 + (trial % 3);
            const M = randomMatrix(n, n, rand);
            expect(determinant(M)).toBeCloseTo(cofactorDeterminant(M), 10);
        }
    });

    it('is multiplicative: det(A*B) = det(A)*det(B)', () => {
        const rand = makeRandom(31337);
        const A = randomMatrix(3, 3, rand);
        const B = randomMatrix(3, 3, rand);
        expect(determinant(multiplyAB(A, B)))
            .toBeCloseTo(determinant(A) * determinant(B), 10);
    });

    it('throws for a nonsquare matrix', () => {
        expect(() => determinant(new Matrix(2, 3)))
            .toThrow('Matrix must be square.');
    });
});

describe('Matrix inverse', () => {
    it('has its known value for a 2x2 matrix', () => {
        const M = Matrix.fromArray(2, 2, [1, 2, 3, 4]);
        const result = inverse(M);
        expect(result.invertible).toBe(true);
        expectMatrixClose(result.inverse,
            Matrix.fromArray(2, 2, [-2, 1, 1.5, -0.5]), 1e-12);
    });

    it('round-trips: M * inverse(M) = inverse(M) * M = I', () => {
        const rand = makeRandom(8675309);
        for (let trial = 0; trial < 10; ++trial) {
            const n = 2 + (trial % 4);
            const M = randomMatrix(n, n, rand);
            const result = inverse(M);
            if (!result.invertible) {
                continue;
            }
            expectMatrixClose(multiplyAB(M, result.inverse),
                Matrix.identity(n, n), 1e-9);
            expectMatrixClose(multiplyAB(result.inverse, M),
                Matrix.identity(n, n), 1e-9);
        }
    });

    it('reports noninvertibility and returns the zero matrix', () => {
        const result = inverse(Matrix.fromArray(2, 2, [1, 2, 2, 4]));
        expect(result.invertible).toBe(false);
        expect(result.inverse.values).toEqual([0, 0, 0, 0]);
    });

    it('the inverse of the identity is the identity', () => {
        const result = inverse(Matrix.identity(3, 3));
        expect(result.invertible).toBe(true);
        expectMatrixClose(result.inverse, Matrix.identity(3, 3), 1e-15);
    });

    it('throws for a nonsquare matrix', () => {
        expect(() => inverse(new Matrix(2, 3)))
            .toThrow('Matrix must be square.');
    });
});

describe('Matrix makeDiagonal, hlift and hproject', () => {
    it('makeDiagonal fills the diagonal, even when nonsquare', () => {
        const M = Matrix.fromArray(2, 3, [1, 2, 3, 4, 5, 6]);
        makeDiagonal(Vector.fromArray([7, 8]), M);
        expect(M.values).toEqual([7, 0, 0, 0, 8, 0]);

        const S = new Matrix(3, 3);
        makeDiagonal(Vector.fromArray([1, 2, 3]), S);
        expect(S.values).toEqual([1, 0, 0, 0, 2, 0, 0, 0, 3]);
    });

    it('hliftMatrix embeds the matrix in the upper block', () => {
        const M = Matrix.fromArray(2, 2, [1, 2, 3, 4]);
        const H = hliftMatrix(M);
        expect(H.numRows).toBe(3);
        expect(H.numCols).toBe(3);
        expect(H.values).toEqual([1, 2, 0, 3, 4, 0, 0, 0, 1]);
    });

    it('hprojectMatrix extracts the upper block; hproject(hlift(M)) = M', () => {
        const M = Matrix.fromArray(2, 2, [1, 2, 3, 4]);
        expectMatrixClose(hprojectMatrix(hliftMatrix(M)), M, 0);
        expect(hprojectMatrix(Matrix.fromArray(3, 3,
            [1, 2, 3, 4, 5, 6, 7, 8, 9])).values).toEqual([1, 2, 4, 5]);
    });

    it('hliftMatrix of an affine transform composes as the transforms do',
        () => {
            const rand = makeRandom(5150);
            const A = randomMatrix(3, 3, rand);
            const B = randomMatrix(3, 3, rand);
            expectMatrixClose(hliftMatrix(multiplyAB(A, B)),
                multiplyAB(hliftMatrix(A), hliftMatrix(B)), 1e-12);
        });

    it('hprojectMatrix throws for dimension 1', () => {
        expect(() => hprojectMatrix(new Matrix(1, 1)))
            .toThrow('Invalid matrix dimension.');
    });

    it('hliftMatrix and hprojectMatrix require square matrices', () => {
        expect(() => hliftMatrix(new Matrix(2, 3)))
            .toThrow('Matrix must be square.');
        expect(() => hprojectMatrix(new Matrix(2, 3)))
            .toThrow('Matrix must be square.');
    });
});

// ---------------------------------------------------------------------------
// Verification wave: property-based checks against upstream Matrix.h.
// ---------------------------------------------------------------------------

describe('Matrix verification', () => {
    it('storage is row major and every accessor agrees with it', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 4 }),
            fc.integer({ min: 1, max: 4 })).chain(([r, c]) =>
            fc.tuple(fc.constant(r), fc.constant(c), matrix(r, c))),
            ([numRows, numCols, M]) => {
                for (let r = 0; r < numRows; ++r) {
                    for (let c = 0; c < numCols; ++c) {
                        const i = c + numCols * r;
                        // M(r, c) and the storage-order-independent M[i] are
                        // the same element under the port's row-major choice.
                        expect(M.get(r, c)).toBe(M.values[i]);
                        expect(M.getFlat(i)).toBe(M.values[i]);
                        expect(M.getRow(r).get(c)).toBe(M.get(r, c));
                        expect(M.getCol(c).get(r)).toBe(M.get(r, c));
                    }
                }
                expect(M.numElements).toBe(numRows * numCols);
                expect(M.getSize()).toEqual({ numRows, numCols });
            });
    });

    it('fromArray truncates and zero-fills like the initializer-list ctor',
        () => {
            check(fc.tuple(fc.integer({ min: 1, max: 3 }),
                fc.integer({ min: 1, max: 3 }),
                fc.array(finite(), { minLength: 0, maxLength: 12 })),
                ([numRows, numCols, vals]) => {
                    const M = Matrix.fromArray(numRows, numCols, vals);
                    for (let i = 0; i < numRows * numCols; ++i) {
                        expect(M.values[i]).toBe(i < vals.length ? vals[i] : 0);
                    }
                });
        });

    it('setRow/setCol invert getRow/getCol', () => {
        check(fc.tuple(matrix(3, 4), vector(4), vector(3),
            fc.integer({ min: 0, max: 2 }), fc.integer({ min: 0, max: 3 })),
            ([M0, row, col, r, c]) => {
                const M = M0.clone();
                M.setRow(r, row);
                expect(M.getRow(r).values).toEqual(row.values);
                M.setCol(c, col);
                expect(M.getCol(c).values).toEqual(col.values);
                // setCol overwrote one element of the row just written.
                expect(M.get(r, c)).toBe(col.get(r));
            });
    });

    it('the transposed products equal transpose + multiplyAB exactly', () => {
        // The accumulation order of each entry is identical, so these agree
        // bit for bit and no tolerance is needed.
        check(fc.tuple(matrix(3, 2), matrix(4, 2), matrix(3, 4)),
            ([A32, B42, C34]) => {
                expect(multiplyABT(A32, B42).values)
                    .toEqual(multiplyAB(A32, transpose(B42)).values);
                expect(multiplyATB(A32, C34).values)
                    .toEqual(multiplyAB(transpose(A32), C34).values);
                expect(multiplyATBT(A32, transpose(C34)).values)
                    .toEqual(multiplyAB(transpose(A32), C34).values);
            });
    });

    it('multiplyMD/multiplyDM/outerProduct match their definitions', () => {
        check(fc.tuple(matrix(3, 4), vector(4), vector(3)),
            ([M, dCols, dRows]) => {
                const MD = multiplyMD(M, dCols);
                const DM = multiplyDM(dRows, M);
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 4; ++c) {
                        expect(MD.get(r, c)).toBe(M.get(r, c) * dCols.get(c));
                        expect(DM.get(r, c)).toBe(dRows.get(r) * M.get(r, c));
                    }
                }
                // MD = M * diag(d): the inner sums have a single nonzero term,
                // so the agreement is exact.
                const diag = new Matrix(4, 4);
                makeDiagonal(dCols, diag);
                expect(MD.values.map(x => x + 0))
                    .toEqual(multiplyAB(M, diag).values.map(x => x + 0));
                const outer = outerProduct(dRows, dCols);
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 4; ++c) {
                        expect(outer.get(r, c))
                            .toBe(dRows.get(r) * dCols.get(c));
                    }
                }
            });
    });

    it('M*V is the linear combination of columns, V^T*M of rows', () => {
        check(fc.tuple(matrix(3, 4), vector(4), vector(3)), ([M, x, y]) => {
            const Mx = mulMatrix(M, x) as Vector;
            const yM = mulMatrix(y, M) as Vector;
            expect(Mx.size).toBe(3);
            expect(yM.size).toBe(4);
            // Upstream M*V accumulates from an explicit zero while upstream
            // Dot seeds with v0[0]*v1[0]. Adding a leading exact zero changes
            // nothing but the sign of a zero result, which '+ 0' normalizes.
            for (let r = 0; r < 3; ++r) {
                expect(Mx.get(r) + 0).toBe(dot(M.getRow(r), x) + 0);
            }
            for (let c = 0; c < 4; ++c) {
                expect(yM.get(c) + 0).toBe(dot(y, M.getCol(c)) + 0);
            }
            // <y, M x> = <M^T y, x>: the same products summed in a different
            // order, so only up to rounding of the accumulation.
            expectClose(dot(y, Mx), dot(yM, x), 1e-9, 1e-11);
        });
    });

    it('transpose is an involution and reverses products', () => {
        check(fc.tuple(matrix(2, 3), matrix(3, 4)), ([A, B]) => {
            expect(transpose(transpose(A)).values).toEqual(A.values);
            expect(transpose(A).numRows).toBe(A.numCols);
            // (A*B)^T = B^T * A^T. The inner sums run over the same index in
            // the same order, but the two factors are swapped, so the products
            // are identical and the agreement is exact.
            expect(transpose(multiplyAB(A, B)).values)
                .toEqual(multiplyAB(transpose(B), transpose(A)).values);
        });
    });

    it('determinant agrees with cofactor expansion', () => {
        // Gaussian elimination with full pivoting versus a completely
        // different algorithm; the tolerance is relative to the Hadamard
        // bound (the product of the row norms), which bounds |det|.
        check(fc.constantFrom(2, 3, 4).chain(n =>
            fc.tuple(fc.constant(n), wellScaledMatrix(n, n))), ([n, M]) => {
            let bound = 1;
            for (let r = 0; r < n; ++r) { bound *= length(M.getRow(r)); }
            expectClose(determinant(M), cofactorDeterminant(M),
                1e-11 * bound, 0);
        });
    });

    it('determinant is multiplicative and inverse round-trips', () => {
        check(wellScaledMatrix(3, 3), M => {
            const d = determinant(M);
            if (Math.abs(d) < 1e-2) {
                return;   // ill conditioned: the inverse amplifies rounding
            }
            const res = inverse(M);
            expect(res.invertible).toBe(true);
            expectMatrixClose(multiplyAB(M, res.inverse),
                Matrix.identity(3, 3), 1e-8);
            expectMatrixClose(multiplyAB(res.inverse, M),
                Matrix.identity(3, 3), 1e-8);
            // det(A*M) = det(A)*det(M) with det(A) = 1 for this shear.
            const A = Matrix.identity(3, 3);
            A.set(0, 1, 2);
            A.set(2, 0, -3);
            expectClose(determinant(multiplyAB(A, M)), d, 1e-9,
                1e-9 * Math.abs(d));
        });
    });

    it('the norms satisfy their ordering and homogeneity', () => {
        check(fc.tuple(matrix(3, 3), finite()), ([M, s]) => {
            const l1 = l1Norm(M), l2 = l2Norm(M), li = lInfinityNorm(M);
            expect(li).toBeLessThanOrEqual(l2 + 1e-12);
            expect(l2).toBeLessThanOrEqual(l1 + 1e-12);
            // Absolute homogeneity, up to the rounding of the products.
            expectClose(l1Norm(mulMatrix(M, s) as Matrix), Math.abs(s) * l1,
                1e-12, 1e-12);
            expectClose(lInfinityNorm(mulMatrix(M, s) as Matrix),
                Math.abs(s) * li, 1e-12, 1e-12);
            // The L2 norm is the Frobenius norm, accumulated in index order.
            let sum = 0;
            for (const x of M.values) { sum += x * x; }
            expect(l2).toBe(Math.sqrt(sum));
        });
    });

    it('lInfinityNorm follows GMatrix.h, which ignores a leading NaN', () => {
        // Documented deviation. Upstream Matrix.h seeds the maximum with
        // fabs(M[0]), so a NaN in element 0 propagates while a NaN anywhere
        // else is silently ignored ('absElement > maxAbsElement' is false for
        // NaN). Upstream GMatrix.h seeds with 0 and therefore ignores every
        // NaN. The port has one shared implementation and follows GMatrix.h,
        // which also keeps the norm of the empty 0-by-0 GMatrix equal to 0.
        expect(lInfinityNorm(Matrix.fromArray(1, 2, [NaN, 5]))).toBe(5);
        expect(lInfinityNorm(Matrix.fromArray(1, 2, [5, NaN]))).toBe(5);
        expect(lInfinityNorm(new Matrix(0, 0))).toBe(0);
    });

    it('hliftMatrix and hprojectMatrix round trip at every dimension', () => {
        check(fc.constantFrom(1, 2, 3, 4).chain(n =>
            fc.tuple(fc.constant(n), matrix(n, n))), ([n, M]) => {
            const H = hliftMatrix(M);
            expect(H.numRows).toBe(n + 1);
            expect(H.get(n, n)).toBe(1);
            for (let i = 0; i < n; ++i) {
                expect(H.get(i, n)).toBe(0);
                expect(H.get(n, i)).toBe(0);
            }
            expect(hprojectMatrix(H).values).toEqual(M.values);
        });
    });

    it('makeDiagonal zeroes the off-diagonal, even when nonsquare', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 4 }),
            fc.integer({ min: 1, max: 4 })).chain(([r, c]) =>
            fc.tuple(fc.constant(r), fc.constant(c), matrix(r, c),
                vector(Math.min(r, c)))), ([numRows, numCols, M0, d]) => {
            const M = M0.clone();
            makeDiagonal(d, M);
            for (let r = 0; r < numRows; ++r) {
                for (let c = 0; c < numCols; ++c) {
                    expect(M.get(r, c)).toBe(r === c ? d.get(r) : 0);
                }
            }
        });
    });

    it('the arithmetic operators satisfy the vector-space identities', () => {
        check(fc.tuple(matrix(2, 3), matrix(2, 3), nonzero()),
            ([A, B, s]) => {
                const Z = new Matrix(2, 3);
                // '+ 0' normalizes the signed zeros of x + (-x).
                expect(addMatrix(A, negateMatrix(A)).values.map(x => x + 0))
                    .toEqual(Z.values);
                expect(subMatrix(A, A).values.map(x => x + 0))
                    .toEqual(Z.values);
                for (let i = 0; i < A.numElements; ++i) {
                    expect(addMatrix(A, B).values[i])
                        .toBe(A.values[i] + B.values[i]);
                    expect(subMatrix(A, B).values[i])
                        .toBe(A.values[i] - B.values[i]);
                }
                // Division is multiplication by the reciprocal, as upstream.
                expect((divMatrix(A, s)).values)
                    .toEqual((mulMatrix(A, 1 / s) as Matrix).values);
                expect(divMatrix(A, 0).values).toEqual(Z.values);
            });
    });

    it('comparisons are the lexicographic std::array ordering', () => {
        check(fc.tuple(matrix(2, 2), matrix(2, 2)), ([A, B]) => {
            let expected = 0;
            for (let i = 0; i < 4 && expected === 0; ++i) {
                expected = A.values[i] < B.values[i] ? -1
                    : (A.values[i] > B.values[i] ? 1 : 0);
            }
            expect(A.lessThan(B)).toBe(expected < 0);
            expect(A.greaterThan(B)).toBe(expected > 0);
            expect(A.lessThanOrEqual(B)).toBe(expected <= 0);
            expect(A.greaterThanOrEqual(B)).toBe(expected >= 0);
            expect(A.equals(B)).toBe(expected === 0);
            expect(A.notEquals(B)).toBe(expected !== 0);
        });
    });

    it('equals follows std::array: a NaN element breaks self-equality', () => {
        // Regression: the port compared lexicographically, which treats NaN
        // as equivalent to itself; C++ 'NaN == NaN' is false, so upstream's
        // 'mTable.mStorage == mat.mTable.mStorage' is false.
        const M = Matrix.fromArray(1, 2, [1, NaN]);
        expect(M.equals(M)).toBe(false);
        expect(M.notEquals(M)).toBe(true);
        expect(M.equals(M.clone())).toBe(false);
        // The ordering operators keep treating NaN as equivalent, exactly as
        // std::lexicographical_compare does.
        expect(M.lessThan(M.clone())).toBe(false);
        expect(M.lessThanOrEqual(M.clone())).toBe(true);
        expect(Matrix.fromArray(1, 2, [1, 2])
            .equals(Matrix.fromArray(1, 2, [1, 2]))).toBe(true);
    });

    it('clone and the free functions never alias their inputs', () => {
        check(matrix(2, 2), M => {
            const c = M.clone();
            c.set(0, 0, c.get(0, 0) + 1);
            expect(M.values).not.toEqual(c.values);
            for (const R of [addMatrix(M, M), subMatrix(M, M),
                negateMatrix(M), transpose(M), multiplyAB(M, M),
                mulMatrix(M, 2) as Matrix, divMatrix(M, 2)]) {
                expect(R.values).not.toBe(M.values);
            }
        });
    });
});
