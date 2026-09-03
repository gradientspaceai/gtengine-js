import { describe, it, expect } from 'vitest';
import { check, matrix, fc } from './helpers/arbitraries.js';
import { GMatrix } from '../src/GMatrix.js';
import {
    Matrix, addMatrix, divMatrix, multiplyAB, determinant, inverse,
    transpose, mulMatrix
} from '../src/Matrix.js';
import { Vector } from '../src/Vector.js';

describe('GMatrix construction and sizing', () => {
    it('the default constructor produces the 0-by-0 matrix', () => {
        const M = new GMatrix();
        expect(M.numRows).toBe(0);
        expect(M.numCols).toBe(0);
        expect(M.numElements).toBe(0);
        expect(M.values).toEqual([]);
    });

    it('GMatrix(numRows, numCols) zero-fills', () => {
        const M = new GMatrix(2, 3);
        expect(M.getSize()).toEqual({ numRows: 2, numCols: 3 });
        expect(M.values).toEqual([0, 0, 0, 0, 0, 0]);
    });

    it('a nonpositive dimension produces the 0-by-0 matrix', () => {
        const M = new GMatrix(0, 5);
        expect(M.numRows).toBe(0);
        expect(M.numCols).toBe(0);
        expect(M.values).toEqual([]);
    });

    it('GMatrix(numRows, numCols, r, c) makes a Euclidean basis matrix', () => {
        expect(new GMatrix(2, 3, 1, 2).values).toEqual([0, 0, 0, 0, 0, 1]);
    });

    it('GMatrix(numRows, numCols, r, c) throws for an invalid index, unlike '
        + 'the fixed-size Matrix', () => {
        expect(() => new GMatrix(2, 2, 2, 0)).toThrow('Invalid index.');
        expect(() => new GMatrix(2, 2, 0, -1)).toThrow('Invalid index.');
    });

    it('is a Matrix (upstream API parity via subclassing)', () => {
        expect(new GMatrix(2, 2)).toBeInstanceOf(Matrix);
    });

    it('the static factories return GMatrix instances', () => {
        expect(GMatrix.zero(2, 2)).toBeInstanceOf(GMatrix);
        expect(GMatrix.unit(2, 2, 0, 1)).toBeInstanceOf(GMatrix);
        expect(GMatrix.identity(2, 3)).toBeInstanceOf(GMatrix);
        expect(GMatrix.fromArray(2, 2, [1, 2, 3, 4])).toBeInstanceOf(GMatrix);
        expect(GMatrix.identity(2, 3).values).toEqual([1, 0, 0, 0, 1, 0]);
        expect(GMatrix.unit(2, 2, 0, 1).values).toEqual([0, 1, 0, 0]);
    });

    it('fromMatrix copies a plain Matrix', () => {
        const M = Matrix.fromArray(2, 2, [1, 2, 3, 4]);
        const G = GMatrix.fromMatrix(M);
        expect(G).toBeInstanceOf(GMatrix);
        expect(G.values).toEqual([1, 2, 3, 4]);
        G.set(0, 0, 9);
        expect(M.get(0, 0)).toBe(1);
    });

    it('clone is a deep copy that is a GMatrix', () => {
        const M = GMatrix.fromArray(2, 2, [1, 2, 3, 4]);
        const N = M.clone();
        expect(N).toBeInstanceOf(GMatrix);
        N.set(1, 1, 7);
        expect(M.get(1, 1)).toBe(4);
    });
});

describe('GMatrix setSize', () => {
    it('grows the table, zeroing the new elements and keeping the leading '
        + 'ones (std::vector::resize semantics)', () => {
        const M = GMatrix.fromArray(2, 2, [1, 2, 3, 4]);
        M.setSize(2, 3);
        expect(M.numRows).toBe(2);
        expect(M.numCols).toBe(3);
        expect(M.values).toEqual([1, 2, 3, 4, 0, 0]);
        // The leading elements are preserved as a flat sequence, so their
        // (row, column) meaning changes with the number of columns.
        expect(M.get(0, 2)).toBe(3);
    });

    it('shrinks the table', () => {
        const M = GMatrix.fromArray(2, 3, [1, 2, 3, 4, 5, 6]);
        M.setSize(2, 2);
        expect(M.values).toEqual([1, 2, 3, 4]);
    });

    it('a nonpositive dimension clears the matrix', () => {
        const M = GMatrix.fromArray(2, 2, [1, 2, 3, 4]);
        M.setSize(0, 2);
        expect(M.numRows).toBe(0);
        expect(M.numCols).toBe(0);
        expect(M.values).toEqual([]);
        M.setSize(3, -1);
        expect(M.numElements).toBe(0);
    });

    it('resizing then filling behaves like a fresh matrix', () => {
        const M = new GMatrix();
        M.setSize(2, 2);
        M.makeIdentity();
        expect(M.values).toEqual([1, 0, 0, 1]);
    });
});

describe('GMatrix range-checked access', () => {
    const M = GMatrix.fromArray(2, 2, [1, 2, 3, 4]);

    it('get and set throw for an invalid index', () => {
        expect(M.get(1, 1)).toBe(4);
        expect(() => M.get(2, 0)).toThrow('Invalid index.');
        expect(() => M.get(0, -1)).toThrow('Invalid index.');
        expect(() => M.clone().set(0, 2, 1)).toThrow('Invalid index.');
    });

    it('makeUnit throws for an invalid index', () => {
        expect(() => M.clone().makeUnit(0, 5)).toThrow('Invalid index.');
        const N = M.clone();
        N.makeUnit(1, 0);
        expect(N.values).toEqual([0, 0, 1, 0]);
    });

    it('setRow and setCol validate the index and the vector size', () => {
        const N = new GMatrix(2, 3);
        N.setRow(1, Vector.fromArray([1, 2, 3]));
        expect(N.values).toEqual([0, 0, 0, 1, 2, 3]);
        N.setCol(0, Vector.fromArray([4, 5]));
        expect(N.values).toEqual([4, 0, 0, 5, 2, 3]);
        expect(() => N.setRow(2, Vector.fromArray([1, 2, 3])))
            .toThrow('Invalid index.');
        expect(() => N.setRow(0, Vector.fromArray([1, 2])))
            .toThrow('Mismatched sizes.');
        expect(() => N.setCol(3, Vector.fromArray([1, 2])))
            .toThrow('Invalid index.');
        expect(() => N.setCol(0, Vector.fromArray([1, 2, 3])))
            .toThrow('Mismatched sizes.');
    });

    it('getRow and getCol validate the index', () => {
        const N = GMatrix.fromArray(2, 3, [1, 2, 3, 4, 5, 6]);
        expect(N.getRow(0).values).toEqual([1, 2, 3]);
        expect(N.getCol(1).values).toEqual([2, 5]);
        expect(() => N.getRow(-1)).toThrow('Invalid index.');
        expect(() => N.getCol(3)).toThrow('Invalid index.');
    });
});

describe('GMatrix comparisons', () => {
    const A = GMatrix.fromArray(2, 2, [1, 2, 3, 4]);
    const B = GMatrix.fromArray(2, 2, [1, 2, 3, 5]);
    const C = GMatrix.fromArray(2, 3, [1, 2, 3, 4, 5, 6]);

    it('are lexicographic when the dimensions agree', () => {
        expect(A.equals(A.clone())).toBe(true);
        expect(A.notEquals(B)).toBe(true);
        expect(A.lessThan(B)).toBe(true);
        expect(A.lessThanOrEqual(B)).toBe(true);
        expect(B.greaterThan(A)).toBe(true);
        expect(B.greaterThanOrEqual(A)).toBe(true);
        expect(A.greaterThanOrEqual(A.clone())).toBe(true);
    });

    it('do not throw on mismatched dimensions; instead nothing holds but '
        + 'the inequality (upstream GMatrix semantics)', () => {
        expect(A.equals(C)).toBe(false);
        expect(A.notEquals(C)).toBe(true);
        expect(A.lessThan(C)).toBe(false);
        expect(C.lessThan(A)).toBe(false);
        expect(A.greaterThan(C)).toBe(false);
        expect(C.greaterThan(A)).toBe(false);
        expect(A.lessThanOrEqual(C)).toBe(false);
        expect(A.greaterThanOrEqual(C)).toBe(false);
    });
});

describe('GMatrix works with the shared Matrix free functions', () => {
    it('arithmetic and products', () => {
        const A = GMatrix.fromArray(2, 3, [1, 2, 3, 4, 5, 6]);
        const B = GMatrix.fromArray(3, 2, [7, 8, 9, 10, 11, 12]);
        expect(multiplyAB(A, B).values).toEqual([58, 64, 139, 154]);
        expect(addMatrix(A, A).values).toEqual([2, 4, 6, 8, 10, 12]);
        expect(transpose(A).values).toEqual([1, 4, 2, 5, 3, 6]);
        expect(mulMatrix(A, Vector.fromArray([1, 0, -1])).values)
            .toEqual([-2, -2]);
    });

    it('determinant and inverse', () => {
        const M = GMatrix.fromArray(3, 3, [6, 1, 1, 4, -2, 5, 2, 8, 7]);
        expect(determinant(M)).toBeCloseTo(-306, 10);
        const result = inverse(M);
        expect(result.invertible).toBe(true);
        const product = multiplyAB(M, result.inverse);
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 3; ++c) {
                expect(product.get(r, c)).toBeCloseTo(r === c ? 1 : 0, 10);
            }
        }
        expect(() => determinant(new GMatrix(2, 3)))
            .toThrow('Matrix must be square.');
    });

    it('division by zero returns the zero matrix (documented deviation from '
        + 'upstream GMatrix, which throws)', () => {
        const M = GMatrix.fromArray(2, 2, [1, 2, 3, 4]);
        expect(divMatrix(M, 0).values).toEqual([0, 0, 0, 0]);
        expect(divMatrix(M, 2).values).toEqual([0.5, 1, 1.5, 2]);
    });

    it('the shared free functions return plain Matrix objects that can be '
        + 'converted back', () => {
        const A = GMatrix.fromArray(2, 2, [1, 2, 3, 4]);
        const sum = addMatrix(A, A);
        expect(sum).toBeInstanceOf(Matrix);
        expect(sum).not.toBeInstanceOf(GMatrix);
        const G = GMatrix.fromMatrix(sum);
        G.setSize(2, 3);
        expect(G.numCols).toBe(3);
    });
});

// ---------------------------------------------------------------------------
// Verification wave: property-based checks against upstream GMatrix.h.
// ---------------------------------------------------------------------------

describe('GMatrix verification', () => {
    it('setSize follows std::vector::resize on the flat table', () => {
        // Upstream SetSize resizes mElements, which keeps the leading
        // elements as a flat sequence (their (row, column) meaning changes
        // with the number of columns) and value-initializes the new ones.
        check(fc.tuple(fc.integer({ min: 1, max: 3 }),
            fc.integer({ min: 1, max: 3 }), fc.integer({ min: -1, max: 4 }),
            fc.integer({ min: -1, max: 4 })).chain(([r0, c0, r1, c1]) =>
            fc.tuple(fc.constant(r1), fc.constant(c1), matrix(r0, c0))),
            ([numRows, numCols, M0]) => {
                const M = GMatrix.fromMatrix(M0);
                const old = [...M.values];
                M.setSize(numRows, numCols);
                if (numRows > 0 && numCols > 0) {
                    expect(M.numRows).toBe(numRows);
                    expect(M.numCols).toBe(numCols);
                    expect(M.numElements).toBe(numRows * numCols);
                    for (let i = 0; i < M.numElements; ++i) {
                        expect(M.values[i]).toBe(i < old.length ? old[i] : 0);
                    }
                } else {
                    // A nonpositive dimension clears the table to 0-by-0.
                    expect(M.numRows).toBe(0);
                    expect(M.numCols).toBe(0);
                    expect(M.values).toEqual([]);
                }
            });
    });

    it('the constructor rejects nonpositive dimensions the same way', () => {
        check(fc.tuple(fc.integer({ min: -2, max: 3 }),
            fc.integer({ min: -2, max: 3 })), ([numRows, numCols]) => {
            const M = new GMatrix(numRows, numCols);
            const valid = numRows > 0 && numCols > 0;
            expect(M.numRows).toBe(valid ? numRows : 0);
            expect(M.numCols).toBe(valid ? numCols : 0);
            expect(M.values.every(x => x === 0)).toBe(true);
        });
    });

    it('element access is range checked in both directions', () => {
        check(fc.tuple(matrix(2, 3), fc.integer({ min: -2, max: 3 }),
            fc.integer({ min: -2, max: 4 })), ([M0, r, c]) => {
            const M = GMatrix.fromMatrix(M0);
            const inRange = 0 <= r && r < 2 && 0 <= c && c < 3;
            if (inRange) {
                expect(M.get(r, c)).toBe(M.values[c + 3 * r]);
                M.set(r, c, 7);
                expect(M.get(r, c)).toBe(7);
                M.makeUnit(r, c);
                expect(M.values.reduce((a, b) => a + b, 0)).toBe(1);
                expect(M.get(r, c)).toBe(1);
            } else {
                expect(() => M.get(r, c)).toThrow('Invalid index.');
                expect(() => M.set(r, c, 7)).toThrow('Invalid index.');
                // Unlike the fixed-size Matrix, which silently zeroes.
                expect(() => M.makeUnit(r, c)).toThrow('Invalid index.');
            }
        });
    });

    it('getRow/getCol and setRow/setCol validate index and size', () => {
        check(fc.tuple(matrix(2, 3), fc.integer({ min: -1, max: 3 }),
            fc.integer({ min: 1, max: 4 })), ([M0, index, size]) => {
            const M = GMatrix.fromMatrix(M0);
            const v = Vector.fromArray(Array.from({ length: size },
                (_, i) => i + 1));
            if (0 <= index && index < 2) {
                expect(M.getRow(index).values)
                    .toEqual([0, 1, 2].map(c => M.get(index, c)));
                if (size === 3) {
                    M.setRow(index, v);
                    expect(M.getRow(index).values).toEqual(v.values);
                } else {
                    expect(() => M.setRow(index, v))
                        .toThrow('Mismatched sizes.');
                }
            } else {
                expect(() => M.getRow(index)).toThrow('Invalid index.');
                expect(() => M.setRow(index, v)).toThrow('Invalid index.');
            }
            if (0 <= index && index < 3) {
                expect(M.getCol(index).values)
                    .toEqual([0, 1].map(r => M.get(r, index)));
                if (size === 2) {
                    M.setCol(index, v);
                    expect(M.getCol(index).values).toEqual(v.values);
                } else {
                    expect(() => M.setCol(index, v))
                        .toThrow('Mismatched sizes.');
                }
            } else {
                expect(() => M.getCol(index)).toThrow('Invalid index.');
                expect(() => M.setCol(index, v)).toThrow('Invalid index.');
            }
        });
    });

    it('comparisons require equal dimensions (upstream GMatrix)', () => {
        // Upstream ANDs 'mNumRows == && mNumCols ==' into every relational
        // operator, so matrices of different dimensions satisfy none of them.
        // That is not a strict weak ordering; the quirk is preserved.
        check(fc.tuple(matrix(2, 2), matrix(2, 2), matrix(1, 4)),
            ([A0, B0, C0]) => {
                const A = GMatrix.fromMatrix(A0);
                const B = GMatrix.fromMatrix(B0);
                const C = GMatrix.fromMatrix(C0);
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
                // Same element count, different dimensions: nothing holds.
                expect(A.equals(C)).toBe(false);
                expect(A.notEquals(C)).toBe(true);
                expect(A.lessThan(C)).toBe(false);
                expect(C.lessThan(A)).toBe(false);
                expect(A.lessThanOrEqual(C)).toBe(false);
                expect(A.greaterThanOrEqual(C)).toBe(false);
            });
    });

    it('equals follows std::vector: a NaN element breaks self-equality', () => {
        // Regression: the port compared lexicographically, which treats NaN
        // as equivalent to itself; C++ 'NaN == NaN' is false, so upstream's
        // 'mElements == mat.mElements' is false.
        const M = GMatrix.fromArray(1, 2, [1, NaN]);
        expect(M.equals(M)).toBe(false);
        expect(M.notEquals(M)).toBe(true);
        expect(M.equals(M.clone())).toBe(false);
        // The ordering operators keep treating NaN as equivalent.
        expect(M.lessThan(M.clone())).toBe(false);
        expect(M.lessThanOrEqual(M.clone())).toBe(true);
        expect(GMatrix.fromArray(1, 2, [1, 2])
            .equals(GMatrix.fromArray(1, 2, [1, 2]))).toBe(true);
    });

    it('clone is deep, resizable and keeps the GMatrix type', () => {
        check(matrix(2, 3), M0 => {
            const M = GMatrix.fromMatrix(M0);
            const c = M.clone();
            expect(c).toBeInstanceOf(GMatrix);
            expect(c.values).toEqual(M.values);
            c.setSize(1, 2);
            expect(M.numRows).toBe(2);
            expect(M.numCols).toBe(3);
            expect(M.values).toEqual(M0.values);
        });
    });

    it('the shared Matrix free functions agree on GMatrix and Matrix', () => {
        // Upstream GMatrix.h re-implements each operation with the same
        // algorithm, so a GMatrix must give bit-identical results.
        check(fc.tuple(matrix(3, 2), matrix(2, 3)), ([A, B]) => {
            const GA = GMatrix.fromMatrix(A);
            const GB = GMatrix.fromMatrix(B);
            expect(multiplyAB(GA, GB).values).toEqual(multiplyAB(A, B).values);
            expect(addMatrix(GA, GA).values).toEqual(addMatrix(A, A).values);
            expect(transpose(GA).values).toEqual(transpose(A).values);
            expect(determinant(GMatrix.fromMatrix(multiplyAB(A, B))))
                .toBe(determinant(multiplyAB(A, B)));
            expect((mulMatrix(GA, Vector.fromArray([1, -2])) as Vector).values)
                .toEqual((mulMatrix(A, Vector.fromArray([1, -2])) as Vector)
                    .values);
        });
    });

    it('makeIdentity and makeZero work through the checked accessors', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 4 }),
            fc.integer({ min: 1, max: 4 })), ([numRows, numCols]) => {
            const M = GMatrix.identity(numRows, numCols);
            for (let r = 0; r < numRows; ++r) {
                for (let c = 0; c < numCols; ++c) {
                    expect(M.get(r, c)).toBe(r === c ? 1 : 0);
                }
            }
            M.makeZero();
            expect(M.values.every(x => x === 0)).toBe(true);
        });
    });
});
