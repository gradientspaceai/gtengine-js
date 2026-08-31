import { describe, it, expect } from 'vitest';
import { GMatrix } from '../src/GMatrix';
import {
    Matrix, addMatrix, divMatrix, multiplyAB, determinant, inverse,
    transpose, mulMatrix
} from '../src/Matrix';
import { Vector } from '../src/Vector';

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
