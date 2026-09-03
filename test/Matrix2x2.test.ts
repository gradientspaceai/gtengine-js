import { describe, it, expect } from 'vitest';
import {
    makeRotation2x2, getRotationAngle2x2, inverse2x2, adjoint2x2,
    determinant2x2, trace2x2, doTransform2x2, setBasis2x2, getBasis2x2
} from '../src/Matrix2x2.js';
import {
    Matrix, inverse, determinant, multiplyAB, mulMatrix, subMatrix,
    lInfinityNorm
} from '../src/Matrix.js';
import { Vector } from '../src/Vector.js';

function expectMatrixClose(actual: Matrix, expected: Matrix,
    tolerance: number = 1e-12): void {
    expect(actual.numRows).toBe(expected.numRows);
    expect(actual.numCols).toBe(expected.numCols);
    for (let i = 0; i < expected.numElements; ++i) {
        expect(Math.abs(actual.values[i] - expected.values[i]))
            .toBeLessThanOrEqual(tolerance);
    }
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296 - 0.5;
    };
}

function randomMatrix2x2(rand: () => number): Matrix {
    const M = new Matrix(2, 2);
    for (let i = 0; i < 4; ++i) {
        M.values[i] = 4 * rand();
    }
    return M;
}

describe('Matrix2x2', () => {
    it('makeRotation2x2 produces {{c,-s},{s,c}}', () => {
        const R = new Matrix(2, 2);
        makeRotation2x2(Math.PI / 2, R);
        expectMatrixClose(R, Matrix.fromArray(2, 2, [0, -1, 1, 0]), 1e-15);

        makeRotation2x2(0, R);
        expectMatrixClose(R, Matrix.identity(2, 2));

        const angle = 0.3;
        makeRotation2x2(angle, R);
        expect(R.get(0, 0)).toBeCloseTo(Math.cos(angle), 15);
        expect(R.get(0, 1)).toBeCloseTo(-Math.sin(angle), 15);
        expect(R.get(1, 0)).toBeCloseTo(Math.sin(angle), 15);
        expect(R.get(1, 1)).toBeCloseTo(Math.cos(angle), 15);
    });

    it('rotates a vector counterclockwise', () => {
        const R = new Matrix(2, 2);
        makeRotation2x2(Math.PI / 2, R);
        const v = doTransform2x2(R, Vector.fromArray([1, 0]));
        expect(v.values[0]).toBeCloseTo(0, 14);
        expect(v.values[1]).toBeCloseTo(1, 14);
    });

    it('getRotationAngle2x2 inverts makeRotation2x2', () => {
        const R = new Matrix(2, 2);
        for (const angle of [-3.0, -1.0, -0.25, 0, 0.25, 1.0, 3.0]) {
            makeRotation2x2(angle, R);
            expect(getRotationAngle2x2(R)).toBeCloseTo(angle, 14);
        }
    });

    it('rotation matrices are orthonormal with determinant 1', () => {
        const R = new Matrix(2, 2);
        makeRotation2x2(0.7, R);
        expect(determinant2x2(R)).toBeCloseTo(1, 15);
        expect(trace2x2(R)).toBeCloseTo(2 * Math.cos(0.7), 15);
        const inv = inverse2x2(R).inverse;
        // The inverse of a rotation is its transpose.
        expect(inv.get(0, 0)).toBeCloseTo(R.get(0, 0), 14);
        expect(inv.get(0, 1)).toBeCloseTo(R.get(1, 0), 14);
        expect(inv.get(1, 0)).toBeCloseTo(R.get(0, 1), 14);
        expect(inv.get(1, 1)).toBeCloseTo(R.get(1, 1), 14);
    });

    it('determinant2x2 matches a known value and the generic determinant',
        () => {
            const M = Matrix.fromArray(2, 2, [1, 2, 3, 4]);
            expect(determinant2x2(M)).toBe(-2);

            const rand = makeRandom(12345);
            for (let trial = 0; trial < 50; ++trial) {
                const R = randomMatrix2x2(rand);
                expect(determinant2x2(R)).toBeCloseTo(determinant(R), 12);
            }
        });

    it('trace2x2 is the sum of the diagonal', () => {
        expect(trace2x2(Matrix.fromArray(2, 2, [1, 2, 3, 4]))).toBe(5);
        expect(trace2x2(Matrix.identity(2, 2))).toBe(2);
    });

    it('inverse2x2 matches a known value and the generic inverse', () => {
        const M = Matrix.fromArray(2, 2, [1, 2, 3, 4]);
        const result = inverse2x2(M);
        expect(result.invertible).toBe(true);
        expectMatrixClose(result.inverse,
            Matrix.fromArray(2, 2, [-2, 1, 1.5, -0.5]), 1e-15);

        const rand = makeRandom(999);
        for (let trial = 0; trial < 50; ++trial) {
            const R = randomMatrix2x2(rand);
            if (Math.abs(determinant2x2(R)) < 1e-3) {
                continue;
            }
            const fast = inverse2x2(R);
            const generic = inverse(R);
            expect(fast.invertible).toBe(generic.invertible);
            expectMatrixClose(fast.inverse, generic.inverse, 1e-9);
            expectMatrixClose(multiplyAB(R, fast.inverse),
                Matrix.identity(2, 2), 1e-12);
        }
    });

    it('inverse2x2 of a singular matrix is zero and reports false', () => {
        const M = Matrix.fromArray(2, 2, [1, 2, 2, 4]);
        const result = inverse2x2(M);
        expect(result.invertible).toBe(false);
        expectMatrixClose(result.inverse, new Matrix(2, 2));

        const zero = inverse2x2(new Matrix(2, 2));
        expect(zero.invertible).toBe(false);
        expectMatrixClose(zero.inverse, new Matrix(2, 2));
    });

    it('adjoint2x2 satisfies M*adj(M) = det(M)*I', () => {
        const rand = makeRandom(2024);
        for (let trial = 0; trial < 40; ++trial) {
            const M = randomMatrix2x2(rand);
            const adj = adjoint2x2(M);
            const det = determinant2x2(M);
            expectMatrixClose(multiplyAB(M, adj),
                mulMatrix(Matrix.identity(2, 2), det), 1e-12);
            expectMatrixClose(multiplyAB(adj, M),
                mulMatrix(Matrix.identity(2, 2), det), 1e-12);
        }
    });

    it('adjoint2x2 of a singular matrix is well defined', () => {
        // adj is defined even where the inverse is not.
        const M = Matrix.fromArray(2, 2, [1, 2, 2, 4]);
        expectMatrixClose(adjoint2x2(M), Matrix.fromArray(2, 2, [4, -2, -2, 1]));
        expect(lInfinityNorm(subMatrix(multiplyAB(M, adjoint2x2(M)),
            new Matrix(2, 2)))).toBe(0);
    });

    it('doTransform2x2 is M*V and A*B (GTE_USE_MAT_VEC)', () => {
        const A = Matrix.fromArray(2, 2, [1, 2, 3, 4]);
        const B = Matrix.fromArray(2, 2, [5, 6, 7, 8]);
        expectMatrixClose(doTransform2x2(A, B), multiplyAB(A, B));

        const V = Vector.fromArray([1, -1]);
        const w = doTransform2x2(A, V);
        expect(w.values).toEqual([-1, -1]);
    });

    it('setBasis2x2/getBasis2x2 access the columns', () => {
        const M = new Matrix(2, 2);
        setBasis2x2(M, 0, Vector.fromArray([1, 2]));
        setBasis2x2(M, 1, Vector.fromArray([3, 4]));
        expectMatrixClose(M, Matrix.fromArray(2, 2, [1, 3, 2, 4]));
        expect(getBasis2x2(M, 0).values).toEqual([1, 2]);
        expect(getBasis2x2(M, 1).values).toEqual([3, 4]);
    });

    it('rejects inputs of the wrong size', () => {
        expect(() => determinant2x2(new Matrix(3, 3))).toThrow();
        expect(() => trace2x2(new Matrix(2, 3))).toThrow();
        expect(() => inverse2x2(new Matrix(3, 3))).toThrow();
        expect(() => adjoint2x2(new Matrix(1, 1))).toThrow();
        expect(() => makeRotation2x2(1, new Matrix(3, 3))).toThrow();
        expect(() => getRotationAngle2x2(new Matrix(3, 3))).toThrow();
        expect(() => setBasis2x2(new Matrix(2, 2), 0, new Vector(3))).toThrow();
        expect(() => getBasis2x2(new Matrix(4, 4), 0)).toThrow();
        expect(() => doTransform2x2(new Matrix(2, 2), new Vector(3))).toThrow();
        expect(() => doTransform2x2(new Matrix(2, 2), new Matrix(3, 3)))
            .toThrow();
    });
});
