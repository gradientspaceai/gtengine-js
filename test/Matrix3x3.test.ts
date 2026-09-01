import { describe, it, expect } from 'vitest';
import {
    inverse3x3, adjoint3x3, determinant3x3, trace3x3, doTransform3x3,
    setBasis3x3, getBasis3x3
} from '../src/Matrix3x3';
import {
    Matrix, inverse, determinant, multiplyAB, mulMatrix, transpose
} from '../src/Matrix';
import { Vector, dot } from '../src/Vector';
import { cross } from '../src/Vector3';

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

function randomMatrix3x3(rand: () => number): Matrix {
    const M = new Matrix(3, 3);
    for (let i = 0; i < 9; ++i) {
        M.values[i] = 4 * rand();
    }
    return M;
}

// A rotation about a unit-length axis by 'angle', built with the Rodrigues
// formula, used as an independent construction of an orthogonal matrix.
function rodrigues(axis: Vector, angle: number): Matrix {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    const R = new Matrix(3, 3);
    const a = axis.values;
    for (let r = 0; r < 3; ++r) {
        for (let col = 0; col < 3; ++col) {
            R.set(r, col, (1 - c) * a[r] * a[col] + (r === col ? c : 0));
        }
    }
    R.set(0, 1, R.get(0, 1) - s * a[2]);
    R.set(0, 2, R.get(0, 2) + s * a[1]);
    R.set(1, 0, R.get(1, 0) + s * a[2]);
    R.set(1, 2, R.get(1, 2) - s * a[0]);
    R.set(2, 0, R.get(2, 0) - s * a[1]);
    R.set(2, 1, R.get(2, 1) + s * a[0]);
    return R;
}

describe('Matrix3x3', () => {
    it('determinant3x3 matches known values and the generic determinant',
        () => {
            expect(determinant3x3(Matrix.identity(3, 3))).toBe(1);
            expect(determinant3x3(new Matrix(3, 3))).toBe(0);
            // A matrix with linearly dependent rows.
            expect(determinant3x3(Matrix.fromArray(3, 3,
                [1, 2, 3, 4, 5, 6, 7, 8, 9]))).toBe(0);
            expect(determinant3x3(Matrix.fromArray(3, 3,
                [2, 0, 0, 0, 3, 0, 0, 0, 4]))).toBe(24);

            const rand = makeRandom(4242);
            for (let trial = 0; trial < 50; ++trial) {
                const M = randomMatrix3x3(rand);
                expect(determinant3x3(M)).toBeCloseTo(determinant(M), 11);
            }
        });

    it('determinant3x3 is the triple scalar product of the rows', () => {
        const rand = makeRandom(77);
        for (let trial = 0; trial < 30; ++trial) {
            const M = randomMatrix3x3(rand);
            const triple = dot(M.getRow(0), cross(M.getRow(1), M.getRow(2)));
            expect(determinant3x3(M)).toBeCloseTo(triple, 12);
        }
    });

    it('trace3x3 is the sum of the diagonal', () => {
        expect(trace3x3(Matrix.identity(3, 3))).toBe(3);
        expect(trace3x3(Matrix.fromArray(3, 3,
            [1, 2, 3, 4, 5, 6, 7, 8, 9]))).toBe(15);
    });

    it('inverse3x3 matches a known value and the generic inverse', () => {
        const D = Matrix.fromArray(3, 3, [2, 0, 0, 0, 4, 0, 0, 0, 5]);
        const dInv = inverse3x3(D);
        expect(dInv.invertible).toBe(true);
        expectMatrixClose(dInv.inverse,
            Matrix.fromArray(3, 3, [0.5, 0, 0, 0, 0.25, 0, 0, 0, 0.2]), 1e-16);

        const rand = makeRandom(31337);
        for (let trial = 0; trial < 50; ++trial) {
            const M = randomMatrix3x3(rand);
            if (Math.abs(determinant3x3(M)) < 1e-2) {
                continue;
            }
            const fast = inverse3x3(M);
            const generic = inverse(M);
            expect(fast.invertible).toBe(true);
            expect(generic.invertible).toBe(true);
            expectMatrixClose(fast.inverse, generic.inverse, 1e-8);
            expectMatrixClose(multiplyAB(M, fast.inverse),
                Matrix.identity(3, 3), 1e-11);
            expectMatrixClose(multiplyAB(fast.inverse, M),
                Matrix.identity(3, 3), 1e-11);
        }
    });

    it('inverse3x3 of a singular matrix is zero and reports false', () => {
        const result = inverse3x3(Matrix.fromArray(3, 3,
            [1, 2, 3, 4, 5, 6, 7, 8, 9]));
        expect(result.invertible).toBe(false);
        expectMatrixClose(result.inverse, new Matrix(3, 3));
    });

    it('the inverse of a rotation is its transpose', () => {
        const axis = Vector.fromArray([1 / 3, 2 / 3, 2 / 3]);
        const R = rodrigues(axis, 0.9);
        expect(determinant3x3(R)).toBeCloseTo(1, 13);
        expectMatrixClose(inverse3x3(R).inverse, transpose(R), 1e-14);
    });

    it('adjoint3x3 satisfies M*adj(M) = det(M)*I and equals det*inverse',
        () => {
            const rand = makeRandom(555);
            for (let trial = 0; trial < 40; ++trial) {
                const M = randomMatrix3x3(rand);
                const adj = adjoint3x3(M);
                const det = determinant3x3(M);
                const detI = mulMatrix(Matrix.identity(3, 3), det);
                expectMatrixClose(multiplyAB(M, adj), detI, 1e-11);
                expectMatrixClose(multiplyAB(adj, M), detI, 1e-11);
                if (Math.abs(det) > 1e-2) {
                    expectMatrixClose(mulMatrix(inverse3x3(M).inverse, det),
                        adj, 1e-11);
                }
            }
        });

    it('adjoint3x3 of a singular matrix is well defined', () => {
        const M = Matrix.fromArray(3, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9]);
        expectMatrixClose(adjoint3x3(M),
            Matrix.fromArray(3, 3, [-3, 6, -3, 6, -12, 6, -3, 6, -3]), 1e-13);
    });

    it('doTransform3x3 is M*V and A*B (GTE_USE_MAT_VEC)', () => {
        const A = Matrix.fromArray(3, 3, [1, 2, 3, 4, 5, 6, 7, 8, 10]);
        const B = Matrix.identity(3, 3);
        expectMatrixClose(doTransform3x3(A, B), A);
        const V = Vector.fromArray([1, 0, -1]);
        expect(doTransform3x3(A, V).values).toEqual([-2, -2, -3]);
    });

    it('setBasis3x3/getBasis3x3 access the columns', () => {
        const M = new Matrix(3, 3);
        setBasis3x3(M, 0, Vector.fromArray([1, 0, 0]));
        setBasis3x3(M, 1, Vector.fromArray([0, 0, 1]));
        setBasis3x3(M, 2, Vector.fromArray([0, 1, 0]));
        expect(getBasis3x3(M, 1).values).toEqual([0, 0, 1]);
        expect(determinant3x3(M)).toBe(-1);
    });

    it('rejects inputs of the wrong size', () => {
        expect(() => determinant3x3(new Matrix(2, 2))).toThrow();
        expect(() => trace3x3(new Matrix(4, 4))).toThrow();
        expect(() => inverse3x3(new Matrix(3, 4))).toThrow();
        expect(() => adjoint3x3(new Matrix(2, 3))).toThrow();
        expect(() => setBasis3x3(new Matrix(3, 3), 0, new Vector(2))).toThrow();
        expect(() => getBasis3x3(new Matrix(2, 2), 0)).toThrow();
        expect(() => doTransform3x3(new Matrix(3, 3), new Vector(4))).toThrow();
        expect(() => doTransform3x3(new Matrix(3, 3), new Matrix(2, 2)))
            .toThrow();
    });
});
