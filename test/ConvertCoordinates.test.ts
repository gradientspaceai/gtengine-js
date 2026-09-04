import { describe, it, expect } from 'vitest';
import { ConvertCoordinates } from '../src/ConvertCoordinates.js';
import {
    Matrix, multiplyAB, mulMatrix, determinant, transpose
} from '../src/Matrix.js';
import { Vector } from '../src/Vector.js';
import { determinant3x3 } from '../src/Matrix3x3.js';
import {
    check, expectVectorClose as expectVectorsClose, fc, invertibleMatrix,
    vector
} from './helpers/arbitraries.js';

function expectVectorClose(actual: Vector, expected: readonly number[],
    tolerance: number = 1e-12): void {
    expect(actual.size).toBe(expected.length);
    for (let i = 0; i < expected.length; ++i) {
        expect(Math.abs(actual.values[i] - expected[i]))
            .toBeLessThanOrEqual(tolerance);
    }
}

function expectMatrixClose(actual: Matrix, expected: Matrix,
    tolerance: number = 1e-12): void {
    expect(actual.numRows).toBe(expected.numRows);
    expect(actual.numCols).toBe(expected.numCols);
    for (let i = 0; i < expected.numElements; ++i) {
        expect(Math.abs(actual.values[i] - expected.values[i]))
            .toBeLessThanOrEqual(tolerance);
    }
}

// Build a matrix from its columns (the basis vectors of a coordinate
// system), as the upstream documentation does with SetCol.
function fromColumns(cols: readonly (readonly number[])[]): Matrix {
    const n = cols.length;
    const M = new Matrix(n, n);
    for (let c = 0; c < n; ++c) {
        M.setCol(c, Vector.fromArray(cols[c]));
    }
    return M;
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296 - 0.5;
    };
}

function randomMatrix(n: number, rand: () => number): Matrix {
    const M = new Matrix(n, n);
    for (let i = 0; i < M.numElements; ++i) {
        M.values[i] = 4 * rand();
    }
    return M;
}

describe('ConvertCoordinates', () => {
    it('starts as the identity conversion', () => {
        const convert = new ConvertCoordinates(3);
        expect(convert.getDimension()).toBe(3);
        expectMatrixClose(convert.getC(), Matrix.identity(3, 3));
        expectMatrixClose(convert.getInverseC(), Matrix.identity(3, 3));
        expect(convert.isVectorOnRightU()).toBe(true);
        expect(convert.isVectorOnRightV()).toBe(true);
        expect(convert.isRightHandedU()).toBe(true);
        expect(convert.isRightHandedV()).toBe(true);
        expectVectorClose(convert.uToV(Vector.fromArray([1, 2, 3])),
            [1, 2, 3]);
    });

    describe('the documented 3D linear change of basis', () => {
        const U = fromColumns([[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
        const V = fromColumns([[1, 0, 0], [0, 0, 1], [0, 1, 0]]);
        const convert = new ConvertCoordinates(3);
        const ok = convert.compute(U, true, V, true);

        it('reports the handedness of each system', () => {
            expect(ok).toBe(true);
            expect(convert.isRightHandedU()).toBe(true);
            expect(convert.isRightHandedV()).toBe(false);
            expect(convert.isVectorOnRightU()).toBe(true);
            expect(convert.isVectorOnRightV()).toBe(true);
        });

        it('computes C = U^{-1}*V and its inverse', () => {
            expectMatrixClose(convert.getC(), V);
            expectMatrixClose(multiplyAB(convert.getC(),
                convert.getInverseC()), Matrix.identity(3, 3));
        });

        it('converts the documented vectors', () => {
            const X = Vector.fromArray([1, 2, 3]);
            const Y = convert.uToV(X);
            expectVectorClose(Y, [1, 3, 2]);
            // Equal Cartesian coordinates.
            expectVectorClose(mulMatrix(U, X), mulMatrix(V, Y).values);

            const Y2 = Vector.fromArray([0, 1, 2]);
            const X2 = convert.vToU(Y2);
            expectVectorClose(X2, [0, 2, 1]);
            expectVectorClose(mulMatrix(U, X2), mulMatrix(V, Y2).values);
        });

        it('converts the documented transformation', () => {
            const c = 0.6, s = 0.8;
            const A = fromColumns([[c, s, 0], [-s, c, 0], [0, 0, 1]]);
            const B = convert.uToV(A);
            expectVectorClose(B.getCol(0), [c, 0, s]);
            expectVectorClose(B.getCol(1), [0, 1, 0]);
            expectVectorClose(B.getCol(2), [-s, 0, c]);

            // Applying the transformations in each system gives the same
            // Cartesian result (both systems are vector-on-right).
            const X = convert.vToU(Vector.fromArray([0, 1, 2]));
            const Y = convert.uToV(X);
            const X1 = mulMatrix(A, X);
            const Y1 = mulMatrix(B, Y);
            expectVectorClose(mulMatrix(U, X1), mulMatrix(V, Y1).values);

            // vToU inverts uToV for transformations.
            expectMatrixClose(convert.vToU(B), A);
        });
    });

    describe('the documented 3D affine change of basis', () => {
        const U = fromColumns([[-1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0],
            [1, 2, 3, 1]]);
        const V = fromColumns([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 1, 0],
            [4, 5, 6, 1]]);
        const convert = new ConvertCoordinates(4);
        const ok = convert.compute(U, true, V, false);

        it('reports the handedness of each system', () => {
            expect(ok).toBe(true);
            expect(convert.isRightHandedU()).toBe(false);
            expect(convert.isRightHandedV()).toBe(true);
            expect(convert.isVectorOnRightU()).toBe(true);
            expect(convert.isVectorOnRightV()).toBe(false);
        });

        it('converts the documented points', () => {
            const X = Vector.fromArray([-1, 4, -3, 1]);
            const Y = convert.uToV(X);
            expectVectorClose(Y, [0, 2, 1, 1]);
            expectVectorClose(mulMatrix(U, X), mulMatrix(V, Y).values);
            // Both represent the Cartesian point (2,5,7,1).
            expectVectorClose(mulMatrix(U, X), [2, 5, 7, 1]);

            const Y2 = Vector.fromArray([1, 2, 3, 1]);
            const X2 = convert.vToU(Y2);
            expectVectorClose(X2, [-1, 6, -4, 1]);
            expectVectorClose(mulMatrix(U, X2), mulMatrix(V, Y2).values);
        });

        it('converts the documented transformation', () => {
            const c = 0.6, s = 0.8;
            const A = fromColumns([[c, s, 0, 0], [-s, c, 0, 0], [0, 0, 1, 0],
                [0.3, 1, -2, 1]]);
            const B = convert.uToV(A);
            // Upstream's header comment lists these four tuples as
            // "B.GetCol(i)", but because the V system is vector-on-left the
            // output is transposed, so they are the ROWS of B. See the
            // "Upstream bug suspects" section of the port PR.
            expectVectorClose(B.getRow(0), [1, 0, 0, 0], 1e-12);
            expectVectorClose(B.getRow(1), [0, c, s, 0], 1e-12);
            expectVectorClose(B.getRow(2), [0, -s, c, 0], 1e-12);
            expectVectorClose(B.getRow(3), [2.0, -0.9, -2.6, 1], 1e-12);

            // U is vector-on-right and V is vector-on-left, so the
            // transformations are applied as X' = A*X and Y' = Y*B.
            const X = Vector.fromArray([-1, 4, -3, 1]);
            const Y = convert.uToV(X);
            const X1 = mulMatrix(A, X);
            const Y1 = mulMatrix(Y, B);
            expectVectorClose(mulMatrix(U, X1), mulMatrix(V, Y1).values);

            expectMatrixClose(convert.vToU(B), A, 1e-12);
        });
    });

    it('round-trips vectors and matrices for all convention combinations',
        () => {
            const rand = makeRandom(778899);
            for (const n of [2, 3, 4]) {
                for (const vorU of [true, false]) {
                    for (const vorV of [true, false]) {
                        let U = randomMatrix(n, rand);
                        let V = randomMatrix(n, rand);
                        while (Math.abs(determinant(U)) < 1e-2) {
                            U = randomMatrix(n, rand);
                        }
                        while (Math.abs(determinant(V)) < 1e-2) {
                            V = randomMatrix(n, rand);
                        }
                        const convert = new ConvertCoordinates(n);
                        expect(convert.compute(U, vorU, V, vorV)).toBe(true);

                        const X = new Vector(n);
                        for (let i = 0; i < n; ++i) {
                            X.values[i] = 4 * rand();
                        }
                        const Y = convert.uToV(X);
                        expectVectorClose(convert.vToU(Y), X.values, 1e-9);
                        // Equal Cartesian coordinates.
                        expectVectorClose(mulMatrix(U, X),
                            mulMatrix(V, Y).values, 1e-9);

                        const A = randomMatrix(n, rand);
                        const B = convert.uToV(A);
                        expectMatrixClose(convert.vToU(B), A, 1e-9);

                        // The transformation of X by A, expressed in the U
                        // convention, has the same Cartesian coordinates as
                        // the transformation of Y by B in the V convention.
                        const X1 = (vorU ? mulMatrix(A, X) : mulMatrix(X, A));
                        const Y1 = (vorV ? mulMatrix(B, Y) : mulMatrix(Y, B));
                        expectVectorClose(mulMatrix(U, X1),
                            mulMatrix(V, Y1).values, 1e-8);
                    }
                }
            }
        });

    it('implements the documented transformation tables', () => {
        const rand = makeRandom(5150);
        const n = 3;
        const U = fromColumns([[1, 0, 0], [0, 2, 0], [0, 0, 3]]);
        const V = fromColumns([[0, 1, 0], [1, 0, 0], [0, 0, 1]]);
        for (const vorU of [true, false]) {
            for (const vorV of [true, false]) {
                const convert = new ConvertCoordinates(n);
                expect(convert.compute(U, vorU, V, vorV)).toBe(true);
                const C = convert.getC();
                const invC = convert.getInverseC();
                const A = randomMatrix(n, rand);

                let expectedUToV = multiplyAB(multiplyAB(invC,
                    vorU ? A : transpose(A)), C);
                if (!vorV) {
                    expectedUToV = transpose(expectedUToV);
                }
                expectMatrixClose(convert.uToV(A), expectedUToV, 1e-12);

                const B = randomMatrix(n, rand);
                let expectedVToU = multiplyAB(multiplyAB(C,
                    vorV ? B : transpose(B)), invC);
                if (!vorU) {
                    expectedVToU = transpose(expectedVToU);
                }
                expectMatrixClose(convert.vToU(B), expectedVToU, 1e-12);
            }
        }
    });

    it('fails and resets when U or V is singular', () => {
        const convert = new ConvertCoordinates(3);
        const good = fromColumns([[1, 0, 0], [0, 0, 1], [0, 1, 0]]);
        expect(convert.compute(good, false, good, false)).toBe(true);
        expect(convert.isVectorOnRightU()).toBe(false);

        const singular = Matrix.fromArray(3, 3, [1, 2, 3, 4, 5, 6, 7, 8, 9]);
        expect(convert.compute(singular, true, good, true)).toBe(false);
        // The state is reset to the identity conversion.
        expectMatrixClose(convert.getC(), Matrix.identity(3, 3));
        expectMatrixClose(convert.getInverseC(), Matrix.identity(3, 3));
        expect(convert.isVectorOnRightU()).toBe(true);
        expect(convert.isVectorOnRightV()).toBe(true);
        expect(convert.isRightHandedU()).toBe(true);
        expect(convert.isRightHandedV()).toBe(true);

        expect(convert.compute(good, true, singular, true)).toBe(false);
        expectMatrixClose(convert.getC(), Matrix.identity(3, 3));
    });

    it('rejects mismatched sizes', () => {
        expect(() => new ConvertCoordinates(0)).toThrow();
        const convert = new ConvertCoordinates(3);
        expect(() => convert.compute(new Matrix(2, 2), true,
            Matrix.identity(3, 3), true)).toThrow();
        expect(() => convert.compute(Matrix.identity(3, 3), true,
            new Matrix(4, 4), true)).toThrow();
        expect(() => convert.uToV(new Vector(4))).toThrow();
        expect(() => convert.vToU(new Vector(2))).toThrow();
        expect(() => convert.uToV(new Matrix(2, 2))).toThrow();
        expect(() => convert.vToU(new Matrix(4, 4))).toThrow();
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V02): property-based re-checks against
// ConvertCoordinates.h.
// ---------------------------------------------------------------------------

// The Cartesian point represented by the coordinates X of the system whose
// basis vectors are the columns of the basis matrix.
function cartesian(basis: Matrix, X: Vector): Vector {
    return mulMatrix(basis, X) as Vector;
}

// Apply the transformation matrix A to the coordinate vector X under the
// stated multiplication convention.
function apply(A: Matrix, X: Vector, vectorOnRight: boolean): Vector {
    return (vectorOnRight
        ? mulMatrix(A, X) as Vector
        : mulMatrix(X, A) as Vector);
}

describe('ConvertCoordinates verification', () => {
    it('uToV and vToU preserve the Cartesian point and invert each other',
        () => {
            check(fc.tuple(fc.constantFrom(2, 3, 4), fc.boolean(), fc.boolean())
                .chain(([n, vorU, vorV]) => fc.tuple(fc.constant(n),
                    invertibleMatrix(n), invertibleMatrix(n), vector(n),
                    fc.constant(vorU), fc.constant(vorV))),
                ([n, U, V, X, vorU, vorV]) => {
                    const convert = new ConvertCoordinates(n);
                    expect(convert.compute(U, vorU, V, vorV)).toBe(true);

                    // The defining property: X and Y = C^{-1}X denote the
                    // same Cartesian point.
                    // invertibleMatrix allows relative determinants down to
                    // 1e-3, so C = V^{-1} U can have condition number ~1e3-1e4
                    // and the round trips lose up to ~1e-6 (observed 1.03e-6);
                    // 1e-5 keeps the check well below any structural error.
                    const Y = convert.uToV(X);
                    expectVectorsClose(cartesian(V, Y), cartesian(U, X), 1e-5,
                        1e-5);

                    // Round trip.
                    expectVectorsClose(convert.vToU(Y), X, 1e-5, 1e-5);
                    expectVectorsClose(convert.uToV(convert.vToU(X)), X, 1e-5,
                        1e-5);

                    // C and its inverse really are inverses of each other.
                    expectMatrixClose(
                        multiplyAB(convert.getC(), convert.getInverseC()),
                        Matrix.identity(n, n), 1e-5);
                }, 40);
        });

    it('uToV of a transformation acts on the V-coordinates as the original '
        + 'acts on the U-coordinates, for all four conventions', () => {
            check(fc.tuple(fc.constantFrom(2, 3), fc.boolean(), fc.boolean())
                .chain(([n, vorU, vorV]) => fc.tuple(fc.constant(n),
                    invertibleMatrix(n), invertibleMatrix(n),
                    invertibleMatrix(n), vector(n), fc.constant(vorU),
                    fc.constant(vorV))),
                ([n, U, V, A, X, vorU, vorV]) => {
                    const convert = new ConvertCoordinates(n);
                    expect(convert.compute(U, vorU, V, vorV)).toBe(true);
                    expect(convert.isVectorOnRightU()).toBe(vorU);
                    expect(convert.isVectorOnRightV()).toBe(vorV);

                    const B = convert.uToV(A);
                    const Y = convert.uToV(X);

                    // Transform in each system, then compare the Cartesian
                    // points. This is the whole contract of the matrix
                    // conversion table, including its transposes.
                    const XT = apply(A, X, vorU);
                    const YT = apply(B, Y, vorV);
                    expectVectorsClose(cartesian(V, YT), cartesian(U, XT), 1e-5,
                        1e-5);

                    // vToU undoes uToV on transformations too.
                    expectMatrixClose(convert.vToU(B), A, 1e-6);
                }, 40);
        });

    it('the handedness flags report the signs of the determinants', () => {
        check(fc.tuple(invertibleMatrix(3), invertibleMatrix(3)), ([U, V]) => {
            const convert = new ConvertCoordinates(3);
            expect(convert.compute(U, true, V, true)).toBe(true);
            expect(convert.isRightHandedU()).toBe(determinant(U) > 0);
            expect(convert.isRightHandedV()).toBe(determinant(V) > 0);

            // Swapping two columns of U reverses its handedness.
            const swapped = U.clone();
            swapped.setCol(0, U.getCol(1));
            swapped.setCol(1, U.getCol(0));
            expect(convert.compute(swapped, true, V, true)).toBe(true);
            expect(convert.isRightHandedU()).toBe(!(determinant(U) > 0));
        });
    });

    it('a singular U or V leaves the converter in its identity state', () => {
        check(fc.tuple(invertibleMatrix(3), vector(3), vector(3)),
            ([V, u, w]) => {
                // A rank-deficient matrix whose third basis vector is zero.
                // A zero column makes the elimination pivot exactly zero, so
                // the singularity is detected without depending on round-off
                // (duplicate columns need not cancel exactly in floating
                // point).
                const U = new Matrix(3, 3);
                U.setCol(0, u);
                U.setCol(1, w);
                expect(Math.abs(determinant3x3(U))).toBe(0);

                const convert = new ConvertCoordinates(3);
                // Seed the converter with a valid conversion first, so the
                // test observes the documented reset rather than the
                // constructor's state.
                expect(convert.compute(V, false, V, false)).toBe(true);
                expect(convert.compute(U, false, V, false)).toBe(false);
                expectMatrixClose(convert.getC(), Matrix.identity(3, 3), 0);
                expectMatrixClose(convert.getInverseC(), Matrix.identity(3, 3),
                    0);
                expect(convert.isVectorOnRightU()).toBe(true);
                expect(convert.isVectorOnRightV()).toBe(true);
                expect(convert.isRightHandedU()).toBe(true);
                expect(convert.isRightHandedV()).toBe(true);

                // A singular V is rejected as well.
                expect(convert.compute(V, false, U, false)).toBe(false);
            });
    });
});
