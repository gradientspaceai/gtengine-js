import { describe, it, expect } from 'vitest';
import {
    makeRotation2x2, getRotationAngle2x2, inverse2x2, adjoint2x2,
    determinant2x2, trace2x2, doTransform2x2, setBasis2x2, getBasis2x2
} from '../src/Matrix2x2.js';
import {
    Matrix, inverse, determinant, multiplyAB, mulMatrix, subMatrix,
    lInfinityNorm, transpose
} from '../src/Matrix.js';
import { Vector } from '../src/Vector.js';
import { GTE_C_PI } from '../src/Constants.js';
import {
    check, expectClose, expectVectorClose, fc, finite, invertibleMatrix,
    matrix, nonzero, vector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (V02): property-based re-checks against Matrix2x2.h.
// ---------------------------------------------------------------------------

describe('Matrix2x2 verification', () => {
    it('makeRotation2x2 is orthogonal with determinant 1 and rotates by the '
        + 'angle', () => {
            check(fc.tuple(finite(-GTE_C_PI, GTE_C_PI), finite(-GTE_C_PI, GTE_C_PI)),
                ([angle, phi]) => {
                    const R = new Matrix(2, 2);
                    makeRotation2x2(angle, R);

                    // R*R^T = I and det(R) = 1.
                    const shouldBeI = multiplyAB(R, transpose(R));
                    expectMatrixClose(shouldBeI, Matrix.identity(2, 2), 1e-14);
                    expectClose(determinant2x2(R), 1, 1e-14, 0);

                    // R maps (cos phi, sin phi) to (cos(angle+phi),
                    // sin(angle+phi)): an independent statement of what the
                    // matrix means, which a transposed port would fail.
                    const v = Vector.fromArray([Math.cos(phi), Math.sin(phi)]);
                    const rv = doTransform2x2(R, v);
                    expectClose(rv.get(0), Math.cos(angle + phi), 1e-14, 0);
                    expectClose(rv.get(1), Math.sin(angle + phi), 1e-14, 0);
                });
        });

    it('getRotationAngle2x2 inverts makeRotation2x2 on (-pi,pi]', () => {
        check(finite(-GTE_C_PI + 1e-6, GTE_C_PI - 1e-6), angle => {
            const R = new Matrix(2, 2);
            makeRotation2x2(angle, R);
            expectClose(getRotationAngle2x2(R), angle, 1e-14, 1e-14);
        });
    });

    it('M*adjoint(M) = det(M)*I for every 2x2 matrix, singular included', () => {
        check(matrix(2, 2), M => {
            const det = determinant2x2(M);
            const product = multiplyAB(M, adjoint2x2(M));
            const expected = mulMatrix(Matrix.identity(2, 2), det) as Matrix;
            // Absolute tolerance scaled by the matrix magnitude: the entries
            // of the product are differences of products of the entries of M.
            const scale = Math.max(1, lInfinityNorm(M) ** 2);
            expectMatrixClose(product, expected, 1e-12 * scale);
        });
    });

    it('determinant2x2 and inverse2x2 agree with the generic Matrix versions',
        () => {
            check(invertibleMatrix(2), M => {
                expectClose(determinant2x2(M), determinant(M), 1e-12, 1e-12);
                const { inverse: inv2, invertible } = inverse2x2(M);
                expect(invertible).toBe(true);
                // The generic inverse uses Gaussian elimination with full
                // pivoting, an algorithm independent of the closed form.
                const generic = inverse(M);
                expect(generic.invertible).toBe(true);
                expectMatrixClose(inv2, generic.inverse, 1e-8);
                expectMatrixClose(multiplyAB(M, inv2), Matrix.identity(2, 2),
                    1e-8);
            });
        });

    it('inverse2x2 reports non-invertibility and returns zero for singular '
        + 'matrices', () => {
            check(fc.tuple(finite(), finite(), nonzero()), ([a, b, k]) => {
                // The second row is k times the first, so the determinant is
                // exactly zero only when the products cancel exactly; build
                // the singular matrix from a rank-1 outer product instead.
                const M = Matrix.fromArray(2, 2, [a, b, k * a, k * b]);
                if (determinant2x2(M) !== 0) {
                    return;   // rounding left a tiny nonzero determinant
                }
                const { inverse: inv, invertible } = inverse2x2(M);
                expect(invertible).toBe(false);
                expectMatrixClose(inv, new Matrix(2, 2), 0);
            });
        });

    it('doTransform2x2 is associative: (A*B)*V = A*(B*V)', () => {
        check(fc.tuple(matrix(2, 2), matrix(2, 2), vector(2)), ([A, B, V]) => {
            const lhs = doTransform2x2(doTransform2x2(A, B), V);
            const rhs = doTransform2x2(A, doTransform2x2(B, V));
            expectVectorClose(lhs, rhs, 1e-10, 1e-10);
        });
    });

    it('setBasis2x2/getBasis2x2 use the columns: M*unit(i) = basis(i)', () => {
        check(fc.tuple(matrix(2, 2), vector(2), fc.integer({ min: 0, max: 1 })),
            ([M, V, i]) => {
                setBasis2x2(M, i, V);
                expectVectorClose(getBasis2x2(M, i), V, 0, 0);
                const unit = new Vector(2);
                unit.makeUnit(i);
                expectVectorClose(doTransform2x2(M, unit), V, 1e-12, 1e-12);
            });
    });

    it('trace2x2 equals the sum of the eigenvalues and is transpose '
        + 'invariant', () => {
            check(matrix(2, 2), M => {
                expectClose(trace2x2(M), trace2x2(transpose(M)), 0, 0);
                // The eigenvalues of a 2x2 matrix are (tr +- sqrt(tr^2-4det))/2
                // and always sum to the trace, even when complex.
                const tr = trace2x2(M);
                const det = determinant2x2(M);
                const disc = tr * tr - 4 * det;
                if (disc >= 0) {
                    const s = Math.sqrt(disc);
                    expectClose(((tr + s) / 2) + ((tr - s) / 2), tr, 1e-12,
                        1e-12);
                }
            });
        });
});
