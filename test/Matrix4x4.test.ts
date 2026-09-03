import { describe, it, expect } from 'vitest';
import {
    inverse4x4, adjoint4x4, determinant4x4, trace4x4, doTransform4x4,
    setBasis4x4, getBasis4x4, makeObliqueProjection4x4,
    makePerspectiveProjection4x4, makeReflection4x4
} from '../src/Matrix4x4.js';
import {
    Matrix, inverse, determinant, multiplyAB, mulMatrix, lInfinityNorm
} from '../src/Matrix.js';
import {
    Vector, add, dot, hproject, length, mul, normalize, sub
} from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, invertibleMatrix,
    matrix, unitVector, vector
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

function randomMatrix4x4(rand: () => number): Matrix {
    const M = new Matrix(4, 4);
    for (let i = 0; i < 16; ++i) {
        M.values[i] = 4 * rand();
    }
    return M;
}

// Apply the homogeneous transform M to the 3D point p and divide by w.
function applyHomogeneous(M: Matrix, p: Vector): Vector {
    const h = Vector.fromArray([p.values[0], p.values[1], p.values[2], 1]);
    const q = mulMatrix(M, h);
    const w = q.values[3];
    return Vector.fromArray([q.values[0] / w, q.values[1] / w,
        q.values[2] / w]);
}

describe('Matrix4x4', () => {
    it('determinant4x4 matches known values and the generic determinant',
        () => {
            expect(determinant4x4(Matrix.identity(4, 4))).toBe(1);
            expect(determinant4x4(new Matrix(4, 4))).toBe(0);
            expect(determinant4x4(Matrix.fromArray(4, 4,
                [2, 0, 0, 0, 0, 3, 0, 0, 0, 0, 5, 0, 0, 0, 0, 7]))).toBe(210);
            // Two identical rows.
            expect(determinant4x4(Matrix.fromArray(4, 4,
                [1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 9, 0, 1, 0, 1]))).toBe(0);

            const rand = makeRandom(8675309);
            for (let trial = 0; trial < 50; ++trial) {
                const M = randomMatrix4x4(rand);
                expect(determinant4x4(M)).toBeCloseTo(determinant(M), 10);
            }
        });

    it('trace4x4 is the sum of the diagonal', () => {
        expect(trace4x4(Matrix.identity(4, 4))).toBe(4);
        expect(trace4x4(Matrix.fromArray(4, 4,
            [1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 3, 0, 0, 0, 0, 4]))).toBe(10);
    });

    it('inverse4x4 matches a known value and the generic inverse', () => {
        const D = Matrix.fromArray(4, 4,
            [2, 0, 0, 0, 0, 4, 0, 0, 0, 0, 5, 0, 0, 0, 0, 8]);
        const dInv = inverse4x4(D);
        expect(dInv.invertible).toBe(true);
        expectMatrixClose(dInv.inverse, Matrix.fromArray(4, 4,
            [0.5, 0, 0, 0, 0, 0.25, 0, 0, 0, 0, 0.2, 0, 0, 0, 0, 0.125]),
            1e-16);

        const rand = makeRandom(20260830);
        let tested = 0;
        for (let trial = 0; trial < 60; ++trial) {
            const M = randomMatrix4x4(rand);
            if (Math.abs(determinant4x4(M)) < 1e-1) {
                continue;
            }
            ++tested;
            const fast = inverse4x4(M);
            const generic = inverse(M);
            expect(fast.invertible).toBe(true);
            expect(generic.invertible).toBe(true);
            expectMatrixClose(fast.inverse, generic.inverse, 1e-7);
            expectMatrixClose(multiplyAB(M, fast.inverse),
                Matrix.identity(4, 4), 1e-10);
            expectMatrixClose(multiplyAB(fast.inverse, M),
                Matrix.identity(4, 4), 1e-10);
        }
        expect(tested).toBeGreaterThan(20);
    });

    it('inverse4x4 of a singular matrix is zero and reports false', () => {
        const M = Matrix.fromArray(4, 4,
            [1, 2, 3, 4, 1, 2, 3, 4, 5, 6, 7, 9, 0, 1, 0, 1]);
        const result = inverse4x4(M);
        expect(result.invertible).toBe(false);
        expectMatrixClose(result.inverse, new Matrix(4, 4));
    });

    it('adjoint4x4 satisfies M*adj(M) = det(M)*I and equals det*inverse',
        () => {
            const rand = makeRandom(13579);
            for (let trial = 0; trial < 40; ++trial) {
                const M = randomMatrix4x4(rand);
                const adj = adjoint4x4(M);
                const det = determinant4x4(M);
                const detI = mulMatrix(Matrix.identity(4, 4), det);
                expectMatrixClose(multiplyAB(M, adj), detI, 1e-10);
                expectMatrixClose(multiplyAB(adj, M), detI, 1e-10);
                if (Math.abs(det) > 1e-1) {
                    expectMatrixClose(mulMatrix(inverse4x4(M).inverse, det),
                        adj, 1e-9);
                }
            }
        });

    it('doTransform4x4 is M*V and A*B (GTE_USE_MAT_VEC)', () => {
        const A = Matrix.fromArray(4, 4,
            [1, 0, 0, 3, 0, 1, 0, 4, 0, 0, 1, 5, 0, 0, 0, 1]);
        expectMatrixClose(doTransform4x4(A, Matrix.identity(4, 4)), A);
        // A is a translation; applying it to a point adds (3,4,5).
        const V = Vector.fromArray([1, 1, 1, 1]);
        expect(doTransform4x4(A, V).values).toEqual([4, 5, 6, 1]);
        // Composition of two translations.
        const B = Matrix.fromArray(4, 4,
            [1, 0, 0, -3, 0, 1, 0, -4, 0, 0, 1, -5, 0, 0, 0, 1]);
        expectMatrixClose(doTransform4x4(A, B), Matrix.identity(4, 4));
    });

    it('setBasis4x4/getBasis4x4 access the columns', () => {
        const M = new Matrix(4, 4);
        setBasis4x4(M, 0, Vector.fromArray([-1, 0, 0, 0]));
        setBasis4x4(M, 1, Vector.fromArray([0, 0, 1, 0]));
        setBasis4x4(M, 2, Vector.fromArray([0, -1, 0, 0]));
        setBasis4x4(M, 3, Vector.fromArray([1, 2, 3, 1]));
        expect(getBasis4x4(M, 3).values).toEqual([1, 2, 3, 1]);
        expect(determinant4x4(M)).toBe(-1);
    });

    describe('makeReflection4x4', () => {
        const normal = Vector.fromArray([0, 0, 1, 0]);
        const origin = Vector.fromArray([0, 0, 2, 1]);
        const M = makeReflection4x4(origin, normal);

        it('is the documented matrix I-2*N*N^T with 2*Dot(N,P)*N', () => {
            expectMatrixClose(M, Matrix.fromArray(4, 4,
                [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 4, 0, 0, 0, 1]), 1e-15);
        });

        it('fixes points on the plane and mirrors others', () => {
            const onPlane = Vector.fromArray([3, -7, 2]);
            const r = applyHomogeneous(M, onPlane);
            expect(r.values[0]).toBeCloseTo(3, 14);
            expect(r.values[1]).toBeCloseTo(-7, 14);
            expect(r.values[2]).toBeCloseTo(2, 14);

            const off = applyHomogeneous(M, Vector.fromArray([1, 2, 5]));
            expect(off.values).toEqual([1, 2, -1]);
        });

        it('is an involution for a general plane', () => {
            const n = Vector.fromArray([1, 2, -2, 0]);
            normalize(n);
            const p = Vector.fromArray([0.5, -1, 3, 1]);
            const R = makeReflection4x4(p, n);
            expectMatrixClose(multiplyAB(R, R), Matrix.identity(4, 4), 1e-13);
            // The reflection reverses orientation and its linear part
            // I-2*N*N^T is symmetric.
            expect(determinant4x4(R)).toBeCloseTo(-1, 12);
            for (let i = 0; i < 3; ++i) {
                for (let j = 0; j < 3; ++j) {
                    expect(R.get(i, j)).toBeCloseTo(R.get(j, i), 15);
                }
            }

            // The midpoint of a point and its reflection lies on the plane.
            const u = Vector.fromArray([4, -3, 1]);
            const v = applyHomogeneous(R, u);
            const mid = Vector.fromArray([
                0.5 * (u.values[0] + v.values[0]),
                0.5 * (u.values[1] + v.values[1]),
                0.5 * (u.values[2] + v.values[2]), 1]);
            expect(dot(n, sub(mid, p))).toBeCloseTo(0, 13);
        });
    });

    describe('makeObliqueProjection4x4', () => {
        it('projects along the direction onto the plane', () => {
            const n = Vector.fromArray([0.3, -0.5, 1, 0]);
            normalize(n);
            const p = Vector.fromArray([1, 1, 1, 1]);
            const d = Vector.fromArray([0.2, 0.1, -1, 0]);
            const M = makeObliqueProjection4x4(p, n, d);

            // M(3,3) > 0 whenever Dot(N,D) < 0.
            expect(dot(n, d)).toBeLessThan(0);
            expect(M.get(3, 3)).toBeGreaterThan(0);

            for (const u of [Vector.fromArray([2, 3, 4]),
                Vector.fromArray([-1, 0, 6]), Vector.fromArray([0, 0, 0])]) {
                const q = applyHomogeneous(M, u);
                const q4 = Vector.fromArray([q.values[0], q.values[1],
                    q.values[2], 1]);
                // The image lies on the plane.
                expect(dot(n, sub(q4, p))).toBeCloseTo(0, 12);
                // The displacement is parallel to the direction.
                const t = -dot(n, sub(Vector.fromArray(
                    [u.values[0], u.values[1], u.values[2], 1]), p)) / dot(n, d);
                expect(q.values[0]).toBeCloseTo(u.values[0] + t * d.values[0], 12);
                expect(q.values[1]).toBeCloseTo(u.values[1] + t * d.values[1], 12);
                expect(q.values[2]).toBeCloseTo(u.values[2] + t * d.values[2], 12);
            }
        });

        it('fixes points already on the plane', () => {
            const n = Vector.fromArray([0, 1, 0, 0]);
            const p = Vector.fromArray([0, 5, 0, 1]);
            const d = Vector.fromArray([1, -1, 0, 0]);
            const M = makeObliqueProjection4x4(p, n, d);
            const q = applyHomogeneous(M, Vector.fromArray([3, 5, -2]));
            expect(q.values[0]).toBeCloseTo(3, 13);
            expect(q.values[1]).toBeCloseTo(5, 13);
            expect(q.values[2]).toBeCloseTo(-2, 13);
        });

        it('is idempotent up to homogeneous scale', () => {
            const n = Vector.fromArray([0, 0, 1, 0]);
            const p = Vector.fromArray([0, 0, 0, 1]);
            const d = Vector.fromArray([0.5, 0.25, -1, 0]);
            const M = makeObliqueProjection4x4(p, n, d);
            const u = Vector.fromArray([2, -3, 7]);
            const once = applyHomogeneous(M, u);
            const twice = applyHomogeneous(M, once);
            for (let i = 0; i < 3; ++i) {
                expect(twice.values[i]).toBeCloseTo(once.values[i], 12);
            }
        });
    });

    describe('makePerspectiveProjection4x4', () => {
        it('projects points onto the plane along rays through the eye', () => {
            const n = Vector.fromArray([0, 0, 1, 0]);
            const p = Vector.fromArray([0, 0, 0, 1]);
            const e = Vector.fromArray([1, 2, 10, 1]);
            const M = makePerspectiveProjection4x4(p, n, e);

            for (const u of [Vector.fromArray([3, 4, 5]),
                Vector.fromArray([-2, 1, 2]), Vector.fromArray([0, 0, -4])]) {
                const q = applyHomogeneous(M, u);
                const q4 = Vector.fromArray([q.values[0], q.values[1],
                    q.values[2], 1]);
                // On the plane.
                expect(dot(n, sub(q4, p))).toBeCloseTo(0, 12);
                // Collinear with the eye and u: q = e + s*(u - e).
                const s = (q.values[2] - e.values[2])
                    / (u.values[2] - e.values[2]);
                expect(q.values[0]).toBeCloseTo(
                    e.values[0] + s * (u.values[0] - e.values[0]), 12);
                expect(q.values[1]).toBeCloseTo(
                    e.values[1] + s * (u.values[1] - e.values[1]), 12);
            }
        });

        it('fixes points already on the plane', () => {
            const n = Vector.fromArray([0, 0, 1, 0]);
            const p = Vector.fromArray([0, 0, 0, 1]);
            const e = Vector.fromArray([1, 2, 10, 1]);
            const M = makePerspectiveProjection4x4(p, n, e);
            const q = applyHomogeneous(M, Vector.fromArray([7, -3, 0]));
            expect(q.values[0]).toBeCloseTo(7, 12);
            expect(q.values[1]).toBeCloseTo(-3, 12);
            expect(q.values[2]).toBeCloseTo(0, 12);
        });

        it('has the documented block structure', () => {
            const n = Vector.fromArray([0, 0, 1, 0]);
            const p = Vector.fromArray([0, 0, 0, 1]);
            const e = Vector.fromArray([0, 0, 4, 1]);
            const M = makePerspectiveProjection4x4(p, n, e);
            // dotND = Dot(N, E-P) = 4, E*N^T has only the (2,2) entry 4.
            expectMatrixClose(M, Matrix.fromArray(4, 4, [
                4, 0, 0, 0,
                0, 4, 0, 0,
                0, 0, 0, 0,
                0, 0, -1, 4
            ]), 1e-14);
        });
    });

    it('rejects inputs of the wrong size', () => {
        expect(() => determinant4x4(new Matrix(3, 3))).toThrow();
        expect(() => trace4x4(new Matrix(2, 2))).toThrow();
        expect(() => inverse4x4(new Matrix(3, 3))).toThrow();
        expect(() => adjoint4x4(new Matrix(4, 3))).toThrow();
        expect(() => setBasis4x4(new Matrix(4, 4), 0, new Vector(3))).toThrow();
        expect(() => getBasis4x4(new Matrix(3, 3), 0)).toThrow();
        expect(() => doTransform4x4(new Matrix(4, 4), new Vector(3))).toThrow();
        expect(() => doTransform4x4(new Matrix(4, 4), new Matrix(3, 3)))
            .toThrow();
        expect(() => makeReflection4x4(new Vector(3), new Vector(4))).toThrow();
        expect(() => makeObliqueProjection4x4(new Vector(4), new Vector(4),
            new Vector(3))).toThrow();
        expect(() => makePerspectiveProjection4x4(new Vector(4),
            new Vector(2), new Vector(4))).toThrow();
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V02): property-based re-checks against Matrix4x4.h.
// ---------------------------------------------------------------------------

function vec4(x: number, y: number, z: number, w: number): Vector {
    return Vector.fromArray([x, y, z, w]);
}

describe('Matrix4x4 verification', () => {
    it('M*adjoint(M) = adjoint(M)*M = det(M)*I for every 4x4 matrix', () => {
        check(matrix(4, 4), M => {
            const det = determinant4x4(M);
            const adj = adjoint4x4(M);
            const expected = mulMatrix(Matrix.identity(4, 4), det) as Matrix;
            // The cofactors are degree-3 polynomials in the entries, so the
            // absolute error scales with the fourth power of the norm.
            const scale = Math.max(1, lInfinityNorm(M) ** 4);
            expectMatrixClose(multiplyAB(M, adj), expected, 1e-12 * scale);
            expectMatrixClose(multiplyAB(adj, M), expected, 1e-12 * scale);
        }, 100);
    });

    it('determinant4x4 and inverse4x4 agree with the generic Matrix versions',
        () => {
            check(invertibleMatrix(4), M => {
                const scale = Math.max(1, lInfinityNorm(M) ** 4);
                expectClose(determinant4x4(M), determinant(M), 1e-10 * scale,
                    1e-10);
                const { inverse: inv4, invertible } = inverse4x4(M);
                expect(invertible).toBe(true);
                const generic = inverse(M);
                expect(generic.invertible).toBe(true);
                expectMatrixClose(inv4, generic.inverse, 1e-6);
                expectMatrixClose(multiplyAB(M, inv4), Matrix.identity(4, 4),
                    1e-8);
                expectMatrixClose(multiplyAB(inv4, M), Matrix.identity(4, 4),
                    1e-8);
            }, 60);
        });

    it('inverse4x4 of a rank-deficient matrix reports non-invertibility',
        () => {
            check(fc.tuple(matrix(4, 4), fc.integer({ min: 0, max: 3 })),
                ([M, r]) => {
                    // Duplicate a row to force a zero determinant.
                    const other = (r + 1) % 4;
                    for (let c = 0; c < 4; ++c) {
                        M.set(r, c, M.get(other, c));
                    }
                    if (determinant4x4(M) !== 0) {
                        return;   // rounding left a tiny nonzero determinant
                    }
                    const { inverse: inv, invertible } = inverse4x4(M);
                    expect(invertible).toBe(false);
                    expectMatrixClose(inv, new Matrix(4, 4), 0);
                });
        });

    it('trace4x4 is transpose invariant; doTransform4x4 is associative and '
        + 'uses the columns as basis vectors', () => {
            check(fc.tuple(matrix(4, 4), matrix(4, 4), vector(4),
                fc.integer({ min: 0, max: 3 })), ([A, B, V, i]) => {
                    const At = new Matrix(4, 4);
                    for (let r = 0; r < 4; ++r) {
                        for (let c = 0; c < 4; ++c) {
                            At.set(c, r, A.get(r, c));
                        }
                    }
                    expectClose(trace4x4(A), trace4x4(At), 0, 0);

                    expectVectorClose(doTransform4x4(doTransform4x4(A, B), V),
                        doTransform4x4(A, doTransform4x4(B, V)), 1e-9, 1e-9);

                    setBasis4x4(A, i, V);
                    expectVectorClose(getBasis4x4(A, i), V, 0, 0);
                    const unit = new Vector(4);
                    unit.makeUnit(i);
                    expectVectorClose(doTransform4x4(A, unit), V, 1e-12, 1e-12);
                });
        });

    it('makeReflection4x4 is an involution fixing the plane and negating the '
        + 'normal offset', () => {
            check(fc.tuple(unitVector(3), vector(3, -5, 5), vector(3, -5, 5),
                finite(-5, 5)), ([n3, p3, u3, s]) => {
                    const N = vec4(n3.get(0), n3.get(1), n3.get(2), 0);
                    const P = vec4(p3.get(0), p3.get(1), p3.get(2), 1);
                    const M = makeReflection4x4(P, N);

                    // An involution: M*M = I, hence det(M) = -1 (a reflection
                    // reverses orientation).
                    expectMatrixClose(multiplyAB(M, M), Matrix.identity(4, 4),
                        1e-12);
                    expectClose(determinant4x4(M), -1, 1e-12, 1e-12);

                    // Points of the plane are fixed: X = P + (a tangent).
                    const tangent = sub(u3, mul(n3, dot(n3, u3)));
                    const X = vec4(p3.get(0) + tangent.get(0),
                        p3.get(1) + tangent.get(1),
                        p3.get(2) + tangent.get(2), 1);
                    expectVectorClose(doTransform4x4(M, X), X, 1e-10, 1e-10);

                    // Moving off the plane along the normal flips the offset.
                    const Y = add(X, mul(N, s));
                    const reflected = doTransform4x4(M, Y);
                    expectVectorClose(reflected, sub(X, mul(N, s)), 1e-10,
                        1e-10);
                });
        });

    it('makeObliqueProjection4x4 maps every point of a projector line to the '
        + 'same point of the plane', () => {
            check(fc.tuple(unitVector(3), vector(3, -5, 5), unitVector(3),
                vector(3, -5, 5), finite(-4, 4)),
                ([n3, p3, d3, u3, s]) => {
                    const dotND = dot(n3, d3);
                    if (Math.abs(dotND) < 0.2) {
                        return;   // the projection is ill-conditioned
                    }
                    const N = vec4(n3.get(0), n3.get(1), n3.get(2), 0);
                    const P = vec4(p3.get(0), p3.get(1), p3.get(2), 1);
                    const D = vec4(d3.get(0), d3.get(1), d3.get(2), 0);
                    const M = makeObliqueProjection4x4(P, N, D);

                    const U = vec4(u3.get(0), u3.get(1), u3.get(2), 1);
                    const proj = (X: Vector): Vector => {
                        const h = doTransform4x4(M, X);
                        return mul(h, 1 / h.get(3));
                    };

                    // The image lies on the plane Dot(N, X - P) = 0.
                    const X0 = proj(U);
                    expectClose(dot(n3, sub(hproject(X0), p3)), 0, 1e-9, 1e-9);
                    expectClose(X0.get(3), 1, 1e-12, 1e-12);

                    // U and U + s*D project to the same point, since the
                    // projection is along D.
                    const X1 = proj(add(U, mul(D, s)));
                    expectVectorClose(X1, X0, 1e-8, 1e-8);

                    // The projection is idempotent on the plane.
                    expectVectorClose(proj(X0), X0, 1e-8, 1e-8);
                });
        });

    it('makePerspectiveProjection4x4 projects along rays through the eye onto '
        + 'the plane', () => {
            check(fc.tuple(unitVector(3), finite(-4, 4), vector(3, -5, 5),
                vector(3, -5, 5)), ([n3, d, e3, u3]) => {
                    // The plane is Dot(N,X) = d with origin point P = d*N.
                    const p3 = mul(n3, d);
                    const N = vec4(n3.get(0), n3.get(1), n3.get(2), 0);
                    const P = vec4(p3.get(0), p3.get(1), p3.get(2), 1);
                    const E = vec4(e3.get(0), e3.get(1), e3.get(2), 1);
                    const V = sub(u3, e3);
                    const denom = dot(n3, V);
                    if (Math.abs(denom) < 0.25
                        || Math.abs(d - dot(n3, e3)) < 0.25) {
                        return;   // eye on the plane or ray parallel to it
                    }

                    const M = makePerspectiveProjection4x4(P, N, E);
                    const U = vec4(u3.get(0), u3.get(1), u3.get(2), 1);
                    const h = doTransform4x4(M, U);
                    if (Math.abs(h.get(3)) < 1e-6) {
                        return;
                    }
                    const X = hproject(mul(h, 1 / h.get(3)));

                    // The image is on the plane ...
                    expectClose(dot(n3, X), d, 1e-8, 1e-8);
                    // ... and on the line through E and U.
                    const cr = cross(sub(X, e3), V);
                    expectClose(length(cr), 0, 1e-7, 1e-7);
                });
        });
});
