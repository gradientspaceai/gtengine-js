import { describe, expect, it } from 'vitest';
import {
    LieSE2, LieSE3, LieSO2, LieSO3
} from '../src/LieGroupsAlgebras';
import { GTE_C_PI } from '../src/Constants';
import {
    Matrix, addMatrix, mulMatrix, multiplyAB, multiplyABT, subMatrix,
    transpose, lInfinityNorm
} from '../src/Matrix';
import { inverse3x3 } from '../src/Matrix3x3';
import { inverse4x4 } from '../src/Matrix4x4';
import { Vector, dot, sub } from '../src/Vector';

function maxDiff(A: Matrix, B: Matrix): number {
    return lInfinityNorm(subMatrix(A, B));
}

function vecMaxDiff(u: Vector, v: Vector): number {
    const d = sub(u, v);
    let m = 0;
    for (const value of d.values) {
        m = Math.max(m, Math.abs(value));
    }
    return m;
}

// A simple deterministic pseudorandom generator so the randomized checks are
// reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function isRotation(M: Matrix): number {
    const n = M.numRows;
    return maxDiff(multiplyABT(M, M), Matrix.identity(n, n));
}

// Build an SE(2) group element from an angle and a translation.
function makeSE2(angle: number, t0: number, t1: number): Matrix {
    const cs = Math.cos(angle);
    const sn = Math.sin(angle);
    return Matrix.fromArray(3, 3, [
        cs, -sn, t0,
        sn, cs, t1,
        0, 0, 1
    ]);
}

// Build an SE(3) group element from a rotation matrix and a translation.
function makeSE3(R: Matrix, t: Vector): Matrix {
    const M = Matrix.identity(4, 4);
    for (let r = 0; r < 3; ++r) {
        for (let c = 0; c < 3; ++c) {
            M.set(r, c, R.get(r, c));
        }
        M.set(r, 3, t.values[r]);
    }
    return M;
}

describe('LieSO2', () => {
    it('toGroup and toAlgebra are inverses', () => {
        const X = LieSO2.toGroup(0.75);
        expect(X.get(0, 0)).toBe(0);
        expect(X.get(1, 1)).toBe(0);
        expect(X.get(0, 1)).toBe(-0.75);
        expect(X.get(1, 0)).toBe(0.75);
        expect(LieSO2.toAlgebra(X)).toBe(0.75);
    });

    it('exp produces the standard 2D rotation matrix', () => {
        const angle = GTE_C_PI / 6;
        const Y = LieSO2.exp(angle);
        expect(Y.get(0, 0)).toBeCloseTo(Math.cos(angle), 15);
        expect(Y.get(0, 1)).toBeCloseTo(-Math.sin(angle), 15);
        expect(Y.get(1, 0)).toBeCloseTo(Math.sin(angle), 15);
        expect(Y.get(1, 1)).toBeCloseTo(Math.cos(angle), 15);
        expect(isRotation(Y)).toBeLessThan(1e-15);

        // exp(0) is the identity.
        expect(maxDiff(LieSO2.exp(0), Matrix.identity(2, 2))).toBe(0);
    });

    it('log(exp(x)) = x for generic and near-identity angles', () => {
        const angles = [
            0, 1e-14, 1e-8, 1e-4, 0.0625, 0.0626, 0.5, 1.5, 3.14,
            -1e-8, -0.5, -2.5, -3.14
        ];
        for (const angle of angles) {
            expect(LieSO2.log(LieSO2.exp(angle))).toBeCloseTo(angle, 12);
        }
    });

    it('exp(log(Y)) = Y for rotation matrices', () => {
        const rand = makeRandom(12345);
        let worst = 0;
        for (let i = 0; i < 200; ++i) {
            const angle = (2 * rand() - 1) * GTE_C_PI;
            const Y = LieSO2.exp(angle);
            worst = Math.max(worst, maxDiff(LieSO2.exp(LieSO2.log(Y)), Y));
        }
        expect(worst).toBeLessThan(1e-14);
    });

    it('the adjoint is the 1x1 identity', () => {
        expect(LieSO2.adjoint(LieSO2.exp(1.2))).toBe(1);
    });

    it('logM1M0Inv is the relative angle', () => {
        const M0 = LieSO2.exp(0.3);
        const M1 = LieSO2.exp(1.1);
        expect(LieSO2.logM1M0Inv(M0, M1)).toBeCloseTo(0.8, 14);
    });

    it('geodesicPath interpolates rotations', () => {
        const M0 = LieSO2.exp(0.3);
        const M1 = LieSO2.exp(1.1);
        expect(maxDiff(LieSO2.geodesicPath(0, M0, M1), M0)).toBeLessThan(1e-15);
        expect(maxDiff(LieSO2.geodesicPath(1, M0, M1), M1)).toBeLessThan(1e-14);
        expect(maxDiff(LieSO2.geodesicPath(0.5, M0, M1), LieSO2.exp(0.7)))
            .toBeLessThan(1e-14);

        // The precomputed-logarithm overload agrees with the other one.
        const logValue = LieSO2.logM1M0Inv(M0, M1);
        for (const t of [0, 0.25, 0.5, 0.75, 1]) {
            expect(maxDiff(LieSO2.geodesicPath(t, M0, logValue),
                LieSO2.geodesicPath(t, M0, M1))).toBeLessThan(1e-15);
        }
    });
});

describe('LieSE2', () => {
    it('toGroup and toAlgebra are inverses', () => {
        const x = Vector.fromArray([0.4, -1.5, 2.5]);
        const X = LieSE2.toGroup(x);
        expect(X.get(0, 1)).toBe(-0.4);
        expect(X.get(1, 0)).toBe(0.4);
        expect(X.get(0, 2)).toBe(-1.5);
        expect(X.get(1, 2)).toBe(2.5);
        expect(X.getRow(2).values).toEqual([0, 0, 0]);
        expect(vecMaxDiff(LieSE2.toAlgebra(X), x)).toBe(0);
    });

    it('exp of a pure translation is the translation', () => {
        const Y = LieSE2.exp(Vector.fromArray([0, 3, -4]));
        expect(maxDiff(Y, makeSE2(0, 3, -4))).toBeLessThan(1e-15);
    });

    it('exp produces a rigid motion with the expected rotation block', () => {
        const angle = 0.9;
        const Y = LieSE2.exp(Vector.fromArray([angle, 1, 2]));
        expect(Y.get(0, 0)).toBeCloseTo(Math.cos(angle), 15);
        expect(Y.get(1, 0)).toBeCloseTo(Math.sin(angle), 15);
        expect(Y.getRow(2).values).toEqual([0, 0, 1]);

        // The translation is V*u with V = [[a0,-a1],[a1,a0]], where
        // a0 = sin(t)/t and a1 = (1-cos(t))/t.
        const a0 = Math.sin(angle) / angle;
        const a1 = (1 - Math.cos(angle)) / angle;
        expect(Y.get(0, 2)).toBeCloseTo(a0 * 1 - a1 * 2, 14);
        expect(Y.get(1, 2)).toBeCloseTo(a1 * 1 + a0 * 2, 14);
    });

    it('log(exp(x)) = x across generic and near-identity inputs', () => {
        const rand = makeRandom(777);
        const angles = [
            0, 1e-12, 1e-6, 0.0624, 0.0625, 0.0626, 0.5, 2.0, 3.0,
            -1e-6, -0.5, -3.0
        ];
        let worst = 0;
        for (const angle of angles) {
            for (let i = 0; i < 10; ++i) {
                const x = Vector.fromArray([
                    angle,
                    4 * rand() - 2,
                    4 * rand() - 2
                ]);
                const y = LieSE2.log(LieSE2.exp(x));
                worst = Math.max(worst, vecMaxDiff(x, y));
            }
        }
        expect(worst).toBeLessThan(1e-11);
    });

    it('exp(log(M)) = M for rigid motions', () => {
        const rand = makeRandom(2468);
        let worst = 0;
        for (let i = 0; i < 200; ++i) {
            const angle = (2 * rand() - 1) * 3.14;
            const M = makeSE2(angle, 6 * rand() - 3, 6 * rand() - 3);
            worst = Math.max(worst, maxDiff(LieSE2.exp(LieSE2.log(M)), M));
        }
        expect(worst).toBeLessThan(1e-12);
    });

    it('the adjoint satisfies L(A(M)*x) = M*L(x)*Inverse(M)', () => {
        const rand = makeRandom(99);
        let worst = 0;
        for (let i = 0; i < 100; ++i) {
            const M = makeSE2((2 * rand() - 1) * 3, 4 * rand() - 2,
                4 * rand() - 2);
            const x = Vector.fromArray([
                2 * rand() - 1, 4 * rand() - 2, 4 * rand() - 2
            ]);
            const A = LieSE2.adjoint(M);
            const lhs = LieSE2.toGroup(mulMatrix(A, x));
            const rhs = multiplyAB(multiplyAB(M, LieSE2.toGroup(x)),
                inverse3x3(M).inverse);
            worst = Math.max(worst, maxDiff(lhs, rhs));
        }
        expect(worst).toBeLessThan(1e-12);

        // The adjoint of the identity is the identity.
        expect(maxDiff(LieSE2.adjoint(Matrix.identity(3, 3)),
            Matrix.identity(3, 3))).toBe(0);
    });

    it('geodesicPath has the expected endpoints and midpoint', () => {
        const M0 = makeSE2(0.2, 1, -1);
        const M1 = makeSE2(1.4, -2, 3);
        expect(maxDiff(LieSE2.geodesicPath(0, M0, M1), M0)).toBeLessThan(1e-13);
        expect(maxDiff(LieSE2.geodesicPath(1, M0, M1), M1)).toBeLessThan(1e-13);

        // The midpoint is a rigid motion whose rotation is the half angle of
        // the relative rotation applied to M0.
        const mid = LieSE2.geodesicPath(0.5, M0, M1);
        const midAngle = Math.atan2(mid.get(1, 0), mid.get(0, 0));
        expect(midAngle).toBeCloseTo(0.2 + 0.6, 12);

        const logValue = LieSE2.logM1M0Inv(M0, M1);
        for (const t of [0, 0.3, 0.5, 1]) {
            expect(maxDiff(LieSE2.geodesicPath(t, M0, logValue),
                LieSE2.geodesicPath(t, M0, M1))).toBeLessThan(1e-14);
        }
    });
});

describe('LieSO3', () => {
    it('toGroup builds the skew-symmetric matrix and toAlgebra inverts it',
        () => {
            const x = Vector.fromArray([1, 2, 3]);
            const X = LieSO3.toGroup(x);
            expect(maxDiff(X, Matrix.fromArray(3, 3, [
                0, -3, 2,
                3, 0, -1,
                -2, 1, 0
            ]))).toBe(0);
            expect(maxDiff(addMatrix(X, transpose(X)), Matrix.zero(3, 3)))
                .toBe(0);
            expect(vecMaxDiff(LieSO3.toAlgebra(X), x)).toBe(0);

            // L(x)*v = Cross(x, v).
            const v = Vector.fromArray([-4, 5, 6]);
            const Xv = mulMatrix(X, v);
            expect(Xv.values[0]).toBeCloseTo(2 * 6 - 3 * 5, 12);
            expect(Xv.values[1]).toBeCloseTo(3 * -4 - 1 * 6, 12);
            expect(Xv.values[2]).toBeCloseTo(1 * 5 - 2 * -4, 12);
        });

    it('exp gives the known rotation about the z axis', () => {
        const Y = LieSO3.exp(Vector.fromArray([0, 0, GTE_C_PI / 2]));
        expect(maxDiff(Y, Matrix.fromArray(3, 3, [
            0, -1, 0,
            1, 0, 0,
            0, 0, 1
        ]))).toBeLessThan(1e-15);

        // exp of the zero algebra element is the identity.
        expect(maxDiff(LieSO3.exp(Vector.zero(3)), Matrix.identity(3, 3)))
            .toBe(0);
        expect(vecMaxDiff(LieSO3.log(Matrix.identity(3, 3)), Vector.zero(3)))
            .toBe(0);
    });

    it('exp produces rotation matrices with determinant 1', () => {
        const rand = makeRandom(31415);
        let worstOrtho = 0;
        for (let i = 0; i < 200; ++i) {
            const x = Vector.fromArray([
                6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3
            ]);
            const Y = LieSO3.exp(x);
            worstOrtho = Math.max(worstOrtho, isRotation(Y));

            // The rotation axis is fixed by the rotation.
            expect(vecMaxDiff(mulMatrix(Y, x), x)).toBeLessThan(1e-13);
        }
        expect(worstOrtho).toBeLessThan(1e-14);
    });

    it('log(exp(x)) = x for generic and near-identity inputs', () => {
        const rand = makeRandom(5150);
        const angles = [
            1e-12, 1e-7, 1e-3, 0.0624, 0.0625, 0.0626, 0.5, 1.5, 2.5, 3.0,
            3.14
        ];
        let worst = 0;
        for (const angle of angles) {
            for (let i = 0; i < 10; ++i) {
                const axis = Vector.fromArray([
                    2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1
                ]);
                const len = Math.sqrt(dot(axis, axis));
                const x = Vector.fromArray(axis.values.map(
                    (v) => v * angle / len));
                const y = LieSO3.log(LieSO3.exp(x));
                worst = Math.max(worst, vecMaxDiff(x, y));
            }
        }
        expect(worst).toBeLessThan(1e-9);

        // Also check log(0) and the exact zero-angle path.
        expect(vecMaxDiff(LieSO3.log(LieSO3.exp(Vector.zero(3))),
            Vector.zero(3))).toBe(0);
    });

    it('log handles the angle = pi branch for each dominant diagonal', () => {
        // Rotations by pi about a unit axis n are R = 2*n*n^T - I. The
        // matrices below are exactly representable and have trace exactly
        // -1, so LieSO3.log takes its arg = -1 branch. The three cases
        // exercise the r00, r11 and r22 selections of the maximum diagonal
        // term.
        const cases: Array<{ R: Matrix, axis: Vector }> = [
            {   // r00 is the maximum diagonal term
                R: Matrix.fromArray(3, 3, [1, 0, 0, 0, -1, 0, 0, 0, -1]),
                axis: Vector.fromArray([1, 0, 0])
            },
            {   // r11 is the maximum diagonal term
                R: Matrix.fromArray(3, 3, [-1, 0, 0, 0, 1, 0, 0, 0, -1]),
                axis: Vector.fromArray([0, 1, 0])
            },
            {   // r22 is the maximum diagonal term
                R: Matrix.fromArray(3, 3, [-1, 0, 0, 0, -1, 0, 0, 0, 1]),
                axis: Vector.fromArray([0, 0, 1])
            },
            {   // r00 selected, axis (1,1,0)/sqrt(2)
                R: Matrix.fromArray(3, 3, [0, 1, 0, 1, 0, 0, 0, 0, -1]),
                axis: Vector.fromArray([1, 1, 0])
            },
            {   // r11 selected, axis (0,1,1)/sqrt(2)
                R: Matrix.fromArray(3, 3, [-1, 0, 0, 0, 0, 1, 0, 1, 0]),
                axis: Vector.fromArray([0, 1, 1])
            },
            {   // r00 selected, axis (1,0,1)/sqrt(2)
                R: Matrix.fromArray(3, 3, [0, 0, 1, 0, -1, 0, 1, 0, 0]),
                axis: Vector.fromArray([1, 0, 1])
            }
        ];

        for (const { R, axis } of cases) {
            expect(trace(R)).toBe(-1);
            expect(isRotation(R)).toBeLessThan(1e-15);

            const y = LieSO3.log(R);

            // The recovered element has length pi. Upstream scales the
            // normalized row by pi/sqrt(2), which would produce length
            // 2.2214...; the port uses pi (see the comment in LieSO3.log).
            expect(Math.sqrt(dot(y, y))).toBeCloseTo(GTE_C_PI, 12);

            // It is parallel to the rotation axis, up to sign, since x and
            // -x produce the same rotation matrix.
            const len = Math.sqrt(dot(axis, axis));
            const unit = Vector.fromArray(axis.values.map((v) => v / len));
            const sign = dot(y, unit) >= 0 ? 1 : -1;
            expect(vecMaxDiff(y, Vector.fromArray(
                unit.values.map((v) => sign * v * GTE_C_PI))))
                .toBeLessThan(1e-14);

            // The round trip reproduces the rotation matrix.
            expect(maxDiff(LieSO3.exp(y), R)).toBeLessThan(1e-14);
        }
    });

    it('log is well behaved just below the angle = pi knife edge', () => {
        // At an angle of exactly pi, Y - Y^T is zero and the arg > -1 branch
        // is ill conditioned; upstream (and this port) rely on the trace
        // being exactly -1 to reach the dedicated branch. Just below pi the
        // generic branch recovers the algebra element accurately.
        for (const angle of [3.0, 3.1, 3.14, 3.1415]) {
            const axis = Vector.fromArray([1, -2, 3]);
            const len = Math.sqrt(dot(axis, axis));
            const x = Vector.fromArray(axis.values.map(
                (v) => v * angle / len));
            const y = LieSO3.log(LieSO3.exp(x));
            expect(vecMaxDiff(x, y)).toBeLessThan(1e-8);
        }
    });

    it('the adjoint is the group element and satisfies the conjugation '
        + 'identity', () => {
            const rand = makeRandom(60606);
            let worst = 0;
            for (let i = 0; i < 100; ++i) {
                const M = LieSO3.exp(Vector.fromArray([
                    4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2
                ]));
                const x = Vector.fromArray([
                    2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1
                ]);
                const A = LieSO3.adjoint(M);
                expect(maxDiff(A, M)).toBe(0);

                // exp(L(M*x)) = M * exp(L(x)) * M^T.
                const lhs = LieSO3.exp(mulMatrix(A, x));
                const rhs = multiplyABT(multiplyAB(M, LieSO3.exp(x)), M);
                worst = Math.max(worst, maxDiff(lhs, rhs));
            }
            expect(worst).toBeLessThan(1e-13);
        });

    it('geodesicPath has the expected endpoints and midpoint', () => {
        const axis = Vector.fromArray([0, 0, 1]);
        const M0 = LieSO3.exp(Vector.fromArray([0, 0, 0.4]));
        const M1 = LieSO3.exp(Vector.fromArray([0, 0, 1.6]));
        expect(maxDiff(LieSO3.geodesicPath(0, M0, M1), M0)).toBeLessThan(1e-14);
        expect(maxDiff(LieSO3.geodesicPath(1, M0, M1), M1)).toBeLessThan(1e-13);
        expect(maxDiff(LieSO3.geodesicPath(0.5, M0, M1),
            LieSO3.exp(Vector.fromArray([0, 0, 1.0])))).toBeLessThan(1e-13);
        expect(axis.size).toBe(3);

        // A general pair of rotations: the path stays in SO(3) and the two
        // overloads agree.
        const N0 = LieSO3.exp(Vector.fromArray([0.3, -0.7, 1.1]));
        const N1 = LieSO3.exp(Vector.fromArray([-1.2, 0.5, 0.2]));
        const logValue = LieSO3.logM1M0Inv(N0, N1);
        for (const t of [0, 0.25, 0.5, 0.75, 1]) {
            const P = LieSO3.geodesicPath(t, N0, N1);
            expect(isRotation(P)).toBeLessThan(1e-13);
            expect(maxDiff(LieSO3.geodesicPath(t, N0, logValue), P))
                .toBeLessThan(1e-14);
        }
        expect(maxDiff(LieSO3.geodesicPath(1, N0, N1), N1)).toBeLessThan(1e-13);
    });
});

describe('LieSE3', () => {
    it('toGroup and toAlgebra are inverses', () => {
        const x = Vector.fromArray([1, 2, 3, 4, 5, 6]);
        const X = LieSE3.toGroup(x);
        expect(maxDiff(X, Matrix.fromArray(4, 4, [
            0, -3, 2, 4,
            3, 0, -1, 5,
            -2, 1, 0, 6,
            0, 0, 0, 0
        ]))).toBe(0);
        expect(vecMaxDiff(LieSE3.toAlgebra(X), x)).toBe(0);
    });

    it('exp of a pure translation is the translation', () => {
        const Y = LieSE3.exp(Vector.fromArray([0, 0, 0, 7, -8, 9]));
        expect(maxDiff(Y, makeSE3(Matrix.identity(3, 3),
            Vector.fromArray([7, -8, 9])))).toBeLessThan(1e-15);
    });

    it('the rotation block of exp agrees with LieSO3.exp', () => {
        const s = Vector.fromArray([0.3, -0.7, 1.1]);
        const Y = LieSE3.exp(Vector.fromArray([
            s.values[0], s.values[1], s.values[2], 1, -2, 3
        ]));
        const R = LieSO3.exp(s);
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 3; ++c) {
                expect(Y.get(r, c)).toBeCloseTo(R.get(r, c), 14);
            }
        }
        expect(Y.getRow(3).values).toEqual([0, 0, 0, 1]);
    });

    it('log(exp(x)) = x across generic and near-identity inputs', () => {
        const rand = makeRandom(9091);
        const angles = [
            0, 1e-12, 1e-6, 0.0624, 0.0625, 0.0626, 0.5, 1.5, 2.5, 3.0
        ];
        let worst = 0;
        for (const angle of angles) {
            for (let i = 0; i < 10; ++i) {
                const axis = [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1];
                const len = Math.hypot(axis[0], axis[1], axis[2]);
                const x = Vector.fromArray([
                    axis[0] * angle / len,
                    axis[1] * angle / len,
                    axis[2] * angle / len,
                    4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2
                ]);
                worst = Math.max(worst,
                    vecMaxDiff(LieSE3.log(LieSE3.exp(x)), x));
            }
        }
        expect(worst).toBeLessThan(1e-9);
    });

    it('exp(log(M)) = M for rigid motions', () => {
        const rand = makeRandom(4242);
        let worst = 0;
        for (let i = 0; i < 100; ++i) {
            const R = LieSO3.exp(Vector.fromArray([
                4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2
            ]));
            const M = makeSE3(R, Vector.fromArray([
                6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3
            ]));
            worst = Math.max(worst, maxDiff(LieSE3.exp(LieSE3.log(M)), M));
        }
        expect(worst).toBeLessThan(1e-11);
    });

    it('the adjoint satisfies L(A(M)*x) = M*L(x)*Inverse(M)', () => {
        const rand = makeRandom(1717);
        let worst = 0;
        for (let i = 0; i < 100; ++i) {
            const R = LieSO3.exp(Vector.fromArray([
                4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2
            ]));
            const M = makeSE3(R, Vector.fromArray([
                4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2
            ]));
            const x = Vector.fromArray([
                2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1,
                2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1
            ]);
            const A = LieSE3.adjoint(M);
            expect(A.numRows).toBe(6);
            expect(A.numCols).toBe(6);
            const lhs = LieSE3.toGroup(mulMatrix(A, x));
            const rhs = multiplyAB(multiplyAB(M, LieSE3.toGroup(x)),
                inverse4x4(M).inverse);
            worst = Math.max(worst, maxDiff(lhs, rhs));
        }
        expect(worst).toBeLessThan(1e-12);

        // The adjoint of the identity is the 6x6 identity.
        expect(maxDiff(LieSE3.adjoint(Matrix.identity(4, 4)),
            Matrix.identity(6, 6))).toBe(0);
    });

    it('geodesicPath has the expected endpoints and stays rigid', () => {
        const M0 = makeSE3(LieSO3.exp(Vector.fromArray([0.1, 0.2, 0.3])),
            Vector.fromArray([1, 2, 3]));
        const M1 = makeSE3(LieSO3.exp(Vector.fromArray([-0.8, 0.4, 1.2])),
            Vector.fromArray([-2, 0.5, 4]));
        expect(maxDiff(LieSE3.geodesicPath(0, M0, M1), M0)).toBeLessThan(1e-12);
        expect(maxDiff(LieSE3.geodesicPath(1, M0, M1), M1)).toBeLessThan(1e-11);

        const logValue = LieSE3.logM1M0Inv(M0, M1);
        for (const t of [0, 0.25, 0.5, 0.75, 1]) {
            const P = LieSE3.geodesicPath(t, M0, logValue);
            expect(maxDiff(P, LieSE3.geodesicPath(t, M0, M1)))
                .toBeLessThan(1e-13);
            expect(P.getRow(3).values[3]).toBeCloseTo(1, 14);

            // The upper 3x3 block is a rotation matrix.
            const R = Matrix.fromArray(3, 3, [
                P.get(0, 0), P.get(0, 1), P.get(0, 2),
                P.get(1, 0), P.get(1, 1), P.get(1, 2),
                P.get(2, 0), P.get(2, 1), P.get(2, 2)
            ]);
            expect(isRotation(R)).toBeLessThan(1e-12);
        }
    });
});

describe('LieGroupsAlgebras degenerate inputs', () => {
    it('rejects mismatched dimensions', () => {
        expect(() => LieSO2.log(Matrix.identity(3, 3))).toThrow();
        expect(() => LieSE2.exp(Vector.zero(2))).toThrow();
        expect(() => LieSO3.exp(Vector.zero(4))).toThrow();
        expect(() => LieSE3.exp(Vector.zero(3))).toThrow();
        expect(() => LieSE3.log(Matrix.identity(3, 3))).toThrow();
    });
});

function trace(M: Matrix): number {
    let sum = 0;
    for (let i = 0; i < M.numRows; ++i) {
        sum += M.get(i, i);
    }
    return sum;
}
