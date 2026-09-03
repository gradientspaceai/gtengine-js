import { describe, it, expect } from 'vitest';
import { AxisAngle } from '../src/AxisAngle.js';
import { GTE_C_HALF_PI, GTE_C_PI } from '../src/Constants.js';
import { EulerAngles, EulerResult } from '../src/EulerAngles.js';
import {
    Matrix, determinant, mulMatrix, multiplyAB, transpose
} from '../src/Matrix.js';
import { Quaternion, rotate } from '../src/Quaternion.js';
import { Rotation } from '../src/Rotation.js';
import { Vector, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, unitVector, vector
} from './helpers/arbitraries.js';

function expectMatrixClose(actual: Matrix, expected: Matrix,
    tolerance: number = 1e-12): void {
    expect(actual.numRows).toBe(expected.numRows);
    expect(actual.numCols).toBe(expected.numCols);
    for (let i = 0; i < expected.numElements; ++i) {
        expect(Math.abs(actual.values[i] - expected.values[i]),
            `element ${i}`).toBeLessThanOrEqual(tolerance);
    }
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// A uniformly-distributed unit quaternion (Shoemake's method).
function randomUnitQuaternion(rand: () => number): Quaternion {
    const u0 = rand(), u1 = rand(), u2 = rand();
    const s0 = Math.sqrt(1 - u0), s1 = Math.sqrt(u0);
    return new Quaternion(
        s0 * Math.sin(2 * GTE_C_PI * u1),
        s0 * Math.cos(2 * GTE_C_PI * u1),
        s1 * Math.sin(2 * GTE_C_PI * u2),
        s1 * Math.cos(2 * GTE_C_PI * u2));
}

function randomUnitAxis(rand: () => number): Vector {
    const z = 2 * rand() - 1;
    const phi = 2 * GTE_C_PI * rand();
    const r = Math.sqrt(Math.max(0, 1 - z * z));
    return Vector.fromArray([r * Math.cos(phi), r * Math.sin(phi), z]);
}

// Quaternions q and -q represent the same rotation; compare up to sign.
function quaternionDifferenceUpToSign(q0: Quaternion, q1: Quaternion): number {
    let same = 0, flipped = 0;
    for (let i = 0; i < 4; ++i) {
        same = Math.max(same, Math.abs(q0.values[i] - q1.values[i]));
        flipped = Math.max(flipped, Math.abs(q0.values[i] + q1.values[i]));
    }
    return Math.min(same, flipped);
}

// The 12 axis triples upstream supports: the six "Tait-Bryan" orders with
// distinct axes and the six "proper Euler" orders that repeat the first axis.
const DISTINCT_ORDERS: ReadonlyArray<[number, number, number]> = [
    [0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]
];
const REPEATED_ORDERS: ReadonlyArray<[number, number, number]> = [
    [0, 1, 0], [0, 2, 0], [1, 0, 1], [1, 2, 1], [2, 0, 2], [2, 1, 2]
];

describe('Rotation: axis-angle and matrix', () => {
    it('builds the known coordinate-axis rotations (GTE_USE_MAT_VEC)', () => {
        const t = 0.7;
        const c = Math.cos(t), s = Math.sin(t);
        const expected = [
            Matrix.fromArray(3, 3, [1, 0, 0, 0, c, -s, 0, s, c]),
            Matrix.fromArray(3, 3, [c, 0, s, 0, 1, 0, -s, 0, c]),
            Matrix.fromArray(3, 3, [c, -s, 0, s, c, 0, 0, 0, 1])
        ];
        for (let d = 0; d < 3; ++d) {
            const R = Rotation.fromAxisAngle(
                new AxisAngle(Vector.unit(3, d), t)).toMatrix();
            expectMatrixClose(R, expected[d], 1e-15);
        }
    });

    it('round-trips random axis-angle pairs through the matrix', () => {
        const rand = makeRandom(11111);
        for (let trial = 0; trial < 200; ++trial) {
            const axis = randomUnitAxis(rand);
            // Keep the angle in (0,pi), the range the inverse produces.
            const angle = 0.05 + rand() * (GTE_C_PI - 0.1);
            const a = new AxisAngle(axis, angle);
            const R = Rotation.fromAxisAngle(a).toMatrix();
            const b = Rotation.fromMatrix(R).toAxisAngle();
            expect(Math.abs(b.angle - angle)).toBeLessThanOrEqual(1e-9);
            for (let d = 0; d < 3; ++d) {
                expect(Math.abs(b.axis.values[d] - axis.values[d]))
                    .toBeLessThanOrEqual(1e-8);
            }
            // The matrix rebuilt from the recovered pair is the same.
            expectMatrixClose(Rotation.fromAxisAngle(b).toMatrix(), R, 1e-12);
        }
    });

    it('handles the identity: any axis works, upstream picks Unit(0)', () => {
        const a = Rotation.fromMatrix(Matrix.identity(3, 3)).toAxisAngle();
        expect(a.angle).toBe(0);
        expect(a.axis.values).toEqual([1, 0, 0]);
    });

    it('handles the angle-pi case, where R is symmetric', () => {
        // R(U,pi) = 2*U*U^T - I. These matrices have trace exactly -1, so
        // acos returns exactly pi and the symmetric branch (which normalizes
        // the largest row of R+I) is the one exercised. Each entry of the
        // list is [axis, R].
        const third = 1 / 3, twoThirds = 2 / 3;
        const cases: Array<[Vector, Matrix]> = [
            [Vector.fromArray([1, 0, 0]),
                Matrix.fromArray(3, 3, [1, 0, 0, 0, -1, 0, 0, 0, -1])],
            [Vector.fromArray([0, 1, 0]),
                Matrix.fromArray(3, 3, [-1, 0, 0, 0, 1, 0, 0, 0, -1])],
            [Vector.fromArray([0, 0, 1]),
                Matrix.fromArray(3, 3, [-1, 0, 0, 0, -1, 0, 0, 0, 1])],
            [Vector.fromArray([Math.SQRT1_2, Math.SQRT1_2, 0]),
                Matrix.fromArray(3, 3, [0, 1, 0, 1, 0, 0, 0, 0, -1])],
            [Vector.fromArray([0, Math.SQRT1_2, Math.SQRT1_2]),
                Matrix.fromArray(3, 3, [-1, 0, 0, 0, 0, 1, 0, 1, 0])],
            [Vector.fromArray([Math.SQRT1_2, 0, -Math.SQRT1_2]),
                Matrix.fromArray(3, 3, [0, 0, -1, 0, -1, 0, -1, 0, 0])],
            [Vector.fromArray([Math.sqrt(third), Math.sqrt(third),
                Math.sqrt(third)]),
                Matrix.fromArray(3, 3, [
                    -third, twoThirds, twoThirds,
                    twoThirds, -third, twoThirds,
                    twoThirds, twoThirds, -third])]
        ];
        for (const [U, R] of cases) {
            const a = Rotation.fromMatrix(R).toAxisAngle();
            expect(a.angle).toBe(GTE_C_PI);
            // R(U,pi) = R(-U,pi), so only the direction up to sign matters.
            expect(Math.abs(Math.abs(dot(a.axis, U)) - 1))
                .toBeLessThanOrEqual(1e-15);
            expectMatrixClose(Rotation.fromAxisAngle(a).toMatrix(), R, 1e-15);
        }
    });

    it('is ill-conditioned but stable near angle pi, as upstream is', () => {
        // Just below pi the antisymmetric extraction is used and acos loses
        // about half the digits (the derivative of acos is unbounded at -1).
        // The recovered pair still rebuilds the input matrix.
        const rand = makeRandom(60606);
        for (let trial = 0; trial < 50; ++trial) {
            const axis = randomUnitAxis(rand);
            const angle = GTE_C_PI - 1e-7 * rand();
            const R = Rotation.fromAxisAngle(
                new AxisAngle(axis, angle)).toMatrix();
            const a = Rotation.fromMatrix(R).toAxisAngle();
            expect(Math.abs(a.angle - angle)).toBeLessThanOrEqual(1e-7);
            expectMatrixClose(Rotation.fromAxisAngle(a).toMatrix(), R, 1e-7);
        }
    });
});

describe('Rotation: quaternion and matrix', () => {
    it('maps the identity quaternion to the identity matrix', () => {
        expectMatrixClose(Rotation.fromQuaternion(Quaternion.identity())
            .toMatrix(), Matrix.identity(3, 3), 0);
    });

    it('round-trips random quaternions through the matrix (up to sign)',
        () => {
            const rand = makeRandom(2468);
            for (let trial = 0; trial < 300; ++trial) {
                const q = randomUnitQuaternion(rand);
                const R = Rotation.fromQuaternion(q).toMatrix();
                // The matrix must be a rotation.
                expectMatrixClose(multiplyAB(transpose(R), R),
                    Matrix.identity(3, 3), 1e-14);
                const q2 = Rotation.fromMatrix(R).toQuaternion();
                expect(quaternionDifferenceUpToSign(q, q2))
                    .toBeLessThanOrEqual(1e-14);
                expectMatrixClose(Rotation.fromQuaternion(q2).toMatrix(), R,
                    1e-14);
            }
        });

    it('exercises all four branches of the matrix-to-quaternion extraction',
        () => {
            // Each branch is selected by the signs of r22 and of r11 -/+ r00;
            // 180-degree rotations about x, y, z and the identity hit them.
            const cases: Quaternion[] = [
                new Quaternion(1, 0, 0, 0),                       // x branch
                new Quaternion(0, 1, 0, 0),                       // y branch
                new Quaternion(0, 0, 1, 0),                       // z branch
                new Quaternion(0, 0, 0, 1)                        // w branch
            ];
            for (const q of cases) {
                const R = Rotation.fromQuaternion(q).toMatrix();
                const q2 = Rotation.fromMatrix(R).toQuaternion();
                expect(quaternionDifferenceUpToSign(q, q2))
                    .toBeLessThanOrEqual(1e-15);
            }
        });

    it('agrees with the axis-angle conversion on random rotations', () => {
        const rand = makeRandom(13579);
        for (let trial = 0; trial < 200; ++trial) {
            const axis = randomUnitAxis(rand);
            const angle = 0.05 + rand() * (GTE_C_PI - 0.1);
            const a = new AxisAngle(axis, angle);
            const q = Rotation.fromAxisAngle(a).toQuaternion();
            // q = sin(angle/2)*axis + cos(angle/2)
            const sn = Math.sin(0.5 * angle);
            expect(Math.abs(q.values[0] - sn * axis.values[0]))
                .toBeLessThanOrEqual(1e-15);
            expect(Math.abs(q.values[3] - Math.cos(0.5 * angle)))
                .toBeLessThanOrEqual(1e-15);
            const b = Rotation.fromQuaternion(q).toAxisAngle();
            expect(Math.abs(b.angle - angle)).toBeLessThanOrEqual(1e-8);
            expectMatrixClose(Rotation.fromQuaternion(q).toMatrix(),
                Rotation.fromAxisAngle(a).toMatrix(), 1e-14);
        }
    });

    it('maps the zero-axis quaternion to the Unit(0) axis and angle 0', () => {
        const a = Rotation.fromQuaternion(new Quaternion(0, 0, 0, 1))
            .toAxisAngle();
        expect(a.angle).toBe(0);
        expect(a.axis.values).toEqual([1, 0, 0]);
    });
});

describe('Rotation: Euler angles', () => {
    it('composes as R(N2,a2)*R(N1,a1)*R(N0,a0) (GTE_USE_MAT_VEC)', () => {
        const e = new EulerAngles(0, 1, 2, 0.3, -0.4, 0.5);
        const R = Rotation.fromEulerAngles(e).toMatrix();
        const r0 = Rotation.fromAxisAngle(
            new AxisAngle(Vector.unit(3, 0), 0.3)).toMatrix();
        const r1 = Rotation.fromAxisAngle(
            new AxisAngle(Vector.unit(3, 1), -0.4)).toMatrix();
        const r2 = Rotation.fromAxisAngle(
            new AxisAngle(Vector.unit(3, 2), 0.5)).toMatrix();
        expectMatrixClose(R, multiplyAB(multiplyAB(r2, r1), r0), 1e-15);
    });

    it('round-trips the six distinct-axis orders', () => {
        const rand = makeRandom(97531);
        for (const [i0, i1, i2] of DISTINCT_ORDERS) {
            for (let trial = 0; trial < 40; ++trial) {
                // For distinct axes, angle1 is in (-pi/2,pi/2) and the outer
                // angles are in (-pi,pi); stay away from the gimbal edges.
                const a0 = (2 * rand() - 1) * (GTE_C_PI - 0.05);
                const a1 = (2 * rand() - 1) * (GTE_C_HALF_PI - 0.1);
                const a2 = (2 * rand() - 1) * (GTE_C_PI - 0.05);
                const e = new EulerAngles(i0, i1, i2, a0, a1, a2);
                const R = Rotation.fromEulerAngles(e).toMatrix();
                const f = Rotation.fromMatrix(R).toEulerAngles(i0, i1, i2);
                expect(f.result).toBe(EulerResult.UNIQUE);
                expect(f.axis).toEqual([i0, i1, i2]);
                expect(Math.abs(f.angle[0] - a0)).toBeLessThanOrEqual(1e-9);
                expect(Math.abs(f.angle[1] - a1)).toBeLessThanOrEqual(1e-9);
                expect(Math.abs(f.angle[2] - a2)).toBeLessThanOrEqual(1e-9);
                expectMatrixClose(Rotation.fromEulerAngles(f).toMatrix(), R,
                    1e-12);
            }
        }
    });

    it('round-trips the six repeated-axis orders', () => {
        const rand = makeRandom(86420);
        for (const [i0, i1, i2] of REPEATED_ORDERS) {
            for (let trial = 0; trial < 40; ++trial) {
                // For repeated axes, angle1 is in (0,pi).
                const a0 = (2 * rand() - 1) * (GTE_C_PI - 0.05);
                const a1 = 0.05 + rand() * (GTE_C_PI - 0.1);
                const a2 = (2 * rand() - 1) * (GTE_C_PI - 0.05);
                const e = new EulerAngles(i0, i1, i2, a0, a1, a2);
                const R = Rotation.fromEulerAngles(e).toMatrix();
                const f = Rotation.fromMatrix(R).toEulerAngles(i0, i1, i2);
                expect(f.result).toBe(EulerResult.UNIQUE);
                expect(Math.abs(f.angle[0] - a0)).toBeLessThanOrEqual(1e-9);
                expect(Math.abs(f.angle[1] - a1)).toBeLessThanOrEqual(1e-9);
                expect(Math.abs(f.angle[2] - a2)).toBeLessThanOrEqual(1e-9);
                expectMatrixClose(Rotation.fromEulerAngles(f).toMatrix(), R,
                    1e-12);
            }
        }
    });

    it('reports the gimbal cases for distinct axes', () => {
        // Rotation by +90 degrees about y has r(2,0) = -1 exactly.
        const yPlus = Matrix.fromArray(3, 3, [0, 0, 1, 0, 1, 0, -1, 0, 0]);
        const dif = Rotation.fromMatrix(yPlus).toEulerAngles(0, 1, 2);
        expect(dif.result).toBe(EulerResult.NOT_UNIQUE_DIF);
        // In the nonunique cases upstream returns angleN0 = 0. For the
        // MAT_VEC branch that is e.angle[2].
        expect(dif.angle[2]).toBe(0);
        expect(Math.abs(dif.angle[1] - GTE_C_HALF_PI))
            .toBeLessThanOrEqual(1e-15);
        expectMatrixClose(Rotation.fromEulerAngles(dif).toMatrix(), yPlus,
            1e-15);

        // Rotation by -90 degrees about y has r(2,0) = +1 exactly.
        const yMinus = Matrix.fromArray(3, 3, [0, 0, -1, 0, 1, 0, 1, 0, 0]);
        const sum = Rotation.fromMatrix(yMinus).toEulerAngles(0, 1, 2);
        expect(sum.result).toBe(EulerResult.NOT_UNIQUE_SUM);
        expect(sum.angle[2]).toBe(0);
        expect(Math.abs(sum.angle[1] + GTE_C_HALF_PI))
            .toBeLessThanOrEqual(1e-15);
        expectMatrixClose(Rotation.fromEulerAngles(sum).toMatrix(), yMinus,
            1e-15);
    });

    it('reports the gimbal cases for repeated axes', () => {
        // For order (0,1,0) the branch is selected by r(0,0).
        const identity = Matrix.identity(3, 3);
        const sum = Rotation.fromMatrix(identity).toEulerAngles(0, 1, 0);
        expect(sum.result).toBe(EulerResult.NOT_UNIQUE_SUM);
        expect(sum.angle[1]).toBe(0);
        expect(sum.angle[2]).toBe(0);
        expectMatrixClose(Rotation.fromEulerAngles(sum).toMatrix(), identity,
            1e-15);

        // Rotation by 180 degrees about z has r(0,0) = -1 exactly.
        const zPi = Matrix.fromArray(3, 3, [-1, 0, 0, 0, -1, 0, 0, 0, 1]);
        const dif = Rotation.fromMatrix(zPi).toEulerAngles(0, 1, 0);
        expect(dif.result).toBe(EulerResult.NOT_UNIQUE_DIF);
        expect(dif.angle[1]).toBe(GTE_C_PI);
        expect(dif.angle[2]).toBe(0);
        expectMatrixClose(Rotation.fromEulerAngles(dif).toMatrix(), zPi,
            1e-15);
    });

    it('reports INVALID for bad axis triples and yields the identity', () => {
        const R = Rotation.fromQuaternion(
            randomUnitQuaternion(makeRandom(5))).toMatrix();
        for (const [i0, i1, i2] of [[0, 0, 1], [1, 1, 1], [0, 1, 1],
            [-1, 1, 2], [0, 3, 2], [0, 1, 5]] as Array<
                [number, number, number]>) {
            const e = Rotation.fromMatrix(R).toEulerAngles(i0, i1, i2);
            expect(e.result).toBe(EulerResult.INVALID);
            expect(e.angle).toEqual([0, 0, 0]);
            // Converting invalid Euler angles back gives the identity.
            expectMatrixClose(Rotation.fromEulerAngles(e).toMatrix(),
                Matrix.identity(3, 3), 0);
        }
    });

    it('reaches Euler angles from quaternions and axis-angle pairs', () => {
        const rand = makeRandom(24680);
        for (let trial = 0; trial < 100; ++trial) {
            const q = randomUnitQuaternion(rand);
            const R = Rotation.fromQuaternion(q).toMatrix();
            const fromMatrix = Rotation.fromMatrix(R).toEulerAngles(2, 0, 1);
            const fromQuat =
                Rotation.fromQuaternion(q).toEulerAngles(2, 0, 1);
            expect(fromQuat.result).toBe(fromMatrix.result);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(fromQuat.angle[i] - fromMatrix.angle[i]))
                    .toBeLessThanOrEqual(1e-12);
            }

            const a = Rotation.fromQuaternion(q).toAxisAngle();
            const fromAxisAngle =
                Rotation.fromAxisAngle(a).toEulerAngles(2, 0, 1);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(fromAxisAngle.angle[i] - fromMatrix.angle[i]))
                    .toBeLessThanOrEqual(1e-7);
            }
        }
    });

    it('reaches quaternions and axis-angle pairs from Euler angles', () => {
        const rand = makeRandom(1122334);
        for (let trial = 0; trial < 100; ++trial) {
            const a0 = (2 * rand() - 1) * 3;
            const a1 = (2 * rand() - 1) * 1.4;
            const a2 = (2 * rand() - 1) * 3;
            const e = new EulerAngles(1, 2, 0, a0, a1, a2);
            const R = Rotation.fromEulerAngles(e).toMatrix();
            const q = Rotation.fromEulerAngles(e).toQuaternion();
            expectMatrixClose(Rotation.fromQuaternion(q).toMatrix(), R, 1e-14);
            const aa = Rotation.fromEulerAngles(e).toAxisAngle();
            expectMatrixClose(Rotation.fromAxisAngle(aa).toMatrix(), R, 1e-12);
        }
    });

    it('returns the stored factorization unchanged when the axes match',
        () => {
            const e = new EulerAngles(2, 0, 1, 0.1, 0.2, 0.3);
            const f = Rotation.fromEulerAngles(e).toEulerAngles(2, 0, 1);
            expect(f.axis).toEqual([2, 0, 1]);
            expect(f.angle).toEqual([0.1, 0.2, 0.3]);
            expect(f.result).toBe(EulerResult.UNIQUE);
        });

    it('recomputes when a different factorization is requested (upstream bug '
        + 'fix)', () => {
            // Upstream's operator()(i0,i1,i2) relabels the cached angles with
            // the requested axes without recomputing them, which describes a
            // different rotation. The port converts through the matrix.
            const e = new EulerAngles(0, 1, 2, 0.3, -0.4, 0.5);
            const R = Rotation.fromEulerAngles(e).toMatrix();
            const f = Rotation.fromEulerAngles(e).toEulerAngles(2, 1, 0);
            expect(f.axis).toEqual([2, 1, 0]);
            expect(f.result).toBe(EulerResult.UNIQUE);
            expectMatrixClose(Rotation.fromEulerAngles(f).toMatrix(), R,
                1e-12);
            // The upstream behavior would have kept (0.3,-0.4,0.5).
            expect(f.angle).not.toEqual([0.3, -0.4, 0.5]);
        });
});

describe('Rotation: dimension 4 (affine)', () => {
    it('produces 4x4 matrices whose upper-left 3x3 block is the rotation',
        () => {
            const rand = makeRandom(444);
            for (let trial = 0; trial < 50; ++trial) {
                const q = randomUnitQuaternion(rand);
                const R3 = Rotation.fromQuaternion(q, 3).toMatrix();
                const R4 = Rotation.fromQuaternion(q, 4).toMatrix();
                expect(R4.numRows).toBe(4);
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 3; ++c) {
                        expect(R4.get(r, c)).toBe(R3.get(r, c));
                    }
                }
                for (let i = 0; i < 3; ++i) {
                    expect(R4.get(i, 3)).toBe(0);
                    expect(R4.get(3, i)).toBe(0);
                }
                expect(R4.get(3, 3)).toBe(1);
            }
        });

    it('round-trips 4D axis-angle pairs, keeping the w-component zero', () => {
        const rand = makeRandom(555);
        for (let trial = 0; trial < 50; ++trial) {
            const axis3 = randomUnitAxis(rand);
            const axis4 = Vector.fromArray([axis3.values[0], axis3.values[1],
                axis3.values[2], 0]);
            const angle = 0.05 + rand() * (GTE_C_PI - 0.1);
            const R4 = Rotation.fromAxisAngle(
                new AxisAngle(axis4, angle)).toMatrix();
            const a = Rotation.fromMatrix(R4).toAxisAngle();
            expect(a.axis.size).toBe(4);
            expect(a.axis.values[3]).toBe(0);
            expect(Math.abs(a.angle - angle)).toBeLessThanOrEqual(1e-9);
            expectMatrixClose(Rotation.fromAxisAngle(a).toMatrix(), R4, 1e-12);
        }
    });

    it('reaches Euler angles in 4D', () => {
        const e = new EulerAngles(0, 1, 2, 0.25, -0.5, 1.1);
        const R4 = Rotation.fromEulerAngles(e, 4).toMatrix();
        expect(R4.numRows).toBe(4);
        const f = Rotation.fromMatrix(R4).toEulerAngles(0, 1, 2);
        expect(f.result).toBe(EulerResult.UNIQUE);
        expect(Math.abs(f.angle[0] - 0.25)).toBeLessThanOrEqual(1e-12);
        expect(Math.abs(f.angle[1] + 0.5)).toBeLessThanOrEqual(1e-12);
        expect(Math.abs(f.angle[2] - 1.1)).toBeLessThanOrEqual(1e-12);
    });
});

describe('Rotation: construction and copying', () => {
    it('rejects dimensions other than 3 and 4', () => {
        expect(() => Rotation.fromQuaternion(Quaternion.identity(), 2))
            .toThrow('Dimension must be 3 or 4.');
        expect(() => Rotation.fromMatrix(Matrix.identity(2, 2)))
            .toThrow('Dimension must be 3 or 4.');
        expect(() => Rotation.fromMatrix(new Matrix(3, 4)))
            .toThrow('Rotation: expecting a square matrix.');
    });

    it('copies its inputs and its outputs', () => {
        const q = new Quaternion(0, 0, Math.SQRT1_2, Math.SQRT1_2);
        const rotation = Rotation.fromQuaternion(q);
        q.values[0] = 99;  // mutating the input must not affect the rotation
        const R = rotation.toMatrix();
        R.set(0, 0, 12345);  // mutating the output must not affect the cache
        expectMatrixClose(rotation.toMatrix(),
            Matrix.fromArray(3, 3, [0, -1, 0, 1, 0, 0, 0, 0, 1]), 1e-15);

        const a = new AxisAngle(Vector.unit(3, 2), 1.0);
        const rotation2 = Rotation.fromAxisAngle(a);
        a.axis.values[2] = 0;
        a.angle = 0;
        expectMatrixClose(rotation2.toMatrix(),
            Rotation.fromAxisAngle(
                new AxisAngle(Vector.unit(3, 2), 1.0)).toMatrix(), 0);

        const e = new EulerAngles(0, 1, 2, 0.1, 0.2, 0.3);
        const rotation3 = Rotation.fromEulerAngles(e);
        e.angle[0] = 9;
        const f = rotation3.toEulerAngles(0, 1, 2);
        expect(f.angle[0]).toBe(0.1);
    });

    it('is idempotent: converting to the source representation is a no-op',
        () => {
            const R = Rotation.fromQuaternion(
                randomUnitQuaternion(makeRandom(31))).toMatrix();
            expectMatrixClose(Rotation.fromMatrix(R).toMatrix(), R, 0);

            const q = new Quaternion(0.5, 0.5, 0.5, 0.5);
            expect(Rotation.fromQuaternion(q).toQuaternion().values)
                .toEqual(q.values);

            const a = new AxisAngle(Vector.unit(3, 1), 0.9);
            const b = Rotation.fromAxisAngle(a).toAxisAngle();
            expect(b.angle).toBe(0.9);
            expect(b.axis.values).toEqual([0, 1, 0]);
        });
});

// ---------------------------------------------------------------------------
// Verification wave (V02): property-based re-checks against Rotation.h.
// ---------------------------------------------------------------------------

// The elementary rotation matrices of the upstream comment block, written out
// by hand so that they are an independent reference for the conversions:
//   R(unit(0),t) = {{ 1, 0, 0}, { 0, c,-s}, { 0, s, c}}
//   R(unit(1),t) = {{ c, 0, s}, { 0, 1, 0}, {-s, 0, c}}
//   R(unit(2),t) = {{ c,-s, 0}, { s, c, 0}, { 0, 0, 1}}
function elementaryRotation(axis: number, t: number): Matrix {
    const c = Math.cos(t), s = Math.sin(t);
    switch (axis) {
        case 0: return Matrix.fromArray(3, 3, [1, 0, 0, 0, c, -s, 0, s, c]);
        case 1: return Matrix.fromArray(3, 3, [c, 0, s, 0, 1, 0, -s, 0, c]);
        default: return Matrix.fromArray(3, 3, [c, -s, 0, s, c, 0, 0, 0, 1]);
    }
}

const ALL_ORDERS: ReadonlyArray<[number, number, number]> =
    [...DISTINCT_ORDERS, ...REPEATED_ORDERS];

const arbAxis = () => unitVector(3);
const arbAngle = () => finite(-GTE_C_PI, GTE_C_PI);
const arbUnitQuat = () => unitVector(4).map(v => Quaternion.fromArray(v.values));
const arbOrder = () => fc.constantFrom(...ALL_ORDERS);

function expectRotationMatrix(R: Matrix, tolerance = 1e-12): void {
    expectMatrixClose(multiplyAB(R, transpose(R)),
        Matrix.identity(R.numRows, R.numCols), tolerance);
    expectClose(determinant(R), 1, tolerance, tolerance);
}

describe('Rotation verification', () => {
    it('axis-angle to matrix is a rotation that fixes the axis and turns by '
        + 'the angle', () => {
            check(fc.tuple(arbAxis(), arbAngle()), ([axis, angle]) => {
                const R = Rotation.fromAxisAngle(new AxisAngle(axis, angle))
                    .toMatrix();
                expectRotationMatrix(R, 1e-12);

                // The axis is the eigenvector for eigenvalue 1.
                expectVectorClose(mulMatrix(R, axis) as Vector, axis, 1e-12,
                    1e-12);

                // trace(R) = 1 + 2*cos(angle) for every rotation.
                expectClose(R.get(0, 0) + R.get(1, 1) + R.get(2, 2),
                    1 + 2 * Math.cos(angle), 1e-12, 1e-12);

                // A vector perpendicular to the axis turns by exactly the
                // angle, with the sign fixed by the right-hand rule.
                const seed = Vector.fromArray([axis.get(1), -axis.get(2),
                    axis.get(0) + 0.5]);
                const perp = sub(seed, mul(axis, dot(axis, seed)));
                if (length(perp) > 1e-2) {
                    normalize(perp);
                    const turned = mulMatrix(R, perp) as Vector;
                    expectClose(dot(turned, perp), Math.cos(angle), 1e-10,
                        1e-10);
                    expectClose(dot(cross(perp, turned), axis),
                        Math.sin(angle), 1e-10, 1e-10);
                }
            });
        });

    it('the axis-angle round trip through the matrix preserves the rotation',
        () => {
            check(fc.tuple(arbAxis(), finite(1e-3, GTE_C_PI - 1e-3)),
                ([axis, angle]) => {
                    const R = Rotation.fromAxisAngle(new AxisAngle(axis, angle))
                        .toMatrix();
                    const aa = Rotation.fromMatrix(R).toAxisAngle();
                    // The recovered angle is in [0,pi]; for a positive input
                    // angle in (0,pi) the pair is recovered exactly.
                    expectClose(aa.angle, angle, 1e-8, 1e-8);
                    expectVectorClose(aa.axis, axis, 1e-7, 1e-7);
                    expectMatrixClose(Rotation.fromAxisAngle(aa).toMatrix(), R,
                        1e-9);
                });
        });

    it('quaternion and matrix conversions are mutually inverse and agree with '
        + 'Quaternion.rotate', () => {
            check(fc.tuple(arbUnitQuat(), vector(3, -5, 5)), ([q, u]) => {
                const R = Rotation.fromQuaternion(q).toMatrix();
                expectRotationMatrix(R, 1e-12);

                // The matrix is the one Quaternion.rotate implements; that
                // function is an independent (Eisele) derivation.
                expectVectorClose(mulMatrix(R, u) as Vector, rotate(q, u),
                    1e-10, 1e-10);

                // Matrix -> quaternion recovers q up to sign.
                const q2 = Rotation.fromMatrix(R).toQuaternion();
                expect(quaternionDifferenceUpToSign(q, q2))
                    .toBeLessThanOrEqual(1e-8);
                expectMatrixClose(Rotation.fromQuaternion(q2).toMatrix(), R,
                    1e-10);
            });
        });

    it('quaternion and axis-angle conversions agree through the matrix', () => {
        check(arbUnitQuat(), q => {
            const vectorPart = Math.hypot(q.get(0), q.get(1), q.get(2));
            if (vectorPart < 1e-6) {
                // The rotation is within acos resolution of the identity, so
                // its axis is not determined: acos(1 - 1e-16) already returns
                // 1.5e-8 rather than 0.
                return;
            }
            const viaQuaternion = Rotation.fromQuaternion(q).toAxisAngle();
            const R = Rotation.fromQuaternion(q).toMatrix();
            const viaMatrix = Rotation.fromMatrix(R).toAxisAngle();
            // Both paths describe the same rotation matrix. The tolerance is
            // loose because extracting an angle near pi from a matrix is
            // ill-conditioned: acos'(x) blows up at x = -1, so the recovered
            // angle carries an error of order sqrt(eps).
            expectMatrixClose(Rotation.fromAxisAngle(viaQuaternion).toMatrix(),
                Rotation.fromAxisAngle(viaMatrix).toMatrix(), 1e-6);

            // axis-angle -> quaternion inverts quaternion -> axis-angle up to
            // the sign of the quaternion.
            const q2 = Rotation.fromAxisAngle(viaQuaternion).toQuaternion();
            // Loose because recovering an angle from acos loses half the
            // digits when the angle is near 0 or pi.
            expect(quaternionDifferenceUpToSign(q, q2))
                .toBeLessThanOrEqual(1e-6);
        });
    });

    it('Euler angles compose as R2*R1*R0 for all 12 axis orders', () => {
        check(fc.tuple(arbOrder(), arbAngle(), arbAngle(), arbAngle()),
            ([order, a0, a1, a2]) => {
                const e = new EulerAngles(order[0], order[1], order[2],
                    a0, a1, a2);
                const R = Rotation.fromEulerAngles(e).toMatrix();
                const expected = multiplyAB(
                    multiplyAB(elementaryRotation(order[2], a2),
                        elementaryRotation(order[1], a1)),
                    elementaryRotation(order[0], a0));
                expectMatrixClose(R, expected, 1e-13);
                expectRotationMatrix(R, 1e-12);
            });
    });

    it('the Euler factorization of a matrix reproduces that matrix, for all '
        + '12 axis orders and both uniqueness cases', () => {
            check(fc.tuple(arbOrder(), arbUnitQuat()), ([order, q]) => {
                const R = Rotation.fromQuaternion(q).toMatrix();
                const e = Rotation.fromMatrix(R).toEulerAngles(order[0],
                    order[1], order[2]);
                expect(e.result).not.toBe(EulerResult.INVALID);
                expect(e.axis).toEqual([order[0], order[1], order[2]]);

                // Upstream switches to the non-unique branches only when the
                // deciding entry reaches exactly +-1. Within an ulp of that
                // boundary the unique branch is taken with atan2 arguments
                // that are pure round-off, and the factorization it returns
                // need not describe R at all (reported in the PR as an
                // upstream suspect). Skip that neighborhood.
                const deciding = (order[0] !== order[2]
                    ? R.get(order[2], order[0])
                    : R.get(order[2], order[2]));
                if (Math.abs(deciding) > 1 - 1e-9) {
                    return;
                }
                // Recomposition must return the original rotation whether the
                // factorization is unique or not (in the non-unique cases the
                // angles differ but the product does not). The tolerance is
                // loose because a factorization close to gimbal lock has an
                // ill-conditioned atan2 pair: the angles carry a large
                // relative error that cancels in the product.
                expectMatrixClose(Rotation.fromEulerAngles(e).toMatrix(), R,
                    1e-6);
            });
        });

    it('the Euler factorization at gimbal lock is flagged and still '
        + 'recomposes', () => {
            check(fc.tuple(fc.constantFrom(...DISTINCT_ORDERS), arbAngle(),
                fc.constantFrom(-1, 1)), ([order, a2, sign]) => {
                    // Middle angle +-pi/2 collapses the outer two angles.
                    const e = new EulerAngles(order[0], order[1], order[2],
                        0, sign * GTE_C_HALF_PI, a2);
                    const R = Rotation.fromEulerAngles(e).toMatrix();
                    const f = Rotation.fromMatrix(R).toEulerAngles(order[0],
                        order[1], order[2]);

                    // Upstream reports non-uniqueness from the test
                    // |r(i2,i0)| >= 1. At gimbal lock that entry is +-1 in
                    // exact arithmetic, but building R from cos/sin can leave
                    // it one ulp inside the interval, in which case the
                    // UNIQUE branch is taken with an ill-conditioned atan2
                    // pair. Both outcomes must recompose to R.
                    const locked = Math.abs(R.get(order[2], order[0])) >= 1;
                    if (locked) {
                        expect(f.result === EulerResult.NOT_UNIQUE_SUM
                            || f.result === EulerResult.NOT_UNIQUE_DIF)
                            .toBe(true);
                        expect(f.angle[2]).toBe(0);
                    }
                    expectMatrixClose(Rotation.fromEulerAngles(f).toMatrix(), R,
                        1e-7);
                });
        });

    it('an invalid axis order yields INVALID angles and the identity matrix',
        () => {
            check(fc.tuple(fc.integer({ min: -2, max: 4 }),
                fc.integer({ min: -2, max: 4 }), fc.integer({ min: -2, max: 4 }),
                arbUnitQuat()), ([i0, i1, i2, q]) => {
                    const valid = 0 <= i0 && i0 < 3 && 0 <= i1 && i1 < 3
                        && 0 <= i2 && i2 < 3 && i1 !== i0 && i1 !== i2;
                    if (valid) {
                        return;
                    }
                    const R = Rotation.fromQuaternion(q).toMatrix();
                    const e = Rotation.fromMatrix(R).toEulerAngles(i0, i1, i2);
                    expect(e.result).toBe(EulerResult.INVALID);
                    expect(e.angle).toEqual([0, 0, 0]);

                    const bad = new EulerAngles(i0, i1, i2, 1, 2, 3);
                    expectMatrixClose(Rotation.fromEulerAngles(bad).toMatrix(),
                        Matrix.identity(3, 3), 0);
                });
        });

    // Regression test for the documented deviation from upstream: upstream's
    // operator()(i0,i1,i2) relabels the cached Euler angles when the source
    // representation is IS_EULER_ANGLES, returning the original angles under
    // the requested axis order -- a different rotation. The port refactors.
    it('re-factoring an Euler-angle rotation into another axis order gives '
        + 'the same rotation (upstream returns the relabeled angles)', () => {
            check(fc.tuple(arbOrder(), arbOrder(), arbAngle(), arbAngle(),
                arbAngle()), ([from, to, a0, a1, a2]) => {
                    const e = new EulerAngles(from[0], from[1], from[2],
                        a0, a1, a2);
                    const rotation = Rotation.fromEulerAngles(e);
                    const R = rotation.toMatrix();

                    // Skip the gimbal-lock boundary for the target order; see
                    // the note in the factorization property above.
                    const deciding = (to[0] !== to[2]
                        ? R.get(to[2], to[0])
                        : R.get(to[2], to[2]));
                    if (Math.abs(deciding) > 1 - 1e-9) {
                        return;
                    }

                    const f = rotation.toEulerAngles(to[0], to[1], to[2]);
                    expect(f.axis).toEqual([to[0], to[1], to[2]]);
                    expectMatrixClose(Rotation.fromEulerAngles(f).toMatrix(), R,
                        1e-8);

                    // Asking for the original order returns the original
                    // angles untouched, as upstream does.
                    const same = Rotation.fromEulerAngles(e)
                        .toEulerAngles(from[0], from[1], from[2]);
                    expect(same.angle).toEqual([a0, a1, a2]);
                });
        });

    it('every accessor returns a copy that does not alias the cache', () => {
        check(fc.tuple(arbAxis(), arbAngle(), arbOrder()),
            ([axis, angle, order]) => {
                const rotation = Rotation.fromAxisAngle(
                    new AxisAngle(axis, angle));

                const m0 = rotation.toMatrix();
                m0.set(0, 0, 12345);
                expect(rotation.toMatrix().get(0, 0)).not.toBe(12345);

                const aa = rotation.toAxisAngle();
                aa.axis.set(0, 12345);
                expect(rotation.toAxisAngle().axis.get(0)).not.toBe(12345);

                const q = rotation.toQuaternion();
                q.set(0, 12345);
                expect(rotation.toQuaternion().get(0)).not.toBe(12345);

                const e = rotation.toEulerAngles(order[0], order[1], order[2]);
                e.angle[0] = 12345;
                expect(rotation.toEulerAngles(order[0], order[1], order[2])
                    .angle[0]).not.toBe(12345);

                // The constructor copies its argument too.
                const axisCopy = axis.clone();
                const fromAA = Rotation.fromAxisAngle(
                    new AxisAngle(axisCopy, angle));
                axisCopy.set(0, 12345);
                expectMatrixClose(fromAA.toMatrix(), rotation.toMatrix(), 0);
            });
    });

    // Regression test for the second documented deviation from upstream: a
    // rotation matrix whose angle is numerically just below pi has an
    // antisymmetric part that underflows to exactly zero, and upstream's
    // 'angle < pi' branch then returns a zero axis. Before the fix, the axis
    // was (0,0,0) and the round trip produced -I.
    it('extracts a unit axis from a matrix whose angle is numerically just '
        + 'below pi', () => {
            const q = new Quaternion(0, 0, -0.9999999999999999, 0);
            const R = Rotation.fromQuaternion(q).toMatrix();
            // The premise of the regression: acos reports an angle strictly
            // below pi while R is symmetric to the last bit.
            const trace = R.get(0, 0) + R.get(1, 1) + R.get(2, 2);
            expect(0.5 * (trace - 1)).toBeGreaterThan(-1);
            expect(R.get(2, 1) - R.get(1, 2)).toBe(0);
            expect(R.get(0, 2) - R.get(2, 0)).toBe(0);
            expect(R.get(1, 0) - R.get(0, 1)).toBe(0);

            const aa = Rotation.fromMatrix(R).toAxisAngle();
            expectClose(dot(aa.axis, aa.axis), 1, 1e-12, 1e-12);
            expectClose(Math.abs(aa.axis.get(2)), 1, 1e-12, 1e-12);
            expectMatrixClose(Rotation.fromAxisAngle(aa).toMatrix(), R, 1e-7);
        });

    // The same defect with a merely underflowed (not zero) antisymmetric
    // part: before the fix the extracted axis had length 0.707.
    it('extracts a unit axis when the antisymmetric part underflows to '
        + 'denormals', () => {
            const q = new Quaternion(0, 0.9999999720664132,
                -0.00023636237563083313, 3.929319571516562e-163);
            const R = Rotation.fromQuaternion(q).toMatrix();
            const aa = Rotation.fromMatrix(R).toAxisAngle();
            expectClose(length(aa.axis), 1, 1e-12, 1e-12);
            expectMatrixClose(Rotation.fromAxisAngle(aa).toMatrix(), R, 1e-6);
        });

    it('never returns a degenerate axis for a rotation matrix', () => {
        check(arbUnitQuat(), q => {
            const R = Rotation.fromQuaternion(q).toMatrix();
            const aa = Rotation.fromMatrix(R).toAxisAngle();
            // The axis is always unit length: the zero-angle case is covered
            // by upstream's explicit Unit(0) fallback and every other path
            // normalizes (see the port's guard in convertMatrixToAxisAngle).
            expectClose(length(aa.axis), 1, 1e-11, 1e-11);

            if (Math.abs(q.get(3)) > 1e-3) {
                // The round trip is only checked away from a rotation angle
                // of pi. There, R is symmetric to within round-off and the
                // entries of R - Transpose(R) that carry the axis are
                // absorbed by the larger entries they are added to, so the
                // extracted direction can be wrong by a finite amount. That
                // is an upstream limitation reported in the PR, not a port
                // defect; the port only guarantees a valid (unit) axis.
                expectMatrixClose(Rotation.fromAxisAngle(aa).toMatrix(), R,
                    1e-6);
            }
        });
    });

    it('dimension 4 embeds the 3D rotation with an unchanged w channel', () => {
        check(fc.tuple(arbUnitQuat(), arbOrder(), arbAngle(), arbAngle(),
            arbAngle()), ([q, order, a0, a1, a2]) => {
                const R3 = Rotation.fromQuaternion(q, 3).toMatrix();
                const R4 = Rotation.fromQuaternion(q, 4).toMatrix();
                expect(R4.numRows).toBe(4);
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 3; ++c) {
                        expectClose(R4.get(r, c), R3.get(r, c), 0, 0);
                    }
                    expectClose(R4.get(r, 3), 0, 0, 0);
                    expectClose(R4.get(3, r), 0, 0, 0);
                }
                expectClose(R4.get(3, 3), 1, 0, 0);

                const e = new EulerAngles(order[0], order[1], order[2],
                    a0, a1, a2);
                const E3 = Rotation.fromEulerAngles(e, 3).toMatrix();
                const E4 = Rotation.fromEulerAngles(e, 4).toMatrix();
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 3; ++c) {
                        expectClose(E4.get(r, c), E3.get(r, c), 1e-14, 1e-14);
                    }
                }
            });
    });
});
