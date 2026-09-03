import { describe, it, expect } from 'vitest';
import { Transform, mulTransform } from '../src/Transform.js';
import { Matrix, mulMatrix, multiplyAB, transpose } from '../src/Matrix.js';
import { inverse4x4 } from '../src/Matrix4x4.js';
import { Vector, hlift, normalize } from '../src/Vector.js';
import { GTE_C_PI } from '../src/Constants.js';
import {
    check, expectClose, expectVectorClose as expectVectorsClose, fc, finite,
    invertibleMatrix, matrix, unitVector, vector
} from './helpers/arbitraries.js';

import { AxisAngle } from '../src/AxisAngle.js';
import { EulerAngles } from '../src/EulerAngles.js';
import { Rotation } from '../src/Rotation.js';

function expectMatrixClose(actual: Matrix, expected: Matrix,
    tolerance: number = 1e-12): void {
    expect(actual.numRows).toBe(expected.numRows);
    expect(actual.numCols).toBe(expected.numCols);
    for (let i = 0; i < expected.numElements; ++i) {
        expect(Math.abs(actual.values[i] - expected.values[i]))
            .toBeLessThanOrEqual(tolerance);
    }
}

function expectVectorClose(actual: Vector, expected: readonly number[],
    tolerance: number = 1e-12): void {
    expect(actual.size).toBe(expected.length);
    for (let i = 0; i < expected.length; ++i) {
        expect(Math.abs(actual.values[i] - (expected[i] as number)))
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

// Rotation about the z-axis by 'angle' as a 4x4 matrix (MAT_VEC convention).
function rotateZ(angle: number): Matrix {
    const cs = Math.cos(angle);
    const sn = Math.sin(angle);
    return Matrix.fromArray(4, 4, [
        cs, -sn, 0, 0,
        sn, cs, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    ]);
}

function randomRotation4x4(rand: () => number): Matrix {
    const axis = Vector.fromArray([rand(), rand(), rand(), 0]);
    if (normalize(axis) === 0) {
        axis.values[0] = 1;
    }
    return Rotation.fromAxisAngle(new AxisAngle(axis, 3 * rand())).toMatrix();
}

describe('Transform', () => {
    it('default-constructs the identity transformation', () => {
        const t = new Transform();
        expect(t.isIdentity()).toBe(true);
        expect(t.isRSMatrix()).toBe(true);
        expect(t.isUniformScale()).toBe(true);
        expectMatrixClose(t.getHMatrix(), Matrix.identity(4, 4));
        expectMatrixClose(t.getHInverse(), Matrix.identity(4, 4));
        expectMatrixClose(t.getMatrix(), Matrix.identity(4, 4));
        expectVectorClose(t.getTranslation(), [0, 0, 0]);
        expectVectorClose(t.getScale(), [1, 1, 1]);
        expect(t.getUniformScale()).toBe(1);
        expect(t.getNorm()).toBe(1);
        expectMatrixClose(Transform.identity().getHMatrix(),
            Matrix.identity(4, 4));
    });

    it('builds the homogeneous matrix for a pure translation', () => {
        const t = new Transform();
        t.setTranslation(2, -3, 5);
        expect(t.isIdentity()).toBe(false);
        expect(t.isRSMatrix()).toBe(true);
        expectMatrixClose(t.getHMatrix(), Matrix.fromArray(4, 4, [
            1, 0, 0, 2,
            0, 1, 0, -3,
            0, 0, 1, 5,
            0, 0, 0, 1
        ]));
        expectMatrixClose(t.getHInverse(), Matrix.fromArray(4, 4, [
            1, 0, 0, -2,
            0, 1, 0, 3,
            0, 0, 1, -5,
            0, 0, 0, 1
        ]));
        expectVectorClose(t.getTranslationW0(), [2, -3, 5, 0]);
        expectVectorClose(t.getTranslationW1(), [2, -3, 5, 1]);

        // The Vector overloads agree with the scalar overload.
        const t3 = new Transform();
        t3.setTranslation(Vector.fromArray([2, -3, 5]));
        expectMatrixClose(t3.getHMatrix(), t.getHMatrix());
        const t4 = new Transform();
        t4.setTranslation(Vector.fromArray([2, -3, 5, 1]));
        expectMatrixClose(t4.getHMatrix(), t.getHMatrix());
    });

    it('composes rotation, uniform scale and translation by hand', () => {
        // A 90-degree rotation about z, uniform scale 2, translation
        // (1,2,3). H = {{2*Rz, T},{0,1}} with Rz = [[0,-1,0],[1,0,0],[0,0,1]].
        const t = new Transform();
        t.setRotation(rotateZ(Math.PI / 2));
        t.setUniformScale(2);
        t.setTranslation(1, 2, 3);

        expect(t.isUniformScale()).toBe(true);
        expect(t.getUniformScale()).toBe(2);
        expect(t.getNorm()).toBe(2);
        expectMatrixClose(t.getHMatrix(), Matrix.fromArray(4, 4, [
            0, -2, 0, 1,
            2, 0, 0, 2,
            0, 0, 2, 3,
            0, 0, 0, 1
        ]), 1e-15);

        // The forward transform of (1,0,0) is 2*(0,1,0) + (1,2,3).
        const y = mulTransform(t, Vector.fromArray([1, 0, 0, 1]));
        expectVectorClose(y, [1, 4, 3, 1], 1e-15);

        // The inverse transform takes it back.
        const x = mulMatrix(t.getHInverse(), y);
        expectVectorClose(x, [1, 0, 0, 1], 1e-15);
    });

    it('applies a nonuniform scale on the right of the rotation', () => {
        // M = R*S scales column j of R by s_j.
        const t = new Transform();
        t.setRotation(rotateZ(Math.PI / 2));
        t.setScale(2, 3, 4);
        t.setTranslation(0, 0, 0);

        expect(t.isRSMatrix()).toBe(true);
        expect(t.isUniformScale()).toBe(false);
        expect(t.getNorm()).toBe(4);
        expectVectorClose(t.getScale(), [2, 3, 4]);
        expectVectorClose(t.getScaleW1(), [2, 3, 4, 1]);
        expectMatrixClose(t.getHMatrix(), Matrix.fromArray(4, 4, [
            0, -3, 0, 0,
            2, 0, 0, 0,
            0, 0, 4, 0,
            0, 0, 0, 1
        ]), 1e-15);

        // H^{-1} = S^{-1}*R^t.
        expectMatrixClose(multiplyAB(t.getHInverse(), t.getHMatrix()),
            Matrix.identity(4, 4), 1e-15);

        // The scale vector overload agrees.
        const t2 = new Transform();
        t2.setRotation(rotateZ(Math.PI / 2));
        t2.setScale(Vector.fromArray([2, 3, 4]));
        expectMatrixClose(t2.getHMatrix(), t.getHMatrix());
    });

    it('supports a general (non-RS) matrix channel', () => {
        const m = Matrix.fromArray(4, 4, [
            1, 2, 0, 0,
            0, 1, 3, 0,
            4, 0, 1, 0,
            0, 0, 0, 1
        ]);
        const t = new Transform();
        t.setMatrix(m);
        t.setTranslation(1, -1, 2);

        expect(t.isRSMatrix()).toBe(false);
        expect(t.isUniformScale()).toBe(false);
        expectMatrixClose(t.getMatrix(), m);
        // Max-row-sum norm: rows sum to 3, 4, 5.
        expect(t.getNorm()).toBe(5);
        expectMatrixClose(t.getHMatrix(), Matrix.fromArray(4, 4, [
            1, 2, 0, 1,
            0, 1, 3, -1,
            4, 0, 1, 2,
            0, 0, 0, 1
        ]));
        expectMatrixClose(t.getHInverse(),
            inverse4x4(t.getHMatrix()).inverse, 1e-14);
        expectMatrixClose(multiplyAB(t.getHInverse(), t.getHMatrix()),
            Matrix.identity(4, 4), 1e-14);

        // The rotation/scale accessors are unavailable for a general matrix.
        expect(() => t.getRotation()).toThrow();
        expect(() => t.getScale()).toThrow();
        expect(() => t.getUniformScale()).toThrow();
        expect(() => t.setScale(1, 2, 3)).toThrow();
        expect(() => t.setUniformScale(2)).toThrow();
        expect(() => t.makeUnitScale()).toThrow();
    });

    it('rejects zero scales and validates matrix sizes', () => {
        const t = new Transform();
        expect(() => t.setScale(0, 1, 1)).toThrow();
        expect(() => t.setScale(1, 0, 1)).toThrow();
        expect(() => t.setScale(1, 1, 0)).toThrow();
        expect(() => t.setUniformScale(0)).toThrow();
        expect(() => t.setMatrix(Matrix.identity(3, 3))).toThrow();
        expect(() => t.setRotation(Matrix.identity(2, 2))).toThrow();
    });

    it('makes identity and unit scale', () => {
        const t = new Transform();
        t.setRotation(rotateZ(0.3));
        t.setScale(2, 3, 4);
        t.setTranslation(5, 6, 7);

        t.makeUnitScale();
        expect(t.isUniformScale()).toBe(true);
        expectVectorClose(t.getScale(), [1, 1, 1]);
        expectMatrixClose(t.getHMatrix(), (() => {
            const h = rotateZ(0.3);
            h.set(0, 3, 5);
            h.set(1, 3, 6);
            h.set(2, 3, 7);
            return h;
        })());

        t.makeIdentity();
        expect(t.isIdentity()).toBe(true);
        expectMatrixClose(t.getHMatrix(), Matrix.identity(4, 4));
        expectMatrixClose(t.getHInverse(), Matrix.identity(4, 4));
    });

    it('sets the rotation from 3x3, quaternion, axis-angle, Euler', () => {
        const angle = 0.7;
        const axis = normalizedAxis();
        const axis4 = Vector.fromArray([axis.values[0], axis.values[1],
            axis.values[2], 0]);
        const expected = Rotation.fromAxisAngle(
            new AxisAngle(axis4, angle)).toMatrix();

        // 4x4 matrix.
        const t4 = new Transform();
        t4.setRotation(expected);
        expectMatrixClose(t4.getHMatrix(), expected, 1e-15);

        // 3x3 matrix.
        const r3 = new Matrix(3, 3);
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 3; ++c) {
                r3.set(r, c, expected.get(r, c));
            }
        }
        const t3 = new Transform();
        t3.setRotation(r3);
        expectMatrixClose(t3.getHMatrix(), expected, 1e-15);
        expectMatrixClose(t3.getRotationMatrix3x3(), r3, 1e-15);

        // Quaternion.
        const q = Rotation.fromMatrix(expected).toQuaternion();
        const tq = new Transform();
        tq.setRotation(q);
        expectMatrixClose(tq.getHMatrix(), expected, 1e-14);
        const q2 = tq.getRotationQuaternion();
        // The quaternion round-trip is unique up to sign.
        const sign = q.values[3] * q2.values[3] >= 0 ? 1 : -1;
        for (let i = 0; i < 4; ++i) {
            expect(Math.abs(q.values[i] - sign * q2.values[i]))
                .toBeLessThanOrEqual(1e-14);
        }

        // Axis-angle with a 3-tuple axis.
        const taa3 = new Transform();
        taa3.setRotation(new AxisAngle(axis, angle));
        expectMatrixClose(taa3.getHMatrix(), expected, 1e-15);
        const aa3 = taa3.getRotationAxisAngle3();
        expect(aa3.axis.size).toBe(3);
        expectMatrixClose(
            Rotation.fromAxisAngle(new AxisAngle(aa3.axis, aa3.angle))
                .toMatrix(),
            Matrix.fromArray(3, 3, [
                expected.get(0, 0), expected.get(0, 1), expected.get(0, 2),
                expected.get(1, 0), expected.get(1, 1), expected.get(1, 2),
                expected.get(2, 0), expected.get(2, 1), expected.get(2, 2)
            ]), 1e-14);

        // Axis-angle with a 4-tuple axis.
        const taa4 = new Transform();
        taa4.setRotation(new AxisAngle(axis4, angle));
        expectMatrixClose(taa4.getHMatrix(), expected, 1e-15);
        const aa4 = taa4.getRotationAxisAngle4();
        expect(aa4.axis.size).toBe(4);

        // Euler angles.
        const euler = new EulerAngles(0, 1, 2, 0.3, -0.4, 1.1);
        const te = new Transform();
        te.setRotation(euler);
        expectMatrixClose(te.getHMatrix(),
            Rotation.fromEulerAngles(euler, 4).toMatrix(), 1e-15);
        const back = te.getRotationEulerAngles(0, 1, 2);
        for (let i = 0; i < 3; ++i) {
            expect(Math.abs(back.angle[i] - euler.angle[i]))
                .toBeLessThanOrEqual(1e-14);
        }
    });

    it('inverts a rotation-uniform-scale-translation transform', () => {
        const t = new Transform();
        t.setRotation(rotateZ(0.9));
        t.setUniformScale(3);
        t.setTranslation(1, -2, 4);

        const inv = t.inverse();
        expect(inv.isRSMatrix()).toBe(true);
        expect(inv.isUniformScale()).toBe(true);
        expect(Math.abs(inv.getUniformScale() - 1 / 3))
            .toBeLessThanOrEqual(1e-15);
        expectMatrixClose(inv.getRotation(), transpose(t.getRotation()),
            1e-15);

        expectMatrixClose(multiplyAB(inv.getHMatrix(), t.getHMatrix()),
            Matrix.identity(4, 4), 1e-14);
        expectMatrixClose(inv.getHMatrix(), t.getHInverse(), 1e-14);
    });

    it('inverts nonuniform-scale and general transforms', () => {
        const t = new Transform();
        t.setRotation(rotateZ(0.4));
        t.setScale(2, -3, 5);
        t.setTranslation(7, 8, 9);

        const inv = t.inverse();
        expect(inv.isRSMatrix()).toBe(false);
        expectMatrixClose(multiplyAB(inv.getHMatrix(), t.getHMatrix()),
            Matrix.identity(4, 4), 1e-14);

        // The fixed 'M' channel of the inverse has a zero last column and
        // row; only the translation channel carries the translation.
        const m = inv.getMatrix();
        for (let i = 0; i < 3; ++i) {
            expect(m.get(i, 3)).toBe(0);
            expect(m.get(3, i)).toBe(0);
        }
        expect(m.get(3, 3)).toBe(1);
        expectVectorClose(inv.getTranslation(),
            [inv.getHMatrix().get(0, 3), inv.getHMatrix().get(1, 3),
                inv.getHMatrix().get(2, 3)]);

        // The identity inverts to the identity.
        expect(new Transform().inverse().isIdentity()).toBe(true);
    });

    it('multiplies a transform with vectors and 4x4 matrices', () => {
        const t = new Transform();
        t.setRotation(rotateZ(0.25));
        t.setScale(1, 2, 3);
        t.setTranslation(-1, 0, 4);
        const h = t.getHMatrix();

        const v = Vector.fromArray([0.5, -1.5, 2.5, 1]);
        expectVectorClose(mulTransform(t, v), mulMatrix(h, v).values, 1e-15);
        expectVectorClose(mulTransform(v, t), mulMatrix(v, h).values, 1e-15);

        const a = Matrix.fromArray(4, 4, [
            1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7
        ]);
        expectMatrixClose(mulTransform(a, t), mulMatrix(a, h), 1e-13);
        expectMatrixClose(mulTransform(t, a), mulMatrix(h, a), 1e-13);
    });

    it('short-circuits products with the identity transform', () => {
        const t = new Transform();
        t.setRotation(rotateZ(0.6));
        t.setTranslation(1, 2, 3);
        const id = new Transform();

        const p0 = mulTransform(id, t);
        expectMatrixClose(p0.getHMatrix(), t.getHMatrix());
        const p1 = mulTransform(t, id);
        expectMatrixClose(p1.getHMatrix(), t.getHMatrix());

        // The results are copies, not aliases.
        p0.setTranslation(9, 9, 9);
        expectVectorClose(t.getTranslation(), [1, 2, 3]);
    });

    it('composes transforms through the RS fast path', () => {
        // A has a uniform scale, so A*B stays an RS-matrix.
        const a = new Transform();
        a.setRotation(rotateZ(Math.PI / 2));
        a.setUniformScale(2);
        a.setTranslation(1, 0, 0);

        const b = new Transform();
        b.setRotation(rotateZ(Math.PI / 2));
        b.setUniformScale(3);
        b.setTranslation(0, 1, 0);

        const p = mulTransform(a, b);
        expect(p.isRSMatrix()).toBe(true);
        expect(p.isUniformScale()).toBe(true);
        expect(Math.abs(p.getUniformScale() - 6)).toBeLessThanOrEqual(1e-15);
        expectMatrixClose(p.getRotation(), rotateZ(Math.PI), 1e-15);
        // T = 2*Rz(90)*(0,1,0) + (1,0,0) = 2*(-1,0,0) + (1,0,0) = (-1,0,0).
        expectVectorClose(p.getTranslation(), [-1, 0, 0], 1e-15);
        expectMatrixClose(p.getHMatrix(),
            multiplyAB(a.getHMatrix(), b.getHMatrix()), 1e-14);

        // A uniform, B nonuniform: still RS, with scale a*S_B.
        const b2 = new Transform();
        b2.setRotation(rotateZ(0.3));
        b2.setScale(1, 2, 3);
        b2.setTranslation(4, 5, 6);
        const p2 = mulTransform(a, b2);
        expect(p2.isRSMatrix()).toBe(true);
        expect(p2.isUniformScale()).toBe(false);
        expectVectorClose(p2.getScale(), [2, 4, 6], 1e-15);
        expectMatrixClose(p2.getHMatrix(),
            multiplyAB(a.getHMatrix(), b2.getHMatrix()), 1e-14);
    });

    it('composes transforms through the general path', () => {
        // A has a nonuniform scale, so the product is a general matrix.
        const a = new Transform();
        a.setRotation(rotateZ(0.35));
        a.setScale(1, 2, 3);
        a.setTranslation(1, 2, 3);

        const b = new Transform();
        b.setRotation(rotateZ(-0.8));
        b.setScale(4, 5, 6);
        b.setTranslation(-1, 0, 2);

        const p = mulTransform(a, b);
        expect(p.isRSMatrix()).toBe(false);
        expectMatrixClose(p.getHMatrix(),
            multiplyAB(a.getHMatrix(), b.getHMatrix()), 1e-14);

        // With a general (non-RS) A.
        const g = new Transform();
        g.setMatrix(Matrix.fromArray(4, 4, [
            1, 2, 0, 0, 0, 1, 3, 0, 4, 0, 1, 0, 0, 0, 0, 1
        ]));
        g.setTranslation(0, -1, 1);
        const pg = mulTransform(g, b);
        expectMatrixClose(pg.getHMatrix(),
            multiplyAB(g.getHMatrix(), b.getHMatrix()), 1e-14);
        const pg2 = mulTransform(b, g);
        expectMatrixClose(pg2.getHMatrix(),
            multiplyAB(b.getHMatrix(), g.getHMatrix()), 1e-14);
    });

    it('composes a three-level hierarchy associatively', () => {
        const a = new Transform();
        a.setRotation(rotateZ(0.2));
        a.setUniformScale(1.5);
        a.setTranslation(1, 0, -1);

        const b = new Transform();
        b.setRotation(rotateZ(-1.1));
        b.setScale(2, 0.5, 3);
        b.setTranslation(0, 2, 0);

        const c = new Transform();
        c.setRotation(rotateZ(0.75));
        c.setUniformScale(0.25);
        c.setTranslation(3, -2, 1);

        const left = mulTransform(mulTransform(a, b), c);
        const right = mulTransform(a, mulTransform(b, c));
        expectMatrixClose(left.getHMatrix(), right.getHMatrix(), 1e-13);
        expectMatrixClose(left.getHMatrix(), multiplyAB(a.getHMatrix(),
            multiplyAB(b.getHMatrix(), c.getHMatrix())), 1e-13);

        // A world point mapped through the chain.
        const p = Vector.fromArray([0.3, -0.7, 1.9, 1]);
        expectVectorClose(mulTransform(left, p),
            mulTransform(a, mulTransform(b, mulTransform(c, p))).values,
            1e-13);
    });

    it('keeps getters from aliasing internal state', () => {
        const t = new Transform();
        t.setRotation(rotateZ(0.5));
        t.setTranslation(1, 2, 3);
        const h = t.getHMatrix();
        h.set(0, 0, 1234);
        expect(t.getHMatrix().get(0, 0)).not.toBe(1234);
        const m = t.getMatrix();
        m.set(1, 1, 4321);
        expect(t.getMatrix().get(1, 1)).not.toBe(4321);
        const hi = t.getHInverse();
        hi.set(2, 2, 999);
        expect(t.getHInverse().get(2, 2)).not.toBe(999);
    });

    it('recomputes the cached inverse after each channel change', () => {
        const t = new Transform();
        t.setRotation(rotateZ(0.5));
        expectMatrixClose(multiplyAB(t.getHInverse(), t.getHMatrix()),
            Matrix.identity(4, 4), 1e-15);
        t.setUniformScale(4);
        expectMatrixClose(multiplyAB(t.getHInverse(), t.getHMatrix()),
            Matrix.identity(4, 4), 1e-15);
        t.setTranslation(1, 2, 3);
        expectMatrixClose(multiplyAB(t.getHInverse(), t.getHMatrix()),
            Matrix.identity(4, 4), 1e-14);
        t.setScale(1, 2, 3);
        expectMatrixClose(multiplyAB(t.getHInverse(), t.getHMatrix()),
            Matrix.identity(4, 4), 1e-14);
        t.setMatrix(Matrix.fromArray(4, 4, [
            2, 1, 0, 0, 0, 3, 1, 0, 1, 0, 4, 0, 0, 0, 0, 1
        ]));
        expectMatrixClose(multiplyAB(t.getHInverse(), t.getHMatrix()),
            Matrix.identity(4, 4), 1e-14);
        t.makeIdentity();
        expectMatrixClose(t.getHInverse(), Matrix.identity(4, 4));
    });

    it('restores a valid last row after a singular general matrix', () => {
        // Port fix: the RS branch of getHInverse writes the last row.
        const t = new Transform();
        t.setMatrix(Matrix.fromArray(4, 4, [
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        ]));
        const singular = t.getHInverse();
        expect(singular.get(3, 3)).toBe(0);

        t.setRotation(rotateZ(0.3));
        t.setUniformScale(2);
        t.setTranslation(1, 2, 3);
        const hinv = t.getHInverse();
        expect(hinv.get(3, 0)).toBe(0);
        expect(hinv.get(3, 1)).toBe(0);
        expect(hinv.get(3, 2)).toBe(0);
        expect(hinv.get(3, 3)).toBe(1);
        expectMatrixClose(multiplyAB(hinv, t.getHMatrix()),
            Matrix.identity(4, 4), 1e-14);
    });

    it('matches independent matrix algebra over random transforms', () => {
        const rand = makeRandom(20260901);
        for (let trial = 0; trial < 200; ++trial) {
            const a = randomTransform(rand, trial % 3);
            const b = randomTransform(rand, (trial + 1) % 3);

            const ha = a.getHMatrix();
            const hb = b.getHMatrix();

            // The cached inverse agrees with the general 4x4 inverse.
            expectMatrixClose(multiplyAB(a.getHInverse(), ha),
                Matrix.identity(4, 4), 1e-9);
            expectMatrixClose(a.getHInverse(), inverse4x4(ha).inverse, 1e-9);

            // The transform product agrees with the matrix product.
            const p = mulTransform(a, b);
            expectMatrixClose(p.getHMatrix(), multiplyAB(ha, hb), 1e-9);

            // The transform inverse agrees with the cached inverse.
            expectMatrixClose(a.inverse().getHMatrix(), a.getHInverse(),
                1e-9);

            // Forward then inverse maps a point back to itself.
            const x = Vector.fromArray([rand(), rand(), rand(), 1]);
            const y = mulTransform(a, x);
            expectVectorClose(mulMatrix(a.getHInverse(), y), x.values, 1e-9);

            // getNorm is an upper bound on the scale for RS transforms.
            if (a.isRSMatrix()) {
                const s = a.getScale();
                expect(a.getNorm()).toBeGreaterThanOrEqual(
                    Math.max(Math.abs(s.values[0]),
                        Math.abs(s.values[1]), Math.abs(s.values[2]))
                    - 1e-15);
            }
        }
    });
});

function normalizedAxis(): Vector {
    const axis = Vector.fromArray([1, -2, 3]);
    normalize(axis);
    return axis;
}

// kind 0: rotation + uniform scale; 1: rotation + nonuniform scale;
// 2: general invertible matrix.
function randomTransform(rand: () => number, kind: number): Transform {
    const t = new Transform();
    if (kind === 2) {
        // A rotation times a nonsingular diagonal plus a small perturbation
        // keeps the matrix well conditioned.
        const r = randomRotation4x4(rand);
        const m = Matrix.identity(4, 4);
        for (let i = 0; i < 3; ++i) {
            for (let j = 0; j < 3; ++j) {
                m.set(i, j, r.get(i, j) * (1 + 0.5 * rand())
                    + 0.1 * rand());
            }
        }
        t.setMatrix(m);
    } else {
        t.setRotation(randomRotation4x4(rand));
        if (kind === 0) {
            t.setUniformScale(1 + 0.5 * rand());
        } else {
            t.setScale(1 + 0.5 * rand(), 1 + 0.5 * rand(),
                1 + 0.5 * rand());
        }
    }
    t.setTranslation(4 * rand(), 4 * rand(), 4 * rand());
    return t;
}

// ---------------------------------------------------------------------------
// Verification wave (V02): property-based re-checks against Transform.h.
// ---------------------------------------------------------------------------

const arbRotation4 = () => fc.tuple(unitVector(3), finite(-GTE_C_PI, GTE_C_PI))
    .map(([axis, angle]) => Rotation.fromAxisAngle(
        new AxisAngle(hlift(axis, 0), angle)).toMatrix());

const arbNonzeroScale = () => finite(-4, 4).filter(s => Math.abs(s) > 0.1);

// Transforms of each structural kind: identity, rotation with uniform scale,
// rotation with nonuniform scale, and a general (non-RS) matrix.
const arbTransform = () => fc.oneof(
    fc.constant(null).map(() => new Transform()),
    fc.tuple(arbRotation4(), arbNonzeroScale(), vector(3, -5, 5))
        .map(([R, s, t]) => {
            const x = new Transform();
            x.setRotation(R);
            x.setUniformScale(s);
            x.setTranslation(t);
            return x;
        }),
    fc.tuple(arbRotation4(), arbNonzeroScale(), arbNonzeroScale(),
        arbNonzeroScale(), vector(3, -5, 5))
        .map(([R, s0, s1, s2, t]) => {
            const x = new Transform();
            x.setRotation(R);
            x.setScale(s0, s1, s2);
            x.setTranslation(t);
            return x;
        }),
    fc.tuple(invertibleMatrix(3), vector(3, -5, 5)).map(([M3, t]) => {
        const M = Matrix.identity(4, 4);
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 3; ++c) {
                M.set(r, c, M3.get(r, c));
            }
        }
        const x = new Transform();
        x.setMatrix(M);
        x.setTranslation(t);
        return x;
    }));

describe('Transform verification', () => {
    it('the homogeneous matrix is {{M*S, T}, {0, 1}}', () => {
        check(arbTransform(), x => {
            const H = x.getHMatrix();
            // The affine structure is invariant.
            expectClose(H.get(3, 0), 0, 0, 0);
            expectClose(H.get(3, 1), 0, 0, 0);
            expectClose(H.get(3, 2), 0, 0, 0);
            expectClose(H.get(3, 3), 1, 0, 0);
            expectVectorsClose(x.getTranslation(),
                Vector.fromArray([H.get(0, 3), H.get(1, 3), H.get(2, 3)]),
                0, 0);

            const M = x.getMatrix();
            for (let r = 0; r < 3; ++r) {
                for (let c = 0; c < 3; ++c) {
                    // For an RS transform column c of M is scaled by S[c]
                    // (M = R*S under MAT_VEC); otherwise M is copied.
                    const s = (x.isRSMatrix() ? x.getScale().get(c) : 1);
                    expectClose(H.get(r, c), M.get(r, c) * s, 1e-12, 1e-12);
                }
            }
        });
    });

    it('getHInverse inverts the homogeneous matrix and equals the inverse '
        + 'transform', () => {
            check(arbTransform(), x => {
                const H = x.getHMatrix();
                const Hinv = x.getHInverse();
                expectMatrixClose(multiplyAB(Hinv, H), Matrix.identity(4, 4),
                    1e-8);
                expectMatrixClose(multiplyAB(H, Hinv), Matrix.identity(4, 4),
                    1e-8);
                // The last row is always (0,0,0,1): the inverse of an affine
                // transform is affine. The identity and rotation-scale
                // branches write those entries directly, so they are exact;
                // the general branch obtains them from the cofactor formula
                // of inverse4x4, which is exact only to round-off (and can
                // produce a signed zero).
                const exact = x.isIdentity() || x.isRSMatrix();
                const tol = (exact ? 0 : 1e-12);
                expectClose(Hinv.get(3, 0), 0, tol, 0);
                expectClose(Hinv.get(3, 1), 0, tol, 0);
                expectClose(Hinv.get(3, 2), 0, tol, 0);
                expectClose(Hinv.get(3, 3), 1, tol, tol);

                expectMatrixClose(x.inverse().getHMatrix(), Hinv, 1e-8);
            }, 60);
        });

    it('the inverse transform undoes the forward transform on points', () => {
        check(fc.tuple(arbTransform(), vector(3, -5, 5)), ([x, p]) => {
            const P = hlift(p, 1);
            const forward = mulTransform(x, P);
            const back = mulTransform(x.inverse(), forward);
            expectVectorsClose(back, P, 1e-7, 1e-7);

            // Composing a transform with its inverse gives the identity.
            const composed = mulTransform(x, x.inverse());
            expectMatrixClose(composed.getHMatrix(), Matrix.identity(4, 4),
                1e-7);
        });
    });

    it('composition of transforms is composition of homogeneous matrices',
        () => {
            check(fc.tuple(arbTransform(), arbTransform()), ([a, b]) => {
                const product = mulTransform(a, b);
                expectMatrixClose(product.getHMatrix(),
                    multiplyAB(a.getHMatrix(), b.getHMatrix()), 1e-9);

                // The RS fast path (A uniform scale, both RS) must agree with
                // the general path, and must keep the RS channels.
                if (a.isUniformScale() && a.isRSMatrix() && b.isRSMatrix()
                    && !a.isIdentity() && !b.isIdentity()) {
                    expect(product.isRSMatrix()).toBe(true);
                }
            });
        });

    it('composition is associative through the homogeneous matrices', () => {
        check(fc.tuple(arbTransform(), arbTransform(), arbTransform()),
            ([a, b, c]) => {
                const left = mulTransform(mulTransform(a, b), c);
                const right = mulTransform(a, mulTransform(b, c));
                expectMatrixClose(left.getHMatrix(), right.getHMatrix(), 1e-8);
            });
    });

    it('the mixed products with vectors and 4x4 matrices use the homogeneous '
        + 'matrix', () => {
            check(fc.tuple(arbTransform(), vector(4, -5, 5), matrix(4, 4)),
                ([x, v, M]) => {
                    const H = x.getHMatrix();
                    expectVectorsClose(mulTransform(x, v),
                        mulMatrix(H, v) as Vector, 1e-12, 1e-12);
                    expectVectorsClose(mulTransform(v, x),
                        mulMatrix(v, H) as Vector, 1e-12, 1e-12);
                    expectMatrixClose(mulTransform(x, M) as Matrix,
                        multiplyAB(H, M), 1e-12);
                    expectMatrixClose(mulTransform(M, x) as Matrix,
                        multiplyAB(M, H), 1e-12);
                }, 80);
        });

    it('a product with the identity transform returns a copy of the other '
        + 'operand', () => {
            check(arbTransform(), x => {
                const id = new Transform();
                for (const product of [mulTransform(id, x),
                    mulTransform(x, id)]) {
                    expectMatrixClose(product.getHMatrix(), x.getHMatrix(), 0);
                    expect(product.isRSMatrix()).toBe(x.isRSMatrix());
                    // A copy, not the same object: mutating the result must
                    // not touch the operand.
                    product.setTranslation(99, 99, 99);
                    expect(x.getTranslation().get(0)).not.toBe(99);
                }
            }, 100);
        });

    it('every way of setting the rotation produces the same transform', () => {
        check(fc.tuple(unitVector(3), finite(-GTE_C_PI + 1e-3,
            GTE_C_PI - 1e-3), vector(3, -5, 5)), ([axis, angle, t]) => {
                const R4 = Rotation.fromAxisAngle(
                    new AxisAngle(hlift(axis, 0), angle)).toMatrix();
                const R3 = new Matrix(3, 3);
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 3; ++c) {
                        R3.set(r, c, R4.get(r, c));
                    }
                }
                const q = Rotation.fromMatrix(R4).toQuaternion();
                const aa3 = new AxisAngle(axis, angle);
                const aa4 = new AxisAngle(hlift(axis, 0), angle);
                const euler = Rotation.fromMatrix(R4).toEulerAngles(0, 1, 2);

                const reference = new Transform();
                reference.setRotation(R4);
                reference.setTranslation(t);
                const H = reference.getHMatrix();

                for (const rotation of [R3, q, aa3, aa4, euler]) {
                    const x = new Transform();
                    x.setRotation(rotation);
                    x.setTranslation(t);
                    expect(x.isRSMatrix()).toBe(true);
                    expect(x.isIdentity()).toBe(false);
                    expectMatrixClose(x.getHMatrix(), H, 1e-7);
                }
            }, 60);
    });

    it('the rotation getters recover the rotation that was set', () => {
        check(fc.tuple(unitVector(3), finite(1e-2, GTE_C_PI - 1e-2)),
            ([axis, angle]) => {
                const R4 = Rotation.fromAxisAngle(
                    new AxisAngle(hlift(axis, 0), angle)).toMatrix();
                const x = new Transform();
                x.setRotation(R4);

                expectMatrixClose(x.getRotation(), R4, 0);

                const R3 = x.getRotationMatrix3x3();
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 3; ++c) {
                        expectClose(R3.get(r, c), R4.get(r, c), 0, 0);
                    }
                }

                const y = new Transform();
                y.setRotation(x.getRotationQuaternion());
                expectMatrixClose(y.getHMatrix(), x.getHMatrix(), 1e-8);

                const z = new Transform();
                z.setRotation(x.getRotationAxisAngle3());
                expectMatrixClose(z.getHMatrix(), x.getHMatrix(), 1e-7);

                const w = new Transform();
                w.setRotation(x.getRotationAxisAngle4());
                expectMatrixClose(w.getHMatrix(), x.getHMatrix(), 1e-7);

                const e = new Transform();
                e.setRotation(x.getRotationEulerAngles(2, 0, 1));
                expectMatrixClose(e.getHMatrix(), x.getHMatrix(), 1e-7);
            }, 60);
    });

    it('getNorm is the largest |scale| for RS transforms and the max row sum '
        + 'otherwise', () => {
            check(arbTransform(), x => {
                if (x.isRSMatrix()) {
                    const s = x.getScale();
                    expectClose(x.getNorm(), Math.max(Math.abs(s.get(0)),
                        Math.abs(s.get(1)), Math.abs(s.get(2))), 0, 0);
                } else {
                    const M = x.getMatrix();
                    let expected = 0;
                    for (let r = 0; r < 3; ++r) {
                        expected = Math.max(expected, Math.abs(M.get(r, 0))
                            + Math.abs(M.get(r, 1)) + Math.abs(M.get(r, 2)));
                    }
                    expectClose(x.getNorm(), expected, 0, 0);
                }
            });
        });

    it('the structural flags follow the documented setter rules', () => {
        check(fc.tuple(arbRotation4(), arbNonzeroScale(), vector(3, -5, 5)),
            ([R, s, t]) => {
                const x = new Transform();
                expect(x.isIdentity()).toBe(true);
                expect(x.isRSMatrix()).toBe(true);
                expect(x.isUniformScale()).toBe(true);

                x.setTranslation(t);
                expect(x.isIdentity()).toBe(false);
                expect(x.isRSMatrix()).toBe(true);

                x.setRotation(R);
                expect(x.isRSMatrix()).toBe(true);

                x.setScale(s, s, s);
                // setScale always clears the uniform-scale hint, even when
                // the three scales are equal.
                expect(x.isUniformScale()).toBe(false);

                x.setUniformScale(s);
                expect(x.isUniformScale()).toBe(true);

                x.makeUnitScale();
                expect(x.isUniformScale()).toBe(true);
                expectVectorsClose(x.getScale(), Vector.fromArray([1, 1, 1]),
                    0, 0);

                const general = Matrix.identity(4, 4);
                general.set(0, 1, 2);
                x.setMatrix(general);
                expect(x.isRSMatrix()).toBe(false);
                expect(x.isUniformScale()).toBe(false);
                // The rotation-scale getters are unavailable in this state.
                expect(() => x.getRotation()).toThrow();
                expect(() => x.getScale()).toThrow();
                expect(() => x.getUniformScale()).toThrow();
                expect(() => x.makeUnitScale()).toThrow();
                expect(() => x.setScale(1, 1, 1)).toThrow();
                expect(() => x.setUniformScale(1)).toThrow();

                x.makeIdentity();
                expect(x.isIdentity()).toBe(true);
                expect(x.isRSMatrix()).toBe(true);
                expectMatrixClose(x.getHMatrix(), Matrix.identity(4, 4), 0);
            }, 60);
    });

    it('no accessor aliases the internal state', () => {
        check(arbTransform(), x => {
            const before = x.getHMatrix();

            x.getHMatrix().set(0, 0, 12345);
            x.getMatrix().set(0, 0, 12345);
            x.getHInverse().set(0, 0, 12345);
            x.getTranslation().set(0, 12345);
            x.getTranslationW0().set(0, 12345);
            x.getTranslationW1().set(0, 12345);
            if (x.isRSMatrix()) {
                x.getRotation().set(0, 0, 12345);
                x.getScale().set(0, 12345);
                x.getScaleW1().set(0, 12345);
            }

            expectMatrixClose(x.getHMatrix(), before, 0);

            // clone() is a deep copy.
            const copy = x.clone();
            copy.setTranslation(99, 99, 99);
            expect(x.getTranslation().get(0)).not.toBe(99);
        });
    });
});
