// Replays oracle/cpp/cases/v02-algebra.cpp. Keep the two files in the same
// order. Cases whose upstream path is arithmetic only are declared exact;
// cases that call sin, cos, acos, asin or atan2 carry the default 1e-12
// tolerance and say so. Where a conversion mixes the two, the arithmetic-only
// outputs use io.outRealExact so they are still required to be bit-identical.
import { describe } from 'vitest';
import { AxisAngle } from '../../src/AxisAngle.js';
import { ConvertCoordinates } from '../../src/ConvertCoordinates.js';
import { EulerAngles } from '../../src/EulerAngles.js';
import { Hyperellipsoid } from '../../src/Hyperellipsoid.js';
import { Hyperplane } from '../../src/Hyperplane.js';
import { Line } from '../../src/Line.js';
import { Matrix } from '../../src/Matrix.js';
import {
    adjoint2x2, determinant2x2, doTransform2x2, getBasis2x2,
    getRotationAngle2x2, inverse2x2, makeRotation2x2, setBasis2x2, trace2x2
} from '../../src/Matrix2x2.js';
import {
    adjoint3x3, determinant3x3, doTransform3x3, getBasis3x3, inverse3x3,
    setBasis3x3, trace3x3
} from '../../src/Matrix3x3.js';
import {
    adjoint4x4, determinant4x4, doTransform4x4, getBasis4x4, inverse4x4,
    makeObliqueProjection4x4, makePerspectiveProjection4x4, makeReflection4x4,
    setBasis4x4, trace4x4
} from '../../src/Matrix4x4.js';
import { perspectiveProject, projectEllipse2, projectEllipsoid3 } from '../../src/Projection.js';
import {
    Quaternion, addQuaternion, conjugate, divQuaternion, inverseQuaternion,
    mulQuaternion, negateQuaternion, rotate, slerpQuaternion, subQuaternion
} from '../../src/Quaternion.js';
import { Rotation } from '../../src/Rotation.js';
import {
    getRotC0EstimateMaxError, getRotC1EstimateMaxError,
    getRotC2EstimateMaxError, getRotC3EstimateMaxError,
    getRotC4EstimateMaxError, rotC0Estimate, rotC1Estimate, rotC2Estimate,
    rotC3Estimate, rotC4Estimate, rotationAndDerivativeEstimate,
    rotationDerivativeEstimate, rotationEstimate
} from '../../src/RotationEstimate.js';
import { slerpUsingCosAngle, slerpUsingMidpoint } from '../../src/Slerp.js';
import { Transform, mulTransform } from '../../src/Transform.js';
import { Vector, dot, length, normalize } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---- input readers, mirroring the C++ generators ----

function mat(io: OracleIO, n: number): Matrix {
    const m = new Matrix(n, n);
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) { m.set(r, c, io.real()); }
    }
    return m;
}

function quat(io: OracleIO): Quaternion {
    return new Quaternion(io.real(), io.real(), io.real(), io.real());
}

// The port of the C++ Lift3To4 helper: a 3x3 rotation in the upper-left block
// of the 4x4 identity.
function lift3To4(r3: Matrix): Matrix {
    const r4 = Matrix.identity(4, 4);
    for (let r = 0; r < 3; ++r) {
        for (let c = 0; c < 3; ++c) { r4.set(r, c, r3.get(r, c)); }
    }
    return r4;
}

// The port of the C++ EulerAxes helper.
const EULER_AXES: readonly (readonly [number, number, number])[] = [
    [0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0],
    [0, 1, 0], [0, 2, 0], [1, 0, 1], [1, 2, 1], [2, 0, 2], [2, 1, 2],
    [1, 1, 2]
];

// ---- output writers, mirroring the C++ recorders ----

// The axis is arithmetic only (a difference of matrix entries or a quaternion
// component divided by a sqrt, then a Normalize), the angle comes from acos.
function outAxisAngle(io: OracleIO, a: AxisAngle): void {
    io.outVecExact(a.axis);
    io.outReal(a.angle);
}

// The same, for a conversion whose source is itself a libm result (a matrix
// or quaternion built from sin and cos): the axis is then no longer a
// function of recorded inputs alone, so it carries the case tolerance too.
function outAxisAngleTol(io: OracleIO, a: AxisAngle): void {
    io.outVec(a.axis);
    io.outReal(a.angle);
}

function outEulerAngles(io: OracleIO, e: EulerAngles): void {
    for (let i = 0; i < 3; ++i) { io.outInt(e.axis[i]); }
    for (let i = 0; i < 3; ++i) { io.outReal(e.angle[i]); }
    io.outInt(e.result);
}

// The channels that every Transform defines.
function outTransform(io: OracleIO, t: Transform): void {
    io.outBool(t.isIdentity());
    io.outBool(t.isRSMatrix());
    io.outBool(t.isUniformScale());
    io.outMat(t.getHMatrix());
    io.outMat(t.getHInverse());
    io.outReal(t.getNorm());
}

// The full accessor set. The rotation-only accessors are emitted only when
// the is-rsmatrix hint holds, since otherwise they throw.
function outTransformFull(io: OracleIO, t: Transform): void {
    outTransform(io, t);
    io.outMat(t.getMatrix());
    io.outVec(t.getTranslation());
    io.outVec(t.getTranslationW0());
    io.outVec(t.getTranslationW1());
    if (t.isRSMatrix()) {
        io.outMat(t.getRotation());
        io.outMat(t.getRotationMatrix3x3());
        io.outVec(t.getScale());
        io.outVec(t.getScaleW1());
        io.outVec(t.getRotationQuaternion());
    }
    if (t.isUniformScale()) {
        io.outReal(t.getUniformScale());
    }
}

// The port of the C++ ApplyTranslation helper.
function applyTranslation(io: OracleIO, t: Transform, tr: Vector): void {
    const which = io.index % 3;
    if (which === 0) {
        t.setTranslation(tr.values[0], tr.values[1], tr.values[2]);
    } else if (which === 1) {
        t.setTranslation(tr);
    } else {
        t.setTranslation(Vector.fromArray([tr.values[0], tr.values[1],
            tr.values[2], 1]));
    }
}

// The port of the C++ ApplyScale helper.
function applyScale(io: OracleIO, t: Transform, s0: number, s1: number,
    s2: number, us: number): void {
    const which = io.index % 4;
    if (which === 1) {
        t.setScale(s0, s1, s2);
    } else if (which === 2) {
        t.setUniformScale(us);
    } else if (which === 3) {
        t.setScale(Vector.fromArray([s0, s1, s2]));
        t.makeUnitScale();
    }
}

// The port of the C++ MakeTransform helper.
function makeTransform(kind: number, rot3: Matrix, general: Matrix,
    tr: Vector, s0: number, s1: number, s2: number, us: number): Transform {
    const t = new Transform();
    if (kind === 0) {
        return t;
    }
    if (kind === 4) {
        t.setMatrix(general);
    } else {
        t.setRotation(rot3);
    }
    t.setTranslation(tr.values[0], tr.values[1], tr.values[2]);
    if (kind === 2) {
        t.setUniformScale(us);
    } else if (kind === 3) {
        t.setScale(s0, s1, s2);
    }
    return t;
}

interface ProjectionSetup {
    ellipsoid: Hyperellipsoid;
    E: Vector;
    N: Vector;
    U: Vector;
    V: Vector;
    constant: number;
    n: number;
}

// The port of the C++ MakeProjectionSetup helper. Every value it needs is a
// recorded input, so nothing is recomputed here.
function projectionSetup(io: OracleIO): ProjectionSetup {
    const N = io.vec(3);
    const U = io.vec(3);
    const V = io.vec(3);
    const ellipsoid = new Hyperellipsoid(3);
    ellipsoid.axis[0] = io.vec(3);
    ellipsoid.axis[1] = io.vec(3);
    ellipsoid.axis[2] = io.vec(3);
    ellipsoid.extent = io.vec(3);
    const E = io.vec(3);
    ellipsoid.center = io.vec(3);
    const n = io.real();
    const constant = io.real();
    return { ellipsoid, E, N, U, V, constant, n };
}

function outEllipse2(io: OracleIO, e: Hyperellipsoid): void {
    io.outVec(e.center);
    io.outVec(e.axis[0]);
    io.outVec(e.axis[1]);
    io.outVec(e.extent);
}

describe('oracle: v02-algebra', () => {
    const family = new OracleFamily('v02-algebra');

    // ---------------- Matrix2x2.h ----------------

    family.case('Matrix2x2.inverse', (io) => {
        const M = mat(io, 2);
        const r = inverse2x2(M);
        io.outBool(r.invertible);
        io.outMat(r.inverse);
        io.outMat(adjoint2x2(M));
        io.outReal(determinant2x2(M));
        io.outReal(trace2x2(M));
        io.outMat(inverse2x2(M).inverse);
    }, { exact: true });

    // MakeRotation calls cos and sin, GetRotationAngle calls atan2.
    family.case('Matrix2x2.rotation', (io) => {
        const angle = io.real();
        const rotation = new Matrix(2, 2);
        makeRotation2x2(angle, rotation);
        io.outMat(rotation);
        io.outReal(getRotationAngle2x2(rotation));
        const other = mat(io, 2);
        io.outReal(getRotationAngle2x2(other));
    });

    family.case('Matrix2x2.doTransform', (io) => {
        const A = mat(io, 2);
        const B = mat(io, 2);
        const V = io.vec(2);
        const W = io.vec(2);
        const i = io.integer();
        io.outVec(doTransform2x2(A, V));
        io.outMat(doTransform2x2(A, B));
        io.outVec(getBasis2x2(A, i));
        const C = A.clone();
        setBasis2x2(C, i, W);
        io.outMat(C);
    }, { exact: true });

    // ---------------- Matrix3x3.h ----------------

    family.case('Matrix3x3.inverse', (io) => {
        const M = mat(io, 3);
        const r = inverse3x3(M);
        io.outBool(r.invertible);
        io.outMat(r.inverse);
        io.outMat(adjoint3x3(M));
        io.outReal(determinant3x3(M));
        io.outReal(trace3x3(M));
        io.outMat(inverse3x3(M).inverse);
    }, { exact: true });

    family.case('Matrix3x3.doTransform', (io) => {
        const A = mat(io, 3);
        const B = mat(io, 3);
        const V = io.vec(3);
        const W = io.vec(3);
        const i = io.integer();
        io.outVec(doTransform3x3(A, V));
        io.outMat(doTransform3x3(A, B));
        io.outVec(getBasis3x3(A, i));
        const C = A.clone();
        setBasis3x3(C, i, W);
        io.outMat(C);
    }, { exact: true });

    // ---------------- Matrix4x4.h ----------------

    family.case('Matrix4x4.inverse', (io) => {
        const M = mat(io, 4);
        const r = inverse4x4(M);
        io.outBool(r.invertible);
        io.outMat(r.inverse);
        io.outMat(adjoint4x4(M));
        io.outReal(determinant4x4(M));
        io.outReal(trace4x4(M));
        io.outMat(inverse4x4(M).inverse);
    }, { exact: true });

    family.case('Matrix4x4.doTransform', (io) => {
        const A = mat(io, 4);
        const B = mat(io, 4);
        const V = io.vec(4);
        const W = io.vec(4);
        const i = io.integer();
        io.outVec(doTransform4x4(A, V));
        io.outMat(doTransform4x4(A, B));
        io.outVec(getBasis4x4(A, i));
        const C = A.clone();
        setBasis4x4(C, i, W);
        io.outMat(C);
    }, { exact: true });

    family.case('Matrix4x4.makeObliqueProjection', (io) => {
        const origin = io.vec(4);
        const normal = io.vec(4);
        const direction = io.vec(4);
        io.outMat(makeObliqueProjection4x4(origin, normal, direction));
    }, { exact: true });

    family.case('Matrix4x4.makePerspectiveProjection', (io) => {
        const origin = io.vec(4);
        const normal = io.vec(4);
        const eye = io.vec(4);
        io.outMat(makePerspectiveProjection4x4(origin, normal, eye));
    }, { exact: true });

    family.case('Matrix4x4.makeReflection', (io) => {
        const origin = io.vec(4);
        const normal = io.vec(4);
        io.outMat(makeReflection4x4(origin, normal));
    }, { exact: true });

    // ---------------- Quaternion.h ----------------

    family.case('Quaternion.arithmetic', (io) => {
        const q0 = quat(io);
        const q1 = quat(io);
        const s = io.real();
        io.outVec(negateQuaternion(q0));
        // Upstream's unary operator+ returns its argument unchanged; the port
        // has no separate function for it.
        io.outVec(q0);
        io.outVec(addQuaternion(q0, q1));
        io.outVec(subQuaternion(q0, q1));
        io.outVec(mulQuaternion(q0, s));
        io.outVec(mulQuaternion(s, q0));
        io.outVec(divQuaternion(q0, s));
        io.outReal(dot(q0, q1));
        io.outReal(length(q0));
        const unit = q0.clone();
        io.outReal(normalize(unit));
        io.outVec(unit);
        io.outBool(q0.equals(q1));
        io.outBool(q0.notEquals(q1));
        io.outBool(q0.lessThan(q1));
        io.outBool(q0.lessThanOrEqual(q1));
        io.outBool(q0.greaterThan(q1));
        io.outBool(q0.greaterThanOrEqual(q1));
        io.outVec(Quaternion.zero());
        io.outVec(Quaternion.i());
        io.outVec(Quaternion.j());
        io.outVec(Quaternion.k());
        io.outVec(Quaternion.identity());
    }, { exact: true });

    family.case('Quaternion.multiply', (io) => {
        const q0 = quat(io);
        const q1 = quat(io);
        io.outVec(mulQuaternion(q0, q1));
        io.outVec(mulQuaternion(q1, q0));
        io.outVec(conjugate(q0));
        io.outVec(inverseQuaternion(q0));
        io.outVec(inverseQuaternion(Quaternion.zero()));
    }, { exact: true });

    family.case('Quaternion.rotate.3d', (io) => {
        const q = quat(io);
        const u = io.vec(3);
        io.outVec(rotate(q, u));
    }, { exact: true });

    family.case('Quaternion.rotate.4d', (io) => {
        const q = quat(io);
        const u = io.vec(4);
        io.outVec(rotate(q, u));
    }, { exact: true });

    // ChebyshevRatiosUsingCosAngle calls acos and sin unless the cosine is
    // at least 1, which is the exact arm.
    family.case('Quaternion.slerp', (io) => {
        const q0 = quat(io);
        const q1 = quat(io);
        const t = io.real();
        io.outVec(slerpQuaternion(t, q0, q1));
    });

    // Upstream's SlerpR is Slerp.ts's slerpUsingCosAngle with the dot product
    // of the (preprocessed) quaternions as the cosine; see the port note in
    // src/Quaternion.ts.
    family.case('Quaternion.slerpR', (io) => {
        const q0 = quat(io);
        const q1 = quat(io);
        const t = io.real();
        io.outReals(slerpUsingCosAngle(t, q0.values, q1.values, dot(q0, q1)));
    });

    family.case('Quaternion.slerpRP', (io) => {
        const q0 = quat(io);
        const q1 = quat(io);
        const cosA = io.real();
        const t = io.real();
        io.outReals(slerpUsingCosAngle(t, q0.values, q1.values, cosA));
    });

    family.case('Quaternion.slerpRPH', (io) => {
        const q0 = quat(io);
        const q1 = quat(io);
        const qh = quat(io);
        const cosAH = io.real();
        const t = io.real();
        io.outReals(slerpUsingMidpoint(t, q0.values, q1.values, qh.values,
            cosAH));
    });

    // Throw parity: chebyshevRatiosUsingCosAngle calls logError for an angle
    // of pi (cosine -1).
    family.case('Quaternion.slerpRP.invalidAngle', (io) => {
        const q0 = quat(io);
        const q1 = quat(io);
        const cosA = io.real();
        const t = io.real();
        io.outReals(slerpUsingCosAngle(t, q0.values, q1.values, cosA));
    });

    // ---------------- Rotation.h ----------------

    family.case('Rotation.quaternionToMatrix', (io) => {
        const q = quat(io);
        io.outMat(Rotation.fromQuaternion(q, 3).toMatrix());
        io.outMat(Rotation.fromQuaternion(q, 4).toMatrix());
    }, { exact: true });

    family.case('Rotation.matrixToQuaternion', (io) => {
        const m3 = mat(io, 3);
        io.outVec(Rotation.fromMatrix(m3).toQuaternion());
        io.outVec(Rotation.fromMatrix(lift3To4(m3)).toQuaternion());
    }, { exact: true });

    // Calls cos and sin.
    family.case('Rotation.axisAngleToMatrix', (io) => {
        const axis = io.vec(3);
        const angle = io.real();
        io.outMat(Rotation.fromAxisAngle(new AxisAngle(axis, angle)).toMatrix());
        const axis4 = Vector.fromArray([axis.values[0], axis.values[1],
            axis.values[2], 0]);
        io.outMat(Rotation.fromAxisAngle(new AxisAngle(axis4, angle)).toMatrix());
    });

    // The axis is arithmetic only, the angle comes from acos.
    family.case('Rotation.matrixToAxisAngle', (io) => {
        const m3 = mat(io, 3);
        outAxisAngle(io, Rotation.fromMatrix(m3).toAxisAngle());
        outAxisAngle(io, Rotation.fromMatrix(lift3To4(m3)).toAxisAngle());
    });

    // Calls sin and cos.
    family.case('Rotation.axisAngleToQuaternion', (io) => {
        const axis = io.vec(3);
        const angle = io.real();
        io.outVec(Rotation.fromAxisAngle(new AxisAngle(axis, angle))
            .toQuaternion());
        const axis4 = Vector.fromArray([axis.values[0], axis.values[1],
            axis.values[2], 0]);
        io.outVec(Rotation.fromAxisAngle(new AxisAngle(axis4, angle))
            .toQuaternion());
    });

    // The axis is arithmetic only, the angle comes from acos.
    family.case('Rotation.quaternionToAxisAngle', (io) => {
        const q = quat(io);
        outAxisAngle(io, Rotation.fromQuaternion(q, 3).toAxisAngle());
        outAxisAngle(io, Rotation.fromQuaternion(q, 4).toAxisAngle());
    });

    // Calls atan2, asin and acos. The input matrix is recorded, so the
    // gimbal-lock tests see identical doubles on both sides.
    family.case('Rotation.matrixToEulerAngles', (io) => {
        const which = io.integer();
        const [i0, i1, i2] = EULER_AXES[which];
        const m3 = mat(io, 3);
        outEulerAngles(io, Rotation.fromMatrix(m3).toEulerAngles(i0, i1, i2));
        outEulerAngles(io,
            Rotation.fromMatrix(lift3To4(m3)).toEulerAngles(i0, i1, i2));
    });

    // Calls cos and sin.
    family.case('Rotation.eulerAnglesToMatrix', (io) => {
        const which = io.integer();
        const [i0, i1, i2] = EULER_AXES[which];
        const e = new EulerAngles(i0, i1, i2, io.real(), io.real(), io.real());
        io.outMat(Rotation.fromEulerAngles(e, 3).toMatrix());
        io.outMat(Rotation.fromEulerAngles(e, 4).toMatrix());
    });

    family.case('Rotation.quaternionToEulerAngles', (io) => {
        const which = io.integer();
        const [i0, i1, i2] = EULER_AXES[which];
        const q = quat(io);
        outEulerAngles(io,
            Rotation.fromQuaternion(q, 3).toEulerAngles(i0, i1, i2));
        outEulerAngles(io,
            Rotation.fromQuaternion(q, 4).toEulerAngles(i0, i1, i2));
    });

    family.case('Rotation.eulerAnglesToQuaternion', (io) => {
        const which = io.integer();
        const [i0, i1, i2] = EULER_AXES[which];
        const e = new EulerAngles(i0, i1, i2, io.real(), io.real(), io.real());
        io.outVec(Rotation.fromEulerAngles(e, 3).toQuaternion());
        io.outVec(Rotation.fromEulerAngles(e, 4).toQuaternion());
    });

    family.case('Rotation.axisAngleToEulerAngles', (io) => {
        const which = io.integer();
        const [i0, i1, i2] = EULER_AXES[which];
        const axis = io.vec(3);
        const angle = io.real();
        outEulerAngles(io, Rotation.fromAxisAngle(new AxisAngle(axis, angle))
            .toEulerAngles(i0, i1, i2));
        const axis4 = Vector.fromArray([axis.values[0], axis.values[1],
            axis.values[2], 0]);
        outEulerAngles(io, Rotation.fromAxisAngle(new AxisAngle(axis4, angle))
            .toEulerAngles(i0, i1, i2));
    });

    family.case('Rotation.eulerAnglesToAxisAngle', (io) => {
        const which = io.integer();
        const [i0, i1, i2] = EULER_AXES[which];
        const e = new EulerAngles(i0, i1, i2, io.real(), io.real(), io.real());
        outAxisAngleTol(io, Rotation.fromEulerAngles(e, 3).toAxisAngle());
        outAxisAngleTol(io, Rotation.fromEulerAngles(e, 4).toAxisAngle());
    });

    family.case('Rotation.passthrough', (io) => {
        const which = io.integer();
        const [i0, i1, i2] = EULER_AXES[which];
        const m3 = mat(io, 3);
        const q = quat(io);
        const axis = io.vec(3);
        const angle = io.real();
        const e = new EulerAngles(i0, i1, i2, io.real(), io.real(), io.real());
        io.outMat(Rotation.fromMatrix(m3).toMatrix());
        io.outVec(Rotation.fromQuaternion(q, 3).toQuaternion());
        outAxisAngle(io,
            Rotation.fromAxisAngle(new AxisAngle(axis, angle)).toAxisAngle());
        outEulerAngles(io,
            Rotation.fromEulerAngles(e, 3).toEulerAngles(i0, i1, i2));
    }, { exact: true });

    // Deliberate port deviation: issue #374, the axis extracted from
    // R - Transpose(R) near angle pi.
    family.case('Rotation.matrixToAxisAngle.nearPi', (io) => {
        const m3 = mat(io, 3);
        outAxisAngle(io, Rotation.fromMatrix(m3).toAxisAngle());
        outAxisAngle(io, Rotation.fromMatrix(lift3To4(m3)).toAxisAngle());
    }, { deviation: 'issue #374 (Rotation.h Convert(Matrix, AxisAngle))' });

    // Deliberate port deviation: issue #225, operator()(i0,i1,i2) on an
    // Euler-sourced Rotation asked for a different factorization.
    family.case('Rotation.eulerAngles.refactorization', (io) => {
        const which0 = io.integer();
        const step = io.integer();
        const [i0, i1, i2] = EULER_AXES[which0];
        const [j0, j1, j2] = EULER_AXES[(which0 + step) % 12];
        const e = new EulerAngles(i0, i1, i2, io.real(), io.real(), io.real());
        outEulerAngles(io,
            Rotation.fromEulerAngles(e, 3).toEulerAngles(j0, j1, j2));
    }, { deviation: 'issue #225 (Rotation.h operator()(i0,i1,i2))' });

    // ---------------- RotationEstimate.h ----------------

    family.case('RotationEstimate.rotCEstimate', (io) => {
        const degree = 4 + 2 * io.integer();
        const t = io.real();
        io.outReal(rotC0Estimate(t, degree));
        io.outReal(rotC1Estimate(t, degree));
        io.outReal(rotC2Estimate(t, degree));
        io.outReal(rotC3Estimate(t, degree));
        io.outReal(rotC4Estimate(t, degree));
        io.outReal(getRotC0EstimateMaxError(degree));
        io.outReal(getRotC1EstimateMaxError(degree));
        io.outReal(getRotC2EstimateMaxError(degree));
        io.outReal(getRotC3EstimateMaxError(degree));
        io.outReal(getRotC4EstimateMaxError(degree));
    }, { exact: true });

    family.case('RotationEstimate.rotationEstimate', (io) => {
        const degree = 4 + 2 * io.integer();
        const p = io.vec(3);
        io.outMat(rotationEstimate(p, degree));
        for (const m of rotationDerivativeEstimate(p, degree)) { io.outMat(m); }
        const both = rotationAndDerivativeEstimate(p, degree);
        io.outMat(both.R);
        for (const m of both.Rder) { io.outMat(m); }
    }, { exact: true });

    // ---------------- Transform.h ----------------

    family.case('Transform.rotationMatrix', (io) => {
        const rot3 = mat(io, 3);
        const tr = io.vec(3);
        const s0 = io.real();
        const s1 = io.real();
        const s2 = io.real();
        const us = io.real();

        const t = new Transform();
        outTransform(io, t);
        t.setRotation(io.index % 2 === 0 ? lift3To4(rot3) : rot3);
        applyTranslation(io, t, tr);
        applyScale(io, t, s0, s1, s2, us);
        outTransformFull(io, t);
        t.makeIdentity();
        outTransform(io, t);
    }, { exact: true });

    family.case('Transform.rotationQuaternion', (io) => {
        const q = quat(io);
        const tr = io.vec(3);
        const us = io.real();
        const t = new Transform();
        t.setRotation(q);
        applyTranslation(io, t, tr);
        t.setUniformScale(us);
        outTransformFull(io, t);
    }, { exact: true });

    family.case('Transform.setMatrix', (io) => {
        const m = mat(io, 4);
        const tr = io.vec(3);
        const t = new Transform();
        t.setMatrix(m);
        applyTranslation(io, t, tr);
        outTransformFull(io, t);
    }, { exact: true });

    // Calls sin, cos and acos.
    family.case('Transform.rotationAxisAngle', (io) => {
        const axis = io.vec(3);
        const magnitude = io.real();
        const sign = io.real();
        const angle = magnitude * sign;
        const tr = io.vec(3);

        const t = new Transform();
        if (io.index % 2 === 0) {
            t.setRotation(new AxisAngle(axis, angle));
        } else {
            const axis4 = Vector.fromArray([axis.values[0], axis.values[1],
                axis.values[2], 0]);
            t.setRotation(new AxisAngle(axis4, angle));
        }
        applyTranslation(io, t, tr);
        outTransform(io, t);
        outAxisAngleTol(io, t.getRotationAxisAngle3());
        outAxisAngleTol(io, t.getRotationAxisAngle4());
    });

    // Calls sin, cos, atan2, asin and acos.
    family.case('Transform.rotationEulerAngles', (io) => {
        const which = io.integer();
        const [i0, i1, i2] = EULER_AXES[which];
        const e = new EulerAngles(i0, i1, i2, io.real(), io.real(), io.real());
        const tr = io.vec(3);
        const t = new Transform();
        t.setRotation(e);
        applyTranslation(io, t, tr);
        outTransform(io, t);
        outEulerAngles(io, t.getRotationEulerAngles(i0, i1, i2));
    });

    family.case('Transform.inverse', (io) => {
        const kind = io.integer();
        const rot3 = mat(io, 3);
        const general = mat(io, 4);
        const tr = io.vec(3);
        const s0 = io.real();
        const s1 = io.real();
        const s2 = io.real();
        const us = io.real();

        const t = makeTransform(kind, rot3, general, tr, s0, s1, s2, us);
        io.outMat(t.getHInverse());
        const inv = t.inverse();
        outTransform(io, inv);
        io.outVec(inv.getTranslation());
        if (inv.isRSMatrix()) {
            io.outMat(inv.getRotation());
            io.outVec(inv.getScale());
        }
        if (inv.isUniformScale()) {
            io.outReal(inv.getUniformScale());
        }
    }, { exact: true });

    // Deliberate port deviation: issue #265, Transform::Inverse stores the
    // translation in the M channel.
    family.case('Transform.inverse.matrixChannel', (io) => {
        const general = mat(io, 4);
        const tr = io.vec(3);
        const t = new Transform();
        t.setMatrix(general);
        t.setTranslation(tr.values[0], tr.values[1], tr.values[2]);
        io.outMat(t.inverse().getMatrix());
    }, { deviation: 'issue #265 (Transform.h Inverse M channel)' });

    // Deliberate port deviation: issue #265, GetHInverse leaves the last row
    // of a previously cached singular general inverse in place.
    family.case('Transform.getHInverse.lastRow', (io) => {
        const singular = mat(io, 4);
        const rot3 = mat(io, 3);
        const tr = io.vec(3);
        const t = new Transform();
        t.setMatrix(singular);
        io.outMat(t.getHInverse());
        t.setRotation(rot3);
        t.setTranslation(tr.values[0], tr.values[1], tr.values[2]);
        io.outMat(t.getHInverse());
    }, { deviation: 'issue #265 (Transform.h GetHInverse last row)' });

    family.case('Transform.multiply', (io) => {
        const kindA = io.integer();
        const kindB = io.integer();
        const rotA = mat(io, 3);
        const rotB = mat(io, 3);
        const generalA = mat(io, 4);
        const generalB = mat(io, 4);
        const trA = io.vec(3);
        const trB = io.vec(3);
        const s0 = io.real();
        const s1 = io.real();
        const s2 = io.real();
        const us = io.real();
        const V = io.vec(4);
        const K = mat(io, 4);

        const A = makeTransform(kindA, rotA, generalA, trA, s0, s1, s2, us);
        const B = makeTransform(kindB, rotB, generalB, trB, s2, s0, s1, us);
        // Formed first: the product asserts when the two scales multiply to
        // zero, and a throwing record must not have emitted anything.
        const product = mulTransform(A, B);
        io.outVec(mulTransform(A, V));
        io.outVec(mulTransform(V, A));
        outTransform(io, product);
        io.outVec(product.getTranslation());
        io.outMat(mulTransform(K, A));
        io.outMat(mulTransform(A, K));
    }, { exact: true });

    // Throw parity for the logAssert preconditions.
    family.case('Transform.assertions', (io) => {
        const which = io.index % 4;
        const rot3 = mat(io, 3);
        const general = mat(io, 4);
        const s = io.real();

        const t = new Transform();
        t.setRotation(rot3);
        if (which === 0) {
            t.setScale(0, s, s);
        } else if (which === 1) {
            t.setMatrix(general);
            t.setScale(s, s, s);
        } else if (which === 2) {
            t.setScale(s, s + 1, s + 2);
            io.outReal(t.getUniformScale());
        } else {
            t.setUniformScale(s);
        }
        outTransform(io, t);
    }, { exact: true });

    family.case('Transform.identity', (io) => {
        outTransformFull(io, new Transform());
        const identity = Transform.identity();
        outTransformFull(io, identity);
        io.outMat(identity.getHMatrix());
    }, { exact: true });

    // ---------------- ConvertCoordinates.h ----------------

    family.case('ConvertCoordinates.convert', (io) => {
        const n = io.integer();
        const vU = io.boolean();
        const vV = io.boolean();
        const U = mat(io, n);
        const V = mat(io, n);
        const X = io.vec(n);
        const Y = io.vec(n);
        const A = mat(io, n);
        const B = mat(io, n);

        const convert = new ConvertCoordinates(n);
        io.outBool(convert.compute(U, vU, V, vV));
        io.outMat(convert.getC());
        io.outMat(convert.getInverseC());
        io.outBool(convert.isVectorOnRightU());
        io.outBool(convert.isVectorOnRightV());
        io.outBool(convert.isRightHandedU());
        io.outBool(convert.isRightHandedV());
        io.outVec(convert.uToV(X));
        io.outVec(convert.vToU(Y));
        io.outMat(convert.uToV(A));
        io.outMat(convert.vToU(B));
    }, { exact: true });

    // ---------------- Projection.h ----------------

    family.case('Projection.projectEllipse2', (io) => {
        const ellipse = new Hyperellipsoid(2);
        ellipse.center = io.vec(2);
        ellipse.axis[0] = io.vec(2);
        ellipse.axis[1] = io.vec(2);
        ellipse.extent = io.vec(2);
        const origin = io.vec(2);
        const direction = io.vec(2);
        const r = projectEllipse2(ellipse, Line.fromOriginDirection(origin,
            direction));
        io.outReal(r.smin);
        io.outReal(r.smax);
    }, { exact: true });

    family.case('Projection.projectEllipsoid3', (io) => {
        const ellipsoid = new Hyperellipsoid(3);
        ellipsoid.center = io.vec(3);
        ellipsoid.axis[0] = io.vec(3);
        ellipsoid.axis[1] = io.vec(3);
        ellipsoid.axis[2] = io.vec(3);
        ellipsoid.extent = io.vec(3);
        const origin = io.vec(3);
        const direction = io.vec(3);
        const r = projectEllipsoid3(ellipsoid, Line.fromOriginDirection(origin,
            direction));
        io.outReal(r.smin);
        io.outReal(r.smax);
    }, { exact: true });

    // Ellipse2::FromCoefficients runs a SymmetricEigensolver, but that solver
    // uses only sqrt and fabs, so the whole path is arithmetic and the
    // comparison is exact. The 2x2 Inverse it takes is the closed form of
    // Matrix2x2.h, which this translation unit makes visible; the port agrees
    // since PR #508.
    family.case('Projection.perspectiveProject.basis', (io) => {
        const s = projectionSetup(io);
        const r = perspectiveProject(s.ellipsoid, s.E, s.N, s.U, s.V, s.n);
        if (!r.isEllipse) {
            throw new Error('fromCoefficientsABC returned false');
        }
        outEllipse2(io, r.ellipse);
    }, { exact: true });

    family.case('Projection.perspectiveProject.plane', (io) => {
        const s = projectionSetup(io);
        const plane = Hyperplane.fromNormalConstant(s.N, s.constant);
        const r = perspectiveProject(s.ellipsoid, s.E, plane);
        if (!r.isEllipse) {
            throw new Error('fromCoefficientsABC returned false');
        }
        outEllipse2(io, r.ellipse);
    }, { exact: true });

    family.finish();
});
