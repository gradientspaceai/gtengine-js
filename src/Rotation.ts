// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Rotation.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Conversions among various representations of rotations. The value of N must
// be 3 or 4. The latter case supports affine algebra when you use 4-tuple
// vectors (w-component is 1 for points and 0 for vectors) and 4x4 matrices
// for affine transformations. Rotation axes must be unit length. The angles
// are in radians. The Euler angles are in world coordinates; upstream has not
// yet added support for body coordinates.
//
// Port notes:
// - Upstream's compile-time dimension 'Rotation<N, Real>' becomes a runtime
//   dimension. The matrix and axis-angle factories infer N from their
//   argument; the quaternion and Euler-angle factories take N explicitly
//   (defaulting to 3), since neither representation carries a dimension.
// - The C++ conversion operators become methods: 'operator Matrix<N,N,Real>'
//   is toMatrix(), 'operator Quaternion<Real>' is toQuaternion(),
//   'operator AxisAngle<N,Real>' is toAxisAngle(), and
//   'operator()(i0,i1,i2)' is toEulerAngles(i0, i1, i2).
// - Upstream caches the converted representation in a 'mutable' member and
//   returns it. The port caches likewise, but every accessor returns a fresh
//   copy so that callers cannot mutate the cache through the returned object
//   (C++ returns by value for the matrix/quaternion/axis-angle operators).
// - Only the GTE_USE_MAT_VEC branches are ported, matching the rest of the
//   library (see Matrix3x3.ts, Quaternion.ts).
// - toEulerAngles fixes an upstream defect for Euler-angle-sourced rotations;
//   see the comment in that method.

import { AxisAngle } from './AxisAngle.js';
import { GTE_C_HALF_PI, GTE_C_PI } from './Constants.js';
import { EulerAngles, EulerResult } from './EulerAngles.js';
import { logAssert } from './Logger.js';
import { Matrix, multiplyAB } from './Matrix.js';
import { Quaternion } from './Quaternion.js';
import { Vector, normalize } from './Vector.js';

function assertDimension(n: number): void {
    logAssert(n === 3 || n === 4, 'Dimension must be 3 or 4.');
}

// Convert a rotation matrix to a quaternion.
//
// x^2 = (+r00 - r11 - r22 + 1)/4
// y^2 = (-r00 + r11 - r22 + 1)/4
// z^2 = (-r00 - r11 + r22 + 1)/4
// w^2 = (+r00 + r11 + r22 + 1)/4
// x^2 + y^2 = (1 - r22)/2
// z^2 + w^2 = (1 + r22)/2
// y^2 - x^2 = (r11 - r00)/2
// w^2 - z^2 = (r11 + r00)/2
// x*y = (r01 + r10)/4
// x*z = (r02 + r20)/4
// y*z = (r12 + r21)/4
// x*w = (r21 - r12)/4
// y*w = (r02 - r20)/4
// z*w = (r10 - r01)/4
//
// If Q is the 4x1 column vector (x,y,z,w), the previous equations give us
//         +-                  -+
//         | x*x  x*y  x*z  x*w |
// Q*Q^T = | y*x  y*y  y*z  y*w |
//         | z*x  z*y  z*z  z*w |
//         | w*x  w*y  w*z  w*w |
//         +-                  -+
// The code extracts the row of maximum length, normalizing it to obtain the
// result q.
function convertMatrixToQuaternion(r: Matrix): Quaternion {
    const q = new Quaternion();

    const r22 = r.get(2, 2);
    if (r22 <= 0) {  // x^2 + y^2 >= z^2 + w^2
        const dif10 = r.get(1, 1) - r.get(0, 0);
        const omr22 = 1 - r22;
        if (dif10 <= 0) {  // x^2 >= y^2
            const fourXSqr = omr22 - dif10;
            const inv4x = 0.5 / Math.sqrt(fourXSqr);
            q.values[0] = fourXSqr * inv4x;
            q.values[1] = (r.get(0, 1) + r.get(1, 0)) * inv4x;
            q.values[2] = (r.get(0, 2) + r.get(2, 0)) * inv4x;
            q.values[3] = (r.get(2, 1) - r.get(1, 2)) * inv4x;
        } else {  // y^2 >= x^2
            const fourYSqr = omr22 + dif10;
            const inv4y = 0.5 / Math.sqrt(fourYSqr);
            q.values[0] = (r.get(0, 1) + r.get(1, 0)) * inv4y;
            q.values[1] = fourYSqr * inv4y;
            q.values[2] = (r.get(1, 2) + r.get(2, 1)) * inv4y;
            q.values[3] = (r.get(0, 2) - r.get(2, 0)) * inv4y;
        }
    } else {  // z^2 + w^2 >= x^2 + y^2
        const sum10 = r.get(1, 1) + r.get(0, 0);
        const opr22 = 1 + r22;
        if (sum10 <= 0) {  // z^2 >= w^2
            const fourZSqr = opr22 - sum10;
            const inv4z = 0.5 / Math.sqrt(fourZSqr);
            q.values[0] = (r.get(0, 2) + r.get(2, 0)) * inv4z;
            q.values[1] = (r.get(1, 2) + r.get(2, 1)) * inv4z;
            q.values[2] = fourZSqr * inv4z;
            q.values[3] = (r.get(1, 0) - r.get(0, 1)) * inv4z;
        } else {  // w^2 >= z^2
            const fourWSqr = opr22 + sum10;
            const inv4w = 0.5 / Math.sqrt(fourWSqr);
            q.values[0] = (r.get(2, 1) - r.get(1, 2)) * inv4w;
            q.values[1] = (r.get(0, 2) - r.get(2, 0)) * inv4w;
            q.values[2] = (r.get(1, 0) - r.get(0, 1)) * inv4w;
            q.values[3] = fourWSqr * inv4w;
        }
    }

    return q;
}

// Convert a quaternion q = x*i + y*j + z*k + w to a rotation matrix.
//     +-           -+   +-                                     -+
// R = | r00 r01 r02 | = | 1-2y^2-2z^2  2(xy-zw)     2(xz+yw)    |
//     | r10 r11 r12 |   | 2(xy+zw)     1-2x^2-2z^2  2(yz-xw)    |
//     | r20 r21 r22 |   | 2(xz-yw)     2(yz+xw)     1-2x^2-2y^2 |
//     +-           -+   +-                                     -+
function convertQuaternionToMatrix(q: Quaternion, n: number): Matrix {
    assertDimension(n);
    const r = Matrix.identity(n, n);

    const twoX = 2 * q.values[0];
    const twoY = 2 * q.values[1];
    const twoZ = 2 * q.values[2];
    const twoXX = twoX * q.values[0];
    const twoXY = twoX * q.values[1];
    const twoXZ = twoX * q.values[2];
    const twoXW = twoX * q.values[3];
    const twoYY = twoY * q.values[1];
    const twoYZ = twoY * q.values[2];
    const twoYW = twoY * q.values[3];
    const twoZZ = twoZ * q.values[2];
    const twoZW = twoZ * q.values[3];

    r.set(0, 0, 1 - twoYY - twoZZ);
    r.set(0, 1, twoXY - twoZW);
    r.set(0, 2, twoXZ + twoYW);
    r.set(1, 0, twoXY + twoZW);
    r.set(1, 1, 1 - twoXX - twoZZ);
    r.set(1, 2, twoYZ - twoXW);
    r.set(2, 0, twoXZ - twoYW);
    r.set(2, 1, twoYZ + twoXW);
    r.set(2, 2, 1 - twoXX - twoYY);

    return r;
}

// Convert a rotation matrix to an axis-angle pair. Let (x0,x1,x2) be the axis
// and let t be an angle of rotation. The rotation matrix is
//   R = I + sin(t)*S + (1-cos(t))*S^2
// where I is the identity and S = {{0,-x2,x1},{x2,0,-x0},{-x1,x0,0}} where
// the inner-brace triples are the rows of the matrix. If t > 0, R represents
// a counterclockwise rotation. It may be shown that cos(t) = (trace(R)-1)/2
// and R - Transpose(R) = 2*sin(t)*S. As long as sin(t) is not zero, we may
// solve for S in the second equation, which produces the axis direction
// U = (S21,S02,S10). When t = 0, the rotation is the identity, in which case
// any axis direction is valid; we choose (1,0,0). When t = pi, it must be
// that R - Transpose(R) = 0, which prevents us from extracting the axis.
// Instead, note that (R+I)/2 = I+S^2 = U*U^T, where U is a unit-length axis
// direction.
function convertMatrixToAxisAngle(r: Matrix, n: number): AxisAngle {
    assertDimension(n);

    const trace = r.get(0, 0) + r.get(1, 1) + r.get(2, 2);
    const half = 0.5;
    let cs = half * (trace - 1);
    cs = Math.max(Math.min(cs, 1), -1);
    // The angle is in [0,pi]. The axis starts out zero (upstream calls
    // a.axis.MakeZero()).
    const a = new AxisAngle(new Vector(n), Math.acos(cs));

    if (a.angle > 0) {
        if (a.angle < GTE_C_PI) {
            // The angle is in (0,pi).
            a.axis.values[0] = r.get(2, 1) - r.get(1, 2);
            a.axis.values[1] = r.get(0, 2) - r.get(2, 0);
            a.axis.values[2] = r.get(1, 0) - r.get(0, 1);
            normalize(a.axis);
        } else {
            // The angle is pi, in which case R is symmetric and
            // R+I = 2*(I+S^2) = 2*U*U^T, where U = (u0,u1,u2) is the
            // unit-length direction of the rotation axis. Determine the
            // largest diagonal entry of R+I and normalize the corresponding
            // row to produce U. It does not matter the sign on u[d] for the
            // chosen diagonal d, because R(U,pi) = R(-U,pi).
            const one = 1;
            if (r.get(0, 0) >= r.get(1, 1)) {
                if (r.get(0, 0) >= r.get(2, 2)) {
                    // r00 is maximum diagonal term
                    a.axis.values[0] = r.get(0, 0) + one;
                    a.axis.values[1] = half * (r.get(0, 1) + r.get(1, 0));
                    a.axis.values[2] = half * (r.get(0, 2) + r.get(2, 0));
                } else {
                    // r22 is maximum diagonal term
                    a.axis.values[0] = half * (r.get(2, 0) + r.get(0, 2));
                    a.axis.values[1] = half * (r.get(2, 1) + r.get(1, 2));
                    a.axis.values[2] = r.get(2, 2) + one;
                }
            } else {
                if (r.get(1, 1) >= r.get(2, 2)) {
                    // r11 is maximum diagonal term
                    a.axis.values[0] = half * (r.get(1, 0) + r.get(0, 1));
                    a.axis.values[1] = r.get(1, 1) + one;
                    a.axis.values[2] = half * (r.get(1, 2) + r.get(2, 1));
                } else {
                    // r22 is maximum diagonal term
                    a.axis.values[0] = half * (r.get(2, 0) + r.get(0, 2));
                    a.axis.values[1] = half * (r.get(2, 1) + r.get(1, 2));
                    a.axis.values[2] = r.get(2, 2) + one;
                }
            }
            normalize(a.axis);
        }
    } else {
        // The angle is 0 and the matrix is the identity. Any axis will work,
        // so choose the Unit(0) axis.
        a.axis.values[0] = 1;
    }

    return a;
}

// Convert an axis-angle pair to a rotation matrix. Assuming (x0,x1,x2) is for
// a right-handed world (x0 to right, x1 up, x2 out of plane of page), a
// positive angle corresponds to a counterclockwise rotation from the
// perspective of an observer looking at the origin of the plane of rotation
// and having view direction the negative of the axis direction. The
// coordinate-axis rotations are, where unit(0) = (1,0,0), unit(1) = (0,1,0),
// unit(2) = (0,0,1),
//   R(unit(0),t) = {{ 1, 0, 0}, { 0, c,-s}, { 0, s, c}}
//   R(unit(1),t) = {{ c, 0, s}, { 0, 1, 0}, {-s, 0, c}}
//   R(unit(2),t) = {{ c,-s, 0}, { s, c, 0}, { 0, 0, 1}}
// where c = cos(t), s = sin(t), and the inner-brace triples are rows of the
// matrix. The general matrix is
//      +-                                                          -+
//  R = | (1-c)*x0^2  + c     (1-c)*x0*x1 - s*x2  (1-c)*x0*x2 + s*x1 |
//      | (1-c)*x0*x1 + s*x2  (1-c)*x1^2  + c     (1-c)*x1*x2 - s*x0 |
//      | (1-c)*x0*x2 - s*x1  (1-c)*x1*x2 + s*x0  (1-c)*x2^2  + c    |
//      +-                                                          -+
function convertAxisAngleToMatrix(a: AxisAngle, n: number): Matrix {
    assertDimension(n);
    const r = Matrix.identity(n, n);

    const cs = Math.cos(a.angle);
    const sn = Math.sin(a.angle);
    const oneMinusCos = 1 - cs;
    const x0sqr = a.axis.values[0] * a.axis.values[0];
    const x1sqr = a.axis.values[1] * a.axis.values[1];
    const x2sqr = a.axis.values[2] * a.axis.values[2];
    const x0x1m = a.axis.values[0] * a.axis.values[1] * oneMinusCos;
    const x0x2m = a.axis.values[0] * a.axis.values[2] * oneMinusCos;
    const x1x2m = a.axis.values[1] * a.axis.values[2] * oneMinusCos;
    const x0Sin = a.axis.values[0] * sn;
    const x1Sin = a.axis.values[1] * sn;
    const x2Sin = a.axis.values[2] * sn;

    r.set(0, 0, x0sqr * oneMinusCos + cs);
    r.set(0, 1, x0x1m - x2Sin);
    r.set(0, 2, x0x2m + x1Sin);
    r.set(1, 0, x0x1m + x2Sin);
    r.set(1, 1, x1sqr * oneMinusCos + cs);
    r.set(1, 2, x1x2m - x0Sin);
    r.set(2, 0, x0x2m - x1Sin);
    r.set(2, 1, x1x2m + x0Sin);
    r.set(2, 2, x2sqr * oneMinusCos + cs);

    return r;
}

// Convert a rotation matrix to Euler angles. Factorization into Euler angles
// is not necessarily unique. If the result is NOT_UNIQUE_SUM, then the
// multiple solutions occur because angleN2+angleN0 is constant. If the result
// is NOT_UNIQUE_DIF, then the multiple solutions occur because
// angleN2-angleN0 is constant. In either type of nonuniqueness, the function
// returns angleN0=0.
function convertMatrixToEulerAngles(r: Matrix, i0: number, i1: number,
    i2: number): EulerAngles {
    const e = new EulerAngles();
    e.axis = [i0, i1, i2];

    if (0 <= i0 && i0 < 3 && 0 <= i1 && i1 < 3 && 0 <= i2 && i2 < 3
        && i1 !== i0 && i1 !== i2) {
        if (i0 !== i2) {
            // Map (0,1,2), (1,2,0), and (2,0,1) to +1.
            // Map (0,2,1), (2,1,0), and (1,0,2) to -1.
            const parity = (((i2 | (i1 << 2)) >> i0) & 1);
            const sgn = ((parity & 1) !== 0 ? -1 : +1);

            if (r.get(i2, i0) < 1) {
                if (r.get(i2, i0) > -1) {
                    e.angle[2] = Math.atan2(sgn * r.get(i1, i0),
                        r.get(i0, i0));
                    e.angle[1] = Math.asin(-sgn * r.get(i2, i0));
                    e.angle[0] = Math.atan2(sgn * r.get(i2, i1),
                        r.get(i2, i2));
                    e.result = EulerResult.UNIQUE;
                } else {
                    e.angle[2] = 0;
                    e.angle[1] = sgn * GTE_C_HALF_PI;
                    e.angle[0] = Math.atan2(-sgn * r.get(i1, i2),
                        r.get(i1, i1));
                    e.result = EulerResult.NOT_UNIQUE_DIF;
                }
            } else {
                e.angle[2] = 0;
                e.angle[1] = -sgn * GTE_C_HALF_PI;
                e.angle[0] = Math.atan2(-sgn * r.get(i1, i2), r.get(i1, i1));
                e.result = EulerResult.NOT_UNIQUE_SUM;
            }
        } else {
            // Map (0,2,0), (1,0,1), and (2,1,2) to +1.
            // Map (0,1,0), (1,2,1), and (2,0,2) to -1.
            const b0 = 3 - i1 - i2;
            const parity = (((b0 | (i1 << 2)) >> i2) & 1);
            const sgn = ((parity & 1) !== 0 ? +1 : -1);

            if (r.get(i2, i2) < 1) {
                if (r.get(i2, i2) > -1) {
                    e.angle[2] = Math.atan2(r.get(i1, i2),
                        sgn * r.get(b0, i2));
                    e.angle[1] = Math.acos(r.get(i2, i2));
                    e.angle[0] = Math.atan2(r.get(i2, i1),
                        -sgn * r.get(i2, b0));
                    e.result = EulerResult.UNIQUE;
                } else {
                    e.angle[2] = 0;
                    e.angle[1] = GTE_C_PI;
                    e.angle[0] = Math.atan2(sgn * r.get(i1, b0),
                        r.get(i1, i1));
                    e.result = EulerResult.NOT_UNIQUE_DIF;
                }
            } else {
                e.angle[2] = 0;
                e.angle[1] = 0;
                e.angle[0] = Math.atan2(sgn * r.get(i1, b0), r.get(i1, i1));
                e.result = EulerResult.NOT_UNIQUE_SUM;
            }
        }
    } else {
        // Invalid angles.
        e.angle[0] = 0;
        e.angle[1] = 0;
        e.angle[2] = 0;
        e.result = EulerResult.INVALID;
    }

    return e;
}

// Convert Euler angles to a rotation matrix. The three integer inputs are in
// {0,1,2} and correspond to world directions unit(0) = (1,0,0),
// unit(1) = (0,1,0), or unit(2) = (0,0,1). The triples (N0,N1,N2) must be in
//   {(0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0),
//    (0,1,0),(0,2,0),(1,0,1),(1,2,1),(2,0,2),(2,1,2)}
// The rotation matrix is
//   R(unit(N2),angleN2)*R(unit(N1),angleN1)*R(unit(N0),angleN0)
// The conventions of convertAxisAngleToMatrix apply here as well.
//
// NOTE: Upstream's reversal of order is chosen so that a rotation matrix
// built with one multiplication convention is the transpose of the rotation
// matrix built with the other multiplication convention. With GTE_USE_MAT_VEC
// (the branch ported here), U = R*V is (u0,u1,u2) = R2*R1*R0*V.
function convertEulerAnglesToMatrix(e: EulerAngles, n: number): Matrix {
    assertDimension(n);

    if (0 <= e.axis[0] && e.axis[0] < 3
        && 0 <= e.axis[1] && e.axis[1] < 3
        && 0 <= e.axis[2] && e.axis[2] < 3
        && e.axis[1] !== e.axis[0]
        && e.axis[1] !== e.axis[2]) {
        const r0 = convertAxisAngleToMatrix(
            new AxisAngle(Vector.unit(n, e.axis[0]), e.angle[0]), n);
        const r1 = convertAxisAngleToMatrix(
            new AxisAngle(Vector.unit(n, e.axis[1]), e.angle[1]), n);
        const r2 = convertAxisAngleToMatrix(
            new AxisAngle(Vector.unit(n, e.axis[2]), e.angle[2]), n);
        return multiplyAB(multiplyAB(r2, r1), r0);
    } else {
        // Invalid angles.
        return Matrix.identity(n, n);
    }
}

// Convert a quaternion to an axis-angle pair, where
//   q = sin(angle/2)*(axis[0]*i+axis[1]*j+axis[2]*k)+cos(angle/2)
function convertQuaternionToAxisAngle(q: Quaternion, n: number): AxisAngle {
    assertDimension(n);
    const a = new AxisAngle(new Vector(n), 0);

    const axisSqrLen = q.values[0] * q.values[0] + q.values[1] * q.values[1]
        + q.values[2] * q.values[2];
    if (axisSqrLen > 0) {
        const adjust = 1 / Math.sqrt(axisSqrLen);
        a.axis.values[0] = q.values[0] * adjust;
        a.axis.values[1] = q.values[1] * adjust;
        a.axis.values[2] = q.values[2] * adjust;
        const cs = Math.max(Math.min(q.values[3], 1), -1);
        a.angle = 2 * Math.acos(cs);
    } else {
        // The angle is 0 (modulo 2*pi). Any axis will work, so choose the
        // Unit(0) axis.
        a.axis.values[0] = 1;
        a.angle = 0;
    }

    return a;
}

// Convert an axis-angle pair to a quaternion, where
//   q = sin(angle/2)*(axis[0]*i+axis[1]*j+axis[2]*k)+cos(angle/2)
function convertAxisAngleToQuaternion(a: AxisAngle): Quaternion {
    const halfAngle = 0.5 * a.angle;
    const sn = Math.sin(halfAngle);
    return new Quaternion(sn * a.axis.values[0], sn * a.axis.values[1],
        sn * a.axis.values[2], Math.cos(halfAngle));
}

function cloneEulerAngles(e: EulerAngles): EulerAngles {
    const copy = new EulerAngles();
    copy.axis = [e.axis[0], e.axis[1], e.axis[2]];
    copy.angle = [e.angle[0], e.angle[1], e.angle[2]];
    copy.result = e.result;
    return copy;
}

// The port of the private 'Rotation::Type' enumeration.
enum RotationType {
    IS_MATRIX,
    IS_QUATERNION,
    IS_AXIS_ANGLE,
    IS_EULER_ANGLES
}

export class Rotation {
    // The dimension N, which must be 3 or 4.
    readonly n: number;

    private mType: RotationType;
    private mMatrix: Matrix;
    private mQuaternion: Quaternion;
    private mAxisAngle: AxisAngle;
    private mEulerAngles: EulerAngles;

    private constructor(n: number, type: RotationType) {
        assertDimension(n);
        this.n = n;
        this.mType = type;
        this.mMatrix = Matrix.zero(n, n);
        this.mQuaternion = new Quaternion();
        this.mAxisAngle = new AxisAngle(new Vector(n), 0);
        this.mEulerAngles = new EulerAngles();
    }

    // Create rotations from various representations. The arguments are
    // copied, matching C++ value semantics.
    static fromMatrix(matrix: Matrix): Rotation {
        logAssert(matrix.numRows === matrix.numCols,
            'Rotation: expecting a square matrix.');
        const rotation = new Rotation(matrix.numRows, RotationType.IS_MATRIX);
        rotation.mMatrix = matrix.clone();
        return rotation;
    }

    static fromQuaternion(quaternion: Quaternion, n: number = 3): Rotation {
        const rotation = new Rotation(n, RotationType.IS_QUATERNION);
        rotation.mQuaternion = quaternion.clone();
        return rotation;
    }

    static fromAxisAngle(axisAngle: AxisAngle): Rotation {
        const rotation = new Rotation(axisAngle.axis.size,
            RotationType.IS_AXIS_ANGLE);
        rotation.mAxisAngle = new AxisAngle(axisAngle.axis, axisAngle.angle);
        return rotation;
    }

    static fromEulerAngles(eulerAngles: EulerAngles,
        n: number = 3): Rotation {
        const rotation = new Rotation(n, RotationType.IS_EULER_ANGLES);
        rotation.mEulerAngles = cloneEulerAngles(eulerAngles);
        return rotation;
    }

    // Convert one representation to another.
    toMatrix(): Matrix {
        switch (this.mType) {
            case RotationType.IS_MATRIX:
                break;
            case RotationType.IS_QUATERNION:
                this.mMatrix = convertQuaternionToMatrix(this.mQuaternion,
                    this.n);
                break;
            case RotationType.IS_AXIS_ANGLE:
                this.mMatrix = convertAxisAngleToMatrix(this.mAxisAngle,
                    this.n);
                break;
            case RotationType.IS_EULER_ANGLES:
                this.mMatrix = convertEulerAnglesToMatrix(this.mEulerAngles,
                    this.n);
                break;
        }

        return this.mMatrix.clone();
    }

    toQuaternion(): Quaternion {
        switch (this.mType) {
            case RotationType.IS_MATRIX:
                this.mQuaternion = convertMatrixToQuaternion(this.mMatrix);
                break;
            case RotationType.IS_QUATERNION:
                break;
            case RotationType.IS_AXIS_ANGLE:
                this.mQuaternion =
                    convertAxisAngleToQuaternion(this.mAxisAngle);
                break;
            case RotationType.IS_EULER_ANGLES:
                // Upstream converts the Euler angles to a matrix which is
                // then converted to a quaternion.
                this.mQuaternion = convertMatrixToQuaternion(
                    convertEulerAnglesToMatrix(this.mEulerAngles, this.n));
                break;
        }

        return this.mQuaternion.clone();
    }

    toAxisAngle(): AxisAngle {
        switch (this.mType) {
            case RotationType.IS_MATRIX:
                this.mAxisAngle = convertMatrixToAxisAngle(this.mMatrix,
                    this.n);
                break;
            case RotationType.IS_QUATERNION:
                this.mAxisAngle = convertQuaternionToAxisAngle(
                    this.mQuaternion, this.n);
                break;
            case RotationType.IS_AXIS_ANGLE:
                break;
            case RotationType.IS_EULER_ANGLES:
                // Upstream converts the Euler angles to a quaternion which is
                // then converted to an axis-angle pair.
                this.mAxisAngle = convertQuaternionToAxisAngle(
                    convertMatrixToQuaternion(convertEulerAnglesToMatrix(
                        this.mEulerAngles, this.n)), this.n);
                break;
        }

        return new AxisAngle(this.mAxisAngle.axis, this.mAxisAngle.angle);
    }

    // The port of 'operator()(i0, i1, i2)'. The axis indices select the
    // Euler-angle factorization.
    toEulerAngles(i0: number, i1: number, i2: number): EulerAngles {
        switch (this.mType) {
            case RotationType.IS_MATRIX:
                this.mEulerAngles = convertMatrixToEulerAngles(this.mMatrix,
                    i0, i1, i2);
                break;
            case RotationType.IS_QUATERNION:
                // Upstream converts the quaternion to a matrix which is then
                // converted to Euler angles.
                this.mEulerAngles = convertMatrixToEulerAngles(
                    convertQuaternionToMatrix(this.mQuaternion, this.n),
                    i0, i1, i2);
                break;
            case RotationType.IS_AXIS_ANGLE:
                // Upstream converts the axis-angle pair to a quaternion which
                // is then converted to Euler angles (and the quaternion to
                // Euler-angle conversion goes through a matrix).
                this.mEulerAngles = convertMatrixToEulerAngles(
                    convertQuaternionToMatrix(
                        convertAxisAngleToQuaternion(this.mAxisAngle),
                        this.n),
                    i0, i1, i2);
                break;
            case RotationType.IS_EULER_ANGLES:
                if (this.mEulerAngles.axis[0] === i0
                    && this.mEulerAngles.axis[1] === i1
                    && this.mEulerAngles.axis[2] === i2) {
                    // The requested factorization is the one the object was
                    // constructed with; return it unchanged, as upstream.
                    break;
                }
                // UPSTREAM BUG (Rotation.h, 'operator()(i0,i1,i2)'): upstream
                // assigns the requested axis indices to the cached Euler
                // angles and then does nothing in the IS_EULER_ANGLES case,
                // so a request for a *different* factorization silently
                // returns the original angles relabeled with the new axes -
                // a different rotation. The port recomputes through the
                // rotation matrix, which is what every other source type
                // does.
                this.mEulerAngles = convertMatrixToEulerAngles(
                    convertEulerAnglesToMatrix(this.mEulerAngles, this.n),
                    i0, i1, i2);
                break;
        }

        return cloneEulerAngles(this.mEulerAngles);
    }
}
