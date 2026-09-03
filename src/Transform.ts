// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Transform.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The transform is Y = M*X+T, where M is a 3-by-3 matrix and T is a 3x1
// translation. In most cases, M = R, a rotation matrix, or M = R*S, where R
// is a rotation matrix and S is a diagonal matrix whose diagonal entries are
// positive scales. To support modeling packages that allow general affine
// transforms, M can be any invertible 3x3 matrix. The vector X is transformed
// in the "forward" direction to Y. The "inverse" direction transforms Y to X,
// namely X = M^{-1}*(Y-T) in the general case. In the special case of
// M = R*S, the inverse direction is X = S^{-1}*R^t*(Y-T), where S^{-1} is the
// diagonal matrix whose diagonal entries are the reciprocals of those of S
// and where R^t is the transpose of R. A homogeneous matrix
// H = {{M,T},{0,1}} is stored by this class. The forward transform is
// {Y,1} = H*{X,1} and the inverse transform is {X,1} = H^{-1}*{Y,1}.
//
// A matrix M = R*S is referred to as an "RS-matrix". The class does not
// provide a member function to compute the inverse of a transform as a
// transform channel-by-channel identity, because the inverse of an RS-matrix
// is not generally an RS-matrix; that is, the inverse of R*S is S^{-1}*R^t
// which cannot always be factored as S^{-1}*R^t = R'*S'. The 'inverse'
// member handles this by falling back to a general matrix.
//
// Port notes:
// - Only the GTE_USE_MAT_VEC branches are ported, matching the precedent set
//   by src/Matrix4x4.ts (B56).
// - Upstream's overloaded 'SetRotation' becomes a single 'setRotation' that
//   dispatches on the runtime type of its argument (a 3x3 or 4x4 Matrix, a
//   Quaternion, an AxisAngle of dimension 3 or 4, or EulerAngles). The
//   overloaded 'GetRotation' out-parameter forms become separate named
//   getters that return their results.
// - Accessors that upstream returns by const reference return clones here so
//   that callers cannot corrupt the transform's internal state.
// - Upstream's 'operator Matrix4x4<Real> const&' implicit conversion has no
//   TypeScript equivalent; use 'getHMatrix()'.

import { logAssert } from './Logger.js';
import {
    Matrix, mulMatrix, multiplyMD, transpose
} from './Matrix.js';
import { inverse4x4 } from './Matrix4x4.js';
import { Quaternion } from './Quaternion.js';
import { AxisAngle } from './AxisAngle.js';
import { EulerAngles } from './EulerAngles.js';
import { Rotation } from './Rotation.js';
import { Vector, add, mul, hlift, hproject } from './Vector.js';

function vector4(x0: number, x1: number, x2: number, x3: number): Vector {
    return Vector.fromArray([x0, x1, x2, x3]);
}

export class Transform {
    // The full 4x4 homogeneous matrix H and its inverse H^{-1}. The inverse
    // is computed only on demand.
    private mHMatrix: Matrix;
    private mInvHMatrix: Matrix;

    private mMatrix: Matrix;    // M (general) or R (rotation)
    private mTranslate: Vector; // T, stored as (t0,t1,t2,1)
    private mScale: Vector;     // S, stored as (s0,s1,s2,1)
    private mIsIdentity: boolean;
    private mIsRSMatrix: boolean;
    private mIsUniformScale: boolean;
    private mInverseNeedsUpdate: boolean;

    // The default constructor produces the identity transformation.
    constructor() {
        this.mHMatrix = Matrix.identity(4, 4);
        this.mInvHMatrix = Matrix.identity(4, 4);
        this.mMatrix = Matrix.identity(4, 4);
        this.mTranslate = vector4(0, 0, 0, 1);
        this.mScale = vector4(1, 1, 1, 1);
        this.mIsIdentity = true;
        this.mIsRSMatrix = true;
        this.mIsUniformScale = true;
        this.mInverseNeedsUpdate = false;
    }

    // The port of the compiler-generated copy constructor / assignment
    // operator; C++ has value semantics where TypeScript objects alias.
    clone(): Transform {
        const copy = new Transform();
        copy.mHMatrix = this.mHMatrix.clone();
        copy.mInvHMatrix = this.mInvHMatrix.clone();
        copy.mMatrix = this.mMatrix.clone();
        copy.mTranslate = this.mTranslate.clone();
        copy.mScale = this.mScale.clone();
        copy.mIsIdentity = this.mIsIdentity;
        copy.mIsRSMatrix = this.mIsRSMatrix;
        copy.mIsUniformScale = this.mIsUniformScale;
        copy.mInverseNeedsUpdate = this.mInverseNeedsUpdate;
        return copy;
    }

    // Set the transformation to the identity matrix.
    makeIdentity(): void {
        this.mMatrix.makeIdentity();
        this.mTranslate = vector4(0, 0, 0, 1);
        this.mScale = vector4(1, 1, 1, 1);
        this.mIsIdentity = true;
        this.mIsRSMatrix = true;
        this.mIsUniformScale = true;
        this.updateHMatrix();
    }

    // Set the transformation to have scales of 1.
    makeUnitScale(): void {
        logAssert(this.mIsRSMatrix, 'Transform is not rotation-scale.');
        this.mScale = vector4(1, 1, 1, 1);
        this.mIsUniformScale = true;
        this.updateHMatrix();
    }

    // Hints about the structure of the transformation.

    // M = I
    isIdentity(): boolean {
        return this.mIsIdentity;
    }

    // M = R*S
    isRSMatrix(): boolean {
        return this.mIsRSMatrix;
    }

    // RS-matrix with S = c*I
    isUniformScale(): boolean {
        return this.mIsRSMatrix && this.mIsUniformScale;
    }

    // Member access.
    // (1) The set* functions set the is-identity hint to false.
    // (2) The setRotation function sets the is-rsmatrix hint to true. If this
    //     hint is false, the get*-rotation functions throw.
    // (3) The setMatrix function sets the is-rsmatrix and is-uniform-scale
    //     hints to false.
    // (4) The setScale function sets the is-uniform-scale hint to false. The
    //     setUniformScale function sets the is-uniform-scale hint to true. If
    //     this hint is false, getUniformScale throws.
    // (5) All set* functions set the inverse-needs-update to true. When
    //     getHInverse is called, the inverse must be computed in this case
    //     and the inverse-needs-update is reset to false.

    // {{R,0},{0,1}}. The argument is a 4x4 or 3x3 rotation matrix, a unit
    // quaternion, an axis-angle pair of dimension 3 or 4 (the axis is unit
    // length and the angle is in radians), or Euler angles (in radians).
    setRotation(rotate: Matrix | Quaternion | AxisAngle | EulerAngles): void {
        if (rotate instanceof Matrix) {
            if (rotate.numRows === 4 && rotate.numCols === 4) {
                this.mMatrix = rotate.clone();
            } else {
                logAssert(rotate.numRows === 3 && rotate.numCols === 3,
                    'The rotation matrix must be 3x3 or 4x4.');
                this.mMatrix.makeIdentity();
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 3; ++c) {
                        this.mMatrix.set(r, c, rotate.get(r, c));
                    }
                }
            }
        } else if (rotate instanceof Quaternion) {
            this.mMatrix = Rotation.fromQuaternion(rotate, 4).toMatrix();
        } else if (rotate instanceof AxisAngle) {
            if (rotate.axis.size === 3) {
                // Upstream lifts the 3D axis with a last component of 1 even
                // though AxisAngle<4> is documented to have axis (x,y,z,0).
                // The conversion to a matrix ignores the last component, so
                // the resulting rotation is the same; see the PR notes.
                const aa4 = new AxisAngle(hlift(rotate.axis, 1), rotate.angle);
                this.mMatrix = Rotation.fromAxisAngle(aa4).toMatrix();
            } else {
                this.mMatrix = Rotation.fromAxisAngle(rotate).toMatrix();
            }
        } else {
            this.mMatrix = Rotation.fromEulerAngles(rotate, 4).toMatrix();
        }

        this.mIsIdentity = false;
        this.mIsRSMatrix = true;
        this.updateHMatrix();
    }

    // {{M,0},{0,1}}
    setMatrix(matrix: Matrix): void {
        logAssert(matrix.numRows === 4 && matrix.numCols === 4,
            'The matrix must be 4x4.');
        this.mMatrix = matrix.clone();
        this.mIsIdentity = false;
        this.mIsRSMatrix = false;
        this.mIsUniformScale = false;
        this.updateHMatrix();
    }

    setTranslation(x0: number, x1: number, x2: number): void;
    setTranslation(translate: Vector): void;
    setTranslation(arg0: number | Vector, x1?: number, x2?: number): void {
        if (arg0 instanceof Vector) {
            logAssert(arg0.size === 3 || arg0.size === 4,
                'The translation must be a 3-tuple or 4-tuple.');
            this.mTranslate = vector4(arg0.values[0], arg0.values[1],
                arg0.values[2], 1);
        } else {
            this.mTranslate = vector4(arg0, x1 as number, x2 as number, 1);
        }
        this.mIsIdentity = false;
        this.updateHMatrix();
    }

    setScale(s0: number, s1: number, s2: number): void;
    setScale(scale: Vector): void;
    setScale(arg0: number | Vector, s1?: number, s2?: number): void {
        logAssert(this.mIsRSMatrix, 'Transform is not rotation-scale.');
        let v0: number, v1: number, v2: number;
        if (arg0 instanceof Vector) {
            logAssert(arg0.size === 3 || arg0.size === 4,
                'The scale must be a 3-tuple or 4-tuple.');
            v0 = arg0.values[0];
            v1 = arg0.values[1];
            v2 = arg0.values[2];
        } else {
            v0 = arg0;
            v1 = s1 as number;
            v2 = s2 as number;
        }
        logAssert(v0 !== 0 && v1 !== 0 && v2 !== 0, 'Scales must be nonzero.');
        this.mScale = vector4(v0, v1, v2, 1);
        this.mIsIdentity = false;
        this.mIsUniformScale = false;
        this.updateHMatrix();
    }

    setUniformScale(scale: number): void {
        logAssert(this.mIsRSMatrix, 'Transform is not rotation-scale.');
        logAssert(scale !== 0, 'Scale must be nonzero.');
        this.mScale = vector4(scale, scale, scale, 1);
        this.mIsIdentity = false;
        this.mIsUniformScale = true;
        this.updateHMatrix();
    }

    // {{R,0},{0,1}}
    getRotation(): Matrix {
        logAssert(this.mIsRSMatrix, 'Transform is not rotation-scale.');
        return this.mMatrix.clone();
    }

    // {{M,0},{0,1}}
    getMatrix(): Matrix {
        return this.mMatrix.clone();
    }

    // (x,y,z)
    getTranslation(): Vector {
        return Vector.fromArray([this.mTranslate.values[0],
            this.mTranslate.values[1], this.mTranslate.values[2]]);
    }

    // (x,y,z,0)
    getTranslationW0(): Vector {
        return vector4(this.mTranslate.values[0], this.mTranslate.values[1],
            this.mTranslate.values[2], 0);
    }

    // (x,y,z,1)
    getTranslationW1(): Vector {
        return vector4(this.mTranslate.values[0], this.mTranslate.values[1],
            this.mTranslate.values[2], 1);
    }

    // (s0,s1,s2)
    getScale(): Vector {
        logAssert(this.mIsRSMatrix, 'Transform is not rotation-scale.');
        return Vector.fromArray([this.mScale.values[0],
            this.mScale.values[1], this.mScale.values[2]]);
    }

    // (s0,s1,s2,1)
    getScaleW1(): Vector {
        logAssert(this.mIsRSMatrix, 'Transform is not rotation-scale.');
        return vector4(this.mScale.values[0], this.mScale.values[1],
            this.mScale.values[2], 1);
    }

    getUniformScale(): number {
        logAssert(this.mIsRSMatrix, 'Transform is not rotation-scale.');
        logAssert(this.mIsUniformScale, 'Transform is not uniform scale.');
        return this.mScale.values[0];
    }

    // Alternate representations to get the rotation.

    // The upper-left 3x3 block of the rotation matrix.
    getRotationMatrix3x3(): Matrix {
        logAssert(this.mIsRSMatrix, 'Transform is not rotation-scale.');
        const rotate = new Matrix(3, 3);
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 3; ++c) {
                rotate.set(r, c, this.mMatrix.get(r, c));
            }
        }
        return rotate;
    }

    // The quaternion is unit length.
    getRotationQuaternion(): Quaternion {
        logAssert(this.mIsRSMatrix, 'Transform is not rotation-scale.');
        return Rotation.fromMatrix(this.mMatrix).toQuaternion();
    }

    // The axis is unit length and the angle is in radians. The axis is a
    // 3-tuple.
    getRotationAxisAngle3(): AxisAngle {
        logAssert(this.mIsRSMatrix, 'Transform is not rotation-scale.');
        const aa4 = Rotation.fromMatrix(this.mMatrix).toAxisAngle();
        return new AxisAngle(hproject(aa4.axis), aa4.angle);
    }

    // The axis is unit length and the angle is in radians. The axis is a
    // 4-tuple (x,y,z,0).
    getRotationAxisAngle4(): AxisAngle {
        logAssert(this.mIsRSMatrix, 'Transform is not rotation-scale.');
        return Rotation.fromMatrix(this.mMatrix).toAxisAngle();
    }

    // The Euler angles are in radians. The (i0,i1,i2) triple is the axis
    // order of the desired factorization. Upstream passes an EulerAngles
    // object whose 'axis' member selects the order and whose 'angle' member
    // is filled in; the port takes the order directly and returns the result.
    getRotationEulerAngles(i0: number, i1: number, i2: number): EulerAngles {
        logAssert(this.mIsRSMatrix, 'Transform is not rotation-scale.');
        return Rotation.fromMatrix(this.mMatrix).toEulerAngles(i0, i1, i2);
    }

    // For M = R*S, the largest value of S in absolute value is returned. For
    // general M, the max-row-sum norm is returned, which is a reasonable
    // measure of maximum scale of the transformation.
    getNorm(): number {
        let sum0: number, sum1: number, sum2: number;

        if (this.mIsRSMatrix) {
            sum0 = Math.abs(this.mScale.values[0]);
            sum1 = Math.abs(this.mScale.values[1]);
            sum2 = Math.abs(this.mScale.values[2]);
        } else {
            // The spectral norm (the maximum absolute value of the
            // eigenvalues) is smaller than or equal to this norm. Therefore,
            // this function returns an approximation to the maximum scale.
            const m = this.mMatrix;
            sum0 = Math.abs(m.get(0, 0)) + Math.abs(m.get(0, 1))
                + Math.abs(m.get(0, 2));
            sum1 = Math.abs(m.get(1, 0)) + Math.abs(m.get(1, 1))
                + Math.abs(m.get(1, 2));
            sum2 = Math.abs(m.get(2, 0)) + Math.abs(m.get(2, 1))
                + Math.abs(m.get(2, 2));
        }

        return Math.max(Math.max(sum0, sum1), sum2);
    }

    // Get the homogeneous matrix (composite of all channels).
    getHMatrix(): Matrix {
        return this.mHMatrix.clone();
    }

    // Get the inverse homogeneous matrix, recomputing it when necessary. For
    // H = {{M,T},{0,1}}, then H^{-1} = {{M^{-1},-M^{-1}*T},{0,1}}.
    getHInverse(): Matrix {
        if (this.mInverseNeedsUpdate) {
            if (this.mIsIdentity) {
                this.mInvHMatrix.makeIdentity();
            } else {
                const m = this.mMatrix;
                const inv = this.mInvHMatrix;
                if (this.mIsRSMatrix) {
                    if (this.mIsUniformScale) {
                        const invScale = 1 / this.mScale.values[0];
                        inv.set(0, 0, invScale * m.get(0, 0));
                        inv.set(0, 1, invScale * m.get(1, 0));
                        inv.set(0, 2, invScale * m.get(2, 0));
                        inv.set(1, 0, invScale * m.get(0, 1));
                        inv.set(1, 1, invScale * m.get(1, 1));
                        inv.set(1, 2, invScale * m.get(2, 1));
                        inv.set(2, 0, invScale * m.get(0, 2));
                        inv.set(2, 1, invScale * m.get(1, 2));
                        inv.set(2, 2, invScale * m.get(2, 2));
                    } else {
                        // Replace 3 reciprocals by 6 multiplies and
                        // 1 reciprocal.
                        const s01 = this.mScale.values[0]
                            * this.mScale.values[1];
                        const s02 = this.mScale.values[0]
                            * this.mScale.values[2];
                        const s12 = this.mScale.values[1]
                            * this.mScale.values[2];
                        const invs012 = 1 / (s01 * this.mScale.values[2]);
                        const invS0 = s12 * invs012;
                        const invS1 = s02 * invs012;
                        const invS2 = s01 * invs012;
                        inv.set(0, 0, invS0 * m.get(0, 0));
                        inv.set(0, 1, invS0 * m.get(1, 0));
                        inv.set(0, 2, invS0 * m.get(2, 0));
                        inv.set(1, 0, invS1 * m.get(0, 1));
                        inv.set(1, 1, invS1 * m.get(1, 1));
                        inv.set(1, 2, invS1 * m.get(2, 1));
                        inv.set(2, 0, invS2 * m.get(0, 2));
                        inv.set(2, 1, invS2 * m.get(1, 2));
                        inv.set(2, 2, invS2 * m.get(2, 2));
                    }

                    // Compute the translation part of the inverse from the
                    // 3x3 block that was just set.
                    inv.set(0, 3, -(
                        inv.get(0, 0) * this.mTranslate.values[0] +
                        inv.get(0, 1) * this.mTranslate.values[1] +
                        inv.get(0, 2) * this.mTranslate.values[2]));
                    inv.set(1, 3, -(
                        inv.get(1, 0) * this.mTranslate.values[0] +
                        inv.get(1, 1) * this.mTranslate.values[1] +
                        inv.get(1, 2) * this.mTranslate.values[2]));
                    inv.set(2, 3, -(
                        inv.get(2, 0) * this.mTranslate.values[0] +
                        inv.get(2, 1) * this.mTranslate.values[1] +
                        inv.get(2, 2) * this.mTranslate.values[2]));

                    // UPSTREAM BUG (Transform.h, GetHInverse): the RS branch
                    // writes only the upper-left 3x3 block and the last
                    // column, relying on the last row of mInvHMatrix still
                    // being (0,0,0,1) from the constructor. That assumption
                    // is violated once the non-RS branch has assigned
                    // 'Inverse(mHMatrix)' for a singular mHMatrix, which
                    // yields the zero matrix and leaves the last row all
                    // zeros. A later switch back to an RS transform then
                    // produces an H^{-1} with last row (0,0,0,0). The port
                    // sets the last row explicitly.
                    inv.set(3, 0, 0);
                    inv.set(3, 1, 0);
                    inv.set(3, 2, 0);
                    inv.set(3, 3, 1);
                } else {
                    // gte::Inverse produces the correct inverse including the
                    // translation, so no overwrite of the last column is
                    // needed here.
                    this.mInvHMatrix = inverse4x4(this.mHMatrix).inverse;
                }
            }

            this.mInverseNeedsUpdate = false;
        }

        return this.mInvHMatrix.clone();
    }

    // Invert the transform. If possible, the channels are properly assigned.
    // For example, if the input has is-rsmatrix equal to true and the scale
    // is uniform, then the inverse also has is-rsmatrix equal to true, its
    // matrix is a rotation matrix and its scale is set accordingly.
    inverse(): Transform {
        const inverse = new Transform(); // = the identity

        if (!this.mIsIdentity) {
            if (this.mIsRSMatrix && this.mIsUniformScale) {
                const invRotate = transpose(this.getRotation());
                const invScale = 1 / this.getUniformScale();
                const invTranslate = mul(-invScale,
                    mulMatrix(invRotate, this.getTranslationW1()));
                inverse.setRotation(invRotate);
                inverse.setUniformScale(invScale);
                inverse.setTranslation(invTranslate);
            } else {
                const invMatrix = inverse4x4(this.mHMatrix).inverse;
                const invTranslate = invMatrix.getCol(3);
                // UPSTREAM BUG (Transform.h, Inverse): upstream passes the
                // full 4x4 inverse to SetMatrix, so the inverse transform's
                // 'M' channel carries the translation in its last column even
                // though GetMatrix is documented to return {{M,0},{0,1}}. The
                // homogeneous matrix is unaffected (UpdateHMatrix copies only
                // the 3x3 block), but the accessor is wrong. The port clears
                // the last column and row before storing M.
                for (let i = 0; i < 3; ++i) {
                    invMatrix.set(i, 3, 0);
                    invMatrix.set(3, i, 0);
                }
                invMatrix.set(3, 3, 1);
                inverse.setMatrix(invMatrix);
                inverse.setTranslation(invTranslate);
            }
        }

        return inverse;
    }

    // The identity transformation.
    static identity(): Transform {
        return new Transform();
    }

    // Fill in the entries of mHMatrix whenever one of the components
    // mMatrix, mTranslate or mScale changes.
    private updateHMatrix(): void {
        if (this.mIsIdentity) {
            this.mHMatrix.makeIdentity();
        } else {
            const h = this.mHMatrix;
            const m = this.mMatrix;
            if (this.mIsRSMatrix) {
                h.set(0, 0, m.get(0, 0) * this.mScale.values[0]);
                h.set(0, 1, m.get(0, 1) * this.mScale.values[1]);
                h.set(0, 2, m.get(0, 2) * this.mScale.values[2]);
                h.set(1, 0, m.get(1, 0) * this.mScale.values[0]);
                h.set(1, 1, m.get(1, 1) * this.mScale.values[1]);
                h.set(1, 2, m.get(1, 2) * this.mScale.values[2]);
                h.set(2, 0, m.get(2, 0) * this.mScale.values[0]);
                h.set(2, 1, m.get(2, 1) * this.mScale.values[1]);
                h.set(2, 2, m.get(2, 2) * this.mScale.values[2]);
            } else {
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 3; ++c) {
                        h.set(r, c, m.get(r, c));
                    }
                }
            }

            h.set(0, 3, this.mTranslate.values[0]);
            h.set(1, 3, this.mTranslate.values[1]);
            h.set(2, 3, this.mTranslate.values[2]);

            // The last row of mHMatrix is always (0,0,0,1) for an affine
            // transformation, so it is set once in the constructor. It is not
            // necessary to reset it here.
        }

        this.mInverseNeedsUpdate = true;
    }
}

// Compute M*V, V^T*M, A*B for transforms A and B, and the mixed products of
// a transform with a 4x4 matrix. The port of upstream's 'operator*'
// overloads.
export function mulTransform(M: Transform, V: Vector): Vector;
export function mulTransform(V: Vector, M: Transform): Vector;
export function mulTransform(A: Transform, B: Transform): Transform;
export function mulTransform(A: Matrix, B: Transform): Matrix;
export function mulTransform(A: Transform, B: Matrix): Matrix;
export function mulTransform(arg0: Transform | Vector | Matrix,
    arg1: Transform | Vector | Matrix): Transform | Vector | Matrix {
    if (arg0 instanceof Transform && arg1 instanceof Transform) {
        return mulTransformTransform(arg0, arg1);
    }

    if (arg0 instanceof Transform) {
        const h = arg0.getHMatrix();
        return arg1 instanceof Matrix
            ? mulMatrix(h, arg1)
            : mulMatrix(h, arg1 as Vector);
    }

    const h = (arg1 as Transform).getHMatrix();
    return arg0 instanceof Matrix
        ? mulMatrix(arg0, h)
        : mulMatrix(arg0 as Vector, h);
}

function mulTransformTransform(A: Transform, B: Transform): Transform {
    if (A.isIdentity()) {
        return B.clone();
    }

    if (B.isIdentity()) {
        return A.clone();
    }

    const product = new Transform();

    if (A.isRSMatrix() && B.isRSMatrix()) {
        if (A.isUniformScale()) {
            // A*B = (a*R_A)*(R_B*S_B)*X + (a*R_A*T_B + T_A), which is an
            // RS-matrix with rotation R_A*R_B and scale a*S_B.
            product.setRotation(mulMatrix(A.getRotation(), B.getRotation()));

            product.setTranslation(add(
                mul(A.getUniformScale(),
                    mulMatrix(A.getRotation(), B.getTranslationW0())),
                A.getTranslationW1()));

            if (B.isUniformScale()) {
                product.setUniformScale(
                    A.getUniformScale() * B.getUniformScale());
            } else {
                product.setScale(mul(A.getUniformScale(), B.getScale()));
            }

            return product;
        }
    }

    // In all remaining cases, the matrix cannot be written as R*S*X+T.
    const matMA = A.isRSMatrix()
        ? multiplyMD(A.getRotation(), A.getScaleW1())
        : A.getMatrix();

    const matMB = B.isRSMatrix()
        ? multiplyMD(B.getRotation(), B.getScaleW1())
        : B.getMatrix();

    product.setMatrix(mulMatrix(matMA, matMB));
    product.setTranslation(add(mulMatrix(matMA, B.getTranslationW0()),
        A.getTranslationW1()));
    return product;
}
