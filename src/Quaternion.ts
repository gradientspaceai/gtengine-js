// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Quaternion.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A quaternion is of the form
//   q = x * i + y * j + z * k + w * 1 = x * i + y * j + z * k + w
// where w, x, y, and z are real numbers. The scalar and vector parts are
//   Vector(q) = x * i + y * j + z * k
//   Scalar(q) = w
//   q = Vector(q) + Scalar(q)
// See https://www.geometrictools.com/Documentation/Quaternions.pdf for the
// arithmetic and algebraic properties of quaternions.
//
// The rotate() function requires fewer arithmetic operations than the
// original implementation using rotatedU = q * (0,u) * conjugate(q). The new
// implementation is based on Robert Eisele's derivation at
// https://raw.org/proof/vector-rotation-using-quaternions/
//
// Port notes:
// - Following PORTING.md's "reuse via subclassing" precedent, Quaternion
//   extends Vector as a fixed size-4 tuple in the order (x,y,z,w). Vector's
//   free functions therefore apply unchanged and are NOT duplicated here:
//   dot(q0,q1), length(q) and normalize(q) are exactly upstream's Dot,
//   Length and Normalize for quaternions (normalize mutates in place and
//   returns the length, and zeroes the tuple when the length is 0, as
//   upstream). The element comparisons (equals, lessThan, ...) are
//   inherited too.
// - The upstream operators that must keep the Quaternion type are ported as
//   free functions suffixed with the type context, matching Matrix.ts's
//   negateMatrix/addMatrix/... scheme: negateQuaternion, addQuaternion,
//   subQuaternion, mulQuaternion, divQuaternion. mulQuaternion is the only
//   one whose mathematics differs from the Vector version: for two
//   quaternion arguments it is the (noncommutative) Hamilton product.
//   Compound assignments (+=, -=, ...) have no in-place ports; unary
//   'operator+' is the identity and has no port.
// - Upstream's Inverse and Conjugate become inverseQuaternion (the name
//   'inverse' is owned by Matrix.ts) and conjugate.
// - Slerp: upstream's SlerpR, SlerpRP and SlerpRPH are algorithmically
//   identical to the already-ported Slerp.h functions slerp,
//   slerpUsingCosAngle and slerpUsingMidpoint (which take number[] tuples,
//   so pass q.values), and are not duplicated. Upstream's Quaternion::Slerp
//   is NOT the same as Slerp.h's slerp: it additionally handles an obtuse
//   angle by negating the contribution of q1, so it is ported here as
//   slerpQuaternion.
// - The port uses GTE_USE_MAT_VEC (the GTE default), so rotate(q,u) is the
//   vector R*u where R is the rotation matrix corresponding to q.

import { chebyshevRatiosUsingCosAngle } from './ChebyshevRatio';
import { logAssert } from './Logger';
import { Vector, add, mul } from './Vector';
import { cross } from './Vector3';

export class Quaternion extends Vector {
    // The quaternions are of the form q = x*i + y*j + z*k + w. In tuple
    // form, q = (x,y,z,w). Upstream's default constructor zeroes the tuple.
    constructor(x: number = 0, y: number = 0, z: number = 0, w: number = 0) {
        super(4);
        this.values[0] = x;
        this.values[1] = y;
        this.values[2] = z;
        this.values[3] = w;
    }

    // Create a quaternion that copies the input values, which must be the
    // four components (x,y,z,w).
    static override fromArray(values: readonly number[]): Quaternion {
        logAssert(values.length === 4, 'Quaternion: expecting 4 components.');
        return new Quaternion(values[0], values[1], values[2], values[3]);
    }

    override clone(): Quaternion {
        return new Quaternion(this.values[0], this.values[1], this.values[2],
            this.values[3]);
    }

    // Special quaternions.

    // z = 0*i + 0*j + 0*k + 0
    static override zero(): Quaternion {
        return new Quaternion(0, 0, 0, 0);
    }

    // i = 1*i + 0*j + 0*k + 0
    static i(): Quaternion {
        return new Quaternion(1, 0, 0, 0);
    }

    // j = 0*i + 1*j + 0*k + 0
    static j(): Quaternion {
        return new Quaternion(0, 1, 0, 0);
    }

    // k = 0*i + 0*j + 1*k + 0
    static k(): Quaternion {
        return new Quaternion(0, 0, 1, 0);
    }

    // 1 = 0*i + 0*j + 0*k + 1
    static identity(): Quaternion {
        return new Quaternion(0, 0, 0, 1);
    }
}

// Unary operations. The port of unary 'operator-'.
export function negateQuaternion(q: Quaternion): Quaternion {
    return new Quaternion(-q.values[0], -q.values[1], -q.values[2],
        -q.values[3]);
}

// Linear algebraic operations.
export function addQuaternion(q0: Quaternion, q1: Quaternion): Quaternion {
    return new Quaternion(
        q0.values[0] + q1.values[0],
        q0.values[1] + q1.values[1],
        q0.values[2] + q1.values[2],
        q0.values[3] + q1.values[3]);
}

export function subQuaternion(q0: Quaternion, q1: Quaternion): Quaternion {
    return new Quaternion(
        q0.values[0] - q1.values[0],
        q0.values[1] - q1.values[1],
        q0.values[2] - q1.values[2],
        q0.values[3] - q1.values[3]);
}

// The port of upstream's 'operator*' overloads: 'q * scalar',
// 'scalar * q' and the quaternion product 'q0 * q1'.
//
// Multiplication of quaternions is not generally commutative; that is,
// q0*q1 and q1*q0 are not usually the same value.
// (x0*i + y0*j + z0*k + w0)*(x1*i + y1*j + z1*k + w1)
// =
// i*(+x0*w1 + y0*z1 - z0*y1 + w0*x1) +
// j*(-x0*z1 + y0*w1 + z0*x1 + w0*y1) +
// k*(+x0*y1 - y0*x1 + z0*w1 + w0*z1) +
// 1*(-x0*x1 - y0*y1 - z0*z1 + w0*w1)
export function mulQuaternion(q: Quaternion, scalar: number): Quaternion;
export function mulQuaternion(scalar: number, q: Quaternion): Quaternion;
export function mulQuaternion(q0: Quaternion, q1: Quaternion): Quaternion;
export function mulQuaternion(arg0: Quaternion | number,
    arg1: Quaternion | number): Quaternion {
    if (typeof arg0 === 'number' || typeof arg1 === 'number') {
        const q = (typeof arg0 === 'number' ? arg1 : arg0) as Quaternion;
        const scalar = (typeof arg0 === 'number' ? arg0 : arg1) as number;
        return new Quaternion(q.values[0] * scalar, q.values[1] * scalar,
            q.values[2] * scalar, q.values[3] * scalar);
    }

    const q0 = arg0.values;
    const q1 = arg1.values;
    return new Quaternion(
        +q0[0] * q1[3] + q0[1] * q1[2] - q0[2] * q1[1] + q0[3] * q1[0],
        -q0[0] * q1[2] + q0[1] * q1[3] + q0[2] * q1[0] + q0[3] * q1[1],
        +q0[0] * q1[1] - q0[1] * q1[0] + q0[2] * q1[3] + q0[3] * q1[2],
        -q0[0] * q1[0] - q0[1] * q1[1] - q0[2] * q1[2] + q0[3] * q1[3]);
}

// The port of 'q / scalar'. As upstream's operator/=, division by zero
// produces the zero quaternion.
export function divQuaternion(q: Quaternion, scalar: number): Quaternion {
    if (scalar !== 0) {
        return new Quaternion(q.values[0] / scalar, q.values[1] / scalar,
            q.values[2] / scalar, q.values[3] / scalar);
    }
    return Quaternion.zero();
}

// The conjugate of q = (x,y,z,w) is conj(q) = (-x,-y,-z,w).
export function conjugate(q: Quaternion): Quaternion {
    return new Quaternion(-q.values[0], -q.values[1], -q.values[2],
        +q.values[3]);
}

// For a nonzero quaternion q = (x,y,z,w), inv(q) = (-x,-y,-z,w)/|q|^2, where
// |q| is the length of the quaternion. When q is zero, the function returns
// zero, which is considered to be an improbable case.
export function inverseQuaternion(q: Quaternion): Quaternion {
    const sqrLen = q.values[0] * q.values[0] + q.values[1] * q.values[1]
        + q.values[2] * q.values[2] + q.values[3] * q.values[3];
    if (sqrLen > 0) {
        return divQuaternion(conjugate(q), sqrLen);
    }
    return Quaternion.zero();
}

// Rotate a 3D vector u = (u0,u1,u2), or a 3D vector represented as a
// homogeneous 4D vector u = (u0,u1,u2,0), using quaternion multiplication.
// The input quaternion must be unit length. If R is the rotation matrix
// corresponding to the quaternion q, the rotated vector is R*u (the port
// uses GTE_USE_MAT_VEC).
export function rotate(q: Quaternion, u: Vector): Vector {
    logAssert(u.size === 3 || u.size === 4, 'Dimension must be 3 or 4.');

    const v = new Vector(u.size);
    v.values[0] = q.values[0];
    v.values[1] = q.values[1];
    v.values[2] = q.values[2];

    const t = mul(cross(v, u), 2);
    return add(add(u, mul(t, q.values[3])), cross(v, t));
}

// The spherical linear interpolation (slerp) of unit-length quaternions q0
// and q1 for t in [0,1] is
//     slerp(t,q0,q1) = [sin((1-t)*theta)*q0 + sin(t*theta)*q1]/sin(theta)
// where theta is the angle between q0 and q1 [cos(theta) = Dot(q0,q1)].
// This function is a parameterization of the great spherical arc between q0
// and q1 on the unit hypersphere. Moreover, the parameterization is one of
// normalized arclength: a particle traveling along the arc through time t
// does so with constant speed.
//
// The angle between q0 and q1 is in [0,pi). There are no angle restrictions
// and nothing is precomputed; when the angle is obtuse, q1 is negated first
// so that the interpolation takes the shorter arc (q and -q represent the
// same rotation). This is what distinguishes it from Slerp.ts's slerp.
//
// Upstream's SlerpR, SlerpRP and SlerpRPH (which require the angle to be in
// [0,pi/2] after preprocessing) are the already-ported slerp,
// slerpUsingCosAngle and slerpUsingMidpoint from Slerp.ts; call them with
// q.values. See SlerpEstimate.ts for the approximations.
export function slerpQuaternion(t: number, q0: Quaternion,
    q1: Quaternion): Quaternion {
    let cosA = q0.values[0] * q1.values[0] + q0.values[1] * q1.values[1]
        + q0.values[2] * q1.values[2] + q0.values[3] * q1.values[3];
    let sign: number;
    if (cosA >= 0) {
        sign = 1;
    } else {
        cosA = -cosA;
        sign = -1;
    }

    const f = chebyshevRatiosUsingCosAngle(t, cosA);
    return addQuaternion(mulQuaternion(q0, f[0]),
        mulQuaternion(q1, sign * f[1]));
}
