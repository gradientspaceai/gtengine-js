// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RigidBody.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The classes in this file support an implementation of collision response
// for acceleration-based constrained motion using impulse functions. The
// description is in Chapter 6 of "Game Physics, 2nd edition". For a
// description of the construction of impulse forces, see
// https://www.geometrictools.com/Documentation/ComputingImpulsiveForces.pdf
//
// Port notes:
// - 'std::function<Vector3<T>(T, RigidBodyState<T> const&)>' becomes the
//   exported type alias RigidBodyFunction, following OdeSolver.ts's
//   OdeFunction precedent. Upstream default-constructs the Force and Torque
//   functionals (calling one before it is set throws std::bad_function_call);
//   the port initializes them to null and asserts in update().
// - 'std::shared_ptr<RigidBody<T>>' becomes a plain reference; the members A
//   and B of RigidBodyContact default to null as upstream's default-
//   constructed shared pointers do, and applyImpulse() asserts they are set.
// - Upstream's accessors return 'const&' into the state. The port returns
//   copies (as Transform.ts does) so that callers cannot mutate the state
//   through a returned object; the setters likewise copy their arguments,
//   matching C++ value semantics.
// - Only the GTE_USE_MAT_VEC branches are ported, so dq/dt = 0.5 * w * q and
//   world vectors are M*v.
// - Upstream bug fixed: RigidBodyState::SetQOrientation computes
//   mROrientation from the *input* quaternion rather than from the member
//   mQOrientation that was just normalized. The Runge-Kutta integrator in
//   Update() always passes normalize = true with a non-unit quaternion
//   (q + h*dq/dt), so upstream builds the rotation matrix (and therefore the
//   world inertia tensors) from a non-unit quaternion, which scales them by
//   |q|^2 and |q|^4 respectively and corrupts the integration. The port uses
//   the normalized member. See the comment in setQOrientation.

import { logAssert } from './Logger.js';
import {
    Matrix, mulMatrix, multiplyABT
} from './Matrix.js';
import { inverse3x3 } from './Matrix3x3.js';
import { LinearSystem } from './LinearSystem.js';
import {
    Quaternion, addQuaternion, mulQuaternion
} from './Quaternion.js';
import { Rotation } from './Rotation.js';
import { Vector, add, dot, mul, normalize, sub } from './Vector.js';
import { cross } from './Vector3.js';

// The force and torque functionals. The first input is the simulation time.
// The second input is the rigid body state.
export type RigidBodyFunction =
    (t: number, state: RigidBodyState) => Vector;

// The rigid body state is stored in a separate class so that the force and
// torque functionals can be passed a single object to avoid a large number
// of parameters that would otherwise have to be passed to the functionals.
// This makes the Runge-Kutta ODE solver easier to read. The RigidBody class
// provides wrappers around the state accessors to avoid exposing a public
// state member.
export class RigidBodyState {
    // Constant quantities during the simulation.
    private mMass: number;
    private mInvMass: number;
    private mBodyInertia: Matrix;
    private mBodyInverseInertia: Matrix;

    // State variables in the differential equations of motion.
    private mPosition: Vector;
    private mQOrientation: Quaternion;
    private mLinearMomentum: Vector;
    private mAngularMomentum: Vector;

    // Quantities derived from the state variables.
    private mWorldInertia: Matrix;
    private mWorldInverseInertia: Matrix;
    private mROrientation: Matrix;
    private mLinearVelocity: Vector;
    private mAngularVelocity: Vector;
    private mQAngularVelocity: Quaternion;

    constructor() {
        this.mMass = 0;
        this.mInvMass = 0;
        this.mBodyInertia = Matrix.zero(3, 3);
        this.mBodyInverseInertia = Matrix.zero(3, 3);
        this.mPosition = Vector.zero(3);
        this.mQOrientation = Quaternion.identity();
        this.mLinearMomentum = Vector.zero(3);
        this.mAngularMomentum = Vector.zero(3);
        this.mWorldInertia = Matrix.zero(3, 3);
        this.mWorldInverseInertia = Matrix.zero(3, 3);
        this.mROrientation = Matrix.identity(3, 3);
        this.mLinearVelocity = Vector.zero(3);
        this.mAngularVelocity = Vector.zero(3);
        this.mQAngularVelocity = new Quaternion();
    }

    // Set the mass to a positive number for movable bodies. Set the mass to
    // zero for immovable objects. A body is immovable in the physics
    // simulation, but you can position and orient the immovable body
    // manually, typically during the creation of the physics objects.
    setMass(mass: number): void {
        if (mass > 0) {
            this.mMass = mass;
            this.mInvMass = 1 / mass;
        } else {
            this.mMass = 0;
            this.mInvMass = 0;
        }
    }

    // Set the body inertia to a positive definite matrix for movable bodies.
    // Set the body inertia to the zero matrix for immovable objects, but you
    // can position and orient the immovable body manually, typically during
    // the creation of the physics objects.
    setBodyInertia(bodyInertia: Matrix): void {
        logAssert(bodyInertia.numRows === 3 && bodyInertia.numCols === 3,
            'RigidBodyState: expecting a 3x3 body inertia.');
        const zero = Matrix.zero(3, 3);

        if (bodyInertia.notEquals(zero)) {
            this.mBodyInertia = bodyInertia.clone();
            this.mBodyInverseInertia = inverse3x3(bodyInertia).inverse;
            this.updateWorldInertialQuantities();
        } else {
            this.mBodyInertia = zero;
            this.mBodyInverseInertia = Matrix.zero(3, 3);
            this.mWorldInertia = Matrix.zero(3, 3);
            this.mWorldInverseInertia = Matrix.zero(3, 3);
        }
    }

    isMovable(): boolean {
        return this.mMass > 0;
    }

    isImmovable(): boolean {
        return this.mMass === 0;
    }

    setPosition(position: Vector): void {
        logAssert(position.size === 3, 'RigidBodyState: expecting a 3-tuple.');
        this.mPosition = Vector.fromArray(position.values);
    }

    setQOrientation(qOrientation: Quaternion,
        normalizeIt: boolean = false): void {
        this.mQOrientation = qOrientation.clone();
        if (normalizeIt) {
            normalize(this.mQOrientation);
        }

        // Upstream passes the *input* qOrientation here rather than the
        // normalized member, which produces a non-rotation matrix (and
        // therefore incorrect world inertia tensors) whenever the caller
        // passes a non-unit quaternion together with normalizeIt = true, as
        // RigidBody.update() does at every Runge-Kutta stage. The port uses
        // the normalized member.
        this.mROrientation = Rotation.fromQuaternion(this.mQOrientation, 3)
            .toMatrix();
        if (this.isMovable()) {
            this.updateWorldInertialQuantities();
        }
    }

    setLinearMomentum(linearMomentum: Vector): void {
        logAssert(linearMomentum.size === 3,
            'RigidBodyState: expecting a 3-tuple.');
        if (this.isMovable()) {
            this.mLinearMomentum = Vector.fromArray(linearMomentum.values);
            this.mLinearVelocity = mul(linearMomentum, this.mInvMass);
        }
    }

    setAngularMomentum(angularMomentum: Vector): void {
        logAssert(angularMomentum.size === 3,
            'RigidBodyState: expecting a 3-tuple.');
        if (this.isMovable()) {
            this.mAngularMomentum = Vector.fromArray(angularMomentum.values);
            this.mAngularVelocity =
                mulMatrix(this.mWorldInverseInertia, angularMomentum);
            this.copyAngularVelocityToQuaternion();
        }
    }

    setROrientation(rOrientation: Matrix): void {
        logAssert(rOrientation.numRows === 3 && rOrientation.numCols === 3,
            'RigidBodyState: expecting a 3x3 orientation.');
        this.mROrientation = rOrientation.clone();
        this.mQOrientation = Rotation.fromMatrix(rOrientation).toQuaternion();
        if (this.isMovable()) {
            this.updateWorldInertialQuantities();
        }
    }

    setLinearVelocity(linearVelocity: Vector): void {
        logAssert(linearVelocity.size === 3,
            'RigidBodyState: expecting a 3-tuple.');
        if (this.isMovable()) {
            this.mLinearVelocity = Vector.fromArray(linearVelocity.values);
            this.mLinearMomentum = mul(linearVelocity, this.mMass);
        }
    }

    setAngularVelocity(angularVelocity: Vector): void {
        logAssert(angularVelocity.size === 3,
            'RigidBodyState: expecting a 3-tuple.');
        if (this.isMovable()) {
            this.mAngularVelocity = Vector.fromArray(angularVelocity.values);
            this.mAngularMomentum =
                mulMatrix(this.mWorldInertia, angularVelocity);
            this.copyAngularVelocityToQuaternion();
        }
    }

    getMass(): number {
        return this.mMass;
    }

    getInverseMass(): number {
        return this.mInvMass;
    }

    getBodyInertia(): Matrix {
        return this.mBodyInertia.clone();
    }

    getBodyInverseInertia(): Matrix {
        return this.mBodyInverseInertia.clone();
    }

    getWorldInertia(): Matrix {
        return this.mWorldInertia.clone();
    }

    getWorldInverseInertia(): Matrix {
        return this.mWorldInverseInertia.clone();
    }

    getPosition(): Vector {
        return Vector.fromArray(this.mPosition.values);
    }

    getQOrientation(): Quaternion {
        return this.mQOrientation.clone();
    }

    getLinearMomentum(): Vector {
        return Vector.fromArray(this.mLinearMomentum.values);
    }

    getAngularMomentum(): Vector {
        return Vector.fromArray(this.mAngularMomentum.values);
    }

    getROrientation(): Matrix {
        return this.mROrientation.clone();
    }

    getLinearVelocity(): Vector {
        return Vector.fromArray(this.mLinearVelocity.values);
    }

    getAngularVelocity(): Vector {
        return Vector.fromArray(this.mAngularVelocity.values);
    }

    getQAngularVelocity(): Quaternion {
        return this.mQAngularVelocity.clone();
    }

    private copyAngularVelocityToQuaternion(): void {
        this.mQAngularVelocity = new Quaternion(
            this.mAngularVelocity.values[0],
            this.mAngularVelocity.values[1],
            this.mAngularVelocity.values[2],
            0);
    }

    private updateWorldInertialQuantities(): void {
        this.mWorldInertia = multiplyABT(
            mulMatrix(this.mROrientation, this.mBodyInertia),
            this.mROrientation);

        this.mWorldInverseInertia = multiplyABT(
            mulMatrix(this.mROrientation, this.mBodyInverseInertia),
            this.mROrientation);
    }
}

export class RigidBody {
    // Force and torque functions. The first input is the simulation time.
    // The second input is the rigid body state. These functions must be set
    // before starting the simulation.
    force: RigidBodyFunction | null;
    torque: RigidBodyFunction | null;

    private mState: RigidBodyState;

    // The rigid body state is initialized to zero values. Set the members
    // before starting the simulation. For immovable objects, set mass to
    // zero.
    constructor() {
        this.force = null;
        this.torque = null;
        this.mState = new RigidBodyState();
    }

    // Set the mass to a positive number for movable bodies. Set the mass to
    // zero for immovable objects.
    setMass(mass: number): void {
        this.mState.setMass(mass);
    }

    // Set the body inertia to a positive definite matrix for movable bodies.
    // Set the body inertia to the zero matrix for immovable objects.
    setBodyInertia(bodyInertia: Matrix): void {
        this.mState.setBodyInertia(bodyInertia);
    }

    isMovable(): boolean {
        return this.mState.isMovable();
    }

    isImmovable(): boolean {
        return this.mState.isImmovable();
    }

    setPosition(position: Vector): void {
        this.mState.setPosition(position);
    }

    setQOrientation(qOrientation: Quaternion,
        normalizeIt: boolean = false): void {
        this.mState.setQOrientation(qOrientation, normalizeIt);
    }

    setLinearMomentum(linearMomentum: Vector): void {
        this.mState.setLinearMomentum(linearMomentum);
    }

    setAngularMomentum(angularMomentum: Vector): void {
        this.mState.setAngularMomentum(angularMomentum);
    }

    setROrientation(rOrientation: Matrix): void {
        this.mState.setROrientation(rOrientation);
    }

    setLinearVelocity(linearVelocity: Vector): void {
        this.mState.setLinearVelocity(linearVelocity);
    }

    setAngularVelocity(angularVelocity: Vector): void {
        this.mState.setAngularVelocity(angularVelocity);
    }

    getMass(): number {
        return this.mState.getMass();
    }

    getInverseMass(): number {
        return this.mState.getInverseMass();
    }

    getBodyInertia(): Matrix {
        return this.mState.getBodyInertia();
    }

    getBodyInverseInertia(): Matrix {
        return this.mState.getBodyInverseInertia();
    }

    getWorldInertia(): Matrix {
        return this.mState.getWorldInertia();
    }

    getWorldInverseInertia(): Matrix {
        return this.mState.getWorldInverseInertia();
    }

    getPosition(): Vector {
        return this.mState.getPosition();
    }

    getQOrientation(): Quaternion {
        return this.mState.getQOrientation();
    }

    getLinearMomentum(): Vector {
        return this.mState.getLinearMomentum();
    }

    getAngularMomentum(): Vector {
        return this.mState.getAngularMomentum();
    }

    getROrientation(): Matrix {
        return this.mState.getROrientation();
    }

    getLinearVelocity(): Vector {
        return this.mState.getLinearVelocity();
    }

    getAngularVelocity(): Vector {
        return this.mState.getAngularVelocity();
    }

    getQAngularVelocity(): Quaternion {
        return this.mState.getQAngularVelocity();
    }

    // Runge-Kutta fourth-order differential equation solver.
    update(t: number, dt: number): void {
        const force = this.force;
        const torque = this.torque;
        logAssert(force !== null && torque !== null,
            'RigidBody: the force and torque functions must be set.');

        const halfDT = 0.5 * dt;
        const sixthDT = dt / 6;
        const tpHalfDT = t + halfDT;
        const tpDT = t + dt;

        const newState = new RigidBodyState();
        newState.setMass(this.getMass());
        newState.setBodyInertia(this.getBodyInertia());

        // A1 = G(T,S0), B1 = S0 + (DT/2)*A1
        const a1DXDT = this.getLinearVelocity();
        let w = this.getQAngularVelocity();
        const a1DQDT = mulQuaternion(mulQuaternion(0.5, w),
            this.getQOrientation());
        const a1DPDT = force(t, this.mState);
        const a1DLDT = torque(t, this.mState);
        newState.setPosition(add(this.getPosition(), mul(a1DXDT, halfDT)));
        newState.setQOrientation(addQuaternion(this.getQOrientation(),
            mulQuaternion(a1DQDT, halfDT)), true);
        newState.setLinearMomentum(add(this.getLinearMomentum(),
            mul(a1DPDT, halfDT)));
        newState.setAngularMomentum(add(this.getAngularMomentum(),
            mul(a1DLDT, halfDT)));

        // A2 = G(T+DT/2,B1), B2 = S0 + (DT/2)*A2
        const a2DXDT = newState.getLinearVelocity();
        w = newState.getQAngularVelocity();
        const a2DQDT = mulQuaternion(mulQuaternion(0.5, w),
            newState.getQOrientation());
        const a2DPDT = force(tpHalfDT, newState);
        const a2DLDT = torque(tpHalfDT, newState);
        newState.setPosition(add(this.getPosition(), mul(a2DXDT, halfDT)));
        newState.setQOrientation(addQuaternion(this.getQOrientation(),
            mulQuaternion(a2DQDT, halfDT)), true);
        newState.setLinearMomentum(add(this.getLinearMomentum(),
            mul(a2DPDT, halfDT)));
        newState.setAngularMomentum(add(this.getAngularMomentum(),
            mul(a2DLDT, halfDT)));

        // A3 = G(T+DT/2,B2), B3 = S0 + DT*A3
        const a3DXDT = newState.getLinearVelocity();
        w = newState.getQAngularVelocity();
        const a3DQDT = mulQuaternion(mulQuaternion(0.5, w),
            newState.getQOrientation());
        const a3DPDT = force(tpHalfDT, newState);
        const a3DLDT = torque(tpHalfDT, newState);
        newState.setPosition(add(this.getPosition(), mul(a3DXDT, dt)));
        newState.setQOrientation(addQuaternion(this.getQOrientation(),
            mulQuaternion(a3DQDT, dt)), true);
        newState.setLinearMomentum(add(this.getLinearMomentum(),
            mul(a3DPDT, dt)));
        newState.setAngularMomentum(add(this.getAngularMomentum(),
            mul(a3DLDT, dt)));

        // A4 = G(T+DT,B3), S1 = S0 + (DT/6)*(A1+2*(A2+A3)+A4)
        const a4DXDT = newState.getLinearVelocity();
        w = newState.getQAngularVelocity();
        const a4DQDT = mulQuaternion(mulQuaternion(0.5, w),
            newState.getQOrientation());
        const a4DPDT = force(tpDT, newState);
        const a4DLDT = torque(tpDT, newState);

        this.setPosition(add(this.getPosition(), mul(
            rk4Combine(a1DXDT, a2DXDT, a3DXDT, a4DXDT), sixthDT)));

        this.setQOrientation(addQuaternion(this.getQOrientation(),
            mulQuaternion(rk4CombineQuaternion(a1DQDT, a2DQDT, a3DQDT,
                a4DQDT), sixthDT)), true);

        this.setLinearMomentum(add(this.getLinearMomentum(), mul(
            rk4Combine(a1DPDT, a2DPDT, a3DPDT, a4DPDT), sixthDT)));

        this.setAngularMomentum(add(this.getAngularMomentum(), mul(
            rk4Combine(a1DLDT, a2DLDT, a3DLDT, a4DLDT), sixthDT)));
    }
}

// The weighted sum a1 + 2*(a2 + a3) + a4 of the Runge-Kutta stages.
function rk4Combine(a1: Vector, a2: Vector, a3: Vector, a4: Vector): Vector {
    return add(add(a1, mul(add(a2, a3), 2)), a4);
}

function rk4CombineQuaternion(a1: Quaternion, a2: Quaternion, a3: Quaternion,
    a4: Quaternion): Quaternion {
    return addQuaternion(
        addQuaternion(a1, mulQuaternion(addQuaternion(a2, a3), 2)), a4);
}

// The rigid body contact stores basic information. The class can be extended
// by derivation to allow for additional information that is specific to a
// simulation.
export class RigidBodyContact {
    // Body A has the vertex in a vertex-face contact or edge-edge contact.
    A: RigidBody | null;

    // Body B has the face in a vertex-face contact, and the normal N is for
    // that face. If there is instead an edge-edge contact, the normal N is
    // the cross product of the edges.
    B: RigidBody | null;

    // The intersection point at contact.
    P: Vector;

    // The outward unit-length normal to the face at the contact point.
    N: Vector;

    // The coefficient of restitution which is in [0,1]. This allows for the
    // loss of kinetic energy at a contact point.
    restitution: number;

    constructor() {
        this.A = null;
        this.B = null;
        this.P = Vector.zero(3);
        this.N = Vector.zero(3);
        this.restitution = 0;
    }

    // Call applyImpulse for each rigid body in your simulation at the
    // current simulation time. After all such calls are made, then iterate
    // over all your rigid bodies and call their update(time, deltaTime)
    // functions.
    applyImpulse(): void {
        const A = this.A;
        const B = this.B;
        logAssert(A !== null && B !== null,
            'RigidBodyContact: the bodies A and B must be set.');

        // The positions of the centers of mass.
        const XA = A.getPosition();
        const XB = B.getPosition();

        // The location of the contact points relative to the centers of
        // mass.
        const rA = sub(this.P, XA);
        const rB = sub(this.P, XB);

        // The preimpulse linear velocities of the centers of mass.
        const linvelANeg = A.getLinearVelocity();
        const linvelBNeg = B.getLinearVelocity();

        // The preimpulse angular velocities about the centers of mass.
        const angvelANeg = A.getAngularVelocity();
        const angvelBNeg = B.getAngularVelocity();

        // The preimpulse velocities of P0.
        const velANeg = add(linvelANeg, cross(angvelANeg, rA));
        const velBNeg = add(linvelBNeg, cross(angvelBNeg, rB));
        const velDiffNeg = sub(velANeg, velBNeg);

        // The preimpulse linear momenta of the centers of mass.
        const linmomANeg = A.getLinearMomentum();
        const linmomBNeg = B.getLinearMomentum();

        // The preimpulse angular momenta about the centers of mass.
        const angmomANeg = A.getAngularMomentum();
        const angmomBNeg = B.getAngularMomentum();

        // The inverse masses, inverse world inertia tensors and quadratic
        // forms associated with these tensors.
        const sumInvMasses = A.getInverseMass() + B.getInverseMass();
        const invJA = A.getWorldInverseInertia();
        const invJB = B.getWorldInverseInertia();

        const T0 = sub(velDiffNeg,
            mul(this.N, dot(this.N, velDiffNeg)));
        // normalize() zeroes T0 and returns 0 when T0 has zero length, so a
        // positive return value is upstream's 'T0 != Zero' test.
        const t0Length = normalize(T0);
        if (t0Length > 0) {
            // T0 is tangent at P, unit length and perpendicular to N.
            const T1 = cross(this.N, T0);
            const rAxN = cross(rA, this.N);
            const rAxT0 = cross(rA, T0);
            const rAxT1 = cross(rA, T1);
            const rBxN = cross(rB, this.N);
            const rBxT0 = cross(rB, T0);
            const rBxT1 = cross(rB, T1);

            // The matrix constructed here is positive definite. This ensures
            // the linear system always has a solution, so the invertible
            // flag from LinearSystem.solve3x3 is ignored.
            const sysMatrix = new Matrix(3, 3);
            sysMatrix.set(0, 0, sumInvMasses
                + dot(rAxN, mulMatrix(invJA, rAxN))
                + dot(rBxN, mulMatrix(invJB, rBxN)));
            sysMatrix.set(1, 1, sumInvMasses
                + dot(rAxT0, mulMatrix(invJA, rAxT0))
                + dot(rBxT0, mulMatrix(invJB, rBxT0)));
            sysMatrix.set(2, 2, sumInvMasses
                + dot(rAxT1, mulMatrix(invJA, rAxT1))
                + dot(rBxT1, mulMatrix(invJB, rBxT1)));
            sysMatrix.set(0, 1, dot(rAxN, mulMatrix(invJA, rAxT0))
                + dot(rBxN, mulMatrix(invJB, rBxT0)));
            sysMatrix.set(0, 2, dot(rAxN, mulMatrix(invJA, rAxT1))
                + dot(rBxN, mulMatrix(invJB, rBxT1)));
            sysMatrix.set(1, 2, dot(rAxT0, mulMatrix(invJA, rAxT1))
                + dot(rBxT0, mulMatrix(invJB, rBxT1)));
            sysMatrix.set(1, 0, sysMatrix.get(0, 1));
            sysMatrix.set(2, 0, sysMatrix.get(0, 2));
            sysMatrix.set(2, 1, sysMatrix.get(1, 2));
            const sysInput = Vector.fromArray([
                -(1 + this.restitution) * dot(this.N, velDiffNeg), 0, 0
            ]);
            const sysOutput = LinearSystem.solve3x3(sysMatrix, sysInput).X;

            // Apply the impulsive force to the bodies to change linear and
            // angular momentum.
            const impulse = add(add(
                mul(this.N, sysOutput.values[0]),
                mul(T0, sysOutput.values[1])),
                mul(T1, sysOutput.values[2]));
            A.setLinearMomentum(add(linmomANeg, impulse));
            B.setLinearMomentum(sub(linmomBNeg, impulse));
            A.setAngularMomentum(add(angmomANeg, cross(rA, impulse)));
            B.setAngularMomentum(sub(angmomBNeg, cross(rB, impulse)));
        } else {
            // Fall back to the impulse force f*N0 when the relative velocity
            // at the contact P0 is parallel to N0.
            const rAxN = cross(rA, this.N);
            const rBxN = cross(rB, this.N);
            const quadformA = dot(rAxN, mulMatrix(invJA, rAxN));
            const quadformB = dot(rBxN, mulMatrix(invJB, rBxN));

            // The magnitude of the impulse force.
            const numer = -(1 + this.restitution) * dot(this.N, velDiffNeg);
            const denom = sumInvMasses + quadformA + quadformB;
            const f = numer / denom;

            // Apply the impulsive force to the bodies to change linear and
            // angular momentum.
            const impulse = mul(this.N, f);
            A.setLinearMomentum(add(linmomANeg, impulse));
            B.setLinearMomentum(sub(linmomBNeg, impulse));
            A.setAngularMomentum(add(angmomANeg, cross(rA, impulse)));
            B.setAngularMomentum(sub(angmomBNeg, cross(rB, impulse)));
        }
    }
}
