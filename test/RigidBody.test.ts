import { describe, expect, it } from 'vitest';
import {
    RigidBody, RigidBodyContact, RigidBodyState
} from '../src/RigidBody.js';
import { AxisAngle } from '../src/AxisAngle.js';
import {
    Matrix, lInfinityNorm, mulMatrix, multiplyAB, multiplyABT, subMatrix
} from '../src/Matrix.js';
import { Quaternion } from '../src/Quaternion.js';
import { Rotation } from '../src/Rotation.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { check, expectClose, fc, scaled } from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function diag3(a: number, b: number, c: number): Matrix {
    return Matrix.fromArray(3, 3, [a, 0, 0, 0, b, 0, 0, 0, c]);
}

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

function len(v: Vector): number {
    return Math.sqrt(dot(v, v));
}

const zeroForce = () => Vector.zero(3);

// A body with the given mass and body inertia, at the origin with the
// identity orientation and no motion, and with zero force and torque.
function makeFreeBody(mass: number, inertia: Matrix): RigidBody {
    const body = new RigidBody();
    body.setMass(mass);
    body.setBodyInertia(inertia);
    body.setPosition(Vector.zero(3));
    body.setQOrientation(Quaternion.identity());
    body.setLinearMomentum(Vector.zero(3));
    body.setAngularMomentum(Vector.zero(3));
    body.force = zeroForce;
    body.torque = zeroForce;
    return body;
}

describe('RigidBodyState', () => {
    it('handles movable and immovable bodies', () => {
        const state = new RigidBodyState();
        expect(state.isImmovable()).toBe(true);
        expect(state.isMovable()).toBe(false);
        expect(state.getMass()).toBe(0);
        expect(state.getInverseMass()).toBe(0);

        // An immovable body ignores the momentum and velocity setters.
        state.setLinearMomentum(v3(1, 2, 3));
        state.setAngularMomentum(v3(1, 2, 3));
        state.setLinearVelocity(v3(1, 2, 3));
        state.setAngularVelocity(v3(1, 2, 3));
        expect(state.getLinearMomentum().values).toEqual([0, 0, 0]);
        expect(state.getAngularMomentum().values).toEqual([0, 0, 0]);
        expect(state.getLinearVelocity().values).toEqual([0, 0, 0]);
        expect(state.getAngularVelocity().values).toEqual([0, 0, 0]);

        // A position may still be set for an immovable body.
        state.setPosition(v3(4, 5, 6));
        expect(state.getPosition().values).toEqual([4, 5, 6]);

        // A negative mass is clamped to the immovable case.
        state.setMass(-3);
        expect(state.isImmovable()).toBe(true);

        state.setMass(4);
        expect(state.isMovable()).toBe(true);
        expect(state.getInverseMass()).toBe(0.25);
        state.setLinearMomentum(v3(8, 0, -4));
        expect(state.getLinearVelocity().values).toEqual([2, 0, -1]);
        state.setLinearVelocity(v3(1, 2, 3));
        expect(state.getLinearMomentum().values).toEqual([4, 8, 12]);
    });

    it('inverts the body inertia and zeroes it for the zero matrix', () => {
        const state = new RigidBodyState();
        state.setMass(1);
        state.setBodyInertia(diag3(2, 4, 8));
        expect(maxDiff(state.getBodyInverseInertia(), diag3(0.5, 0.25, 0.125)))
            .toBeLessThan(1e-15);

        // The orientation is the identity, so the world tensors equal the
        // body tensors.
        expect(maxDiff(state.getWorldInertia(), diag3(2, 4, 8)))
            .toBeLessThan(1e-15);
        expect(maxDiff(state.getWorldInverseInertia(),
            diag3(0.5, 0.25, 0.125))).toBeLessThan(1e-15);

        state.setBodyInertia(Matrix.zero(3, 3));
        expect(maxDiff(state.getBodyInertia(), Matrix.zero(3, 3))).toBe(0);
        expect(maxDiff(state.getBodyInverseInertia(), Matrix.zero(3, 3)))
            .toBe(0);
        expect(maxDiff(state.getWorldInertia(), Matrix.zero(3, 3))).toBe(0);
        expect(maxDiff(state.getWorldInverseInertia(), Matrix.zero(3, 3)))
            .toBe(0);
    });

    it('keeps the world inertia equal to R*J*R^T', () => {
        const state = new RigidBodyState();
        state.setMass(3);
        const J = diag3(2, 5, 7);
        state.setBodyInertia(J);

        const axis = v3(1, -2, 3);
        normalize(axis);
        const R = Rotation.fromAxisAngle(new AxisAngle(axis, 0.9)).toMatrix();
        state.setROrientation(R);

        expect(maxDiff(state.getWorldInertia(),
            multiplyABT(multiplyAB(R, J), R))).toBeLessThan(1e-14);
        expect(maxDiff(multiplyAB(state.getWorldInertia(),
            state.getWorldInverseInertia()), Matrix.identity(3, 3)))
            .toBeLessThan(1e-14);

        // The quaternion and matrix orientations agree.
        const q = state.getQOrientation();
        expect(len(q)).toBeCloseTo(1, 14);
        expect(maxDiff(Rotation.fromQuaternion(q, 3).toMatrix(), R))
            .toBeLessThan(1e-14);
    });

    it('builds the rotation matrix from the normalized quaternion', () => {
        // Upstream builds mROrientation from the un-normalized input
        // quaternion, which yields a matrix scaled by |q|^2. The port uses
        // the normalized member, so the matrix is always a rotation.
        const state = new RigidBodyState();
        state.setMass(1);
        state.setBodyInertia(diag3(1, 1, 1));

        const scale = 3;
        const unit = new Quaternion(0.5, 0.5, 0.5, 0.5);
        const scaled = new Quaternion(0.5 * scale, 0.5 * scale, 0.5 * scale,
            0.5 * scale);
        state.setQOrientation(scaled, true);

        const R = state.getROrientation();
        expect(maxDiff(multiplyABT(R, R), Matrix.identity(3, 3)))
            .toBeLessThan(1e-14);
        expect(maxDiff(R, Rotation.fromQuaternion(unit, 3).toMatrix()))
            .toBeLessThan(1e-14);
        expect(len(state.getQOrientation())).toBeCloseTo(1, 15);

        // The world inertia of an isotropic body is unchanged by any
        // rotation; upstream's |q|^4 scaling would break this.
        expect(maxDiff(state.getWorldInertia(), diag3(1, 1, 1)))
            .toBeLessThan(1e-14);
    });

    it('copies the angular velocity into the quaternion form', () => {
        const state = new RigidBodyState();
        state.setMass(1);
        state.setBodyInertia(diag3(2, 2, 2));
        state.setAngularVelocity(v3(1, -2, 3));
        expect(state.getQAngularVelocity().values).toEqual([1, -2, 3, 0]);
        expect(state.getAngularMomentum().values).toEqual([2, -4, 6]);

        state.setAngularMomentum(v3(2, -4, 6));
        expect(state.getAngularVelocity().values).toEqual([1, -2, 3]);
        expect(state.getQAngularVelocity().values).toEqual([1, -2, 3, 0]);
    });

    it('returns copies so callers cannot mutate the state', () => {
        const state = new RigidBodyState();
        state.setMass(1);
        state.setPosition(v3(1, 2, 3));
        const p = state.getPosition();
        p.values[0] = 99;
        expect(state.getPosition().values[0]).toBe(1);

        const input = v3(4, 5, 6);
        state.setPosition(input);
        input.values[1] = 99;
        expect(state.getPosition().values[1]).toBe(5);
    });
});

describe('RigidBody integration', () => {
    it('requires the force and torque functions', () => {
        const body = new RigidBody();
        body.setMass(1);
        body.setBodyInertia(diag3(1, 1, 1));
        expect(() => body.update(0, 0.01)).toThrow();
    });

    it('reproduces the closed form for a constant force', () => {
        // Runge-Kutta 4 integrates a quadratic trajectory exactly (up to
        // round-off), so the numerical solution matches
        //   x(t) = x0 + v0*t + (1/2)*(F/m)*t^2
        //   p(t) = p0 + F*t
        const mass = 2.5;
        const F = v3(0, -9.81 * mass, 1.5);
        const body = makeFreeBody(mass, diag3(1, 1, 1));
        const x0 = v3(1, 20, -3);
        const v0 = v3(2, 5, -1);
        body.setPosition(x0);
        body.setLinearVelocity(v0);
        body.force = () => Vector.fromArray(F.values);

        const dt = 0.01;
        const numSteps = 400;
        for (let i = 0; i < numSteps; ++i) {
            body.update(i * dt, dt);
        }

        const t = numSteps * dt;
        const a = mul(F, 1 / mass);
        const expected = add(add(x0, mul(v0, t)), mul(a, 0.5 * t * t));
        expect(vecMaxDiff(body.getPosition(), expected)).toBeLessThan(1e-10);
        expect(vecMaxDiff(body.getLinearMomentum(),
            add(mul(v0, mass), mul(F, t)))).toBeLessThan(1e-10);

        // With no torque the orientation never changes.
        expect(maxDiff(body.getROrientation(), Matrix.identity(3, 3)))
            .toBeLessThan(1e-14);
    });

    it('rotates an isotropic body at a constant rate about a fixed axis',
        () => {
            // For an isotropic inertia tensor the world inverse inertia is
            // constant, so the angular velocity is constant and the
            // orientation is the rotation about w by |w|*t.
            const inertia = 3;
            const body = makeFreeBody(1, diag3(inertia, inertia, inertia));
            const w = v3(0.7, -1.3, 0.4);
            body.setAngularVelocity(w);

            const dt = 0.002;
            const numSteps = 1000;
            for (let i = 0; i < numSteps; ++i) {
                body.update(i * dt, dt);
            }

            const t = numSteps * dt;
            const speed = len(w);
            const axis = mul(w, 1 / speed);
            const expected = Rotation.fromAxisAngle(
                new AxisAngle(axis, speed * t)).toMatrix();
            expect(maxDiff(body.getROrientation(), expected))
                .toBeLessThan(1e-10);

            // The angular velocity and momentum are unchanged.
            expect(vecMaxDiff(body.getAngularVelocity(), w))
                .toBeLessThan(1e-12);
            expect(vecMaxDiff(body.getAngularMomentum(),
                mul(w, inertia))).toBeLessThan(1e-12);
            expect(len(body.getQOrientation())).toBeCloseTo(1, 12);
        });

    it('conserves angular momentum and kinetic energy for a free '
        + 'asymmetric top', () => {
            const body = makeFreeBody(1, diag3(1, 2, 3));
            const w0 = v3(1.2, -0.8, 0.5);
            body.setAngularVelocity(w0);

            const L0 = body.getAngularMomentum();
            const E0 = 0.5 * dot(body.getAngularVelocity(), L0);
            expect(E0).toBeGreaterThan(0);

            const dt = 0.005;
            const numSteps = 1000;
            let worstL = 0;
            let worstE = 0;
            let worstNorm = 0;
            let worstOrtho = 0;
            for (let i = 0; i < numSteps; ++i) {
                body.update(i * dt, dt);
                worstL = Math.max(worstL,
                    vecMaxDiff(body.getAngularMomentum(), L0));
                const E = 0.5 * dot(body.getAngularVelocity(),
                    body.getAngularMomentum());
                worstE = Math.max(worstE, Math.abs(E - E0) / E0);
                worstNorm = Math.max(worstNorm,
                    Math.abs(len(body.getQOrientation()) - 1));
                const R = body.getROrientation();
                worstOrtho = Math.max(worstOrtho,
                    maxDiff(multiplyABT(R, R), Matrix.identity(3, 3)));
            }

            // The torque is zero, so the angular momentum is conserved by
            // construction. The kinetic energy is conserved only if the
            // orientation (and hence the world inertia tensor) is
            // integrated correctly.
            expect(worstL).toBeLessThan(1e-12);
            expect(worstE).toBeLessThan(1e-8);
            expect(worstNorm).toBeLessThan(1e-12);
            expect(worstOrtho).toBeLessThan(1e-12);

            // The body-frame angular velocity traces the polhode: the two
            // Euler invariants 2*E and |L|^2 are unchanged in body
            // coordinates as well.
            const R = body.getROrientation();
            const bodyW = mulMatrix(transpose3(R), body.getAngularVelocity());
            const bodyL = mulMatrix(transpose3(R), body.getAngularMomentum());
            expect(0.5 * dot(bodyW, bodyL)).toBeCloseTo(E0, 6);
            expect(len(bodyL)).toBeCloseTo(len(L0), 8);
        });

    it('precesses a torque-free symmetric top with the analytic period',
        () => {
            // For a free symmetric top with body inertia diag(I1,I1,I3), the
            // symmetry axis precesses uniformly about the (fixed) angular
            // momentum with rate |L|/I1, so it returns to its initial world
            // direction after T = 2*pi*I1/|L|.
            const I1 = 1;
            const I3 = 2;
            const body = makeFreeBody(1, diag3(I1, I1, I3));
            const w0 = v3(0.6, 0, 0.8);
            body.setAngularVelocity(w0);

            const L = body.getAngularMomentum();
            expect(vecMaxDiff(L, v3(0.6, 0, 1.6))).toBeLessThan(1e-15);
            const period = 2 * Math.PI * I1 / len(L);

            const e3 = v3(0, 0, 1);
            const axis0 = mulMatrix(body.getROrientation(), e3);
            const cos0 = dot(axis0, L) / (len(axis0) * len(L));

            const numSteps = 2000;
            const dt = period / numSteps;
            let worstCos = 0;
            let worstSpeed = 0;
            const speed0 = len(body.getAngularVelocity());
            for (let i = 0; i < numSteps; ++i) {
                body.update(i * dt, dt);
                const axis = mulMatrix(body.getROrientation(), e3);
                worstCos = Math.max(worstCos, Math.abs(
                    dot(axis, L) / (len(axis) * len(L)) - cos0));
                worstSpeed = Math.max(worstSpeed,
                    Math.abs(len(body.getAngularVelocity()) - speed0));
            }

            // The angle between the symmetry axis and L is constant, as is
            // the magnitude of the angular velocity.
            expect(worstCos).toBeLessThan(1e-8);
            expect(worstSpeed).toBeLessThan(1e-8);

            // After one precession period the symmetry axis is back where it
            // started.
            const axisT = mulMatrix(body.getROrientation(), e3);
            expect(vecMaxDiff(axisT, axis0)).toBeLessThan(1e-6);

            // The angular velocity has also returned (it precesses about L
            // at the same rate).
            expect(vecMaxDiff(body.getAngularVelocity(), w0))
                .toBeLessThan(1e-6);
        });
});

describe('RigidBodyContact', () => {
    it('requires both bodies', () => {
        const contact = new RigidBodyContact();
        expect(() => contact.applyImpulse()).toThrow();
    });

    it('exchanges the velocities of equal masses in a head-on collision',
        () => {
            // The relative velocity at the contact is parallel to N, so the
            // fallback branch (a single impulse magnitude along N) is used.
            const A = makeFreeBody(1, diag3(1, 1, 1));
            const B = makeFreeBody(1, diag3(1, 1, 1));
            A.setPosition(v3(-1, 0, 0));
            B.setPosition(v3(1, 0, 0));
            A.setLinearVelocity(v3(1, 0, 0));
            B.setLinearVelocity(v3(0, 0, 0));

            const contact = new RigidBodyContact();
            contact.A = A;
            contact.B = B;
            contact.P = v3(0, 0, 0);
            contact.N = v3(-1, 0, 0);
            contact.restitution = 1;
            contact.applyImpulse();

            expect(vecMaxDiff(A.getLinearVelocity(), v3(0, 0, 0)))
                .toBeLessThan(1e-15);
            expect(vecMaxDiff(B.getLinearVelocity(), v3(1, 0, 0)))
                .toBeLessThan(1e-15);
            expect(vecMaxDiff(A.getAngularVelocity(), v3(0, 0, 0)))
                .toBeLessThan(1e-15);
            expect(vecMaxDiff(B.getAngularVelocity(), v3(0, 0, 0)))
                .toBeLessThan(1e-15);

            // A perfectly inelastic contact removes the normal relative
            // velocity instead.
            const A2 = makeFreeBody(1, diag3(1, 1, 1));
            const B2 = makeFreeBody(1, diag3(1, 1, 1));
            A2.setPosition(v3(-1, 0, 0));
            B2.setPosition(v3(1, 0, 0));
            A2.setLinearVelocity(v3(1, 0, 0));
            const contact2 = new RigidBodyContact();
            contact2.A = A2;
            contact2.B = B2;
            contact2.P = v3(0, 0, 0);
            contact2.N = v3(-1, 0, 0);
            contact2.restitution = 0;
            contact2.applyImpulse();
            expect(vecMaxDiff(A2.getLinearVelocity(), v3(0.5, 0, 0)))
                .toBeLessThan(1e-15);
            expect(vecMaxDiff(B2.getLinearVelocity(), v3(0.5, 0, 0)))
                .toBeLessThan(1e-15);
        });

    it('applies the three-component impulse with the expected relative '
        + 'velocity change', () => {
            for (const restitution of [0, 0.4, 1]) {
                const A = makeFreeBody(2, diag3(2, 3, 4));
                const B = makeFreeBody(5, diag3(4, 4, 6));
                A.setPosition(v3(-1, 0, 0));
                B.setPosition(v3(1, 0.2, -0.1));
                A.setLinearVelocity(v3(1.3, 0.5, -0.2));
                B.setLinearVelocity(v3(-0.4, 0.1, 0.3));
                A.setAngularVelocity(v3(0.2, -0.6, 0.9));
                B.setAngularVelocity(v3(-0.5, 0.3, 0.1));

                const P = v3(0, 0.3, 0.15);
                const N = v3(-1, 0, 0);

                const contact = new RigidBodyContact();
                contact.A = A;
                contact.B = B;
                contact.P = P;
                contact.N = N;
                contact.restitution = restitution;

                const rA = sub(P, A.getPosition());
                const rB = sub(P, B.getPosition());
                const contactVelocity = (body: RigidBody, r: Vector) =>
                    add(body.getLinearVelocity(),
                        cross(body.getAngularVelocity(), r));
                const before = sub(contactVelocity(A, rA),
                    contactVelocity(B, rB));

                // The tangent basis used by the implementation.
                const T0 = sub(before, mul(N, dot(N, before)));
                normalize(T0);
                const T1 = cross(N, T0);

                const totalLinearBefore = add(A.getLinearMomentum(),
                    B.getLinearMomentum());
                // The angular momentum about the contact point P is
                // L_i + (X_i - P) x p_i = L_i - r_i x p_i.
                const angularAboutP = () => add(
                    sub(A.getAngularMomentum(),
                        cross(rA, A.getLinearMomentum())),
                    sub(B.getAngularMomentum(),
                        cross(rB, B.getLinearMomentum())));
                const totalAngularBefore = angularAboutP();

                contact.applyImpulse();

                const after = sub(contactVelocity(A, rA),
                    contactVelocity(B, rB));

                // The impulse reverses the normal relative velocity scaled
                // by the coefficient of restitution and leaves the tangential
                // relative velocity unchanged.
                expect(dot(N, after)).toBeCloseTo(
                    -restitution * dot(N, before), 10);
                expect(dot(T0, after)).toBeCloseTo(dot(T0, before), 10);
                expect(dot(T1, after)).toBeCloseTo(dot(T1, before), 10);

                // The impulse is internal, so the total linear momentum and
                // the total angular momentum about the contact point are
                // conserved.
                expect(vecMaxDiff(add(A.getLinearMomentum(),
                    B.getLinearMomentum()), totalLinearBefore))
                    .toBeLessThan(1e-12);
                expect(vecMaxDiff(angularAboutP(), totalAngularBefore))
                    .toBeLessThan(1e-12);
            }
        });
});

function transpose3(M: Matrix): Matrix {
    const R = new Matrix(3, 3);
    for (let r = 0; r < 3; ++r) {
        for (let c = 0; c < 3; ++c) {
            R.set(r, c, M.get(c, r));
        }
    }
    return R;
}

// ---------------------------------------------------------------------------
// Verification (V41): properties of the state bookkeeping, the Runge-Kutta
// integrator (including the unit-quaternion regression of upstream issue
// #313) and the contact impulse.
// ---------------------------------------------------------------------------

const rbCoord = (max = 3) => scaled(-max, max, 512);

const unitQuaternion = (): fc.Arbitrary<Quaternion> =>
    fc.array(scaled(-1, 1, 512), { minLength: 4, maxLength: 4 })
        .filter(a => a.reduce((s, x) => s + x * x, 0) > 0.05)
        .map(a => {
            const q = new Quaternion(a[0], a[1], a[2], a[3]);
            normalize(q);
            return q;
        });

const unitVector3 = (): fc.Arbitrary<Vector> =>
    fc.array(scaled(-1, 1, 512), { minLength: 3, maxLength: 3 })
        .filter(a => a[0] ** 2 + a[1] ** 2 + a[2] ** 2 > 0.05)
        .map(a => {
            const v = Vector.fromArray(a);
            normalize(v);
            return v;
        });

const positiveInertia = (): fc.Arbitrary<Matrix> =>
    fc.tuple(scaled(0.3, 4, 128), scaled(0.3, 4, 128), scaled(0.3, 4, 128))
        .map(([a, b, c]) => diag3(a, b, c));

const vector3 = (max = 3) =>
    fc.array(rbCoord(max), { minLength: 3, maxLength: 3 })
        .map(a => Vector.fromArray(a));

// A movable body with a random state and no force or torque.
const freeBody = (): fc.Arbitrary<RigidBody> =>
    fc.tuple(scaled(0.4, 4, 128), positiveInertia(), unitQuaternion(),
        vector3(), vector3(), vector3())
        .map(([mass, inertia, q, position, linear, angular]) => {
            const body = new RigidBody();
            body.setMass(mass);
            body.setBodyInertia(inertia);
            body.setPosition(position);
            body.setQOrientation(q);
            body.setLinearMomentum(linear);
            body.setAngularMomentum(angular);
            body.force = zeroForce;
            body.torque = zeroForce;
            return body;
        });

function kineticEnergy(body: RigidBody): number {
    const v = body.getLinearVelocity();
    const w = body.getAngularVelocity();
    return 0.5 * body.getMass() * dot(v, v)
        + 0.5 * dot(w, body.getAngularMomentum());
}

describe('RigidBody verification', () => {
    it('keeps the world inertia consistent with R*J*R^T and its inverse',
        () => {
            check(fc.tuple(scaled(0.4, 4, 128), positiveInertia(),
                unitQuaternion()), ([mass, inertia, q]) => {
                    const state = new RigidBodyState();
                    state.setMass(mass);
                    state.setBodyInertia(inertia);
                    state.setQOrientation(q);
                    const R = state.getROrientation();
                    // R is a rotation because the port builds it from the
                    // normalized quaternion (upstream #313 uses the raw
                    // input, which is only a rotation when the caller passes
                    // a unit quaternion).
                    expect(maxDiff(multiplyABT(R, R), Matrix.identity(3, 3)))
                        .toBeLessThan(1e-13);
                    expect(maxDiff(state.getWorldInertia(),
                        multiplyABT(mulMatrix(R, inertia), R)))
                        .toBeLessThan(1e-12);
                    expect(maxDiff(multiplyAB(state.getWorldInertia(),
                        state.getWorldInverseInertia()),
                        Matrix.identity(3, 3))).toBeLessThan(1e-11);
                });
        });

    it('setQOrientation normalizes only when asked, and the rotation always ' +
        'comes from the stored quaternion', () => {
            check(fc.tuple(unitQuaternion(), scaled(0.2, 5, 128),
                positiveInertia()), ([q, factor, inertia]) => {
                    const scaledQ = new Quaternion(q.values[0] * factor,
                        q.values[1] * factor, q.values[2] * factor,
                        q.values[3] * factor);

                    const normalized = new RigidBodyState();
                    normalized.setMass(1);
                    normalized.setBodyInertia(inertia);
                    normalized.setQOrientation(scaledQ, true);
                    expectClose(Math.sqrt(dot(normalized.getQOrientation(),
                        normalized.getQOrientation())), 1, 1e-14, 1e-14);
                    // The regression for upstream #313: the matrix is built
                    // from the normalized member, so it stays a rotation.
                    const R = normalized.getROrientation();
                    expect(maxDiff(multiplyABT(R, R), Matrix.identity(3, 3)))
                        .toBeLessThan(1e-13);
                    expect(maxDiff(R, Rotation.fromQuaternion(q, 3).toMatrix()))
                        .toBeLessThan(1e-13);

                    const raw = new RigidBodyState();
                    raw.setMass(1);
                    raw.setBodyInertia(inertia);
                    raw.setQOrientation(scaledQ);
                    expect(vecMaxDiff(raw.getQOrientation(), scaledQ)).toBe(0);
                });
        });

    it('setMass, setBodyInertia and the momentum setters honour immovability',
        () => {
            check(fc.tuple(scaled(-3, 3, 256), positiveInertia(), vector3()),
                ([mass, inertia, momentum]) => {
                    const state = new RigidBodyState();
                    state.setMass(mass);
                    state.setBodyInertia(inertia);
                    state.setLinearMomentum(momentum);
                    state.setAngularMomentum(momentum);
                    if (mass > 0) {
                        expect(state.isMovable()).toBe(true);
                        expect(state.isImmovable()).toBe(false);
                        expectClose(state.getInverseMass(), 1 / mass, 0, 1e-15);
                        expect(vecMaxDiff(state.getLinearMomentum(), momentum))
                            .toBe(0);
                    }
                    else {
                        expect(state.getMass()).toBe(0);
                        expect(state.getInverseMass()).toBe(0);
                        expect(state.isImmovable()).toBe(true);
                        // The setters are no-ops for an immovable body.
                        expect(state.getLinearMomentum().values
                            .every(x => x === 0)).toBe(true);
                        expect(state.getAngularMomentum().values
                            .every(x => x === 0)).toBe(true);
                    }
                });
        });

    it('free motion keeps the quaternion unit, the rotation orthonormal and ' +
        'the momenta exactly constant', () => {
            check(fc.tuple(freeBody(), scaled(0.005, 0.05, 64)),
                ([body, dt]) => {
                    const p0 = body.getLinearMomentum();
                    const l0 = body.getAngularMomentum();
                    const x0 = body.getPosition();
                    for (let k = 0; k < 25; ++k) {
                        body.update(k * dt, dt);
                        const q = body.getQOrientation();
                        // Every stage renormalizes, so the stored quaternion
                        // is unit after every step (upstream #313).
                        expectClose(Math.sqrt(dot(q, q)), 1, 1e-14, 1e-14);
                        const R = body.getROrientation();
                        expect(maxDiff(multiplyABT(R, R),
                            Matrix.identity(3, 3))).toBeLessThan(1e-12);
                    }
                    // Zero force and torque leave both momenta untouched.
                    expect(vecMaxDiff(body.getLinearMomentum(), p0)).toBe(0);
                    expect(vecMaxDiff(body.getAngularMomentum(), l0)).toBe(0);
                    // The center of mass travels in a straight line at the
                    // constant linear velocity.
                    const expected = add(x0, mul(body.getLinearVelocity(),
                        25 * dt));
                    expect(vecMaxDiff(body.getPosition(), expected))
                        .toBeLessThan(1e-12);
                });
        });

    it('conserves the kinetic energy of a torque-free body to the accuracy ' +
        'of the integrator', () => {
            // Free rigid-body motion conserves 0.5*m*|v|^2 + 0.5*w.L exactly;
            // the numerical drift is the RK4 truncation error, which grows
            // like (|w|*dt)^4 per step. The angular momentum and the step are
            // kept small so that the bound below is far tighter than the
            // error any coefficient mistake in the integrator would produce.
            check(fc.tuple(scaled(0.4, 4, 128), positiveInertia(),
                unitQuaternion(), vector3(0.5), scaled(0.002, 0.008, 64)),
                ([mass, inertia, q, angular, dt]) => {
                    const body = new RigidBody();
                    body.setMass(mass);
                    body.setBodyInertia(inertia);
                    body.setPosition(Vector.zero(3));
                    body.setQOrientation(q);
                    body.setLinearMomentum(Vector.zero(3));
                    body.setAngularMomentum(angular);
                    body.force = zeroForce;
                    body.torque = zeroForce;
                    const e0 = kineticEnergy(body);
                    for (let k = 0; k < 25; ++k) { body.update(k * dt, dt); }
                    expectClose(kineticEnergy(body), e0, 1e-12, 1e-9);
                }, 100);
        });

    it('integrates a constant force exactly (RK4 is exact for a quadratic)',
        () => {
            check(fc.tuple(scaled(0.4, 4, 128), positiveInertia(), vector3(),
                vector3(), vector3(), scaled(0.005, 0.05, 64)),
                ([mass, inertia, position, velocity, force, dt]) => {
                    const body = new RigidBody();
                    body.setMass(mass);
                    body.setBodyInertia(inertia);
                    body.setPosition(position);
                    body.setQOrientation(Quaternion.identity());
                    body.setLinearVelocity(velocity);
                    body.setAngularMomentum(Vector.zero(3));
                    body.force = () => force.clone();
                    body.torque = zeroForce;
                    const numSteps = 8;
                    for (let k = 0; k < numSteps; ++k) {
                        body.update(k * dt, dt);
                    }
                    const t = numSteps * dt;
                    const expected = add(add(position, mul(velocity, t)),
                        mul(force, (0.5 * t * t) / mass));
                    expect(vecMaxDiff(body.getPosition(), expected))
                        .toBeLessThan(1e-12);
                    expect(vecMaxDiff(body.getLinearMomentum(),
                        add(mul(velocity, mass), mul(force, t))))
                        .toBeLessThan(1e-12);
                });
        });

    it('rotates an isotropic torque-free body about a fixed axis at a ' +
        'constant rate', () => {
            check(fc.tuple(scaled(0.4, 4, 128), scaled(0.3, 3, 128),
                vector3(2), scaled(0.002, 0.02, 64)),
                ([mass, moment, angular, dt]) => {
                    const body = new RigidBody();
                    body.setMass(mass);
                    body.setBodyInertia(diag3(moment, moment, moment));
                    body.setPosition(Vector.zero(3));
                    body.setQOrientation(Quaternion.identity());
                    body.setAngularVelocity(angular);
                    body.force = zeroForce;
                    body.torque = zeroForce;
                    const w0 = body.getAngularVelocity();
                    const numSteps = 30;
                    for (let k = 0; k < numSteps; ++k) {
                        body.update(k * dt, dt);
                    }
                    // An isotropic body has no precession: the angular
                    // velocity is constant in the world frame.
                    expect(vecMaxDiff(body.getAngularVelocity(), w0))
                        .toBeLessThan(1e-11);
                    // The orientation is the rotation by |w|*t about w.
                    const speed = len(w0);
                    if (speed > 1e-3) {
                        const axis = mul(w0, 1 / speed);
                        const expected = Rotation.fromAxisAngle(
                            new AxisAngle(axis, speed * numSteps * dt))
                            .toMatrix();
                        // The residual is the RK4 truncation error over the
                        // run; the convergence-order property below pins the
                        // stage weights that produce it.
                        expect(maxDiff(body.getROrientation(), expected))
                            .toBeLessThan(5e-8);
                    }
                });
        });

    it('converges at fourth order in the step size', () => {
        // Integrating a torque-free isotropic body to a fixed final time with
        // n and 2n steps must reduce the error against the analytic rotation
        // by about 2^4; this is what distinguishes the RK4 stage weights from
        // any lower-order combination.
        check(fc.tuple(scaled(0.4, 4, 128), scaled(0.3, 3, 128), vector3(2)),
            ([mass, moment, angular]) => {
                const speed = len(angular);
                fc.pre(speed > 0.5);
                const finalTime = 1;
                const axis = mul(angular, 1 / speed);
                const exact = Rotation.fromAxisAngle(
                    new AxisAngle(axis, speed * finalTime)).toMatrix();
                const errorFor = (numSteps: number): number => {
                    const body = new RigidBody();
                    body.setMass(mass);
                    body.setBodyInertia(diag3(moment, moment, moment));
                    body.setPosition(Vector.zero(3));
                    body.setQOrientation(Quaternion.identity());
                    body.setAngularVelocity(angular);
                    body.force = zeroForce;
                    body.torque = zeroForce;
                    const dt = finalTime / numSteps;
                    for (let k = 0; k < numSteps; ++k) {
                        body.update(k * dt, dt);
                    }
                    return maxDiff(body.getROrientation(), exact);
                };
                const coarse = errorFor(20);
                const fine = errorFor(40);
                // Skip draws where the coarse error is already at round-off,
                // where the ratio carries no information.
                fc.pre(coarse > 1e-11);
                expect(fine).toBeLessThan(coarse / 8);
            }, 60);
    });

    it('the contact impulse conserves total linear and angular momentum',
        () => {
            check(fc.tuple(freeBody(), freeBody(), vector3(), unitVector3(),
                scaled(0, 1, 64)), ([A, B, P, N, restitution]) => {
                    const contact = new RigidBodyContact();
                    contact.A = A;
                    contact.B = B;
                    contact.P = P;
                    contact.N = N;
                    contact.restitution = restitution;

                    const pBefore = add(A.getLinearMomentum(),
                        B.getLinearMomentum());
                    // Total angular momentum about the origin includes the
                    // orbital term x_i x p_i.
                    const lBefore = add(
                        add(A.getAngularMomentum(),
                            cross(A.getPosition(), A.getLinearMomentum())),
                        add(B.getAngularMomentum(),
                            cross(B.getPosition(), B.getLinearMomentum())));

                    contact.applyImpulse();

                    const pAfter = add(A.getLinearMomentum(),
                        B.getLinearMomentum());
                    const lAfter = add(
                        add(A.getAngularMomentum(),
                            cross(A.getPosition(), A.getLinearMomentum())),
                        add(B.getAngularMomentum(),
                            cross(B.getPosition(), B.getLinearMomentum())));

                    const scale = Math.max(1, len(pBefore), len(lBefore));
                    expect(vecMaxDiff(pAfter, pBefore))
                        .toBeLessThan(1e-10 * scale);
                    expect(vecMaxDiff(lAfter, lBefore))
                        .toBeLessThan(1e-10 * scale);
                });
        });

    it('the contact impulse reverses the normal relative velocity by the ' +
        'restitution and preserves the tangential components', () => {
            check(fc.tuple(freeBody(), freeBody(), vector3(), unitVector3(),
                scaled(0, 1, 64)), ([A, B, P, N, restitution]) => {
                    const contact = new RigidBodyContact();
                    contact.A = A;
                    contact.B = B;
                    contact.P = P;
                    contact.N = N;
                    contact.restitution = restitution;

                    const rA = sub(P, A.getPosition());
                    const rB = sub(P, B.getPosition());
                    const relative = (): Vector => sub(
                        add(A.getLinearVelocity(),
                            cross(A.getAngularVelocity(), rA)),
                        add(B.getLinearVelocity(),
                            cross(B.getAngularVelocity(), rB)));

                    const before = relative();
                    const normalBefore = dot(N, before);
                    // T0 is the tangent direction the query builds; when the
                    // relative velocity is parallel to N the fallback branch
                    // runs and there is no tangential constraint.
                    const T0 = sub(before, mul(N, normalBefore));
                    const tangentLength = normalize(T0);

                    contact.applyImpulse();

                    const after = relative();
                    // Both branches solve for the normal component:
                    // N.v+ = -restitution * N.v-.
                    const scale = Math.max(1, len(before));
                    expectClose(dot(N, after), -restitution * normalBefore,
                        1e-10 * scale, 1e-9);
                    if (tangentLength > 0) {
                        // The 3x3 system also pins both tangential
                        // components of the relative velocity.
                        const T1 = cross(N, T0);
                        expectClose(dot(T0, after), dot(T0, before),
                            1e-10 * scale, 1e-9);
                        expectClose(dot(T1, after), dot(T1, before),
                            1e-10 * scale, 1e-9);
                    }
                });
        });

    it('rejects an unset force, torque or contact body', () => {
        const body = new RigidBody();
        body.setMass(1);
        expect(() => body.update(0, 0.1)).toThrow();
        body.force = zeroForce;
        expect(() => body.update(0, 0.1)).toThrow();
        body.torque = zeroForce;
        expect(() => body.update(0, 0.1)).not.toThrow();

        const contact = new RigidBodyContact();
        expect(() => contact.applyImpulse()).toThrow();
        contact.A = body;
        expect(() => contact.applyImpulse()).toThrow();
    });
});
