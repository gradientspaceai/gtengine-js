import { describe, it, expect } from 'vitest';
import { ParticleSystem } from '../src/ParticleSystem';
import { Vector, add, sub, mul, dot } from '../src/Vector';

// A system whose particles all feel the same constant acceleration (for
// example gravity).
class ConstantAccelSystem extends ParticleSystem {
    private mAccel: Vector;

    constructor(dimension: number, numParticles: number, step: number,
        accel: Vector) {
        super(dimension, numParticles, step);
        this.mAccel = accel.clone();
    }

    protected acceleration(_i: number, _time: number,
        _position: readonly Vector[], _velocity: readonly Vector[]): Vector {
        return this.mAccel.clone();
    }
}

// A system driven by a caller-supplied acceleration function; the calls are
// recorded so the tests can inspect the Runge-Kutta stage evaluations.
interface AccelCall {
    i: number;
    time: number;
    position: Vector[];
    velocity: Vector[];
}

class CallbackSystem extends ParticleSystem {
    readonly calls: AccelCall[] = [];
    private mF: (i: number, time: number, position: readonly Vector[],
        velocity: readonly Vector[]) => Vector;

    constructor(dimension: number, numParticles: number, step: number,
        F: (i: number, time: number, position: readonly Vector[],
            velocity: readonly Vector[]) => Vector) {
        super(dimension, numParticles, step);
        this.mF = F;
    }

    protected acceleration(i: number, time: number,
        position: readonly Vector[], velocity: readonly Vector[]): Vector {
        this.calls.push({
            i,
            time,
            position: position.map((p) => p.clone()),
            velocity: velocity.map((v) => v.clone())
        });
        return this.mF(i, time, position, velocity);
    }
}

function vec(...values: number[]): Vector {
    return Vector.fromArray(values);
}

describe('ParticleSystem', () => {
    it('constructs with zero masses, positions and velocities', () => {
        const system = new ConstantAccelSystem(3, 4, 0.25, vec(0, 0, 0));
        expect(system.getDimension()).toBe(3);
        expect(system.getNumParticles()).toBe(4);
        expect(system.getStep()).toBe(0.25);
        for (let i = 0; i < 4; ++i) {
            expect(system.getMass(i)).toBe(0);
            expect(system.getPosition(i).values).toEqual([0, 0, 0]);
            expect(system.getVelocity(i).values).toEqual([0, 0, 0]);
        }
    });

    it('setMass stores finite positive masses and treats everything else as immovable', () => {
        const system = new ConstantAccelSystem(2, 6, 0.1, vec(0, -1));

        system.setMass(0, 2);
        expect(system.getMass(0)).toBe(2);

        system.setMass(1, 0);
        expect(system.getMass(1)).toBe(Number.MAX_VALUE);

        system.setMass(2, -3);
        expect(system.getMass(2)).toBe(Number.MAX_VALUE);

        system.setMass(3, Number.MAX_VALUE);
        expect(system.getMass(3)).toBe(Number.MAX_VALUE);

        system.setMass(4, Number.POSITIVE_INFINITY);
        expect(system.getMass(4)).toBe(Number.MAX_VALUE);

        system.setMass(5, Number.NaN);
        expect(system.getMass(5)).toBe(Number.MAX_VALUE);
    });

    it('setPosition and setVelocity copy their inputs (C++ value semantics)', () => {
        const system = new ConstantAccelSystem(3, 1, 0.1, vec(0, 0, 0));
        const p = vec(1, 2, 3);
        const v = vec(4, 5, 6);
        system.setPosition(0, p);
        system.setVelocity(0, v);

        p.set(0, 100);
        v.set(0, 200);

        expect(system.getPosition(0).values).toEqual([1, 2, 3]);
        expect(system.getVelocity(0).values).toEqual([4, 5, 6]);
    });

    it('setStep replaces the step (and the derived half/sixth steps)', () => {
        const accel = vec(0, -9.81, 0);
        const system = new ConstantAccelSystem(3, 1, 1, accel);
        system.setMass(0, 1);
        system.setStep(0.5);
        expect(system.getStep()).toBe(0.5);

        // A single update with the new step must advance by 0.5, not 1.
        system.update(0);
        const h = 0.5;
        expect(system.getPosition(0).get(1)).toBeCloseTo(0.5 * -9.81 * h * h, 12);
        expect(system.getVelocity(0).get(1)).toBeCloseTo(-9.81 * h, 12);
    });

    it('reproduces the analytic constant-acceleration trajectory', () => {
        // RK4 is exact for x" = a with constant a, so the numerical result
        // matches p0 + t*v0 + t^2*a/2 to round-off.
        const h = 1 / 64;
        const accel = vec(0.5, -9.81, 0.25);
        const system = new ConstantAccelSystem(3, 1, h, accel);
        system.setMass(0, 3);
        const p0 = vec(1, 2, -3);
        const v0 = vec(-0.5, 4, 1.5);
        system.setPosition(0, p0);
        system.setVelocity(0, v0);

        const numSteps = 64;
        for (let k = 0; k < numSteps; ++k) {
            system.update(k * h);
        }

        const t = numSteps * h;
        const expectedP = add(add(p0, mul(t, v0)), mul(0.5 * t * t, accel));
        const expectedV = add(v0, mul(t, accel));
        for (let d = 0; d < 3; ++d) {
            expect(system.getPosition(0).get(d)).toBeCloseTo(expectedP.get(d), 12);
            expect(system.getVelocity(0).get(d)).toBeCloseTo(expectedV.get(d), 12);
        }
    });

    it('moves a particle in a straight line when the acceleration is zero', () => {
        const h = 0.1;
        const system = new ConstantAccelSystem(2, 1, h, vec(0, 0));
        system.setMass(0, 1);
        system.setPosition(0, vec(-1, 5));
        system.setVelocity(0, vec(2, -3));

        for (let k = 0; k < 20; ++k) {
            system.update(k * h);
        }

        expect(system.getPosition(0).get(0)).toBeCloseTo(-1 + 2 * 2, 12);
        expect(system.getPosition(0).get(1)).toBeCloseTo(5 - 3 * 2, 12);
        expect(system.getVelocity(0).values).toEqual([2, -3]);
    });

    it('keeps immovable particles fixed and lets movable ones fall', () => {
        const h = 0.05;
        const system = new ConstantAccelSystem(3, 3, h, vec(0, -10, 0));
        // Particle 0: pinned by the documented "infinite" mass.
        system.setMass(0, Number.MAX_VALUE);
        // Particle 1: movable.
        system.setMass(1, 1);
        // Particle 2: mass never assigned, so it is pinned as well (the
        // zero-initialized mass has zero inverse mass).
        system.setPosition(0, vec(1, 1, 1));
        system.setVelocity(0, vec(7, 7, 7));
        system.setPosition(1, vec(2, 2, 2));
        system.setPosition(2, vec(3, 3, 3));

        for (let k = 0; k < 10; ++k) {
            system.update(k * h);
        }

        // The pinned particles keep both their positions and their
        // velocities: upstream skips them in the integration entirely.
        expect(system.getPosition(0).values).toEqual([1, 1, 1]);
        expect(system.getVelocity(0).values).toEqual([7, 7, 7]);
        expect(system.getPosition(2).values).toEqual([3, 3, 3]);

        const t = 10 * h;
        expect(system.getPosition(1).get(1)).toBeCloseTo(2 - 5 * t * t, 12);
        expect(system.getVelocity(1).get(1)).toBeCloseTo(-10 * t, 12);
    });

    it('evaluates the acceleration at the four Runge-Kutta stage times', () => {
        const h = 0.2;
        const system = new CallbackSystem(2, 1, h, () => vec(0, 0));
        system.setMass(0, 1);
        system.update(5);

        expect(system.calls.length).toBe(4);
        expect(system.calls.map((c) => c.time)).toEqual([5, 5.1, 5.1, 5.2]);
        expect(system.calls.every((c) => c.i === 0)).toBe(true);
    });

    it('never calls the acceleration hook for immovable particles', () => {
        const h = 0.1;
        const system = new CallbackSystem(2, 3, h, () => vec(1, 1));
        system.setMass(0, Number.MAX_VALUE);
        system.setMass(1, 2);
        system.setMass(2, 4);
        system.update(0);

        expect(system.calls.length).toBe(8);
        expect(system.calls.some((c) => c.i === 0)).toBe(false);
        expect(system.calls.filter((c) => c.i === 1).length).toBe(4);
        expect(system.calls.filter((c) => c.i === 2).length).toBe(4);
    });

    it('passes the current state to the first stage and intermediate states later', () => {
        const h = 0.4;
        const accel = vec(0, -2);
        const system = new CallbackSystem(2, 1, h, () => accel.clone());
        system.setMass(0, 1);
        system.setPosition(0, vec(3, 7));
        system.setVelocity(0, vec(1, 0));
        system.update(0);

        // Stage 1 sees the current state exactly.
        expect(system.calls[0].position[0].values).toEqual([3, 7]);
        expect(system.calls[0].velocity[0].values).toEqual([1, 0]);

        // Stage 2 sees p + (h/2)*v and v + (h/2)*a.
        expect(system.calls[1].position[0].get(0)).toBeCloseTo(3 + 0.2 * 1, 12);
        expect(system.calls[1].position[0].get(1)).toBeCloseTo(7, 12);
        expect(system.calls[1].velocity[0].get(1)).toBeCloseTo(0 + 0.2 * -2, 12);

        // Stage 4 sees p + h*(stage-3 velocity) and v + h*a.
        expect(system.calls[3].velocity[0].get(1)).toBeCloseTo(0 + 0.4 * -2, 12);

        // The pinned-particle branch aside, the recorded states must never be
        // the live internal vectors.
        expect(system.calls[0].position[0]).not.toBe(system.getPosition(0));
    });

    it('solves the harmonic oscillator with fourth-order accuracy', () => {
        // x" = -x with x(0) = 1, x'(0) = 0 has the solution x(t) = cos(t).
        const solveError = (numSteps: number): number => {
            const h = 1 / numSteps;
            const system = new CallbackSystem(1, 1, h,
                (i, _t, position) => vec(-position[i].get(0)));
            system.setMass(0, 1);
            system.setPosition(0, vec(1));
            system.setVelocity(0, vec(0));
            for (let k = 0; k < numSteps; ++k) {
                system.update(k * h);
            }
            return Math.abs(system.getPosition(0).get(0) - Math.cos(1));
        };

        const e1 = solveError(8);
        const e2 = solveError(16);
        expect(e1).toBeLessThan(1e-5);
        // Halving the step must reduce the error by roughly 2^4 = 16.
        expect(e2).toBeLessThan(e1 / 10);
        expect(solveError(256)).toBeLessThan(1e-11);
    });

    it('handles a velocity-dependent acceleration (linear drag)', () => {
        // v' = -k*v has the solution v(t) = v0*exp(-k*t) and, with x(0) = 0,
        // x(t) = (v0/k)*(1 - exp(-k*t)).
        const k = 1.5;
        const v0 = 4;
        const numSteps = 100;
        const h = 2 / numSteps;
        const system = new CallbackSystem(1, 1, h,
            (i, _t, _p, velocity) => vec(-k * velocity[i].get(0)));
        system.setMass(0, 1);
        system.setVelocity(0, vec(v0));

        for (let j = 0; j < numSteps; ++j) {
            system.update(j * h);
        }

        const t = 2;
        expect(Math.abs(system.getVelocity(0).get(0)
            - v0 * Math.exp(-k * t))).toBeLessThan(1e-7);
        expect(Math.abs(system.getPosition(0).get(0)
            - (v0 / k) * (1 - Math.exp(-k * t)))).toBeLessThan(1e-7);
    });

    it('handles a time-dependent acceleration', () => {
        // x" = t with x(0) = x'(0) = 0 has x(t) = t^3/6 and x'(t) = t^2/2.
        const numSteps = 50;
        const h = 1 / numSteps;
        const system = new CallbackSystem(1, 1, h, (_i, time) => vec(time));
        system.setMass(0, 1);

        for (let j = 0; j < numSteps; ++j) {
            system.update(j * h);
        }

        expect(system.getVelocity(0).get(0)).toBeCloseTo(0.5, 10);
        expect(system.getPosition(0).get(0)).toBeCloseTo(1 / 6, 10);
    });

    it('integrates interacting particles and conserves linear momentum', () => {
        // A spring between two particles produces equal and opposite forces,
        // so the total momentum m0*v0 + m1*v1 is invariant.
        const m0 = 2;
        const m1 = 5;
        const stiffness = 3;
        const masses = [m0, m1];
        const h = 0.01;
        const system = new CallbackSystem(2, 2, h, (i, _t, position) => {
            const other = 1 - i;
            const force = mul(stiffness, sub(position[other], position[i]));
            return mul(1 / masses[i], force);
        });
        system.setMass(0, m0);
        system.setMass(1, m1);
        system.setPosition(0, vec(-1, 0));
        system.setPosition(1, vec(1, 0));
        system.setVelocity(0, vec(0, 1));
        system.setVelocity(1, vec(0, -2));

        const momentum = (): Vector => add(mul(m0, system.getVelocity(0)),
            mul(m1, system.getVelocity(1)));
        const p0 = momentum();

        for (let j = 0; j < 200; ++j) {
            system.update(j * h);
        }

        const p1 = momentum();
        expect(p1.get(0)).toBeCloseTo(p0.get(0), 10);
        expect(p1.get(1)).toBeCloseTo(p0.get(1), 10);

        // Cross-check against the analytic two-body solution. The center of
        // mass X moves at constant velocity and the relative coordinate
        // r = p1 - p0 satisfies r" = -w^2*r with w^2 = k*(1/m0 + 1/m1).
        const t = 200 * h;
        const total = m0 + m1;
        const w = Math.sqrt(stiffness * (1 / m0 + 1 / m1));
        const r0 = vec(2, 0);
        const dr0 = vec(0, -3);
        const r = add(mul(Math.cos(w * t), r0), mul(Math.sin(w * t) / w, dr0));
        const x0 = vec(3 / 7, 0);
        const dx0 = vec(0, -8 / 7);
        const com = add(x0, mul(t, dx0));

        const expected0 = sub(com, mul(m1 / total, r));
        const expected1 = add(com, mul(m0 / total, r));
        for (let d = 0; d < 2; ++d) {
            expect(Math.abs(system.getPosition(0).get(d)
                - expected0.get(d))).toBeLessThan(1e-8);
            expect(Math.abs(system.getPosition(1).get(d)
                - expected1.get(d))).toBeLessThan(1e-8);
        }

        // The particles actually moved apart from their initial separation.
        const separation = Math.sqrt(dot(sub(system.getPosition(1),
            system.getPosition(0)), sub(system.getPosition(1),
                system.getPosition(0))));
        expect(Math.abs(separation - 2)).toBeGreaterThan(1e-3);
    });

    it('integrates several independent particles at once', () => {
        const h = 0.25;
        const accel = vec(0, -8, 0);
        const system = new ConstantAccelSystem(3, 3, h, accel);
        for (let i = 0; i < 3; ++i) {
            system.setMass(i, i + 1);
            system.setPosition(i, vec(i, 0, 0));
            system.setVelocity(i, vec(0, i, 0));
        }

        const numSteps = 8;
        for (let k = 0; k < numSteps; ++k) {
            system.update(k * h);
        }

        const t = numSteps * h;
        for (let i = 0; i < 3; ++i) {
            expect(system.getPosition(i).get(0)).toBeCloseTo(i, 12);
            expect(system.getPosition(i).get(1)).toBeCloseTo(i * t - 4 * t * t, 12);
            expect(system.getVelocity(i).get(1)).toBeCloseTo(i - 8 * t, 12);
        }
    });

    it('supports a derived class that overrides update for pre/post semantics', () => {
        const events: string[] = [];
        class BoundedSystem extends ConstantAccelSystem {
            override update(time: number): void {
                events.push(`pre:${time}`);
                super.update(time);
                // Post-update: clamp the height to a floor at y = 0.
                if (this.mPosition[0].get(1) < 0) {
                    this.mPosition[0].set(1, 0);
                    this.mVelocity[0].set(1, 0);
                }
                events.push(`post:${time}`);
            }
        }

        const h = 0.5;
        const system = new BoundedSystem(2, 1, h, vec(0, -10));
        system.setMass(0, 1);
        system.setPosition(0, vec(0, 1));

        for (let k = 0; k < 4; ++k) {
            system.update(k * h);
        }

        expect(events).toEqual(['pre:0', 'post:0', 'pre:0.5', 'post:0.5',
            'pre:1', 'post:1', 'pre:1.5', 'post:1.5']);
        expect(system.getPosition(0).get(1)).toBe(0);
        expect(system.getVelocity(0).get(1)).toBe(0);
    });
});
