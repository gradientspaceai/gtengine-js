import { describe, expect, it } from 'vitest';
import {
    MassSpringArbitrary,
    MassSpringArbitrarySpring
} from '../src/MassSpringArbitrary.js';
import { MassSpringCurve } from '../src/MassSpringCurve.js';
import { Vector, add, dot, length as vectorLength, mul, sub } from '../src/Vector.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// Exposes the protected acceleration(...) callback for testing.
class TestMassSpringArbitrary extends MassSpringArbitrary {
    accelerationAt(i: number, time: number, position: readonly Vector[],
        velocity: readonly Vector[]): Vector {
        return this.acceleration(i, time, position, velocity);
    }
}

// A gravity-like constant external acceleration.
class GravityMassSpringArbitrary extends MassSpringArbitrary {
    constructor(dimension: number, numParticles: number, numSprings: number,
        step: number, private g: Vector) {
        super(dimension, numParticles, numSprings, step);
    }

    override externalAcceleration(_i: number, _time: number,
        _position: readonly Vector[], _velocity: readonly Vector[]): Vector {
        return this.g.clone();
    }

    accelerationAt(i: number, time: number, position: readonly Vector[],
        velocity: readonly Vector[]): Vector {
        return this.acceleration(i, time, position, velocity);
    }
}

function spring(p0: number, p1: number, constant: number,
    length: number): MassSpringArbitrarySpring {
    return new MassSpringArbitrarySpring(p0, p1, constant, length);
}

// An equilateral triangle of side 'side', centered at the origin.
function triangle(side: number, restLength: number, constant: number,
    mass: number, step: number): TestMassSpringArbitrary {
    const system = new TestMassSpringArbitrary(2, 3, 3, step);
    const radius = side / Math.sqrt(3);
    for (let i = 0; i < 3; ++i) {
        const angle = (2 * Math.PI * i) / 3 + Math.PI / 2;
        system.setMass(i, mass);
        system.setPosition(i, v2(radius * Math.cos(angle), radius * Math.sin(angle)));
        system.setVelocity(i, v2(0, 0));
    }
    system.setSpring(0, spring(0, 1, constant, restLength));
    system.setSpring(1, spring(1, 2, constant, restLength));
    system.setSpring(2, spring(2, 0, constant, restLength));
    return system;
}

describe('MassSpringArbitrary construction and member access', () => {
    it('reports the particle and spring counts', () => {
        const system = new MassSpringArbitrary(3, 5, 7, 0.01);
        expect(system.getDimension()).toBe(3);
        expect(system.getNumParticles()).toBe(5);
        expect(system.getNumSprings()).toBe(7);
        expect(system.getStep()).toBe(0.01);
    });

    it('default-constructs the springs to zero', () => {
        const system = new MassSpringArbitrary(2, 4, 3, 0.1);
        for (let i = 0; i < system.getNumSprings(); ++i) {
            const s = system.getSpring(i);
            expect(s.particle0).toBe(0);
            expect(s.particle1).toBe(0);
            expect(s.constant).toBe(0);
            expect(s.length).toBe(0);
        }
    });

    it('stores springs by value, so the caller may reuse the input', () => {
        const system = new MassSpringArbitrary(2, 4, 2, 0.1);
        const s = spring(0, 1, 10, 0.5);
        system.setSpring(0, s);
        s.particle1 = 3;
        s.constant = 99;
        expect(system.getSpring(0).particle1).toBe(1);
        expect(system.getSpring(0).constant).toBe(10);
        expect(system.getSpring(0).length).toBe(0.5);
    });

    it('the default external acceleration is zero', () => {
        const system = new MassSpringArbitrary(3, 3, 2, 0.1);
        expect(system.externalAcceleration(1, 0, [], []).values).toEqual([0, 0, 0]);
    });
});

describe('MassSpringArbitrary acceleration', () => {
    it('is zero when every spring is at its resting length', () => {
        const side = 2;
        const system = triangle(side, side, 5, 1, 0.01);
        for (let i = 0; i < 3; ++i) {
            const a = system.accelerationAt(i, 0, [
                system.getPosition(0), system.getPosition(1), system.getPosition(2)
            ], []);
            expect(vectorLength(a)).toBeLessThan(1e-12);
        }
    });

    it('is the Hooke force divided by the mass for a single spring', () => {
        const constant = 7;
        const restLength = 1;
        const mass = 2;
        const system = new TestMassSpringArbitrary(2, 2, 1, 0.01);
        system.setMass(0, mass);
        system.setMass(1, mass);
        const p0 = v2(0, 0);
        const p1 = v2(1.5, 0);
        system.setPosition(0, p0);
        system.setPosition(1, p1);
        system.setSpring(0, spring(0, 1, constant, restLength));

        // The stretched spring pulls each particle toward the other with
        // magnitude constant * (d - restLength) = 7 * 0.5 = 3.5, so the
        // accelerations are +/- 1.75 along the x-axis.
        const a0 = system.accelerationAt(0, 0, [p0, p1], []);
        const a1 = system.accelerationAt(1, 0, [p0, p1], []);
        expect(a0.values[0]).toBeCloseTo(1.75, 12);
        expect(a0.values[1]).toBeCloseTo(0, 12);
        expect(a1.values[0]).toBeCloseTo(-1.75, 12);
        expect(a1.values[1]).toBeCloseTo(0, 12);
    });

    it('sums the contributions of all adjacent springs', () => {
        // Particle 1 is connected to particles 0 and 2, both springs
        // stretched by 0.5 with constant 4 and unit mass. The forces are
        // 2 in -x and 2 in +y, so the accelerations add.
        const system = new TestMassSpringArbitrary(2, 3, 2, 0.01);
        for (let i = 0; i < 3; ++i) {
            system.setMass(i, 1);
        }
        const p = [v2(0, 0), v2(1.5, 0), v2(1.5, 1.5)];
        system.setSpring(0, spring(0, 1, 4, 1));
        system.setSpring(1, spring(1, 2, 4, 1));
        const a1 = system.accelerationAt(1, 0, p, []);
        expect(a1.values[0]).toBeCloseTo(-2, 12);
        expect(a1.values[1]).toBeCloseTo(2, 12);
    });

    it('adds the external acceleration', () => {
        const g = v2(0, -9.81);
        const system = new GravityMassSpringArbitrary(2, 2, 1, 0.01, g);
        system.setMass(0, 1);
        system.setMass(1, 1);
        const p = [v2(0, 0), v2(1, 0)];
        system.setSpring(0, spring(0, 1, 3, 1));
        // The spring is at rest, so only gravity remains.
        const a = system.accelerationAt(0, 0, p, []);
        expect(a.values[0]).toBeCloseTo(0, 12);
        expect(a.values[1]).toBeCloseTo(-9.81, 12);
    });
});

describe('MassSpringArbitrary dynamics', () => {
    it('keeps an equilateral triangle at rest in equilibrium', () => {
        const side = 2;
        const system = triangle(side, side, 5, 1, 0.005);
        const initial = [0, 1, 2].map(i => system.getPosition(i).clone());
        let time = 0;
        for (let k = 0; k < 200; ++k, time += 0.005) {
            system.update(time);
        }
        for (let i = 0; i < 3; ++i) {
            expect(vectorLength(sub(system.getPosition(i), initial[i])))
                .toBeLessThan(1e-10);
        }
    });

    it('expands a compressed triangle symmetrically, preserving the centroid',
        () => {
            const restLength = 2;
            const system = triangle(1.5, restLength, 5, 1, 0.002);
            let time = 0;
            for (let k = 0; k < 300; ++k, time += 0.002) {
                system.update(time);
            }
            // The three sides remain equal (symmetry) and are longer than the
            // initial 1.5.
            const p = [0, 1, 2].map(i => system.getPosition(i));
            const d01 = vectorLength(sub(p[0], p[1]));
            const d12 = vectorLength(sub(p[1], p[2]));
            const d20 = vectorLength(sub(p[2], p[0]));
            expect(d12).toBeCloseTo(d01, 9);
            expect(d20).toBeCloseTo(d01, 9);
            expect(d01).toBeGreaterThan(1.5);
            expect(d01).toBeLessThan(2 * restLength);

            // No external force acts, so the centroid does not move (it is at
            // the origin by construction).
            let centroid = new Vector(2);
            for (let i = 0; i < 3; ++i) {
                centroid = add(centroid, p[i]);
            }
            expect(vectorLength(centroid)).toBeLessThan(1e-10);
        });

    it('oscillates a fixed-mass/free-mass pair with the analytic period', () => {
        // Particle 0 is immovable (infinite mass); particle 1 has mass m and
        // is attached by a spring of constant k. In one dimension the
        // restoring force is exactly -k*(d - L), so the motion is simple
        // harmonic with period 2*pi*sqrt(m/k).
        const k = 8;
        const m = 2;
        const restLength = 1;
        const amplitude = 0.25;
        const period = 2 * Math.PI * Math.sqrt(m / k);
        const numSteps = 2000;
        const step = period / numSteps;

        const system = new MassSpringArbitrary(2, 2, 1, step);
        system.setMass(0, Number.MAX_VALUE);
        system.setMass(1, m);
        system.setPosition(0, v2(0, 0));
        system.setPosition(1, v2(restLength + amplitude, 0));
        system.setSpring(0, spring(0, 1, k, restLength));

        // After a quarter period the displaced mass is at the resting length
        // with maximum speed sqrt(k/m)*amplitude.
        let time = 0;
        for (let s = 0; s < numSteps / 4; ++s, time += step) {
            system.update(time);
        }
        expect(system.getPosition(1).values[0]).toBeCloseTo(restLength, 6);
        expect(system.getVelocity(1).values[0])
            .toBeCloseTo(-Math.sqrt(k / m) * amplitude, 6);

        // After a half period it is on the opposite side of the equilibrium.
        for (let s = 0; s < numSteps / 4; ++s, time += step) {
            system.update(time);
        }
        expect(system.getPosition(1).values[0])
            .toBeCloseTo(restLength - amplitude, 6);

        // After a full period it is back where it started, at rest, and the
        // immovable particle never moved.
        for (let s = 0; s < numSteps / 2; ++s, time += step) {
            system.update(time);
        }
        expect(system.getPosition(1).values[0])
            .toBeCloseTo(restLength + amplitude, 6);
        expect(system.getVelocity(1).values[0]).toBeCloseTo(0, 6);
        expect(system.getPosition(0).values).toEqual([0, 0]);
        expect(system.getVelocity(0).values).toEqual([0, 0]);
    });

    it('conserves energy for a free two-particle spring', () => {
        const k = 6;
        const m = 1.5;
        const restLength = 1;
        const step = 1e-3;
        const system = new MassSpringArbitrary(2, 2, 1, step);
        system.setMass(0, m);
        system.setMass(1, m);
        system.setPosition(0, v2(0, 0));
        system.setPosition(1, v2(1.4, 0));
        system.setSpring(0, spring(0, 1, k, restLength));

        const energy = (): number => {
            const d = vectorLength(sub(system.getPosition(1),
                system.getPosition(0))) - restLength;
            let kinetic = 0;
            for (let i = 0; i < 2; ++i) {
                const v = system.getVelocity(i);
                kinetic += 0.5 * m * dot(v, v);
            }
            return kinetic + 0.5 * k * d * d;
        };

        const initialEnergy = energy();
        let time = 0;
        for (let s = 0; s < 3000; ++s, time += step) {
            system.update(time);
            expect(Math.abs(energy() - initialEnergy)).toBeLessThan(1e-8);
        }

        // The center of mass of the isolated pair does not move.
        const center = mul(0.5, add(system.getPosition(0), system.getPosition(1)));
        expect(center.values[0]).toBeCloseTo(0.7, 10);
        expect(center.values[1]).toBeCloseTo(0, 10);
    });

    it('leaves all-immovable particles fixed', () => {
        const system = triangle(1.5, 2, 5, Number.MAX_VALUE, 0.01);
        const initial = [0, 1, 2].map(i => system.getPosition(i).clone());
        let time = 0;
        for (let s = 0; s < 50; ++s, time += 0.01) {
            system.update(time);
        }
        for (let i = 0; i < 3; ++i) {
            expect(system.getPosition(i).values).toEqual(initial[i].values);
            expect(system.getVelocity(i).values).toEqual([0, 0]);
        }
    });

    it('reproduces MassSpringCurve when the graph is a chain', () => {
        const numParticles = 6;
        const step = 1e-3;
        const constant = 12;
        const restLength = 1;
        const mass = 0.7;

        const curve = new MassSpringCurve(2, numParticles, step);
        const graph = new MassSpringArbitrary(2, numParticles,
            numParticles - 1, step);
        for (let i = 0; i < numParticles; ++i) {
            // A perturbed chain, so the dynamics are nontrivial.
            const p = v2(1.1 * i, 0.05 * Math.sin(i));
            curve.setMass(i, mass);
            curve.setPosition(i, p);
            graph.setMass(i, mass);
            graph.setPosition(i, p);
        }
        for (let i = 0; i + 1 < numParticles; ++i) {
            curve.setConstant(i, constant);
            curve.setLength(i, restLength);
            graph.setSpring(i, spring(i, i + 1, constant, restLength));
        }

        let time = 0;
        for (let s = 0; s < 200; ++s, time += step) {
            curve.update(time);
            graph.update(time);
        }
        for (let i = 0; i < numParticles; ++i) {
            expect(vectorLength(sub(curve.getPosition(i), graph.getPosition(i))))
                .toBeLessThan(1e-12);
            expect(vectorLength(sub(curve.getVelocity(i), graph.getVelocity(i))))
                .toBeLessThan(1e-12);
        }
    });
});
