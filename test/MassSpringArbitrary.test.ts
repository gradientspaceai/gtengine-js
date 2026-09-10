import { describe, expect, it } from 'vitest';
import {
    MassSpringArbitrary,
    MassSpringArbitrarySpring
} from '../src/MassSpringArbitrary.js';
import { MassSpringCurve } from '../src/MassSpringCurve.js';
import { Vector, add, dot, length as vectorLength, mul, sub } from '../src/Vector.js';
import { check, expectClose, fc, scaled } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V41): properties of the arbitrary-topology mass-spring system.
// ---------------------------------------------------------------------------

// A random graph on numParticles nodes: distinct unordered pairs, positions,
// masses, spring constants and resting lengths.
interface Topology {
    numParticles: number;
    edges: [number, number][];
    positions: number[][];
    velocities: number[][];
    masses: number[];
    constants: number[];
    lengths: number[];
}

const coord = () => scaled(-4, 4, 512);

const topology = (): fc.Arbitrary<Topology> =>
    fc.integer({ min: 2, max: 5 }).chain(numParticles => {
        const pairs: [number, number][] = [];
        for (let i = 0; i < numParticles; ++i) {
            for (let j = i + 1; j < numParticles; ++j) { pairs.push([i, j]); }
        }
        return fc.record({
            numParticles: fc.constant(numParticles),
            edges: fc.uniqueArray(fc.constantFrom(...pairs),
                { minLength: 1, maxLength: pairs.length }),
            positions: fc.array(fc.array(coord(),
                { minLength: 2, maxLength: 2 }),
                { minLength: numParticles, maxLength: numParticles }),
            velocities: fc.array(fc.array(scaled(-1, 1, 256),
                { minLength: 2, maxLength: 2 }),
                { minLength: numParticles, maxLength: numParticles }),
            masses: fc.array(scaled(0.5, 3, 64),
                { minLength: numParticles, maxLength: numParticles }),
            constants: fc.array(scaled(0.1, 4, 64),
                { minLength: 1, maxLength: pairs.length }),
            lengths: fc.array(scaled(0.1, 3, 64),
                { minLength: 1, maxLength: pairs.length })
        });
    })
        // Reject configurations with coincident endpoints: the upstream
        // acceleration divides by the spring length, which is 0 there.
        .filter(t => t.edges.every(([i, j]) => {
            const dx = t.positions[j][0] - t.positions[i][0];
            const dy = t.positions[j][1] - t.positions[i][1];
            return dx * dx + dy * dy > 1e-4;
        }));

function buildSystem(t: Topology, step: number): TestMassSpringArbitrary {
    const system = new TestMassSpringArbitrary(2, t.numParticles,
        t.edges.length, step);
    for (let i = 0; i < t.numParticles; ++i) {
        system.setMass(i, t.masses[i]);
        system.setPosition(i, v2(t.positions[i][0], t.positions[i][1]));
        system.setVelocity(i, v2(t.velocities[i][0], t.velocities[i][1]));
    }
    for (let e = 0; e < t.edges.length; ++e) {
        system.setSpring(e, spring(t.edges[e][0], t.edges[e][1],
            t.constants[e % t.constants.length],
            t.lengths[e % t.lengths.length]));
    }
    return system;
}

function currentPositions(system: MassSpringArbitrary): Vector[] {
    const out: Vector[] = [];
    for (let i = 0; i < system.getNumParticles(); ++i) {
        out.push(system.getPosition(i).clone());
    }
    return out;
}

function currentVelocities(system: MassSpringArbitrary): Vector[] {
    const out: Vector[] = [];
    for (let i = 0; i < system.getNumParticles(); ++i) {
        out.push(system.getVelocity(i).clone());
    }
    return out;
}

// Total linear momentum sum_i m_i v_i.
function momentum(system: MassSpringArbitrary): Vector {
    let p = Vector.zero(2);
    for (let i = 0; i < system.getNumParticles(); ++i) {
        p = add(p, mul(system.getVelocity(i), system.getMass(i)));
    }
    return p;
}

// Kinetic plus spring potential energy.
function energy(system: MassSpringArbitrary): number {
    let e = 0;
    for (let i = 0; i < system.getNumParticles(); ++i) {
        const v = system.getVelocity(i);
        e += 0.5 * system.getMass(i) * dot(v, v);
    }
    for (let s = 0; s < system.getNumSprings(); ++s) {
        const sp = system.getSpring(s);
        const d = vectorLength(sub(system.getPosition(sp.particle1),
            system.getPosition(sp.particle0)));
        e += 0.5 * sp.constant * (d - sp.length) * (d - sp.length);
    }
    return e;
}

describe('MassSpringArbitrary verification', () => {
    it('acceleration equals the brute-force Hooke sum over adjacent springs',
        () => {
            check(topology(), (t) => {
                const system = buildSystem(t, 0.01);
                const position = currentPositions(system);
                const velocity = currentVelocities(system);
                for (let i = 0; i < t.numParticles; ++i) {
                    // Independent computation: sum over the springs that
                    // touch i of k*(1 - L/|d|)*d, divided by the mass.
                    let force = Vector.zero(2);
                    for (let e = 0; e < t.edges.length; ++e) {
                        const sp = system.getSpring(e);
                        let other = -1;
                        if (sp.particle0 === i) { other = sp.particle1; }
                        else if (sp.particle1 === i) { other = sp.particle0; }
                        if (other < 0) { continue; }
                        const d = sub(position[other], position[i]);
                        const ratio = sp.length / vectorLength(d);
                        force = add(force, mul(d, sp.constant * (1 - ratio)));
                    }
                    const expected = mul(force, 1 / t.masses[i]);
                    const actual = system.accelerationAt(i, 0, position,
                        velocity);
                    for (let k = 0; k < 2; ++k) {
                        expectClose(actual.get(k), expected.get(k), 1e-12,
                            1e-11);
                    }
                }
            }, 150);
        });

    it('obeys Newton third law: the spring forces on the two endpoints are ' +
        'equal and opposite', () => {
            const arb = fc.tuple(coord(), coord(), coord(), coord(),
                scaled(0.1, 4, 64), scaled(0.1, 3, 64), scaled(0.5, 3, 64),
                scaled(0.5, 3, 64))
                .filter(([x0, y0, x1, y1]) =>
                    (x1 - x0) ** 2 + (y1 - y0) ** 2 > 1e-4);
            check(arb, ([x0, y0, x1, y1, k, restLength, m0, m1]) => {
                const system = new TestMassSpringArbitrary(2, 2, 1, 0.01);
                system.setMass(0, m0);
                system.setMass(1, m1);
                system.setPosition(0, v2(x0, y0));
                system.setPosition(1, v2(x1, y1));
                system.setSpring(0, spring(0, 1, k, restLength));
                const position = currentPositions(system);
                const velocity = currentVelocities(system);
                const f0 = mul(system.accelerationAt(0, 0, position, velocity),
                    m0);
                const f1 = mul(system.accelerationAt(1, 0, position, velocity),
                    m1);
                for (let i = 0; i < 2; ++i) {
                    expectClose(f0.get(i), -f1.get(i), 1e-12, 1e-11);
                }
            });
        });

    it('a configuration at the resting lengths with zero velocity is an ' +
        'exact fixed point', () => {
            check(topology(), (t) => {
                const system = new TestMassSpringArbitrary(2, t.numParticles,
                    t.edges.length, 0.05);
                for (let i = 0; i < t.numParticles; ++i) {
                    system.setMass(i, t.masses[i]);
                    system.setPosition(i,
                        v2(t.positions[i][0], t.positions[i][1]));
                }
                for (let e = 0; e < t.edges.length; ++e) {
                    const [i, j] = t.edges[e];
                    const rest = vectorLength(sub(system.getPosition(j),
                        system.getPosition(i)));
                    system.setSpring(e, spring(i, j,
                        t.constants[e % t.constants.length], rest));
                }
                const position = currentPositions(system);
                const velocity = currentVelocities(system);
                for (let i = 0; i < t.numParticles; ++i) {
                    const a = system.accelerationAt(i, 0, position, velocity);
                    expect(a.get(0) + 0).toBe(0);
                    expect(a.get(1) + 0).toBe(0);
                }
                system.update(0);
                for (let i = 0; i < t.numParticles; ++i) {
                    expect(system.getPosition(i).get(0)).toBe(t.positions[i][0]);
                    expect(system.getPosition(i).get(1)).toBe(t.positions[i][1]);
                    expect(system.getVelocity(i).get(0) + 0).toBe(0);
                    expect(system.getVelocity(i).get(1) + 0).toBe(0);
                }
            }, 100);
        });

    it('conserves total linear momentum with no external force', () => {
        check(fc.tuple(topology(), scaled(0.002, 0.02, 32)), ([t, step]) => {
            const system = buildSystem(t, step);
            const before = momentum(system);
            for (let k = 0; k < 8; ++k) { system.update(k * step); }
            const after = momentum(system);
            let scale = 0;
            for (let i = 0; i < t.numParticles; ++i) {
                scale += t.masses[i]
                    * vectorLength(v2(t.velocities[i][0], t.velocities[i][1]));
            }
            // Newton third law makes sum_i m_i a_i vanish at every stage, so
            // the momentum drift is pure floating-point round-off.
            const tol = 1e-11 * Math.max(1, scale);
            expect(Math.abs(after.get(0) - before.get(0)))
                .toBeLessThanOrEqual(tol);
            expect(Math.abs(after.get(1) - before.get(1)))
                .toBeLessThanOrEqual(tol);
        }, 100);
    });

    it('keeps the total energy of an undamped system bounded over a run',
        () => {
            check(fc.tuple(topology(), scaled(0.001, 0.006, 32)),
                ([t, step]) => {
                    const system = buildSystem(t, step);
                    const e0 = energy(system);
                    let worst = 0;
                    for (let k = 0; k < 40; ++k) {
                        system.update(k * step);
                        worst = Math.max(worst,
                            Math.abs(energy(system) - e0));
                    }
                    // RK4 is not symplectic, but with these step sizes the
                    // energy drift over 40 steps stays far below 1% of the
                    // initial energy plus the spring scale.
                    const scale = Math.max(e0, 1);
                    expect(worst).toBeLessThanOrEqual(0.01 * scale);
                }, 60);
        });

    it('never moves immovable particles', () => {
        check(fc.tuple(topology(), scaled(0.002, 0.02, 32),
            fc.integer({ min: 0, max: 4 })), ([t, step, pinned]) => {
                const system = buildSystem(t, step);
                const fixedIndex = pinned % t.numParticles;
                system.setMass(fixedIndex, Number.MAX_VALUE);
                const p0 = system.getPosition(fixedIndex).clone();
                const v0 = system.getVelocity(fixedIndex).clone();
                for (let k = 0; k < 5; ++k) { system.update(k * step); }
                expect(system.getPosition(fixedIndex).get(0)).toBe(p0.get(0));
                expect(system.getPosition(fixedIndex).get(1)).toBe(p0.get(1));
                expect(system.getVelocity(fixedIndex).get(0)).toBe(v0.get(0));
                expect(system.getVelocity(fixedIndex).get(1)).toBe(v0.get(1));
            }, 100);
    });

    it('setSpring copies the spring and records both endpoints as adjacent',
        () => {
            check(fc.tuple(fc.integer({ min: 0, max: 3 }),
                fc.integer({ min: 0, max: 3 }), scaled(0.1, 3, 32),
                scaled(0.1, 3, 32)).filter(([i, j]) => i !== j),
                ([i, j, k, restLength]) => {
                    const system = new TestMassSpringArbitrary(2, 4, 1, 0.01);
                    const input = spring(i, j, k, restLength);
                    system.setSpring(0, input);
                    input.constant = -12345;
                    expect(system.getSpring(0).constant).toBe(k);
                    // Both endpoints see a nonzero force when the spring is
                    // stretched away from its resting length.
                    for (let n = 0; n < 4; ++n) {
                        system.setMass(n, 1);
                        system.setPosition(n, v2(n, 0));
                    }
                    // Place the endpoints 4 apart, which the generated
                    // resting length in (0.1, 3] can never equal, so both
                    // endpoints feel a nonzero force.
                    system.setPosition(i, v2(0, 0));
                    system.setPosition(j, v2(4, 0));
                    const position = currentPositions(system);
                    const velocity = currentVelocities(system);
                    const ai = system.accelerationAt(i, 0, position, velocity);
                    const aj = system.accelerationAt(j, 0, position, velocity);
                    expect(vectorLength(ai)).toBeGreaterThan(0);
                    expect(vectorLength(aj)).toBeGreaterThan(0);
                });
        });
});
