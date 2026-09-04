import { describe, expect, it } from 'vitest';
import { MassSpringCurve } from '../src/MassSpringCurve.js';
import { Vector, length as vectorLength, sub } from '../src/Vector.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// A uniform chain of 'numParticles' masses spaced 'spacing' apart along the
// x-axis, with every spring of constant 'constant' and resting length
// 'restLength'.
function makeChain(numParticles: number, spacing: number, restLength: number,
    constant: number, mass: number, step: number): MassSpringCurve {
    const chain = new MassSpringCurve(2, numParticles, step);
    for (let i = 0; i < numParticles; ++i) {
        chain.setMass(i, mass);
        chain.setPosition(i, v2(i * spacing, 0));
        chain.setVelocity(i, v2(0, 0));
    }
    for (let i = 0; i < numParticles - 1; ++i) {
        chain.setConstant(i, constant);
        chain.setLength(i, restLength);
    }
    return chain;
}

describe('MassSpringCurve construction and member access', () => {
    it('has one fewer spring than particles', () => {
        const chain = new MassSpringCurve(3, 6, 0.01);
        expect(chain.getNumParticles()).toBe(6);
        expect(chain.getNumSprings()).toBe(5);
        expect(chain.getDimension()).toBe(3);
        expect(chain.getStep()).toBe(0.01);
    });

    it('spring constants and resting lengths default to zero', () => {
        const chain = new MassSpringCurve(2, 4, 0.1);
        for (let i = 0; i < chain.getNumSprings(); ++i) {
            expect(chain.getConstant(i)).toBe(0);
            expect(chain.getLength(i)).toBe(0);
        }
    });

    it('stores the spring constants and resting lengths', () => {
        const chain = new MassSpringCurve(2, 4, 0.1);
        for (let i = 0; i < chain.getNumSprings(); ++i) {
            chain.setConstant(i, 2 * i + 1);
            chain.setLength(i, 0.5 * (i + 1));
        }
        for (let i = 0; i < chain.getNumSprings(); ++i) {
            expect(chain.getConstant(i)).toBe(2 * i + 1);
            expect(chain.getLength(i)).toBe(0.5 * (i + 1));
        }
    });

    it('the default external acceleration is zero', () => {
        const chain = new MassSpringCurve(3, 3, 0.1);
        const a = chain.externalAcceleration(1, 0, [], []);
        expect(a.values).toEqual([0, 0, 0]);
    });
});

describe('MassSpringCurve equilibrium', () => {
    it('a chain at its resting lengths does not move', () => {
        const spacing = 1.5;
        const chain = makeChain(6, spacing, spacing, 10, 1, 0.01);
        for (let n = 0; n < 25; ++n) {
            chain.update(n * chain.getStep());
        }
        for (let i = 0; i < chain.getNumParticles(); ++i) {
            expect(chain.getPosition(i).get(0)).toBeCloseTo(i * spacing, 12);
            expect(chain.getPosition(i).get(1)).toBeCloseTo(0, 12);
            expect(vectorLength(chain.getVelocity(i))).toBeCloseTo(0, 12);
        }
    });

    it('a chain with nonuniform resting lengths is also at equilibrium', () => {
        // Springs at their individual resting lengths exert no force even when
        // the lengths differ.
        const rest = [1, 2, 0.5];
        const chain = new MassSpringCurve(2, 4, 0.01);
        let x = 0;
        for (let i = 0; i < 4; ++i) {
            chain.setMass(i, 1 + i);
            chain.setPosition(i, v2(x, 0));
            if (i < 3) {
                chain.setConstant(i, 3 + i);
                chain.setLength(i, rest[i]);
                x += rest[i];
            }
        }
        const initial = [0, 1, 3, 3.5];
        for (let n = 0; n < 10; ++n) {
            chain.update(n * 0.01);
        }
        for (let i = 0; i < 4; ++i) {
            expect(chain.getPosition(i).get(0)).toBeCloseTo(initial[i], 12);
        }
    });
});

describe('MassSpringCurve dynamics', () => {
    it('a stretched chain contracts symmetrically', () => {
        // Three particles at x = 0, 2, 4 with resting lengths 1: the ends pull
        // inward and, by symmetry, the middle particle does not move.
        const chain = makeChain(3, 2, 1, 4, 1, 0.005);
        for (let n = 0; n < 20; ++n) {
            chain.update(n * chain.getStep());
        }
        const p0 = chain.getPosition(0).get(0);
        const p1 = chain.getPosition(1).get(0);
        const p2 = chain.getPosition(2).get(0);
        expect(p0).toBeGreaterThan(0);
        expect(p2).toBeLessThan(4);
        expect(p1).toBeCloseTo(2, 10);
        // Symmetry about x = 2.
        expect(p0 + p2).toBeCloseTo(4, 10);
        // Also symmetric in the velocities.
        expect(chain.getVelocity(0).get(0) + chain.getVelocity(2).get(0))
            .toBeCloseTo(0, 10);
    });

    it('a compressed chain expands', () => {
        const chain = makeChain(3, 1, 2, 4, 1, 0.005);
        chain.update(0);
        expect(chain.getPosition(0).get(0)).toBeLessThan(0);
        expect(chain.getPosition(2).get(0)).toBeGreaterThan(2);
    });

    it('immovable endpoints stay fixed while the interior moves', () => {
        const chain = makeChain(5, 2, 1, 4, 1, 0.005);
        chain.setMass(0, Number.MAX_VALUE);
        chain.setMass(4, Number.MAX_VALUE);
        // A uniformly stretched chain is already an equilibrium (every
        // interior particle feels equal and opposite spring forces), so
        // displace one particle to make the interior move.
        chain.setPosition(2, v2(4, 0.5));
        for (let n = 0; n < 30; ++n) {
            chain.update(n * chain.getStep());
        }
        expect(chain.getPosition(0).values).toEqual([0, 0]);
        expect(chain.getPosition(4).values).toEqual([8, 0]);
        // The displaced particle is pulled back toward the axis.
        expect(chain.getPosition(2).get(1)).toBeLessThan(0.5);
        expect(chain.getPosition(2).get(1)).toBeGreaterThan(0);
        // Its neighbors are pulled off the axis.
        expect(chain.getPosition(1).get(1)).toBeGreaterThan(0);
        expect(chain.getPosition(3).get(1)).toBeGreaterThan(0);
    });

    it('a single spring is a harmonic oscillator with one fixed end', () => {
        // Particle 0 immovable at the origin, particle 1 of mass m attached by
        // a spring of constant c and resting length L. With x(0) = L + A and
        // zero initial velocity, x(t) = L + A*cos(w*t), w = sqrt(c/m).
        const m = 2;
        const c = 8;
        const L = 1;
        const A = 0.25;
        const w = Math.sqrt(c / m);
        const step = 1e-4;
        const chain = new MassSpringCurve(2, 2, step);
        chain.setMass(0, Number.MAX_VALUE);
        chain.setMass(1, m);
        chain.setPosition(0, v2(0, 0));
        chain.setPosition(1, v2(L + A, 0));
        chain.setConstant(0, c);
        chain.setLength(0, L);

        const numSteps = 2000;
        for (let n = 0; n < numSteps; ++n) {
            chain.update(n * step);
        }
        const t = numSteps * step;
        expect(chain.getPosition(1).get(0))
            .toBeCloseTo(L + A * Math.cos(w * t), 8);
        expect(chain.getVelocity(1).get(0))
            .toBeCloseTo(-A * w * Math.sin(w * t), 8);
    });
});

describe('MassSpringCurve external acceleration', () => {
    it('constant gravity gives the exact free-fall displacement', () => {
        // With every spring at its resting length the only acceleration is the
        // external one; RK4 is exact for constant acceleration, so after one
        // step of size h the displacement is h^2*a/2 and the speed is h*a.
        const g = -9.8;
        const h = 0.05;

        class FallingChain extends MassSpringCurve {
            override externalAcceleration(): Vector {
                return v2(0, g);
            }
        }

        const chain = new FallingChain(2, 4, h);
        for (let i = 0; i < 4; ++i) {
            chain.setMass(i, 1);
            chain.setPosition(i, v2(i, 0));
        }
        for (let i = 0; i < 3; ++i) {
            chain.setConstant(i, 5);
            chain.setLength(i, 1);
        }

        chain.update(0);
        for (let i = 0; i < 4; ++i) {
            expect(chain.getPosition(i).get(0)).toBeCloseTo(i, 12);
            expect(chain.getPosition(i).get(1)).toBeCloseTo(0.5 * h * h * g, 14);
            expect(chain.getVelocity(i).get(1)).toBeCloseTo(h * g, 14);
        }
    });

    it('external acceleration is scaled by nothing (it is already F/m)', () => {
        // The spring force is divided by the mass, but the external term is an
        // acceleration already, so two chains with different masses fall the
        // same way.
        const g = -1;

        class FallingChain extends MassSpringCurve {
            override externalAcceleration(): Vector {
                return v2(0, g);
            }
        }

        const light = new FallingChain(2, 2, 0.1);
        const heavy = new FallingChain(2, 2, 0.1);
        for (const [chain, mass] of [[light, 1], [heavy, 1000]] as const) {
            chain.setMass(0, mass);
            chain.setMass(1, mass);
            chain.setPosition(0, v2(0, 0));
            chain.setPosition(1, v2(1, 0));
            chain.setConstant(0, 7);
            chain.setLength(0, 1);
            chain.update(0);
        }
        expect(light.getPosition(0).get(1))
            .toBeCloseTo(heavy.getPosition(0).get(1), 14);
    });
});

describe('MassSpringCurve energy', () => {
    it('total energy is nearly conserved over many small steps', () => {
        const m = 1;
        const c = 20;
        const L = 1;
        const step = 1e-4;
        const numParticles = 5;
        const chain = new MassSpringCurve(2, numParticles, step);
        for (let i = 0; i < numParticles; ++i) {
            chain.setMass(i, m);
            chain.setPosition(i, v2(i * L, i === 2 ? 0.3 : 0));
            chain.setVelocity(i, v2(0, 0));
        }
        for (let i = 0; i < numParticles - 1; ++i) {
            chain.setConstant(i, c);
            chain.setLength(i, L);
        }

        const energy = (): number => {
            let e = 0;
            for (let i = 0; i < numParticles; ++i) {
                const v = chain.getVelocity(i);
                e += 0.5 * m * (v.get(0) * v.get(0) + v.get(1) * v.get(1));
            }
            for (let i = 0; i + 1 < numParticles; ++i) {
                const d = vectorLength(sub(chain.getPosition(i + 1),
                    chain.getPosition(i)));
                e += 0.5 * c * (d - L) * (d - L);
            }
            return e;
        };

        const e0 = energy();
        for (let n = 0; n < 3000; ++n) {
            chain.update(n * step);
        }
        expect(energy()).toBeCloseTo(e0, 6);
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). MassSpringCurve.h was read
// line by line against src/MassSpringCurve.ts. The whole file is the
// Acceleration callback, so the properties below check it against the physics
// it is supposed to implement: the force is the negative gradient of the
// Hooke potential, the forces on the two ends of a spring cancel, and the
// result is invariant under rigid motions of the configuration.
import {
    check, fc, expectClose, wellScaledVector, rotationFrame
} from './helpers/arbitraries.js';

// Exposes the protected acceleration callback.
class ProbeChain extends MassSpringCurve {
    accelerationAt(i: number, time: number, position: readonly Vector[],
        velocity: readonly Vector[]): Vector {
        return this.acceleration(i, time, position, velocity);
    }
}

interface ChainData {
    chain: ProbeChain;
    position: Vector[];
    mass: number[];
    constant: number[];
    restLength: number[];
}

// Random chains whose particles are well separated: Length(diff) appears in a
// denominator, so nearly coincident particles make the force arbitrarily
// large and the finite-difference cross-check meaningless.
const chainData = (dim: number): fc.Arbitrary<ChainData> =>
    fc.integer({ min: 2, max: 6 }).chain(n => fc.tuple(
        fc.array(wellScaledVector(dim, -4, 4), { minLength: n, maxLength: n }),
        fc.array(fc.double({ min: 0.5, max: 4, noNaN: true,
            noDefaultInfinity: true }), { minLength: n, maxLength: n }),
        fc.array(fc.double({ min: 0.5, max: 5, noNaN: true,
            noDefaultInfinity: true }), { minLength: n - 1, maxLength: n - 1 }),
        fc.array(fc.double({ min: 0.2, max: 3, noNaN: true,
            noDefaultInfinity: true }), { minLength: n - 1, maxLength: n - 1 })))
        .map(([position, mass, constant, restLength]) => {
            const chain = new ProbeChain(dim, position.length, 0.01);
            for (let i = 0; i < position.length; ++i) {
                chain.setMass(i, mass[i]);
                chain.setPosition(i, position[i]);
                chain.setVelocity(i, new Vector(dim));
            }
            for (let i = 0; i + 1 < position.length; ++i) {
                chain.setConstant(i, constant[i]);
                chain.setLength(i, restLength[i]);
            }
            return { chain, position, mass, constant, restLength };
        })
        .filter(d => {
            for (let i = 0; i + 1 < d.position.length; ++i) {
                if (vectorLength(sub(d.position[i + 1], d.position[i])) < 0.5) {
                    return false;
                }
            }
            return true;
        });

// The Hooke potential of the whole chain, U = sum_j c_j (|d_j| - L_j)^2 / 2.
function potential(d: ChainData, position: readonly Vector[]): number {
    let u = 0;
    for (let j = 0; j + 1 < position.length; ++j) {
        const len = vectorLength(sub(position[j + 1], position[j]));
        u += 0.5 * d.constant[j] * (len - d.restLength[j]) ** 2;
    }
    return u;
}

describe('MassSpringCurve verification', () => {
    it('acceleration is minus the potential gradient over the mass', () => {
        // m_i * a_i = -dU/dP_i for the Hooke potential; the gradient is taken
        // by centered differences of an independently written U.
        check(chainData(3), d => {
            const dim = 3;
            const h = 1e-5;
            for (let i = 0; i < d.position.length; ++i) {
                const a = d.chain.accelerationAt(i, 0, d.position, d.position);
                for (let k = 0; k < dim; ++k) {
                    const plus = d.position.map(p => p.clone());
                    const minus = d.position.map(p => p.clone());
                    plus[i].values[k] += h;
                    minus[i].values[k] -= h;
                    const grad = (potential(d, plus) - potential(d, minus)) /
                        (2 * h);
                    expectClose(d.mass[i] * a.values[k], -grad, 1e-6, 1e-6);
                }
            }
        });
    });

    it('the spring forces of the whole chain sum to zero', () => {
        // Newton's third law: with no external force the total force on the
        // chain vanishes, so sum_i m_i a_i = 0. This catches a sign or index
        // error in either branch of Acceleration.
        check(chainData(3), d => {
            const total = new Vector(3);
            for (let i = 0; i < d.position.length; ++i) {
                const a = d.chain.accelerationAt(i, 0, d.position, d.position);
                for (let k = 0; k < 3; ++k) {
                    total.values[k] += d.mass[i] * a.values[k];
                }
            }
            const scale = Math.max(1, ...d.constant) *
                Math.max(1, ...d.restLength);
            for (let k = 0; k < 3; ++k) {
                expect(Math.abs(total.values[k])).toBeLessThan(1e-9 * scale);
            }
        });
    });

    it('a chain at its resting lengths has zero acceleration', () => {
        // Place the particles on a ray so that consecutive distances are
        // exactly the resting lengths.
        check(fc.tuple(fc.integer({ min: 2, max: 6 }),
            fc.array(fc.double({ min: 0.3, max: 3, noNaN: true,
                noDefaultInfinity: true }), { minLength: 5, maxLength: 5 }),
            fc.array(fc.double({ min: 0.5, max: 4, noNaN: true,
                noDefaultInfinity: true }), { minLength: 6, maxLength: 6 })),
        ([n, lengths, constants]) => {
            const chain = new ProbeChain(2, n, 0.01);
            const position: Vector[] = [];
            let x = 0;
            for (let i = 0; i < n; ++i) {
                position.push(Vector.fromArray([x, 0]));
                chain.setMass(i, 1 + i);
                chain.setPosition(i, position[i]);
                chain.setVelocity(i, new Vector(2));
                if (i + 1 < n) { x += lengths[i]; }
            }
            for (let i = 0; i + 1 < n; ++i) {
                chain.setConstant(i, constants[i]);
                chain.setLength(i, lengths[i]);
            }
            for (let i = 0; i < n; ++i) {
                const a = chain.accelerationAt(i, 0, position, position);
                expectClose(a.values[0], 0, 1e-12, 1e-12);
                expectClose(a.values[1], 0, 1e-12, 1e-12);
            }
        });
    });

    it('is invariant under translation and equivariant under rotation', () => {
        check(fc.tuple(chainData(3), wellScaledVector(3, -6, 6),
            rotationFrame(3)), ([d, shift, frame]) => {
            const moved: Vector[] = d.position.map(p => {
                const q = Vector.fromArray([0, 1, 2].map(k =>
                    frame[0].values[k] * p.values[0] +
                    frame[1].values[k] * p.values[1] +
                    frame[2].values[k] * p.values[2] + shift.values[k]));
                return q;
            });
            for (let i = 0; i < d.position.length; ++i) {
                const a = d.chain.accelerationAt(i, 0, d.position, d.position);
                const b = d.chain.accelerationAt(i, 0, moved, moved);
                const rotated = [0, 1, 2].map(k =>
                    frame[0].values[k] * a.values[0] +
                    frame[1].values[k] * a.values[1] +
                    frame[2].values[k] * a.values[2]);
                for (let k = 0; k < 3; ++k) {
                    expectClose(b.values[k], rotated[k], 1e-9, 1e-9);
                }
            }
        });
    });

    it('ignores the velocities and the time', () => {
        // The spring model is position-only; a translation-in-time or a
        // different velocity array must not change the acceleration.
        check(fc.tuple(chainData(2), fc.double({ min: -10, max: 10,
            noNaN: true, noDefaultInfinity: true })), ([d, time]) => {
            const other = d.position.map(() => Vector.fromArray([3, -7]));
            for (let i = 0; i < d.position.length; ++i) {
                const a = d.chain.accelerationAt(i, 0, d.position, d.position);
                const b = d.chain.accelerationAt(i, time, d.position, other);
                expect(b.values).toEqual(a.values);
            }
        });
    });

    it('an immovable particle stays put through an update', () => {
        // setMass(i, 0) gives an infinite mass; ParticleSystem skips those in
        // every stage of the RK4 update.
        check(chainData(2), d => {
            const n = d.position.length;
            d.chain.setMass(0, 0);
            d.chain.update(0);
            expect(d.chain.getPosition(0).values)
                .toEqual(d.position[0].values);
            for (let i = 0; i < n; ++i) {
                for (const x of d.chain.getPosition(i).values) {
                    expect(Number.isFinite(x)).toBe(true);
                }
            }
        });
    });
});
