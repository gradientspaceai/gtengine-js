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
