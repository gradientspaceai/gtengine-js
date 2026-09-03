import { describe, expect, it } from 'vitest';
import { MassSpringSurface } from '../src/MassSpringSurface.js';
import { MassSpringVolume } from '../src/MassSpringVolume.js';
import { MassSpringCurve } from '../src/MassSpringCurve.js';
import { Vector, add, dot, length as vectorLength, sub } from '../src/Vector.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// Exposes the protected members for testing.
class TestMassSpringSurface extends MassSpringSurface {
    accelerationAt(i: number, time: number, position: readonly Vector[],
        velocity: readonly Vector[]): Vector {
        return this.acceleration(i, time, position, velocity);
    }

    indexOf(r: number, c: number): number {
        return this.getIndex(r, c);
    }

    coordinatesOf(i: number): { r: number, c: number } {
        return this.getCoordinates(i);
    }
}

// A gravity-like constant external acceleration.
class GravityMassSpringSurface extends TestMassSpringSurface {
    constructor(dimension: number, numRows: number, numCols: number,
        step: number, private g: Vector) {
        super(dimension, numRows, numCols, step);
    }

    override externalAcceleration(_i: number, _time: number,
        _position: readonly Vector[], _velocity: readonly Vector[]): Vector {
        return this.g.clone();
    }
}

// A uniform lattice of masses at spacing 'spacing' in each direction, with
// every spring of constant 'constant' and resting length 'restLength'.
function makeLattice(numRows: number, numCols: number, spacing: number,
    restLength: number, constant: number, mass: number,
    step: number): TestMassSpringSurface {
    const surface = new TestMassSpringSurface(3, numRows, numCols, step);
    for (let r = 0; r < numRows; ++r) {
        for (let c = 0; c < numCols; ++c) {
            surface.setMassAt(r, c, mass);
            surface.setPositionAt(r, c, v3(c * spacing, r * spacing, 0));
            surface.setVelocityAt(r, c, v3(0, 0, 0));
            surface.setConstantR(r, c, constant);
            surface.setLengthR(r, c, restLength);
            surface.setConstantC(r, c, constant);
            surface.setLengthC(r, c, restLength);
        }
    }
    return surface;
}

describe('MassSpringSurface construction and member access', () => {
    it('reports the lattice dimensions and particle count', () => {
        const surface = new MassSpringSurface(3, 3, 4, 0.01);
        expect(surface.getNumRows()).toBe(3);
        expect(surface.getNumCols()).toBe(4);
        expect(surface.getNumParticles()).toBe(12);
        expect(surface.getDimension()).toBe(3);
        expect(surface.getStep()).toBe(0.01);
    });

    it('uses row-major indexing, index = c + C*r', () => {
        const surface = new TestMassSpringSurface(3, 3, 4, 0.01);
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 4; ++c) {
                const i = surface.indexOf(r, c);
                expect(i).toBe(c + 4 * r);
                expect(surface.coordinatesOf(i)).toEqual({ r, c });
            }
        }
    });

    it('constants and lengths default to zero', () => {
        const surface = new MassSpringSurface(3, 2, 2, 0.1);
        for (let r = 0; r < 2; ++r) {
            for (let c = 0; c < 2; ++c) {
                expect(surface.getConstantR(r, c)).toBe(0);
                expect(surface.getLengthR(r, c)).toBe(0);
                expect(surface.getConstantC(r, c)).toBe(0);
                expect(surface.getLengthC(r, c)).toBe(0);
            }
        }
    });

    it('stores the per-direction constants and lengths independently', () => {
        const surface = new MassSpringSurface(3, 2, 2, 0.1);
        surface.setConstantR(1, 0, 4);
        surface.setLengthR(1, 0, 0.6);
        surface.setConstantC(1, 0, 5);
        surface.setLengthC(1, 0, 0.7);
        expect(surface.getConstantR(1, 0)).toBe(4);
        expect(surface.getLengthR(1, 0)).toBe(0.6);
        expect(surface.getConstantC(1, 0)).toBe(5);
        expect(surface.getLengthC(1, 0)).toBe(0.7);
        // No other cell was touched.
        expect(surface.getConstantR(0, 0)).toBe(0);
        expect(surface.getConstantC(1, 1)).toBe(0);
    });

    it('the lattice accessors agree with the linear-index accessors', () => {
        const surface = new TestMassSpringSurface(3, 2, 3, 0.1);
        surface.setMassAt(1, 2, 2.5);
        surface.setPositionAt(1, 2, v3(1, 2, 3));
        surface.setVelocityAt(1, 2, v3(4, 5, 6));
        const i = surface.indexOf(1, 2);
        expect(surface.getMass(i)).toBe(2.5);
        expect(surface.getMassAt(1, 2)).toBe(2.5);
        expect(surface.getPosition(i).values).toEqual([1, 2, 3]);
        expect(surface.getPositionAt(1, 2).values).toEqual([1, 2, 3]);
        expect(surface.getVelocityAt(1, 2).values).toEqual([4, 5, 6]);
    });

    it('copies the position and velocity inputs', () => {
        const surface = new MassSpringSurface(3, 1, 2, 0.1);
        const p = v3(1, 2, 3);
        surface.setPositionAt(0, 1, p);
        p.set(0, 99);
        expect(surface.getPositionAt(0, 1).values).toEqual([1, 2, 3]);
    });

    it('the default external acceleration is zero', () => {
        const surface = new MassSpringSurface(3, 2, 2, 0.1);
        expect(surface.externalAcceleration(0, 0, [], []).values)
            .toEqual([0, 0, 0]);
    });

    it('supports dimensions other than 3', () => {
        const surface = new MassSpringSurface(2, 2, 2, 0.1);
        expect(surface.getDimension()).toBe(2);
        expect(surface.externalAcceleration(0, 0, [], []).values).toEqual([0, 0]);
    });
});

describe('MassSpringSurface acceleration', () => {
    it('is zero when every spring is at its resting length', () => {
        const spacing = 1.5;
        const surface = makeLattice(3, 4, spacing, spacing, 7, 1, 0.01);
        const position: Vector[] = [];
        for (let i = 0; i < surface.getNumParticles(); ++i) {
            position.push(surface.getPosition(i));
        }
        for (let i = 0; i < surface.getNumParticles(); ++i) {
            const a = surface.accelerationAt(i, 0, position, []);
            expect(vectorLength(a)).toBeLessThan(1e-12);
        }
    });

    it('uses only the springs that exist for corner, edge and interior masses',
        () => {
            // A lattice compressed uniformly to below the resting length: each
            // existing spring pushes the mass away from its neighbor with
            // magnitude constant * (restLength - spacing). Counting the
            // springs at a site gives the acceleration exactly.
            const spacing = 1;
            const restLength = 2;
            const constant = 3;
            const surface = makeLattice(3, 3, spacing, restLength, constant, 1, 0.01);
            const position: Vector[] = [];
            for (let i = 0; i < surface.getNumParticles(); ++i) {
                position.push(surface.getPosition(i));
            }
            const magnitude = constant * (restLength - spacing);

            // Corner (0,0): neighbors at +r and +c, so it is pushed toward
            // (-1,-1) (x is the c axis, y is the r axis).
            const corner = surface.accelerationAt(surface.indexOf(0, 0), 0,
                position, []);
            expect(corner.values[0]).toBeCloseTo(-magnitude, 12);
            expect(corner.values[1]).toBeCloseTo(-magnitude, 12);
            expect(corner.values[2]).toBeCloseTo(0, 12);

            // Opposite corner (2,2): pushed toward (+1,+1).
            const corner22 = surface.accelerationAt(surface.indexOf(2, 2), 0,
                position, []);
            expect(corner22.values[0]).toBeCloseTo(magnitude, 12);
            expect(corner22.values[1]).toBeCloseTo(magnitude, 12);

            // Interior (1,1): every direction is balanced.
            const interior = surface.accelerationAt(surface.indexOf(1, 1), 0,
                position, []);
            expect(vectorLength(interior)).toBeLessThan(1e-12);

            // Edge center (0,1): only the +r spring is unbalanced.
            const edge = surface.accelerationAt(surface.indexOf(0, 1), 0,
                position, []);
            expect(edge.values[0]).toBeCloseTo(0, 12);
            expect(edge.values[1]).toBeCloseTo(-magnitude, 12);

            // Edge center (1,0): only the +c spring is unbalanced.
            const edge10 = surface.accelerationAt(surface.indexOf(1, 0), 0,
                position, []);
            expect(edge10.values[0]).toBeCloseTo(-magnitude, 12);
            expect(edge10.values[1]).toBeCloseTo(0, 12);
        });

    it('uses the spring constants and lengths of the correct cell', () => {
        // The mass at (r,c) owns the springs to (r+1,c) and (r,c+1); the
        // springs to (r-1,c) and (r,c-1) are owned by those neighbors.
        const surface = new TestMassSpringSurface(3, 2, 2, 0.01);
        for (let r = 0; r < 2; ++r) {
            for (let c = 0; c < 2; ++c) {
                surface.setMassAt(r, c, 1);
                surface.setPositionAt(r, c, v3(c, r, 0));
            }
        }
        // Only the spring from (0,0) to (0,1) is active, with rest length 2.
        surface.setConstantC(0, 0, 5);
        surface.setLengthC(0, 0, 2);

        const position = [surface.getPosition(0), surface.getPosition(1),
            surface.getPosition(2), surface.getPosition(3)];
        const a00 = surface.accelerationAt(surface.indexOf(0, 0), 0, position, []);
        const a01 = surface.accelerationAt(surface.indexOf(0, 1), 0, position, []);
        // 5 * (1 - 2/1) * (+1) = -5 for the mass at (0,0).
        expect(a00.values[0]).toBeCloseTo(-5, 12);
        expect(a01.values[0]).toBeCloseTo(5, 12);
        // The masses in row 1 feel nothing.
        const a10 = surface.accelerationAt(surface.indexOf(1, 0), 0, position, []);
        expect(vectorLength(a10)).toBe(0);
    });

    it('divides the spring force by the mass', () => {
        const surface = makeLattice(1, 2, 1, 2, 3, 4, 0.01);
        const position = [surface.getPosition(0), surface.getPosition(1)];
        const a = surface.accelerationAt(0, 0, position, []);
        // constant * (spacing - restLength) / mass = 3 * (1 - 2) / 4.
        expect(a.values[0]).toBeCloseTo(-3 / 4, 12);
    });

    it('adds the external acceleration', () => {
        const g = v3(0, 0, -9.81);
        const surface = new GravityMassSpringSurface(3, 1, 2, 0.01, g);
        surface.setMassAt(0, 0, 1);
        surface.setMassAt(0, 1, 1);
        surface.setPositionAt(0, 0, v3(0, 0, 0));
        surface.setPositionAt(0, 1, v3(1, 0, 0));
        surface.setConstantC(0, 0, 5);
        surface.setLengthC(0, 0, 1);
        const a = surface.accelerationAt(0, 0,
            [surface.getPosition(0), surface.getPosition(1)], []);
        expect(a.values[0]).toBeCloseTo(0, 12);
        expect(a.values[2]).toBeCloseTo(-9.81, 12);
    });

    it('is evaluated at the supplied positions, not the stored ones', () => {
        const surface = makeLattice(1, 2, 1, 1, 4, 1, 0.01);
        // Stretch the pair to length 2 through the argument only.
        const position = [v3(0, 0, 0), v3(2, 0, 0)];
        const a = surface.accelerationAt(0, 0, position, []);
        // 4 * (1 - 1/2) * (+2) = 4 pulling toward the neighbor.
        expect(a.values[0]).toBeCloseTo(4, 12);
        // The stored positions are unchanged.
        expect(surface.getPositionAt(0, 1).values).toEqual([1, 0, 0]);
    });
});

describe('MassSpringSurface dynamics', () => {
    it('keeps a lattice at its resting lengths in equilibrium', () => {
        const spacing = 1.25;
        const surface = makeLattice(3, 3, spacing, spacing, 9, 1.5, 0.005);
        const initial: Vector[] = [];
        for (let i = 0; i < surface.getNumParticles(); ++i) {
            initial.push(surface.getPosition(i).clone());
        }
        let time = 0;
        for (let k = 0; k < 200; ++k, time += 0.005) {
            surface.update(time);
        }
        for (let i = 0; i < surface.getNumParticles(); ++i) {
            expect(vectorLength(sub(surface.getPosition(i), initial[i])))
                .toBeLessThan(1e-10);
        }
    });

    it('expands a uniformly compressed lattice while preserving the center of '
        + 'mass', () => {
            const restLength = 1;
            const surface = makeLattice(2, 2, 0.7, restLength, 8, 1, 0.002);
            let center = new Vector(3);
            for (let i = 0; i < surface.getNumParticles(); ++i) {
                center = add(center, surface.getPosition(i));
            }
            let time = 0;
            for (let k = 0; k < 250; ++k, time += 0.002) {
                surface.update(time);
            }
            let newCenter = new Vector(3);
            for (let i = 0; i < surface.getNumParticles(); ++i) {
                newCenter = add(newCenter, surface.getPosition(i));
            }
            expect(vectorLength(sub(newCenter, center))).toBeLessThan(1e-10);

            // The square stays a square (by symmetry) and has expanded.
            const edge = vectorLength(sub(surface.getPositionAt(0, 1),
                surface.getPositionAt(0, 0)));
            expect(edge).toBeGreaterThan(0.7);
            expect(edge).toBeLessThan(2 * restLength - 0.7 + 1e-9);
            const other = vectorLength(sub(surface.getPositionAt(1, 1),
                surface.getPositionAt(0, 1)));
            expect(other).toBeCloseTo(edge, 9);
        });

    it('holds immovable boundary masses in place while the interior moves',
        () => {
            const restLength = 1;
            const surface = makeLattice(1, 3, 1, restLength, 10, 1, 0.002);
            // Pin the two end masses and displace the middle one.
            surface.setMassAt(0, 0, Number.MAX_VALUE);
            surface.setMassAt(0, 2, Number.MAX_VALUE);
            surface.setPositionAt(0, 1, v3(1, 0.3, 0));
            const end0 = surface.getPositionAt(0, 0).clone();
            const end2 = surface.getPositionAt(0, 2).clone();

            let time = 0;
            let maxOffset = 0;
            for (let k = 0; k < 500; ++k, time += 0.002) {
                surface.update(time);
                maxOffset = Math.max(maxOffset,
                    Math.abs(surface.getPositionAt(0, 1).values[1]));
            }
            expect(surface.getPositionAt(0, 0).values).toEqual(end0.values);
            expect(surface.getPositionAt(0, 2).values).toEqual(end2.values);
            expect(surface.getVelocityAt(0, 0).values).toEqual([0, 0, 0]);
            expect(maxOffset).toBeGreaterThan(0.05);
            expect(surface.getPositionAt(0, 1).values[0])
                .toBeGreaterThan(end0.values[0]);
            expect(surface.getPositionAt(0, 1).values[0])
                .toBeLessThan(end2.values[0]);
        });

    it('holds a pinned boundary row of a sagging cloth in place', () => {
        // A cloth pinned along its top row and pulled down by gravity: the
        // pinned row stays put, the free rows sag in -z.
        const numRows = 4, numCols = 4;
        const step = 1e-3;
        const surface = new GravityMassSpringSurface(3, numRows, numCols, step,
            v3(0, 0, -9.81));
        for (let r = 0; r < numRows; ++r) {
            for (let c = 0; c < numCols; ++c) {
                surface.setMassAt(r, c, r === 0 ? Number.MAX_VALUE : 1);
                surface.setPositionAt(r, c, v3(c, r, 0));
                surface.setConstantR(r, c, 500);
                surface.setLengthR(r, c, 1);
                surface.setConstantC(r, c, 500);
                surface.setLengthC(r, c, 1);
            }
        }
        let time = 0;
        for (let k = 0; k < 400; ++k, time += step) {
            surface.update(time);
        }
        for (let c = 0; c < numCols; ++c) {
            expect(surface.getPositionAt(0, c).values).toEqual([c, 0, 0]);
        }
        // Every free row has moved downward, more so farther from the pins.
        let previous = 0;
        for (let r = 1; r < numRows; ++r) {
            const z = surface.getPositionAt(r, 1).values[2];
            expect(z).toBeLessThan(previous);
            previous = z;
        }
    });

    it('oscillates a pinned/free pair with the analytic period', () => {
        // A 1x2 lattice: mass (0,0) is immovable and mass (0,1) has mass m. In
        // one dimension the restoring force is exactly -k*(d - L), so the
        // motion is simple harmonic with period 2*pi*sqrt(m/k).
        const k = 5;
        const m = 1.25;
        const restLength = 1;
        const amplitude = 0.2;
        const period = 2 * Math.PI * Math.sqrt(m / k);
        const numSteps = 2000;
        const step = period / numSteps;

        const surface = new MassSpringSurface(3, 1, 2, step);
        surface.setMassAt(0, 0, Number.MAX_VALUE);
        surface.setMassAt(0, 1, m);
        surface.setPositionAt(0, 0, v3(0, 0, 0));
        surface.setPositionAt(0, 1, v3(restLength + amplitude, 0, 0));
        surface.setConstantC(0, 0, k);
        surface.setLengthC(0, 0, restLength);

        let time = 0;
        for (let s = 0; s < numSteps / 2; ++s, time += step) {
            surface.update(time);
        }
        expect(surface.getPositionAt(0, 1).values[0])
            .toBeCloseTo(restLength - amplitude, 6);
        for (let s = 0; s < numSteps / 2; ++s, time += step) {
            surface.update(time);
        }
        expect(surface.getPositionAt(0, 1).values[0])
            .toBeCloseTo(restLength + amplitude, 6);
        expect(surface.getVelocityAt(0, 1).values[0]).toBeCloseTo(0, 6);
        expect(surface.getPositionAt(0, 0).values).toEqual([0, 0, 0]);
    });

    it('conserves energy for a free two-mass spring', () => {
        const k = 6;
        const m = 1.5;
        const restLength = 1;
        const step = 1e-3;
        const surface = new MassSpringSurface(3, 1, 2, step);
        surface.setMassAt(0, 0, m);
        surface.setMassAt(0, 1, m);
        surface.setPositionAt(0, 0, v3(0, 0, 0));
        surface.setPositionAt(0, 1, v3(1.4, 0, 0));
        surface.setConstantC(0, 0, k);
        surface.setLengthC(0, 0, restLength);

        const energy = (): number => {
            const d = vectorLength(sub(surface.getPositionAt(0, 1),
                surface.getPositionAt(0, 0))) - restLength;
            let kinetic = 0;
            for (let i = 0; i < 2; ++i) {
                const v = surface.getVelocity(i);
                kinetic += 0.5 * m * dot(v, v);
            }
            return kinetic + 0.5 * k * d * d;
        };

        const initialEnergy = energy();
        let time = 0;
        for (let s = 0; s < 3000; ++s, time += step) {
            surface.update(time);
            expect(Math.abs(energy() - initialEnergy)).toBeLessThan(1e-8);
        }
    });
});

describe('MassSpringSurface cross-checks against the curve and volume', () => {
    it('reproduces MassSpringCurve for a 1xC lattice', () => {
        const numCols = 6;
        const step = 1e-3;
        const constant = 12;
        const restLength = 1;
        const mass = 0.7;

        const curve = new MassSpringCurve(3, numCols, step);
        const surface = new MassSpringSurface(3, 1, numCols, step);
        for (let c = 0; c < numCols; ++c) {
            const p = v3(1.1 * c, 0.05 * Math.sin(c), 0);
            curve.setMass(c, mass);
            curve.setPosition(c, p);
            surface.setMassAt(0, c, mass);
            surface.setPositionAt(0, c, p);
        }
        for (let c = 0; c + 1 < numCols; ++c) {
            curve.setConstant(c, constant);
            curve.setLength(c, restLength);
            surface.setConstantC(0, c, constant);
            surface.setLengthC(0, c, restLength);
        }

        let time = 0;
        for (let s = 0; s < 300; ++s, time += step) {
            curve.update(time);
            surface.update(time);
        }
        for (let c = 0; c < numCols; ++c) {
            expect(vectorLength(sub(surface.getPositionAt(0, c),
                curve.getPosition(c)))).toBeLessThan(1e-12);
        }
    });

    it('reproduces MassSpringVolume for a 1xRxC lattice', () => {
        const numRows = 3, numCols = 4;
        const step = 1e-3;
        const constant = 9;
        const restLength = 1;
        const mass = 1.3;

        const surface = new MassSpringSurface(3, numRows, numCols, step);
        const volume = new MassSpringVolume(3, 1, numRows, numCols, step);
        for (let r = 0; r < numRows; ++r) {
            for (let c = 0; c < numCols; ++c) {
                const p = v3(0.9 * c, 1.05 * r, 0.1 * Math.cos(r + c));
                surface.setMassAt(r, c, mass);
                surface.setPositionAt(r, c, p);
                surface.setConstantR(r, c, constant);
                surface.setLengthR(r, c, restLength);
                surface.setConstantC(r, c, constant);
                surface.setLengthC(r, c, restLength);
                volume.setMassAt(0, r, c, mass);
                volume.setPositionAt(0, r, c, p);
                volume.setConstantR(0, r, c, constant);
                volume.setLengthR(0, r, c, restLength);
                volume.setConstantC(0, r, c, constant);
                volume.setLengthC(0, r, c, restLength);
            }
        }

        let time = 0;
        for (let s = 0; s < 250; ++s, time += step) {
            surface.update(time);
            volume.update(time);
        }
        for (let r = 0; r < numRows; ++r) {
            for (let c = 0; c < numCols; ++c) {
                expect(vectorLength(sub(surface.getPositionAt(r, c),
                    volume.getPositionAt(0, r, c)))).toBeLessThan(1e-12);
            }
        }
    });
});
