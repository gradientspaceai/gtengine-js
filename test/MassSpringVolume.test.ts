import { describe, expect, it } from 'vitest';
import { MassSpringVolume } from '../src/MassSpringVolume.js';
import { MassSpringCurve } from '../src/MassSpringCurve.js';
import { Vector, add, dot, length as vectorLength, sub } from '../src/Vector.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// Exposes the protected members for testing.
class TestMassSpringVolume extends MassSpringVolume {
    accelerationAt(i: number, time: number, position: readonly Vector[],
        velocity: readonly Vector[]): Vector {
        return this.acceleration(i, time, position, velocity);
    }

    indexOf(s: number, r: number, c: number): number {
        return this.getIndex(s, r, c);
    }

    coordinatesOf(i: number): { s: number, r: number, c: number } {
        return this.getCoordinates(i);
    }
}

// A gravity-like constant external acceleration.
class GravityMassSpringVolume extends TestMassSpringVolume {
    constructor(dimension: number, numSlices: number, numRows: number,
        numCols: number, step: number, private g: Vector) {
        super(dimension, numSlices, numRows, numCols, step);
    }

    override externalAcceleration(_i: number, _time: number,
        _position: readonly Vector[], _velocity: readonly Vector[]): Vector {
        return this.g.clone();
    }
}

// A uniform lattice of masses at spacing 'spacing' in each direction, with
// every spring of constant 'constant' and resting length 'restLength'.
function makeLattice(numSlices: number, numRows: number, numCols: number,
    spacing: number, restLength: number, constant: number, mass: number,
    step: number): TestMassSpringVolume {
    const volume = new TestMassSpringVolume(3, numSlices, numRows, numCols, step);
    for (let s = 0; s < numSlices; ++s) {
        for (let r = 0; r < numRows; ++r) {
            for (let c = 0; c < numCols; ++c) {
                volume.setMassAt(s, r, c, mass);
                volume.setPositionAt(s, r, c,
                    v3(c * spacing, r * spacing, s * spacing));
                volume.setVelocityAt(s, r, c, v3(0, 0, 0));
                volume.setConstantS(s, r, c, constant);
                volume.setLengthS(s, r, c, restLength);
                volume.setConstantR(s, r, c, constant);
                volume.setLengthR(s, r, c, restLength);
                volume.setConstantC(s, r, c, constant);
                volume.setLengthC(s, r, c, restLength);
            }
        }
    }
    return volume;
}

describe('MassSpringVolume construction and member access', () => {
    it('reports the lattice dimensions and particle count', () => {
        const volume = new MassSpringVolume(3, 2, 3, 4, 0.01);
        expect(volume.getNumSlices()).toBe(2);
        expect(volume.getNumRows()).toBe(3);
        expect(volume.getNumCols()).toBe(4);
        expect(volume.getNumParticles()).toBe(24);
        expect(volume.getDimension()).toBe(3);
        expect(volume.getStep()).toBe(0.01);
    });

    it('uses lexicographical indexing, index = c + C*(r + R*s)', () => {
        const volume = new TestMassSpringVolume(3, 2, 3, 4, 0.01);
        for (let s = 0; s < 2; ++s) {
            for (let r = 0; r < 3; ++r) {
                for (let c = 0; c < 4; ++c) {
                    const i = volume.indexOf(s, r, c);
                    expect(i).toBe(c + 4 * (r + 3 * s));
                    expect(volume.coordinatesOf(i)).toEqual({ s, r, c });
                }
            }
        }
    });

    it('constants and lengths default to zero', () => {
        const volume = new MassSpringVolume(3, 2, 2, 2, 0.1);
        for (let s = 0; s < 2; ++s) {
            for (let r = 0; r < 2; ++r) {
                for (let c = 0; c < 2; ++c) {
                    expect(volume.getConstantS(s, r, c)).toBe(0);
                    expect(volume.getLengthS(s, r, c)).toBe(0);
                    expect(volume.getConstantR(s, r, c)).toBe(0);
                    expect(volume.getLengthR(s, r, c)).toBe(0);
                    expect(volume.getConstantC(s, r, c)).toBe(0);
                    expect(volume.getLengthC(s, r, c)).toBe(0);
                }
            }
        }
    });

    it('stores the per-direction constants and lengths independently', () => {
        const volume = new MassSpringVolume(3, 2, 2, 2, 0.1);
        volume.setConstantS(1, 0, 1, 3);
        volume.setLengthS(1, 0, 1, 0.5);
        volume.setConstantR(1, 0, 1, 4);
        volume.setLengthR(1, 0, 1, 0.6);
        volume.setConstantC(1, 0, 1, 5);
        volume.setLengthC(1, 0, 1, 0.7);
        expect(volume.getConstantS(1, 0, 1)).toBe(3);
        expect(volume.getLengthS(1, 0, 1)).toBe(0.5);
        expect(volume.getConstantR(1, 0, 1)).toBe(4);
        expect(volume.getLengthR(1, 0, 1)).toBe(0.6);
        expect(volume.getConstantC(1, 0, 1)).toBe(5);
        expect(volume.getLengthC(1, 0, 1)).toBe(0.7);
        // No other cell was touched.
        expect(volume.getConstantS(0, 0, 0)).toBe(0);
        expect(volume.getConstantC(1, 1, 1)).toBe(0);
    });

    it('the lattice accessors agree with the linear-index accessors', () => {
        const volume = new TestMassSpringVolume(3, 2, 2, 3, 0.1);
        volume.setMassAt(1, 1, 2, 2.5);
        volume.setPositionAt(1, 1, 2, v3(1, 2, 3));
        volume.setVelocityAt(1, 1, 2, v3(4, 5, 6));
        const i = volume.indexOf(1, 1, 2);
        expect(volume.getMass(i)).toBe(2.5);
        expect(volume.getMassAt(1, 1, 2)).toBe(2.5);
        expect(volume.getPosition(i).values).toEqual([1, 2, 3]);
        expect(volume.getPositionAt(1, 1, 2).values).toEqual([1, 2, 3]);
        expect(volume.getVelocityAt(1, 1, 2).values).toEqual([4, 5, 6]);
    });

    it('copies the position and velocity inputs', () => {
        const volume = new MassSpringVolume(3, 1, 1, 2, 0.1);
        const p = v3(1, 2, 3);
        volume.setPositionAt(0, 0, 1, p);
        p.set(0, 99);
        expect(volume.getPositionAt(0, 0, 1).values).toEqual([1, 2, 3]);
    });

    it('the default external acceleration is zero', () => {
        const volume = new MassSpringVolume(3, 2, 2, 2, 0.1);
        expect(volume.externalAcceleration(0, 0, [], []).values).toEqual([0, 0, 0]);
    });
});

describe('MassSpringVolume acceleration', () => {
    it('is zero when every spring is at its resting length', () => {
        const spacing = 1.5;
        const volume = makeLattice(3, 3, 3, spacing, spacing, 7, 1, 0.01);
        const position: Vector[] = [];
        for (let i = 0; i < volume.getNumParticles(); ++i) {
            position.push(volume.getPosition(i));
        }
        for (let i = 0; i < volume.getNumParticles(); ++i) {
            const a = volume.accelerationAt(i, 0, position, []);
            expect(vectorLength(a)).toBeLessThan(1e-12);
        }
    });

    it('uses only the springs that exist for corner, edge, face and interior '
        + 'masses', () => {
            // A lattice compressed uniformly by half the resting length: each
            // existing spring pushes the mass away from its neighbor with
            // magnitude constant * (spacing - restLength). Counting the
            // springs at a site gives the acceleration exactly.
            const spacing = 1;
            const restLength = 2;
            const constant = 3;
            const volume = makeLattice(3, 3, 3, spacing, restLength, constant,
                1, 0.01);
            const position: Vector[] = [];
            for (let i = 0; i < volume.getNumParticles(); ++i) {
                position.push(volume.getPosition(i));
            }
            // Each spring contributes constant*(1 - restLength/spacing)*diff
            // = 3*(1-2)*(+/-1) = -/+3 along its axis; the push is away from
            // the neighbor.
            const magnitude = constant * (restLength - spacing);

            // Corner (0,0,0): neighbors at +s, +r, +c, so it is pushed toward
            // (-1,-1,-1) (x is the c axis, y is the r axis, z is the s axis).
            const corner = volume.accelerationAt(volume.indexOf(0, 0, 0), 0,
                position, []);
            expect(corner.values[0]).toBeCloseTo(-magnitude, 12);
            expect(corner.values[1]).toBeCloseTo(-magnitude, 12);
            expect(corner.values[2]).toBeCloseTo(-magnitude, 12);

            // Interior (1,1,1): every direction is balanced.
            const interior = volume.accelerationAt(volume.indexOf(1, 1, 1), 0,
                position, []);
            expect(vectorLength(interior)).toBeLessThan(1e-12);

            // Face center (0,1,1): only the +s spring is unbalanced.
            const face = volume.accelerationAt(volume.indexOf(0, 1, 1), 0,
                position, []);
            expect(face.values[0]).toBeCloseTo(0, 12);
            expect(face.values[1]).toBeCloseTo(0, 12);
            expect(face.values[2]).toBeCloseTo(-magnitude, 12);

            // Edge center (0, 0, 1): the +s and +r springs are unbalanced.
            const edge = volume.accelerationAt(volume.indexOf(0, 0, 1), 0,
                position, []);
            expect(edge.values[0]).toBeCloseTo(0, 12);
            expect(edge.values[1]).toBeCloseTo(-magnitude, 12);
            expect(edge.values[2]).toBeCloseTo(-magnitude, 12);
        });

    it('divides the spring force by the mass', () => {
        const volume = makeLattice(1, 1, 2, 1, 2, 3, 4, 0.01);
        const position = [volume.getPosition(0), volume.getPosition(1)];
        const a = volume.accelerationAt(0, 0, position, []);
        // constant * (spacing - restLength) / mass = 3 * (1 - 2) / 4.
        expect(a.values[0]).toBeCloseTo(-3 / 4, 12);
    });

    it('adds the external acceleration', () => {
        const g = v3(0, 0, -9.81);
        const volume = new GravityMassSpringVolume(3, 1, 1, 2, 0.01, g);
        volume.setMassAt(0, 0, 0, 1);
        volume.setMassAt(0, 0, 1, 1);
        volume.setPositionAt(0, 0, 0, v3(0, 0, 0));
        volume.setPositionAt(0, 0, 1, v3(1, 0, 0));
        volume.setConstantC(0, 0, 0, 5);
        volume.setLengthC(0, 0, 0, 1);
        const a = volume.accelerationAt(0, 0,
            [volume.getPosition(0), volume.getPosition(1)], []);
        expect(a.values[0]).toBeCloseTo(0, 12);
        expect(a.values[2]).toBeCloseTo(-9.81, 12);
    });
});

describe('MassSpringVolume dynamics', () => {
    it('keeps a lattice at its resting lengths in equilibrium', () => {
        const spacing = 1.25;
        const volume = makeLattice(2, 3, 3, spacing, spacing, 9, 1.5, 0.005);
        const initial: Vector[] = [];
        for (let i = 0; i < volume.getNumParticles(); ++i) {
            initial.push(volume.getPosition(i).clone());
        }
        let time = 0;
        for (let k = 0; k < 200; ++k, time += 0.005) {
            volume.update(time);
        }
        for (let i = 0; i < volume.getNumParticles(); ++i) {
            expect(vectorLength(sub(volume.getPosition(i), initial[i])))
                .toBeLessThan(1e-10);
        }
    });

    it('expands a uniformly compressed lattice while preserving the center of '
        + 'mass', () => {
            const restLength = 1;
            const volume = makeLattice(2, 2, 2, 0.7, restLength, 8, 1, 0.002);
            let center = new Vector(3);
            for (let i = 0; i < volume.getNumParticles(); ++i) {
                center = add(center, volume.getPosition(i));
            }
            let time = 0;
            for (let k = 0; k < 250; ++k, time += 0.002) {
                volume.update(time);
            }
            let newCenter = new Vector(3);
            for (let i = 0; i < volume.getNumParticles(); ++i) {
                newCenter = add(newCenter, volume.getPosition(i));
            }
            expect(vectorLength(sub(newCenter, center))).toBeLessThan(1e-10);

            // The cube stays a cube (by symmetry) and has expanded.
            const edge = vectorLength(sub(volume.getPositionAt(0, 0, 1),
                volume.getPositionAt(0, 0, 0)));
            expect(edge).toBeGreaterThan(0.7);
            expect(edge).toBeLessThan(2 * restLength - 0.7 + 1e-9);
            for (const [a, b] of [[[0, 0, 0], [0, 1, 0]], [[0, 0, 0], [1, 0, 0]],
                [[1, 1, 0], [1, 1, 1]]] as number[][][]) {
                const d = vectorLength(sub(
                    volume.getPositionAt(a[0], a[1], a[2]),
                    volume.getPositionAt(b[0], b[1], b[2])));
                expect(d).toBeCloseTo(edge, 9);
            }
        });

    it('holds immovable boundary masses in place while the interior moves',
        () => {
            const restLength = 1;
            const volume = makeLattice(1, 1, 3, 1, restLength, 10, 1, 0.002);
            // Pin the two end masses and displace the middle one.
            volume.setMassAt(0, 0, 0, Number.MAX_VALUE);
            volume.setMassAt(0, 0, 2, Number.MAX_VALUE);
            volume.setPositionAt(0, 0, 1, v3(1, 0.3, 0));
            const end0 = volume.getPositionAt(0, 0, 0).clone();
            const end2 = volume.getPositionAt(0, 0, 2).clone();

            let time = 0;
            let maxOffset = 0;
            for (let k = 0; k < 500; ++k, time += 0.002) {
                volume.update(time);
                maxOffset = Math.max(maxOffset,
                    Math.abs(volume.getPositionAt(0, 0, 1).values[1]));
            }
            expect(volume.getPositionAt(0, 0, 0).values).toEqual(end0.values);
            expect(volume.getPositionAt(0, 0, 2).values).toEqual(end2.values);
            expect(volume.getVelocityAt(0, 0, 0).values).toEqual([0, 0, 0]);
            expect(maxOffset).toBeGreaterThan(0.05);
            // The middle mass stays between the pinned ends.
            expect(volume.getPositionAt(0, 0, 1).values[0])
                .toBeGreaterThan(end0.values[0]);
            expect(volume.getPositionAt(0, 0, 1).values[0])
                .toBeLessThan(end2.values[0]);
        });

    it('oscillates a pinned/free pair with the analytic period', () => {
        // A 1x1x2 lattice: mass (0,0,0) is immovable and mass (0,0,1) has
        // mass m. In one dimension the restoring force is exactly
        // -k*(d - L), so the motion is simple harmonic with period
        // 2*pi*sqrt(m/k).
        const k = 5;
        const m = 1.25;
        const restLength = 1;
        const amplitude = 0.2;
        const period = 2 * Math.PI * Math.sqrt(m / k);
        const numSteps = 2000;
        const step = period / numSteps;

        const volume = new MassSpringVolume(3, 1, 1, 2, step);
        volume.setMassAt(0, 0, 0, Number.MAX_VALUE);
        volume.setMassAt(0, 0, 1, m);
        volume.setPositionAt(0, 0, 0, v3(0, 0, 0));
        volume.setPositionAt(0, 0, 1, v3(restLength + amplitude, 0, 0));
        volume.setConstantC(0, 0, 0, k);
        volume.setLengthC(0, 0, 0, restLength);

        let time = 0;
        for (let s = 0; s < numSteps / 2; ++s, time += step) {
            volume.update(time);
        }
        expect(volume.getPositionAt(0, 0, 1).values[0])
            .toBeCloseTo(restLength - amplitude, 6);
        for (let s = 0; s < numSteps / 2; ++s, time += step) {
            volume.update(time);
        }
        expect(volume.getPositionAt(0, 0, 1).values[0])
            .toBeCloseTo(restLength + amplitude, 6);
        expect(volume.getVelocityAt(0, 0, 1).values[0]).toBeCloseTo(0, 6);
        expect(volume.getPositionAt(0, 0, 0).values).toEqual([0, 0, 0]);
    });

    it('conserves energy for a free two-mass spring', () => {
        const k = 6;
        const m = 1.5;
        const restLength = 1;
        const step = 1e-3;
        const volume = new MassSpringVolume(3, 1, 1, 2, step);
        volume.setMassAt(0, 0, 0, m);
        volume.setMassAt(0, 0, 1, m);
        volume.setPositionAt(0, 0, 0, v3(0, 0, 0));
        volume.setPositionAt(0, 0, 1, v3(1.4, 0, 0));
        volume.setConstantC(0, 0, 0, k);
        volume.setLengthC(0, 0, 0, restLength);

        const energy = (): number => {
            const d = vectorLength(sub(volume.getPositionAt(0, 0, 1),
                volume.getPositionAt(0, 0, 0))) - restLength;
            let kinetic = 0;
            for (let i = 0; i < 2; ++i) {
                const v = volume.getVelocity(i);
                kinetic += 0.5 * m * dot(v, v);
            }
            return kinetic + 0.5 * k * d * d;
        };

        const initialEnergy = energy();
        let time = 0;
        for (let s = 0; s < 3000; ++s, time += step) {
            volume.update(time);
            expect(Math.abs(energy() - initialEnergy)).toBeLessThan(1e-8);
        }
    });

    it('reproduces MassSpringCurve for a 1x1xC lattice', () => {
        const numCols = 6;
        const step = 1e-3;
        const constant = 12;
        const restLength = 1;
        const mass = 0.7;

        const curve = new MassSpringCurve(3, numCols, step);
        const volume = new MassSpringVolume(3, 1, 1, numCols, step);
        for (let c = 0; c < numCols; ++c) {
            const p = v3(1.1 * c, 0.05 * Math.sin(c), 0);
            curve.setMass(c, mass);
            curve.setPosition(c, p);
            volume.setMassAt(0, 0, c, mass);
            volume.setPositionAt(0, 0, c, p);
        }
        for (let c = 0; c + 1 < numCols; ++c) {
            curve.setConstant(c, constant);
            curve.setLength(c, restLength);
            volume.setConstantC(0, 0, c, constant);
            volume.setLengthC(0, 0, c, restLength);
        }

        let time = 0;
        for (let s = 0; s < 200; ++s, time += step) {
            curve.update(time);
            volume.update(time);
        }
        for (let c = 0; c < numCols; ++c) {
            expect(vectorLength(sub(curve.getPosition(c),
                volume.getPositionAt(0, 0, c)))).toBeLessThan(1e-12);
            expect(vectorLength(sub(curve.getVelocity(c),
                volume.getVelocityAt(0, 0, c)))).toBeLessThan(1e-12);
        }
    });
});
