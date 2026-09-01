import { describe, it, expect } from 'vitest';
import { Sector2 } from '../src/Sector2';
import { GTE_C_PI } from '../src/Constants';
import { Vector } from '../src/Vector';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

describe('Sector2 construction', () => {
    it('the default constructor is the unit disk', () => {
        const sector = new Sector2();
        expect(sector.vertex.values).toEqual([0, 0]);
        expect(sector.radius).toBe(1);
        expect(sector.direction.values).toEqual([1, 0]);
        expect(sector.angle).toBe(GTE_C_PI);
        expect(sector.cosAngle).toBe(-1);
        expect(sector.sinAngle).toBe(0);
    });

    it('the factory copies the vectors and sets cos/sin from the angle', () => {
        const vertex = v2(1, 2);
        const direction = v2(0, 1);
        const sector = Sector2.fromVertexRadiusDirectionAngle(vertex, 3,
            direction, Math.PI / 3);
        vertex.set(0, 99);
        direction.set(1, 99);
        expect(sector.vertex.values).toEqual([1, 2]);
        expect(sector.direction.values).toEqual([0, 1]);
        expect(sector.radius).toBe(3);
        expect(sector.angle).toBe(Math.PI / 3);
        expect(sector.cosAngle).toBeCloseTo(0.5, 15);
        expect(sector.sinAngle).toBeCloseTo(Math.sqrt(3) / 2, 15);
    });

    it('rejects vectors that are not 2D', () => {
        expect(() => Sector2.fromVertexRadiusDirectionAngle(
            Vector.fromArray([0, 0, 0]), 1, v2(1, 0), 1)).toThrow();
    });

    it('setAngle keeps angle, cos and sin consistent', () => {
        const sector = new Sector2();
        sector.setAngle(Math.PI / 2);
        expect(sector.angle).toBe(Math.PI / 2);
        expect(sector.cosAngle).toBeCloseTo(0, 15);
        expect(sector.sinAngle).toBe(1);
    });

    it('clone is a deep copy that preserves the stored cos/sin', () => {
        const sector = new Sector2();
        // Deliberately inconsistent state, as C++ copying would preserve.
        sector.angle = 0.25;
        const copy = sector.clone();
        expect(copy.angle).toBe(0.25);
        expect(copy.cosAngle).toBe(-1);
        copy.vertex.set(0, 5);
        expect(sector.vertex.values).toEqual([0, 0]);
    });
});

describe('Sector2 containment', () => {
    it('the default sector is the full unit disk', () => {
        const sector = new Sector2();
        expect(sector.contains(v2(0, 0))).toBe(true);
        expect(sector.contains(v2(0.5, 0.5))).toBe(true);
        expect(sector.contains(v2(-0.9, 0))).toBe(true);
        expect(sector.contains(v2(1, 0))).toBe(true);
        expect(sector.contains(v2(1.0001, 0))).toBe(false);
    });

    it('a quarter-plane sector of angle pi/4 accepts only its wedge', () => {
        // Vertex (0,0), radius 2, direction (1,0), half-angle pi/4.
        const sector = Sector2.fromVertexRadiusDirectionAngle(v2(0, 0), 2,
            v2(1, 0), Math.PI / 4);
        // Inside the wedge and the disk.
        expect(sector.contains(v2(1, 0))).toBe(true);
        expect(sector.contains(v2(1, 0.5))).toBe(true);
        expect(sector.contains(v2(1, -0.5))).toBe(true);
        // Inside the wedge but outside the disk.
        expect(sector.contains(v2(3, 0))).toBe(false);
        // Inside the disk but outside the wedge.
        expect(sector.contains(v2(0.5, 1))).toBe(false);
        expect(sector.contains(v2(-1, 0))).toBe(false);
    });

    it('points exactly on the wedge boundary are contained', () => {
        const angle = Math.PI / 3;
        const sector = Sector2.fromVertexRadiusDirectionAngle(v2(0, 0), 2,
            v2(1, 0), angle);
        const boundary = v2(Math.cos(angle), Math.sin(angle));
        expect(sector.contains(boundary)).toBe(true);
        // Just outside the wedge.
        const outside = v2(Math.cos(angle + 1e-3), Math.sin(angle + 1e-3));
        expect(sector.contains(outside)).toBe(false);
    });

    it('the vertex itself is always contained (length 0)', () => {
        const sector = Sector2.fromVertexRadiusDirectionAngle(v2(4, -3), 1,
            v2(0, 1), 0.1);
        expect(sector.contains(v2(4, -3))).toBe(true);
    });

    it('agrees with the angle/radius criterion over a sampled grid', () => {
        const vertex = v2(1, -2);
        const radius = 2.5;
        const dirAngle = 0.7;
        const halfAngle = 0.9;
        const sector = Sector2.fromVertexRadiusDirectionAngle(vertex, radius,
            v2(Math.cos(dirAngle), Math.sin(dirAngle)), halfAngle);
        for (let i = 0; i < 20; ++i) {
            for (let j = 0; j < 20; ++j) {
                const r = (i + 0.5) * (radius * 1.4 / 20);
                const t = -Math.PI + (j + 0.5) * (2 * Math.PI / 20);
                const p = v2(vertex.get(0) + r * Math.cos(t),
                    vertex.get(1) + r * Math.sin(t));
                let delta = Math.abs(t - dirAngle);
                while (delta > Math.PI) {
                    delta = 2 * Math.PI - delta;
                }
                const expected = r <= radius && delta <= halfAngle;
                expect(sector.contains(p)).toBe(expected);
            }
        }
    });
});

describe('Sector2 comparisons', () => {
    const base = new Sector2();

    it('equals compares vertex, radius, direction and angle', () => {
        expect(base.equals(new Sector2())).toBe(true);
        expect(base.notEquals(new Sector2())).toBe(false);

        const other = base.clone();
        other.setAngle(1);
        expect(base.equals(other)).toBe(false);
        expect(base.notEquals(other)).toBe(true);
    });

    it('equals ignores the derived cos/sin members, as upstream does', () => {
        const other = base.clone();
        other.cosAngle = 0.5;
        other.sinAngle = 0.5;
        expect(base.equals(other)).toBe(true);
    });

    it('lessThan orders by vertex, then radius, then direction, then angle', () => {
        const smallVertex = base.clone();
        smallVertex.vertex = v2(-1, 0);
        expect(smallVertex.lessThan(base)).toBe(true);

        const smallRadius = base.clone();
        smallRadius.radius = 0.5;
        expect(smallRadius.lessThan(base)).toBe(true);

        const smallDirection = base.clone();
        smallDirection.direction = v2(0.5, 0);
        expect(smallDirection.lessThan(base)).toBe(true);

        const smallAngle = base.clone();
        smallAngle.setAngle(1);
        expect(smallAngle.lessThan(base)).toBe(true);
        expect(base.lessThan(smallAngle)).toBe(false);
    });

    it('the derived comparisons are consistent', () => {
        const smaller = base.clone();
        smaller.setAngle(1);
        expect(smaller.lessThanOrEqual(base)).toBe(true);
        expect(base.lessThanOrEqual(base.clone())).toBe(true);
        expect(base.greaterThan(smaller)).toBe(true);
        expect(base.greaterThanOrEqual(smaller)).toBe(true);
        expect(base.greaterThan(base.clone())).toBe(false);
        expect(base.greaterThanOrEqual(base.clone())).toBe(true);
    });
});
