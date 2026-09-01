import { describe, it, expect } from 'vitest';
import { Parallelepiped3 } from '../src/Parallelepiped3';
import { Vector } from '../src/Vector';
import { dotCross } from '../src/Vector3';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

describe('Parallelepiped3 construction', () => {
    it('the default constructor is the [-1,1]^3 cube', () => {
        const p = new Parallelepiped3();
        expect(p.center.values).toEqual([0, 0, 0]);
        expect(p.axis.length).toBe(3);
        expect(p.axis[0].values).toEqual([1, 0, 0]);
        expect(p.axis[1].values).toEqual([0, 1, 0]);
        expect(p.axis[2].values).toEqual([0, 0, 1]);
    });

    it('fromCenterAxis copies the inputs', () => {
        const center = v3(1, 2, 3);
        const axis = [v3(2, 0, 0), v3(0, 3, 0), v3(0, 0, 4)];
        const p = Parallelepiped3.fromCenterAxis(center, axis);
        center.set(0, 99);
        axis[2].set(2, 99);
        expect(p.center.values).toEqual([1, 2, 3]);
        expect(p.axis[2].values).toEqual([0, 0, 4]);
    });

    it('requires a right-handed basis', () => {
        // Swapping two axes makes the basis left-handed.
        expect(() => Parallelepiped3.fromCenterAxis(v3(0, 0, 0),
            [v3(0, 1, 0), v3(1, 0, 0), v3(0, 0, 1)])).toThrow(
            'The axes must form a right-handed basis.');
        // Coplanar axes are degenerate: DotCross = 0.
        expect(() => Parallelepiped3.fromCenterAxis(v3(0, 0, 0),
            [v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, 0)])).toThrow();
        // A right-handed non-orthogonal basis is accepted.
        const sheared = [v3(1, 0, 0), v3(1, 1, 0), v3(1, 1, 1)];
        expect(dotCross(sheared[0], sheared[1], sheared[2])).toBe(1);
        expect(() => Parallelepiped3.fromCenterAxis(v3(0, 0, 0),
            sheared)).not.toThrow();
    });

    it('rejects vectors that are not 3D', () => {
        expect(() => Parallelepiped3.fromCenterAxis(Vector.fromArray([0, 0]),
            [v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)])).toThrow();
    });

    it('clone is a deep copy', () => {
        const p = new Parallelepiped3();
        const copy = p.clone();
        copy.center.set(0, 5);
        copy.axis[2].set(2, 7);
        expect(p.center.values).toEqual([0, 0, 0]);
        expect(p.axis[2].values).toEqual([0, 0, 1]);
    });
});

describe('Parallelepiped3 getVertices', () => {
    it('gives the bit-pattern-ordered corners of the default cube', () => {
        const vertices = new Parallelepiped3().getVertices();
        expect(vertices.length).toBe(8);
        expect(vertices[0].values).toEqual([-1, -1, -1]);
        expect(vertices[1].values).toEqual([1, -1, -1]);
        expect(vertices[2].values).toEqual([-1, 1, -1]);
        expect(vertices[3].values).toEqual([1, 1, -1]);
        expect(vertices[4].values).toEqual([-1, -1, 1]);
        expect(vertices[5].values).toEqual([1, -1, 1]);
        expect(vertices[6].values).toEqual([-1, 1, 1]);
        expect(vertices[7].values).toEqual([1, 1, 1]);
    });

    it('handles a scaled, translated box', () => {
        const p = Parallelepiped3.fromCenterAxis(v3(1, 2, 3),
            [v3(2, 0, 0), v3(0, 3, 0), v3(0, 0, 4)]);
        const vertices = p.getVertices();
        expect(vertices[0].values).toEqual([-1, -1, -1]);
        expect(vertices[7].values).toEqual([3, 5, 7]);
        expect(vertices[3].values).toEqual([3, 5, -1]);
        expect(vertices[4].values).toEqual([-1, -1, 7]);
    });

    it('the vertices are the center plus the eight sign combinations', () => {
        const p = Parallelepiped3.fromCenterAxis(v3(-3, 7, 2),
            [v3(1, 0, 0), v3(1, 1, 0), v3(1, 1, 1)]);
        const vertices = p.getVertices();
        for (let i = 0; i < 8; ++i) {
            const s = [2 * (i & 1) - 1, 2 * ((i >> 1) & 1) - 1,
                2 * ((i >> 2) & 1) - 1];
            for (let c = 0; c < 3; ++c) {
                const expected = p.center.get(c) + s[0] * p.axis[0].get(c)
                    + s[1] * p.axis[1].get(c) + s[2] * p.axis[2].get(c);
                expect(vertices[i].get(c)).toBe(expected);
            }
        }
    });

    it('the volume is 8*DotCross(A0,A1,A2)', () => {
        // For the axis-aligned case the box is 4-by-6-by-8 = 192, and
        // 8 * DotCross((2,0,0),(0,3,0),(0,0,4)) = 8 * 24 = 192.
        const axis = [v3(2, 0, 0), v3(0, 3, 0), v3(0, 0, 4)];
        expect(8 * dotCross(axis[0], axis[1], axis[2])).toBe(192);
        const p = Parallelepiped3.fromCenterAxis(v3(1, 2, 3), axis);
        const vertices = p.getVertices();
        // Extents along each axis direction are twice the axis lengths.
        expect(vertices[1].get(0) - vertices[0].get(0)).toBe(4);
        expect(vertices[2].get(1) - vertices[0].get(1)).toBe(6);
        expect(vertices[4].get(2) - vertices[0].get(2)).toBe(8);
    });
});

describe('Parallelepiped3 comparisons', () => {
    const base = new Parallelepiped3();

    it('equals compares the center and all three axes', () => {
        expect(base.equals(new Parallelepiped3())).toBe(true);
        expect(base.notEquals(new Parallelepiped3())).toBe(false);

        const other = base.clone();
        other.axis[2] = v3(0, 0, 2);
        expect(base.equals(other)).toBe(false);
        expect(base.notEquals(other)).toBe(true);
    });

    it('lessThan orders by center, then the axis array lexicographically', () => {
        const smallCenter = base.clone();
        smallCenter.center = v3(-1, 0, 0);
        expect(smallCenter.lessThan(base)).toBe(true);

        const smallAxis0 = base.clone();
        smallAxis0.axis[0] = v3(0.5, 0, 0);
        expect(smallAxis0.lessThan(base)).toBe(true);

        const smallAxis2 = base.clone();
        smallAxis2.axis[2] = v3(0, 0, 0.5);
        expect(smallAxis2.lessThan(base)).toBe(true);
        expect(base.lessThan(smallAxis2)).toBe(false);
    });

    it('the derived comparisons are consistent', () => {
        const bigger = base.clone();
        bigger.axis[2] = v3(0, 0, 2);
        expect(base.lessThanOrEqual(bigger)).toBe(true);
        expect(base.lessThanOrEqual(base.clone())).toBe(true);
        expect(bigger.greaterThan(base)).toBe(true);
        expect(bigger.greaterThanOrEqual(base)).toBe(true);
        expect(base.greaterThan(base.clone())).toBe(false);
        expect(base.greaterThanOrEqual(base.clone())).toBe(true);
    });
});
