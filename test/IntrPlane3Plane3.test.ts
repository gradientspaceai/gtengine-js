import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane.js';
import {
    IntrPlane3Plane3TI,
    IntrPlane3Plane3FI,
    defaultIntrPlane3Plane3FIResult
} from '../src/IntrPlane3Plane3.js';
import { Vector, add, dot, mul, normalize } from '../src/Vector.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function plane(normal: number[], constant: number): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalConstant(n, constant);
}

describe('IntrPlane3Plane3', () => {
    const ti = new IntrPlane3Plane3TI();
    const fi = new IntrPlane3Plane3FI();

    it('has the documented default result', () => {
        const result = defaultIntrPlane3Plane3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.isLine).toBe(false);
        expect(result.line.origin.values).toEqual([0, 0, 0]);
        expect(result.plane.normal.values).toEqual([0, 0, 0]);
    });

    it('intersects two orthogonal planes in the expected line', () => {
        const p0 = plane([1, 0, 0], 0);
        const p1 = plane([0, 1, 0], 0);
        expect(ti.test(p0, p1).intersect).toBe(true);

        const result = fi.find(p0, p1);
        expect(result.intersect).toBe(true);
        expect(result.isLine).toBe(true);
        expect(result.line.origin.values[0]).toBeCloseTo(0, 12);
        expect(result.line.origin.values[1]).toBeCloseTo(0, 12);
        expect(result.line.origin.values[2]).toBeCloseTo(0, 12);
        // Cross((1,0,0),(0,1,0)) = (0,0,1).
        expect(result.line.direction.values[0]).toBeCloseTo(0, 12);
        expect(result.line.direction.values[1]).toBeCloseTo(0, 12);
        expect(result.line.direction.values[2]).toBeCloseTo(1, 12);
    });

    it('offsets the line origin by the plane constants', () => {
        const p0 = plane([1, 0, 0], 2);
        const p1 = plane([0, 1, 0], 3);
        const result = fi.find(p0, p1);
        expect(result.isLine).toBe(true);
        expect(result.line.origin.values[0]).toBeCloseTo(2, 12);
        expect(result.line.origin.values[1]).toBeCloseTo(3, 12);
        expect(result.line.origin.values[2]).toBeCloseTo(0, 12);
    });

    it('reports coplanar planes with the same normal direction', () => {
        const p0 = plane([0, 0, 1], 4);
        const p1 = plane([0, 0, 1], 4);
        expect(ti.test(p0, p1).intersect).toBe(true);
        const result = fi.find(p0, p1);
        expect(result.intersect).toBe(true);
        expect(result.isLine).toBe(false);
        expect(result.plane.constant).toBe(4);
        expect(result.plane.normal.values).toEqual([0, 0, 1]);
    });

    it('reports coplanar planes with opposite normal directions', () => {
        const p0 = plane([0, 0, 1], 4);
        const p1 = plane([0, 0, -1], -4);
        expect(ti.test(p0, p1).intersect).toBe(true);
        const result = fi.find(p0, p1);
        expect(result.intersect).toBe(true);
        expect(result.isLine).toBe(false);
    });

    it('reports parallel but distinct planes as disjoint', () => {
        const p0 = plane([0, 0, 1], 4);
        const p1 = plane([0, 0, 1], 5);
        expect(ti.test(p0, p1).intersect).toBe(false);
        const result = fi.find(p0, p1);
        expect(result.intersect).toBe(false);
        expect(result.isLine).toBe(false);

        const p2 = plane([0, 0, -1], 5);
        expect(ti.test(p0, p2).intersect).toBe(false);
        expect(fi.find(p0, p2).intersect).toBe(false);
    });

    it('rejects non-3D planes', () => {
        const p2 = Hyperplane.fromNormalConstant(
            Vector.fromArray([1, 0]), 0);
        expect(() => ti.test(p2, p2)).toThrow();
        expect(() => fi.find(p2, p2)).toThrow();
    });

    it('produces a line lying in both planes for random pairs', () => {
        let seed = 777333;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 300; ++trial) {
            const p0 = plane([rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                rand() * 6 - 3);
            const p1 = plane([rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                rand() * 6 - 3);
            const result = fi.find(p0, p1);
            expect(ti.test(p0, p1).intersect).toBe(result.intersect);
            if (!result.isLine) {
                continue;
            }

            // The direction is unit length and orthogonal to both normals.
            expect(dot(result.line.direction, result.line.direction))
                .toBeCloseTo(1, 9);
            expect(dot(result.line.direction, p0.normal)).toBeCloseTo(0, 9);
            expect(dot(result.line.direction, p1.normal)).toBeCloseTo(0, 9);

            // Several points of the line satisfy both plane equations.
            for (const t of [-3, 0, 2.5]) {
                const X = add(result.line.origin, mul(t, result.line.direction));
                expect(dot(p0.normal, X) - p0.constant).toBeCloseTo(0, 8);
                expect(dot(p1.normal, X) - p1.constant).toBeCloseTo(0, 8);
            }
        }
    });
});
