import { describe, it, expect } from 'vitest';
import {
    getContainerSphere3,
    inContainerSphere3,
    mergeContainersSphere3
} from '../src/ContSphere3';
import { Hypersphere, type Sphere3 } from '../src/Hypersphere';
import { Vector, length, sub } from '../src/Vector';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function sphere(x: number, y: number, z: number, radius: number): Sphere3 {
    return Hypersphere.fromCenterRadius(v(x, y, z), radius);
}

describe('getContainerSphere3', () => {
    it('computes the average-center sphere of a cube', () => {
        const points: Vector[] = [];
        for (const x of [-1, 1]) {
            for (const y of [-1, 1]) {
                for (const z of [-1, 1]) {
                    points.push(v(x, y, z));
                }
            }
        }
        const s = getContainerSphere3(points);
        expect(s.dimension).toBe(3);
        expect(s.center.values).toEqual([0, 0, 0]);
        expect(s.radius).toBeCloseTo(Math.sqrt(3), 14);
    });

    it('computes the average-center sphere of the unit octahedron vertices', () => {
        const points = [
            v(1, 0, 0), v(-1, 0, 0), v(0, 1, 0),
            v(0, -1, 0), v(0, 0, 1), v(0, 0, -1)
        ];
        const s = getContainerSphere3(points);
        expect(s.center.values[0]).toBeCloseTo(0, 14);
        expect(s.center.values[1]).toBeCloseTo(0, 14);
        expect(s.center.values[2]).toBeCloseTo(0, 14);
        expect(s.radius).toBeCloseTo(1, 14);
    });

    it('degenerates to a zero-radius sphere for one point', () => {
        const s = getContainerSphere3([v(3, -4, 5)]);
        expect(s.center.values).toEqual([3, -4, 5]);
        expect(s.radius).toBe(0);
    });

    it('uses the average of the points as the center, not the optimal center', () => {
        const points = [v(0, 0, 0), v(0, 0, 0), v(0, 0, 0), v(4, 0, 0)];
        const s = getContainerSphere3(points);
        expect(s.center.values[0]).toBeCloseTo(1, 14);
        expect(s.radius).toBeCloseTo(3, 14);
        // The minimum-volume sphere would have radius 2.
        expect(s.radius).toBeGreaterThan(2);
    });

    it('contains every input point, and the radius is attained (randomized)', () => {
        let seed = 20240813;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff - 0.5;
        };

        for (let trial = 0; trial < 50; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 12; ++i) {
                points.push(v(20 * rand(), 20 * rand(), 20 * rand()));
            }
            const s = getContainerSphere3(points);

            const mean = [0, 0, 0];
            for (const p of points) {
                for (let d = 0; d < 3; ++d) {
                    mean[d] += p.values[d];
                }
            }
            for (let d = 0; d < 3; ++d) {
                expect(s.center.values[d]).toBeCloseTo(mean[d] / points.length, 12);
            }

            let maxDist = 0;
            for (const p of points) {
                const dist = length(sub(p, s.center));
                maxDist = Math.max(maxDist, dist);
                expect(inContainerSphere3(p, s)).toBe(true);
            }
            expect(s.radius).toBeCloseTo(maxDist, 12);
        }
    });

    it('throws on an empty point set', () => {
        expect(() => getContainerSphere3([]))
            .toThrow('getContainerSphere3: no points.');
    });

    it('throws when the points are not 3D', () => {
        expect(() => getContainerSphere3([Vector.fromArray([0, 0])]))
            .toThrow('getContainerSphere3: points must be 3D.');
        expect(() => getContainerSphere3([v(0, 0, 0), Vector.fromArray([1, 1])]))
            .toThrow('getContainerSphere3: points must be 3D.');
    });
});

describe('inContainerSphere3', () => {
    const s = sphere(1, 2, -1, 3);

    it('accepts the center and interior points', () => {
        expect(inContainerSphere3(v(1, 2, -1), s)).toBe(true);
        expect(inContainerSphere3(v(2, 3, 0), s)).toBe(true);
    });

    it('accepts boundary points (the boundary is part of the sphere)', () => {
        expect(inContainerSphere3(v(4, 2, -1), s)).toBe(true);
        expect(inContainerSphere3(v(1, 2, 2), s)).toBe(true);
    });

    it('rejects exterior points', () => {
        expect(inContainerSphere3(v(4.0001, 2, -1), s)).toBe(false);
        expect(inContainerSphere3(v(5, 6, 7), s)).toBe(false);
    });

    it('handles a zero-radius sphere', () => {
        const point = sphere(0, 0, 0, 0);
        expect(inContainerSphere3(v(0, 0, 0), point)).toBe(true);
        expect(inContainerSphere3(v(0, 1e-12, 0), point)).toBe(false);
    });

    it('throws when the inputs are not 3D', () => {
        expect(() => inContainerSphere3(Vector.fromArray([0, 0]), s))
            .toThrow('inContainerSphere3: inputs must be 3D.');
    });
});

describe('mergeContainersSphere3', () => {
    it('merges two separated spheres of equal radius', () => {
        const merge = mergeContainersSphere3(
            sphere(0, 0, 0, 1), sphere(0, 4, 0, 1));
        expect(merge.center.values[0]).toBeCloseTo(0, 14);
        expect(merge.center.values[1]).toBeCloseTo(2, 14);
        expect(merge.center.values[2]).toBeCloseTo(0, 14);
        expect(merge.radius).toBeCloseTo(3, 14);
    });

    it('merges two separated spheres of different radii', () => {
        // Centers 3 apart with radii 1 and 2: the merged sphere spans from
        // (-1,0,0) to (5,0,0), so its center is (2,0,0) and its radius is 3.
        const merge = mergeContainersSphere3(
            sphere(0, 0, 0, 1), sphere(3, 0, 0, 2));
        expect(merge.center.values[0]).toBeCloseTo(2, 14);
        expect(merge.radius).toBeCloseTo(3, 14);
    });

    it('returns the containing sphere when one sphere contains the other', () => {
        const outer = sphere(0, 0, 0, 10);
        const inner = sphere(1, 1, 1, 2);
        expect(mergeContainersSphere3(outer, inner).equals(outer)).toBe(true);
        expect(mergeContainersSphere3(inner, outer).equals(outer)).toBe(true);
    });

    it('returns the larger sphere for concentric spheres', () => {
        const merge = mergeContainersSphere3(
            sphere(2, 3, 4, 1), sphere(2, 3, 4, 5));
        expect(merge.center.values).toEqual([2, 3, 4]);
        expect(merge.radius).toBe(5);
    });

    it('is idempotent for identical spheres', () => {
        const s = sphere(-1, 4, 2, 2.5);
        expect(mergeContainersSphere3(s, s).equals(s)).toBe(true);
    });

    it('does not alias its inputs when one sphere contains the other', () => {
        const outer = sphere(0, 0, 0, 10);
        const merge = mergeContainersSphere3(outer, sphere(0, 0, 0, 1));
        merge.center.values[0] = 100;
        merge.radius = 0;
        expect(outer.center.values[0]).toBe(0);
        expect(outer.radius).toBe(10);
    });

    it('is commutative and contains both inputs (randomized)', () => {
        let seed = 55555;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };

        for (let trial = 0; trial < 200; ++trial) {
            const s0 = sphere(10 * (rand() - 0.5), 10 * (rand() - 0.5),
                10 * (rand() - 0.5), 5 * rand());
            const s1 = sphere(10 * (rand() - 0.5), 10 * (rand() - 0.5),
                10 * (rand() - 0.5), 5 * rand());
            const m01 = mergeContainersSphere3(s0, s1);
            const m10 = mergeContainersSphere3(s1, s0);

            for (let d = 0; d < 3; ++d) {
                expect(m01.center.values[d]).toBeCloseTo(m10.center.values[d], 12);
            }
            expect(m01.radius).toBeCloseTo(m10.radius, 12);

            for (const s of [s0, s1]) {
                const dist = length(sub(s.center, m01.center));
                expect(dist + s.radius).toBeLessThanOrEqual(
                    m01.radius * (1 + 1e-12) + 1e-12);
            }
        }
    });

    it('throws when the inputs are not 3D', () => {
        expect(() => mergeContainersSphere3(sphere(0, 0, 0, 1), new Hypersphere(2)))
            .toThrow('mergeContainersSphere3: inputs must be 3D.');
    });
});
