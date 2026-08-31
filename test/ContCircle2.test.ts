import { describe, it, expect } from 'vitest';
import {
    getContainerCircle2,
    inContainerCircle2,
    mergeContainersCircle2
} from '../src/ContCircle2';
import { Hypersphere, type Circle2 } from '../src/Hypersphere';
import { Vector, length, sub } from '../src/Vector';

function v(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function circle(x: number, y: number, radius: number): Circle2 {
    return Hypersphere.fromCenterRadius(v(x, y), radius);
}

describe('getContainerCircle2', () => {
    it('computes the average-center circle of a square', () => {
        const points = [v(-1, -1), v(1, -1), v(1, 1), v(-1, 1)];
        const c = getContainerCircle2(points);
        expect(c.dimension).toBe(2);
        expect(c.center.values).toEqual([0, 0]);
        expect(c.radius).toBeCloseTo(Math.SQRT2, 14);
    });

    it('degenerates to a zero-radius circle for one point', () => {
        const c = getContainerCircle2([v(3, -4)]);
        expect(c.center.values).toEqual([3, -4]);
        expect(c.radius).toBe(0);
    });

    it('uses the average of the points as the center, not the optimal center', () => {
        // Three points clustered near the origin plus one far away: the
        // average center is pulled toward the cluster, so the circle is
        // larger than the minimum-area circle.
        const points = [v(0, 0), v(0, 0), v(0, 0), v(4, 0)];
        const c = getContainerCircle2(points);
        expect(c.center.values[0]).toBeCloseTo(1, 14);
        expect(c.center.values[1]).toBe(0);
        expect(c.radius).toBeCloseTo(3, 14);
        // The minimum-area circle would have radius 2.
        expect(c.radius).toBeGreaterThan(2);
    });

    it('is invariant to the order of the points', () => {
        const points = [v(1, 2), v(-3, 4), v(0, -5), v(6, 1)];
        const a = getContainerCircle2(points);
        const b = getContainerCircle2([points[2], points[0], points[3], points[1]]);
        expect(a.center.values[0]).toBeCloseTo(b.center.values[0], 14);
        expect(a.center.values[1]).toBeCloseTo(b.center.values[1], 14);
        expect(a.radius).toBeCloseTo(b.radius, 14);
    });

    it('contains every input point, and the radius is attained (randomized)', () => {
        let seed = 987654321;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff - 0.5;
        };

        for (let trial = 0; trial < 50; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 12; ++i) {
                points.push(v(20 * rand(), 20 * rand()));
            }
            const c = getContainerCircle2(points);

            // Center is the arithmetic mean.
            const mean = [0, 0];
            for (const p of points) {
                mean[0] += p.values[0];
                mean[1] += p.values[1];
            }
            expect(c.center.values[0]).toBeCloseTo(mean[0] / points.length, 12);
            expect(c.center.values[1]).toBeCloseTo(mean[1] / points.length, 12);

            // Radius is the largest distance from the center.
            let maxDist = 0;
            for (const p of points) {
                maxDist = Math.max(maxDist, length(sub(p, c.center)));
                expect(length(sub(p, c.center))).toBeLessThanOrEqual(
                    c.radius * (1 + 1e-12));
            }
            expect(c.radius).toBeCloseTo(maxDist, 12);
        }
    });

    it('throws on an empty point set', () => {
        expect(() => getContainerCircle2([]))
            .toThrow('getContainerCircle2: no points.');
    });

    it('throws when the points are not 2D', () => {
        expect(() => getContainerCircle2([Vector.fromArray([0, 0, 0])]))
            .toThrow('getContainerCircle2: points must be 2D.');
        expect(() => getContainerCircle2([v(0, 0), Vector.fromArray([1, 1, 1])]))
            .toThrow('getContainerCircle2: points must be 2D.');
    });
});

describe('inContainerCircle2', () => {
    const c = circle(1, 2, 3);

    it('accepts the center and interior points', () => {
        expect(inContainerCircle2(v(1, 2), c)).toBe(true);
        expect(inContainerCircle2(v(3, 3), c)).toBe(true);
    });

    it('accepts boundary points (the boundary is part of the circle)', () => {
        expect(inContainerCircle2(v(4, 2), c)).toBe(true);
        expect(inContainerCircle2(v(1, -1), c)).toBe(true);
    });

    it('rejects exterior points', () => {
        expect(inContainerCircle2(v(4.0001, 2), c)).toBe(false);
        expect(inContainerCircle2(v(5, 6), c)).toBe(false);
    });

    it('handles a zero-radius circle', () => {
        const point = circle(0, 0, 0);
        expect(inContainerCircle2(v(0, 0), point)).toBe(true);
        expect(inContainerCircle2(v(1e-12, 0), point)).toBe(false);
    });

    it('throws when the inputs are not 2D', () => {
        expect(() => inContainerCircle2(Vector.fromArray([0, 0, 0]), c))
            .toThrow('inContainerCircle2: inputs must be 2D.');
    });
});

describe('mergeContainersCircle2', () => {
    it('merges two separated circles of equal radius', () => {
        const merge = mergeContainersCircle2(circle(0, 0, 1), circle(4, 0, 1));
        expect(merge.center.values[0]).toBeCloseTo(2, 14);
        expect(merge.center.values[1]).toBeCloseTo(0, 14);
        expect(merge.radius).toBeCloseTo(3, 14);
    });

    it('merges two separated circles of different radii', () => {
        // Centers 3 apart with radii 1 and 2: the merged circle spans from
        // (-1,0) to (5,0), so its center is (2,0) and its radius is 3.
        const merge = mergeContainersCircle2(circle(0, 0, 1), circle(3, 0, 2));
        expect(merge.center.values[0]).toBeCloseTo(2, 14);
        expect(merge.center.values[1]).toBeCloseTo(0, 14);
        expect(merge.radius).toBeCloseTo(3, 14);
    });

    it('returns the containing circle when one circle contains the other', () => {
        const outer = circle(0, 0, 10);
        const inner = circle(1, 1, 2);
        expect(mergeContainersCircle2(outer, inner).equals(outer)).toBe(true);
        expect(mergeContainersCircle2(inner, outer).equals(outer)).toBe(true);
    });

    it('returns the larger circle for concentric circles', () => {
        const merge = mergeContainersCircle2(circle(2, 3, 1), circle(2, 3, 5));
        expect(merge.center.values).toEqual([2, 3]);
        expect(merge.radius).toBe(5);
    });

    it('is idempotent for identical circles', () => {
        const c = circle(-1, 4, 2.5);
        expect(mergeContainersCircle2(c, c).equals(c)).toBe(true);
    });

    it('does not alias its inputs when one circle contains the other', () => {
        const outer = circle(0, 0, 10);
        const merge = mergeContainersCircle2(outer, circle(0, 0, 1));
        merge.center.values[0] = 100;
        merge.radius = 0;
        expect(outer.center.values[0]).toBe(0);
        expect(outer.radius).toBe(10);
    });

    it('is commutative and contains both inputs (randomized)', () => {
        let seed = 424242;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };

        for (let trial = 0; trial < 200; ++trial) {
            const c0 = circle(10 * (rand() - 0.5), 10 * (rand() - 0.5), 5 * rand());
            const c1 = circle(10 * (rand() - 0.5), 10 * (rand() - 0.5), 5 * rand());
            const m01 = mergeContainersCircle2(c0, c1);
            const m10 = mergeContainersCircle2(c1, c0);

            expect(m01.center.values[0]).toBeCloseTo(m10.center.values[0], 12);
            expect(m01.center.values[1]).toBeCloseTo(m10.center.values[1], 12);
            expect(m01.radius).toBeCloseTo(m10.radius, 12);

            // Each input circle is inside the merged circle.
            for (const c of [c0, c1]) {
                const d = length(sub(c.center, m01.center));
                expect(d + c.radius).toBeLessThanOrEqual(m01.radius * (1 + 1e-12) + 1e-12);
            }
        }
    });

    it('throws when the inputs are not 2D', () => {
        expect(() => mergeContainersCircle2(circle(0, 0, 1), new Hypersphere(3)))
            .toThrow('mergeContainersCircle2: inputs must be 2D.');
    });
});
