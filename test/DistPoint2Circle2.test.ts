import { describe, expect, it } from 'vitest';
import { DistPoint2Circle2 } from '../src/DistPoint2Circle2';
import { Hypersphere } from '../src/Hypersphere';
import { Vector, dot, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function circle(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(v(...center), radius);
}

describe('DistPoint2Circle2', () => {
    const query = new DistPoint2Circle2();

    it('measures the distance from a point outside the circle', () => {
        const result = query.compute(v(5, 0), circle([0, 0], 2));
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.sqrDistance).toBeCloseTo(9, 12);
        expect(result.equidistant).toBe(false);
        expect(result.closest[1].values[0]).toBeCloseTo(2, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 12);
    });

    it('measures the distance from a point inside the circle', () => {
        const result = query.compute(v(0, 1), circle([0, 0], 4));
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(0, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(4, 12);
    });

    it('reports zero distance for a point on the circle', () => {
        const result = query.compute(v(3, 4), circle([0, 0], 5));
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.equidistant).toBe(false);
    });

    it('flags the circle center as equidistant', () => {
        const result = query.compute(v(2, -3), circle([2, -3], 7));
        expect(result.equidistant).toBe(true);
        expect(result.distance).toBe(7);
        expect(result.sqrDistance).toBe(49);
        // The reported point is C + r*(1,0).
        expect(result.closest[1].values).toEqual([9, -3]);
    });

    it('handles a zero-radius (degenerate) circle', () => {
        const result = query.compute(v(3, 4), circle([0, 0], 0));
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(0, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 12);
    });

    it('places the closest point on the circle at the reported distance',
        () => {
            let seed = 20250;
            const rand = () => {
                seed = (seed * 1103515245 + 12345) % 2147483648;
                return seed / 2147483648 * 10 - 5;
            };
            const c = circle([1, -2], 2.5);
            for (let trial = 0; trial < 100; ++trial) {
                const point = v(rand(), rand());
                const result = query.compute(point, c);

                // The closest point is on the circle.
                const radial = sub(result.closest[1], c.center);
                expect(Math.sqrt(dot(radial, radial))).toBeCloseTo(c.radius,
                    9);

                // It realizes the reported distance.
                const diff = sub(result.closest[0], result.closest[1]);
                expect(Math.sqrt(dot(diff, diff))).toBeCloseTo(
                    result.distance, 9);

                // No sampled circle point is closer.
                for (let k = 0; k < 720; ++k) {
                    const s = k * Math.PI / 360;
                    const q = v(c.center.values[0] + c.radius * Math.cos(s),
                        c.center.values[1] + c.radius * Math.sin(s));
                    const d = sub(point, q);
                    expect(dot(d, d)).toBeGreaterThanOrEqual(
                        result.sqrDistance - 1e-9);
                }
            }
        });
});
