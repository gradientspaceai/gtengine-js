import { describe, expect, it } from 'vitest';
import { DistPointRectangle } from '../src/DistPointRectangle';
import { Rectangle } from '../src/Rectangle';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

// A rectangle in 3D whose plane is z = 0, centered at the origin, with the
// standard axes and the given extents.
function unitRect(e0: number, e1: number): Rectangle {
    return Rectangle.fromCenterAxisExtent(v(0, 0, 0),
        [v(1, 0, 0), v(0, 1, 0)], v(e0, e1));
}

describe('DistPointRectangle', () => {
    const query = new DistPointRectangle();

    it('measures the plane offset for a point over the rectangle', () => {
        const result = query.compute(v(0.5, -0.25, 3), unitRect(1, 2));
        expect(result.cartesian).toEqual([0.5, -0.25]);
        expect(result.closest[1].values).toEqual([0.5, -0.25, 0]);
        expect(result.distance).toBeCloseTo(3, 12);
    });

    it('reports zero distance for a point on the rectangle', () => {
        const result = query.compute(v(-1, 2, 0), unitRect(1, 2));
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.cartesian).toEqual([-1, 2]);
    });

    it('clamps the coordinates outside the extents', () => {
        const result = query.compute(v(5, -7, 0), unitRect(1, 2));
        expect(result.cartesian).toEqual([1, -2]);
        expect(result.closest[1].values).toEqual([1, -2, 0]);
        expect(result.distance).toBeCloseTo(Math.sqrt(16 + 25), 12);
    });

    it('clamps a single coordinate (edge region)', () => {
        const result = query.compute(v(0.5, 9, 4), unitRect(1, 2));
        expect(result.cartesian).toEqual([0.5, 2]);
        expect(result.distance).toBeCloseTo(Math.sqrt(49 + 16), 12);
    });

    it('handles a rotated rectangle in 2D', () => {
        const a0 = v(1, 1);
        const a1 = v(-1, 1);
        normalize(a0);
        normalize(a1);
        const rect = Rectangle.fromCenterAxisExtent(v(2, 3), [a0, a1],
            v(1, 1));
        // The point is 5 units along a0 from the center, so it is clamped to
        // the extent 1 and the distance is 4.
        const point = add(rect.center, mul(5, a0));
        const result = query.compute(point, rect);
        expect(result.cartesian[0]).toBeCloseTo(1, 12);
        expect(result.cartesian[1]).toBeCloseTo(0, 12);
        expect(result.distance).toBeCloseTo(4, 12);
    });

    it('agrees with a sampled minimum over the rectangle', () => {
        let seed = 777;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 8 - 4;
        };
        const rect = unitRect(1.5, 0.75);
        for (let trial = 0; trial < 40; ++trial) {
            const point = v(rand(), rand(), rand());
            const result = query.compute(point, rect);

            expect(Math.abs(result.cartesian[0])).toBeLessThanOrEqual(1.5);
            expect(Math.abs(result.cartesian[1])).toBeLessThanOrEqual(0.75);

            let best = Number.MAX_VALUE;
            for (let i = 0; i <= 120; ++i) {
                const s0 = -1.5 + 3 * i / 120;
                for (let j = 0; j <= 120; ++j) {
                    const s1 = -0.75 + 1.5 * j / 120;
                    const q = add(rect.center,
                        add(mul(s0, rect.axis[0]), mul(s1, rect.axis[1])));
                    const diff = sub(point, q);
                    best = Math.min(best, dot(diff, diff));
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
        }
    });
});
