import { describe, expect, it } from 'vitest';
import { DistRay3Rectangle3 } from '../src/DistRay3Rectangle3';
import { Ray } from '../src/Ray';
import { Rectangle } from '../src/Rectangle';
import { Vector, add, dot, length, mul, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

function rectangle(center: number[], axis: number[][],
    extent: number[]): Rectangle {
    return Rectangle.fromCenterAxisExtent(v(...center),
        axis.map(a => v(...a)), v(...extent));
}

// The exact squared distance from a point to the solid rectangle, computed
// independently of the library (the axes are orthonormal, so the rectangle
// coordinates are clamped independently).
function pointRectangleSqrDistance(p: Vector, r: Rectangle): number {
    const delta = sub(p, r.center);
    let closest = r.center.clone();
    let s0 = dot(delta, r.axis[0]);
    let s1 = dot(delta, r.axis[1]);
    s0 = Math.max(-r.extent.values[0], Math.min(r.extent.values[0], s0));
    s1 = Math.max(-r.extent.values[1], Math.min(r.extent.values[1], s1));
    closest = add(closest, add(mul(s0, r.axis[0]), mul(s1, r.axis[1])));
    const diff = sub(p, closest);
    return dot(diff, diff);
}

// The distance from a ray point to the solid rectangle is a convex function
// of the ray parameter, so a ternary search over [0,tmax] converges to the
// true minimum independently of the query under test.
function bruteForce(r: Ray, rect: Rectangle, tmax: number): number {
    const f = (t: number): number => pointRectangleSqrDistance(
        add(r.origin, mul(t, r.direction)), rect);
    let lo = 0;
    let hi = tmax;
    for (let i = 0; i < 200; ++i) {
        const m0 = lo + (hi - lo) / 3;
        const m1 = hi - (hi - lo) / 3;
        if (f(m0) < f(m1)) {
            hi = m1;
        }
        else {
            lo = m0;
        }
    }
    return Math.sqrt(f(0.5 * (lo + hi)));
}

function verifyClosest(r: Ray, rect: Rectangle,
    result: {
        distance: number, parameter: number,
        cartesian: [number, number], closest: [Vector, Vector]
    }): void {
    expect(result.parameter).toBeGreaterThanOrEqual(0);
    const onRay = add(r.origin, mul(result.parameter, r.direction));
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onRay.values[i], 9);
    }
    // The reported cartesian coordinates must be in range and reproduce the
    // reported rectangle point.
    expect(Math.abs(result.cartesian[0]))
        .toBeLessThanOrEqual(rect.extent.values[0] + 1e-9);
    expect(Math.abs(result.cartesian[1]))
        .toBeLessThanOrEqual(rect.extent.values[1] + 1e-9);
    const onRect = add(rect.center,
        add(mul(result.cartesian[0], rect.axis[0]),
            mul(result.cartesian[1], rect.axis[1])));
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[1].values[i]).toBeCloseTo(onRect.values[i], 9);
    }
    expect(length(sub(result.closest[0], result.closest[1])))
        .toBeCloseTo(result.distance, 9);
}

describe('DistRay3Rectangle3', () => {
    const query = new DistRay3Rectangle3();
    const unitRect = rectangle([0, 0, 0], [[1, 0, 0], [0, 1, 0]], [1, 1]);

    it('returns the ray origin when the ray points away', () => {
        const r = ray([0, 0, 5], [0, 0, 1]);
        const result = query.compute(r, unitRect);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.sqrDistance).toBeCloseTo(25, 12);
        expect(result.parameter).toBe(0);
        expect(result.cartesian).toEqual([0, 0]);
        expect(result.closest[1].values).toEqual([0, 0, 0]);
    });

    it('reports zero distance when the ray hits the rectangle', () => {
        const r = ray([0.25, -0.5, 5], [0, 0, -1]);
        const result = query.compute(r, unitRect);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.cartesian[0]).toBeCloseTo(0.25, 10);
        expect(result.cartesian[1]).toBeCloseTo(-0.5, 10);
    });

    it('finds the rectangle edge when the ray passes beside it', () => {
        // The ray is parallel to the rectangle plane at height 4 and moves
        // away from the rectangle in x.
        const r = ray([3, 0, 4], [1, 0, 0]);
        const result = query.compute(r, unitRect);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 16), 10);
        expect(result.parameter).toBe(0);
        expect(result.cartesian[0]).toBeCloseTo(1, 10);
        verifyClosest(r, unitRect, result);
    });

    it('handles a ray whose line meets the rectangle behind the origin',
        () => {
            const r = ray([0, 0, 2], [0, 0, 1]);
            const result = query.compute(r, unitRect);
            expect(result.distance).toBeCloseTo(2, 12);
            expect(result.parameter).toBe(0);
        });

    it('handles a degenerate rectangle with zero extents', () => {
        const rect = rectangle([0, 0, 0], [[1, 0, 0], [0, 1, 0]], [0, 0]);
        const r = ray([0, 0, 3], [0, 0, 1]);
        const result = query.compute(r, rect);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values).toEqual([0, 0, 0]);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 20240817;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 150; ++trial) {
            const angle = 2 * Math.PI * rand();
            const c = Math.cos(angle);
            const s = Math.sin(angle);
            const rect = rectangle([2 * rand() - 1, 2 * rand() - 1,
                2 * rand() - 1], [[c, s, 0], [-s, c, 0]],
            [0.2 + rand(), 0.2 + rand()]);
            const r = ray([8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            if (length(r.direction) < 1e-6) {
                continue;
            }
            const result = query.compute(r, rect);
            const brute = bruteForce(r, rect, 1e6);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-6);
            expect(result.distance).toBeCloseTo(brute, 6);
            verifyClosest(r, rect, result);
        }
    });
});
