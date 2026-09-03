import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistRay3AlignedBox3 } from '../src/DistRay3AlignedBox3.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, length, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(v(...min), v(...max));
}

// The exact squared distance from a point to the solid aligned box, computed
// independently of the library.
function pointBoxSqrDistance(p: Vector, b: AlignedBox): number {
    let sqrDistance = 0;
    for (let i = 0; i < p.size; ++i) {
        const delta = Math.max(0, b.min.values[i] - p.values[i],
            p.values[i] - b.max.values[i]);
        sqrDistance += delta * delta;
    }
    return sqrDistance;
}

// Verify that the reported box point is inside the box.
function verifyBoxPoint(b: AlignedBox, q: Vector): void {
    for (let i = 0; i < q.size; ++i) {
        expect(q.values[i]).toBeGreaterThanOrEqual(b.min.values[i] - 1e-9);
        expect(q.values[i]).toBeLessThanOrEqual(b.max.values[i] + 1e-9);
    }
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

// The distance from a ray point to the solid box is a convex function of the
// ray parameter, so a ternary search over [0,tmax] converges to the true
// minimum independently of the query under test.
function bruteForce(r: Ray, b: AlignedBox, tmax: number): number {
    const f = (t: number): number =>
        pointBoxSqrDistance(add(r.origin, mul(t, r.direction)), b);
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

// Verify that the reported closest points are consistent with the reported
// distance and lie on their primitives.
function verifyClosest(r: Ray, b: AlignedBox,
    result: { distance: number, parameter: number, closest: [Vector, Vector] }
): void {
    expect(result.parameter).toBeGreaterThanOrEqual(0);
    const onRay = add(r.origin, mul(result.parameter, r.direction));
    for (let i = 0; i < onRay.size; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onRay.values[i], 9);
    }
    verifyBoxPoint(b, result.closest[1]);
    expect(length(sub(result.closest[0], result.closest[1])))
        .toBeCloseTo(result.distance, 9);
}

describe('DistRay3AlignedBox3', () => {
    const query = new DistRay3AlignedBox3();
    const unitBox = box([-1, -1, -1], [1, 1, 1]);

    it('returns the ray origin when the ray points away from the box', () => {
        const r = ray([4, 0, 0], [1, 0, 0]);
        const result = query.compute(r, unitBox);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.sqrDistance).toBeCloseTo(9, 12);
        expect(result.parameter).toBe(0);
        expect(result.closest[0].values).toEqual([4, 0, 0]);
        expect(result.closest[1].values).toEqual([1, 0, 0]);
    });

    it('reports zero distance when the ray meets the box', () => {
        const r = ray([-5, 0, 0], [1, 0, 0]);
        const result = query.compute(r, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        verifyClosest(r, unitBox, result);
    });

    it('handles a ray parallel to a box face', () => {
        const r = ray([-4, 3, 0], [1, 0, 0]);
        const result = query.compute(r, unitBox);
        expect(result.distance).toBeCloseTo(2, 12);
        verifyClosest(r, unitBox, result);
    });

    it('finds the corner as the closest point', () => {
        const r = ray([5, 5, 5], [1, 1, 1]);
        const result = query.compute(r, unitBox);
        expect(result.distance).toBeCloseTo(Math.sqrt(3) * 4, 10);
        expect(result.parameter).toBe(0);
        expect(result.closest[1].values).toEqual([1, 1, 1]);
    });

    it('handles a degenerate box that is a single point', () => {
        const b = box([2, 0, 0], [2, 0, 0]);
        const r = ray([0, 0, 0], [0, 1, 0]);
        const result = query.compute(r, b);
        expect(result.distance).toBeCloseTo(2, 10);
        expect(result.parameter).toBe(0);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 13572468;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 150; ++trial) {
            const b = box([-1 - rand(), -1 - rand(), -1 - rand()],
                [1 + rand(), 1 + rand(), 1 + rand()]);
            const r = ray([8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            if (length(r.direction) < 1e-6) {
                continue;
            }
            const result = query.compute(r, b);
            const brute = bruteForce(r, b, 1e6);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-6);
            expect(result.distance).toBeCloseTo(brute, 6);
            verifyClosest(r, b, result);
        }
    });
});
