import { describe, expect, it } from 'vitest';
import { DistRayRay, type DistRayRayResult } from '../src/DistRayRay.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

// The minimum of |X0 - X1|^2 over a dense sampling of both rays, with the
// parameters restricted to [0,maxT]. This is an upper bound for the true
// squared distance whenever the true minimum is attained in that range.
function sampledSqrDistance(ray0: Ray, ray1: Ray, maxT: number,
    samples: number): number {
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= samples; ++i) {
        const p = add(ray0.origin, mul(maxT * i / samples, ray0.direction));
        for (let j = 0; j <= samples; ++j) {
            const q = add(ray1.origin,
                mul(maxT * j / samples, ray1.direction));
            const d = sub(p, q);
            best = Math.min(best, dot(d, d));
        }
    }
    return best;
}

describe('DistRayRay', () => {
    const query = new DistRayRay();

    // The reported closest points must lie on the rays (nonnegative
    // parameters) and realize the reported distance.
    function verify(result: DistRayRayResult, ray0: Ray, ray1: Ray): void {
        expect(result.parameter[0]).toBeGreaterThanOrEqual(0);
        expect(result.parameter[1]).toBeGreaterThanOrEqual(0);
        const p = add(ray0.origin, mul(result.parameter[0], ray0.direction));
        const q = add(ray1.origin, mul(result.parameter[1], ray1.direction));
        for (let i = 0; i < p.size; ++i) {
            expect(result.closest[0].values[i]).toBeCloseTo(p.values[i], 10);
            expect(result.closest[1].values[i]).toBeCloseTo(q.values[i], 10);
        }
        const diff = sub(result.closest[0], result.closest[1]);
        expect(result.sqrDistance).toBeCloseTo(dot(diff, diff), 10);
        expect(result.distance).toBeCloseTo(Math.sqrt(result.sqrDistance), 12);
    }

    it('finds interior closest points of skew rays (region 0)', () => {
        // The x-axis ray and a ray along +y starting below the x-axis at
        // x = 2, z = 3.
        const r0 = ray([0, 0, 0], [1, 0, 0]);
        const r1 = ray([2, -1, 3], [0, 1, 0]);
        const result = query.compute(r0, r1);
        expect(result.parameter[0]).toBeCloseTo(2, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.distance).toBeCloseTo(3, 12);
        verify(result, r0, r1);
    });

    it('reports zero distance for intersecting rays', () => {
        const r0 = ray([0, 0, 0], [1, 1, 0]);
        const r1 = ray([4, 0, 0], [-1, 1, 0]);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.closest[0].values[0]).toBeCloseTo(2, 12);
        expect(result.closest[0].values[1]).toBeCloseTo(2, 12);
        verify(result, r0, r1);
    });

    it('clamps ray1 to its origin (region 3)', () => {
        // The perpendicular foot on ray1 has a negative parameter, so the
        // closest ray1 point is its origin.
        const r0 = ray([0, 0, 0], [1, 0, 0]);
        const r1 = ray([2, 1, 3], [0, 1, 0]);
        const result = query.compute(r0, r1);
        expect(result.parameter[1]).toBe(0);
        expect(result.parameter[0]).toBeCloseTo(2, 12);
        expect(result.distance).toBeCloseTo(Math.sqrt(1 + 9), 12);
        verify(result, r0, r1);
    });

    it('clamps ray0 to its origin (region 1)', () => {
        const r0 = ray([0, 0, 0], [1, 0, 0]);
        const r1 = ray([-2, -1, 3], [0, 1, 0]);
        const result = query.compute(r0, r1);
        expect(result.parameter[0]).toBe(0);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 9), 12);
        verify(result, r0, r1);
    });

    it('clamps both rays to their origins (region 2)', () => {
        // Both rays point away from each other.
        const r0 = ray([0, 0, 0], [1, 0, 0]);
        const r1 = ray([-2, 1, 3], [0, 1, 0]);
        const result = query.compute(r0, r1);
        expect(result.parameter[0]).toBe(0);
        expect(result.parameter[1]).toBe(0);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 1 + 9), 12);
        verify(result, r0, r1);
    });

    it('handles parallel rays with the same direction, ray1 ahead', () => {
        const r0 = ray([0, 0, 0], [1, 0, 0]);
        const r1 = ray([5, 3, 0], [2, 0, 0]);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.parameter[1]).toBe(0);
        expect(result.parameter[0]).toBeCloseTo(5, 12);
        verify(result, r0, r1);
    });

    it('handles parallel rays with the same direction, ray1 behind', () => {
        const r0 = ray([0, 0, 0], [1, 0, 0]);
        const r1 = ray([-5, 3, 0], [1, 0, 0]);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.parameter[0]).toBe(0);
        expect(result.parameter[1]).toBeCloseTo(5, 12);
        verify(result, r0, r1);
    });

    it('handles antiparallel rays that overlap', () => {
        const r0 = ray([0, 0, 0], [1, 0, 0]);
        const r1 = ray([5, 3, 0], [-1, 0, 0]);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.parameter[0]).toBeCloseTo(5, 12);
        expect(result.parameter[1]).toBe(0);
        verify(result, r0, r1);
    });

    it('handles antiparallel rays that diverge', () => {
        const r0 = ray([0, 0, 0], [1, 0, 0]);
        const r1 = ray([-5, 3, 0], [-1, 0, 0]);
        const result = query.compute(r0, r1);
        // The closest pair is the two origins.
        expect(result.parameter[0]).toBe(0);
        expect(result.parameter[1]).toBe(0);
        expect(result.distance).toBeCloseTo(Math.sqrt(25 + 9), 12);
        verify(result, r0, r1);
    });

    it('handles collinear rays that overlap (zero distance)', () => {
        const r0 = ray([0, 0, 0], [1, 0, 0]);
        const r1 = ray([4, 0, 0], [-1, 0, 0]);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(0, 12);
        verify(result, r0, r1);
    });

    it('is symmetric in its arguments', () => {
        const cases: Array<[Ray, Ray]> = [
            [ray([0, 0, 0], [1, 0, 0]), ray([2, -1, 3], [0, 1, 0])],
            [ray([0, 0, 0], [1, 0, 0]), ray([-2, 1, 3], [0, 1, 0])],
            [ray([0, 0, 0], [1, 0, 0]), ray([5, 3, 0], [2, 0, 0])],
            [ray([1, 2, 3], [1, 1, 1]), ray([-4, 0, 2], [0, -1, 2])]
        ];
        for (const [r0, r1] of cases) {
            const forward = query.compute(r0, r1);
            const reverse = query.compute(r1, r0);
            expect(reverse.distance).toBeCloseTo(forward.distance, 12);
            verify(forward, r0, r1);
            verify(reverse, r1, r0);
        }
    });

    it('works in 2D', () => {
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([3, 4], [1, 0]);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(4, 12);
        verify(result, r0, r1);
    });

    it('agrees with a dense sampling of both rays', () => {
        let seed = 4242;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 40; ++trial) {
            const r0 = ray([rand(), rand(), rand()],
                [rand(), rand(), rand()]);
            const r1 = ray([rand(), rand(), rand()],
                [rand(), rand(), rand()]);
            const result = query.compute(r0, r1);
            verify(result, r0, r1);

            const best = sampledSqrDistance(r0, r1, 20, 300);
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);

            // The squared distance is convex on the (nonnegative) parameter
            // domain, so a local minimum is the global one. Perturb both
            // parameters and confirm no nearby feasible pair is closer.
            const h = 1e-3;
            for (const d0 of [-h, 0, h]) {
                for (const d1 of [-h, 0, h]) {
                    const t0 = Math.max(result.parameter[0] + d0, 0);
                    const t1 = Math.max(result.parameter[1] + d1, 0);
                    const p = add(r0.origin, mul(t0, r0.direction));
                    const q = add(r1.origin, mul(t1, r1.direction));
                    const d = sub(p, q);
                    expect(dot(d, d)).toBeGreaterThanOrEqual(
                        result.sqrDistance - 1e-8);
                }
            }
        }
    });
});
