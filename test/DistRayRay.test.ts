import { describe, expect, it } from 'vitest';
import { DistRayRay, type DistRayRayResult } from '../src/DistRayRay.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { DistLineLine } from '../src/DistLineLine.js';
import { DistPointRay } from '../src/DistPointRay.js';
import { Line } from '../src/Line.js';
import { check, expectClose, expectVectorClose, fc, rotationFrame, seededRandom, wellScaledVector } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of DistRayRay.ts
// against the upstream header DistRayRay.h.
// ---------------------------------------------------------------------------

function rot(R: readonly Vector[], p: Vector): Vector {
    let q = mul(p.values[0], R[0]);
    for (let i = 1; i < R.length; ++i) {
        q = add(q, mul(p.values[i], R[i]));
    }
    return q;
}

const nonUnitRay3 = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -3, 3))
    .filter(([, d]) => length(d) > 0.25)
    .map(([o, d]) => Ray.fromOriginDirection(o, d));

// Near-parallel pairs make the upstream determinant and its numerators cancel
// to rounding noise, so the reported parameters lose all significance there.
const wellConditioned = fc.tuple(nonUnitRay3, nonUnitRay3)
    .filter(([r0, r1]) => length(cross(r0.direction, r1.direction))
        > 0.2 * length(r0.direction) * length(r1.direction));

describe('DistRayRay verification', () => {
    const query = new DistRayRay();
    const lineLine = new DistLineLine();
    const pointRay = new DistPointRay();

    it('result is self consistent and both parameters are nonnegative', () => {
        check(fc.tuple(nonUnitRay3, nonUnitRay3), ([r0, r1]) => {
            const r = query.compute(r0, r1);
            expect(r.parameter[0]).toBeGreaterThanOrEqual(0);
            expect(r.parameter[1]).toBeGreaterThanOrEqual(0);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            const diff = sub(r.closest[0], r.closest[1]);
            expectClose(r.sqrDistance, dot(diff, diff), 1e-12, 1e-12);
            expectVectorClose(r.closest[0],
                add(r0.origin, mul(r.parameter[0], r0.direction)), 1e-12,
                1e-12);
            expectVectorClose(r.closest[1],
                add(r1.origin, mul(r.parameter[1], r1.direction)), 1e-12,
                1e-12);
        });
    });

    it('matches the exact minimum over the domain [0,inf)^2', () => {
        check(wellConditioned, ([r0, r1]) => {
            // Convex quadratic on a quadrant: the minimum is the
            // unconstrained line/line critical point when both parameters are
            // nonnegative, and otherwise lies on one of the two boundary
            // faces s0 = 0 or s1 = 0, each a point-ray distance.
            const unconstrained = lineLine.compute(
                Line.fromOriginDirection(r0.origin, r0.direction),
                Line.fromOriginDirection(r1.origin, r1.direction));
            let ref = Math.min(pointRay.compute(r0.origin, r1).distance,
                pointRay.compute(r1.origin, r0).distance);
            if (unconstrained.parameter[0] >= 0
                && unconstrained.parameter[1] >= 0) {
                ref = Math.min(ref, unconstrained.distance);
            }
            expectClose(query.compute(r0, r1).distance, ref, 1e-8, 1e-8);
        });
    });

    it('is symmetric under argument swap', () => {
        check(wellConditioned, ([r0, r1]) => {
            expectClose(query.compute(r0, r1).distance,
                query.compute(r1, r0).distance, 1e-8, 1e-8);
        });
    });

    it('is minimal over sampled point pairs', () => {
        const rand = seededRandom(0x51d8);
        check(wellConditioned, ([r0, r1]) => {
            const r = query.compute(r0, r1);
            for (let k = 0; k < 20; ++k) {
                const p = add(r0.origin, mul(12 * rand(), r0.direction));
                const q = add(r1.origin, mul(12 * rand(), r1.direction));
                const diff = sub(p, q);
                const sqr = dot(diff, diff);
                expect(r.sqrDistance).toBeLessThanOrEqual(sqr + 1e-9 * (1 + sqr));
            }
        }, 60);
    });

    it('is invariant under rigid motions', () => {
        check(fc.tuple(wellConditioned, rotationFrame(3),
            wellScaledVector(3, -5, 5)), ([[r0, r1], R, tr]) => {
            const move = (ray: Ray): Ray => Ray.fromOriginDirection(
                add(rot(R, ray.origin), tr), rot(R, ray.direction));
            expectClose(query.compute(r0, r1).distance,
                query.compute(move(r0), move(r1)).distance, 1e-8, 1e-8);
        });
    });

    it('handles exactly parallel rays', () => {
        // A power-of-two scale keeps the determinant exactly zero, so the
        // parallel branch is exercised deterministically.
        check(fc.tuple(nonUnitRay3, wellScaledVector(3, -8, 8),
            fc.constantFrom(1, -1, 2, -2, 0.5, -0.5, 4, -4)),
        ([r0, o1, k]) => {
            const r1 = Ray.fromOriginDirection(o1, mul(k, r0.direction));
            const r = query.compute(r0, r1);
            expect(r.parameter[0]).toBeGreaterThanOrEqual(0);
            expect(r.parameter[1]).toBeGreaterThanOrEqual(0);
            // The parallel branch always places one closest point at a ray
            // origin, so the answer is one of the two point-ray distances.
            const ref = Math.min(pointRay.compute(r0.origin, r1).distance,
                pointRay.compute(r1.origin, r0).distance);
            expectClose(r.distance, ref, 1e-8, 1e-8);
        });
    });
});
