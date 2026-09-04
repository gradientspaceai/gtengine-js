import { describe, expect, it } from 'vitest';
import {
    DistRaySegment, type DistRaySegmentResult
} from '../src/DistRaySegment.js';
import { Ray } from '../src/Ray.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { DistLineLine } from '../src/DistLineLine.js';
import { DistPointRay } from '../src/DistPointRay.js';
import { DistPointSegment } from '../src/DistPointSegment.js';
import { Line } from '../src/Line.js';
import { check, expectClose, expectVectorClose, fc, rotationFrame, seededRandom, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

// A near-exact reference minimum: sample the segment parameter densely and,
// for each sample, solve for the optimal ray parameter in closed form
// (clamped to be nonnegative).
function referenceSqrDistance(r: Ray, s: Segment, samples: number): number {
    const segDir = sub(s.p[1], s.p[0]);
    const a00 = dot(r.direction, r.direction);
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= samples; ++i) {
        const q = add(s.p[0], mul(i / samples, segDir));
        const t = Math.max(dot(r.direction, sub(q, r.origin)) / a00, 0);
        const p = add(r.origin, mul(t, r.direction));
        const d = sub(p, q);
        best = Math.min(best, dot(d, d));
    }
    return best;
}

describe('DistRaySegment', () => {
    const query = new DistRaySegment();

    // The reported closest points must lie on the ray (nonnegative
    // parameter) and on the segment (parameter in [0,1]), and realize the
    // reported distance.
    function verify(result: DistRaySegmentResult, r: Ray, s: Segment): void {
        expect(result.parameter[0]).toBeGreaterThanOrEqual(0);
        expect(result.parameter[1]).toBeGreaterThanOrEqual(0);
        expect(result.parameter[1]).toBeLessThanOrEqual(1);
        const p = add(r.origin, mul(result.parameter[0], r.direction));
        const q = add(s.p[0], mul(result.parameter[1], sub(s.p[1], s.p[0])));
        for (let i = 0; i < p.size; ++i) {
            expect(result.closest[0].values[i]).toBeCloseTo(p.values[i], 10);
            expect(result.closest[1].values[i]).toBeCloseTo(q.values[i], 10);
        }
        const diff = sub(result.closest[0], result.closest[1]);
        expect(result.sqrDistance).toBeCloseTo(dot(diff, diff), 10);
        expect(result.distance).toBeCloseTo(Math.sqrt(result.sqrDistance), 12);
    }

    it('finds interior closest points (region 0)', () => {
        const r = ray([0, 0, 0], [1, 0, 0]);
        const s = seg([2, 3, 4], [2, 3, -4]);
        const result = query.compute(r, s);
        expect(result.parameter[0]).toBeCloseTo(2, 12);
        expect(result.parameter[1]).toBeCloseTo(0.5, 12);
        expect(result.distance).toBeCloseTo(3, 12);
        verify(result, r, s);
    });

    it('clamps to the segment endpoint Q1 (region 1)', () => {
        const r = ray([0, 0, 0], [1, 0, 0]);
        const s = seg([2, 3, 5], [2, 3, 1]);
        const result = query.compute(r, s);
        expect(result.parameter[1]).toBe(1);
        expect(result.closest[1].values).toEqual([2, 3, 1]);
        expect(result.distance).toBeCloseTo(Math.sqrt(9 + 1), 12);
        verify(result, r, s);
    });

    it('clamps to the segment endpoint Q0 (region 5)', () => {
        const r = ray([0, 0, 0], [1, 0, 0]);
        const s = seg([2, 3, 1], [2, 3, 5]);
        const result = query.compute(r, s);
        expect(result.parameter[1]).toBe(0);
        expect(result.closest[1].values).toEqual([2, 3, 1]);
        expect(result.distance).toBeCloseTo(Math.sqrt(9 + 1), 12);
        verify(result, r, s);
    });

    it('clamps to the ray origin and a segment interior point', () => {
        // The segment straddles the plane x = -2, behind the ray origin, so
        // the closest ray point is the origin.
        const r = ray([0, 0, 0], [1, 0, 0]);
        const s = seg([-2, 3, -4], [-2, 3, 4]);
        const result = query.compute(r, s);
        expect(result.parameter[0]).toBe(0);
        expect(result.parameter[1]).toBeCloseTo(0.5, 12);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 9), 12);
        verify(result, r, s);
    });

    it('clamps to the ray origin and the segment endpoint Q0', () => {
        const r = ray([0, 0, 0], [1, 0, 0]);
        const s = seg([-2, 3, 1], [-2, 3, 5]);
        const result = query.compute(r, s);
        expect(result.parameter[0]).toBe(0);
        expect(result.parameter[1]).toBe(0);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 9 + 1), 12);
        verify(result, r, s);
    });

    it('clamps to the ray origin and the segment endpoint Q1', () => {
        const r = ray([0, 0, 0], [1, 0, 0]);
        const s = seg([-2, 3, 5], [-2, 3, 1]);
        const result = query.compute(r, s);
        expect(result.parameter[0]).toBe(0);
        expect(result.parameter[1]).toBe(1);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 9 + 1), 12);
        verify(result, r, s);
    });

    it('reports zero distance when the segment crosses the ray', () => {
        const r = ray([0, 0, 0], [1, 0, 0]);
        const s = seg([3, -1, 0], [3, 1, 0]);
        const result = query.compute(r, s);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter[0]).toBeCloseTo(3, 12);
        expect(result.parameter[1]).toBeCloseTo(0.5, 12);
        verify(result, r, s);
    });

    it('handles a parallel segment in front of the ray origin', () => {
        const r = ray([0, 0, 0], [1, 0, 0]);
        const s = seg([2, 4, 0], [7, 4, 0]);
        const result = query.compute(r, s);
        expect(result.distance).toBeCloseTo(4, 12);
        verify(result, r, s);
    });

    it('handles a parallel segment straddling the ray origin', () => {
        const r = ray([0, 0, 0], [1, 0, 0]);
        const s = seg([-3, 4, 0], [3, 4, 0]);
        const result = query.compute(r, s);
        expect(result.distance).toBeCloseTo(4, 12);
        verify(result, r, s);
    });

    it('handles a parallel segment entirely behind the ray origin', () => {
        // Upstream's parallel branch returns a negative ray parameter here,
        // which reports a closest "ray" point off the ray and a distance of
        // 1 instead of sqrt(10). The port clamps the ray parameter to zero.
        const r = ray([0, 0, 0], [1, 0, 0]);
        const forward = query.compute(r, seg([-5, 1, 0], [-3, 1, 0]));
        expect(forward.parameter[0]).toBe(0);
        expect(forward.parameter[1]).toBe(1);
        expect(forward.distance).toBeCloseTo(Math.sqrt(9 + 1), 12);
        verify(forward, r, seg([-5, 1, 0], [-3, 1, 0]));

        // The same segment with the endpoints swapped (antiparallel).
        const reverse = query.compute(r, seg([-3, 1, 0], [-5, 1, 0]));
        expect(reverse.parameter[0]).toBe(0);
        expect(reverse.parameter[1]).toBe(0);
        expect(reverse.distance).toBeCloseTo(Math.sqrt(9 + 1), 12);
        verify(reverse, r, seg([-3, 1, 0], [-5, 1, 0]));
    });

    it('handles a degenerate (zero-length) segment', () => {
        const r = ray([0, 0, 0], [1, 0, 0]);
        const ahead = query.compute(r, seg([2, 5, 0], [2, 5, 0]));
        expect(ahead.distance).toBeCloseTo(5, 12);
        expect(ahead.parameter[0]).toBeCloseTo(2, 12);
        verify(ahead, r, seg([2, 5, 0], [2, 5, 0]));

        // A degenerate segment behind the ray origin: the closest ray point
        // is the origin.
        const behind = query.compute(r, seg([-3, 4, 0], [-3, 4, 0]));
        expect(behind.parameter[0]).toBe(0);
        expect(behind.distance).toBeCloseTo(5, 12);
        verify(behind, r, seg([-3, 4, 0], [-3, 4, 0]));
    });

    it('works in 2D', () => {
        const r = ray([0, 0], [0, 1]);
        const s = seg([3, 4], [3, 9]);
        const result = query.compute(r, s);
        expect(result.distance).toBeCloseTo(3, 12);
        verify(result, r, s);
    });

    it('agrees with a dense sampling of the ray and segment', () => {
        let seed = 777;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 40; ++trial) {
            const r = ray([rand(), rand(), rand()], [rand(), rand(), rand()]);
            const s = seg([rand(), rand(), rand()], [rand(), rand(), rand()]);
            const result = query.compute(r, s);
            verify(result, r, s);

            const best = referenceSqrDistance(r, s, 20000);
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
            expect(result.sqrDistance).toBeGreaterThan(best - 1e-6);

            // The squared distance is convex on the feasible domain, so a
            // local minimum is the global one.
            const segDir = sub(s.p[1], s.p[0]);
            const h = 1e-3;
            for (const d0 of [-h, 0, h]) {
                for (const d1 of [-h, 0, h]) {
                    const t0 = Math.max(result.parameter[0] + d0, 0);
                    const t1 = Math.min(
                        Math.max(result.parameter[1] + d1, 0), 1);
                    const p = add(r.origin, mul(t0, r.direction));
                    const q = add(s.p[0], mul(t1, segDir));
                    const d = sub(p, q);
                    expect(dot(d, d)).toBeGreaterThanOrEqual(
                        result.sqrDistance - 1e-8);
                }
            }
        }
    });

    it('agrees with a dense sampling for parallel configurations', () => {
        // The coordinates are small integers and the segment direction is an
        // integer multiple of the ray direction, so a00*a11 - a01*a01 is
        // exactly zero and the parallel branch is taken. (Constructing the
        // parallel case from arbitrary floating-point coordinates leaves a
        // tiny positive determinant, which sends the query down the
        // nonparallel branch where the upstream formulation is known to lose
        // accuracy; that is what upstream's robust segment/segment query
        // addresses, and DistRaySegment has no robust variant.)
        let seed = 31337;
        const randInt = (lo: number, hi: number) => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return lo + (seed % (hi - lo + 1));
        };
        for (let trial = 0; trial < 60; ++trial) {
            const dir = [randInt(-3, 3), randInt(-3, 3), randInt(-3, 3)];
            if (dir[0] === 0 && dir[1] === 0 && dir[2] === 0) {
                dir[0] = 1;
            }
            const r = ray([randInt(-4, 4), randInt(-4, 4), randInt(-4, 4)],
                dir);
            // The segment direction is an integer multiple (possibly zero,
            // making the segment degenerate) of the ray direction.
            const scale = randInt(-3, 3);
            const q0 = [randInt(-4, 4), randInt(-4, 4), randInt(-4, 4)];
            const s = seg(q0, [q0[0] + scale * dir[0],
                q0[1] + scale * dir[1], q0[2] + scale * dir[2]]);
            const result = query.compute(r, s);
            verify(result, r, s);

            const best = referenceSqrDistance(r, s, 20000);
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
            expect(result.sqrDistance).toBeGreaterThan(best - 1e-6);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of DistRaySegment.ts
// against the upstream header DistRaySegment.h, including regression coverage
// for the ray-parameter clamp the port adds (upstream issue #126).
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

const segment3 = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -8, 8))
    .filter(([p0, p1]) => length(sub(p1, p0)) > 0.25)
    .map(([p0, p1]) => Segment.fromEndpoints(p0, p1));

// Near-parallel pairs make the upstream determinant and its numerators cancel
// to rounding noise, so the reported parameters lose all significance there.
const wellConditioned = fc.tuple(nonUnitRay3, segment3)
    .filter(([ray, seg]) => {
        const sd = sub(seg.p[1], seg.p[0]);
        return length(cross(ray.direction, sd))
            > 0.2 * length(ray.direction) * length(sd);
    });

describe('DistRaySegment verification', () => {
    const query = new DistRaySegment();
    const lineLine = new DistLineLine();
    const pointRay = new DistPointRay();
    const pointSegment = new DistPointSegment();

    it('result is self consistent and both parameters are in their domains',
        () => {
            check(fc.tuple(nonUnitRay3, segment3), ([ray, seg]) => {
                const r = query.compute(ray, seg);
                expect(r.parameter[0]).toBeGreaterThanOrEqual(0);
                expect(r.parameter[1]).toBeGreaterThanOrEqual(0);
                expect(r.parameter[1]).toBeLessThanOrEqual(1);
                expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
                const diff = sub(r.closest[0], r.closest[1]);
                expectClose(r.sqrDistance, dot(diff, diff), 1e-12, 1e-12);
                expectVectorClose(r.closest[0],
                    add(ray.origin, mul(r.parameter[0], ray.direction)), 1e-12,
                    1e-12);
                expectVectorClose(r.closest[1],
                    add(seg.p[0], mul(r.parameter[1], sub(seg.p[1], seg.p[0]))),
                    1e-12, 1e-12);
            });
        });

    it('matches the exact minimum over the domain [0,inf) x [0,1]', () => {
        check(wellConditioned, ([ray, seg]) => {
            // Convex quadratic: the minimum is the unconstrained line/line
            // critical point when it lies in the domain, and otherwise on one
            // of the three boundary faces s0 = 0, s1 = 0 and s1 = 1.
            const unconstrained = lineLine.compute(
                Line.fromOriginDirection(ray.origin, ray.direction),
                Line.fromOriginDirection(seg.p[0], sub(seg.p[1], seg.p[0])));
            let ref = pointSegment.compute(ray.origin, seg).distance;
            ref = Math.min(ref, pointRay.compute(seg.p[0], ray).distance);
            ref = Math.min(ref, pointRay.compute(seg.p[1], ray).distance);
            const [s0, s1] = unconstrained.parameter;
            if (s0 >= 0 && s1 >= 0 && s1 <= 1) {
                ref = Math.min(ref, unconstrained.distance);
            }
            expectClose(query.compute(ray, seg).distance, ref, 1e-8, 1e-8);
        });
    });

    // Restricted to wellConditioned: when the two directions are parallel,
    // det = max(a00*a11 - a01*a01, 0) is a difference of two products that
    // cancel only in exact arithmetic, so it can come out one ulp positive.
    // Upstream then takes the nonparallel branch with numerators that are
    // pure rounding noise and reports a point pair that is on the ray and the
    // segment but far from the minimum. This affects DistLineRay,
    // DistLineSegment, DistRayRay and DistRaySegment alike and is inherent in
    // the upstream parallel test, not a port defect.
    it('is minimal over sampled point pairs', () => {
        const rand = seededRandom(0x51d9);
        check(wellConditioned, ([ray, seg]) => {
            const r = query.compute(ray, seg);
            const sd = sub(seg.p[1], seg.p[0]);
            for (let k = 0; k < 20; ++k) {
                const p = add(ray.origin, mul(12 * rand(), ray.direction));
                const q = add(seg.p[0], mul(rand(), sd));
                const diff = sub(p, q);
                const sqr = dot(diff, diff);
                expect(r.sqrDistance).toBeLessThanOrEqual(sqr + 1e-9 * (1 + sqr));
            }
        }, 60);
    });

    it('is invariant under rigid motions', () => {
        check(fc.tuple(wellConditioned, rotationFrame(3),
            wellScaledVector(3, -5, 5)), ([[ray, seg], R, tr]) => {
            const movedRay = Ray.fromOriginDirection(
                add(rot(R, ray.origin), tr), rot(R, ray.direction));
            const movedSeg = Segment.fromEndpoints(
                add(rot(R, seg.p[0]), tr), add(rot(R, seg.p[1]), tr));
            expectClose(query.compute(ray, seg).distance,
                query.compute(movedRay, movedSeg).distance, 1e-8, 1e-8);
        });
    });

    it('clamps the ray parameter in region 1 (upstream issue #126)', () => {
        // The witness recorded in the issue: a well-conditioned nonparallel
        // configuration (det about 2.79) that reaches region 1, where
        // upstream computes s0 = -(a01+b0)/a00 = -0.377 without clamping and
        // reports a closest "ray" point behind the ray origin together with
        // an under-reported distance (sqrDistance 2.675 instead of 4.136).
        const r = ray([-1.9852731227874756, -2.484701693058014,
            -1.566794514656067],
        [-2.7053112387657166, -0.4477686882019043, -1.6644186973571777]);
        const s = seg([1.403347671031952, -2.980962038040161,
            2.2564467592164874],
        [-1.6737277507781982, -2.9632678627967834, 0.3852156177163124]);
        const result = query.compute(r, s);
        expect(result.parameter[0]).toBe(0);
        expect(result.parameter[1]).toBeGreaterThanOrEqual(0);
        expect(result.parameter[1]).toBeLessThanOrEqual(1);
        // The true minimum is on the face s0 = 0, that is the distance from
        // the ray origin to the segment.
        expectClose(result.sqrDistance,
            pointSegment.compute(r.origin, s).sqrDistance, 1e-12, 1e-12);
        expect(result.sqrDistance).toBeGreaterThan(4);
        expectVectorClose(result.closest[0], r.origin, 0, 0);
    });

    it('never reports a distance below the true minimum', () => {
        // The clamp must not be reachable from a configuration whose true
        // minimum has a positive ray parameter: whenever the port clamps, the
        // face s0 = 0 really does carry the minimum. (Restricted to
        // wellConditioned for the reason given above.)
        check(wellConditioned, ([ray, seg]) => {
            const r = query.compute(ray, seg);
            if (r.parameter[0] === 0) {
                expectClose(r.sqrDistance,
                    pointSegment.compute(ray.origin, seg).sqrDistance, 1e-9,
                    1e-9);
            }
        });
    });
});
