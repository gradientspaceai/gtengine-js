import { describe, expect, it } from 'vitest';
import { DistLineRay } from '../src/DistLineRay.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { DistLineLine } from '../src/DistLineLine.js';
import { DistLineSegment } from '../src/DistLineSegment.js';
import { DistPointLine } from '../src/DistPointLine.js';
import { Segment } from '../src/Segment.js';
import { check, expectClose, expectVectorClose, fc, rotationFrame, seededRandom, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

describe('DistLineRay', () => {
    const query = new DistLineRay();

    it('uses interior points when the ray parameter is nonnegative', () => {
        // The x-axis line and the ray (0,3,0)+s*(0,-1,1). A ray point is
        // (0,3-s,s), whose distance to the x-axis is sqrt((3-s)^2+s^2). The
        // minimum is at s = 1.5, giving the point (0,1.5,1.5) and the
        // distance sqrt(4.5).
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            ray([0, 3, 0], [0, -1, 1]));
        expect(result.parameter[1]).toBeCloseTo(1.5, 10);
        expect(result.distance).toBeCloseTo(Math.sqrt(4.5), 10);
        expect(result.closest[1].values[1]).toBeCloseTo(1.5, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(1.5, 10);
    });

    it('clamps to the ray origin when the unconstrained s1 is negative',
        () => {
            // The ray points away from the line.
            const result = query.compute(line([0, 0, 0], [1, 0, 0]),
                ray([0, 3, 0], [0, 1, 0]));
            expect(result.parameter[1]).toBe(0);
            expect(result.distance).toBeCloseTo(3, 12);
            expect(result.closest[1].values).toEqual([0, 3, 0]);
        });

    it('handles a line and ray that are parallel', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            ray([5, 4, 0], [2, 0, 0]));
        expect(result.parameter[1]).toBe(0);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values).toEqual([5, 4, 0]);
    });

    it('reports zero distance when the ray meets the line', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            ray([2, 3, 0], [0, -1, 0]));
        expect(result.distance).toBeCloseTo(0, 10);
        expect(result.parameter[1]).toBeCloseTo(3, 10);
    });

    it('agrees with a sampled minimum', () => {
        let seed = 606;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 40; ++trial) {
            const l = line([rand(), rand(), rand()],
                [rand() + 4, rand(), rand()]);
            const r = ray([rand(), rand(), rand()],
                [rand(), rand() + 4, rand()]);
            const result = query.compute(l, r);

            expect(result.parameter[1]).toBeGreaterThanOrEqual(0);

            // For each sampled ray point, the exact line distance is at least
            // the reported distance.
            const a00 = dot(l.direction, l.direction);
            let best = Number.MAX_VALUE;
            for (let k = 0; k <= 1500; ++k) {
                const s1 = k * 0.01;
                const q = add(r.origin, mul(s1, r.direction));
                const w = sub(q, l.origin);
                const t = dot(l.direction, w) / a00;
                const d = sub(w, mul(t, l.direction));
                best = Math.min(best, dot(d, d));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of DistLineRay.ts
// against the upstream header DistLineRay.h.
// ---------------------------------------------------------------------------

function rot(R: readonly Vector[], p: Vector): Vector {
    let q = mul(p.values[0], R[0]);
    for (let i = 1; i < R.length; ++i) {
        q = add(q, mul(p.values[i], R[i]));
    }
    return q;
}

const nonUnitLine3 = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -3, 3))
    .filter(([, d]) => length(d) > 0.25)
    .map(([o, d]) => Line.fromOriginDirection(o, d));

const nonUnitRay3 = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -3, 3))
    .filter(([, d]) => length(d) > 0.25)
    .map(([o, d]) => Ray.fromOriginDirection(o, d));

// Pairs that are well away from parallel; near-parallel configurations make
// the upstream determinant and its numerators cancel to rounding noise, so
// the parameters (not the algorithm) lose all significance there.
const wellConditioned = fc.tuple(nonUnitLine3, nonUnitRay3)
    .filter(([l, ray]) => length(cross(l.direction, ray.direction))
        > 0.2 * length(l.direction) * length(ray.direction));

describe('DistLineRay verification', () => {
    const query = new DistLineRay();
    const lineLine = new DistLineLine();
    const pointLine = new DistPointLine();
    const lineSegment = new DistLineSegment();

    it('result is self consistent and the ray parameter is nonnegative', () => {
        check(fc.tuple(nonUnitLine3, nonUnitRay3), ([l, ray]) => {
            const r = query.compute(l, ray);
            expect(r.parameter[1]).toBeGreaterThanOrEqual(0);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            const diff = sub(r.closest[0], r.closest[1]);
            expectClose(r.sqrDistance, dot(diff, diff), 1e-12, 1e-12);
            expectVectorClose(r.closest[0],
                add(l.origin, mul(r.parameter[0], l.direction)), 1e-12, 1e-12);
            expectVectorClose(r.closest[1],
                add(ray.origin, mul(r.parameter[1], ray.direction)), 1e-12,
                1e-12);
        });
    });

    it('matches the exact minimum over the domain R x [0,inf)', () => {
        check(wellConditioned, ([l, ray]) => {
            // The squared distance is a convex quadratic on the domain. Its
            // minimum is either the unconstrained line/line critical point,
            // when that point has a nonnegative ray parameter, or the value
            // on the boundary face s1 = 0, which is a point-line distance.
            const unconstrained = lineLine.compute(l,
                Line.fromOriginDirection(ray.origin, ray.direction));
            let ref = pointLine.compute(ray.origin, l).distance;
            if (unconstrained.parameter[1] >= 0) {
                ref = Math.min(ref, unconstrained.distance);
            }
            expectClose(query.compute(l, ray).distance, ref, 1e-8, 1e-8);
        });
    });

    it('is minimal over sampled point pairs', () => {
        const rand = seededRandom(0x51d6);
        check(wellConditioned, ([l, ray]) => {
            const r = query.compute(l, ray);
            for (let k = 0; k < 20; ++k) {
                const p = add(l.origin, mul(20 * (rand() - 0.5), l.direction));
                const q = add(ray.origin, mul(12 * rand(), ray.direction));
                const diff = sub(p, q);
                const sqr = dot(diff, diff);
                expect(r.sqrDistance).toBeLessThanOrEqual(sqr + 1e-9 * (1 + sqr));
            }
        }, 60);
    });

    it('agrees with the line-segment query for a long segment', () => {
        check(wellConditioned, ([l, ray]) => {
            const r = query.compute(l, ray);
            const seg = Segment.fromEndpoints(ray.origin,
                add(ray.origin, mul(1e4, ray.direction)));
            const rs = lineSegment.compute(l, seg);
            if (r.parameter[1] < 1e3) {
                expectClose(r.distance, rs.distance, 1e-7, 1e-7);
            }
        });
    });

    it('is invariant under rigid motions', () => {
        check(fc.tuple(wellConditioned, rotationFrame(3),
            wellScaledVector(3, -5, 5)), ([[l, ray], R, tr]) => {
            const movedLine = Line.fromOriginDirection(
                add(rot(R, l.origin), tr), rot(R, l.direction));
            const movedRay = Ray.fromOriginDirection(
                add(rot(R, ray.origin), tr), rot(R, ray.direction));
            expectClose(query.compute(l, ray).distance,
                query.compute(movedLine, movedRay).distance, 1e-8, 1e-8);
        });
    });

    it('clamps to the ray origin for parallel configurations', () => {
        // A power-of-two scale keeps a01 = -k*a00 and a11 = k^2*a00 exact, so
        // the determinant is exactly zero and the parallel branch is taken.
        check(fc.tuple(nonUnitLine3, wellScaledVector(3, -8, 8),
            fc.constantFrom(1, -1, 2, -2, 0.5, -0.5, 4, -4)),
        ([l, o1, k]) => {
            const ray = Ray.fromOriginDirection(o1, mul(k, l.direction));
            const r = query.compute(l, ray);
            expect(r.parameter[1]).toBe(0);
            expectVectorClose(r.closest[1], o1, 0, 0);
            expectClose(r.distance, pointLine.compute(o1, l).distance, 1e-9,
                1e-9);
        });
    });
});
