import { describe, expect, it } from 'vitest';
import { Arc2 } from '../src/Arc2.js';
import { DistRay2Arc2 } from '../src/DistRay2Arc2.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { DistRay2Circle2 } from '../src/DistRay2Circle2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import {
    check, expectClose, expectVectorClose, fc, positive,
    rotationFrame, unitVector, wellScaled, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

// Build an arc on the circle with the given center and radius, traversing
// counterclockwise from angle a0 to angle a1.
function arc(center: number[], radius: number, a0: number,
    a1: number): Arc2 {
    const c = v(...center);
    const end0 = add(c, v(radius * Math.cos(a0), radius * Math.sin(a0)));
    const end1 = add(c, v(radius * Math.cos(a1), radius * Math.sin(a1)));
    return Arc2.fromCenterRadiusEnds(c, radius, end0, end1);
}

// The exact distance from a point to a ray, computed independently.
function pointRayDistance(p: Vector, r: Ray): number {
    const diff = sub(p, r.origin);
    let t = dot(r.direction, diff) / dot(r.direction, r.direction);
    if (t < 0) {
        t = 0;
    }
    return length(sub(p, add(r.origin, mul(t, r.direction))));
}

// Brute-force minimum distance between the ray and a dense sampling of the
// arc.
function bruteForce(r: Ray, a: Arc2, a0: number, a1: number): number {
    let hi = a1;
    while (hi < a0) {
        hi += 2 * Math.PI;
    }
    const at = (t: number): number => pointRayDistance(add(a.center,
        v(a.radius * Math.cos(t), a.radius * Math.sin(t))), r);
    const n = 4000;
    let best = Number.MAX_VALUE;
    let bt = a0;
    for (let i = 0; i <= n; ++i) {
        const t = a0 + (hi - a0) * i / n;
        const d = at(t);
        if (d < best) {
            best = d;
            bt = t;
        }
    }
    // Refine locally around the best sample; the distance function is smooth
    // and locally convex near its minimum.
    let h = (hi - a0) / n;
    for (let pass = 0; pass < 120; ++pass) {
        for (const sign of [1, -1]) {
            const t = Math.min(hi, Math.max(a0, bt + sign * h));
            const d = at(t);
            if (d < best) {
                best = d;
                bt = t;
            }
        }
        h *= 0.75;
    }
    return best;
}

describe('DistRay2Arc2', () => {
    const query = new DistRay2Arc2();
    // The quarter arc of the unit circle in the first quadrant.
    const quarter = arc([0, 0], 1, 0, Math.PI / 2);

    it('computes the distance when the closest circle point is on the arc',
        () => {
            const r = ray([5, 5], [1, 1]);
            const result = query.compute(r, quarter);
            expect(result.numClosestPairs).toBe(1);
            expect(result.distance).toBeCloseTo(5 * Math.SQRT2 - 1, 10);
            expect(result.closest[0][1].values[0]).toBeCloseTo(Math.SQRT1_2,
                10);
            expect(result.closest[0][1].values[1]).toBeCloseTo(Math.SQRT1_2,
                10);
        });

    it('reports zero distance when the ray meets the arc', () => {
        const r = ray([0, 0], [1, 1]);
        const result = query.compute(r, quarter);
        expect(result.distance).toBeCloseTo(0, 10);
        expect(result.numClosestPairs).toBeGreaterThanOrEqual(1);
    });

    it('falls back to an arc endpoint when the circle point is off the arc',
        () => {
            // The ray is far along the negative x-axis pointing away from
            // the arc; the closest arc point is the endpoint (0,1).
            const r = ray([-4, 1], [-1, 0]);
            const result = query.compute(r, quarter);
            expect(result.numClosestPairs).toBe(1);
            expect(result.distance).toBeCloseTo(4, 10);
            expect(result.closest[0][1].values[0]).toBeCloseTo(0, 10);
            expect(result.closest[0][1].values[1]).toBeCloseTo(1, 10);
        });

    it('reports two closest pairs when the arc endpoints are equidistant',
        () => {
            // The upper half of the unit circle and a ray below it pointing
            // down, so both endpoints (1,0) and (-1,0) are equidistant.
            const upper = arc([0, 0], 1, 0, Math.PI);
            const r = ray([0, -5], [0, -1]);
            const result = query.compute(r, upper);
            expect(result.numClosestPairs).toBe(2);
            expect(result.distance).toBeCloseTo(Math.sqrt(26), 10);
            const xs = [result.closest[0][1].values[0],
                result.closest[1][1].values[0]].sort((p, q) => p - q);
            expect(xs[0]).toBeCloseTo(-1, 10);
            expect(xs[1]).toBeCloseTo(1, 10);
        });

    it('handles a degenerate zero-radius arc', () => {
        const point = Arc2.fromCenterRadiusEnds(v(2, 0), 0, v(2, 0), v(2, 0));
        const r = ray([0, 0], [0, 1]);
        const result = query.compute(r, point);
        expect(result.distance).toBeCloseTo(2, 8);
    });

    it('agrees with a dense sampling of the arc on random inputs', () => {
        let seed = 60708090;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 30; ++trial) {
            const radius = 0.5 + 2 * rand();
            const a0 = 2 * Math.PI * rand();
            const a1 = a0 + 0.3 + 5.5 * rand();
            const center = [4 * rand() - 2, 4 * rand() - 2];
            const a = arc(center, radius, a0, a1);
            const r = ray([10 * rand() - 5, 10 * rand() - 5],
                [2 * rand() - 1, 2 * rand() - 1]);
            if (length(r.direction) < 1e-6) {
                continue;
            }
            const result = query.compute(r, a);
            const brute = bruteForce(r, a, a0, a1);
            expect(result.distance).toBeCloseTo(brute, 7);

            // Every reported pair realizes the reported distance, the ray
            // point is on the ray and the arc point is on the arc.
            for (let j = 0; j < result.numClosestPairs; ++j) {
                const onRay = add(r.origin,
                    mul(result.parameter[j], r.direction));
                expect(result.parameter[j]).toBeGreaterThanOrEqual(0);
                expect(length(sub(result.closest[j][0], onRay)))
                    .toBeLessThan(1e-8);
                expect(length(sub(result.closest[j][1], a.center)))
                    .toBeCloseTo(a.radius, 8);
                expect(length(sub(result.closest[j][0], result.closest[j][1])))
                    .toBeCloseTo(result.distance, 8);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistRay2Arc2.h. Upstream runs the ray-circle query and keeps the
// circle closest points that are on the arc; if none is, it sorts the two arc
// endpoint distances and the ray-origin-to-arc distance and takes the minima.
// ---------------------------------------------------------------------------

// An arc with a well separated pair of endpoints on a circle of moderate
// radius, plus the two angles that bound it (the port has no accessor for
// them).
const v21Arc = fc.tuple(wellScaledVector(2, -6, 6), positive(4, 0.2),
    wellScaled(-Math.PI, Math.PI),
    fc.double({ min: 0.05, max: 2 * Math.PI - 0.05, noNaN: true }))
    .map(([c, r, a0, sweep]) => {
        const a1 = a0 + sweep;
        const e0 = add(c, Vector.fromArray(
            [r * Math.cos(a0), r * Math.sin(a0)]));
        const e1 = add(c, Vector.fromArray(
            [r * Math.cos(a1), r * Math.sin(a1)]));
        return [Arc2.fromCenterRadiusEnds(c, r, e0, e1), a0, a1] as
            [Arc2, number, number];
    });

const v21Ray = fc.tuple(wellScaledVector(2, -8, 8), unitVector(2))
    .map(([o, d]) => Ray.fromOriginDirection(o, d));

// A point of the arc at fractional position u along its counterclockwise
// sweep.
function v21ArcPoint(a: Arc2, a0: number, a1: number, u: number): Vector {
    let hi = a1;
    while (hi < a0) {
        hi += 2 * Math.PI;
    }
    const t = a0 + (hi - a0) * u;
    return add(a.center,
        Vector.fromArray([a.radius * Math.cos(t), a.radius * Math.sin(t)]));
}

describe('DistRay2Arc2 verification', () => {
    const query = new DistRay2Arc2();

    it('every reported pair is on the ray and on the arc', () => {
        check(fc.tuple(v21Ray, v21Arc), ([r, [a, a0, a1]]) => {
            void a0;
            void a1;
            const res = query.compute(r, a);
            expect(res.numClosestPairs === 1 || res.numClosestPairs === 2)
                .toBe(true);
            expectClose(res.distance, Math.sqrt(res.sqrDistance), 1e-9, 1e-9);
            for (let j = 0; j < res.numClosestPairs; ++j) {
                const onRay = res.closest[j][0];
                const onArc = res.closest[j][1];
                // The ray point is at the reported parameter, which is >= 0.
                expect(res.parameter[j]).toBeGreaterThanOrEqual(0);
                expectVectorClose(onRay,
                    add(r.origin, mul(res.parameter[j], r.direction)), 1e-8,
                    1e-8);
                // The arc point is on the circle and inside the arc.
                expectClose(length(sub(onArc, a.center)), a.radius, 1e-8,
                    1e-8);
                expect(a.containsOnCircle(onArc)).toBe(true);
                // The pair realizes the reported distance.
                expectClose(length(sub(onRay, onArc)), res.distance, 1e-8,
                    1e-8);
            }
        });
    });

    it('matches a brute-force minimization over the arc', () => {
        check(fc.tuple(v21Ray, v21Arc), ([r, [a, a0, a1]]) => {
            expectClose(query.compute(r, a).distance, bruteForce(r, a, a0, a1),
                1e-6, 1e-6);
        }, 60);
    }, 30000);

    it('is not larger than the distance to any sampled arc point', () => {
        check(fc.tuple(v21Ray, v21Arc,
            fc.double({ min: 0, max: 1, noNaN: true })),
        ([r, [a, a0, a1], u]) => {
            const q = v21ArcPoint(a, a0, a1, u);
            expect(query.compute(r, a).distance)
                .toBeLessThanOrEqual(pointRayDistance(q, r) + 1e-8);
        });
    });

    it('reports zero distance when the ray starts on the arc', () => {
        check(fc.tuple(v21Arc, fc.double({ min: 0, max: 1, noNaN: true }),
            unitVector(2)), ([[a, a0, a1], u, d]) => {
            const q = v21ArcPoint(a, a0, a1, u);
            const r = Ray.fromOriginDirection(q, d);
            expect(query.compute(r, a).distance).toBeLessThanOrEqual(1e-8);
        });
    });

    it('agrees with the ray-circle query when the arc is the whole circle',
        () => {
            // An arc whose endpoints coincide has DotPerp(P-E0, E1-E0) = 0 for
            // every P, so Contains accepts every circle point and the query
            // degenerates to the ray-circle query.
            check(fc.tuple(v21Ray, wellScaledVector(2, -6, 6),
                positive(4, 0.2), wellScaled(-Math.PI, Math.PI)),
            ([r, c, radius, a0]) => {
                const e = add(c, Vector.fromArray(
                    [radius * Math.cos(a0), radius * Math.sin(a0)]));
                const a = Arc2.fromCenterRadiusEnds(c, radius, e, e.clone());
                const circle = Hypersphere.fromCenterRadius(c, radius);
                const rc = new DistRay2Circle2().compute(r, circle);
                const ra = query.compute(r, a);
                expect(ra.numClosestPairs).toBe(rc.numClosestPairs);
                expectClose(ra.distance, rc.distance, 0, 0);
            });
        });

    it('is equivariant under rigid motions', () => {
        // Rotations preserve the counterclockwise ordering of the arc
        // endpoints, so the moved arc is still a valid Arc2.
        check(fc.tuple(v21Ray, v21Arc, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([r, [a, a0, a1], R, tr]) => {
            void a0;
            void a1;
            const rot = (x: Vector): Vector => add(mul(x.values[0], R[0]),
                mul(x.values[1], R[1]));
            const moved = Arc2.fromCenterRadiusEnds(add(rot(a.center), tr),
                a.radius, add(rot(a.end[0]), tr), add(rot(a.end[1]), tr));
            const movedRay = Ray.fromOriginDirection(add(rot(r.origin), tr),
                rot(r.direction));
            expectClose(query.compute(r, a).distance,
                query.compute(movedRay, moved).distance, 1e-7, 1e-7);
        });
    });

    it('does not mutate its inputs', () => {
        check(fc.tuple(v21Ray, v21Arc), ([r, [a, a0, a1]]) => {
            void a0;
            void a1;
            const o = r.origin.clone();
            const d = r.direction.clone();
            const snapshot = [...a.center.values, a.radius,
                ...a.end[0].values, ...a.end[1].values];
            const res = query.compute(r, a);
            expect(r.origin.values).toEqual(o.values);
            expect(r.direction.values).toEqual(d.values);
            expect([...a.center.values, a.radius, ...a.end[0].values,
                ...a.end[1].values]).toEqual(snapshot);
            for (let j = 0; j < res.numClosestPairs; ++j) {
                res.closest[j][0].values[0] = 777;
                res.closest[j][1].values[0] = 777;
            }
            expect(r.origin.values).toEqual(o.values);
            expect([...a.center.values, a.radius, ...a.end[0].values,
                ...a.end[1].values]).toEqual(snapshot);
        });
    });
});
