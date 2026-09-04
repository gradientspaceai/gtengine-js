import { describe, expect, it } from 'vitest';
import { Arc2 } from '../src/Arc2.js';
import { DistSegment2Arc2 } from '../src/DistSegment2Arc2.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { DistSegment2Circle2 } from '../src/DistSegment2Circle2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import {
    check, expectClose, expectVectorClose, fc, positive,
    rotationFrame, unitVector, wellScaled, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
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

// The exact distance from a point to a segment, computed independently.
function pointSegmentDistance(p: Vector, s: Segment): number {
    const direction = sub(s.p[1], s.p[0]);
    const dd = dot(direction, direction);
    let t = dd > 0 ? dot(direction, sub(p, s.p[0])) / dd : 0;
    t = Math.max(0, Math.min(1, t));
    return length(sub(p, add(s.p[0], mul(t, direction))));
}

// Brute-force minimum distance between the segment and the arc: a coarse
// sampling of the arc followed by a local refinement.
function bruteForce(s: Segment, a: Arc2, a0: number, a1: number): number {
    let hi = a1;
    while (hi < a0) {
        hi += 2 * Math.PI;
    }
    const at = (t: number): number => pointSegmentDistance(add(a.center,
        v(a.radius * Math.cos(t), a.radius * Math.sin(t))), s);
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

describe('DistSegment2Arc2', () => {
    const query = new DistSegment2Arc2();
    // The quarter arc of the unit circle in the first quadrant.
    const quarter = arc([0, 0], 1, 0, Math.PI / 2);

    it('computes the distance when the closest circle point is on the arc',
        () => {
            const s = segment([5, 5], [7, 7]);
            const result = query.compute(s, quarter);
            expect(result.numClosestPairs).toBe(1);
            expect(result.distance).toBeCloseTo(5 * Math.SQRT2 - 1, 10);
            expect(result.closest[0][1].values[0]).toBeCloseTo(Math.SQRT1_2,
                10);
            expect(result.closest[0][1].values[1]).toBeCloseTo(Math.SQRT1_2,
                10);
        });

    it('reports zero distance when the segment crosses the arc', () => {
        const s = segment([0, 0], [2, 2]);
        const result = query.compute(s, quarter);
        expect(result.distance).toBeCloseTo(0, 10);
    });

    it('falls back to an arc endpoint when the circle point is off the arc',
        () => {
            const s = segment([-4, 1], [-6, 1]);
            const result = query.compute(s, quarter);
            expect(result.numClosestPairs).toBe(1);
            expect(result.distance).toBeCloseTo(4, 10);
            expect(result.closest[0][1].values[0]).toBeCloseTo(0, 10);
            expect(result.closest[0][1].values[1]).toBeCloseTo(1, 10);
        });

    it('reports two closest pairs when the arc endpoints are equidistant',
        () => {
            // The upper half of the unit circle and a horizontal segment
            // below it, centered on the y-axis, so both endpoints (1,0) and
            // (-1,0) are equidistant from the segment.
            const upper = arc([0, 0], 1, 0, Math.PI);
            const s = segment([-0.25, -5], [0.25, -5]);
            const result = query.compute(s, upper);
            expect(result.numClosestPairs).toBe(2);
            const xs = [result.closest[0][1].values[0],
                result.closest[1][1].values[0]].sort((p, q) => p - q);
            expect(xs[0]).toBeCloseTo(-1, 10);
            expect(xs[1]).toBeCloseTo(1, 10);
        });

    it('handles a degenerate zero-length segment', () => {
        const s = segment([5, 0], [5, 0]);
        const result = query.compute(s, quarter);
        expect(result.distance).toBeCloseTo(4, 10);
        expect(result.closest[0][1].values[0]).toBeCloseTo(1, 10);
    });

    it('handles a degenerate zero-radius arc', () => {
        const point = Arc2.fromCenterRadiusEnds(v(2, 0), 0, v(2, 0), v(2, 0));
        const s = segment([0, 0], [0, 4]);
        const result = query.compute(s, point);
        expect(result.distance).toBeCloseTo(2, 8);
    });

    it('agrees with a dense sampling of the arc on random inputs', () => {
        let seed = 13141516;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 40; ++trial) {
            const radius = 0.5 + 2 * rand();
            const a0 = 2 * Math.PI * rand();
            const a1 = a0 + 0.3 + 5.5 * rand();
            const center = [4 * rand() - 2, 4 * rand() - 2];
            const a = arc(center, radius, a0, a1);
            const s = segment([10 * rand() - 5, 10 * rand() - 5],
                [10 * rand() - 5, 10 * rand() - 5]);
            const result = query.compute(s, a);
            const brute = bruteForce(s, a, a0, a1);
            expect(result.distance).toBeCloseTo(brute, 7);

            const direction = sub(s.p[1], s.p[0]);
            for (let j = 0; j < result.numClosestPairs; ++j) {
                expect(result.parameter[j]).toBeGreaterThanOrEqual(-1e-12);
                expect(result.parameter[j]).toBeLessThanOrEqual(1 + 1e-12);
                const onSeg = add(s.p[0], mul(result.parameter[j], direction));
                expect(length(sub(result.closest[j][0], onSeg)))
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
// header DistSegment2Arc2.h. Upstream runs the segment-circle query and keeps
// the circle closest points that are on the arc; if none is, it sorts the two
// arc endpoint distances and the two segment endpoint distances to the arc
// and takes the minima.
// ---------------------------------------------------------------------------

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

const v21Segment = fc.tuple(wellScaledVector(2, -8, 8),
    wellScaledVector(2, -8, 8))
    .filter(([a, b]) => length(sub(b, a)) > 1e-2)
    .map(([a, b]) => Segment.fromEndpoints(a, b));

function v21ArcPoint(a: Arc2, a0: number, a1: number, u: number): Vector {
    let hi = a1;
    while (hi < a0) {
        hi += 2 * Math.PI;
    }
    const t = a0 + (hi - a0) * u;
    return add(a.center,
        Vector.fromArray([a.radius * Math.cos(t), a.radius * Math.sin(t)]));
}

describe('DistSegment2Arc2 verification', () => {
    const query = new DistSegment2Arc2();

    it('every reported pair is on the segment and on the arc', () => {
        check(fc.tuple(v21Segment, v21Arc), ([seg, [a, a0, a1]]) => {
            void a0;
            void a1;
            const res = query.compute(seg, a);
            expect(res.numClosestPairs === 1 || res.numClosestPairs === 2)
                .toBe(true);
            expectClose(res.distance, Math.sqrt(res.sqrDistance), 1e-9, 1e-9);
            const dir = sub(seg.p[1], seg.p[0]);
            for (let j = 0; j < res.numClosestPairs; ++j) {
                const onSeg = res.closest[j][0];
                const onArc = res.closest[j][1];
                // Upstream parameterizes the segment as P0 + t*(P1-P0) with
                // 0 <= t <= 1.
                expect(res.parameter[j]).toBeGreaterThanOrEqual(0);
                expect(res.parameter[j]).toBeLessThanOrEqual(1);
                expectVectorClose(onSeg,
                    add(seg.p[0], mul(res.parameter[j], dir)), 1e-8, 1e-8);
                expectClose(length(sub(onArc, a.center)), a.radius, 1e-8,
                    1e-8);
                expect(a.containsOnCircle(onArc)).toBe(true);
                expectClose(length(sub(onSeg, onArc)), res.distance, 1e-8,
                    1e-8);
            }
        });
    });

    it('matches a brute-force minimization over the arc', () => {
        check(fc.tuple(v21Segment, v21Arc), ([seg, [a, a0, a1]]) => {
            expectClose(query.compute(seg, a).distance,
                bruteForce(seg, a, a0, a1), 1e-6, 1e-6);
        }, 60);
    }, 30000);

    it('is not larger than the distance to any sampled arc point', () => {
        check(fc.tuple(v21Segment, v21Arc,
            fc.double({ min: 0, max: 1, noNaN: true })),
        ([seg, [a, a0, a1], u]) => {
            const q = v21ArcPoint(a, a0, a1, u);
            expect(query.compute(seg, a).distance)
                .toBeLessThanOrEqual(pointSegmentDistance(q, seg) + 1e-8);
        });
    });

    it('reports zero distance when an endpoint is on the arc', () => {
        check(fc.tuple(v21Arc, fc.double({ min: 0, max: 1, noNaN: true }),
            wellScaledVector(2, -6, 6)), ([[a, a0, a1], u, other]) => {
            const q = v21ArcPoint(a, a0, a1, u);
            const seg = Segment.fromEndpoints(q, add(q, other));
            expect(query.compute(seg, a).distance).toBeLessThanOrEqual(1e-8);
        });
    });

    it('agrees with the segment-circle query when the arc is a full circle',
        () => {
            check(fc.tuple(v21Segment, wellScaledVector(2, -6, 6),
                positive(4, 0.2), wellScaled(-Math.PI, Math.PI)),
            ([seg, c, radius, a0]) => {
                const e = add(c, Vector.fromArray(
                    [radius * Math.cos(a0), radius * Math.sin(a0)]));
                const a = Arc2.fromCenterRadiusEnds(c, radius, e, e.clone());
                const circle = Hypersphere.fromCenterRadius(c, radius);
                const sc = new DistSegment2Circle2().compute(seg, circle);
                if (sc.numClosestPairs === 0) {
                    // Upstream DistSegment2Circle2 returns a default-
                    // constructed result (no pairs, zero distance) when the
                    // segment is strictly inside the circle; see the API
                    // notes. DistSegment2Arc2 recovers because it falls
                    // through to the endpoint comparison, which is checked by
                    // the brute-force property and by the dedicated test
                    // below.
                    return;
                }
                const sa = query.compute(seg, a);
                expect(sa.numClosestPairs).toBe(sc.numClosestPairs);
                expectClose(sa.distance, sc.distance, 0, 0);
            });
        });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(v21Segment, v21Arc, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([seg, [a, a0, a1], R, tr]) => {
            void a0;
            void a1;
            const rot = (x: Vector): Vector => add(mul(x.values[0], R[0]),
                mul(x.values[1], R[1]));
            const moved = Arc2.fromCenterRadiusEnds(add(rot(a.center), tr),
                a.radius, add(rot(a.end[0]), tr), add(rot(a.end[1]), tr));
            const movedSeg = Segment.fromEndpoints(add(rot(seg.p[0]), tr),
                add(rot(seg.p[1]), tr));
            expectClose(query.compute(seg, a).distance,
                query.compute(movedSeg, moved).distance, 1e-7, 1e-7);
        });
    });

    it('does not mutate its inputs', () => {
        check(fc.tuple(v21Segment, v21Arc), ([seg, [a, a0, a1]]) => {
            void a0;
            void a1;
            const p0 = seg.p[0].clone();
            const p1 = seg.p[1].clone();
            const snapshot = [...a.center.values, a.radius,
                ...a.end[0].values, ...a.end[1].values];
            const res = query.compute(seg, a);
            expect(seg.p[0].values).toEqual(p0.values);
            expect(seg.p[1].values).toEqual(p1.values);
            expect([...a.center.values, a.radius, ...a.end[0].values,
                ...a.end[1].values]).toEqual(snapshot);
            for (let j = 0; j < res.numClosestPairs; ++j) {
                res.closest[j][0].values[0] = 777;
                res.closest[j][1].values[0] = 777;
            }
            expect(seg.p[0].values).toEqual(p0.values);
            expect([...a.center.values, a.radius, ...a.end[0].values,
                ...a.end[1].values]).toEqual(snapshot);
        });
    });
    it('is correct for a segment strictly inside the circle of the arc', () => {
        // DistSegment2Circle2 reports no closest pairs here (upstream
        // defect), so the query falls through to the arc-endpoint and
        // segment-endpoint comparison. The true minimum is r - max_i |Pi-C|
        // when the nearest circle point of the farther endpoint is on the
        // arc, so use the full circle to make that certain.
        check(fc.tuple(wellScaledVector(2, -1, 1), wellScaledVector(2, -1, 1),
            positive(4, 1.5), wellScaled(-Math.PI, Math.PI)),
        ([p0, p1, radius, a0]) => {
            if (length(sub(p1, p0)) < 1e-2) {
                return;
            }
            const c = new Vector(2);
            const e = Vector.fromArray(
                [radius * Math.cos(a0), radius * Math.sin(a0)]);
            const a = Arc2.fromCenterRadiusEnds(c, radius, e, e.clone());
            const seg = Segment.fromEndpoints(p0, p1);
            const expected = radius - Math.max(length(p0), length(p1));
            const res = query.compute(seg, a);
            expectClose(res.distance, expected, 1e-9, 1e-9);
            expect(res.numClosestPairs).toBeGreaterThanOrEqual(1);
        });
    });
});
