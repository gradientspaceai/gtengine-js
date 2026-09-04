import { describe, expect, it } from 'vitest';
import { DistPoint2Circle2 } from '../src/DistPoint2Circle2.js';
import { DistSegment2Circle2 } from '../src/DistSegment2Circle2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

function circle(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(v(...center), radius);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistSegment2Circle2', () => {
    const query = new DistSegment2Circle2();
    const unit = circle([0, 0], 1);

    it('reports two closest pairs when the segment crosses the circle', () => {
        const result = query.compute(segment([-2, 0], [2, 0]), unit);
        expect(result.numClosestPairs).toBe(2);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.closest[0][1].values[0]).toBeCloseTo(-1, 10);
        expect(result.closest[1][1].values[0]).toBeCloseTo(1, 10);
    });

    it('drops the t0-point when the segment starts inside the circle', () => {
        const result = query.compute(segment([0, 0], [3, 0]), unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.closest[0][1].values[0]).toBeCloseTo(1, 10);
    });

    it('drops the t1-point when the segment ends inside the circle', () => {
        const result = query.compute(segment([-3, 0], [0, 0]), unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.closest[0][1].values[0]).toBeCloseTo(-1, 10);
    });

    it('reports no pairs when the segment is strictly inside the circle',
        () => {
            const result = query.compute(segment([-0.5, 0], [0.5, 0]), unit);
            expect(result.numClosestPairs).toBe(0);
            expect(result.distance).toBe(0);
        });

    it('uses the far endpoint when the segment lies before the circle', () => {
        const pcQuery = new DistPoint2Circle2();
        const seg = segment([5, 0], [8, 0]);
        const result = query.compute(seg, unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.parameter[0]).toBe(0);
        const pcResult = pcQuery.compute(seg.p[0], unit);
        expect(result.distance).toBeCloseTo(pcResult.distance, 12);
        expect(result.closest[0][1].values[0]).toBeCloseTo(1, 10);
    });

    it('uses the near endpoint when the segment lies past the circle', () => {
        const seg = segment([-8, 0], [-5, 0]);
        const result = query.compute(seg, unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.parameter[0]).toBe(1);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[0][1].values[0]).toBeCloseTo(-1, 10);
    });

    it('handles a segment whose line misses the circle', () => {
        const result = query.compute(segment([-2, 3], [2, 3]), unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(2, 10);
    });

    it('handles a nearly degenerate (very short) segment', () => {
        // The segment direction P1-P0 must be nonzero: the query divides by
        // Dot(D,D), so a zero-length segment is an invalid input upstream and
        // here.
        const result = query.compute(segment([4, 0], [4 + 1e-9, 0]), unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(3, 8);
        expect(result.closest[0][1].values[0]).toBeCloseTo(1, 8);
    });

    it('handles a zero-radius circle', () => {
        const point = circle([0, 0], 0);
        const result = query.compute(segment([3, 4], [6, 8]), point);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(5, 10);
    });

    it('agrees with a dense sampling of the segment and circle', () => {
        const rnd = makeRandom(4004);
        const center = v(0.25, -0.5);
        const radius = 1.1;
        const c = Hypersphere.fromCenterRadius(center, radius);

        for (let trial = 0; trial < 60; ++trial) {
            const p0 = v(5 * rnd() - 2.5, 5 * rnd() - 2.5);
            const p1 = v(5 * rnd() - 2.5, 5 * rnd() - 2.5);
            const seg = Segment.fromEndpoints(p0, p1);
            const dir = sub(p1, p0);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const result = query.compute(seg, c);

            for (let j = 0; j < result.numClosestPairs; ++j) {
                expect(result.parameter[j]).toBeGreaterThanOrEqual(-1e-12);
                expect(result.parameter[j]).toBeLessThanOrEqual(1 + 1e-12);
                const p = add(p0, mul(result.parameter[j], dir));
                for (let i = 0; i < 2; ++i) {
                    expect(p.values[i]).toBeCloseTo(
                        result.closest[j][0].values[i], 7);
                }
                const d = sub(result.closest[j][1], center);
                expect(Math.sqrt(dot(d, d))).toBeCloseTo(radius, 7);
                const e = sub(result.closest[j][0], result.closest[j][1]);
                expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 6);
            }

            // Brute-force minimum over sampled circle points, with the
            // segment parameter clamped to [0,1].
            const n = 3000;
            const dd = dot(dir, dir);
            let best = Number.MAX_VALUE;
            for (let i = 0; i < n; ++i) {
                const t = 2 * Math.PI * i / n;
                const q = add(center, v(radius * Math.cos(t),
                    radius * Math.sin(t)));
                const w = sub(q, p0);
                let s = dot(w, dir) / dd;
                s = Math.max(0, Math.min(1, s));
                const f = sub(w, mul(s, dir));
                best = Math.min(best, dot(f, f));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistSegment2Circle2.h.
// ---------------------------------------------------------------------------

describe('DistSegment2Circle2 verification', () => {
    const query = new DistSegment2Circle2();

    const circleArb = fc.tuple(wellScaledVector(2, -5, 5), finite(0.25, 4))
        .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));

    const segmentArb = fc.tuple(wellScaledVector(2, -8, 8),
        wellScaledVector(2, -8, 8))
        .filter(([p0, p1]) => length(sub(p1, p0)) > 0.5)
        .map(([p0, p1]) => Segment.fromEndpoints(p0, p1));

    function pointCircleDistance(p: Vector, c: Hypersphere): number {
        return Math.abs(length(sub(p, c.center)) - c.radius);
    }

    // |P(t)-C| is convex on [0,1], so | |P(t)-C| - r | has at most two local
    // minima; a dense scan plus a golden-section refinement finds the global
    // one.
    function bruteForce(seg: Segment, c: Hypersphere): number {
        const d = sub(seg.p[1], seg.p[0]);
        const f = (t: number): number => pointCircleDistance(
            add(seg.p[0], mul(t, d)), c);
        const n = 20000;
        let bestI = 0;
        let best = Number.POSITIVE_INFINITY;
        for (let i = 0; i <= n; ++i) {
            const y = f(i / n);
            if (y < best) { best = y; bestI = i; }
        }
        let lo = Math.max(bestI / n - 1 / n, 0);
        let hi = Math.min(lo + 2 / n, 1);
        const phi = (Math.sqrt(5) - 1) / 2;
        for (let i = 0; i < 200; ++i) {
            const m0 = hi - phi * (hi - lo);
            const m1 = lo + phi * (hi - lo);
            if (f(m0) <= f(m1)) { hi = m1; } else { lo = m0; }
        }
        return Math.min(best, f(0.5 * (lo + hi)));
    }

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(segmentArb, circleArb), ([seg, circle]) => {
            const r = query.compute(seg, circle);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            const d = sub(seg.p[1], seg.p[0]);
            for (let j = 0; j < r.numClosestPairs; ++j) {
                expect(r.parameter[j]).toBeGreaterThanOrEqual(0);
                expect(r.parameter[j]).toBeLessThanOrEqual(1);
                expectVectorClose(r.closest[j][0],
                    add(seg.p[0], mul(r.parameter[j], d)), 1e-8, 1e-8);
                expectClose(length(sub(r.closest[j][1], circle.center)),
                    circle.radius, 1e-8, 1e-8);
                expectClose(length(sub(r.closest[j][0], r.closest[j][1])),
                    r.distance, 1e-8, 1e-8);
            }
        });
    });

    it('matches an independent minimization along the segment', () => {
        check(fc.tuple(segmentArb, circleArb), ([seg, circle]) => {
            const r = query.compute(seg, circle);
            // The "segment strictly inside the circle" case is the upstream
            // quirk documented in src/DistSegment2Circle2.ts: the whole
            // result is reset, so the reported distance is not the segment
            // distance. It is pinned by its own test below.
            if (r.numClosestPairs === 0) { return; }
            expectClose(r.distance, bruteForce(seg, circle), 1e-7, 1e-7);
        }, 60);
    }, 30000);

    it('preserves the upstream zeroed result for an interior segment', () => {
        // Upstream resets the whole Result, so numClosestPairs is 0 (the
        // header documents 1 or 2) and the distance is 0 even though the
        // segment does not touch the circle. The true minimum is attained at
        // the endpoint farthest from the center.
        check(fc.tuple(circleArb, finite(0, 0.4), finite(-Math.PI, Math.PI),
            finite(0, 0.4), finite(-Math.PI, Math.PI)),
            ([circle, f0, a0, f1, a1]) => {
                const at = (f: number, a: number): Vector =>
                    add(circle.center, v(f * circle.radius * Math.cos(a),
                        f * circle.radius * Math.sin(a)));
                const p0 = at(f0, a0);
                const p1 = at(f1, a1);
                if (length(sub(p1, p0)) < 1e-3) { return; }
                const seg = Segment.fromEndpoints(p0, p1);
                const r = query.compute(seg, circle);
                expect(r.numClosestPairs).toBe(0);
                expect(r.distance).toBe(0);
                // The distance the query should have reported.
                const trueDistance = Math.min(
                    pointCircleDistance(p0, circle),
                    pointCircleDistance(p1, circle));
                expect(trueDistance).toBeGreaterThan(0);
                expectClose(trueDistance, bruteForce(seg, circle), 1e-7, 1e-7);
            }, 60);
    }, 30000);

    it('reduces to the point-circle query for a short segment', () => {
        check(fc.tuple(circleArb, wellScaledVector(2, -8, 8), unitVector(2)),
            ([circle, p, dir]) => {
                const seg = Segment.fromEndpoints(p, add(p, mul(1e-6, dir)));
                const r = query.compute(seg, circle);
                if (r.numClosestPairs === 0) { return; }
                expectClose(r.distance, pointCircleDistance(p, circle),
                    1e-5, 1e-5);
            });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(segmentArb, circleArb, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([seg, circle, frame, shift]) => {
            const rot = (p: Vector): Vector =>
                v(frame[0].values[0] * p.values[0]
                    + frame[1].values[0] * p.values[1],
                    frame[0].values[1] * p.values[0]
                    + frame[1].values[1] * p.values[1]);
            const movedSeg = Segment.fromEndpoints(add(shift, rot(seg.p[0])),
                add(shift, rot(seg.p[1])));
            const movedCircle = Hypersphere.fromCenterRadius(
                add(shift, rot(circle.center)), circle.radius);
            const r0 = query.compute(seg, circle);
            const r1 = query.compute(movedSeg, movedCircle);
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
            expect(r0.numClosestPairs).toBe(r1.numClosestPairs);
        });
    });

    it('reverses the parameters when the endpoints are swapped', () => {
        check(fc.tuple(segmentArb, circleArb), ([seg, circle]) => {
            const reversed = Segment.fromEndpoints(seg.p[1], seg.p[0]);
            const r0 = query.compute(seg, circle);
            const r1 = query.compute(reversed, circle);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expect(r0.numClosestPairs).toBe(r1.numClosestPairs);
        });
    });
});
