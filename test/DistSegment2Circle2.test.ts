import { describe, expect, it } from 'vitest';
import { DistPoint2Circle2 } from '../src/DistPoint2Circle2.js';
import { DistSegment2Circle2 } from '../src/DistSegment2Circle2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

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
