import { describe, expect, it } from 'vitest';
import { Arc2 } from '../src/Arc2';
import { DistSegment2Arc2 } from '../src/DistSegment2Arc2';
import { Segment } from '../src/Segment';
import { Vector, add, dot, length, mul, sub } from '../src/Vector';

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
