import { describe, expect, it } from 'vitest';
import { DistLineSegment } from '../src/DistLineSegment.js';
import { Line } from '../src/Line.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { DistLineLine } from '../src/DistLineLine.js';
import { DistPointLine } from '../src/DistPointLine.js';
import { DistSegmentSegment } from '../src/DistSegmentSegment.js';
import { check, expectClose, expectVectorClose, fc, rotationFrame, seededRandom, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

describe('DistLineSegment', () => {
    const query = new DistLineSegment();

    it('uses interior points when the segment parameter is in [0,1]', () => {
        // The x-axis and the segment from (2,3,4) to (2,3,-4). The closest
        // segment point is (2,3,0) at s1 = 0.5, and the distance is 3.
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            seg([2, 3, 4], [2, 3, -4]));
        expect(result.parameter[1]).toBeCloseTo(0.5, 10);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[0].values[0]).toBeCloseTo(2, 10);
    });

    it('clamps to the endpoint Q1', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            seg([0, 3, 5], [0, 3, 1]));
        expect(result.parameter[1]).toBe(1);
        expect(result.closest[1].values).toEqual([0, 3, 1]);
        expect(result.distance).toBeCloseTo(Math.sqrt(9 + 1), 10);
    });

    it('clamps to the endpoint Q0', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            seg([0, 3, 1], [0, 3, 5]));
        expect(result.parameter[1]).toBe(0);
        expect(result.closest[1].values).toEqual([0, 3, 1]);
        expect(result.distance).toBeCloseTo(Math.sqrt(9 + 1), 10);
    });

    it('handles a line parallel to the segment', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            seg([2, 4, 0], [7, 4, 0]));
        expect(result.parameter[1]).toBe(0);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values).toEqual([2, 4, 0]);
    });

    it('handles a degenerate (zero-length) segment', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            seg([2, 5, 0], [2, 5, 0]));
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values).toEqual([2, 5, 0]);
    });

    it('reports zero distance when the segment crosses the line', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            seg([3, -1, 0], [3, 1, 0]));
        expect(result.distance).toBeCloseTo(0, 10);
        expect(result.parameter[1]).toBeCloseTo(0.5, 10);
    });

    it('agrees with a sampled minimum over the segment', () => {
        let seed = 8080;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const l = line([rand(), rand(), rand()],
                [rand() + 4, rand(), rand()]);
            const s = seg([rand(), rand(), rand()],
                [rand(), rand(), rand()]);
            const result = query.compute(l, s);

            expect(result.parameter[1]).toBeGreaterThanOrEqual(0);
            expect(result.parameter[1]).toBeLessThanOrEqual(1);

            const segDir = sub(s.p[1], s.p[0]);
            const a00 = dot(l.direction, l.direction);
            let best = Number.MAX_VALUE;
            for (let k = 0; k <= 2000; ++k) {
                const s1 = k / 2000;
                const q = add(s.p[0], mul(s1, segDir));
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
// Verification wave (V19): property-based cross-checks of DistLineSegment.ts
// against the upstream header DistLineSegment.h.
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

const segment3 = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -8, 8))
    .filter(([p0, p1]) => length(sub(p1, p0)) > 0.25)
    .map(([p0, p1]) => Segment.fromEndpoints(p0, p1));

// Pairs whose supporting directions are well away from parallel; near
// parallel configurations lose all significance in the upstream determinant.
const wellConditioned = fc.tuple(nonUnitLine3, segment3)
    .filter(([l, seg]) => {
        const sd = sub(seg.p[1], seg.p[0]);
        return length(cross(l.direction, sd))
            > 0.2 * length(l.direction) * length(sd);
    });

describe('DistLineSegment verification', () => {
    const query = new DistLineSegment();
    const lineLine = new DistLineLine();
    const pointLine = new DistPointLine();
    const segSeg = new DistSegmentSegment();

    it('result is self consistent and the segment parameter is in [0,1]',
        () => {
            check(fc.tuple(nonUnitLine3, segment3), ([l, seg]) => {
                const r = query.compute(l, seg);
                expect(r.parameter[1]).toBeGreaterThanOrEqual(0);
                expect(r.parameter[1]).toBeLessThanOrEqual(1);
                expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
                const diff = sub(r.closest[0], r.closest[1]);
                expectClose(r.sqrDistance, dot(diff, diff), 1e-12, 1e-12);
                expectVectorClose(r.closest[0],
                    add(l.origin, mul(r.parameter[0], l.direction)), 1e-12,
                    1e-12);
                expectVectorClose(r.closest[1],
                    add(seg.p[0], mul(r.parameter[1], sub(seg.p[1], seg.p[0]))),
                    1e-12, 1e-12);
            });
        });

    it('matches the exact minimum over the domain R x [0,1]', () => {
        check(wellConditioned, ([l, seg]) => {
            // Convex quadratic: the minimum is the unconstrained line/line
            // critical point when its segment parameter is in [0,1], and
            // otherwise lies on one of the two boundary faces s1 = 0, 1.
            const supporting = Line.fromOriginDirection(seg.p[0],
                sub(seg.p[1], seg.p[0]));
            const unconstrained = lineLine.compute(l, supporting);
            let ref = Math.min(pointLine.compute(seg.p[0], l).distance,
                pointLine.compute(seg.p[1], l).distance);
            const s1 = unconstrained.parameter[1];
            if (s1 >= 0 && s1 <= 1) {
                ref = Math.min(ref, unconstrained.distance);
            }
            expectClose(query.compute(l, seg).distance, ref, 1e-8, 1e-8);
        });
    });

    it('is minimal over sampled point pairs', () => {
        const rand = seededRandom(0x51d7);
        check(wellConditioned, ([l, seg]) => {
            const r = query.compute(l, seg);
            const sd = sub(seg.p[1], seg.p[0]);
            for (let k = 0; k < 20; ++k) {
                const p = add(l.origin, mul(20 * (rand() - 0.5), l.direction));
                const q = add(seg.p[0], mul(rand(), sd));
                const diff = sub(p, q);
                const sqr = dot(diff, diff);
                expect(r.sqrDistance).toBeLessThanOrEqual(sqr + 1e-9 * (1 + sqr));
            }
        }, 60);
    });

    it('agrees with the segment-segment query for a long line segment', () => {
        check(wellConditioned, ([l, seg]) => {
            const r = query.compute(l, seg);
            const longSeg = Segment.fromEndpoints(
                add(l.origin, mul(-1e4, l.direction)),
                add(l.origin, mul(1e4, l.direction)));
            if (Math.abs(r.parameter[0]) < 1e3) {
                expectClose(r.distance,
                    segSeg.compute(longSeg, seg).distance, 1e-6, 1e-6);
            }
        });
    });

    it('is invariant under rigid motions', () => {
        check(fc.tuple(wellConditioned, rotationFrame(3),
            wellScaledVector(3, -5, 5)), ([[l, seg], R, tr]) => {
            const movedLine = Line.fromOriginDirection(
                add(rot(R, l.origin), tr), rot(R, l.direction));
            const movedSeg = Segment.fromEndpoints(
                add(rot(R, seg.p[0]), tr), add(rot(R, seg.p[1]), tr));
            expectClose(query.compute(l, seg).distance,
                query.compute(movedLine, movedSeg).distance, 1e-8, 1e-8);
        });
    });

    it('selects the endpoint Q0 for parallel configurations', () => {
        // The line direction is a power-of-two multiple of the segment
        // direction, so a01 = -k*a11 and a00 = k^2*a11 hold exactly in
        // binary64 and the determinant is exactly zero. Building the segment
        // from the line instead would not work: Q0 + k*D does not subtract
        // back to exactly k*D, and the residual leaves a tiny positive
        // determinant.
        check(fc.tuple(segment3, wellScaledVector(3, -8, 8),
            fc.constantFrom(1, -1, 2, -2, 0.5, -0.5, 4, -4)),
        ([seg, origin, k]) => {
            const q0 = seg.p[0];
            const l = Line.fromOriginDirection(origin,
                mul(k, sub(seg.p[1], seg.p[0])));
            const r = query.compute(l, seg);
            expect(r.parameter[1]).toBe(0);
            expectVectorClose(r.closest[1], q0, 0, 0);
            expectClose(r.distance, pointLine.compute(q0, l).distance, 1e-9,
                1e-9);
        });
    });

    it('handles a zero-length segment as a point-line distance', () => {
        check(fc.tuple(nonUnitLine3, wellScaledVector(3, -8, 8)), ([l, q]) => {
            const r = query.compute(l, Segment.fromEndpoints(q, q));
            // a11 = 0 and a01 = 0, so the determinant is zero and upstream
            // takes the parallel branch with s1 = 0.
            expect(r.parameter[1]).toBe(0);
            expectClose(r.distance, pointLine.compute(q, l).distance, 1e-9,
                1e-9);
        });
    });
});
