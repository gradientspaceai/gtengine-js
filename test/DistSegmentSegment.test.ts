import { describe, expect, it } from 'vitest';
import {
    DistSegmentSegment, type DistSegmentSegmentResult
} from '../src/DistSegmentSegment';
import { Segment } from '../src/Segment';
import { Vector, add, dot, mul, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

// A near-exact reference minimum: sample the parameter of segment1 densely
// and, for each sample, solve for the optimal parameter of segment0 in
// closed form (clamped to [0,1]).
function bruteForceSqrDistance(s0: Segment, s1: Segment,
    samples: number): number {
    const dir0 = sub(s0.p[1], s0.p[0]);
    const dir1 = sub(s1.p[1], s1.p[0]);
    const a = dot(dir0, dir0);
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= samples; ++i) {
        const q = add(s1.p[0], mul(i / samples, dir1));
        let t = 0;
        if (a > 0) {
            t = dot(dir0, sub(q, s0.p[0])) / a;
            t = Math.min(Math.max(t, 0), 1);
        }
        const p = add(s0.p[0], mul(t, dir0));
        const d = sub(p, q);
        best = Math.min(best, dot(d, d));
    }
    return best;
}

describe('DistSegmentSegment', () => {
    const query = new DistSegmentSegment();

    // The reported closest points must lie on the segments (parameters in
    // [0,1]) and realize the reported distance.
    function verify(result: DistSegmentSegmentResult, s0: Segment,
        s1: Segment): void {
        for (let i = 0; i < 2; ++i) {
            expect(result.parameter[i]).toBeGreaterThanOrEqual(0);
            expect(result.parameter[i]).toBeLessThanOrEqual(1);
        }
        const p = add(s0.p[0], mul(result.parameter[0],
            sub(s0.p[1], s0.p[0])));
        const q = add(s1.p[0], mul(result.parameter[1],
            sub(s1.p[1], s1.p[0])));
        for (let i = 0; i < p.size; ++i) {
            expect(result.closest[0].values[i]).toBeCloseTo(p.values[i], 10);
            expect(result.closest[1].values[i]).toBeCloseTo(q.values[i], 10);
        }
        const diff = sub(result.closest[0], result.closest[1]);
        expect(result.sqrDistance).toBeCloseTo(dot(diff, diff), 10);
        expect(result.distance).toBeCloseTo(Math.sqrt(result.sqrDistance), 12);
    }

    // Run both the standard and the robust query on the same input and check
    // that they agree, then return the standard result.
    function computeBoth(s0: Segment, s1: Segment): DistSegmentSegmentResult {
        const result = query.compute(s0, s1);
        const robust = query.computeRobust(s0, s1);
        verify(result, s0, s1);
        verify(robust, s0, s1);
        expect(robust.distance).toBeCloseTo(result.distance, 10);
        // The endpoint overloads must match the Segment overloads exactly.
        const byPoints = query.computeEndpoints(s0.p[0], s0.p[1], s1.p[0],
            s1.p[1]);
        expect(byPoints.parameter).toEqual(result.parameter);
        const robustByPoints = query.computeRobustEndpoints(s0.p[0], s0.p[1],
            s1.p[0], s1.p[1]);
        expect(robustByPoints.parameter).toEqual(robust.parameter);
        return result;
    }

    it('finds interior closest points of skew segments (region 0)', () => {
        // The x-axis segment and a segment along y at x = 1, z = 2.
        const s0 = seg([0, 0, 0], [4, 0, 0]);
        const s1 = seg([1, -1, 2], [1, 1, 2]);
        const result = computeBoth(s0, s1);
        expect(result.parameter[0]).toBeCloseTo(0.25, 12);
        expect(result.parameter[1]).toBeCloseTo(0.5, 12);
        expect(result.distance).toBeCloseTo(2, 12);
    });

    it('reports zero distance for crossing segments', () => {
        const s0 = seg([-1, 0, 0], [1, 0, 0]);
        const s1 = seg([0, -1, 0], [0, 1, 0]);
        const result = computeBoth(s0, s1);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter[0]).toBeCloseTo(0.5, 12);
        expect(result.parameter[1]).toBeCloseTo(0.5, 12);
    });

    it('finds the closest pair of endpoints', () => {
        // Two collinear-in-y segments offset in x; the closest pair is
        // (4,0,0) and (6,3,0).
        const s0 = seg([0, 0, 0], [4, 0, 0]);
        const s1 = seg([6, 3, 0], [10, 3, 0]);
        const result = computeBoth(s0, s1);
        expect(result.parameter[0]).toBe(1);
        expect(result.parameter[1]).toBe(0);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 9), 12);
    });

    it('handles parallel segments with a partial overlap', () => {
        const s0 = seg([0, 0, 0], [4, 0, 0]);
        const s1 = seg([2, 5, 0], [8, 5, 0]);
        const result = computeBoth(s0, s1);
        expect(result.distance).toBeCloseTo(5, 12);
    });

    it('handles antiparallel segments with a partial overlap', () => {
        const s0 = seg([0, 0, 0], [4, 0, 0]);
        const s1 = seg([8, 5, 0], [2, 5, 0]);
        const result = computeBoth(s0, s1);
        expect(result.distance).toBeCloseTo(5, 12);
    });

    it('handles collinear overlapping segments (zero distance)', () => {
        const s0 = seg([0, 0, 0], [4, 0, 0]);
        const s1 = seg([2, 0, 0], [8, 0, 0]);
        const result = computeBoth(s0, s1);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('handles collinear disjoint segments', () => {
        const s0 = seg([0, 0, 0], [4, 0, 0]);
        const s1 = seg([7, 0, 0], [9, 0, 0]);
        const result = computeBoth(s0, s1);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.parameter[0]).toBe(1);
        expect(result.parameter[1]).toBe(0);
    });

    it('handles one degenerate (zero-length) segment', () => {
        const s0 = seg([0, 0, 0], [4, 0, 0]);
        const point = seg([1, 3, 0], [1, 3, 0]);
        const result = computeBoth(s0, point);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.parameter[0]).toBeCloseTo(0.25, 12);

        // The other argument order.
        const swapped = computeBoth(point, s0);
        expect(swapped.distance).toBeCloseTo(3, 12);
        expect(swapped.parameter[1]).toBeCloseTo(0.25, 12);
    });

    it('handles a degenerate segment beyond an endpoint', () => {
        const s0 = seg([0, 0, 0], [4, 0, 0]);
        const point = seg([7, 3, 0], [7, 3, 0]);
        const result = computeBoth(s0, point);
        expect(result.parameter[0]).toBe(1);
        expect(result.distance).toBeCloseTo(Math.sqrt(9 + 9), 12);
    });

    it('handles two degenerate segments', () => {
        const p = seg([1, 2, 3], [1, 2, 3]);
        const q = seg([1, 2, 8], [1, 2, 8]);
        const result = computeBoth(p, q);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.parameter[0]).toBe(0);
        expect(result.parameter[1]).toBe(0);
    });

    it('is symmetric in its arguments', () => {
        const cases: Array<[Segment, Segment]> = [
            [seg([0, 0, 0], [4, 0, 0]), seg([1, -1, 2], [1, 1, 2])],
            [seg([0, 0, 0], [4, 0, 0]), seg([6, 3, 0], [10, 3, 0])],
            [seg([0, 0, 0], [4, 0, 0]), seg([2, 5, 0], [8, 5, 0])],
            [seg([1, 2, 3], [-2, 0, 1]), seg([4, -1, 2], [0, 3, -3])]
        ];
        for (const [s0, s1] of cases) {
            const forward = computeBoth(s0, s1);
            const reverse = computeBoth(s1, s0);
            expect(reverse.distance).toBeCloseTo(forward.distance, 10);
        }
    });

    it('works in 2D and 4D', () => {
        const a = computeBoth(seg([0, 0], [4, 0]), seg([1, 6], [3, 6]));
        expect(a.distance).toBeCloseTo(6, 12);

        const b = computeBoth(seg([0, 0, 0, 0], [4, 0, 0, 0]),
            seg([2, 3, 4, 0], [2, 3, 4, 0]));
        expect(b.distance).toBeCloseTo(5, 12);
    });

    it('handles nearly parallel segments robustly', () => {
        // The segments differ in direction by about 1e-9 radians, which is
        // the case the conjugate-gradient (robust) query is designed for.
        const eps = 1e-9;
        const s0 = seg([0, 0, 0], [1, 0, 0]);
        const s1 = seg([0, 1, 0], [1, eps, 0]);
        const result = query.compute(s0, s1);
        const robust = query.computeRobust(s0, s1);
        verify(result, s0, s1);
        verify(robust, s0, s1);
        const best = bruteForceSqrDistance(s0, s1, 200000);
        expect(robust.sqrDistance).toBeLessThanOrEqual(best + 1e-12);
        expect(Math.sqrt(best) - robust.distance).toBeLessThan(1e-6);
        expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-12);
    });

    it('is more accurate than compute for near-parallel segments', () => {
        // Upstream documents computeRobust (ComputeRobust) as the
        // floating-point-robust query: the two segments here are parallel to
        // within rounding, so a*c - b*b is a tiny positive number and the
        // divisions in compute lose most of their significant digits. The
        // robust query matches a dense brute-force minimum; compute is off
        // by a factor of more than two in the squared distance.
        const s0 = seg([-0.2821011543273926, 1.5676589012145996,
            -1.549705982208252],
            [0.2521996389198904, 0.048745550542832916, 0.9485857130375877]);
        const s1 = seg([-0.6344702839851379, 1.122955083847046,
            0.47550272941589355],
            [-1.2362647284939925, 2.833739772316619, -2.3383767968039564]);
        const result = query.compute(s0, s1);
        const robust = query.computeRobust(s0, s1);
        verify(result, s0, s1);
        verify(robust, s0, s1);

        const best = bruteForceSqrDistance(s0, s1, 200000);
        expect(robust.sqrDistance).toBeCloseTo(best, 12);
        expect(robust.sqrDistance).toBeCloseTo(0.9406708459521783, 12);
        // The nonrobust query overestimates here; the port preserves that
        // upstream behavior rather than "fixing" it.
        expect(result.sqrDistance).toBeGreaterThan(best + 1);
        expect(result.sqrDistance).toBeCloseTo(2.61441244705322, 10);
    });

    it('agrees with a brute-force minimum for random segments', () => {
        let seed = 99991;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const s0 = seg([rand(), rand(), rand()],
                [rand(), rand(), rand()]);
            const s1 = seg([rand(), rand(), rand()],
                [rand(), rand(), rand()]);
            const result = query.compute(s0, s1);
            const robust = query.computeRobust(s0, s1);
            verify(result, s0, s1);
            verify(robust, s0, s1);

            const best = bruteForceSqrDistance(s0, s1, 20000);
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
            expect(robust.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
            expect(Math.sqrt(best) - result.distance).toBeLessThan(1e-5);
            expect(Math.sqrt(best) - robust.distance).toBeLessThan(1e-5);
        }
    });

    it('agrees with a brute-force minimum for parallel segments', () => {
        let seed = 5150;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const dir = [rand(), rand(), rand()];
            const p0 = [rand(), rand(), rand()];
            const q0 = [rand(), rand(), rand()];
            const scaleP = rand();
            const scaleQ = rand();
            const s0 = seg(p0, [p0[0] + scaleP * dir[0],
                p0[1] + scaleP * dir[1], p0[2] + scaleP * dir[2]]);
            const s1 = seg(q0, [q0[0] + scaleQ * dir[0],
                q0[1] + scaleQ * dir[1], q0[2] + scaleQ * dir[2]]);
            const result = query.compute(s0, s1);
            const robust = query.computeRobust(s0, s1);
            verify(result, s0, s1);
            verify(robust, s0, s1);

            const best = bruteForceSqrDistance(s0, s1, 20000);
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
            expect(robust.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
            expect(Math.sqrt(best) - result.distance).toBeLessThan(1e-5);
            expect(Math.sqrt(best) - robust.distance).toBeLessThan(1e-5);
        }
    });
});
