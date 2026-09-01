import { describe, it, expect } from 'vitest';
import { Segment } from '../src/Segment';
import { Vector, add, dot, mul, sub } from '../src/Vector';
import { dotPerp } from '../src/Vector2';
import {
    IntrSegment2Segment2TI,
    IntrSegment2Segment2FI
} from '../src/IntrSegment2Segment2';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Independent brute-force test: does any sampled point of segment0 lie within
// tolerance of segment1?
function samplesOverlap(s0: Segment, s1: Segment, tol: number): boolean {
    const d1 = sub(s1.p[1], s1.p[0]);
    const dotD1D1 = dot(d1, d1);
    const samples = 2000;
    for (let k = 0; k <= samples; ++k) {
        const p = add(s0.p[0], mul(k / samples, sub(s0.p[1], s0.p[0])));
        let u = dotD1D1 > 0 ? dot(sub(p, s1.p[0]), d1) / dotD1D1 : 0;
        u = Math.min(1, Math.max(0, u));
        const q = add(s1.p[0], mul(u, d1));
        const diff = sub(p, q);
        if (Math.sqrt(dot(diff, diff)) <= tol) {
            return true;
        }
    }
    return false;
}

describe('IntrSegment2Segment2', () => {
    const ti = new IntrSegment2Segment2TI();
    const fi = new IntrSegment2Segment2FI();

    it('finds a transversal crossing at a single point', () => {
        const s0 = segment([-2, 0], [2, 0]);
        const s1 = segment([1, -3], [1, 3]);
        expect(ti.test(s0, s1).numIntersections).toBe(1);
        const result = fi.find(s0, s1);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
        // The centered form of s0 has center (0,0), direction (1,0), extent 2.
        expect(result.segment0Parameter[0]).toBeCloseTo(1, 12);
        // The centered form of s1 has center (1,0), direction (0,1), extent 3.
        expect(result.segment1Parameter[0]).toBeCloseTo(0, 12);

        // The 'exact' variants use the endpoint parameterization.
        const exact = fi.findExact(s0, s1);
        expect(exact.numIntersections).toBe(1);
        expect(exact.segment0Parameter[0]).toBeCloseTo(0.75, 12);
        expect(exact.segment1Parameter[0]).toBeCloseTo(0.5, 12);
        expect(ti.testExact(s0, s1).numIntersections).toBe(1);
    });

    it('reports a touching endpoint as a single intersection', () => {
        const s0 = segment([0, 0], [2, 0]);
        const s1 = segment([2, 0], [2, 4]);
        expect(ti.test(s0, s1).intersect).toBe(true);
        const result = fi.find(s0, s1);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(2, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
        expect(fi.findExact(s0, s1).segment0Parameter[0]).toBeCloseTo(1, 12);
    });

    it('misses when the crossing lies off both segments', () => {
        const s0 = segment([0, 0], [1, 0]);
        const s1 = segment([5, -1], [5, 1]);
        expect(ti.test(s0, s1).intersect).toBe(false);
        expect(ti.testExact(s0, s1).intersect).toBe(false);
        const result = fi.find(s0, s1);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(fi.findExact(s0, s1).intersect).toBe(false);
    });

    it('reports collinear overlapping segments as a segment', () => {
        const s0 = segment([0, 0], [4, 0]);
        const s1 = segment([2, 0], [6, 0]);
        const tiResult = ti.test(s0, s1);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(2);

        const result = fi.find(s0, s1);
        expect(result.numIntersections).toBe(2);
        // The overlap is x in [2,4]; s0's centered form is center (2,0),
        // direction (1,0), extent 2.
        expect(result.segment0Parameter[0]).toBeCloseTo(0, 12);
        expect(result.segment0Parameter[1]).toBeCloseTo(2, 12);
        expect(result.point[0].values[0]).toBeCloseTo(2, 12);
        expect(result.point[1].values[0]).toBeCloseTo(4, 12);

        const exact = fi.findExact(s0, s1);
        expect(exact.numIntersections).toBe(2);
        expect(exact.segment0Parameter[0]).toBeCloseTo(0.5, 12);
        expect(exact.segment0Parameter[1]).toBeCloseTo(1, 12);
        expect(exact.segment1Parameter[0]).toBeCloseTo(0, 12);
        expect(exact.segment1Parameter[1]).toBeCloseTo(0.5, 12);
    });

    it('reports collinear segments touching at one endpoint as a point', () => {
        const s0 = segment([0, 0], [2, 0]);
        const s1 = segment([2, 0], [5, 0]);
        const tiResult = ti.test(s0, s1);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(1);

        const result = fi.find(s0, s1);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(2, 12);
        expect(result.point[1].values[0]).toBeCloseTo(2, 12);
        expect(result.segment0Parameter[0]).toBe(result.segment0Parameter[1]);
    });

    it('reports disjoint collinear segments as no intersection', () => {
        const s0 = segment([0, 0], [1, 0]);
        const s1 = segment([3, 0], [5, 0]);
        expect(ti.test(s0, s1).intersect).toBe(false);
        expect(fi.find(s0, s1).intersect).toBe(false);
        expect(ti.testExact(s0, s1).intersect).toBe(false);
        expect(fi.findExact(s0, s1).intersect).toBe(false);
    });

    it('reports parallel distinct segments as no intersection', () => {
        const s0 = segment([0, 0], [4, 0]);
        const s1 = segment([0, 1], [4, 1]);
        expect(ti.test(s0, s1).intersect).toBe(false);
        expect(fi.find(s0, s1).intersect).toBe(false);
        expect(ti.testExact(s0, s1).intersect).toBe(false);
        expect(fi.findExact(s0, s1).intersect).toBe(false);
    });

    it('swaps segment1 parameters for a reversed collinear overlap (exact)', () => {
        const s0 = segment([0, 0], [4, 0]);
        // s1 runs in the opposite direction, from x = 6 back to x = 2.
        const s1 = segment([6, 0], [2, 0]);
        const exact = fi.findExact(s0, s1);
        expect(exact.numIntersections).toBe(2);
        expect(exact.segment0Parameter[0]).toBeCloseTo(0.5, 12);
        expect(exact.segment0Parameter[1]).toBeCloseTo(1, 12);
        // Upstream swaps so segment1Parameter is increasing.
        expect(exact.segment1Parameter[0])
            .toBeLessThan(exact.segment1Parameter[1]);
        expect(exact.segment1Parameter[0]).toBeCloseTo(0.5, 12);
        expect(exact.segment1Parameter[1]).toBeCloseTo(1, 12);
    });

    it('handles a degenerate segment on the other segment (exact query)', () => {
        // A zero-length segment cannot be handled by the centered-form query
        // (Normalize produces a zero direction), but the exact query treats it
        // as the point P0 and reports containment through the collinear path.
        const point = segment([1, 0], [1, 0]);
        const s = segment([0, 0], [4, 0]);
        const exact = fi.findExact(s, point);
        expect(exact.intersect).toBe(true);
        expect(exact.numIntersections).toBe(1);
        expect(exact.point[0].values[0]).toBeCloseTo(1, 12);

        // Upstream quirk (preserved): a zero-length segment has a zero
        // direction, so IntrLine2Line2 reports DotPerp(D0,D1) = 0 and
        // DotPerp(Q,D1) = 0, i.e. "the lines are the same", for ANY position
        // of the degenerate segment. The collinear branch then projects the
        // point onto segment0 and reports a spurious intersection at the
        // projection. Callers must reject degenerate segments themselves.
        const off = segment([1, 1], [1, 1]);
        const spurious = fi.findExact(s, off);
        expect(spurious.intersect).toBe(true);
        expect(spurious.point[0].values[0]).toBeCloseTo(1, 12);
        expect(spurious.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('agrees with brute-force sampling on random configurations', () => {
        const rnd = makeRandom(4242);
        let sampleMismatch = 0;
        let tiFiMismatch = 0;
        let exactAgreementMismatch = 0;
        let pointMismatch = 0;
        let hits = 0;

        for (let trial = 0; trial < 400; ++trial) {
            const s0 = Segment.fromEndpoints(
                vec(4 * rnd() - 2, 4 * rnd() - 2),
                vec(4 * rnd() - 2, 4 * rnd() - 2));
            const s1 = Segment.fromEndpoints(
                vec(4 * rnd() - 2, 4 * rnd() - 2),
                vec(4 * rnd() - 2, 4 * rnd() - 2));
            const d0 = sub(s0.p[1], s0.p[0]);
            const d1 = sub(s1.p[1], s1.p[0]);
            if (dot(d0, d0) < 1e-4 || dot(d1, d1) < 1e-4) {
                continue;
            }
            // Skip near-parallel pairs, where the sampling tolerance and the
            // exact zero test for dotPerp disagree.
            if (Math.abs(dotPerp(d0, d1)) < 1e-3) {
                continue;
            }

            const tiResult = ti.test(s0, s1);
            const fiResult = fi.find(s0, s1);
            if (tiResult.intersect !== fiResult.intersect ||
                tiResult.numIntersections !== fiResult.numIntersections) {
                ++tiFiMismatch;
            }

            const exactResult = fi.findExact(s0, s1);
            if (exactResult.intersect !== fiResult.intersect) {
                ++exactAgreementMismatch;
            }

            const sampled = samplesOverlap(s0, s1, 2e-3);
            if (sampled !== fiResult.intersect) {
                // Sampling can only miss an intersection near a shared
                // endpoint; a false positive here would be a real bug.
                if (!fiResult.intersect) {
                    ++sampleMismatch;
                }
            }

            if (fiResult.intersect) {
                ++hits;
                const cf0 = s0.getCenteredForm();
                const cf1 = s1.getCenteredForm();
                for (let i = 0; i < fiResult.numIntersections; ++i) {
                    const q0 = add(cf0.center,
                        mul(fiResult.segment0Parameter[i], cf0.direction));
                    const q1 = add(cf1.center,
                        mul(fiResult.segment1Parameter[i], cf1.direction));
                    const dq = sub(q0, q1);
                    if (Math.sqrt(dot(dq, dq)) > 1e-6) {
                        ++pointMismatch;
                    }
                    const dp = sub(fiResult.point[i], q0);
                    if (Math.sqrt(dot(dp, dp)) > 1e-6) {
                        ++pointMismatch;
                    }
                    // The exact query's point must be the same location.
                    if (exactResult.intersect) {
                        const de = sub(exactResult.point[i], fiResult.point[i]);
                        if (Math.sqrt(dot(de, de)) > 1e-6) {
                            ++pointMismatch;
                        }
                    }
                }
            }
        }

        expect(hits).toBeGreaterThan(20);
        expect([tiFiMismatch, exactAgreementMismatch, sampleMismatch,
            pointMismatch]).toEqual([0, 0, 0, 0]);
    });
});
