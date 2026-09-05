import { describe, it, expect } from 'vitest';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { dotPerp } from '../src/Vector2.js';
import {
    IntrSegment2Segment2TI,
    IntrSegment2Segment2FI
} from '../src/IntrSegment2Segment2.js';
import {
    check, expectClose, expectVectorClose, fc, latticeVector, wellScaledVector
} from './helpers/arbitraries.js';
import { orient2 } from './helpers/exact.js';

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

// ---------------------------------------------------------------------------
// Verification (V32): properties cross-checking the port against upstream
// IntrSegment2Segment2.h.
// ---------------------------------------------------------------------------

describe('IntrSegment2Segment2 verification', () => {
    const ti = new IntrSegment2Segment2TI();
    const fi = new IntrSegment2Segment2FI();

    // Integer-lattice segments: small integers are exact in binary64, so the
    // DotPerp parallelism tests of the 'Exact' queries and the bigint
    // reference predicate below agree bit for bit.
    const arbLatticeSegment = fc.tuple(latticeVector(2, -6, 6),
        latticeVector(2, -6, 6))
        .filter(([a, b]) => !a.equals(b))
        .map(([a, b]) => Segment.fromEndpoints(a, b));
    const arbLatticePair = fc.tuple(arbLatticeSegment, arbLatticeSegment);

    function big(s: Segment, i: number, d: number): bigint {
        return BigInt(s.p[i].get(d));
    }

    // r is on the closed segment pq, given that (p,q,r) are collinear.
    function onSegment(px: bigint, py: bigint, qx: bigint, qy: bigint,
        rx: bigint, ry: bigint): boolean {
        return orient2(px, py, qx, qy, rx, ry) === 0
            && rx >= (px < qx ? px : qx) && rx <= (px > qx ? px : qx)
            && ry >= (py < qy ? py : qy) && ry <= (py > qy ? py : qy);
    }

    // The exact closed-segment intersection predicate on integer coordinates.
    function exactIntersect(s0: Segment, s1: Segment): boolean {
        const ax = big(s0, 0, 0), ay = big(s0, 0, 1);
        const bx = big(s0, 1, 0), by = big(s0, 1, 1);
        const cx = big(s1, 0, 0), cy = big(s1, 0, 1);
        const dx = big(s1, 1, 0), dy = big(s1, 1, 1);
        const o1 = orient2(ax, ay, bx, by, cx, cy);
        const o2 = orient2(ax, ay, bx, by, dx, dy);
        const o3 = orient2(cx, cy, dx, dy, ax, ay);
        const o4 = orient2(cx, cy, dx, dy, bx, by);
        if (o1 !== 0 && o2 !== 0 && o3 !== 0 && o4 !== 0) {
            return o1 !== o2 && o3 !== o4;
        }
        return onSegment(ax, ay, bx, by, cx, cy)
            || onSegment(ax, ay, bx, by, dx, dy)
            || onSegment(cx, cy, dx, dy, ax, ay)
            || onSegment(cx, cy, dx, dy, bx, by);
    }

    it('testExact matches the exact bigint segment-intersection predicate',
        () => {
            check(arbLatticePair, ([s0, s1]) => {
                expect(ti.testExact(s0, s1).intersect)
                    .toBe(exactIntersect(s0, s1));
            });
        });

    it('findExact agrees with testExact and matches the exact predicate', () => {
        check(arbLatticePair, ([s0, s1]) => {
            const res = fi.findExact(s0, s1);
            const tiRes = ti.testExact(s0, s1);
            expect(res.intersect).toBe(exactIntersect(s0, s1));
            expect(res.intersect).toBe(tiRes.intersect);
            expect(res.numIntersections).toBe(tiRes.numIntersections);
        });
    });

    it('findExact points lie on both segments with parameters in [0,1]', () => {
        check(arbLatticePair, ([s0, s1]) => {
            const res = fi.findExact(s0, s1);
            if (!res.intersect) {
                return;
            }
            const d0 = sub(s0.p[1], s0.p[0]);
            const d1 = sub(s1.p[1], s1.p[0]);
            for (let i = 0; i < res.numIntersections; ++i) {
                const t0 = res.segment0Parameter[i];
                const t1 = res.segment1Parameter[i];
                expect(t0).toBeGreaterThanOrEqual(-1e-12);
                expect(t0).toBeLessThanOrEqual(1 + 1e-12);
                expect(t1).toBeGreaterThanOrEqual(-1e-12);
                expect(t1).toBeLessThanOrEqual(1 + 1e-12);
                expectVectorClose(res.point[i], add(s0.p[0], mul(t0, d0)),
                    1e-9, 1e-9);
                // Upstream's 'Exact' query swaps segment1Parameter[0] and
                // [1] when the collinear segments point in opposite
                // directions, in order to keep segment1Parameter ascending.
                // The parameter therefore names one of the two reported
                // points, not necessarily point[i].
                const q = add(s1.p[0], mul(t1, d1));
                const other = res.point[res.numIntersections - 1 - i];
                const matches = (a: Vector, b: Vector): boolean => {
                    const e = sub(a, b);
                    return Math.sqrt(dot(e, e)) <= 1e-9 * (1 + Math.sqrt(
                        dot(a, a)));
                };
                expect(matches(q, res.point[i]) || matches(q, other))
                    .toBe(true);
            }
            if (res.numIntersections === 2) {
                expect(res.segment0Parameter[0])
                    .toBeLessThanOrEqual(res.segment0Parameter[1]);
                expect(res.segment1Parameter[0])
                    .toBeLessThanOrEqual(res.segment1Parameter[1]);
            }
        });
    });

    it('the centered-form find keeps point[i] = C1 + segment1Parameter[i]*D1'
        + ' for antiparallel collinear segments', () => {
        // Upstream computes segment1Parameter[i] = overlap[i] - t, which has
        // the wrong sign when the two centered directions are opposite: for
        // these inputs it reports -2.5 and 0.5, naming (6,0) and (3,0),
        // neither of which is an intersection point. The port negates in the
        // antiparallel case.
        const s0 = segment([0, 0], [4, 0]);
        const s1 = segment([6, 0], [1, 0]);
        const res = fi.find(s0, s1);
        expect(res.numIntersections).toBe(2);
        expectVectorClose(res.point[0], vec(1, 0));
        expectVectorClose(res.point[1], vec(4, 0));
        expect(res.segment1Parameter[0]).toBe(2.5);
        expect(res.segment1Parameter[1]).toBe(-0.5);
        const c1 = s1.getCenteredForm();
        for (let i = 0; i < 2; ++i) {
            expectVectorClose(
                add(c1.center, mul(res.segment1Parameter[i], c1.direction)),
                res.point[i]);
        }
    });

    it('the centered-form find matches upstream for parallel collinear'
        + ' segments with the same direction', () => {
        const s0 = segment([0, 0], [4, 0]);
        const s1 = segment([1, 0], [6, 0]);
        const res = fi.find(s0, s1);
        expect(res.numIntersections).toBe(2);
        const c0 = s0.getCenteredForm(), c1 = s1.getCenteredForm();
        const t = dot(c0.direction, sub(c1.center, c0.center));
        // The unmodified upstream expression.
        expect(res.segment1Parameter[0])
            .toBe(res.segment0Parameter[0] - t);
        expect(res.segment1Parameter[1])
            .toBe(res.segment0Parameter[1] - t);
    });

    it('the queries are symmetric under swapping the two segments', () => {
        check(arbLatticePair, ([s0, s1]) => {
            const a = ti.testExact(s0, s1), b = ti.testExact(s1, s0);
            expect(b.intersect).toBe(a.intersect);
            expect(b.numIntersections).toBe(a.numIntersections);
            const c = fi.findExact(s0, s1), d = fi.findExact(s1, s0);
            expect(d.intersect).toBe(c.intersect);
            expect(d.numIntersections).toBe(c.numIntersections);
            const e = ti.test(s0, s1), f = ti.test(s1, s0);
            expect(f.intersect).toBe(e.intersect);
            expect(f.numIntersections).toBe(e.numIntersections);
        });
    });

    it('the centered-form find reports points on both segments', () => {
        const arbWellScaled = fc.tuple(wellScaledVector(2, -8, 8),
            wellScaledVector(2, -8, 8))
            .filter(([a, b]) => {
                const d = sub(b, a);
                return dot(d, d) > 1e-2;
            })
            .map(([a, b]) => Segment.fromEndpoints(a, b));
        check(fc.tuple(arbWellScaled, arbWellScaled), ([s0, s1]) => {
            const res = fi.find(s0, s1);
            expect(res.intersect).toBe(ti.test(s0, s1).intersect);
            if (!res.intersect) {
                return;
            }
            const c0 = s0.getCenteredForm(), c1 = s1.getCenteredForm();
            for (let i = 0; i < res.numIntersections; ++i) {
                const t0 = res.segment0Parameter[i];
                const t1 = res.segment1Parameter[i];
                expect(Math.abs(t0)).toBeLessThanOrEqual(
                    c0.extent * (1 + 1e-9) + 1e-9);
                expect(Math.abs(t1)).toBeLessThanOrEqual(
                    c1.extent * (1 + 1e-9) + 1e-9);
                // The reported point uses the centered form of segment0.
                expectVectorClose(res.point[i],
                    add(c0.center, mul(t0, c0.direction)), 1e-9, 1e-9);
                expect(Number.isNaN(res.point[i].get(0))).toBe(false);
                expect(Number.isNaN(res.point[i].get(1))).toBe(false);
            }
        });
    });

    it('point[1] is a copy of point[0], not an alias', () => {
        // Upstream copies the Vector2 into point[1]; a plain TS assignment
        // aliases, so a caller mutating point[0] would silently change
        // point[1].
        const s0 = segment([-1, 0], [1, 0]);
        const s1 = segment([0, -1], [0, 1]);
        const res = fi.find(s0, s1);
        expect(res.numIntersections).toBe(1);
        res.point[0].set(0, 42);
        expect(res.point[1].get(0)).toBe(0);

        // Collinear single-point overlap: [0,4] and [4,8] share only x = 4.
        const e0 = segment([0, 0], [4, 0]);
        const e1 = segment([4, 0], [8, 0]);
        const resE = fi.findExact(e0, e1);
        expect(resE.numIntersections).toBe(1);
        resE.point[0].set(0, 42);
        expect(resE.point[1].get(0)).toBe(4);
    });

    it('reports the overlapping sub-segment for collinear segments', () => {
        const s0 = segment([0, 0], [4, 0]);
        const s1 = segment([1, 0], [6, 0]);
        const res = fi.findExact(s0, s1);
        expect(res.intersect).toBe(true);
        expect(res.numIntersections).toBe(2);
        expectVectorClose(res.point[0], vec(1, 0));
        expectVectorClose(res.point[1], vec(4, 0));
        expectClose(res.segment0Parameter[0], 0.25);
        expectClose(res.segment0Parameter[1], 1);
        expectClose(res.segment1Parameter[0], 0);
        expectClose(res.segment1Parameter[1], 0.6);
    });

    it('orders segment1Parameter ascending for an oppositely directed'
        + ' collinear overlap', () => {
        const s0 = segment([0, 0], [4, 0]);
        const s1 = segment([6, 0], [1, 0]);   // reversed direction
        const res = fi.findExact(s0, s1);
        expect(res.intersect).toBe(true);
        expect(res.numIntersections).toBe(2);
        expect(res.segment1Parameter[0])
            .toBeLessThanOrEqual(res.segment1Parameter[1]);
        // Both parameters evaluate to points on segment1.
        const d1 = sub(s1.p[1], s1.p[0]);
        for (let i = 0; i < 2; ++i) {
            const p = add(s1.p[0], mul(res.segment1Parameter[i], d1));
            expect(p.get(1)).toBe(0);
            expect(p.get(0)).toBeGreaterThanOrEqual(1 - 1e-12);
            expect(p.get(0)).toBeLessThanOrEqual(4 + 1e-12);
        }
    });
});
