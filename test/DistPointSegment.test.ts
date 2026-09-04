import { describe, expect, it } from 'vitest';
import { DistPointSegment } from '../src/DistPointSegment.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { DistPointLine } from '../src/DistPointLine.js';
import { Line } from '../src/Line.js';
import { check, expectClose, expectVectorClose, fc, rotationFrame, seededRandom, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

describe('DistPointSegment', () => {
    const query = new DistPointSegment();

    it('uses an interior point when the projection is inside [0,1]', () => {
        const result = query.compute(v(2, 3), seg([0, 0], [4, 0]));
        expect(result.parameter).toBeCloseTo(0.5, 12);
        expect(result.closest[1].values).toEqual([2, 0]);
        expect(result.distance).toBeCloseTo(3, 12);
    });

    it('clamps to the second endpoint', () => {
        const result = query.compute(v(9, 3), seg([0, 0], [4, 0]));
        expect(result.parameter).toBe(1);
        expect(result.closest[1].values).toEqual([4, 0]);
        expect(result.distance).toBeCloseTo(Math.sqrt(25 + 9), 12);
    });

    it('clamps to the first endpoint', () => {
        const result = query.compute(v(-3, 4), seg([0, 0], [4, 0]));
        expect(result.parameter).toBe(0);
        expect(result.closest[1].values).toEqual([0, 0]);
        expect(result.distance).toBeCloseTo(5, 12);
    });

    it('handles a degenerate (zero-length) segment', () => {
        const result = query.compute(v(3, 4), seg([1, 1], [1, 1]));
        // Dot(direction, point - p1) = 0 >= 0, so the query reports t = 1.
        expect(result.parameter).toBe(1);
        expect(result.closest[1].values).toEqual([1, 1]);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 9), 12);
    });

    it('reports zero distance for a point on the segment', () => {
        const s = seg([1, 2, 3], [4, 0, -1]);
        const d = sub(s.p[1], s.p[0]);
        const point = add(s.p[0], mul(0.25, d));
        const result = query.compute(point, s);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter).toBeCloseTo(0.25, 12);
    });

    it('agrees with a sampled minimum over the segment', () => {
        let seed = 4242;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const point = v(rand(), rand(), rand());
            const s = seg([rand(), rand(), rand()], [rand(), rand(), rand()]);
            const result = query.compute(point, s);

            expect(result.parameter).toBeGreaterThanOrEqual(0);
            expect(result.parameter).toBeLessThanOrEqual(1);

            const dir = sub(s.p[1], s.p[0]);
            const onSeg = add(s.p[0], mul(result.parameter, dir));
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[1].values[i]).toBeCloseTo(
                    onSeg.values[i], 9);
            }

            let best = Number.MAX_VALUE;
            for (let k = 0; k <= 2000; ++k) {
                const t = k / 2000;
                const q = add(s.p[0], mul(t, dir));
                const diff = sub(point, q);
                best = Math.min(best, dot(diff, diff));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of DistPointSegment.ts
// against the upstream header DistPointSegment.h.
// ---------------------------------------------------------------------------

function rot(R: readonly Vector[], p: Vector): Vector {
    let q = mul(p.values[0], R[0]);
    for (let i = 1; i < R.length; ++i) {
        q = add(q, mul(p.values[i], R[i]));
    }
    return q;
}

const segment3 = fc.tuple(wellScaledVector(3, -8, 8), wellScaledVector(3, -8, 8))
    .filter(([p0, p1]) => length(sub(p1, p0)) > 0.25)
    .map(([p0, p1]) => Segment.fromEndpoints(p0, p1));

describe('DistPointSegment verification', () => {
    const query = new DistPointSegment();
    const lineQuery = new DistPointLine();

    it('result is self consistent and the parameter is in [0,1]', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), segment3), ([p, seg]) => {
            const r = query.compute(p, seg);
            expect(r.parameter).toBeGreaterThanOrEqual(0);
            expect(r.parameter).toBeLessThanOrEqual(1);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            const diff = sub(r.closest[0], r.closest[1]);
            expectClose(r.distance, Math.sqrt(dot(diff, diff)), 1e-12, 1e-12);
            expectVectorClose(r.closest[0], p, 0, 0);
            expect(r.closest[0]).not.toBe(p);
            expectVectorClose(r.closest[1],
                add(seg.p[0], mul(r.parameter, sub(seg.p[1], seg.p[0]))),
                1e-12, 1e-12);
        });
    });

    it('agrees with the point-line query when the projection is interior',
        () => {
            check(fc.tuple(wellScaledVector(3, -8, 8), segment3), ([p, seg]) => {
                const dir = sub(seg.p[1], seg.p[0]);
                const rl = lineQuery.compute(p,
                    Line.fromOriginDirection(seg.p[0], dir));
                const rs = query.compute(p, seg);
                if (rl.parameter > 1e-6 && rl.parameter < 1 - 1e-6) {
                    expectClose(rs.parameter, rl.parameter, 1e-9, 1e-9);
                    expectClose(rs.distance, rl.distance, 1e-9, 1e-9);
                }
                else {
                    // Outside the domain the segment answer is an endpoint.
                    expect(rs.distance).toBeGreaterThanOrEqual(
                        rl.distance - 1e-9 * (1 + rl.distance));
                }
            });
        });

    it('is minimal over sampled segment points', () => {
        const rand = seededRandom(0x51d2);
        check(fc.tuple(wellScaledVector(3, -8, 8), segment3), ([p, seg]) => {
            const r = query.compute(p, seg);
            const dir = sub(seg.p[1], seg.p[0]);
            for (let k = 0; k < 20; ++k) {
                const t = rand();
                const diff = sub(p, add(seg.p[0], mul(t, dir)));
                const sqr = dot(diff, diff);
                expect(r.sqrDistance).toBeLessThanOrEqual(sqr + 1e-9 * (1 + sqr));
            }
        }, 60);
    });

    it('reverses the parameter when the endpoints are swapped', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), segment3), ([p, seg]) => {
            const r0 = query.compute(p, seg);
            const r1 = query.compute(p,
                Segment.fromEndpoints(seg.p[1], seg.p[0]));
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expectClose(r0.parameter, 1 - r1.parameter, 1e-9, 1e-9);
            expectVectorClose(r0.closest[1], r1.closest[1], 1e-9, 1e-9);
        });
    });

    it('is invariant under rigid motions', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), segment3, rotationFrame(3),
            wellScaledVector(3, -5, 5)), ([p, seg, R, tr]) => {
            const r0 = query.compute(p, seg);
            const r1 = query.compute(add(rot(R, p), tr),
                Segment.fromEndpoints(add(rot(R, seg.p[0]), tr),
                    add(rot(R, seg.p[1]), tr)));
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expectClose(r0.parameter, r1.parameter, 1e-9, 1e-8);
        });
    });

    it('handles a zero-length segment', () => {
        // With P0 = P1 the direction is zero, so Dot(D, P-P1) = 0 >= 0 and
        // upstream takes the first branch: parameter 1, closest point P1.
        check(fc.tuple(wellScaledVector(3, -8, 8), wellScaledVector(3, -8, 8)),
            ([p, q]) => {
                const r = query.compute(p, Segment.fromEndpoints(q, q));
                expect(r.parameter).toBe(1);
                expectVectorClose(r.closest[1], q, 0, 0);
                const diff = sub(p, q);
                expectClose(r.sqrDistance, dot(diff, diff), 1e-12, 1e-12);
            });
    });
});
