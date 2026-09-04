import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistSegment3AlignedBox3 } from '../src/DistSegment3AlignedBox3.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, sub }
    from '../src/Vector.js';
import { DistLine3AlignedBox3 } from '../src/DistLine3AlignedBox3.js';
import { DistPointAlignedBox } from '../src/DistPointAlignedBox.js';
import { Line } from '../src/Line.js';
import {
    check, expectClose, expectVectorClose, fc, positive,
    rotationFrame, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(v(...min), v(...max));
}

// The exact squared distance from a point to the solid aligned box, computed
// independently of the library.
function pointBoxSqrDistance(p: Vector, b: AlignedBox): number {
    let sqrDistance = 0;
    for (let i = 0; i < p.size; ++i) {
        const delta = Math.max(0, b.min.values[i] - p.values[i],
            p.values[i] - b.max.values[i]);
        sqrDistance += delta * delta;
    }
    return sqrDistance;
}

// Verify that the reported box point is inside the box.
function verifyBoxPoint(b: AlignedBox, q: Vector): void {
    for (let i = 0; i < q.size; ++i) {
        expect(q.values[i]).toBeGreaterThanOrEqual(b.min.values[i] - 1e-9);
        expect(q.values[i]).toBeLessThanOrEqual(b.max.values[i] + 1e-9);
    }
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

// The distance from a segment point to the solid box is a convex function of
// the segment parameter, so a ternary search over [0,1] converges to the true
// minimum independently of the query under test.
function bruteForce(s: Segment, b: AlignedBox): number {
    const direction = sub(s.p[1], s.p[0]);
    const f = (t: number): number =>
        pointBoxSqrDistance(add(s.p[0], mul(t, direction)), b);
    let lo = 0;
    let hi = 1;
    for (let i = 0; i < 200; ++i) {
        const m0 = lo + (hi - lo) / 3;
        const m1 = hi - (hi - lo) / 3;
        if (f(m0) < f(m1)) {
            hi = m1;
        }
        else {
            lo = m0;
        }
    }
    return Math.sqrt(f(0.5 * (lo + hi)));
}

// Verify that the reported closest points are consistent with the reported
// distance and lie on their primitives.
function verifyClosest(s: Segment, b: AlignedBox,
    result: { distance: number, parameter: number, closest: [Vector, Vector] }
): void {
    expect(result.parameter).toBeGreaterThanOrEqual(0);
    expect(result.parameter).toBeLessThanOrEqual(1);
    const onSeg = add(s.p[0], mul(result.parameter, sub(s.p[1], s.p[0])));
    for (let i = 0; i < onSeg.size; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onSeg.values[i], 9);
    }
    verifyBoxPoint(b, result.closest[1]);
    expect(length(sub(result.closest[0], result.closest[1])))
        .toBeCloseTo(result.distance, 9);
}

describe('DistSegment3AlignedBox3', () => {
    const query = new DistSegment3AlignedBox3();
    const unitBox = box([-1, -1, -1], [1, 1, 1]);

    it('clamps to the first endpoint', () => {
        const s = segment([4, 0, 0], [6, 0, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.parameter).toBe(0);
        expect(result.closest[1].values).toEqual([1, 0, 0]);
    });

    it('clamps to the second endpoint', () => {
        const s = segment([6, 0, 0], [4, 0, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.parameter).toBe(1);
    });

    it('uses an interior parameter when the projection is inside', () => {
        const s = segment([-3, 3, 0], [3, 3, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.parameter).toBeGreaterThan(0);
        expect(result.parameter).toBeLessThan(1);
        verifyClosest(s, unitBox, result);
    });

    it('reports zero distance when the segment crosses the box', () => {
        const s = segment([-3, 0, 0], [3, 0, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('handles a degenerate zero-length segment', () => {
        const s = segment([4, 4, 4], [4, 4, 4]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(Math.sqrt(27), 10);
        expect(result.closest[1].values).toEqual([1, 1, 1]);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 44332211;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 150; ++trial) {
            const b = box([-1 - rand(), -1 - rand(), -1 - rand()],
                [1 + rand(), 1 + rand(), 1 + rand()]);
            const s = segment(
                [8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4]);
            const result = query.compute(s, b);
            const brute = bruteForce(s, b);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-6);
            expect(result.distance).toBeCloseTo(brute, 6);
            verifyClosest(s, b, result);
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistSegment3AlignedBox3.h. Upstream solves the query on the line containing the
// segment and clamps the result to the segment domain, so the properties below
// check the line/point agreement as well as the geometric invariants.
// ---------------------------------------------------------------------------

const v21Shape = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -8, 8)).map(([a, b]) => {
    const lo = new Vector(3);
    const hi = new Vector(3);
    for (let i = 0; i < 3; ++i) {
        lo.values[i] = Math.min(a.values[i], b.values[i]);
        hi.values[i] = Math.max(a.values[i], b.values[i]);
    }
    return AlignedBox.fromMinMax(lo, hi);
});

// A box together with one of its points.
const v21ShapePoint = fc.tuple(v21Shape,
    fc.array(fc.double({ min: 0, max: 1, noNaN: true }),
        { minLength: 3, maxLength: 3 }))
    .map(([b, u]) => {
        const q = new Vector(3);
        for (let i = 0; i < 3; ++i) {
            q.values[i] = (1 - u[i]) * b.min.values[i] + u[i] * b.max.values[i];
        }
        return [b, q] as [AlignedBox, Vector];
    });

function v21PointDistance(p: Vector, b: AlignedBox): number {
    return new DistPointAlignedBox().compute(p, b).distance;
}

function v21CheckShapePoint(b: AlignedBox,
    res: { closest: [Vector, Vector] }): void {
    for (let i = 0; i < 3; ++i) {
        expect(res.closest[1].values[i])
            .toBeGreaterThanOrEqual(b.min.values[i] - 1e-9);
        expect(res.closest[1].values[i])
            .toBeLessThanOrEqual(b.max.values[i] + 1e-9);
    }
}

// An aligned box is not closed under rotation, so the equivariance property
// uses translations only (rot is ignored).
function v21MoveShape(b: AlignedBox, rot: (x: Vector) => Vector,
    tr: Vector): AlignedBox {
    void rot;
    return AlignedBox.fromMinMax(add(b.min, tr), add(b.max, tr));
}

function v21ShapeSnapshot(b: AlignedBox): number[] {
    return [...b.min.values, ...b.max.values];
}

const v21Segment = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -8, 8))
    .filter(([a, b]) => length(sub(b, a)) > 1e-2)
    .map(([a, b]) => Segment.fromEndpoints(a, b));

// Minimum of a convex function on [lo,hi] by ternary search. The distance
// from a point to a convex set is convex and the segment is an affine image of
// its parameter, so the composition is convex and the search is exact.
function v21MinOnInterval(f: (t: number) => number, lo: number,
    hi: number): number {
    let a = lo;
    let b = hi;
    for (let i = 0; i < 140; ++i) {
        const m1 = a + (b - a) / 3;
        const m2 = b - (b - a) / 3;
        if (f(m1) <= f(m2)) {
            b = m2;
        }
        else {
            a = m1;
        }
    }
    return Math.min(f(a), Math.min(f(b), f(0.5 * (a + b))));
}

describe('DistSegment3AlignedBox3 verification', () => {
    const query = new DistSegment3AlignedBox3();

    it('result is self consistent and the points lie on their primitives',
        () => {
            check(fc.tuple(v21Segment, v21Shape), ([seg, s]) => {
                const res = query.compute(seg, s);
                expectClose(res.distance, Math.sqrt(res.sqrDistance), 1e-12,
                    1e-12);
                const diff = sub(res.closest[0], res.closest[1]);
                expectClose(res.sqrDistance, dot(diff, diff), 1e-9, 1e-9);
                // Upstream parameterizes the segment as P0 + t*(P1-P0) with
                // 0 <= t <= 1 (not the centered form), so the parameter must
                // be in [0,1].
                expect(res.parameter).toBeGreaterThanOrEqual(0);
                expect(res.parameter).toBeLessThanOrEqual(1);
                expectVectorClose(res.closest[0],
                    add(seg.p[0],
                        mul(res.parameter, sub(seg.p[1], seg.p[0]))), 1e-9,
                    1e-9);
                v21CheckShapePoint(s, res);
            });
        });

    it('matches a convex minimization along the segment', () => {
        check(fc.tuple(v21Segment, v21Shape), ([seg, s]) => {
            const dir = sub(seg.p[1], seg.p[0]);
            const f = (t: number): number =>
                v21PointDistance(add(seg.p[0], mul(t, dir)), s);
            expectClose(query.compute(seg, s).distance,
                v21MinOnInterval(f, 0, 1), 1e-7, 1e-7);
        }, 60);
    }, 30000);

    it('agrees with the line query inside [0,1] and the point query outside',
        () => {
            check(fc.tuple(v21Segment, v21Shape), ([seg, s]) => {
                const line = Line.fromOriginDirection(seg.p[0],
                    sub(seg.p[1], seg.p[0]));
                const lr = new DistLine3AlignedBox3().compute(line, s);
                const sr = query.compute(seg, s);
                if (lr.parameter >= 0 && lr.parameter <= 1) {
                    expect(sr.parameter).toBe(lr.parameter);
                    expect(sr.distance).toBe(lr.distance);
                    expectVectorClose(sr.closest[0], lr.closest[0], 0, 0);
                    expectVectorClose(sr.closest[1], lr.closest[1], 0, 0);
                }
                else {
                    const end = lr.parameter < 0 ? 0 : 1;
                    expect(sr.parameter).toBe(end);
                    expectVectorClose(sr.closest[0], seg.p[end], 0, 0);
                    expect(sr.distance)
                        .toBe(v21PointDistance(seg.p[end], s));
                }
            });
        });

    it('reports zero distance when an endpoint is on the shape', () => {
        check(fc.tuple(v21ShapePoint, wellScaledVector(3, -6, 6)),
            ([[s, q], other]) => {
                const seg = Segment.fromEndpoints(q, add(q, other));
                expect(query.compute(seg, s).distance)
                    .toBeLessThanOrEqual(1e-9);
            });
    });

    it('is invariant under a common translation', () => {
        // An aligned box is not closed under rotation, so only translations
        // move the whole configuration rigidly.
        check(fc.tuple(v21Segment, v21Shape, rotationFrame(3),
            wellScaledVector(3, -5, 5)), ([seg, s, R, tr]) => {
            void R;
            const rot = (x: Vector): Vector => x.clone();
            const moved = Segment.fromEndpoints(add(rot(seg.p[0]), tr),
                add(rot(seg.p[1]), tr));
            const r0 = query.compute(seg, s);
            const r1 = query.compute(moved, v21MoveShape(s, rot, tr));
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
        });
    });

    it('a short segment agrees with the point query at its first endpoint',
        () => {
            check(fc.tuple(wellScaledVector(3, -8, 8), unitVector(3),
                v21Shape), ([p0, d, s]) => {
                const seg = Segment.fromEndpoints(p0, add(p0, mul(1e-9, d)));
                expectClose(query.compute(seg, s).distance,
                    v21PointDistance(p0, s), 1e-7, 1e-7);
            });
        });

    it('does not mutate its inputs', () => {
        check(fc.tuple(v21Segment, v21Shape), ([seg, s]) => {
            const p0 = seg.p[0].clone();
            const p1 = seg.p[1].clone();
            const snapshot = v21ShapeSnapshot(s);
            const res = query.compute(seg, s);
            expect(seg.p[0].values).toEqual(p0.values);
            expect(seg.p[1].values).toEqual(p1.values);
            expect(v21ShapeSnapshot(s)).toEqual(snapshot);
            res.closest[0].values[0] = 4242;
            res.closest[1].values[0] = 4242;
            expect(seg.p[0].values).toEqual(p0.values);
            expect(v21ShapeSnapshot(s)).toEqual(snapshot);
        });
    });
});
