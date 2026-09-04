import { describe, expect, it } from 'vitest';
import { DistSegment2OrientedBox2 } from '../src/DistSegment2OrientedBox2.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, sub }
    from '../src/Vector.js';
import { DistLine2OrientedBox2 } from '../src/DistLine2OrientedBox2.js';
import { DistPointOrientedBox } from '../src/DistPointOrientedBox.js';
import { Line } from '../src/Line.js';
import {
    check, expectClose, expectVectorClose, fc, positive,
    rotationFrame, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(center: number[], axis: number[][],
    extent: number[]): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(v(...center),
        axis.map(a => v(...a)), v(...extent));
}

// The exact squared distance from a point to the solid oriented box, computed
// independently of the library.
function pointBoxSqrDistance(p: Vector, b: OrientedBox): number {
    const delta = sub(p, b.center);
    let sqrDistance = 0;
    for (let i = 0; i < p.size; ++i) {
        const y = dot(delta, b.axis[i]);
        const excess = Math.max(0, Math.abs(y) - b.extent.values[i]);
        sqrDistance += excess * excess;
    }
    return sqrDistance;
}

// Verify that the reported box point is inside the box.
function verifyBoxPoint(b: OrientedBox, q: Vector): void {
    const delta = sub(q, b.center);
    for (let i = 0; i < q.size; ++i) {
        expect(Math.abs(dot(delta, b.axis[i])))
            .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
    }
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

// The distance from a segment point to the solid box is a convex function of
// the segment parameter, so a ternary search over [0,1] converges to the true
// minimum independently of the query under test.
function bruteForce(s: Segment, b: OrientedBox): number {
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
function verifyClosest(s: Segment, b: OrientedBox,
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

describe('DistSegment2OrientedBox2', () => {
    const query = new DistSegment2OrientedBox2();
    const axisAligned = box([0, 0], [[1, 0], [0, 1]], [1, 1]);
    const c = Math.SQRT1_2;
    const rotated = box([0, 0], [[c, c], [-c, c]], [1, 1]);

    it('clamps to the first endpoint', () => {
        const s = segment([3, 0], [5, 0]);
        const result = query.compute(s, axisAligned);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.parameter).toBe(0);
    });

    it('clamps to the second endpoint', () => {
        const s = segment([5, 0], [3, 0]);
        const result = query.compute(s, axisAligned);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.parameter).toBe(1);
    });

    it('handles a box rotated by 45 degrees', () => {
        const s = segment([4, 0], [3, 0]);
        const result = query.compute(s, rotated);
        expect(result.distance).toBeCloseTo(3 - Math.SQRT2, 10);
        verifyClosest(s, rotated, result);
    });

    it('handles a degenerate zero-length segment', () => {
        const s = segment([0, 4], [0, 4]);
        const result = query.compute(s, rotated);
        expect(result.distance).toBeCloseTo(4 - Math.SQRT2, 10);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 66554433;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 200; ++trial) {
            const angle = 2 * Math.PI * rand();
            const ca = Math.cos(angle);
            const sa = Math.sin(angle);
            const b = box([2 * rand() - 1, 2 * rand() - 1],
                [[ca, sa], [-sa, ca]],
                [0.2 + rand(), 0.2 + rand()]);
            const s = segment([8 * rand() - 4, 8 * rand() - 4],
                [8 * rand() - 4, 8 * rand() - 4]);
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
// header DistSegment2OrientedBox2.h. Upstream solves the query on the line containing the
// segment and clamps the result to the segment domain, so the properties below
// check the line/point agreement as well as the geometric invariants.
// ---------------------------------------------------------------------------

const v21Shape = fc.tuple(wellScaledVector(2, -8, 8), rotationFrame(2),
    fc.array(positive(4, 0.1), { minLength: 2, maxLength: 2 }))
    .map(([c, axes, e]) => OrientedBox.fromCenterAxisExtent(c, axes,
        Vector.fromArray(e)));

// A box together with one of its points.
const v21ShapePoint = fc.tuple(v21Shape,
    fc.array(fc.double({ min: -1, max: 1, noNaN: true }),
        { minLength: 2, maxLength: 2 }))
    .map(([b, u]) => {
        let q = b.center.clone();
        for (let i = 0; i < 2; ++i) {
            q = add(q, mul(u[i] * b.extent.values[i], b.axis[i]));
        }
        return [b, q] as [OrientedBox, Vector];
    });

function v21PointDistance(p: Vector, b: OrientedBox): number {
    return new DistPointOrientedBox().compute(p, b).distance;
}

function v21CheckShapePoint(b: OrientedBox,
    res: { closest: [Vector, Vector] }): void {
    const d = sub(res.closest[1], b.center);
    for (let i = 0; i < 2; ++i) {
        expect(Math.abs(dot(d, b.axis[i])))
            .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
    }
}

function v21MoveShape(b: OrientedBox, rot: (x: Vector) => Vector,
    tr: Vector): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(add(rot(b.center), tr),
        b.axis.map(a => rot(a)), b.extent);
}

function v21ShapeSnapshot(b: OrientedBox): number[] {
    return [...b.center.values, ...b.extent.values,
        ...b.axis.flatMap(a => [...a.values])];
}

const v21Segment = fc.tuple(wellScaledVector(2, -8, 8),
    wellScaledVector(2, -8, 8))
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

describe('DistSegment2OrientedBox2 verification', () => {
    const query = new DistSegment2OrientedBox2();

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
                const lr = new DistLine2OrientedBox2().compute(line, s);
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
        check(fc.tuple(v21ShapePoint, wellScaledVector(2, -6, 6)),
            ([[s, q], other]) => {
                const seg = Segment.fromEndpoints(q, add(q, other));
                expect(query.compute(seg, s).distance)
                    .toBeLessThanOrEqual(1e-9);
            });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(v21Segment, v21Shape, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([seg, s, R, tr]) => {
            const rot = (x: Vector): Vector => {
                let y = new Vector(2);
                for (let i = 0; i < 2; ++i) {
                    y = add(y, mul(x.values[i], R[i]));
                }
                return y;
            };
            const moved = Segment.fromEndpoints(add(rot(seg.p[0]), tr),
                add(rot(seg.p[1]), tr));
            const r0 = query.compute(seg, s);
            const r1 = query.compute(moved, v21MoveShape(s, rot, tr));
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
        });
    });

    it('a short segment agrees with the point query at its first endpoint',
        () => {
            check(fc.tuple(wellScaledVector(2, -8, 8), unitVector(2),
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
