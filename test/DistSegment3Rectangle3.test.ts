import { describe, expect, it } from 'vitest';
import { DistPointRectangle } from '../src/DistPointRectangle.js';
import { DistSegment3Rectangle3 } from '../src/DistSegment3Rectangle3.js';
import type { DistSegment3Rectangle3Result }
    from '../src/DistSegment3Rectangle3.js';
import { Rectangle } from '../src/Rectangle.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, sub }
    from '../src/Vector.js';
import { DistLine3Rectangle3 } from '../src/DistLine3Rectangle3.js';
import { Line } from '../src/Line.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, positive,
    rotationFrame, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

// The rectangle in the z = 0 plane with the standard x and y axes.
function xyRectangle(e0: number, e1: number, center: number[] = [0, 0, 0]):
    Rectangle {
    return Rectangle.fromCenterAxisExtent(v(...center),
        [v(1, 0, 0), v(0, 1, 0)], v(e0, e1));
}

// An orthonormal frame; the first two vectors are the rectangle axes.
function frame(a: number, b: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    return [
        v(ca, sa, 0),
        v(-sa * cb, ca * cb, sb),
        v(sa * sb, -ca * sb, cb)
    ];
}

// The exact squared distance from a point to a solid rectangle: clamp the
// rectangle coordinates of the point to the extents.
function pointRectangleSqrDistance(p: Vector, r: Rectangle): number {
    const delta = sub(p, r.center);
    let closest = r.center.clone();
    for (let i = 0; i < 2; ++i) {
        const e = r.extent.values[i];
        const s = Math.min(Math.max(dot(r.axis[i], delta), -e), e);
        closest = add(closest, mul(s, r.axis[i]));
    }
    const d = sub(p, closest);
    return dot(d, d);
}

// The squared distance from segment(t) to the solid rectangle is a convex
// function of t on [0,1], so a ternary search finds its minimum. A dense
// sampling is used as a second opinion.
function bruteForceSqrDistance(s: Segment, r: Rectangle): number {
    const f = (u: number) => pointRectangleSqrDistance(
        add(s.p[0], mul(u, sub(s.p[1], s.p[0]))), r);

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
    let best = f(0.5 * (lo + hi));
    for (let i = 0; i <= 2000; ++i) {
        const value = f(i / 2000);
        if (value < best) {
            best = value;
        }
    }
    return best;
}

// Verify the internal consistency of a result: the closest points lie on
// their primitives, the reported W-coordinates describe closest[1], and the
// pair realizes the reported distance.
function expectConsistent(result: DistSegment3Rectangle3Result,
    s: Segment, r: Rectangle): void {
    expect(result.parameter).toBeGreaterThanOrEqual(-1e-12);
    expect(result.parameter).toBeLessThanOrEqual(1 + 1e-12);

    const onSegment = add(s.p[0],
        mul(result.parameter, sub(s.p[1], s.p[0])));
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onSegment.values[i], 6);
    }

    let onRectangle = r.center.clone();
    for (let i = 0; i < 2; ++i) {
        expect(Math.abs(result.cartesian[i]))
            .toBeLessThanOrEqual(r.extent.values[i] + 1e-9);
        onRectangle = add(onRectangle,
            mul(result.cartesian[i], r.axis[i]));
    }
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[1].values[i])
            .toBeCloseTo(onRectangle.values[i], 6);
    }

    const diff = sub(result.closest[0], result.closest[1]);
    expect(Math.sqrt(dot(diff, diff))).toBeCloseTo(result.distance, 6);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 8);
}

// A small deterministic linear congruential generator.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistSegment3Rectangle3', () => {
    const query = new DistSegment3Rectangle3();
    const rect = xyRectangle(2, 1);

    it('measures a segment parallel to the rectangle plane', () => {
        const s = segment([-1, 0, 3], [1, 0, 3]);
        const result = query.compute(s, rect);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(0, 12);
        expectConsistent(result, s, rect);
    });

    it('reports zero distance when the segment crosses the rectangle', () => {
        const s = segment([0.5, 0.25, -2], [0.5, 0.25, 2]);
        const result = query.compute(s, rect);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter).toBeCloseTo(0.5, 10);
        expect(result.cartesian[0]).toBeCloseTo(0.5, 10);
        expect(result.cartesian[1]).toBeCloseTo(0.25, 10);
        expectConsistent(result, s, rect);
    });

    it('reports zero distance when an endpoint lies on the rectangle', () => {
        const s = segment([1, 0.5, 0], [1, 0.5, 4]);
        const result = query.compute(s, rect);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, rect);
    });

    it('misses the rectangle when the crossing point is outside it', () => {
        // The line crosses z = 0 at (5,0,0), well beyond the x extent 2, so
        // the closest rectangle point is on the edge x = 2.
        const s = segment([5, 0, -1], [5, 0, 1]);
        const result = query.compute(s, rect);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.cartesian[0]).toBeCloseTo(2, 10);
        expect(result.cartesian[1]).toBeCloseTo(0, 10);
        expectConsistent(result, s, rect);
    });

    it('clamps to the first endpoint when the line minimum is behind it',
        () => {
            const s = segment([6, 0, 0], [9, 0, 0]);
            const result = query.compute(s, rect);
            expect(result.parameter).toBe(0);
            expect(result.distance).toBeCloseTo(4, 10);
            expect(result.cartesian[0]).toBeCloseTo(2, 10);
            expectConsistent(result, s, rect);
        });

    it('clamps to the second endpoint when the line minimum is beyond it',
        () => {
            const s = segment([9, 0, 0], [6, 0, 0]);
            const result = query.compute(s, rect);
            expect(result.parameter).toBe(1);
            expect(result.distance).toBeCloseTo(4, 10);
            expectConsistent(result, s, rect);
        });

    it('measures the distance to a rectangle corner', () => {
        // The nearest rectangle point to the segment is the corner (2,1,0).
        const s = segment([5, 4, 0], [7, 6, 0]);
        const result = query.compute(s, rect);
        expect(result.parameter).toBe(0);
        expect(result.distance).toBeCloseTo(Math.hypot(3, 3), 10);
        expect(result.cartesian[0]).toBeCloseTo(2, 10);
        expect(result.cartesian[1]).toBeCloseTo(1, 10);
        expectConsistent(result, s, rect);
    });

    it('gives the same distance for both segment orientations', () => {
        const axes = frame(0.6, 1.3);
        const r = Rectangle.fromCenterAxisExtent(v(1, -1, 2),
            [axes[0], axes[1]], v(1.5, 0.75));
        const cases: Array<[number[], number[]]> = [
            [[5, 0, 0], [8, 0, 0]],
            [[-1, 0, 3], [1, 0, 3]],
            [[-3, -3, -3], [3, 3, 3]],
            [[2, 2, -5], [2, 2, 5]],
            [[0.5, 0.5, 0.5], [7, 1, -3]]
        ];
        for (const [p0, p1] of cases) {
            const forward = query.compute(segment(p0, p1), r);
            const backward = query.compute(segment(p1, p0), r);
            expect(forward.distance).toBeCloseTo(backward.distance, 9);
        }
    });

    it('handles a degenerate zero-length segment off the rectangle', () => {
        const s = segment([3, 4, 5], [3, 4, 5]);
        const result = query.compute(s, rect);
        const expected = new DistPointRectangle().compute(v(3, 4, 5), rect);
        expect(result.distance).toBeCloseTo(expected.distance, 10);
        expect(result.closest[0].values[0]).toBeCloseTo(3, 12);
        expectConsistent(result, s, rect);
    });

    it('handles a degenerate zero-length segment on the rectangle', () => {
        const s = segment([0.5, -0.25, 0], [0.5, -0.25, 0]);
        const result = query.compute(s, rect);
        expect(result.distance).toBeCloseTo(0, 10);
        expectConsistent(result, s, rect);
    });

    it('agrees with a brute-force sampling for random configurations', () => {
        const random = makeRandom(87084);
        for (let trial = 0; trial < 250; ++trial) {
            const axes = frame(2 * Math.PI * random(), Math.PI * random());
            const r = Rectangle.fromCenterAxisExtent(
                v(3 * random() - 1.5, 3 * random() - 1.5, 3 * random() - 1.5),
                [axes[0], axes[1]],
                v(0.2 + 2 * random(), 0.2 + 2 * random()));
            const s = segment(
                [6 * random() - 3, 6 * random() - 3, 6 * random() - 3],
                [6 * random() - 3, 6 * random() - 3, 6 * random() - 3]);
            const result = query.compute(s, r);
            const expected = Math.sqrt(bruteForceSqrDistance(s, r));
            expect(Math.abs(result.distance - expected)).toBeLessThan(1e-6);
            expectConsistent(result, s, r);
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistSegment3Rectangle3.h. Upstream solves the query on the line containing the
// segment and clamps the result to the segment domain, so the properties below
// check the line/point agreement as well as the geometric invariants.
// ---------------------------------------------------------------------------

const v21Shape = fc.tuple(wellScaledVector(3, -8, 8), rotationFrame(3),
    positive(4, 0.1), positive(4, 0.1))
    .map(([c, axes, e0, e1]) => Rectangle.fromCenterAxisExtent(c,
        [axes[0], axes[1]], Vector.fromArray([e0, e1])));

// A rectangle together with one of its points.
const v21ShapePoint = fc.tuple(v21Shape,
    fc.double({ min: -1, max: 1, noNaN: true }),
    fc.double({ min: -1, max: 1, noNaN: true }))
    .map(([rect, u0, u1]) => {
        const q = add(rect.center,
            add(mul(u0 * rect.extent.values[0], rect.axis[0]),
                mul(u1 * rect.extent.values[1], rect.axis[1])));
        return [rect, q] as [Rectangle, Vector];
    });

function v21PointDistance(p: Vector, rect: Rectangle): number {
    return new DistPointRectangle().compute(p, rect).distance;
}

function v21CheckShapePoint(rect: Rectangle,
    res: { closest: [Vector, Vector], cartesian: [number, number] }): void {
    const d = sub(res.closest[1], rect.center);
    for (let i = 0; i < 2; ++i) {
        expectClose(dot(d, rect.axis[i]), res.cartesian[i], 1e-8, 1e-8);
        expect(Math.abs(res.cartesian[i]))
            .toBeLessThanOrEqual(rect.extent.values[i] + 1e-9);
    }
    // The closest point is in the plane of the rectangle.
    const normal = cross(rect.axis[0], rect.axis[1]);
    expectClose(dot(d, normal), 0, 1e-8, 1e-8);
    expectVectorClose(res.closest[1], add(rect.center,
        add(mul(res.cartesian[0], rect.axis[0]),
            mul(res.cartesian[1], rect.axis[1]))), 1e-8, 1e-8);
}

function v21MoveShape(rect: Rectangle, rot: (x: Vector) => Vector,
    tr: Vector): Rectangle {
    return Rectangle.fromCenterAxisExtent(add(rot(rect.center), tr),
        rect.axis.map(a => rot(a)), rect.extent);
}

function v21ShapeSnapshot(rect: Rectangle): number[] {
    return [...rect.center.values, ...rect.extent.values,
        ...rect.axis.flatMap(a => [...a.values])];
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

describe('DistSegment3Rectangle3 verification', () => {
    const query = new DistSegment3Rectangle3();

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
            const res = query.compute(seg, s);
            const dir = sub(seg.p[1], seg.p[0]);
            const f = (t: number): number =>
                v21PointDistance(add(seg.p[0], mul(t, dir)), s);
            expectClose(res.distance, v21MinOnInterval(f, 0, 1), 1e-7, 1e-7);
        }, 60);
    }, 30000);

    it('agrees with the line query inside [0,1] and the point query outside',
        () => {
            check(fc.tuple(v21Segment, v21Shape), ([seg, s]) => {
                const line = Line.fromOriginDirection(seg.p[0],
                    sub(seg.p[1], seg.p[0]));
                const lr = new DistLine3Rectangle3().compute(line, s);
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
                const zres = query.compute(seg, s);
                // The segment contains q, so its distance to the shape is at
                // most q's own distance to the shape. That distance is zero
                // in exact arithmetic; the sampled shape point carries
                // rounding error that grows with the coordinate magnitudes
                // and, for a thin triangle, with 1/area, so the near-zero
                // bound is scale relative.
                const pd = v21PointDistance(q, s);
                expect(zres.distance).toBeLessThanOrEqual(pd + 1e-9);
                expect(pd).toBeLessThanOrEqual(1e-7 * (1 + length(q)));
            });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(v21Segment, v21Shape, rotationFrame(3),
            wellScaledVector(3, -5, 5)), ([seg, s, R, tr]) => {
            const rot = (x: Vector): Vector => {
                let y = new Vector(3);
                for (let i = 0; i < 3; ++i) {
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
            check(fc.tuple(wellScaledVector(3, -8, 8), unitVector(3),
                v21Shape), ([p0, d, s]) => {
                const seg = Segment.fromEndpoints(p0, add(p0, mul(1e-9, d)));
                const res = query.compute(seg, s);
                expectClose(res.distance, v21PointDistance(p0, s), 1e-7,
                    1e-7);
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
