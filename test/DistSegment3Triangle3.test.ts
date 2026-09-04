import { describe, expect, it } from 'vitest';
import { DistPointTriangle } from '../src/DistPointTriangle.js';
import { DistSegment3Triangle3 } from '../src/DistSegment3Triangle3.js';
import type { DistSegment3Triangle3Result }
    from '../src/DistSegment3Triangle3.js';
import { Segment } from '../src/Segment.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, length, mul, sub }
    from '../src/Vector.js';
import { DistLine3Triangle3 } from '../src/DistLine3Triangle3.js';
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

function triangle(a: number[], b: number[], c: number[]): Triangle {
    return Triangle.fromVertices(v(...a), v(...b), v(...c));
}

// The exact squared distance from a point to a solid triangle, computed
// independently of the ported query: minimize over the interior candidate
// (the projection onto the plane, if its barycentric coordinates are
// nonnegative) and over the three edges.
function pointTriangleSqrDistance(p: Vector, t: Triangle): number {
    const e1 = sub(t.v[1], t.v[0]);
    const e2 = sub(t.v[2], t.v[0]);
    const d = sub(p, t.v[0]);
    const a = dot(e1, e1), b = dot(e1, e2), c = dot(e2, e2);
    const d1 = dot(e1, d), d2 = dot(e2, d);
    const det = a * c - b * b;

    const sqrLenTo = (q: Vector) => {
        const diff = sub(p, q);
        return dot(diff, diff);
    };

    let best = Number.MAX_VALUE;
    if (det > 0) {
        const s = (c * d1 - b * d2) / det;
        const u = (a * d2 - b * d1) / det;
        if (s >= 0 && u >= 0 && s + u <= 1) {
            best = sqrLenTo(add(t.v[0], add(mul(s, e1), mul(u, e2))));
        }
    }

    // Closest point on each edge segment.
    for (const [i, j] of [[0, 1], [1, 2], [2, 0]]) {
        const dir = sub(t.v[j], t.v[i]);
        const len2 = dot(dir, dir);
        let r = 0;
        if (len2 > 0) {
            r = Math.min(Math.max(dot(sub(p, t.v[i]), dir) / len2, 0), 1);
        }
        const value = sqrLenTo(add(t.v[i], mul(r, dir)));
        if (value < best) {
            best = value;
        }
    }
    return best;
}

// The squared distance from segment(t) to the solid triangle is a convex
// function of t on [0,1], so a ternary search finds its minimum. A dense
// sampling is used as a second opinion.
function bruteForceSqrDistance(s: Segment, t: Triangle): number {
    const f = (u: number) => pointTriangleSqrDistance(
        add(s.p[0], mul(u, sub(s.p[1], s.p[0]))), t);

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
// their primitives, the barycentric coordinates describe closest[1], and the
// pair realizes the reported distance.
function expectConsistent(result: DistSegment3Triangle3Result,
    s: Segment, t: Triangle): void {
    expect(result.parameter).toBeGreaterThanOrEqual(-1e-12);
    expect(result.parameter).toBeLessThanOrEqual(1 + 1e-12);

    const onSegment = add(s.p[0],
        mul(result.parameter, sub(s.p[1], s.p[0])));
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onSegment.values[i], 6);
    }

    let sum = 0;
    let onTriangle = new Vector(3);
    for (let i = 0; i < 3; ++i) {
        expect(result.barycentric[i]).toBeGreaterThanOrEqual(-1e-9);
        sum += result.barycentric[i];
        onTriangle = add(onTriangle, mul(result.barycentric[i], t.v[i]));
    }
    expect(sum).toBeCloseTo(1, 9);
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[1].values[i])
            .toBeCloseTo(onTriangle.values[i], 6);
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

describe('DistSegment3Triangle3', () => {
    const query = new DistSegment3Triangle3();
    // The right triangle with vertices (0,0,0), (1,0,0), (0,1,0).
    const tri = triangle([0, 0, 0], [1, 0, 0], [0, 1, 0]);

    it('measures a segment parallel to the triangle plane', () => {
        const s = segment([0.25, 0.25, 3], [0.3, 0.3, 3]);
        const result = query.compute(s, tri);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(0, 12);
        expectConsistent(result, s, tri);
    });

    it('reports zero distance when the segment crosses the triangle', () => {
        const s = segment([0.25, 0.25, -1], [0.25, 0.25, 1]);
        const result = query.compute(s, tri);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter).toBeCloseTo(0.5, 10);
        expect(result.barycentric[0]).toBeCloseTo(0.5, 10);
        expect(result.barycentric[1]).toBeCloseTo(0.25, 10);
        expect(result.barycentric[2]).toBeCloseTo(0.25, 10);
        expectConsistent(result, s, tri);
    });

    it('reports zero distance when an endpoint touches the triangle', () => {
        const s = segment([0.25, 0.25, 0], [0.25, 0.25, 5]);
        const result = query.compute(s, tri);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, tri);
    });

    it('misses the triangle when the crossing point is outside it', () => {
        // The line crosses z = 0 at (3,0,0), outside the triangle; the
        // closest triangle point is the vertex (1,0,0).
        const s = segment([3, 0, -1], [3, 0, 1]);
        const result = query.compute(s, tri);
        expect(result.distance).toBeCloseTo(2, 10);
        expect(result.barycentric[1]).toBeCloseTo(1, 10);
        expectConsistent(result, s, tri);
    });

    it('measures a segment parallel to an edge', () => {
        // The segment lies at y = -2 in the plane z = 0, spanning x in
        // [0,1]; the closest triangle points are on the edge <V0,V1>.
        const s = segment([0, -2, 0], [1, -2, 0]);
        const result = query.compute(s, tri);
        expect(result.distance).toBeCloseTo(2, 10);
        expect(result.barycentric[2]).toBeCloseTo(0, 10);
        expectConsistent(result, s, tri);
    });

    it('clamps to the first endpoint when the line minimum is behind it',
        () => {
            const s = segment([4, 0, 0], [7, 0, 0]);
            const result = query.compute(s, tri);
            expect(result.parameter).toBe(0);
            expect(result.distance).toBeCloseTo(3, 10);
            expect(result.barycentric[1]).toBeCloseTo(1, 10);
            expectConsistent(result, s, tri);
        });

    it('clamps to the second endpoint when the line minimum is beyond it',
        () => {
            const s = segment([7, 0, 0], [4, 0, 0]);
            const result = query.compute(s, tri);
            expect(result.parameter).toBe(1);
            expect(result.distance).toBeCloseTo(3, 10);
            expectConsistent(result, s, tri);
        });

    it('gives the same distance for both segment orientations', () => {
        const t = triangle([1, 0, 0], [0, 2, 1], [-1, -1, 2]);
        const cases: Array<[number[], number[]]> = [
            [[4, 0, 0], [7, 0, 0]],
            [[0.25, 0.25, 3], [0.3, 0.3, 3]],
            [[-3, -3, -3], [3, 3, 3]],
            [[2, 2, -5], [2, 2, 5]],
            [[0.5, 0.5, 0.5], [7, 1, -3]]
        ];
        for (const [p0, p1] of cases) {
            const forward = query.compute(segment(p0, p1), t);
            const backward = query.compute(segment(p1, p0), t);
            expect(forward.distance).toBeCloseTo(backward.distance, 9);
        }
    });

    it('handles a degenerate zero-length segment off the triangle', () => {
        const s = segment([3, 4, 5], [3, 4, 5]);
        const result = query.compute(s, tri);
        const expected = new DistPointTriangle().compute(v(3, 4, 5), tri);
        expect(result.distance).toBeCloseTo(expected.distance, 10);
        expect(result.closest[0].values[0]).toBeCloseTo(3, 12);
        expectConsistent(result, s, tri);
    });

    it('handles a degenerate zero-length segment on the triangle', () => {
        const s = segment([0.25, 0.25, 0], [0.25, 0.25, 0]);
        const result = query.compute(s, tri);
        expect(result.distance).toBeCloseTo(0, 10);
        expectConsistent(result, s, tri);
    });

    it('agrees with a brute-force sampling for random configurations', () => {
        const random = makeRandom(87082);
        for (let trial = 0; trial < 250; ++trial) {
            const t = triangle(
                [3 * random() - 1.5, 3 * random() - 1.5, 3 * random() - 1.5],
                [3 * random() - 1.5, 3 * random() - 1.5, 3 * random() - 1.5],
                [3 * random() - 1.5, 3 * random() - 1.5, 3 * random() - 1.5]);
            const s = segment(
                [6 * random() - 3, 6 * random() - 3, 6 * random() - 3],
                [6 * random() - 3, 6 * random() - 3, 6 * random() - 3]);
            const result = query.compute(s, t);
            const expected = Math.sqrt(bruteForceSqrDistance(s, t));
            expect(Math.abs(result.distance - expected)).toBeLessThan(1e-6);
            expectConsistent(result, s, t);
        }
    }, 30000);
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistSegment3Triangle3.h. Upstream solves the query on the line containing the
// segment and clamps the result to the segment domain, so the properties below
// check the line/point agreement as well as the geometric invariants.
// ---------------------------------------------------------------------------

const v21Shape = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -8, 8), wellScaledVector(3, -8, 8))
    .filter(([a, b, c]) => {
        const e0 = sub(b, a);
        const e1 = sub(c, a);
        return length(cross(e0, e1)) > 1e-2;
    })
    .map(([a, b, c]) => Triangle.fromVertices(a, b, c));

// A triangle together with one of its points.
const v21ShapePoint = fc.tuple(v21Shape,
    fc.double({ min: 0, max: 1, noNaN: true }),
    fc.double({ min: 0, max: 1, noNaN: true }))
    .map(([tri, u, w]) => {
        // Map the unit square onto the triangle by folding.
        let b1 = u;
        let b2 = w;
        if (b1 + b2 > 1) {
            b1 = 1 - b1;
            b2 = 1 - b2;
        }
        const q = add(tri.v[0], add(mul(b1, sub(tri.v[1], tri.v[0])),
            mul(b2, sub(tri.v[2], tri.v[0]))));
        return [tri, q] as [Triangle, Vector];
    });

function v21PointDistance(p: Vector, tri: Triangle): number {
    return new DistPointTriangle().compute(p, tri).distance;
}

function v21CheckShapePoint(tri: Triangle,
    res: { closest: [Vector, Vector],
        barycentric: [number, number, number] }): void {
    let sum = 0;
    for (let i = 0; i < 3; ++i) {
        expect(res.barycentric[i]).toBeGreaterThanOrEqual(-1e-9);
        expect(res.barycentric[i]).toBeLessThanOrEqual(1 + 1e-9);
        sum += res.barycentric[i];
    }
    expectClose(sum, 1, 1e-9, 1e-9);
    let x = new Vector(3);
    for (let i = 0; i < 3; ++i) {
        x = add(x, mul(res.barycentric[i], tri.v[i]));
    }
    expectVectorClose(res.closest[1], x, 1e-8, 1e-8);
}

function v21MoveShape(tri: Triangle, rot: (x: Vector) => Vector,
    tr: Vector): Triangle {
    return Triangle.fromVertexArray(tri.v.map(x => add(rot(x), tr)));
}

function v21ShapeSnapshot(tri: Triangle): number[] {
    return tri.v.flatMap(x => [...x.values]);
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

describe('DistSegment3Triangle3 verification', () => {
    const query = new DistSegment3Triangle3();

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
                const lr = new DistLine3Triangle3().compute(line, s);
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
