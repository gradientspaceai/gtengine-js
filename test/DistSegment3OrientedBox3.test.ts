import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistPointOrientedBox } from '../src/DistPointOrientedBox.js';
import { DistSegment3CanonicalBox3 } from '../src/DistSegment3CanonicalBox3.js';
import { DistSegment3OrientedBox3 } from '../src/DistSegment3OrientedBox3.js';
import type { DistSegment3OrientedBox3Result }
    from '../src/DistSegment3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, sub }
    from '../src/Vector.js';
import { DistLine3OrientedBox3 } from '../src/DistLine3OrientedBox3.js';
import { Line } from '../src/Line.js';
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

// An axis-aligned oriented box (identity axes) with the given center and
// extents.
function alignedObb(center: number[], extent: number[]): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(v(...center),
        [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], v(...extent));
}

// A rotation about the z axis by angle a, then about the x axis by angle b,
// applied to the standard basis; the result is an orthonormal frame.
function frame(a: number, b: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    // Rz(a) columns, then Rx(b) applied to them.
    const u0 = v(ca, sa, 0);
    const u1 = v(-sa * cb, ca * cb, sb);
    const u2 = v(sa * sb, -ca * sb, cb);
    return [u0, u1, u2];
}

// The exact squared distance from a point to a solid oriented box: clamp the
// box coordinates of the point to the extents.
function pointBoxSqrDistance(p: Vector, b: OrientedBox): number {
    const delta = sub(p, b.center);
    let closest = b.center.clone();
    for (let i = 0; i < 3; ++i) {
        const e = b.extent.values[i];
        const y = Math.min(Math.max(dot(b.axis[i], delta), -e), e);
        closest = add(closest, mul(y, b.axis[i]));
    }
    const d = sub(p, closest);
    return dot(d, d);
}

// The squared distance from segment(t) to the solid box is a convex function
// of t on [0,1], so a ternary search finds its minimum. A dense sampling is
// used as a second opinion.
function bruteForceSqrDistance(s: Segment, b: OrientedBox): number {
    const f = (u: number) => pointBoxSqrDistance(
        add(s.p[0], mul(u, sub(s.p[1], s.p[0]))), b);

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
// their primitives and realize the reported distance.
function expectConsistent(result: DistSegment3OrientedBox3Result,
    s: Segment, b: OrientedBox): void {
    expect(result.parameter).toBeGreaterThanOrEqual(-1e-12);
    expect(result.parameter).toBeLessThanOrEqual(1 + 1e-12);

    const onSegment = add(s.p[0],
        mul(result.parameter, sub(s.p[1], s.p[0])));
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onSegment.values[i], 6);
    }

    // The closest box point must have |y[i]| <= e[i] in the box frame.
    const delta = sub(result.closest[1], b.center);
    for (let i = 0; i < 3; ++i) {
        expect(Math.abs(dot(b.axis[i], delta)))
            .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
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

describe('DistSegment3OrientedBox3', () => {
    const query = new DistSegment3OrientedBox3();
    const unitBox = alignedObb([0, 0, 0], [1, 1, 1]);

    it('measures a segment parallel to a face', () => {
        const s = segment([-2, 0, 4], [2, 0, 4]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 10);
        expectConsistent(result, s, unitBox);
    });

    it('reports zero distance for a segment passing through the box', () => {
        const s = segment([-3, 0, 0], [3, 0, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitBox);
    });

    it('reports zero distance for a segment strictly inside the box', () => {
        const s = segment([-0.5, 0.25, 0.1], [0.5, -0.25, -0.1]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitBox);
    });

    it('reports zero distance when an endpoint touches a face', () => {
        const s = segment([1, 0, 0], [4, 0, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitBox);
    });

    it('clamps to the first endpoint when the line minimum is behind it',
        () => {
            const s = segment([3, 0, 0], [6, 0, 0]);
            const result = query.compute(s, unitBox);
            expect(result.parameter).toBe(0);
            expect(result.distance).toBeCloseTo(2, 10);
            expect(result.closest[1].values[0]).toBeCloseTo(1, 10);
            expectConsistent(result, s, unitBox);
        });

    it('clamps to the second endpoint when the line minimum is beyond it',
        () => {
            const s = segment([6, 0, 0], [3, 0, 0]);
            const result = query.compute(s, unitBox);
            expect(result.parameter).toBe(1);
            expect(result.distance).toBeCloseTo(2, 10);
            expectConsistent(result, s, unitBox);
        });

    it('measures the distance to a box corner', () => {
        const s = segment([3, 3, 3], [4, 5, 6]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(Math.sqrt(12), 10);
        expect(result.parameter).toBe(0);
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i]).toBeCloseTo(1, 10);
        }
        expectConsistent(result, s, unitBox);
    });

    it('measures a rotated box analytically', () => {
        // A box rotated 45 degrees about z, with extents (1,1,1). Its
        // "corner" nearest to +x lies at distance sqrt(2) from the center
        // along x. A segment parallel to the y axis at x = 5 is therefore at
        // distance 5 - sqrt(2).
        const c = Math.SQRT1_2;
        const box45 = OrientedBox.fromCenterAxisExtent(v(0, 0, 0),
            [v(c, c, 0), v(-c, c, 0), v(0, 0, 1)], v(1, 1, 1));
        const s = segment([5, -3, 0], [5, 3, 0]);
        const result = query.compute(s, box45);
        expect(result.distance).toBeCloseTo(5 - Math.SQRT2, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(Math.SQRT2, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 10);
        expectConsistent(result, s, box45);
    });

    it('is translation and rotation invariant', () => {
        const axes = frame(0.7, -1.1);
        const rotated = OrientedBox.fromCenterAxisExtent(v(3, -2, 5), axes,
            v(1, 2, 0.5));
        const cbox = CanonicalBox.fromExtent(v(1, 2, 0.5));
        const cQuery = new DistSegment3CanonicalBox3();
        const random = makeRandom(20260901);
        for (let trial = 0; trial < 60; ++trial) {
            // Build a segment in the box frame, then map it to world space.
            const local: Vector[] = [];
            for (let k = 0; k < 2; ++k) {
                local.push(v(6 * random() - 3, 6 * random() - 3,
                    6 * random() - 3));
            }
            const world = local.map((q) => {
                let w = rotated.center.clone();
                for (let i = 0; i < 3; ++i) {
                    w = add(w, mul(q.values[i], axes[i]));
                }
                return w;
            });
            const expected = cQuery.compute(
                Segment.fromEndpoints(local[0], local[1]), cbox);
            const actual = query.compute(
                Segment.fromEndpoints(world[0], world[1]), rotated);
            expect(actual.distance).toBeCloseTo(expected.distance, 9);
            expect(actual.parameter).toBeCloseTo(expected.parameter, 9);
        }
    });

    it('gives the same distance for both segment orientations', () => {
        const axes = frame(0.35, 0.9);
        const b = OrientedBox.fromCenterAxisExtent(v(1, 1, -1), axes,
            v(2, 0.5, 1));
        const cases: Array<[number[], number[]]> = [
            [[3, 0, 0], [6, 0, 0]],
            [[-2, 0, 4], [2, 0, 4]],
            [[-3, -3, -3], [3, 3, 3]],
            [[2, 2, -5], [2, 2, 5]],
            [[0.5, 0.5, 0.5], [7, 1, -3]]
        ];
        for (const [p0, p1] of cases) {
            const forward = query.compute(segment(p0, p1), b);
            const backward = query.compute(segment(p1, p0), b);
            // Only the distance is orientation independent; when several
            // segment points realize the minimum, the reported parameter
            // depends on the traversal direction.
            expect(forward.distance).toBeCloseTo(backward.distance, 10);
        }
    });

    it('handles a degenerate zero-length segment outside the box', () => {
        const axes = frame(1.2, 0.4);
        const b = OrientedBox.fromCenterAxisExtent(v(0, 1, 0), axes,
            v(1, 2, 3));
        const s = segment([3, 4, 5], [3, 4, 5]);
        const result = query.compute(s, b);
        const expected = new DistPointOrientedBox().compute(v(3, 4, 5), b);
        expect(result.distance).toBeCloseTo(expected.distance, 12);
        expect(result.closest[0].values[0]).toBeCloseTo(3, 12);
        expectConsistent(result, s, b);
    });

    it('handles a degenerate zero-length segment inside the box', () => {
        const s = segment([0.25, -0.5, 0.75], [0.25, -0.5, 0.75]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitBox);
    });

    it('handles a degenerate (flat) box', () => {
        const flat = alignedObb([0, 0, 0], [2, 3, 0]);
        const s = segment([0, 0, 1], [0, 0, 5]);
        const result = query.compute(s, flat);
        expect(result.distance).toBeCloseTo(1, 10);
        expect(result.parameter).toBe(0);
        expectConsistent(result, s, flat);
    });

    it('agrees with a brute-force sampling for random configurations', () => {
        const random = makeRandom(87042);
        for (let trial = 0; trial < 250; ++trial) {
            const axes = frame(2 * Math.PI * random(), Math.PI * random());
            const b = OrientedBox.fromCenterAxisExtent(
                v(4 * random() - 2, 4 * random() - 2, 4 * random() - 2),
                axes,
                v(0.2 + 2 * random(), 0.2 + 2 * random(),
                    0.2 + 2 * random()));
            const s = segment(
                [8 * random() - 4, 8 * random() - 4, 8 * random() - 4],
                [8 * random() - 4, 8 * random() - 4, 8 * random() - 4]);
            const result = query.compute(s, b);
            const expected = Math.sqrt(bruteForceSqrDistance(s, b));
            expect(Math.abs(result.distance - expected)).toBeLessThan(1e-6);
            expectConsistent(result, s, b);
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistSegment3OrientedBox3.h. Upstream solves the query on the line containing the
// segment and clamps the result to the segment domain, so the properties below
// check the line/point agreement as well as the geometric invariants.
// ---------------------------------------------------------------------------

const v21Shape = fc.tuple(wellScaledVector(3, -8, 8), rotationFrame(3),
    fc.array(positive(4, 0.1), { minLength: 3, maxLength: 3 }))
    .map(([c, axes, e]) => OrientedBox.fromCenterAxisExtent(c, axes,
        Vector.fromArray(e)));

// A box together with one of its points.
const v21ShapePoint = fc.tuple(v21Shape,
    fc.array(fc.double({ min: -1, max: 1, noNaN: true }),
        { minLength: 3, maxLength: 3 }))
    .map(([b, u]) => {
        let q = b.center.clone();
        for (let i = 0; i < 3; ++i) {
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
    for (let i = 0; i < 3; ++i) {
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

// Upstream's line-box queries accumulate the squared distance as
// "... + delta * parameter" with parameter = -delta / lenSqr, i.e. as a
// subtraction of two nearly equal quantities. A line that nearly touches the
// box can therefore produce a tiny negative sqrDistance, whose square root is
// NaN. Upstream has the identical expression (DistLine3CanonicalBox3.h and
// its 2D counterpart), so the port inherits it; see the API notes of the V21
// verification. The properties skip results with a non-finite distance rather
// than paper over it.
function v21Usable(res: { distance: number, sqrDistance: number }): boolean {
    return Number.isFinite(res.distance) && res.sqrDistance >= 0;
}

describe('DistSegment3OrientedBox3 verification', () => {
    const query = new DistSegment3OrientedBox3();

    it('result is self consistent and the points lie on their primitives',
        () => {
            check(fc.tuple(v21Segment, v21Shape), ([seg, s]) => {
                const res = query.compute(seg, s);
                if (!v21Usable(res)) {
                    return;
                }
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
            if (!v21Usable(res)) {
                return;
            }
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
                const lr = new DistLine3OrientedBox3().compute(line, s);
                const sr = query.compute(seg, s);
                if (!v21Usable(lr) || !v21Usable(sr)) {
                    return;
                }
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
                if (!v21Usable(zres)) {
                    return;
                }
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
            if (!v21Usable(r0) || !v21Usable(r1)) {
                return;
            }
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
        });
    });

    it('a short segment agrees with the point query at its first endpoint',
        () => {
            check(fc.tuple(wellScaledVector(3, -8, 8), unitVector(3),
                v21Shape), ([p0, d, s]) => {
                const seg = Segment.fromEndpoints(p0, add(p0, mul(1e-9, d)));
                const res = query.compute(seg, s);
                if (!v21Usable(res)) {
                    return;
                }
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
