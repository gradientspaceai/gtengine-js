import { describe, expect, it } from 'vitest';
import { DistPointTriangle } from '../src/DistPointTriangle.js';
import { DistRay3Triangle3 } from '../src/DistRay3Triangle3.js';
import { Ray } from '../src/Ray.js';
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

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

function triangle(v0: number[], v1: number[], v2: number[]): Triangle {
    return Triangle.fromVertices(v(...v0), v(...v1), v(...v2));
}

// Brute-force minimum over a dense sampling of the ray parameter, using the
// independently tested point-triangle query for each sample.
// The distance from a ray point to the (convex) triangle is a convex function
// of the ray parameter, so a ternary search over [0,tmax] converges to the
// true minimum. The per-sample point-triangle distance uses the
// independently tested DistPointTriangle query.
function bruteForce(r: Ray, tri: Triangle, tmax: number): number {
    const ptQuery = new DistPointTriangle();
    const f = (t: number): number => ptQuery.compute(
        add(r.origin, mul(t, r.direction)), tri).distance;
    let lo = 0;
    let hi = tmax;
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
    return f(0.5 * (lo + hi));
}

function verifyClosest(r: Ray, tri: Triangle,
    result: {
        distance: number, parameter: number,
        barycentric: [number, number, number], closest: [Vector, Vector]
    }): void {
    expect(result.parameter).toBeGreaterThanOrEqual(0);
    const onRay = add(r.origin, mul(result.parameter, r.direction));
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onRay.values[i], 9);
    }
    // The barycentric coordinates must be nonnegative, sum to 1 and reproduce
    // the reported triangle point.
    let sum = 0;
    for (let i = 0; i < 3; ++i) {
        expect(result.barycentric[i]).toBeGreaterThanOrEqual(-1e-9);
        sum += result.barycentric[i];
    }
    expect(sum).toBeCloseTo(1, 9);
    let onTri = new Vector(3);
    for (let i = 0; i < 3; ++i) {
        onTri = add(onTri, mul(result.barycentric[i], tri.v[i]));
    }
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[1].values[i]).toBeCloseTo(onTri.values[i], 9);
    }
    expect(length(sub(result.closest[0], result.closest[1])))
        .toBeCloseTo(result.distance, 9);
}

describe('DistRay3Triangle3', () => {
    const query = new DistRay3Triangle3();
    const tri = triangle([0, 0, 0], [1, 0, 0], [0, 1, 0]);

    it('returns the ray origin when the ray points away', () => {
        const r = ray([0, 0, 5], [0, 0, 1]);
        const result = query.compute(r, tri);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.sqrDistance).toBeCloseTo(25, 12);
        expect(result.parameter).toBe(0);
        expect(result.barycentric[0]).toBeCloseTo(1, 10);
        expect(result.closest[1].values).toEqual([0, 0, 0]);
    });

    it('reports zero distance when the ray hits the triangle', () => {
        const r = ray([0.25, 0.25, 5], [0, 0, -1]);
        const result = query.compute(r, tri);
        expect(result.distance).toBeCloseTo(0, 10);
        expect(result.barycentric[1]).toBeCloseTo(0.25, 8);
        expect(result.barycentric[2]).toBeCloseTo(0.25, 8);
    });

    it('finds a vertex as the closest point', () => {
        const r = ray([3, 0, 0], [1, 0, 0]);
        const result = query.compute(r, tri);
        expect(result.distance).toBeCloseTo(2, 10);
        expect(result.parameter).toBe(0);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 10);
    });

    it('finds an edge point as the closest point', () => {
        const r = ray([0.5, -2, 0], [1, 0, 0]);
        const result = query.compute(r, tri);
        // The closest triangle point is on the edge from (0,0,0) to (1,0,0).
        expect(result.distance).toBeCloseTo(2, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 10);
        verifyClosest(r, tri, result);
    });

    it('handles a degenerate triangle collapsed to a segment', () => {
        const degenerate = triangle([0, 0, 0], [1, 0, 0], [2, 0, 0]);
        const r = ray([0, 3, 0], [1, 0, 0]);
        const result = query.compute(r, degenerate);
        expect(result.distance).toBeCloseTo(3, 10);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 31415926;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 120; ++trial) {
            const t = triangle(
                [3 * rand() - 1.5, 3 * rand() - 1.5, 3 * rand() - 1.5],
                [3 * rand() - 1.5, 3 * rand() - 1.5, 3 * rand() - 1.5],
                [3 * rand() - 1.5, 3 * rand() - 1.5, 3 * rand() - 1.5]);
            const r = ray([8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            if (length(r.direction) < 1e-6) {
                continue;
            }
            const result = query.compute(r, t);
            const brute = bruteForce(r, t, 1e6);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-6);
            expect(result.distance).toBeCloseTo(brute, 6);
            verifyClosest(r, t, result);
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistRay3Triangle3.h. Upstream solves the query on the line containing the
// ray and clamps the result to the ray domain, so the properties below
// check the line/point agreement as well as the geometric invariants.
// ---------------------------------------------------------------------------

// Triangles with a bounded aspect ratio: twice the area is at least a fixed
// fraction of the squared longest edge. An unconstrained area bound admits
// slivers (a 16-long, 0.001-thin triangle passes |e0 x e1| > 1e-2), and the
// barycentric coordinates of a sliver carry a 1/area error amplification that
// swamps the tolerances of the cross-checks below.
const v21Shape = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -8, 8), wellScaledVector(3, -8, 8))
    .filter(([a, b, c]) => {
        const e0 = sub(b, a);
        const e1 = sub(c, a);
        const e2 = sub(c, b);
        const maxEdge = Math.max(length(e0), length(e1), length(e2));
        return maxEdge > 1e-1
            && length(cross(e0, e1)) > 0.05 * maxEdge * maxEdge;
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

const v21Ray = fc.tuple(wellScaledVector(3, -8, 8), unitVector(3))
    .map(([o, d]) => Ray.fromOriginDirection(o, d));

// Minimum of a convex function on [lo,hi] by ternary search. The distance
// from a point to a convex set is convex and the ray is an affine image of
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

describe('DistRay3Triangle3 verification', () => {
    const query = new DistRay3Triangle3();

    it('result is self consistent and the points lie on their primitives',
        () => {
            check(fc.tuple(v21Ray, v21Shape), ([r, s]) => {
                const res = query.compute(r, s);
                expectClose(res.distance, Math.sqrt(res.sqrDistance), 1e-12,
                    1e-12);
                const diff = sub(res.closest[0], res.closest[1]);
                expectClose(res.sqrDistance, dot(diff, diff), 1e-9, 1e-9);
                expect(res.parameter).toBeGreaterThanOrEqual(0);
                expectVectorClose(res.closest[0],
                    add(r.origin, mul(res.parameter, r.direction)), 1e-9,
                    1e-9);
                v21CheckShapePoint(s, res);
            });
        });

    it('matches a convex minimization along the ray', () => {
        check(fc.tuple(v21Ray, v21Shape), ([r, s]) => {
            const res = query.compute(r, s);
            const f = (t: number): number =>
                v21PointDistance(add(r.origin, mul(t, r.direction)), s);
            expectClose(res.distance, v21MinOnInterval(f, 0, 100), 1e-7,
                1e-7);
        }, 60);
    }, 30000);

    it('agrees with the line query on the ray and with the point query off it',
        () => {
            check(fc.tuple(v21Ray, v21Shape), ([r, s]) => {
                const line = Line.fromOriginDirection(r.origin, r.direction);
                const lr = new DistLine3Triangle3().compute(line, s);
                const rr = query.compute(r, s);
                if (lr.parameter >= 0) {
                    expect(rr.parameter).toBe(lr.parameter);
                    expect(rr.distance).toBe(lr.distance);
                    expectVectorClose(rr.closest[0], lr.closest[0], 0, 0);
                    expectVectorClose(rr.closest[1], lr.closest[1], 0, 0);
                }
                else {
                    expect(rr.parameter).toBe(0);
                    expectVectorClose(rr.closest[0], r.origin, 0, 0);
                    expect(rr.distance)
                        .toBe(v21PointDistance(r.origin, s));
                }
            });
        });

    it('reports zero distance when the ray starts on the shape', () => {
        check(fc.tuple(v21ShapePoint, unitVector(3)), ([[s, q], d]) => {
            const r = Ray.fromOriginDirection(q, d);
            const zres = query.compute(r, s);
            // The ray contains q, so its distance to the shape is at most q's
            // own distance to the shape. That distance is zero in exact
            // arithmetic; the sampled shape point carries rounding error that
            // grows with the coordinate magnitudes and, for a thin triangle,
            // with 1/area, so the near-zero bound is scale relative.
            const pd = v21PointDistance(q, s);
            expect(zres.distance).toBeLessThanOrEqual(pd + 1e-9);
            expect(pd).toBeLessThanOrEqual(1e-7 * (1 + length(q)));
        });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(v21Ray, v21Shape, rotationFrame(3),
            wellScaledVector(3, -5, 5)), ([r, s, R, tr]) => {
            const rot = (x: Vector): Vector => {
                let y = new Vector(3);
                for (let i = 0; i < 3; ++i) {
                    y = add(y, mul(x.values[i], R[i]));
                }
                return y;
            };
            const movedRay = Ray.fromOriginDirection(add(rot(r.origin), tr),
                rot(r.direction));
            const r0 = query.compute(r, s);
            const r1 = query.compute(movedRay, v21MoveShape(s, rot, tr));
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
        });
    });

    it('does not mutate its inputs', () => {
        check(fc.tuple(v21Ray, v21Shape), ([r, s]) => {
            const o = r.origin.clone();
            const d = r.direction.clone();
            const snapshot = v21ShapeSnapshot(s);
            const res = query.compute(r, s);
            expect(r.origin.values).toEqual(o.values);
            expect(r.direction.values).toEqual(d.values);
            expect(v21ShapeSnapshot(s)).toEqual(snapshot);
            res.closest[0].values[0] = 4242;
            res.closest[1].values[0] = 4242;
            expect(r.origin.values).toEqual(o.values);
            expect(v21ShapeSnapshot(s)).toEqual(snapshot);
        });
    });
});
