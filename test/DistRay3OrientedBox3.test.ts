import { describe, expect, it } from 'vitest';
import { DistRay3OrientedBox3 } from '../src/DistRay3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, sub }
    from '../src/Vector.js';
import { DistLine3OrientedBox3 } from '../src/DistLine3OrientedBox3.js';
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

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

// The distance from a ray point to the solid box is a convex function of the
// ray parameter, so a ternary search over [0,tmax] converges to the true
// minimum independently of the query under test.
function bruteForce(r: Ray, b: OrientedBox, tmax: number): number {
    const f = (t: number): number =>
        pointBoxSqrDistance(add(r.origin, mul(t, r.direction)), b);
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
    return Math.sqrt(f(0.5 * (lo + hi)));
}

// Verify that the reported closest points are consistent with the reported
// distance and lie on their primitives.
function verifyClosest(r: Ray, b: OrientedBox,
    result: { distance: number, parameter: number, closest: [Vector, Vector] }
): void {
    expect(result.parameter).toBeGreaterThanOrEqual(0);
    const onRay = add(r.origin, mul(result.parameter, r.direction));
    for (let i = 0; i < onRay.size; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onRay.values[i], 9);
    }
    verifyBoxPoint(b, result.closest[1]);
    expect(length(sub(result.closest[0], result.closest[1])))
        .toBeCloseTo(result.distance, 9);
}

describe('DistRay3OrientedBox3', () => {
    const query = new DistRay3OrientedBox3();
    const axisAligned = box([0, 0, 0],
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]], [1, 1, 1]);
    const c = Math.SQRT1_2;
    const rotated = box([0, 0, 0],
        [[c, c, 0], [-c, c, 0], [0, 0, 1]], [1, 1, 1]);

    it('matches the aligned result for an axis-aligned oriented box', () => {
        const r = ray([4, 0, 0], [1, 0, 0]);
        const result = query.compute(r, axisAligned);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.parameter).toBe(0);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
    });

    it('handles a box rotated about the z-axis', () => {
        const r = ray([4, 0, 0], [1, 0, 0]);
        const result = query.compute(r, rotated);
        expect(result.distance).toBeCloseTo(4 - Math.SQRT2, 10);
        verifyClosest(r, rotated, result);
    });

    it('reports zero distance when the ray enters the box', () => {
        const r = ray([0, 0, -5], [0, 0, 1]);
        const result = query.compute(r, rotated);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('handles a degenerate box with zero extents', () => {
        const b = box([0, 3, 0], [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            [0, 0, 0]);
        const r = ray([0, 0, 0], [1, 0, 0]);
        const result = query.compute(r, b);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.parameter).toBe(0);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 55667788;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 150; ++trial) {
            const angle = 2 * Math.PI * rand();
            const ca = Math.cos(angle);
            const sa = Math.sin(angle);
            const b = box([2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1],
                [[ca, sa, 0], [-sa, ca, 0], [0, 0, 1]],
                [0.2 + rand(), 0.2 + rand(), 0.2 + rand()]);
            const r = ray([8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            if (length(r.direction) < 1e-6) {
                continue;
            }
            const result = query.compute(r, b);
            const brute = bruteForce(r, b, 1e6);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-6);
            expect(result.distance).toBeCloseTo(brute, 6);
            verifyClosest(r, b, result);
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistRay3OrientedBox3.h. Upstream solves the query on the line containing the
// ray and clamps the result to the ray domain, so the properties below
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

describe('DistRay3OrientedBox3 verification', () => {
    const query = new DistRay3OrientedBox3();

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
                const lr = new DistLine3OrientedBox3().compute(line, s);
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
