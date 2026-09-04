import { describe, expect, it } from 'vitest';
import { DistRay3Rectangle3 } from '../src/DistRay3Rectangle3.js';
import { Ray } from '../src/Ray.js';
import { Rectangle } from '../src/Rectangle.js';
import { Vector, add, dot, length, mul, sub }
    from '../src/Vector.js';
import { DistLine3Rectangle3 } from '../src/DistLine3Rectangle3.js';
import { DistPointRectangle } from '../src/DistPointRectangle.js';
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

function rectangle(center: number[], axis: number[][],
    extent: number[]): Rectangle {
    return Rectangle.fromCenterAxisExtent(v(...center),
        axis.map(a => v(...a)), v(...extent));
}

// The exact squared distance from a point to the solid rectangle, computed
// independently of the library (the axes are orthonormal, so the rectangle
// coordinates are clamped independently).
function pointRectangleSqrDistance(p: Vector, r: Rectangle): number {
    const delta = sub(p, r.center);
    let closest = r.center.clone();
    let s0 = dot(delta, r.axis[0]);
    let s1 = dot(delta, r.axis[1]);
    s0 = Math.max(-r.extent.values[0], Math.min(r.extent.values[0], s0));
    s1 = Math.max(-r.extent.values[1], Math.min(r.extent.values[1], s1));
    closest = add(closest, add(mul(s0, r.axis[0]), mul(s1, r.axis[1])));
    const diff = sub(p, closest);
    return dot(diff, diff);
}

// The distance from a ray point to the solid rectangle is a convex function
// of the ray parameter, so a ternary search over [0,tmax] converges to the
// true minimum independently of the query under test.
function bruteForce(r: Ray, rect: Rectangle, tmax: number): number {
    const f = (t: number): number => pointRectangleSqrDistance(
        add(r.origin, mul(t, r.direction)), rect);
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

function verifyClosest(r: Ray, rect: Rectangle,
    result: {
        distance: number, parameter: number,
        cartesian: [number, number], closest: [Vector, Vector]
    }): void {
    expect(result.parameter).toBeGreaterThanOrEqual(0);
    const onRay = add(r.origin, mul(result.parameter, r.direction));
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onRay.values[i], 9);
    }
    // The reported cartesian coordinates must be in range and reproduce the
    // reported rectangle point.
    expect(Math.abs(result.cartesian[0]))
        .toBeLessThanOrEqual(rect.extent.values[0] + 1e-9);
    expect(Math.abs(result.cartesian[1]))
        .toBeLessThanOrEqual(rect.extent.values[1] + 1e-9);
    const onRect = add(rect.center,
        add(mul(result.cartesian[0], rect.axis[0]),
            mul(result.cartesian[1], rect.axis[1])));
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[1].values[i]).toBeCloseTo(onRect.values[i], 9);
    }
    expect(length(sub(result.closest[0], result.closest[1])))
        .toBeCloseTo(result.distance, 9);
}

describe('DistRay3Rectangle3', () => {
    const query = new DistRay3Rectangle3();
    const unitRect = rectangle([0, 0, 0], [[1, 0, 0], [0, 1, 0]], [1, 1]);

    it('returns the ray origin when the ray points away', () => {
        const r = ray([0, 0, 5], [0, 0, 1]);
        const result = query.compute(r, unitRect);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.sqrDistance).toBeCloseTo(25, 12);
        expect(result.parameter).toBe(0);
        expect(result.cartesian).toEqual([0, 0]);
        expect(result.closest[1].values).toEqual([0, 0, 0]);
    });

    it('reports zero distance when the ray hits the rectangle', () => {
        const r = ray([0.25, -0.5, 5], [0, 0, -1]);
        const result = query.compute(r, unitRect);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.cartesian[0]).toBeCloseTo(0.25, 10);
        expect(result.cartesian[1]).toBeCloseTo(-0.5, 10);
    });

    it('finds the rectangle edge when the ray passes beside it', () => {
        // The ray is parallel to the rectangle plane at height 4 and moves
        // away from the rectangle in x.
        const r = ray([3, 0, 4], [1, 0, 0]);
        const result = query.compute(r, unitRect);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 16), 10);
        expect(result.parameter).toBe(0);
        expect(result.cartesian[0]).toBeCloseTo(1, 10);
        verifyClosest(r, unitRect, result);
    });

    it('handles a ray whose line meets the rectangle behind the origin',
        () => {
            const r = ray([0, 0, 2], [0, 0, 1]);
            const result = query.compute(r, unitRect);
            expect(result.distance).toBeCloseTo(2, 12);
            expect(result.parameter).toBe(0);
        });

    it('handles a degenerate rectangle with zero extents', () => {
        const rect = rectangle([0, 0, 0], [[1, 0, 0], [0, 1, 0]], [0, 0]);
        const r = ray([0, 0, 3], [0, 0, 1]);
        const result = query.compute(r, rect);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values).toEqual([0, 0, 0]);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 20240817;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 150; ++trial) {
            const angle = 2 * Math.PI * rand();
            const c = Math.cos(angle);
            const s = Math.sin(angle);
            const rect = rectangle([2 * rand() - 1, 2 * rand() - 1,
                2 * rand() - 1], [[c, s, 0], [-s, c, 0]],
            [0.2 + rand(), 0.2 + rand()]);
            const r = ray([8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            if (length(r.direction) < 1e-6) {
                continue;
            }
            const result = query.compute(r, rect);
            const brute = bruteForce(r, rect, 1e6);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-6);
            expect(result.distance).toBeCloseTo(brute, 6);
            verifyClosest(r, rect, result);
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistRay3Rectangle3.h. Upstream solves the query on the line containing the
// ray and clamps the result to the ray domain, so the properties below
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

describe('DistRay3Rectangle3 verification', () => {
    const query = new DistRay3Rectangle3();

    it('result is self consistent and the points lie on their primitives',
        () => {
            check(fc.tuple(v21Ray, v21Shape), ([r, s]) => {
                const res = query.compute(r, s);
                if (!v21Usable(res)) {
                    return;
                }
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
            if (!v21Usable(res)) {
                return;
            }
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
                const lr = new DistLine3Rectangle3().compute(line, s);
                const rr = query.compute(r, s);
                if (!v21Usable(lr) || !v21Usable(rr)) {
                    return;
                }
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
            if (!v21Usable(zres)) {
                return;
            }
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
            if (!v21Usable(r0) || !v21Usable(r1)) {
                return;
            }
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
