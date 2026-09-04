import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistRay2AlignedBox2 } from '../src/DistRay2AlignedBox2.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, sub }
    from '../src/Vector.js';
import { DistLine2AlignedBox2 } from '../src/DistLine2AlignedBox2.js';
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

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
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

// The distance from a ray point to the solid box is a convex function of the
// ray parameter, so a ternary search over [0,tmax] converges to the true
// minimum independently of the query under test.
function bruteForce(r: Ray, b: AlignedBox, tmax: number): number {
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
function verifyClosest(r: Ray, b: AlignedBox,
    result: { distance: number, parameter: number, closest: [Vector, Vector] }
): void {
    expect(result.parameter).toBeGreaterThanOrEqual(0);
    const onRay = add(r.origin, mul(result.parameter, r.direction));
    for (let i = 0; i < onRay.size; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onRay.values[i], 10);
        expect(result.closest[1].values[i])
            .toBeGreaterThanOrEqual(b.min.values[i] - 1e-10);
        expect(result.closest[1].values[i])
            .toBeLessThanOrEqual(b.max.values[i] + 1e-10);
    }
    const diff = sub(result.closest[0], result.closest[1]);
    let len = 0;
    for (let i = 0; i < diff.size; ++i) {
        len += diff.values[i] * diff.values[i];
    }
    expect(Math.sqrt(len)).toBeCloseTo(result.distance, 10);
}

describe('DistRay2AlignedBox2', () => {
    const query = new DistRay2AlignedBox2();
    const unitBox = box([-1, -1], [1, 1]);

    it('returns the ray origin when the ray points away from the box', () => {
        const r = ray([3, 0], [1, 0]);
        const result = query.compute(r, unitBox);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.sqrDistance).toBeCloseTo(4, 12);
        expect(result.parameter).toBe(0);
        expect(result.closest[0].values).toEqual([3, 0]);
        expect(result.closest[1].values).toEqual([1, 0]);
    });

    it('reports zero distance when the ray meets the box', () => {
        const r = ray([-5, 0], [1, 0]);
        const result = query.compute(r, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.sqrDistance).toBeCloseTo(0, 12);
        verifyClosest(r, unitBox, result);
    });

    it('reports zero distance when the ray origin is inside the box', () => {
        const r = ray([0.25, -0.5], [1, 2]);
        const result = query.compute(r, unitBox);
        expect(result.distance).toBe(0);
        expect(result.parameter).toBeGreaterThanOrEqual(0);
        verifyClosest(r, unitBox, result);
    });

    it('handles a ray parallel to a box face', () => {
        const r = ray([-4, 3], [1, 0]);
        const result = query.compute(r, unitBox);
        expect(result.distance).toBeCloseTo(2, 12);
        verifyClosest(r, unitBox, result);
    });

    it('matches the analytic distance for a diagonal approach', () => {
        // The line through (4,0) with direction (-1,1) is x + y = 4. The
        // closest box point is the corner (1,1) at distance |1+1-4|/sqrt(2).
        const r = ray([4, 0], [-1, 1]);
        const result = query.compute(r, unitBox);
        expect(result.distance).toBeCloseTo(2 / Math.SQRT2, 10);
        verifyClosest(r, unitBox, result);
    });

    it('handles a degenerate box that is a single point', () => {
        const b = box([2, 2], [2, 2]);
        const r = ray([0, 0], [1, 0]);
        const result = query.compute(r, b);
        expect(result.distance).toBeCloseTo(2, 10);
        expect(result.closest[1].values).toEqual([2, 2]);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 987654321;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 200; ++trial) {
            const b = box([-1 - 2 * rand(), -1 - 2 * rand()],
                [1 + 2 * rand(), 1 + 2 * rand()]);
            const r = ray([8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1]);
            if (r.direction.values[0] === 0 && r.direction.values[1] === 0) {
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
// header DistRay2AlignedBox2.h. Upstream solves the query on the line containing the
// ray and clamps the result to the ray domain, so the properties below
// check the line/point agreement as well as the geometric invariants.
// ---------------------------------------------------------------------------

const v21Shape = fc.tuple(wellScaledVector(2, -8, 8),
    wellScaledVector(2, -8, 8)).map(([a, b]) => {
    const lo = new Vector(2);
    const hi = new Vector(2);
    for (let i = 0; i < 2; ++i) {
        lo.values[i] = Math.min(a.values[i], b.values[i]);
        hi.values[i] = Math.max(a.values[i], b.values[i]);
    }
    return AlignedBox.fromMinMax(lo, hi);
});

// A box together with one of its points.
const v21ShapePoint = fc.tuple(v21Shape,
    fc.array(fc.double({ min: 0, max: 1, noNaN: true }),
        { minLength: 2, maxLength: 2 }))
    .map(([b, u]) => {
        const q = new Vector(2);
        for (let i = 0; i < 2; ++i) {
            q.values[i] = (1 - u[i]) * b.min.values[i] + u[i] * b.max.values[i];
        }
        return [b, q] as [AlignedBox, Vector];
    });

function v21PointDistance(p: Vector, b: AlignedBox): number {
    return new DistPointAlignedBox().compute(p, b).distance;
}

function v21CheckShapePoint(b: AlignedBox,
    res: { closest: [Vector, Vector] }): void {
    for (let i = 0; i < 2; ++i) {
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

const v21Ray = fc.tuple(wellScaledVector(2, -8, 8), unitVector(2))
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

describe('DistRay2AlignedBox2 verification', () => {
    const query = new DistRay2AlignedBox2();

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
                const lr = new DistLine2AlignedBox2().compute(line, s);
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
        check(fc.tuple(v21ShapePoint, unitVector(2)), ([[s, q], d]) => {
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

    it('is invariant under a common translation', () => {
        // An aligned box is not closed under rotation, so only translations
        // move the whole configuration rigidly.
        check(fc.tuple(v21Ray, v21Shape, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([r, s, R, tr]) => {
            void R;
            const rot = (x: Vector): Vector => x.clone();
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
