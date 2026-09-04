import { describe, expect, it } from 'vitest';
import { DistRectangle3Rectangle3 } from '../src/DistRectangle3Rectangle3.js';
import type { DistRectangle3Rectangle3Result }
    from '../src/DistRectangle3Rectangle3.js';
import { Rectangle } from '../src/Rectangle.js';
import { DistSegment3Rectangle3 } from '../src/DistSegment3Rectangle3.js';
import { DistTriangle3Triangle3 } from '../src/DistTriangle3Triangle3.js';
import { Segment } from '../src/Segment.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    seededRandom, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

// An orthonormal frame parameterized by two angles.
function frame(a: number, b: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    return [
        v(ca, sa, 0),
        v(-sa * cb, ca * cb, sb),
        v(sa * sb, -ca * sb, cb)
    ];
}

function rect(center: number[], axis: Vector[], e0: number, e1: number):
    Rectangle {
    return Rectangle.fromCenterAxisExtent(v(...center), [axis[0], axis[1]],
        v(e0, e1));
}

// The exact squared distance from a point to a solid rectangle: clamp the
// rectangle coordinates of the point to the extents.
function pointRectangleSqr(p: Vector, r: Rectangle): number {
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

// A grid sampling of rectangle0 combined with the exact point-rectangle
// distance to rectangle1 gives an upper bound whose error is at most the
// grid spacing.
function sampledDistance(r0: Rectangle, r1: Rectangle, n: number): number {
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= n; ++i) {
        const s0 = r0.extent.values[0] * (-1 + (2 * i) / n);
        const base = add(r0.center, mul(s0, r0.axis[0]));
        for (let j = 0; j <= n; ++j) {
            const s1 = r0.extent.values[1] * (-1 + (2 * j) / n);
            best = Math.min(best,
                pointRectangleSqr(add(base, mul(s1, r0.axis[1])), r1));
        }
    }
    return Math.sqrt(best);
}

function expectConsistent(r0: Rectangle, r1: Rectangle,
    result: DistRectangle3Rectangle3Result): void {
    const delta = sub(result.closest[0], result.closest[1]);
    expect(Math.sqrt(dot(delta, delta))).toBeCloseTo(result.distance, 8);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 8);

    const check = (c: [number, number], r: Rectangle, closest: Vector) => {
        for (let i = 0; i < 2; ++i) {
            expect(Math.abs(c[i]))
                .toBeLessThanOrEqual(r.extent.values[i] + 1e-8);
        }
        const p = add(r.center,
            add(mul(c[0], r.axis[0]), mul(c[1], r.axis[1])));
        for (let i = 0; i < 3; ++i) {
            expect(p.values[i]).toBeCloseTo(closest.values[i], 7);
        }
    };
    check(result.cartesian0, r0, result.closest[0]);
    check(result.cartesian1, r1, result.closest[1]);
}

describe('DistRectangle3Rectangle3', () => {
    const query = new DistRectangle3Rectangle3();

    it('computes the distance between parallel overlapping rectangles', () => {
        const r0 = rect([0, 0, 0], frame(0, 0), 1, 1);
        const r1 = rect([0, 0, 7], frame(0, 0), 1, 1);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(7, 10);
        expect(result.closest[0].values[2]).toBeCloseTo(0, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(7, 10);
        expectConsistent(r0, r1, result);
    });

    it('computes the distance between coplanar separated rectangles', () => {
        const r0 = rect([0, 0, 0], frame(0, 0), 1, 1);
        const r1 = rect([5, 0, 0], frame(0, 0), 1, 1);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.cartesian0[0]).toBeCloseTo(1, 8);
        expect(result.cartesian1[0]).toBeCloseTo(-1, 8);
        expectConsistent(r0, r1, result);
    });

    it('reports zero distance for crossing rectangles', () => {
        const r0 = rect([0, 0, 0], [v(1, 0, 0), v(0, 1, 0)], 2, 2);
        const r1 = rect([0, 0, 0], [v(1, 0, 0), v(0, 0, 1)], 2, 2);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(0, 10);
        expectConsistent(r0, r1, result);
    });

    it('finds a corner-corner closest pair', () => {
        const r0 = rect([0, 0, 0], frame(0, 0), 1, 1);
        const r1 = rect([5, 5, 0], frame(0, 0), 1, 1);
        const result = query.compute(r0, r1);
        expect(result.distance).toBeCloseTo(Math.sqrt(9 + 9), 8);
        expect(result.cartesian0[0]).toBeCloseTo(1, 8);
        expect(result.cartesian0[1]).toBeCloseTo(1, 8);
        expect(result.cartesian1[0]).toBeCloseTo(-1, 8);
        expect(result.cartesian1[1]).toBeCloseTo(-1, 8);
        expectConsistent(r0, r1, result);
    });

    it('is symmetric in its arguments', () => {
        const r0 = rect([0.3, -1.2, 0.4], frame(0.5, 1.1), 1.5, 0.6);
        const r1 = rect([-3.1, 1.4, 2.9], frame(2.0, 0.3), 0.8, 1.7);
        const r01 = query.compute(r0, r1);
        const r10 = query.compute(r1, r0);
        expect(r10.distance).toBeCloseTo(r01.distance, 9);
        expectConsistent(r0, r1, r01);
        expectConsistent(r1, r0, r10);
    });

    it('matches a dense sampling on random inputs', () => {
        let seed = 135792468;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rv = () => [rand() * 6 - 3, rand() * 6 - 3, rand() * 6 - 3];

        for (let trial = 0; trial < 20; ++trial) {
            const r0 = rect(rv(), frame(rand() * Math.PI, rand() * Math.PI),
                rand() * 2 + 0.1, rand() * 2 + 0.1);
            const r1 = rect(rv(), frame(rand() * Math.PI, rand() * Math.PI),
                rand() * 2 + 0.1, rand() * 2 + 0.1);
            const result = query.compute(r0, r1);
            const sampled = Math.min(sampledDistance(r0, r1, 70),
                sampledDistance(r1, r0, 70));
            expect(result.distance).toBeLessThanOrEqual(sampled + 1e-8);
            expect(sampled - result.distance).toBeLessThan(0.2);
            expectConsistent(r0, r1, result);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistRectangle3Rectangle3.h.
// ---------------------------------------------------------------------------

describe('DistRectangle3Rectangle3 verification', () => {
    const query = new DistRectangle3Rectangle3();
    const ttQuery = new DistTriangle3Triangle3();
    const srQuery = new DistSegment3Rectangle3();

    const rectArb = fc.tuple(wellScaledVector(3, -5, 5), rotationFrame(3),
        fc.array(finite(0.05, 3), { minLength: 2, maxLength: 2 }))
        .map(([c, R, e]) =>
            Rectangle.fromCenterAxisExtent(c, [R[0], R[1]], v(e[0], e[1])));

    // The two triangles whose union is the solid rectangle. The vertices come
    // back in bit-pattern order, so <0,1,3> and <0,3,2> tile the rectangle.
    function triangles(r: Rectangle): [Triangle, Triangle] {
        const p = r.getVertices();
        return [Triangle.fromVertices(p[0], p[1], p[3]),
            Triangle.fromVertices(p[0], p[3], p[2])];
    }

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(rectArb, rectArb), ([r0, r1]) => {
            const res = query.compute(r0, r1);
            expectClose(res.sqrDistance, res.distance * res.distance,
                1e-12, 1e-12);
            const d = sub(res.closest[0], res.closest[1]);
            expectClose(Math.sqrt(dot(d, d)), res.distance, 1e-8, 1e-8);

            const verify = (s: [number, number], r: Rectangle, c: Vector) => {
                let rebuilt = r.center.clone();
                for (let i = 0; i < 2; ++i) {
                    expect(Math.abs(s[i]))
                        .toBeLessThanOrEqual(r.extent.values[i] + 1e-8);
                    rebuilt = add(rebuilt, mul(s[i], r.axis[i]));
                }
                expectVectorClose(rebuilt, c, 1e-7, 1e-7);
            };
            verify(res.cartesian0, r0, res.closest[0]);
            verify(res.cartesian1, r1, res.closest[1]);
        });
    });

    it('is symmetric under argument swap', () => {
        check(fc.tuple(rectArb, rectArb), ([r0, r1]) => {
            const a = query.compute(r0, r1);
            const b = query.compute(r1, r0);
            expectClose(a.distance, b.distance, 1e-8, 1e-8);
            const d = sub(b.closest[1], b.closest[0]);
            expectClose(Math.sqrt(dot(d, d)), a.distance, 1e-8, 1e-8);
        });
    });

    // A rectangle is the union of two triangles, so the rectangle-rectangle
    // distance is the minimum over the four triangle pairs. This is an
    // independent code path (DistTriangle3Triangle3 does not use the
    // segment-rectangle query at all).
    it('equals the minimum over the four triangle pairs', () => {
        check(fc.tuple(rectArb, rectArb), ([r0, r1]) => {
            const a = triangles(r0), b = triangles(r1);
            let best = Number.MAX_VALUE;
            for (const t0 of a) {
                for (const t1 of b) {
                    best = Math.min(best, ttQuery.compute(t0, t1).sqrDistance);
                }
            }
            const res = query.compute(r0, r1);
            expectClose(res.sqrDistance, best, 1e-8, 1e-8);
        });
    });

    it('equals the minimum over the eight edge-rectangle queries', () => {
        check(fc.tuple(rectArb, rectArb), ([r0, r1]) => {
            const edges: [number, number][] =
                [[0, 1], [2, 3], [0, 2], [1, 3]];
            const p0 = r0.getVertices(), p1 = r1.getVertices();
            let best = Number.MAX_VALUE;
            for (const [i, j] of edges) {
                best = Math.min(best, srQuery.compute(
                    Segment.fromEndpoints(p0[i], p0[j]), r1).sqrDistance);
                best = Math.min(best, srQuery.compute(
                    Segment.fromEndpoints(p1[i], p1[j]), r0).sqrDistance);
            }
            const res = query.compute(r0, r1);
            expectClose(res.sqrDistance, best, 1e-9, 1e-9);
        });
    });

    it('brackets a grid sampling of the first rectangle', () => {
        const rng = seededRandom(0x7e57ca5e);
        const n = 40;
        for (let k = 0; k < 25; ++k) {
            const mk = () => rect([6 * rng() - 3, 6 * rng() - 3,
                6 * rng() - 3], frame(2 * Math.PI * rng(), 2 * Math.PI * rng()),
                0.5 + rng(), 0.5 + rng());
            const r0 = mk(), r1 = mk();
            const res = query.compute(r0, r1);
            const s = sampledDistance(r0, r1, n);
            expect(res.distance).toBeLessThanOrEqual(s + 1e-9);
            // The sample spacing along axis i is 2*extent[i]/n, so the
            // 1-Lipschitz distance function overshoots by at most the
            // diagonal of one cell.
            const cell = Math.sqrt(
                (2 * r0.extent.values[0] / n) ** 2
                + (2 * r0.extent.values[1] / n) ** 2);
            expect(s - res.distance).toBeLessThanOrEqual(cell + 1e-9);
        }
    }, 30000);

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(rectArb, rectArb, rotationFrame(3),
            wellScaledVector(3, -4, 4)), ([r0, r1, R, t]) => {
                const rot = (q: Vector) => add(add(
                    mul(q.values[0], R[0]), mul(q.values[1], R[1])),
                    mul(q.values[2], R[2]));
                const xr = (r: Rectangle) => Rectangle.fromCenterAxisExtent(
                    add(rot(r.center), t),
                    [rot(r.axis[0]), rot(r.axis[1])], r.extent);
                const a = query.compute(r0, r1);
                const b = query.compute(xr(r0), xr(r1));
                expectClose(a.distance, b.distance, 1e-8, 1e-8);
            });
    });

    it('handles a zero-extent (point) rectangle', () => {
        check(fc.tuple(rectArb, wellScaledVector(3, -5, 5), rotationFrame(3)),
            ([r0, c, R]) => {
                const pt = Rectangle.fromCenterAxisExtent(c, [R[0], R[1]],
                    v(0, 0));
                const res = query.compute(r0, pt);
                expectClose(res.sqrDistance, pointRectangleSqr(c, r0),
                    1e-8, 1e-8);
                expectVectorClose(res.closest[1], c, 1e-8, 1e-8);
            });
    });
});
