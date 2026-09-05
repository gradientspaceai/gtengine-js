import { describe, expect, it } from 'vitest';
import { DistTriangle3Rectangle3 } from '../src/DistTriangle3Rectangle3.js';
import type { DistTriangle3Rectangle3Result }
    from '../src/DistTriangle3Rectangle3.js';
import { Rectangle } from '../src/Rectangle.js';
import { Triangle } from '../src/Triangle.js';
import { DistTriangle3Triangle3 } from '../src/DistTriangle3Triangle3.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, latticeVector,
    rotationFrame, seededRandom, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function tri(v0: number[], v1: number[], v2: number[]): Triangle {
    return Triangle.fromVertices(v(...v0), v(...v1), v(...v2));
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

// The exact distance from a point to a solid rectangle: clamp the rectangle
// coordinates of the point to the extents.
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

// A barycentric grid sampling of the triangle combined with the exact
// point-rectangle distance gives an upper bound whose error is at most the
// grid spacing.
function sampledDistance(t: Triangle, r: Rectangle, n: number): number {
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= n; ++i) {
        for (let j = 0; i + j <= n; ++j) {
            const k = n - i - j;
            let p = mul(i / n, t.v[0]);
            p = add(p, mul(j / n, t.v[1]));
            p = add(p, mul(k / n, t.v[2]));
            best = Math.min(best, pointRectangleSqr(p, r));
        }
    }
    return Math.sqrt(best);
}

function expectConsistent(t: Triangle, r: Rectangle,
    result: DistTriangle3Rectangle3Result): void {
    const delta = sub(result.closest[0], result.closest[1]);
    expect(Math.sqrt(dot(delta, delta))).toBeCloseTo(result.distance, 8);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 8);

    // The barycentric coordinates reproduce closest[0].
    const b = result.barycentric;
    expect(b[0] + b[1] + b[2]).toBeCloseTo(1, 8);
    let reconstructed = new Vector(3);
    for (let i = 0; i < 3; ++i) {
        expect(b[i]).toBeGreaterThanOrEqual(-1e-8);
        reconstructed = add(reconstructed, mul(b[i], t.v[i]));
    }
    for (let i = 0; i < 3; ++i) {
        expect(reconstructed.values[i])
            .toBeCloseTo(result.closest[0].values[i], 7);
    }

    // The cartesian coordinates are within the extents and reproduce
    // closest[1].
    const c = add(r.center, add(mul(result.cartesian[0], r.axis[0]),
        mul(result.cartesian[1], r.axis[1])));
    for (let i = 0; i < 2; ++i) {
        expect(Math.abs(result.cartesian[i]))
            .toBeLessThanOrEqual(r.extent.values[i] + 1e-8);
    }
    for (let i = 0; i < 3; ++i) {
        expect(c.values[i]).toBeCloseTo(result.closest[1].values[i], 7);
    }
}

describe('DistTriangle3Rectangle3', () => {
    const query = new DistTriangle3Rectangle3();

    it('computes the distance between parallel overlapping objects', () => {
        const t = tri([-1, -1, 0], [1, -1, 0], [0, 1, 0]);
        const r = rect([0, 0, 4], frame(0, 0), 2, 2);
        const result = query.compute(t, r);
        expect(result.distance).toBeCloseTo(4, 10);
        expect(result.closest[0].values[2]).toBeCloseTo(0, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(4, 10);
        expectConsistent(t, r, result);
    });

    it('reports zero distance when the objects intersect', () => {
        const t = tri([-2, 0, -2], [2, 0, -2], [0, 0, 3]);
        const r = rect([0, 0, 0], frame(0, 0), 2, 2);
        const result = query.compute(t, r);
        expect(result.distance).toBeCloseTo(0, 10);
        expectConsistent(t, r, result);
    });

    it('computes the distance for coplanar separated objects', () => {
        const t = tri([0, 0, 0], [1, 0, 0], [0, 1, 0]);
        const r = rect([5, 0, 0], frame(0, 0), 1, 1);
        const result = query.compute(t, r);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.cartesian[0]).toBeCloseTo(-1, 8);
        expectConsistent(t, r, result);
    });

    it('finds an edge-vertex closest pair', () => {
        // A triangle far along +x with a vertex nearest the rectangle
        // corner (1,1,0).
        const t = tri([5, 5, 0], [6, 5, 0], [5, 6, 0]);
        const r = rect([0, 0, 0], frame(0, 0), 1, 1);
        const result = query.compute(t, r);
        expect(result.distance).toBeCloseTo(Math.sqrt(16 + 16), 8);
        expect(result.barycentric[0]).toBeCloseTo(1, 8);
        expect(result.cartesian[0]).toBeCloseTo(1, 8);
        expect(result.cartesian[1]).toBeCloseTo(1, 8);
        expectConsistent(t, r, result);
    });

    it('matches a dense sampling on random inputs', () => {
        let seed = 864297531;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rv = () => [rand() * 6 - 3, rand() * 6 - 3, rand() * 6 - 3];

        for (let trial = 0; trial < 20; ++trial) {
            const t = tri(rv(), rv(), rv());
            const r = rect(rv(), frame(rand() * Math.PI, rand() * Math.PI),
                rand() * 2 + 0.1, rand() * 2 + 0.1);
            const result = query.compute(t, r);
            const sampled = sampledDistance(t, r, 64);
            expect(result.distance).toBeLessThanOrEqual(sampled + 1e-8);
            expect(sampled - result.distance).toBeLessThan(0.25);
            expectConsistent(t, r, result);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistTriangle3Rectangle3.h.
// ---------------------------------------------------------------------------

describe('DistTriangle3Rectangle3 verification', () => {
    const query = new DistTriangle3Rectangle3();
    const ttQuery = new DistTriangle3Triangle3();

    const rectArb = fc.tuple(wellScaledVector(3, -5, 5), rotationFrame(3),
        fc.array(finite(0.05, 3), { minLength: 2, maxLength: 2 }))
        .map(([c, R, e]) =>
            Rectangle.fromCenterAxisExtent(c, [R[0], R[1]], v(e[0], e[1])));

    const triArb = fc.tuple(latticeVector(3, -6, 6), latticeVector(3, -6, 6),
        latticeVector(3, -6, 6))
        .filter(([a, b, c]) => {
            const n = cross(sub(b, a), sub(c, a));
            return dot(n, n) > 4;
        })
        .map(([a, b, c]) => Triangle.fromVertices(a, b, c));

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(triArb, rectArb), ([t, r]) => {
            const res = query.compute(t, r);
            // The absolute tolerance is 1e-6: these queries accumulate the
            // squared distance while clamping to faces and edges, so a
            // near-touching configuration loses about half the mantissa and
            // the distance carries an absolute error of order sqrt(eps)
            // times the coordinate scale. A translation or frame error
            // would show up as an O(1) discrepancy.
            expectClose(res.sqrDistance, res.distance * res.distance,
                1e-12, 1e-12);
            const d = sub(res.closest[0], res.closest[1]);
            expectClose(Math.sqrt(dot(d, d)), res.distance, 1e-6, 1e-8);

            const b = res.barycentric;
            expectClose(b[0] + b[1] + b[2], 1, 1e-9, 1e-9);
            let rebuilt = new Vector(3);
            for (let i = 0; i < 3; ++i) {
                expect(b[i]).toBeGreaterThanOrEqual(-1e-8);
                rebuilt = add(rebuilt, mul(b[i], t.v[i]));
            }
            expectVectorClose(rebuilt, res.closest[0], 1e-7, 1e-7);

            let onRect = r.center.clone();
            for (let i = 0; i < 2; ++i) {
                expect(Math.abs(res.cartesian[i]))
                    .toBeLessThanOrEqual(r.extent.values[i] + 1e-8);
                onRect = add(onRect, mul(res.cartesian[i], r.axis[i]));
            }
            expectVectorClose(onRect, res.closest[1], 1e-7, 1e-7);
        });
    });

    // The rectangle is the union of the triangles <0,1,3> and <0,3,2> of its
    // bit-ordered vertices, so the distance equals the minimum over the two
    // triangle-triangle queries.
    it('equals the minimum over the two triangles of the rectangle', () => {
        check(fc.tuple(triArb, rectArb), ([t, r]) => {
            const p = r.getVertices();
            const best = Math.min(
                ttQuery.compute(t,
                    Triangle.fromVertices(p[0], p[1], p[3])).sqrDistance,
                ttQuery.compute(t,
                    Triangle.fromVertices(p[0], p[3], p[2])).sqrDistance);
            const res = query.compute(t, r);
            expectClose(res.sqrDistance, best, 1e-8, 1e-8);
        });
    });

    it('brackets a barycentric sampling of the triangle', () => {
        const rng = seededRandom(0x3141c0de);
        const n = 40;
        for (let k = 0; k < 25; ++k) {
            const p = () => [6 * rng() - 3, 6 * rng() - 3, 6 * rng() - 3];
            const t = tri(p(), p(), p());
            const r = rect(p(), frame(2 * Math.PI * rng(),
                2 * Math.PI * rng()), 0.5 + rng(), 0.5 + rng());
            const res = query.compute(t, r);
            const s = sampledDistance(t, r, n);
            expect(res.distance).toBeLessThanOrEqual(s + 1e-9);
            let maxEdge = 0;
            for (let i = 0; i < 3; ++i) {
                const e = sub(t.v[(i + 1) % 3], t.v[i]);
                maxEdge = Math.max(maxEdge, Math.sqrt(dot(e, e)));
            }
            expect(s - res.distance)
                .toBeLessThanOrEqual(maxEdge / n + 1e-9);
        }
    }, 30000);

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(triArb, rectArb, rotationFrame(3),
            wellScaledVector(3, -4, 4)), ([t, r, R, tr]) => {
                const rot = (q: Vector) => add(add(
                    mul(q.values[0], R[0]), mul(q.values[1], R[1])),
                    mul(q.values[2], R[2]));
                const xf = (q: Vector) => add(rot(q), tr);
                const a = query.compute(t, r);
                const b = query.compute(
                    Triangle.fromVertices(xf(t.v[0]), xf(t.v[1]), xf(t.v[2])),
                    Rectangle.fromCenterAxisExtent(xf(r.center),
                        [rot(r.axis[0]), rot(r.axis[1])], r.extent));
                expectClose(a.distance, b.distance, 1e-8, 1e-8);
            });
    });

    it('handles a zero-extent (point) rectangle', () => {
        check(fc.tuple(triArb, wellScaledVector(3, -5, 5), rotationFrame(3)),
            ([t, c, R]) => {
                const pt = Rectangle.fromCenterAxisExtent(c, [R[0], R[1]],
                    v(0, 0));
                const res = query.compute(t, pt);
                expectVectorClose(res.closest[1], c, 1e-8, 1e-8);
                // '+ 0' normalizes -0, which toBe compares with Object.is.
                expect(res.cartesian[0] + 0).toBe(0);
                expect(res.cartesian[1] + 0).toBe(0);
            });
    });
});
