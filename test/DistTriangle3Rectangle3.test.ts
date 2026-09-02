import { describe, expect, it } from 'vitest';
import { DistTriangle3Rectangle3 } from '../src/DistTriangle3Rectangle3';
import type { DistTriangle3Rectangle3Result }
    from '../src/DistTriangle3Rectangle3';
import { Rectangle } from '../src/Rectangle';
import { Triangle } from '../src/Triangle';
import { Vector, add, dot, mul, sub } from '../src/Vector';

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
