import { describe, expect, it } from 'vitest';
import { DistTriangle3Triangle3 } from '../src/DistTriangle3Triangle3';
import type { DistTriangle3Triangle3Result }
    from '../src/DistTriangle3Triangle3';
import { Triangle } from '../src/Triangle';
import { Vector, add, dot, mul, sub } from '../src/Vector';
import { cross } from '../src/Vector3';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function tri(v0: number[], v1: number[], v2: number[]): Triangle {
    return Triangle.fromVertices(v(...v0), v(...v1), v(...v2));
}

// An independent point-segment squared distance.
function pointSegmentSqr(p: Vector, a: Vector, b: Vector): number {
    const d = sub(b, a);
    const dd = dot(d, d);
    let t = dd > 0 ? dot(sub(p, a), d) / dd : 0;
    t = Math.min(Math.max(t, 0), 1);
    const delta = sub(p, add(a, mul(t, d)));
    return dot(delta, delta);
}

// An independent point-triangle squared distance: project onto the triangle
// plane, and if the projection is inside use the plane distance, otherwise
// take the minimum over the three edges.
function pointTriangleSqr(p: Vector, t: Triangle): number {
    const e10 = sub(t.v[1], t.v[0]);
    const e20 = sub(t.v[2], t.v[0]);
    const n = cross(e10, e20);
    const nn = dot(n, n);
    if (nn > 0) {
        const delta = sub(p, t.v[0]);
        const h = dot(n, delta) / nn;
        const proj = sub(delta, mul(h, n));
        const nxp = cross(n, proj);
        const b1 = dot(e20, nxp) / nn;
        const b2 = -dot(e10, nxp) / nn;
        const b0 = 1 - b1 - b2;
        if (b0 >= 0 && b1 >= 0 && b2 >= 0) {
            return h * h * nn;
        }
    }
    return Math.min(
        pointSegmentSqr(p, t.v[0], t.v[1]),
        pointSegmentSqr(p, t.v[1], t.v[2]),
        pointSegmentSqr(p, t.v[2], t.v[0]));
}

// A barycentric grid sampling of triangle0 combined with the exact
// point-triangle distance to triangle1 gives an upper bound whose error is
// at most the grid spacing.
function sampledDistance(t0: Triangle, t1: Triangle, n: number): number {
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= n; ++i) {
        for (let j = 0; i + j <= n; ++j) {
            const k = n - i - j;
            let p = mul(i / n, t0.v[0]);
            p = add(p, mul(j / n, t0.v[1]));
            p = add(p, mul(k / n, t0.v[2]));
            best = Math.min(best, pointTriangleSqr(p, t1));
        }
    }
    return Math.sqrt(best);
}

function expectConsistent(t0: Triangle, t1: Triangle,
    result: DistTriangle3Triangle3Result): void {
    const delta = sub(result.closest[0], result.closest[1]);
    expect(Math.sqrt(dot(delta, delta))).toBeCloseTo(result.distance, 8);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 8);

    const check = (b: [number, number, number], t: Triangle, c: Vector) => {
        expect(b[0] + b[1] + b[2]).toBeCloseTo(1, 8);
        let reconstructed = new Vector(3);
        for (let i = 0; i < 3; ++i) {
            expect(b[i]).toBeGreaterThanOrEqual(-1e-8);
            reconstructed = add(reconstructed, mul(b[i], t.v[i]));
        }
        for (let i = 0; i < 3; ++i) {
            expect(reconstructed.values[i]).toBeCloseTo(c.values[i], 7);
        }
    };
    check(result.barycentric0, t0, result.closest[0]);
    check(result.barycentric1, t1, result.closest[1]);
}

describe('DistTriangle3Triangle3', () => {
    const query = new DistTriangle3Triangle3();

    it('computes the distance between parallel overlapping triangles', () => {
        const t0 = tri([0, 0, 0], [4, 0, 0], [0, 4, 0]);
        const t1 = tri([0, 0, 3], [4, 0, 3], [0, 4, 3]);
        const result = query.compute(t0, t1);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[0].values[2]).toBeCloseTo(0, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(3, 10);
        expectConsistent(t0, t1, result);
    });

    it('reports zero distance for intersecting triangles', () => {
        const t0 = tri([-2, 0, 0], [2, 0, 0], [0, 3, 0]);
        const t1 = tri([0, 1, -2], [0, 1, 2], [0, 4, 0]);
        const result = query.compute(t0, t1);
        expect(result.distance).toBeCloseTo(0, 10);
        expectConsistent(t0, t1, result);
    });

    it('computes the distance between coplanar separated triangles', () => {
        // Both triangles lie in z = 0 and are separated along x by 3.
        const t0 = tri([0, 0, 0], [1, 0, 0], [0, 1, 0]);
        const t1 = tri([4, 0, 0], [5, 0, 0], [4, 1, 0]);
        const result = query.compute(t0, t1);
        expect(result.distance).toBeCloseTo(3, 10);
        expectConsistent(t0, t1, result);
    });

    it('is symmetric in its arguments', () => {
        const t0 = tri([0.3, -1.2, 0.4], [2.5, 0.7, -0.8], [1.1, 2.2, 1.3]);
        const t1 = tri([-3.1, 1.4, 2.9], [-1.8, 2.6, 4.1], [-2.5, 0.2, 3.3]);
        const r01 = query.compute(t0, t1);
        const r10 = query.compute(t1, t0);
        expect(r10.distance).toBeCloseTo(r01.distance, 9);
        expectConsistent(t0, t1, r01);
        expectConsistent(t1, t0, r10);
    });

    it('handles a degenerate (segment-like) triangle', () => {
        const t0 = tri([0, 0, 0], [2, 0, 0], [1, 0, 0]);
        const t1 = tri([0, 5, 0], [2, 5, 0], [1, 6, 0]);
        const result = query.compute(t0, t1);
        expect(result.distance).toBeCloseTo(5, 8);
    });

    it('matches a dense sampling on random inputs', () => {
        let seed = 192837465;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rv = () => [rand() * 6 - 3, rand() * 6 - 3, rand() * 6 - 3];

        for (let trial = 0; trial < 20; ++trial) {
            const t0 = tri(rv(), rv(), rv());
            const t1 = tri(rv(), rv(), rv());
            const result = query.compute(t0, t1);
            const sampled = Math.min(sampledDistance(t0, t1, 48),
                sampledDistance(t1, t0, 48));
            expect(result.distance).toBeLessThanOrEqual(sampled + 1e-8);
            expect(sampled - result.distance).toBeLessThan(0.3);
            expectConsistent(t0, t1, result);
        }
    });
});
