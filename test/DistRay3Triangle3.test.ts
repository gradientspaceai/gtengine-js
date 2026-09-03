import { describe, expect, it } from 'vitest';
import { DistPointTriangle } from '../src/DistPointTriangle.js';
import { DistRay3Triangle3 } from '../src/DistRay3Triangle3.js';
import { Ray } from '../src/Ray.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, length, mul, sub } from '../src/Vector.js';

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
