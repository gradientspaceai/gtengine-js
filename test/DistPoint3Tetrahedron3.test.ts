import { describe, expect, it } from 'vitest';
import { DistPoint3Tetrahedron3 } from '../src/DistPoint3Tetrahedron3';
import type { DistPoint3Tetrahedron3Result }
    from '../src/DistPoint3Tetrahedron3';
import { Tetrahedron3 } from '../src/Tetrahedron3';
import { Vector, add, dot, mul, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function tetra(v0: number[], v1: number[], v2: number[], v3: number[]):
    Tetrahedron3 {
    return Tetrahedron3.fromVertices(v(...v0), v(...v1), v(...v2), v(...v3));
}

const canonical = tetra([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]);

// A dense barycentric sampling of the solid tetrahedron gives an upper bound
// for the distance that converges to the true value as n grows.
function sampledSqrDistance(p: Vector, t: Tetrahedron3, n: number): number {
    let best = Number.MAX_VALUE;
    for (let i0 = 0; i0 <= n; ++i0) {
        for (let i1 = 0; i0 + i1 <= n; ++i1) {
            for (let i2 = 0; i0 + i1 + i2 <= n; ++i2) {
                const i3 = n - i0 - i1 - i2;
                let q = mul(i0 / n, t.v[0]);
                q = add(q, mul(i1 / n, t.v[1]));
                q = add(q, mul(i2 / n, t.v[2]));
                q = add(q, mul(i3 / n, t.v[3]));
                const d = sub(p, q);
                best = Math.min(best, dot(d, d));
            }
        }
    }
    return best;
}

function expectConsistent(p: Vector, t: Tetrahedron3,
    result: DistPoint3Tetrahedron3Result): void {
    // closest[0] is the query point.
    expect(result.closest[0].equals(p)).toBe(true);

    // The distance is the length of the segment joining the closest points.
    const delta = sub(result.closest[0], result.closest[1]);
    expect(Math.sqrt(dot(delta, delta))).toBeCloseTo(result.distance, 9);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 9);

    // The barycentric coordinates are a partition of unity and reproduce the
    // closest point.
    const sum = result.barycentric[0] + result.barycentric[1]
        + result.barycentric[2] + result.barycentric[3];
    expect(sum).toBeCloseTo(1, 9);
    let reconstructed = new Vector(3);
    for (let i = 0; i < 4; ++i) {
        expect(result.barycentric[i]).toBeGreaterThanOrEqual(-1e-9);
        reconstructed = add(reconstructed, mul(result.barycentric[i], t.v[i]));
    }
    for (let i = 0; i < 3; ++i) {
        expect(reconstructed.values[i])
            .toBeCloseTo(result.closest[1].values[i], 8);
    }
}

describe('DistPoint3Tetrahedron3', () => {
    const query = new DistPoint3Tetrahedron3();

    it('reports zero distance for interior points', () => {
        const points = [v(0.1, 0.1, 0.1), v(0.25, 0.25, 0.25),
            v(0.5, 0.2, 0.2), v(0.01, 0.01, 0.9)];
        for (const p of points) {
            const result = query.compute(p, canonical);
            expect(result.distance).toBe(0);
            expect(result.sqrDistance).toBe(0);
            expect(result.closest[1].equals(p)).toBe(true);
            expectConsistent(p, canonical, result);
        }
    });

    it('reports zero distance for vertices and face centroids', () => {
        for (let i = 0; i < 4; ++i) {
            const result = query.compute(canonical.v[i], canonical);
            expect(result.distance).toBeCloseTo(0, 12);
            expect(result.barycentric[i]).toBeCloseTo(1, 9);
        }

        const centroid = canonical.computeCentroid();
        const result = query.compute(centroid, canonical);
        expect(result.distance).toBe(0);
        for (let i = 0; i < 4; ++i) {
            expect(result.barycentric[i]).toBeCloseTo(0.25, 9);
        }
    });

    it('computes known exterior distances', () => {
        // Beyond the vertex (1,0,0) along +x.
        let result = query.compute(v(3, 0, 0), canonical);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
        expect(result.barycentric[1]).toBeCloseTo(1, 9);

        // Below the z = 0 face.
        result = query.compute(v(0.25, 0.25, -2), canonical);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(0, 12);

        // Outside the slanted face x + y + z = 1; the closest point is the
        // centroid (1/3,1/3,1/3) of that face.
        result = query.compute(v(1, 1, 1), canonical);
        expect(result.distance).toBeCloseTo(2 / Math.sqrt(3), 10);
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i]).toBeCloseTo(1 / 3, 9);
        }
        expect(result.barycentric[0]).toBeCloseTo(0, 9);
    });

    it('is invariant under rigid motion', () => {
        const c = Math.cos(0.7), s = Math.sin(0.7);
        const rotate = (p: Vector) => v(
            c * p.values[0] - s * p.values[1] + 3,
            s * p.values[0] + c * p.values[1] - 1,
            p.values[2] + 5);
        const moved = Tetrahedron3.fromVertices(
            rotate(canonical.v[0]), rotate(canonical.v[1]),
            rotate(canonical.v[2]), rotate(canonical.v[3]));

        const p = v(1.3, -0.4, 0.8);
        const r0 = query.compute(p, canonical);
        const r1 = query.compute(rotate(p), moved);
        expect(r1.distance).toBeCloseTo(r0.distance, 10);
        for (let i = 0; i < 4; ++i) {
            expect(r1.barycentric[i]).toBeCloseTo(r0.barycentric[i], 8);
        }
    });

    it('matches a dense sampling of the solid tetrahedron', () => {
        let seed = 1357911;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const t = tetra([0.2, -0.3, 0.1], [2.1, 0.4, -0.2], [0.3, 1.9, 0.5],
            [-0.1, 0.2, 2.4]);
        for (let trial = 0; trial < 40; ++trial) {
            const p = v(rand() * 6 - 3, rand() * 6 - 3, rand() * 6 - 3);
            const result = query.compute(p, t);
            const sampled = Math.sqrt(sampledSqrDistance(p, t, 12));
            // The query is exact, so it never exceeds a sampled distance,
            // and the sampling is fine enough to bound it from above.
            expect(result.distance).toBeLessThanOrEqual(sampled + 1e-9);
            expect(sampled - result.distance).toBeLessThan(0.25);
            expectConsistent(p, t, result);
        }
    });
});
