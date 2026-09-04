import { describe, expect, it } from 'vitest';
import { DistPoint3Tetrahedron3 } from '../src/DistPoint3Tetrahedron3.js';
import type { DistPoint3Tetrahedron3Result }
    from '../src/DistPoint3Tetrahedron3.js';
import { Tetrahedron3 } from '../src/Tetrahedron3.js';
import { DistPointTriangle } from '../src/DistPointTriangle.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { computeBarycentrics3, cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, latticeVector,
    rotationFrame, seededRandom, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistPoint3Tetrahedron3.h.
// ---------------------------------------------------------------------------

describe('DistPoint3Tetrahedron3 verification', () => {
    const query = new DistPoint3Tetrahedron3();
    const ptQuery = new DistPointTriangle();

    // Lattice vertices keep the tetrahedron well conditioned; the volume
    // filter rejects the near-degenerate draws whose face normals (built by
    // UnitCross in Tetrahedron3.getPlanes) would lose all precision.
    const tetraArb = fc.tuple(latticeVector(3, -5, 5), latticeVector(3, -5, 5),
        latticeVector(3, -5, 5), latticeVector(3, -5, 5))
        .filter(([a, b, c, d]) =>
            Math.abs(dot(sub(b, a), cross(sub(c, a), sub(d, a)))) > 2)
        .map(([a, b, c, d]) => Tetrahedron3.fromVertices(a, b, c, d));

    const pointArb = wellScaledVector(3, -8, 8);

    // The exact squared distance to the solid tetrahedron, computed
    // independently: zero when the point is inside, otherwise the minimum
    // over the four triangular faces.
    function referenceSqrDistance(p: Vector, t: Tetrahedron3): number {
        const bary = computeBarycentrics3(p, t.v[0], t.v[1], t.v[2], t.v[3]);
        if (bary.valid && bary.bary.every(b => b >= 0)) {
            return 0;
        }
        let best = Number.MAX_VALUE;
        for (let f = 0; f < 4; ++f) {
            const idx = Tetrahedron3.getFaceIndices(f);
            const tri = Triangle.fromVertices(t.v[idx[0]], t.v[idx[1]],
                t.v[idx[2]]);
            best = Math.min(best, ptQuery.compute(p, tri).sqrDistance);
        }
        return best;
    }

    it('reports consistent distances and barycentric closest points', () => {
        check(fc.tuple(pointArb, tetraArb), ([p, t]) => {
            const r = query.compute(p, t);
            expect(r.closest[0].equals(p)).toBe(true);
            expectClose(r.sqrDistance, r.distance * r.distance, 1e-12, 1e-12);
            const d = sub(r.closest[0], r.closest[1]);
            expectClose(Math.sqrt(dot(d, d)), r.distance, 1e-9, 1e-9);

            let sum = 0;
            let rebuilt = new Vector(3);
            for (let i = 0; i < 4; ++i) {
                expect(r.barycentric[i]).toBeGreaterThanOrEqual(-1e-9);
                sum += r.barycentric[i];
                rebuilt = add(rebuilt, mul(r.barycentric[i], t.v[i]));
            }
            expectClose(sum, 1, 1e-9, 1e-9);
            expectVectorClose(rebuilt, r.closest[1], 1e-8, 1e-8);
        });
    });

    it('matches the face-minimum reference distance', () => {
        check(fc.tuple(pointArb, tetraArb), ([p, t]) => {
            const r = query.compute(p, t);
            expectClose(r.sqrDistance, referenceSqrDistance(p, t),
                1e-9, 1e-9);
        });
    });

    it('reports zero for points inside the tetrahedron', () => {
        check(fc.tuple(tetraArb, fc.array(finite(0.05, 1),
            { minLength: 4, maxLength: 4 })), ([t, w]) => {
                const s = w[0] + w[1] + w[2] + w[3];
                let p = new Vector(3);
                for (let i = 0; i < 4; ++i) {
                    p = add(p, mul(w[i] / s, t.v[i]));
                }
                const r = query.compute(p, t);
                expect(r.distance).toBe(0);
                expect(r.sqrDistance).toBe(0);
                expectVectorClose(r.closest[1], p, 1e-12, 1e-12);
            });
    });

    it('never exceeds a barycentric sampling of the solid', () => {
        // Deterministic sampling loop; the timeout covers the cost under a
        // loaded machine.
        const rng = seededRandom(0x5eed1234);
        const t = tetra([0, 0, 0], [3, 0, 0], [0, 4, 0], [0, 0, 5]);
        for (let k = 0; k < 60; ++k) {
            const p = v(12 * rng() - 6, 12 * rng() - 6, 12 * rng() - 6);
            const r = query.compute(p, t);
            expect(r.sqrDistance)
                .toBeLessThanOrEqual(sampledSqrDistance(p, t, 12) + 1e-9);
        }
    }, 30000);

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(pointArb, tetraArb, rotationFrame(3),
            wellScaledVector(3, -5, 5)), ([p, t, R, tr]) => {
                const xf = (q: Vector) => add(add(add(
                    mul(q.values[0], R[0]), mul(q.values[1], R[1])),
                    mul(q.values[2], R[2])), tr);
                const r0 = query.compute(p, t);
                const r1 = query.compute(xf(p), Tetrahedron3.fromVertices(
                    xf(t.v[0]), xf(t.v[1]), xf(t.v[2]), xf(t.v[3])));
                expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
                expectVectorClose(xf(r0.closest[1]), r1.closest[1],
                    1e-7, 1e-7);
            });
    });

    it('is unaffected by the winding of the input vertices', () => {
        // Tetrahedron3.getPlanes flips inner-pointing normals, so swapping
        // two vertices must not change the result.
        check(fc.tuple(pointArb, tetraArb), ([p, t]) => {
            const r0 = query.compute(p, t);
            const r1 = query.compute(p, Tetrahedron3.fromVertices(
                t.v[1], t.v[0], t.v[2], t.v[3]));
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expectVectorClose(r0.closest[1], r1.closest[1], 1e-8, 1e-8);
        });
    });
});
