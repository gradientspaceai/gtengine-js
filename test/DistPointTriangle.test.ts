import { describe, expect, it } from 'vitest';
import { DistPointTriangle } from '../src/DistPointTriangle.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function tri(v0: number[], v1: number[], v2: number[]): Triangle {
    return Triangle.fromVertices(v(...v0), v(...v1), v(...v2));
}

// The minimum of |point - X|^2 over a dense sampling of the triangle. This
// is an upper bound for the true squared distance and converges to it as the
// sampling is refined.
function sampledSqrDistance(point: Vector, t: Triangle,
    samples: number): number {
    const edge0 = sub(t.v[1], t.v[0]);
    const edge1 = sub(t.v[2], t.v[0]);
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= samples; ++i) {
        const s = i / samples;
        for (let j = 0; j + i <= samples; ++j) {
            const u = j / samples;
            const x = add(t.v[0], add(mul(s, edge0), mul(u, edge1)));
            const d = sub(point, x);
            best = Math.min(best, dot(d, d));
        }
    }
    return best;
}

describe('DistPointTriangle', () => {
    const query = new DistPointTriangle();
    // The reference triangle used by most of the region tests.
    const unit = tri([0, 0, 0], [1, 0, 0], [0, 1, 0]);

    // Check that the reported closest point really is the barycentric
    // combination of the vertices, that the barycentric coordinates are a
    // valid convex combination, and that the reported distance matches.
    function verify(result: ReturnType<DistPointTriangle['compute']>,
        point: Vector, t: Triangle): void {
        const b = result.barycentric;
        expect(b[0] + b[1] + b[2]).toBeCloseTo(1, 12);
        for (let i = 0; i < 3; ++i) {
            expect(b[i]).toBeGreaterThanOrEqual(-1e-14);
            expect(b[i]).toBeLessThanOrEqual(1 + 1e-14);
        }
        const x = add(add(mul(b[0], t.v[0]), mul(b[1], t.v[1])),
            mul(b[2], t.v[2]));
        for (let i = 0; i < x.size; ++i) {
            expect(result.closest[1].values[i]).toBeCloseTo(x.values[i], 10);
        }
        expect(result.closest[0].values).toEqual(point.values);
        const diff = sub(result.closest[0], result.closest[1]);
        expect(result.sqrDistance).toBeCloseTo(dot(diff, diff), 12);
        expect(result.distance).toBeCloseTo(Math.sqrt(result.sqrDistance), 12);
    }

    it('handles a point above an interior point (region 0)', () => {
        const point = v(0.25, 0.25, 3);
        const result = query.compute(point, unit);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.barycentric[0]).toBeCloseTo(0.5, 12);
        expect(result.barycentric[1]).toBeCloseTo(0.25, 12);
        expect(result.barycentric[2]).toBeCloseTo(0.25, 12);
        expect(result.closest[1].values).toEqual([0.25, 0.25, 0]);
        verify(result, point, unit);
    });

    it('reports zero distance for a point inside the triangle', () => {
        const point = v(0.2, 0.3, 0);
        const result = query.compute(point, unit);
        expect(result.distance).toBeCloseTo(0, 14);
        expect(result.sqrDistance).toBeCloseTo(0, 14);
        verify(result, point, unit);
    });

    it('clamps to the vertex V0 (region 4)', () => {
        const point = v(-1, -1, 0);
        const result = query.compute(point, unit);
        expect(result.closest[1].values).toEqual([0, 0, 0]);
        expect(result.distance).toBeCloseTo(Math.SQRT2, 12);
        expect(result.barycentric).toEqual([1, 0, 0]);
        verify(result, point, unit);
    });

    it('clamps to the edge <V0,V2> (region 3)', () => {
        const point = v(-1, 0.5, 0);
        const result = query.compute(point, unit);
        expect(result.closest[1].values).toEqual([0, 0.5, 0]);
        expect(result.distance).toBeCloseTo(1, 12);
        expect(result.barycentric[1]).toBe(0);
        verify(result, point, unit);
    });

    it('clamps to the edge <V0,V1> (region 5)', () => {
        const point = v(0.5, -1, 0);
        const result = query.compute(point, unit);
        expect(result.closest[1].values).toEqual([0.5, 0, 0]);
        expect(result.distance).toBeCloseTo(1, 12);
        expect(result.barycentric[2]).toBe(0);
        verify(result, point, unit);
    });

    it('clamps to the edge <V1,V2> (region 1)', () => {
        const point = v(1, 1, 0);
        const result = query.compute(point, unit);
        expect(result.closest[1].values[0]).toBeCloseTo(0.5, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0.5, 12);
        expect(result.distance).toBeCloseTo(Math.sqrt(0.5), 12);
        expect(result.barycentric[0]).toBeCloseTo(0, 12);
        verify(result, point, unit);
    });

    it('clamps to the vertex V2 (region 2)', () => {
        const point = v(-0.5, 2, 0);
        const result = query.compute(point, unit);
        expect(result.closest[1].values).toEqual([0, 1, 0]);
        expect(result.distance).toBeCloseTo(Math.sqrt(0.25 + 1), 12);
        expect(result.barycentric).toEqual([0, 0, 1]);
        verify(result, point, unit);
    });

    it('clamps to the vertex V1 (region 6)', () => {
        const point = v(2, -0.5, 0);
        const result = query.compute(point, unit);
        expect(result.closest[1].values).toEqual([1, 0, 0]);
        expect(result.distance).toBeCloseTo(Math.sqrt(1 + 0.25), 12);
        expect(result.barycentric).toEqual([0, 1, 0]);
        verify(result, point, unit);
    });

    it('works in 2D', () => {
        const t = tri([0, 0], [4, 0], [0, 4]);
        const point = v(4, 4);
        const result = query.compute(point, t);
        // The closest point is the midpoint of the hypotenuse, (2,2).
        expect(result.closest[1].values[0]).toBeCloseTo(2, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(2, 12);
        expect(result.distance).toBeCloseTo(Math.sqrt(8), 12);
        verify(result, point, t);
    });

    it('works in 4D', () => {
        const t = tri([0, 0, 0, 0], [2, 0, 0, 0], [0, 2, 0, 0]);
        const point = v(0.5, 0.5, 3, 4);
        const result = query.compute(point, t);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values).toEqual([0.5, 0.5, 0, 0]);
        verify(result, point, t);
    });

    it('agrees with useConjugateGradient on the region tests', () => {
        const points = [
            v(0.25, 0.25, 3), v(0.2, 0.3, 0), v(-1, -1, 0), v(-1, 0.5, 0),
            v(0.5, -1, 0), v(1, 1, 0), v(-0.5, 2, 0), v(2, -0.5, 0)
        ];
        for (const point of points) {
            const r0 = query.compute(point, unit);
            const r1 = query.useConjugateGradient(point, unit);
            expect(r1.distance).toBeCloseTo(r0.distance, 10);
            for (let i = 0; i < 3; ++i) {
                expect(r1.barycentric[i]).toBeCloseTo(r0.barycentric[i], 10);
            }
            verify(r1, point, unit);
        }
    });

    it('useConjugateGradient handles a degenerate triangle', () => {
        // All three vertices are the same point, so the closest triangle
        // point is that vertex. (The compute(...) path divides by a zero
        // determinant for degenerate triangles; upstream documents
        // UseConjugateGradient as the robust alternative.)
        const t = tri([1, 2, 3], [1, 2, 3], [1, 2, 3]);
        const result = query.useConjugateGradient(v(1, 2, 7), t);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values).toEqual([1, 2, 3]);
    });

    it('useConjugateGradient handles a sliver triangle', () => {
        // V2 lies on the segment <V0,V1>, so the triangle is a segment.
        const t = tri([0, 0, 0], [4, 0, 0], [2, 0, 0]);
        const result = query.useConjugateGradient(v(1, 3, 0), t);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 10);
    });

    it('agrees with a dense sampling of the triangle', () => {
        let seed = 20250510;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 8 - 4;
        };
        for (let trial = 0; trial < 40; ++trial) {
            const t = tri([rand(), rand(), rand()],
                [rand(), rand(), rand()], [rand(), rand(), rand()]);
            const point = v(rand(), rand(), rand());
            const result = query.compute(point, t);
            verify(result, point, t);

            const best = sampledSqrDistance(point, t, 150);
            // The query is the true minimum, so it cannot exceed any
            // sampled value, and the sampling is fine enough to be close.
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
            expect(Math.sqrt(best) - result.distance).toBeLessThan(0.15);

            const cg = query.useConjugateGradient(point, t);
            expect(cg.distance).toBeCloseTo(result.distance, 8);
        }
    });
});
