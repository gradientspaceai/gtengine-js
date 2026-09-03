import { describe, expect, it } from 'vitest';
import { DistLine2Triangle2 } from '../src/DistLine2Triangle2.js';
import { Line } from '../src/Line.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function tri(a: number[], b: number[], c: number[]): Triangle {
    return Triangle.fromVertices(v(...a), v(...b), v(...c));
}

// The point of the triangle with the given barycentric coordinates.
function fromBarycentric(t: Triangle, b: readonly number[]): Vector {
    return add(add(mul(b[0], t.v[0]), mul(b[1], t.v[1])), mul(b[2], t.v[2]));
}

// The exact squared distance from a point to the solid triangle, computed by
// dense sampling of the barycentric domain.
function pointTriangleSqrDistance(p: Vector, t: Triangle,
    steps: number): number {
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= steps; ++i) {
        for (let j = 0; i + j <= steps; ++j) {
            const b0 = i / steps;
            const b1 = j / steps;
            const q = fromBarycentric(t, [b0, b1, 1 - b0 - b1]);
            const d = sub(p, q);
            best = Math.min(best, dot(d, d));
        }
    }
    return best;
}

describe('DistLine2Triangle2', () => {
    const query = new DistLine2Triangle2();
    const unitTri = tri([0, 0], [1, 0], [0, 1]);

    it('measures a line strictly on the positive side', () => {
        // The line y = 4 with direction (1,0) has normal Perp(D) = (0,-1),
        // so the triangle is on the negative side; the closest vertex is
        // (0,1) at distance 3.
        const result = query.compute(line([0, 4], [1, 0]), unitTri);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values).toEqual([0, 1]);
        expect(result.barycentric).toEqual([0, 0, 1]);
        expect(result.closest[0].values[1]).toBeCloseTo(4, 10);
    });

    it('measures a line strictly on the other side', () => {
        const result = query.compute(line([0, -2], [1, 0]), unitTri);
        expect(result.distance).toBeCloseTo(2, 10);
        // Two vertices are at y = 0; the first minimum found wins.
        expect(result.closest[1].values[1]).toBe(0);
        expect(result.barycentric[0] + result.barycentric[1]).toBe(1);
    });

    it('reports zero distance for a line crossing the triangle', () => {
        const result = query.compute(line([-1, 0.25], [1, 0]), unitTri);
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
        // The closest points coincide on the triangle boundary.
        expect(result.closest[0].values).toEqual(result.closest[1].values);
    });

    it('reports zero distance when the line contains a vertex', () => {
        // The line through (1,0) with direction (0,1) touches the triangle at
        // that vertex only.
        const result = query.compute(line([1, -5], [0, 1]), unitTri);
        expect(result.distance).toBe(0);
        expect(result.closest[1].values).toEqual([1, 0]);
        expect(result.barycentric).toEqual([0, 1, 0]);
        expect(result.parameter).toBeCloseTo(5, 10);
    });

    it('reports zero distance when the line contains an edge', () => {
        const result = query.compute(line([0, 0], [1, 0]), unitTri);
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
    });

    it('handles a degenerate (single-point) triangle off the line', () => {
        const degenerate = tri([2, 2], [2, 2], [2, 2]);
        const result = query.compute(line([0, 0], [1, 0]), degenerate);
        // All normal components are equal and negative, so the --- branch is
        // taken and the point-to-line distance is reported.
        expect(result.closest[1].values).toEqual([2, 2]);
        expect(result.barycentric).toEqual([1, 0, 0]);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.closest[0].values[0]).toBeCloseTo(2, 12);
        expect(result.closest[0].values[1]).toBeCloseTo(0, 12);
    });

    it('handles a degenerate triangle lying on the line (000 branch)', () => {
        const degenerate = tri([2, 0], [2, 0], [2, 0]);
        const result = query.compute(line([0, 0], [1, 0]), degenerate);
        expect(result.closest[1].values).toEqual([2, 0]);
        expect(result.barycentric).toEqual([1, 0, 0]);
        expect(result.distance).toBe(0);
        expect(result.parameter).toBeCloseTo(2, 12);
    });

    it('places the closest points consistently for separated triangles',
        () => {
            let seed = 1122;
            const rand = () => {
                seed = (seed * 1103515245 + 12345) % 2147483648;
                return seed / 2147483648 * 4 - 2;
            };
            for (let trial = 0; trial < 40; ++trial) {
                // Build a triangle entirely above the line y = 10 so the
                // "no common points" branch is taken for a horizontal line.
                const t = tri([rand(), 11 + Math.abs(rand())],
                    [rand(), 11 + Math.abs(rand())],
                    [rand(), 11 + Math.abs(rand())]);
                const l = line([rand(), 10], [1, 0]);
                const result = query.compute(l, t);

                // The closest triangle point is a vertex.
                const isVertex = t.v.some(vert =>
                    Math.abs(vert.values[0] - result.closest[1].values[0])
                        < 1e-12
                    && Math.abs(vert.values[1] - result.closest[1].values[1])
                        < 1e-12);
                expect(isVertex).toBe(true);

                // The distance equals the minimum vertex height above y = 10.
                const expected = Math.min(...t.v.map(x => x.values[1])) - 10;
                expect(result.distance).toBeCloseTo(expected, 9);

                // The barycentric coordinates reproduce the closest point.
                const q = fromBarycentric(t, result.barycentric);
                expect(q.values[0]).toBeCloseTo(
                    result.closest[1].values[0], 9);
                expect(q.values[1]).toBeCloseTo(
                    result.closest[1].values[1], 9);

                // The line point matches the parameter.
                const onLine = add(l.origin, mul(result.parameter,
                    l.direction));
                expect(result.closest[0].values[0]).toBeCloseTo(
                    onLine.values[0], 9);

                // It agrees with a sampled point-triangle distance.
                const sampled = Math.sqrt(pointTriangleSqrDistance(
                    result.closest[0], t, 60));
                expect(result.distance).toBeLessThanOrEqual(sampled + 1e-9);
            }
        });

    it('always reports a valid barycentric triple', () => {
        let seed = 3344;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 8 - 4;
        };
        for (let trial = 0; trial < 120; ++trial) {
            const t = tri([rand(), rand()], [rand(), rand()],
                [rand(), rand()]);
            const l = line([rand(), rand()], [rand() + 5, rand()]);
            const result = query.compute(l, t);

            const sum = result.barycentric[0] + result.barycentric[1]
                + result.barycentric[2];
            expect(sum).toBeCloseTo(1, 8);
            for (const b of result.barycentric) {
                expect(b).toBeGreaterThanOrEqual(-1e-9);
                expect(b).toBeLessThanOrEqual(1 + 1e-9);
            }

            const q = fromBarycentric(t, result.barycentric);
            expect(q.values[0]).toBeCloseTo(result.closest[1].values[0], 7);
            expect(q.values[1]).toBeCloseTo(result.closest[1].values[1], 7);

            const diff = sub(result.closest[0], result.closest[1]);
            expect(dot(diff, diff)).toBeCloseTo(result.sqrDistance, 8);
        }
    });
});
