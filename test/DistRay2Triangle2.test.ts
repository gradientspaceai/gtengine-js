import { describe, expect, it } from 'vitest';
import { DistLine2Triangle2 } from '../src/DistLine2Triangle2';
import { DistPointTriangle } from '../src/DistPointTriangle';
import { DistRay2Triangle2 } from '../src/DistRay2Triangle2';
import { Line } from '../src/Line';
import { Ray } from '../src/Ray';
import { Triangle } from '../src/Triangle';
import { Vector, add, dot, mul, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistRay2Triangle2', () => {
    const query = new DistRay2Triangle2();
    const tri = Triangle.fromVertices(v(0, 0), v(1, 0), v(0, 1));

    it('reports zero distance for a ray that hits the triangle', () => {
        const result = query.compute(ray([-5, 0.25], [1, 0]), tri);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter).toBeGreaterThanOrEqual(0);
    });

    it('reports zero distance for a ray starting inside the triangle', () => {
        const result = query.compute(ray([0.25, 0.25], [1, 1]), tri);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('uses the ray origin when the ray points away from the triangle',
        () => {
            const result = query.compute(ray([5, 0], [1, 0]), tri);
            expect(result.parameter).toBe(0);
            expect(result.distance).toBeCloseTo(4, 12);
            expect(result.closest[0].values).toEqual([5, 0]);
            expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
            expect(result.closest[1].values[1]).toBeCloseTo(0, 12);
        });

    it('matches the point-triangle query when the ray points away', () => {
        const ptQuery = new DistPointTriangle();
        const origin = v(4, -3);
        const result = query.compute(
            Ray.fromOriginDirection(origin, v(1, -1)), tri);
        const ptResult = ptQuery.compute(origin, tri);
        expect(result.distance).toBeCloseTo(ptResult.distance, 12);
        for (let i = 0; i < 3; ++i) {
            expect(result.barycentric[i]).toBeCloseTo(ptResult.barycentric[i],
                12);
        }
    });

    it('matches the line query when the closest line point is on the ray',
        () => {
            const ltQuery = new DistLine2Triangle2();
            const r = ray([-5, 3], [1, 0]);
            const result = query.compute(r, tri);
            const ltResult = ltQuery.compute(
                Line.fromOriginDirection(r.origin, r.direction), tri);
            expect(ltResult.parameter).toBeGreaterThanOrEqual(0);
            expect(result.distance).toBeCloseTo(ltResult.distance, 12);
        });

    it('handles a sliver triangle', () => {
        // A triangle with a tiny but nonzero area; the closest point is on
        // the long edge.
        const sliver = Triangle.fromVertices(v(0, 0), v(2, 0), v(1, 1e-9));
        const result = query.compute(ray([1, 3], [0, 1]), sliver);
        expect(result.parameter).toBe(0);
        expect(result.distance).toBeCloseTo(3, 8);
    });

    it('handles a degenerate triangle whose edges are traversed', () => {
        // A zero-area triangle collapsed onto a segment. The closest line
        // point has a nonnegative parameter, so the edge-based line query
        // handles it.
        const degenerate = Triangle.fromVertices(v(0, 0), v(2, 0), v(1, 0));
        const result = query.compute(ray([-5, 3], [1, 0]), degenerate);
        expect(result.parameter).toBeGreaterThanOrEqual(0);
        expect(result.distance).toBeCloseTo(3, 10);
    });

    it('agrees with a dense sampling of the ray and triangle', () => {
        const rnd = makeRandom(60013);
        const t = Triangle.fromVertices(v(0.5, -1), v(2, 0.5), v(-1, 1.5));

        for (let trial = 0; trial < 50; ++trial) {
            const origin = v(6 * rnd() - 3, 6 * rnd() - 3);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const r = Ray.fromOriginDirection(origin, dir);
            const result = query.compute(r, t);

            expect(result.parameter).toBeGreaterThanOrEqual(0);
            const onRay = add(r.origin, mul(result.parameter, r.direction));
            for (let i = 0; i < 2; ++i) {
                expect(onRay.values[i]).toBeCloseTo(
                    result.closest[0].values[i], 7);
            }

            const b = result.barycentric;
            expect(b[0] + b[1] + b[2]).toBeCloseTo(1, 8);
            for (const bi of b) {
                expect(bi).toBeGreaterThanOrEqual(-1e-8);
            }

            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 7);

            // No sampled (ray point, triangle point) pair is closer.
            const n = 50;
            const dd = dot(r.direction, r.direction);
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= n; ++i) {
                for (let j = 0; i + j <= n; ++j) {
                    const b1 = i / n, b2 = j / n, b0 = 1 - b1 - b2;
                    const p = add(mul(b0, t.v[0]),
                        add(mul(b1, t.v[1]), mul(b2, t.v[2])));
                    const w = sub(p, r.origin);
                    let s = dot(w, r.direction) / dd;
                    if (s < 0) {
                        s = 0;
                    }
                    const f = sub(w, mul(s, r.direction));
                    best = Math.min(best, dot(f, f));
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});
