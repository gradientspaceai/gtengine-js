import { describe, expect, it } from 'vitest';
import { DistLine2Triangle2 } from '../src/DistLine2Triangle2.js';
import { DistPointTriangle } from '../src/DistPointTriangle.js';
import { DistRay2Triangle2 } from '../src/DistRay2Triangle2.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistRay2Triangle2.h.
// ---------------------------------------------------------------------------

describe('DistRay2Triangle2 verification', () => {
    const query = new DistRay2Triangle2();
    const lineQuery = new DistLine2Triangle2();

    const triArb = fc.tuple(wellScaledVector(2, -5, 5),
        wellScaledVector(2, -5, 5), wellScaledVector(2, -5, 5))
        .filter(([a, b, c]) => Math.abs(
            (b.values[0] - a.values[0]) * (c.values[1] - a.values[1])
            - (b.values[1] - a.values[1]) * (c.values[0] - a.values[0])) > 1)
        .map(([a, b, c]) => Triangle.fromVertices(a, b, c));

    const rayArb = fc.tuple(wellScaledVector(2, -8, 8), unitVector(2))
        .map(([o, d]) => Ray.fromOriginDirection(o, d));

    function pointSegmentDistance(p: Vector, a: Vector, b: Vector): number {
        const d = sub(b, a);
        const dd = dot(d, d);
        const t = dd > 0 ? Math.min(Math.max(dot(sub(p, a), d) / dd, 0), 1) : 0;
        return length(sub(p, add(a, mul(t, d))));
    }

    // Independent distance from a point to a solid triangle in 2D: zero when
    // the point is on the same side of every edge, otherwise the minimum over
    // the three edges.
    function pointTriangleDistance(p: Vector, tri: Triangle): number {
        let positive = 0;
        let negative = 0;
        for (let i = 0; i < 3; ++i) {
            const a = tri.v[i];
            const b = tri.v[(i + 1) % 3];
            const side = (b.values[0] - a.values[0]) * (p.values[1]
                - a.values[1]) - (b.values[1] - a.values[1])
                * (p.values[0] - a.values[0]);
            if (side > 0) { ++positive; }
            else if (side < 0) { ++negative; }
        }
        if (positive === 0 || negative === 0) { return 0; }
        return Math.min(
            pointSegmentDistance(p, tri.v[0], tri.v[1]),
            pointSegmentDistance(p, tri.v[1], tri.v[2]),
            pointSegmentDistance(p, tri.v[2], tri.v[0]));
    }

    // The distance from a point moving along a ray to a convex set is convex
    // in the parameter, so ternary search on [0, tMax] finds the minimum.
    function ternaryMin(f: (t: number) => number, lo: number,
        hi: number): number {
        let a = lo, b = hi;
        for (let i = 0; i < 200; ++i) {
            const m0 = a + (b - a) / 3;
            const m1 = b - (b - a) / 3;
            if (f(m0) <= f(m1)) { b = m1; } else { a = m0; }
        }
        return f(0.5 * (a + b));
    }

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(rayArb, triArb), ([ray, tri]) => {
            const r = query.compute(ray, tri);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            expect(r.parameter).toBeGreaterThanOrEqual(0);
            expectVectorClose(r.closest[0],
                add(ray.origin, mul(r.parameter, ray.direction)), 1e-8, 1e-8);
            expectClose(length(sub(r.closest[0], r.closest[1])), r.distance,
                1e-8, 1e-8);
            const b = r.barycentric;
            expectClose(b[0] + b[1] + b[2], 1, 1e-9, 1e-9);
            for (let i = 0; i < 3; ++i) {
                expect(b[i]).toBeGreaterThanOrEqual(-1e-12);
            }
            expectVectorClose(r.closest[1],
                add(add(mul(b[0], tri.v[0]), mul(b[1], tri.v[1])),
                    mul(b[2], tri.v[2])), 1e-7, 1e-7);
        });
    });

    it('matches an independent convex minimization along the ray', () => {
        check(fc.tuple(rayArb, triArb), ([ray, tri]) => {
            const r = query.compute(ray, tri);
            const best = ternaryMin(t => pointTriangleDistance(
                add(ray.origin, mul(t, ray.direction)), tri), 0, 200);
            expectClose(r.distance, best, 1e-7, 1e-7);
        }, 100);
    });

    it('agrees with the line query when the line parameter is nonnegative',
        () => {
            check(fc.tuple(rayArb, triArb), ([ray, tri]) => {
                const line = Line.fromOriginDirection(ray.origin,
                    ray.direction);
                const rl = lineQuery.compute(line, tri);
                const rr = query.compute(ray, tri);
                if (rl.parameter >= 0) {
                    expectClose(rr.distance, rl.distance, 1e-12, 1e-12);
                    expectClose(rr.parameter, rl.parameter, 1e-12, 1e-12);
                } else {
                    // Otherwise the ray origin is the closest ray point.
                    expect(rr.parameter).toBe(0);
                    expectVectorClose(rr.closest[0], ray.origin, 1e-12,
                        1e-12);
                    expectClose(rr.distance,
                        pointTriangleDistance(ray.origin, tri), 1e-8, 1e-8);
                }
                expect(rr.distance).toBeGreaterThanOrEqual(
                    rl.distance - 1e-12);
            });
        });

    it('reports zero distance for a ray starting inside the triangle', () => {
        check(fc.tuple(triArb, finite(0.05, 0.9), finite(0.05, 0.9),
            unitVector(2)), ([tri, s, t, dir]) => {
            if (s + t > 0.95) { return; }
            const p = add(tri.v[0], add(mul(s, sub(tri.v[1], tri.v[0])),
                mul(t, sub(tri.v[2], tri.v[0]))));
            const r = query.compute(Ray.fromOriginDirection(p, dir), tri);
            expectClose(r.distance, 0, 1e-9, 1e-9);
        });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(rayArb, triArb, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([ray, tri, frame, shift]) => {
            const rot = (p: Vector): Vector =>
                v(frame[0].values[0] * p.values[0]
                    + frame[1].values[0] * p.values[1],
                    frame[0].values[1] * p.values[0]
                    + frame[1].values[1] * p.values[1]);
            const movedRay = Ray.fromOriginDirection(
                add(shift, rot(ray.origin)), rot(ray.direction));
            const movedTri = Triangle.fromVertices(add(shift, rot(tri.v[0])),
                add(shift, rot(tri.v[1])), add(shift, rot(tri.v[2])));
            const r0 = query.compute(ray, tri);
            const r1 = query.compute(movedRay, movedTri);
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
        });
    });
});
