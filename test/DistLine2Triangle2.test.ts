import { describe, expect, it } from 'vitest';
import { DistLine2Triangle2 } from '../src/DistLine2Triangle2.js';
import { Line } from '../src/Line.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { perp } from '../src/Vector2.js';
import { DistPointTriangle } from '../src/DistPointTriangle.js';
import { check, expectClose, expectVectorClose, fc, rotationFrame, seededRandom, wellScaledVector } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of
// DistLine2Triangle2.ts against the upstream header DistLine2Triangle2.h.
// ---------------------------------------------------------------------------

function rot2(R: readonly Vector[], p: Vector): Vector {
    return add(mul(p.values[0], R[0]), mul(p.values[1], R[1]));
}

/**
 * The exact line-triangle distance, computed independently of the query under
 * test. The function t -> distance(P + t*D, triangle) is convex (a norm
 * distance to a convex set composed with an affine map), so a ternary search
 * on a bracket containing the minimizer converges to the global minimum. The
 * bracket is |t - t_c| <= (h + R)/|D| around the foot t_c of the
 * perpendicular from the centroid, with h the perpendicular distance and R
 * the radius of a ball about the centroid containing the triangle.
 */
function referenceLineTriangleDistance(l: Line, tri: Triangle): number {
    const pointTri = new DistPointTriangle();
    const centroid = mul(1 / 3, add(tri.v[0], add(tri.v[1], tri.v[2])));
    let radius = 0;
    for (let i = 0; i < 3; ++i) {
        radius = Math.max(radius, length(sub(tri.v[i], centroid)));
    }
    const dd = dot(l.direction, l.direction);
    const tc = dot(l.direction, sub(centroid, l.origin)) / dd;
    const foot = add(l.origin, mul(tc, l.direction));
    const h = length(sub(centroid, foot));
    const half = (h + radius) / Math.sqrt(dd);
    let lo = tc - half - 1;
    let hi = tc + half + 1;
    const f = (t: number): number => pointTri.compute(
        add(l.origin, mul(t, l.direction)), tri).sqrDistance;
    for (let k = 0; k < 160; ++k) {
        const m0 = lo + (hi - lo) / 3;
        const m1 = hi - (hi - lo) / 3;
        if (f(m0) <= f(m1)) { hi = m1; }
        else { lo = m0; }
    }
    return Math.sqrt(Math.min(f(lo), f(hi), f(0.5 * (lo + hi))));
}

const nonUnitLine2 = fc.tuple(wellScaledVector(2, -8, 8),
    wellScaledVector(2, -3, 3))
    .filter(([, d]) => length(d) > 0.25)
    .map(([o, d]) => Line.fromOriginDirection(o, d));

const triangle2 = fc.tuple(wellScaledVector(2, -8, 8),
    wellScaledVector(2, -8, 8), wellScaledVector(2, -8, 8))
    .filter(([a, b, c]) => {
        const e0 = sub(b, a);
        const e1 = sub(c, a);
        const d00 = dot(e0, e0), d11 = dot(e1, e1), d01 = dot(e0, e1);
        return d00 > 0.25 && d11 > 0.25
            && d00 * d11 - d01 * d01 > 0.05 * d00 * d11;
    })
    .map(([a, b, c]) => Triangle.fromVertices(a, b, c));

describe('DistLine2Triangle2 verification', () => {
    it('does not return NaN when the crossed edge is numerically parallel to the line', () => {
        // Upstream's s = DotPerp(D, P - V0) / DotPerp(D, V1 - V0) rounds the
        // denominator to exactly 0 here although the vertex signs differ
        // (found by the rigid-motion property). Pre-fix: distance NaN.
        const line = Line.fromOriginDirection(
            Vector.fromArray([-1.7086881202421092e-14, 1.9960937490641202]),
            Vector.fromArray([-2.5613545596609416e-14, 2.9921808228827933]));
        const triangle = Triangle.fromVertexArray([
            Vector.fromArray([-6.848127733652727e-14, 7.999999999999985]),
            Vector.fromArray([0, 0]),
            Vector.fromArray([0.49999999999995026, 5.820550528229657])]);
        const result = new DistLine2Triangle2().compute(line, triangle);
        expect(Number.isFinite(result.distance)).toBe(true);
        expect(result.distance).toBe(0);
        expect(Number.isFinite(result.parameter)).toBe(true);
        for (const b of result.barycentric) {
            expect(b).toBeGreaterThanOrEqual(-1e-12);
            expect(b).toBeLessThanOrEqual(1 + 1e-12);
        }
    });

    const query = new DistLine2Triangle2();

    it('result is self consistent with valid barycentrics', () => {
        check(fc.tuple(nonUnitLine2, triangle2), ([l, tri]) => {
            const r = query.compute(l, tri);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            const diff = sub(r.closest[0], r.closest[1]);
            expectClose(r.sqrDistance, dot(diff, diff), 1e-12, 1e-12);
            expectClose(r.barycentric[0] + r.barycentric[1]
                + r.barycentric[2], 1, 1e-9, 1e-9);
            for (let i = 0; i < 3; ++i) {
                expect(r.barycentric[i]).toBeGreaterThanOrEqual(-1e-9);
                expect(r.barycentric[i]).toBeLessThanOrEqual(1 + 1e-9);
            }
            let q = mul(r.barycentric[0], tri.v[0]);
            q = add(q, mul(r.barycentric[1], tri.v[1]));
            q = add(q, mul(r.barycentric[2], tri.v[2]));
            expectVectorClose(r.closest[1], q, 1e-8, 1e-8);
            // closest[0] is on the line at the reported parameter.
            expectVectorClose(r.closest[0],
                add(l.origin, mul(r.parameter, l.direction)), 1e-8, 1e-8);
        });
    });

    it('matches an independent convex minimization along the line', () => {
        check(fc.tuple(nonUnitLine2, triangle2), ([l, tri]) => {
            expectClose(query.compute(l, tri).distance,
                referenceLineTriangleDistance(l, tri), 1e-7, 1e-7);
        }, 100);
    });

    it('reports zero distance exactly when the line meets the triangle',
        () => {
            check(fc.tuple(nonUnitLine2, triangle2), ([l, tri]) => {
                const r = query.compute(l, tri);
                // The sign test that drives the dispatch: the line separates
                // the triangle only when all three normal components share a
                // strict sign.
                const N = perp(l.direction);
                const s = [0, 1, 2].map(i =>
                    Math.sign(dot(N, sub(tri.v[i], l.origin))));
                const separated = (s[0] > 0 && s[1] > 0 && s[2] > 0)
                    || (s[0] < 0 && s[1] < 0 && s[2] < 0);
                if (separated) {
                    expect(r.distance).toBeGreaterThan(0);
                    // The closest triangle point is a vertex.
                    const isVertex = r.barycentric.some(b => b === 1);
                    expect(isVertex).toBe(true);
                }
                else {
                    expect(r.distance).toBe(0);
                }
            });
        });

    it('picks the vertex nearest the line when the triangle is on one side',
        () => {
            check(fc.tuple(nonUnitLine2, triangle2), ([l, tri]) => {
                const r = query.compute(l, tri);
                if (r.distance === 0) { return; }
                const N = perp(l.direction);
                const nl = length(N);
                let best = Infinity;
                for (let i = 0; i < 3; ++i) {
                    best = Math.min(best,
                        Math.abs(dot(N, sub(tri.v[i], l.origin))) / nl);
                }
                expectClose(r.distance, best, 1e-9, 1e-9);
            });
        });

    it('is minimal over sampled line/triangle point pairs', () => {
        const rand = seededRandom(0x51e3);
        check(fc.tuple(nonUnitLine2, triangle2), ([l, tri]) => {
            const r = query.compute(l, tri);
            for (let k = 0; k < 24; ++k) {
                const t = 40 * (rand() - 0.5);
                let b0 = rand();
                let b1 = rand();
                if (b0 + b1 > 1) { b0 = 1 - b0; b1 = 1 - b1; }
                let q = mul(b0, tri.v[0]);
                q = add(q, mul(b1, tri.v[1]));
                q = add(q, mul(1 - b0 - b1, tri.v[2]));
                const gap = length(
                    sub(add(l.origin, mul(t, l.direction)), q));
                expect(r.distance).toBeLessThanOrEqual(gap + 1e-9 * (1 + gap));
            }
        }, 60);
    });

    it('is equivariant under rigid motions of the plane', () => {
        check(fc.tuple(nonUnitLine2, triangle2, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([l, tri, R, tr]) => {
            const movedLine = Line.fromOriginDirection(
                add(rot2(R, l.origin), tr), rot2(R, l.direction));
            const movedTri = Triangle.fromVertices(
                add(rot2(R, tri.v[0]), tr), add(rot2(R, tri.v[1]), tr),
                add(rot2(R, tri.v[2]), tr));
            expectClose(query.compute(l, tri).distance,
                query.compute(movedLine, movedTri).distance, 1e-8, 1e-8);
        });
    });

    it('is permutation equivariant in the triangle vertices', () => {
        check(fc.tuple(nonUnitLine2, triangle2), ([l, tri]) => {
            const base = query.compute(l, tri);
            for (const perm of [[1, 2, 0], [2, 0, 1], [1, 0, 2]]) {
                const permuted = Triangle.fromVertices(tri.v[perm[0]],
                    tri.v[perm[1]], tri.v[perm[2]]);
                expectClose(query.compute(l, permuted).distance,
                    base.distance, 1e-9, 1e-9);
            }
        });
    });

    it('handles a line through a single vertex', () => {
        // Integer coordinates keep Dot(Perp(D), V - P) exactly zero at the
        // chosen vertex, so the "line contains a vertex" branches are the
        // ones exercised.
        check(fc.tuple(fc.integer({ min: -6, max: 6 }),
            fc.integer({ min: -6, max: 6 }), fc.integer({ min: 1, max: 6 }),
            fc.integer({ min: 1, max: 6 }), fc.integer({ min: 0, max: 2 })),
        ([ax, ay, w, h, which]) => {
            const tri = Triangle.fromVertices(v(ax, ay), v(ax + w, ay),
                v(ax, ay + h));
            const vertex = tri.v[which];
            const l = Line.fromOriginDirection(vertex, v(2, 1));
            const r = query.compute(l, tri);
            expect(r.distance).toBe(0);
            expectVectorClose(r.closest[0], r.closest[1], 0, 0);
        });
    });
});
