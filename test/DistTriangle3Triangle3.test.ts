import { describe, expect, it } from 'vitest';
import { DistTriangle3Triangle3 } from '../src/DistTriangle3Triangle3.js';
import type { DistTriangle3Triangle3Result }
    from '../src/DistTriangle3Triangle3.js';
import { Triangle } from '../src/Triangle.js';
import { DistSegment3Triangle3 } from '../src/DistSegment3Triangle3.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, latticeVector, rotationFrame,
    seededRandom, wellScaledVector
} from './helpers/arbitraries.js';
import { cross } from '../src/Vector3.js';

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

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistTriangle3Triangle3.h.
// ---------------------------------------------------------------------------

describe('DistTriangle3Triangle3 verification', () => {
    const query = new DistTriangle3Triangle3();
    const stQuery = new DistSegment3Triangle3();

    const triArb = fc.tuple(latticeVector(3, -6, 6), latticeVector(3, -6, 6),
        latticeVector(3, -6, 6))
        .filter(([a, b, c]) => {
            const n = cross(sub(b, a), sub(c, a));
            return dot(n, n) > 4;
        })
        .map(([a, b, c]) => Triangle.fromVertices(a, b, c));

    function maxEdge(t: Triangle): number {
        let m = 0;
        for (let i = 0; i < 3; ++i) {
            const e = sub(t.v[(i + 1) % 3], t.v[i]);
            m = Math.max(m, Math.sqrt(dot(e, e)));
        }
        return m;
    }

    it('reports consistent distances and barycentric closest points', () => {
        check(fc.tuple(triArb, triArb), ([t0, t1]) => {
            const r = query.compute(t0, t1);
            // The absolute tolerance is 1e-6: these queries accumulate the
            // squared distance while clamping to faces and edges, so a
            // near-touching configuration loses about half the mantissa and
            // the distance carries an absolute error of order sqrt(eps)
            // times the coordinate scale. A translation or frame error
            // would show up as an O(1) discrepancy.
            expectClose(r.sqrDistance, r.distance * r.distance, 1e-12, 1e-12);
            const d = sub(r.closest[0], r.closest[1]);
            expectClose(Math.sqrt(dot(d, d)), r.distance, 1e-6, 1e-8);

            const verify = (b: [number, number, number], t: Triangle,
                c: Vector) => {
                expectClose(b[0] + b[1] + b[2], 1, 1e-9, 1e-9);
                let rebuilt = new Vector(3);
                for (let i = 0; i < 3; ++i) {
                    expect(b[i]).toBeGreaterThanOrEqual(-1e-8);
                    rebuilt = add(rebuilt, mul(b[i], t.v[i]));
                }
                expectVectorClose(rebuilt, c, 1e-7, 1e-7);
            };
            verify(r.barycentric0, t0, r.closest[0]);
            verify(r.barycentric1, t1, r.closest[1]);
        });
    });

    it('is symmetric under argument swap', () => {
        check(fc.tuple(triArb, triArb), ([t0, t1]) => {
            const r0 = query.compute(t0, t1);
            const r1 = query.compute(t1, t0);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            // The distance is symmetric; the reported pair need not be the
            // same representative when several pairs are equidistant, so
            // only check that the swapped pair realizes the same distance.
            const d = sub(r1.closest[1], r1.closest[0]);
            expectClose(Math.sqrt(dot(d, d)), r0.distance, 1e-6, 1e-8);
        });
    });

    // The minimum distance between two solid triangles is attained on the
    // boundary of at least one of them, so it equals the minimum over the six
    // (edge, other triangle) queries.
    it('equals the minimum over the six edge-triangle queries', () => {
        check(fc.tuple(triArb, triArb), ([t0, t1]) => {
            let best = Number.MAX_VALUE;
            for (let i = 0; i < 3; ++i) {
                best = Math.min(best, stQuery.compute(
                    Segment.fromEndpoints(t0.v[i], t0.v[(i + 1) % 3]),
                    t1).sqrDistance);
                best = Math.min(best, stQuery.compute(
                    Segment.fromEndpoints(t1.v[i], t1.v[(i + 1) % 3]),
                    t0).sqrDistance);
            }
            const r = query.compute(t0, t1);
            expectClose(r.sqrDistance, best, 1e-9, 1e-9);
        });
    });

    // Two-sided check against a barycentric sampling of triangle0. The
    // sampling is an upper bound for the true distance, and because the
    // distance function is 1-Lipschitz and the sample spacing is at most
    // maxEdge/n, the upper bound is at most maxEdge/n above the true value.
    it('brackets a barycentric sampling of the first triangle', () => {
        const rng = seededRandom(0x1a2b3c4d);
        const n = 40;
        for (let k = 0; k < 25; ++k) {
            const p = () => v(8 * rng() - 4, 8 * rng() - 4, 8 * rng() - 4);
            const t0 = Triangle.fromVertices(p(), p(), p());
            const t1 = Triangle.fromVertices(p(), p(), p());
            const r = query.compute(t0, t1);
            const s = sampledDistance(t0, t1, n);
            expect(r.distance).toBeLessThanOrEqual(s + 1e-9);
            expect(s - r.distance).toBeLessThanOrEqual(maxEdge(t0) / n + 1e-9);
        }
    }, 30000);

    it('reports zero for intersecting triangles', () => {
        const t0 = tri([0, 0, 0], [4, 0, 0], [0, 4, 0]);
        const t1 = tri([1, 1, -1], [1, 1, 1], [3, 1, 0]);
        const r = query.compute(t0, t1);
        expect(r.distance).toBeCloseTo(0, 9);
        expectVectorClose(r.closest[0], r.closest[1], 1e-8, 1e-8);
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(triArb, triArb, rotationFrame(3),
            wellScaledVector(3, -5, 5)), ([t0, t1, R, tr]) => {
                const xf = (q: Vector) => add(add(add(
                    mul(q.values[0], R[0]), mul(q.values[1], R[1])),
                    mul(q.values[2], R[2])), tr);
                const xt = (t: Triangle) => Triangle.fromVertices(
                    xf(t.v[0]), xf(t.v[1]), xf(t.v[2]));
                const r0 = query.compute(t0, t1);
                const r1 = query.compute(xt(t0), xt(t1));
                expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
            });
    });
});
