import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistTriangle3CanonicalBox3 }
    from '../src/DistTriangle3CanonicalBox3.js';
import type { DistTriangle3CanonicalBox3Result }
    from '../src/DistTriangle3CanonicalBox3.js';
import { Triangle } from '../src/Triangle.js';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistSegment3CanonicalBox3 } from '../src/DistSegment3CanonicalBox3.js';
import { DistTriangle3AlignedBox3 } from '../src/DistTriangle3AlignedBox3.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, latticeVector,
    seededRandom
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function tri(v0: number[], v1: number[], v2: number[]): Triangle {
    return Triangle.fromVertices(v(...v0), v(...v1), v(...v2));
}

// The exact distance from a point to a solid canonical box.
function pointBoxDistance(p: Vector, box: CanonicalBox): number {
    let sqr = 0;
    for (let i = 0; i < 3; ++i) {
        const d = Math.abs(p.values[i]) - box.extent.values[i];
        if (d > 0) {
            sqr += d * d;
        }
    }
    return Math.sqrt(sqr);
}

// A barycentric grid sampling of the solid triangle. Because the distance
// from each sample to the box is exact, the result is an upper bound for the
// triangle-box distance whose error is at most the grid spacing.
function sampledDistance(t: Triangle, box: CanonicalBox, n: number): number {
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= n; ++i) {
        for (let j = 0; i + j <= n; ++j) {
            const k = n - i - j;
            let p = mul(i / n, t.v[0]);
            p = add(p, mul(j / n, t.v[1]));
            p = add(p, mul(k / n, t.v[2]));
            best = Math.min(best, pointBoxDistance(p, box));
        }
    }
    return best;
}

function expectConsistent(t: Triangle, box: CanonicalBox,
    result: DistTriangle3CanonicalBox3Result): void {
    const delta = sub(result.closest[0], result.closest[1]);
    expect(Math.sqrt(dot(delta, delta))).toBeCloseTo(result.distance, 8);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 8);

    // The barycentric coordinates are a partition of unity, are nonnegative
    // and reproduce closest[0].
    const b = result.barycentric;
    expect(b[0] + b[1] + b[2]).toBeCloseTo(1, 8);
    let reconstructed = new Vector(3);
    for (let i = 0; i < 3; ++i) {
        expect(b[i]).toBeGreaterThanOrEqual(-1e-8);
        reconstructed = add(reconstructed, mul(b[i], t.v[i]));
    }
    for (let i = 0; i < 3; ++i) {
        expect(reconstructed.values[i])
            .toBeCloseTo(result.closest[0].values[i], 7);
    }

    // closest[1] is in the box.
    for (let i = 0; i < 3; ++i) {
        expect(Math.abs(result.closest[1].values[i]))
            .toBeLessThanOrEqual(box.extent.values[i] + 1e-8);
    }
}

describe('DistTriangle3CanonicalBox3', () => {
    const query = new DistTriangle3CanonicalBox3();
    const unit = CanonicalBox.fromExtent(v(1, 1, 1));

    it('computes the distance for a triangle parallel to a box face', () => {
        // The triangle lies in z = 5 and its projection covers the origin, so
        // the closest points differ only in z.
        const t = tri([-2, -2, 5], [2, -2, 5], [0, 2, 5]);
        const result = query.compute(t, unit);
        expect(result.distance).toBeCloseTo(4, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 10);
        expectConsistent(t, unit, result);
    });

    it('reports zero distance when the triangle intersects the box', () => {
        const t = tri([-3, 0, 0], [3, 0, 0], [0, 3, 0]);
        const result = query.compute(t, unit);
        expect(result.distance).toBeCloseTo(0, 10);
        expectConsistent(t, unit, result);
    });

    it('computes the distance for a triangle beyond a box vertex', () => {
        // A small triangle centered on the line x = y = z beyond (1,1,1).
        const t = tri([3, 3, 3], [3.5, 3, 3], [3, 3.5, 3]);
        const result = query.compute(t, unit);
        const expected = Math.sqrt(3) * 2;
        expect(result.distance).toBeCloseTo(expected, 10);
        expect(result.barycentric[0]).toBeCloseTo(1, 8);
        expectConsistent(t, unit, result);
    });

    it('handles a degenerate (segment-like) triangle', () => {
        const t = tri([4, 0, 0], [4, 2, 0], [4, 1, 0]);
        const result = query.compute(t, unit);
        expect(result.distance).toBeCloseTo(3, 10);
    });

    it('matches a dense sampling of the triangle on random inputs', () => {
        let seed = 777333111;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const box = CanonicalBox.fromExtent(v(1, 1.5, 0.75));
        for (let trial = 0; trial < 25; ++trial) {
            const t = tri(
                [rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4]);
            const result = query.compute(t, box);
            const sampled = sampledDistance(t, box, 40);
            expect(result.distance).toBeLessThanOrEqual(sampled + 1e-8);
            expect(sampled - result.distance).toBeLessThan(0.35);
            expectConsistent(t, box, result);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistTriangle3CanonicalBox3.h.
// ---------------------------------------------------------------------------

describe('DistTriangle3CanonicalBox3 verification', () => {
    const query = new DistTriangle3CanonicalBox3();
    const abQuery = new DistTriangle3AlignedBox3();
    const sbQuery = new DistSegment3CanonicalBox3();

    const boxArb = fc.array(finite(0.05, 4), { minLength: 3, maxLength: 3 })
        .map(e => CanonicalBox.fromExtent(v(e[0], e[1], e[2])));

    // Lattice vertices keep the triangle normal well conditioned; the area
    // filter rejects the needle triangles whose Normalize(K) underflows.
    const triArb = fc.tuple(latticeVector(3, -6, 6), latticeVector(3, -6, 6),
        latticeVector(3, -6, 6))
        .filter(([a, b, c]) => {
            const n = cross(sub(b, a), sub(c, a));
            return dot(n, n) > 4;
        })
        .map(([a, b, c]) => Triangle.fromVertices(a, b, c));

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(triArb, boxArb), ([t, box]) => {
            const res = query.compute(t, box);
            // The absolute tolerance is 1e-6: these queries accumulate the
            // squared distance while clamping to faces and edges, so a
            // near-touching configuration loses about half the mantissa and
            // the distance carries an absolute error of order sqrt(eps)
            // times the coordinate scale. A translation or frame error
            // would show up as an O(1) discrepancy.
            expectClose(res.sqrDistance, res.distance * res.distance,
                1e-12, 1e-12);
            const d = sub(res.closest[0], res.closest[1]);
            expectClose(Math.sqrt(dot(d, d)), res.distance, 1e-6, 1e-8);

            const b = res.barycentric;
            expectClose(b[0] + b[1] + b[2], 1, 1e-9, 1e-9);
            let rebuilt = new Vector(3);
            for (let i = 0; i < 3; ++i) {
                expect(b[i]).toBeGreaterThanOrEqual(-1e-8);
                rebuilt = add(rebuilt, mul(b[i], t.v[i]));
            }
            expectVectorClose(rebuilt, res.closest[0], 1e-7, 1e-7);

            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(res.closest[1].values[i]))
                    .toBeLessThanOrEqual(box.extent.values[i] + 1e-8);
            }
        });
    });

    it('agrees with the exact point-box distance at closest[0]', () => {
        check(fc.tuple(triArb, boxArb), ([t, box]) => {
            const res = query.compute(t, box);
            expectClose(res.distance, pointBoxDistance(res.closest[0], box),
                1e-7, 1e-7);
        });
    });

    // When the closest plane point is outside the triangle, the upstream
    // fall-back is the minimum over the three edges; the reported distance
    // must then equal that minimum. When the plane point is inside, the
    // distance can only be smaller.
    it('never exceeds the minimum over the three edges', () => {
        check(fc.tuple(triArb, boxArb), ([t, box]) => {
            let best = Number.MAX_VALUE;
            for (let i = 0; i < 3; ++i) {
                const s = Segment.fromEndpoints(t.v[i], t.v[(i + 1) % 3]);
                best = Math.min(best, sbQuery.compute(s, box).sqrDistance);
            }
            const res = query.compute(t, box);
            expect(res.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
        });
    });

    it('never exceeds a barycentric sampling of the triangle', () => {
        const rng = seededRandom(0x0badf00d);
        const box = CanonicalBox.fromExtent(v(1, 2, 0.5));
        for (let k = 0; k < 40; ++k) {
            const p = () => v(10 * rng() - 5, 10 * rng() - 5, 10 * rng() - 5);
            const t = Triangle.fromVertices(p(), p(), p());
            const res = query.compute(t, box);
            expect(res.distance)
                .toBeLessThanOrEqual(sampledDistance(t, box, 30) + 1e-9);
        }
    }, 30000);

    it('agrees with the aligned-box query for an origin-centered box', () => {
        check(fc.tuple(triArb, boxArb), ([t, box]) => {
            const e = box.extent;
            const ab = AlignedBox.fromMinMax(
                v(-e.values[0], -e.values[1], -e.values[2]),
                v(e.values[0], e.values[1], e.values[2]));
            const r0 = query.compute(t, box);
            const r1 = abQuery.compute(t, ab);
            expectClose(r0.distance, r1.distance, 1e-12, 1e-12);
            expectVectorClose(r0.closest[0], r1.closest[0], 1e-12, 1e-12);
            expectVectorClose(r0.closest[1], r1.closest[1], 1e-12, 1e-12);
        });
    });

    it('reports zero when a vertex is inside the box', () => {
        check(fc.tuple(triArb, boxArb), ([t, box]) => {
            // Translate the triangle so that vertex 0 lands at the origin,
            // which is inside every canonical box.
            const shift = Triangle.fromVertices(sub(t.v[0], t.v[0]),
                sub(t.v[1], t.v[0]), sub(t.v[2], t.v[0]));
            const res = query.compute(shift, box);
            expect(res.distance).toBe(0);
            expect(res.sqrDistance).toBe(0);
        });
    });

    it('is invariant under a permutation of the triangle vertices', () => {
        check(fc.tuple(triArb, boxArb), ([t, box]) => {
            const r0 = query.compute(t, box);
            const r1 = query.compute(
                Triangle.fromVertices(t.v[1], t.v[2], t.v[0]), box);
            const r2 = query.compute(
                Triangle.fromVertices(t.v[1], t.v[0], t.v[2]), box);
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
            expectClose(r0.distance, r2.distance, 1e-8, 1e-8);
        });
    });
});
