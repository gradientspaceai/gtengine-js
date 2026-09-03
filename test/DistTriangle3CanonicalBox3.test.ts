import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistTriangle3CanonicalBox3 }
    from '../src/DistTriangle3CanonicalBox3.js';
import type { DistTriangle3CanonicalBox3Result }
    from '../src/DistTriangle3CanonicalBox3.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

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
