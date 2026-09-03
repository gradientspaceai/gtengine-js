import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistTriangle3AlignedBox3 } from '../src/DistTriangle3AlignedBox3.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, length, mul, sub } from '../src/Vector.js';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// The distance from a point to a solid aligned box, computed by clamping.
function distPointBox(p: Vector, box: AlignedBox): number {
    let sqr = 0;
    for (let i = 0; i < 3; ++i) {
        const x = p.values[i];
        const lo = box.min.values[i], hi = box.max.values[i];
        const d = x < lo ? lo - x : (x > hi ? x - hi : 0);
        sqr += d * d;
    }
    return Math.sqrt(sqr);
}

// Brute-force minimum over a barycentric grid on the triangle.
function bruteForce(triangle: Triangle, box: AlignedBox, n: number): number {
    let best = Infinity;
    for (let i = 0; i <= n; ++i) {
        for (let j = 0; i + j <= n; ++j) {
            const b0 = i / n, b1 = j / n, b2 = 1 - b0 - b1;
            const p = add(add(mul(b0, triangle.v[0]), mul(b1, triangle.v[1])),
                mul(b2, triangle.v[2]));
            const d = distPointBox(p, box);
            if (d < best) {
                best = d;
            }
        }
    }
    return best;
}

const unitBox = () => AlignedBox.fromMinMax(v(-1, -1, -1), v(1, 1, 1));

describe('DistTriangle3AlignedBox3', () => {
    it('reports the separation of a triangle above a box face', () => {
        const triangle = Triangle.fromVertices(
            v(0, 0, 4), v(0.5, 0, 4), v(0, 0.5, 4));
        const result = new DistTriangle3AlignedBox3()
            .compute(triangle, unitBox());
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.sqrDistance).toBeCloseTo(9, 12);
        expect(result.closest[0].values[2]).toBeCloseTo(4, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });

    it('finds a vertex as the closest triangle point', () => {
        // The triangle is far from the box except for the vertex (2,0,0),
        // whose closest box point is (1,0,0).
        const triangle = Triangle.fromVertices(
            v(2, 0, 0), v(6, 3, 0), v(6, -3, 0));
        const result = new DistTriangle3AlignedBox3()
            .compute(triangle, unitBox());
        expect(result.distance).toBeCloseTo(1, 12);
        expect(result.closest[0].values).toEqual([2, 0, 0]);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
        expect(result.barycentric[0]).toBeCloseTo(1, 10);
        expect(result.barycentric[1]).toBeCloseTo(0, 10);
        expect(result.barycentric[2]).toBeCloseTo(0, 10);
    });

    it('reports zero distance for a triangle that pierces the box', () => {
        const triangle = Triangle.fromVertices(
            v(-3, -3, 0), v(3, -3, 0), v(0, 3, 0));
        const result = new DistTriangle3AlignedBox3()
            .compute(triangle, unitBox());
        expect(result.distance).toBe(0);
        expect(length(sub(result.closest[0], result.closest[1])))
            .toBeCloseTo(0, 12);
    });

    it('reports zero distance for a triangle inside the box', () => {
        const triangle = Triangle.fromVertices(
            v(-0.5, -0.5, 0), v(0.5, -0.5, 0), v(0, 0.5, 0));
        const result = new DistTriangle3AlignedBox3()
            .compute(triangle, unitBox());
        expect(result.distance).toBe(0);
    });

    it('translates the closest points out of the canonical box frame', () => {
        const box = AlignedBox.fromMinMax(v(10, 20, 30), v(12, 24, 36));
        const triangle = Triangle.fromVertices(
            v(11, 22, 40), v(11.5, 22, 40), v(11, 22.5, 40));
        const result = new DistTriangle3AlignedBox3().compute(triangle, box);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[0].values[2]).toBeCloseTo(40, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(36, 12);
        expect(distPointBox(result.closest[1], box)).toBeCloseTo(0, 12);
    });

    it('agrees with brute-force sampling and reports consistent points', () => {
        let seed = 13572468;
        const rnd = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };
        const query = new DistTriangle3AlignedBox3();
        for (let trial = 0; trial < 40; ++trial) {
            const pt = () => v(rnd() * 6 - 3, rnd() * 6 - 3, rnd() * 6 - 3);
            const triangle = Triangle.fromVertices(pt(), pt(), pt());
            const lo = v(rnd() - 2, rnd() - 2, rnd() - 2);
            const box = AlignedBox.fromMinMax(lo,
                add(lo, v(0.5 + rnd() * 2, 0.5 + rnd() * 2, 0.5 + rnd() * 2)));

            const result = query.compute(triangle, box);
            const brute = bruteForce(triangle, box, 40);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-9);
            expect(result.distance).toBeGreaterThan(brute - 0.15);
            expect(result.sqrDistance)
                .toBeCloseTo(result.distance * result.distance, 9);
            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 8);
            expect(distPointBox(result.closest[1], box)).toBeCloseTo(0, 8);

            // The barycentric coordinates are a convex combination that
            // reproduces the triangle point.
            const b = result.barycentric;
            expect(b[0] + b[1] + b[2]).toBeCloseTo(1, 8);
            for (const bi of b) {
                expect(bi).toBeGreaterThanOrEqual(-1e-9);
            }
            const fromBary = add(add(mul(b[0], triangle.v[0]),
                mul(b[1], triangle.v[1])), mul(b[2], triangle.v[2]));
            expect(length(sub(fromBary, result.closest[0])))
                .toBeCloseTo(0, 8);
        }
    });
});
