import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, dot, normalize } from '../src/Vector.js';
import {
    IntrHalfspace3Triangle3TI,
    IntrHalfspace3Triangle3FI
} from '../src/IntrHalfspace3Triangle3.js';

function halfspace(nx: number, ny: number, nz: number, c: number): Halfspace {
    const n = Vector.fromArray([nx, ny, nz]);
    normalize(n);
    return Halfspace.fromNormalConstant(n, c);
}

function triangle(a: number[], b: number[], c: number[]): Triangle {
    return Triangle.fromVertices(Vector.fromArray(a), Vector.fromArray(b),
        Vector.fromArray(c));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The area of the polygon formed by the first n result points, assumed convex
// and given in order.
function polygonArea(points: Vector[], n: number): number {
    let area = 0;
    for (let i = 1; i + 1 < n; ++i) {
        const u = [
            points[i].values[0] - points[0].values[0],
            points[i].values[1] - points[0].values[1],
            points[i].values[2] - points[0].values[2]
        ];
        const v = [
            points[i + 1].values[0] - points[0].values[0],
            points[i + 1].values[1] - points[0].values[1],
            points[i + 1].values[2] - points[0].values[2]
        ];
        area += 0.5 * Math.hypot(
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0]);
    }
    return area;
}

describe('IntrHalfspace3Triangle3', () => {
    const ti = new IntrHalfspace3Triangle3TI();
    const fi = new IntrHalfspace3Triangle3FI();
    const zUp = halfspace(0, 0, 1, 0);  // z >= 0

    it('returns the original triangle when it is inside the halfspace', () => {
        const tri = triangle([0, 0, 1], [1, 0, 2], [0, 1, 3]);
        expect(ti.test(zUp, tri).intersect).toBe(true);
        const result = fi.find(zUp, tri);
        expect(result.numPoints).toBe(3);
        expect(result.point[0].values).toEqual([0, 0, 1]);
        expect(result.point[2].values).toEqual([0, 1, 3]);
    });

    it('returns nothing when the triangle is strictly outside', () => {
        const tri = triangle([0, 0, -1], [1, 0, -2], [0, 1, -3]);
        expect(ti.test(zUp, tri).intersect).toBe(false);
        const result = fi.find(zUp, tri);
        expect(result.intersect).toBe(false);
        expect(result.numPoints).toBe(0);
    });

    it('clips one vertex away, yielding a quadrilateral', () => {
        // (n,p,z) = (1,2,0). Vertex 0 is below the plane.
        const tri = triangle([0, 0, -1], [2, 0, 1], [0, 2, 1]);
        const result = fi.find(zUp, tri);
        expect(result.intersect).toBe(true);
        expect(result.numPoints).toBe(4);
        // Every reported point is in the closed halfspace.
        for (let i = 0; i < 4; ++i) {
            expect(result.point[i].values[2]).toBeGreaterThanOrEqual(-1e-15);
        }
        // The two positive vertices are kept verbatim.
        expect(result.point[0].values).toEqual([2, 0, 1]);
        expect(result.point[1].values).toEqual([0, 2, 1]);
        // The clipped points are the edge midpoints in this configuration.
        expect(result.point[2].values).toEqual([0, 1, 0]);
        expect(result.point[3].values).toEqual([1, 0, 0]);
        // The quadrilateral has 3/4 the area of the whole triangle.
        const whole = polygonArea([tri.v[0], tri.v[1], tri.v[2]], 3);
        expect(polygonArea(result.point, 4)).toBeCloseTo(0.75 * whole, 10);
    });

    it('clips two vertices away, yielding a triangle', () => {
        // (n,p,z) = (2,1,0). Vertex 0 is above the plane.
        const tri = triangle([0, 0, 1], [2, 0, -1], [0, 2, -1]);
        const result = fi.find(zUp, tri);
        expect(result.numPoints).toBe(3);
        expect(result.point[0].values).toEqual([0, 0, 1]);
        expect(result.point[1].values).toEqual([1, 0, 0]);
        expect(result.point[2].values).toEqual([0, 1, 0]);
        const whole = polygonArea([tri.v[0], tri.v[1], tri.v[2]], 3);
        expect(polygonArea(result.point, 3)).toBeCloseTo(0.25 * whole, 10);
    });

    it('handles one vertex on the plane with one on each side', () => {
        // (n,p,z) = (1,1,1).
        const tri = triangle([0, 0, 0], [1, 0, 2], [0, 1, -2]);
        const result = fi.find(zUp, tri);
        expect(result.numPoints).toBe(3);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values).toEqual([1, 0, 2]);
        expect(result.point[2].values[2]).toBeCloseTo(0, 12);
    });

    it('reports the single edge on the plane when the rest is outside', () => {
        // (n,p,z) = (1,0,2).
        const tri = triangle([0, 0, 0], [1, 0, 0], [0, 1, -1]);
        const result = fi.find(zUp, tri);
        expect(result.intersect).toBe(true);
        expect(result.numPoints).toBe(2);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values).toEqual([1, 0, 0]);
    });

    it('reports the single vertex on the plane when the rest is outside', () => {
        // (n,p,z) = (2,0,1).
        const tri = triangle([0, 0, -1], [1, 0, 0], [0, 1, -1]);
        const result = fi.find(zUp, tri);
        expect(result.intersect).toBe(true);
        expect(result.numPoints).toBe(1);
        expect(result.point[0].values).toEqual([1, 0, 0]);
    });

    it('keeps TI and FI consistent and clips into the halfspace', () => {
        const rand = makeRandom(90210);
        let numHit = 0, numMiss = 0;
        const counts = new Map<number, number>();
        for (let trial = 0; trial < 600; ++trial) {
            const h = halfspace(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1,
                2 * rand() - 1);
            const tri = triangle(
                [4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2],
                [4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2],
                [4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2]);

            const t = ti.test(h, tri).intersect;
            const f = fi.find(h, tri);
            expect(f.intersect).toBe(t);

            const s = [0, 1, 2].map(i => dot(h.normal, tri.v[i]) - h.constant);
            expect(t).toBe(Math.max(s[0], s[1], s[2]) >= 0);

            counts.set(f.numPoints, (counts.get(f.numPoints) ?? 0) + 1);
            if (f.intersect) {
                ++numHit;
                for (let i = 0; i < f.numPoints; ++i) {
                    expect(dot(h.normal, f.point[i]) - h.constant)
                        .toBeGreaterThan(-1e-12);
                }
                // The clipped polygon never exceeds the triangle's area.
                const whole = polygonArea([tri.v[0], tri.v[1], tri.v[2]], 3);
                expect(polygonArea(f.point, f.numPoints))
                    .toBeLessThanOrEqual(whole + 1e-9);
            } else {
                ++numMiss;
                expect(f.numPoints).toBe(0);
            }
        }
        expect(numHit).toBeGreaterThan(50);
        expect(numMiss).toBeGreaterThan(50);
        // Whole triangles, clipped triangles and quadrilaterals all occur.
        expect(counts.get(3) ?? 0).toBeGreaterThan(0);
        expect(counts.get(4) ?? 0).toBeGreaterThan(0);
    });
});
