import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, mul, normalize } from '../src/Vector.js';
import { IntrLine2Triangle2FI } from '../src/IntrLine2Triangle2.js';
import {
    IntrRay2Triangle2TI,
    IntrRay2Triangle2FI
} from '../src/IntrRay2Triangle2.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function ray(p: number[], d: number[]): Ray {
    const dir = vec(d);
    normalize(dir);
    return Ray.fromOriginDirection(vec(p), dir);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The triangle (0,0), (4,0), (0,4).
const tri = Triangle.fromVertices(vec([0, 0]), vec([4, 0]), vec([0, 4]));

function insideTriangle(t: Triangle, x: Vector): boolean {
    // Barycentric sign test with a small tolerance.
    const [v0, v1, v2] = t.v;
    const d = (v1.values[1] - v2.values[1]) * (v0.values[0] - v2.values[0])
        + (v2.values[0] - v1.values[0]) * (v0.values[1] - v2.values[1]);
    const b0 = ((v1.values[1] - v2.values[1]) * (x.values[0] - v2.values[0])
        + (v2.values[0] - v1.values[0]) * (x.values[1] - v2.values[1])) / d;
    const b1 = ((v2.values[1] - v0.values[1]) * (x.values[0] - v2.values[0])
        + (v0.values[0] - v2.values[0]) * (x.values[1] - v2.values[1])) / d;
    const b2 = 1 - b0 - b1;
    const eps = 1e-12;
    return b0 >= -eps && b1 >= -eps && b2 >= -eps;
}

describe('IntrRay2Triangle2', () => {
    const ti = new IntrRay2Triangle2TI();
    const fi = new IntrRay2Triangle2FI();

    it('finds the chord of a ray that crosses the triangle', () => {
        // The ray y = 1, x increasing from -3, crosses the triangle from
        // x = 0 to x = 3.
        const r = ray([-3, 1], [1, 0]);
        expect(ti.test(r, tri).intersect).toBe(true);
        const result = fi.find(r, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(3, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(3, 12);
    });

    it('clips the near end when the ray origin is inside the triangle', () => {
        const r = ray([1, 1], [1, 0]);
        const result = fi.find(r, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 12);
    });

    it('reports no intersection when the ray points away', () => {
        const r = ray([-3, 1], [-1, 0]);
        expect(ti.test(r, tri).intersect).toBe(false);
        const result = fi.find(r, tri);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
        expect(result.point[0].values).toEqual([0, 0]);
    });

    it('reports a single point when the ray touches only a vertex', () => {
        // The line y = 4 supports the triangle at the vertex (0,4); the ray
        // reaches it at t = 3.
        const r = ray([-3, 4], [1, 0]);
        const result = fi.find(r, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(3, 12);
        expect(result.parameter[1]).toBeCloseTo(3, 12);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[1]).toBeCloseTo(4, 12);
    });

    it('handles a ray collinear with an edge, clipped at the origin', () => {
        // The ray starts at the vertex (4,0) and runs along the supporting
        // line of the edge from (0,0) to (4,0), leaving the triangle at once.
        const r = ray([4, 0], [1, 0]);
        const result = fi.find(r, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(0, 12);
    });

    it('misses a triangle that lies behind and to the side', () => {
        const r = ray([-1, -1], [0, -1]);
        expect(ti.test(r, tri).intersect).toBe(false);
        expect(fi.find(r, tri).intersect).toBe(false);
    });

    it('is the line query clipped to t >= 0', () => {
        const rand = makeRandom(2718281);
        const lineFI = new IntrLine2Triangle2FI();
        for (let trial = 0; trial < 400; ++trial) {
            const r = ray([10 * rand() - 5, 10 * rand() - 5],
                [2 * rand() - 1, 2 * rand() - 1]);
            const l = Line.fromOriginDirection(r.origin, r.direction);
            const lineResult = lineFI.find(l, tri);
            const rayResult = fi.find(r, tri);
            expect(ti.test(r, tri).intersect).toBe(rayResult.intersect);

            if (!lineResult.intersect || lineResult.parameter[1] < 0) {
                expect(rayResult.intersect).toBe(false);
            }
            else {
                expect(rayResult.intersect).toBe(true);
                expect(rayResult.parameter[0]).toBeCloseTo(
                    Math.max(lineResult.parameter[0], 0), 12);
                expect(rayResult.parameter[1]).toBeCloseTo(
                    lineResult.parameter[1], 12);
            }
        }
    });

    it('agrees with dense sampling along the ray', () => {
        const rand = makeRandom(161803);
        for (let trial = 0; trial < 80; ++trial) {
            const r = ray([8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1]);
            const result = fi.find(r, tri);

            let tLo = Number.POSITIVE_INFINITY;
            let tHi = Number.NEGATIVE_INFINITY;
            const n = 20000;
            for (let k = 0; k <= n; ++k) {
                const t = (16 * k) / n;
                const x = add(r.origin, mul(t, r.direction));
                if (insideTriangle(tri, x)) {
                    if (t < tLo) { tLo = t; }
                    if (t > tHi) { tHi = t; }
                }
            }

            if (tLo <= tHi) {
                expect(result.intersect).toBe(true);
                expect(result.parameter[0]).toBeLessThanOrEqual(tLo + 1e-9);
                expect(result.parameter[1]).toBeGreaterThanOrEqual(tHi - 1e-9);
                expect(tLo - result.parameter[0]).toBeLessThan(3e-3);
                expect(result.parameter[1] - tHi).toBeLessThan(3e-3);
            }
        }
    });
});
