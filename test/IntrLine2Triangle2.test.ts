import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line';
import { Triangle } from '../src/Triangle';
import { Vector, add, mul, normalize } from '../src/Vector';
import { dotPerp } from '../src/Vector2';
import {
    IntrLine2Triangle2TI,
    IntrLine2Triangle2FI
} from '../src/IntrLine2Triangle2';

function line(px: number, py: number, dx: number, dy: number): Line {
    return Line.fromOriginDirection(Vector.fromArray([px, py]),
        Vector.fromArray([dx, dy]));
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

// Barycentric coordinates of p with respect to the triangle.
function barycentric(p: Vector, t: Triangle): number[] {
    const [v0, v1, v2] = t.v;
    const d = (v1.values[0] - v0.values[0]) * (v2.values[1] - v0.values[1]) -
        (v2.values[0] - v0.values[0]) * (v1.values[1] - v0.values[1]);
    const b1 = ((p.values[0] - v0.values[0]) * (v2.values[1] - v0.values[1]) -
        (v2.values[0] - v0.values[0]) * (p.values[1] - v0.values[1])) / d;
    const b2 = ((v1.values[0] - v0.values[0]) * (p.values[1] - v0.values[1]) -
        (p.values[0] - v0.values[0]) * (v1.values[1] - v0.values[1])) / d;
    return [1 - b1 - b2, b1, b2];
}

describe('IntrLine2Triangle2', () => {
    const ti = new IntrLine2Triangle2TI();
    const fi = new IntrLine2Triangle2FI();
    const tri = triangle([0, 0], [4, 0], [0, 4]);

    it('clips a line crossing the interior, (n,p,z) = (2,1,0)', () => {
        // The line y = 1 crosses the triangle from (0,1) to (3,1).
        const l = line(0, 1, 1, 0);
        expect(ti.test(l, tri).intersect).toBe(true);
        const result = fi.find(l, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(3, 12);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(3, 12);
    });

    it('clips a line through a vertex and the opposite edge, (1,1,1)', () => {
        // The line through (0,0) with direction (1,1) exits at (2,2).
        const l = line(0, 0, 1, 1);
        const result = fi.find(l, tri);
        expect(result.numIntersections).toBe(2);
        // The direction is not unit length, so the parameters are divided by
        // |D|^2 = 2.
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 12);
        expect(result.point[1].values[0]).toBeCloseTo(2, 12);
        expect(result.point[1].values[1]).toBeCloseTo(2, 12);
    });

    it('reports a single point when the line touches a vertex, (2,0,1)', () => {
        const l = line(0, 4, 1, 0);
        expect(ti.test(l, tri).intersect).toBe(true);
        const result = fi.find(l, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBe(result.parameter[1]);
        expect(result.point[0].values).toEqual([0, 4]);
    });

    it('reports the whole edge when the line contains one, (0,1,2)', () => {
        const l = line(0, 0, 1, 0);
        expect(ti.test(l, tri).intersect).toBe(true);
        const result = fi.find(l, tri);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(4, 12);
    });

    it('reports no intersection when all vertices are on one side', () => {
        for (const l of [line(0, -1, 1, 0), line(0, 5, 1, 0)]) {
            expect(ti.test(l, tri).intersect).toBe(false);
            const result = fi.find(l, tri);
            expect(result.intersect).toBe(false);
            expect(result.numIntersections).toBe(0);
        }
    });

    it('treats a degenerate triangle on the line as no intersection, (0,0,3)', () => {
        const collinear = triangle([0, 0], [1, 0], [2, 0]);
        const l = line(0, 0, 1, 0);
        expect(ti.test(l, collinear).intersect).toBe(false);
        expect(fi.find(l, collinear).intersect).toBe(false);
    });

    it('orders the parameters and is direction-reversal consistent', () => {
        const forward = fi.find(line(0, 1, 1, 0), tri);
        const backward = fi.find(line(0, 1, -1, 0), tri);
        expect(forward.parameter[0]).toBeLessThanOrEqual(forward.parameter[1]);
        expect(backward.parameter[0]).toBeLessThanOrEqual(backward.parameter[1]);
        // The intersection segment is the same set of points.
        const fPoints = [forward.point[0].values, forward.point[1].values]
            .sort((a, b) => a[0] - b[0]);
        const bPoints = [backward.point[0].values, backward.point[1].values]
            .sort((a, b) => a[0] - b[0]);
        expect(fPoints[0][0]).toBeCloseTo(bPoints[0][0], 12);
        expect(fPoints[1][0]).toBeCloseTo(bPoints[1][0], 12);
    });

    it('keeps TI and FI consistent and clips inside the triangle', () => {
        const rand = makeRandom(13571113);
        let numHit = 0, numMiss = 0;
        for (let trial = 0; trial < 800; ++trial) {
            const t = triangle(
                [4 * rand() - 2, 4 * rand() - 2],
                [4 * rand() - 2, 4 * rand() - 2],
                [4 * rand() - 2, 4 * rand() - 2]);
            const d = Vector.fromArray([2 * rand() - 1, 2 * rand() - 1]);
            if (normalize(d) < 1e-6) {
                continue;
            }
            const l = Line.fromOriginDirection(
                Vector.fromArray([4 * rand() - 2, 4 * rand() - 2]), d);

            // Independent oracle: the line meets the solid triangle when the
            // vertices are not all strictly on one side.
            const s = [0, 1, 2].map(i =>
                dotPerp(l.direction, Vector.fromArray([
                    t.v[i].values[0] - l.origin.values[0],
                    t.v[i].values[1] - l.origin.values[1]])));
            const oracle = !(s.every(x => x > 0) || s.every(x => x < 0) ||
                s.every(x => x === 0));

            const tiResult = ti.test(l, t).intersect;
            const f = fi.find(l, t);
            expect(tiResult).toBe(oracle);
            expect(f.intersect).toBe(tiResult);

            if (f.intersect) {
                ++numHit;
                expect(f.parameter[0]).toBeLessThanOrEqual(f.parameter[1]);
                // Sampled points of the reported interval are in the triangle.
                for (let k = 0; k <= 4; ++k) {
                    const u = f.parameter[0] +
                        (k / 4) * (f.parameter[1] - f.parameter[0]);
                    const p = add(l.origin, mul(u, l.direction));
                    for (const b of barycentric(p, t)) {
                        expect(b).toBeGreaterThan(-1e-8);
                    }
                }
            } else {
                ++numMiss;
                expect(f.numIntersections).toBe(0);
            }
        }
        expect(numHit).toBeGreaterThan(50);
        expect(numMiss).toBeGreaterThan(50);
    });
});
