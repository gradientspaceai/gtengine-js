import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace';
import { Vector, dot, normalize } from '../src/Vector';
import { IntrHalfspace2Polygon2FI } from '../src/IntrHalfspace2Polygon2';

function halfspace(nx: number, ny: number, c: number): Halfspace {
    const n = Vector.fromArray([nx, ny]);
    normalize(n);
    return Halfspace.fromNormalConstant(n, c);
}

function poly(points: number[][]): Vector[] {
    return points.map(p => Vector.fromArray(p));
}

function area(points: readonly Vector[]): number {
    let a = 0;
    for (let i = 0; i < points.length; ++i) {
        const p = points[i], q = points[(i + 1) % points.length];
        a += p.values[0] * q.values[1] - q.values[0] * p.values[1];
    }
    return 0.5 * a;
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// An independent Sutherland-Hodgman clip of a convex polygon by one halfplane
// dot(n, x) >= c.
function clipOracle(points: readonly Vector[], h: Halfspace): Vector[] {
    const out: Vector[] = [];
    const f = (p: Vector) => dot(h.normal, p) - h.constant;
    for (let i = 0; i < points.length; ++i) {
        const p = points[i], q = points[(i + 1) % points.length];
        const fp = f(p), fq = f(q);
        if (fp >= 0) {
            out.push(p);
        }
        if ((fp > 0 && fq < 0) || (fp < 0 && fq > 0)) {
            const t = fp / (fp - fq);
            out.push(Vector.fromArray([
                p.values[0] + t * (q.values[0] - p.values[0]),
                p.values[1] + t * (q.values[1] - p.values[1])
            ]));
        }
    }
    return out;
}

describe('IntrHalfspace2Polygon2', () => {
    const fi = new IntrHalfspace2Polygon2FI();
    // Counterclockwise unit square.
    const square = poly([[0, 0], [1, 0], [1, 1], [0, 1]]);

    it('reports no intersection when the polygon is strictly outside', () => {
        // x >= 5.
        const result = fi.find(halfspace(1, 0, 5), square);
        expect(result.intersect).toBe(false);
        expect(result.polygon).toHaveLength(0);
    });

    it('reports intersection with an empty polygon when no clipping is needed', () => {
        // x >= -1 contains the square; upstream returns intersect = true with
        // an empty 'polygon' because the caller already has the input.
        const result = fi.find(halfspace(1, 0, -1), square);
        expect(result.intersect).toBe(true);
        expect(result.polygon).toHaveLength(0);
    });

    it('treats a polygon touching the boundary from inside as unclipped', () => {
        // x >= 0: vertices at x = 0 have distance 0, so negative == 0.
        const result = fi.find(halfspace(1, 0, 0), square);
        expect(result.intersect).toBe(true);
        expect(result.polygon).toHaveLength(0);
    });

    it('clips the square with a vertical line, positiveIndex > 0', () => {
        // x <= 0.5, i.e. the halfspace -x >= -0.5. Vertices 0 and 3 (x = 0)
        // are positive, vertices 1 and 2 (x = 1) are negative, so the first
        // positive index is 0. Use the opposite halfspace to exercise the
        // positiveIndex > 0 branch instead: x >= 0.5 keeps vertices 1 and 2.
        const result = fi.find(halfspace(1, 0, 0.5), square);
        expect(result.intersect).toBe(true);
        const values = result.polygon.map(p => [p.values[0], p.values[1]]);
        expect(values).toEqual([[0.5, 0], [1, 0], [1, 1], [0.5, 1]]);
        expect(area(result.polygon)).toBeCloseTo(0.5, 12);
    });

    it('clips the square when the kept run wraps to index 0', () => {
        // -x >= -0.5 keeps vertices 0 and 3, whose indices wrap around.
        const result = fi.find(halfspace(-1, 0, -0.5), square);
        expect(result.intersect).toBe(true);
        expect(area(result.polygon)).toBeCloseTo(0.5, 12);
        for (const p of result.polygon) {
            expect(p.values[0]).toBeLessThanOrEqual(0.5 + 1e-15);
        }
        expect(result.polygon).toHaveLength(4);
    });

    it('clips a triangle to a smaller triangle', () => {
        const triangle = poly([[0, 0], [2, 0], [0, 2]]);
        // y >= 1 keeps only the region near (0,2).
        const result = fi.find(halfspace(0, 1, 1), triangle);
        expect(result.intersect).toBe(true);
        expect(area(result.polygon)).toBeCloseTo(0.25 * area(triangle), 12);
    });

    it('clips exactly through two vertices (zero distances)', () => {
        // The diagonal line x + y = 1 passes through (1,0) and (0,1). Build
        // the constant from the normalized normal so the distances at those
        // vertices are exactly zero.
        const n = Vector.fromArray([1, 1]);
        normalize(n);
        const h = Halfspace.fromNormalConstant(n, n.values[0]);
        const result = fi.find(h, square);
        expect(result.intersect).toBe(true);
        // (n,p,z) = (1,1,2): the kept polygon is the triangle (1,0),(1,1),(0,1).
        const values = result.polygon.map(p => [p.values[0], p.values[1]]);
        expect(values).toEqual([[1, 0], [1, 1], [0, 1]]);
        expect(area(result.polygon)).toBeCloseTo(0.5, 12);
    });

    it('agrees with a Sutherland-Hodgman oracle on random convex polygons', () => {
        const rand = makeRandom(4242);
        let numClipped = 0, numOutside = 0, numWhole = 0;
        for (let trial = 0; trial < 400; ++trial) {
            // A convex polygon: a regular n-gon, randomly scaled, rotated and
            // translated, with vertices in counterclockwise order.
            const n = 3 + Math.floor(6 * rand());
            const cx = 2 * rand() - 1, cy = 2 * rand() - 1;
            const radius = 0.5 + rand();
            const phase = 2 * Math.PI * rand();
            const points: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                const a = phase + (2 * Math.PI * i) / n;
                points.push(Vector.fromArray([
                    cx + radius * Math.cos(a), cy + radius * Math.sin(a)]));
            }

            const h = halfspace(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            const distances = points.map(p => dot(h.normal, p) - h.constant);
            const positive = distances.filter(d => d > 0).length;
            const negative = distances.filter(d => d < 0).length;

            const result = fi.find(h, points);
            expect(result.intersect).toBe(positive > 0);

            if (positive === 0) {
                ++numOutside;
                expect(result.polygon).toHaveLength(0);
            } else if (negative === 0) {
                ++numWhole;
                expect(result.polygon).toHaveLength(0);
            } else {
                ++numClipped;
                const expected = clipOracle(points, h);
                expect(Math.abs(area(result.polygon)))
                    .toBeCloseTo(Math.abs(area(expected)), 10);
                for (const p of result.polygon) {
                    expect(dot(h.normal, p) - h.constant).toBeGreaterThan(-1e-12);
                }
            }
        }
        expect(numClipped).toBeGreaterThan(50);
        expect(numOutside).toBeGreaterThan(10);
        expect(numWhole).toBeGreaterThan(10);
    });
});
