import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, dot, normalize } from '../src/Vector.js';
import {
    IntrHalfspace3Triangle3TI,
    IntrHalfspace3Triangle3FI
} from '../src/IntrHalfspace3Triangle3.js';
import { length, sub } from '../src/Vector.js';
import { check, fc } from './helpers/arbitraries.js';

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

describe('IntrHalfspace3Triangle3 verification', () => {
    const ti = new IntrHalfspace3Triangle3TI();
    const fi = new IntrHalfspace3Triangle3FI();

    // Integer normals, constants and vertices. The query uses only the signs
    // of dot(N,V) - c and ratios s_i/(s_i - s_j), all invariant under a
    // positive scaling of (N,c), so an unnormalized integer normal drives
    // exactly the same branches while making every sign test exact. That is
    // what makes the on-plane cases of the (n,p,z) table common draws.
    const latticeHs = fc.tuple(fc.integer({ min: -3, max: 3 }),
        fc.integer({ min: -3, max: 3 }), fc.integer({ min: -3, max: 3 }),
        fc.integer({ min: -4, max: 4 }))
        .filter(([a, b, c]) => a !== 0 || b !== 0 || c !== 0)
        .map(([a, b, c, d]) =>
            Halfspace.fromNormalConstant(Vector.fromArray([a, b, c]), d));
    const latticeTri = fc.array(
        fc.array(fc.integer({ min: -4, max: 4 }), { minLength: 3, maxLength: 3 }),
        { minLength: 3, maxLength: 3 })
        .map(vs => triangle(vs[0], vs[1], vs[2]));

    // An independent Sutherland-Hodgman clip of the triangle by dot(N,X) >= c.
    function clipOracle(t: Triangle, h: Halfspace): Vector[] {
        const f = (p: Vector): number => dot(h.normal, p) - h.constant;
        const out: Vector[] = [];
        for (let i = 0; i < 3; ++i) {
            const p = t.v[i], q = t.v[(i + 1) % 3];
            const fp = f(p), fq = f(q);
            if (fp >= 0) {
                out.push(p);
            }
            if ((fp > 0 && fq < 0) || (fp < 0 && fq > 0)) {
                const s = fp / (fp - fq);
                out.push(Vector.fromArray([
                    p.values[0] + s * (q.values[0] - p.values[0]),
                    p.values[1] + s * (q.values[1] - p.values[1]),
                    p.values[2] + s * (q.values[2] - p.values[2])]));
            }
        }
        return out;
    }

    function dedupe(points: readonly Vector[], tol: number): Vector[] {
        const out: Vector[] = [];
        for (const p of points) {
            if (!out.some(q => length(sub(p, q)) <= tol)) {
                out.push(p);
            }
        }
        return out;
    }

    it('TI and FI agree and the clipped polygon equals the exact clip', () => {
        check(fc.tuple(latticeHs, latticeTri), ([h, t]) => {
            const r = fi.find(h, t);
            const s = t.v.map(v => dot(h.normal, v) - h.constant);
            const numNeg = s.filter(x => x < 0).length;
            expect(ti.test(h, t).intersect).toBe(Math.max(...s) >= 0);
            expect(r.intersect).toBe(numNeg < 3);

            const oracle = dedupe(clipOracle(t, h), 1e-12);
            const got = dedupe(r.point.slice(0, r.numPoints), 1e-12);
            expect(got.length).toBe(oracle.length);
            // Both polygons are built from the same exact integer data by one
            // division per clipped edge, so matching to 1e-12 is exact up to
            // the final rounding.
            for (const p of got) {
                expect(oracle.some(q => length(sub(p, q)) <= 1e-12)).toBe(true);
            }
            for (const q of oracle) {
                expect(got.some(p => length(sub(p, q)) <= 1e-12)).toBe(true);
            }
        });
    });

    it('reported points lie in the triangle and in the closed halfspace', () => {
        check(fc.tuple(latticeHs, latticeTri), ([h, t]) => {
            const r = fi.find(h, t);
            const e1 = sub(t.v[1], t.v[0]);
            const e2 = sub(t.v[2], t.v[0]);
            const d11 = dot(e1, e1), d12 = dot(e1, e2), d22 = dot(e2, e2);
            const det = d11 * d22 - d12 * d12;
            for (let k = 0; k < r.numPoints; ++k) {
                const p = r.point[k];
                for (let i = 0; i < 3; ++i) {
                    expect(Number.isFinite(p.values[i])).toBe(true);
                }
                expect(dot(h.normal, p) - h.constant)
                    .toBeGreaterThanOrEqual(-1e-11 * (1 + Math.abs(h.constant)));
                if (det <= 1e-9) {
                    continue;    // degenerate triangle: no barycentrics
                }
                const w = sub(p, t.v[0]);
                const dw1 = dot(w, e1), dw2 = dot(w, e2);
                const b1 = (d22 * dw1 - d12 * dw2) / det;
                const b2 = (d11 * dw2 - d12 * dw1) / det;
                expect(b1).toBeGreaterThanOrEqual(-1e-9);
                expect(b2).toBeGreaterThanOrEqual(-1e-9);
                expect(b1 + b2).toBeLessThanOrEqual(1 + 1e-9);
            }
        });
    });

    it('the clipped area never exceeds the triangle area', () => {
        check(fc.tuple(latticeHs, latticeTri), ([h, t]) => {
            const r = fi.find(h, t);
            const full = polygonArea(t.v, 3);
            const clipped = r.numPoints >= 3
                ? polygonArea(r.point.slice(0, r.numPoints), r.numPoints) : 0;
            expect(clipped).toBeLessThanOrEqual(full + 1e-9);
            if (r.numPoints === 3 && r.point.length >= 3) {
                // A fully-inside triangle is returned unchanged.
                const s = t.v.map(v => dot(h.normal, v) - h.constant);
                if (s.every(x => x >= 0)) {
                    expect(clipped).toBeCloseTo(full, 9);
                }
            }
        });
    });

    it('reported vertices are copies, not aliases of the triangle', () => {
        const t = triangle([0, 0, 1], [1, 0, 1], [0, 1, 1]);
        const r = fi.find(halfspace(0, 0, 1, 0), t);
        expect(r.numPoints).toBe(3);
        for (let k = 0; k < 3; ++k) {
            expect(r.point[k]).not.toBe(t.v[k]);
        }
        r.point[0].set(0, 99);
        expect(t.v[0].values[0]).toBe(0);
    });

    it('the (1,0,2) edge case reports exactly the two on-plane vertices', () => {
        const t = triangle([0, 0, 0], [1, 0, 0], [0, 0, -1]);
        const r = fi.find(halfspace(0, 0, 1, 0), t);
        expect(r.intersect).toBe(true);
        expect(r.numPoints).toBe(2);
        expect(r.point[0].values).toEqual([0, 0, 0]);
        expect(r.point[1].values).toEqual([1, 0, 0]);
    });
});
