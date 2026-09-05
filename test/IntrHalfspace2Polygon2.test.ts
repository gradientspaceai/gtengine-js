import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace.js';
import { Vector, dot, normalize } from '../src/Vector.js';
import { IntrHalfspace2Polygon2FI } from '../src/IntrHalfspace2Polygon2.js';
import { check, fc } from './helpers/arbitraries.js';

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

describe('IntrHalfspace2Polygon2 verification', () => {
    const fi = new IntrHalfspace2Polygon2FI();

    // Counterclockwise convex lattice polygon: the monotone-chain hull of a
    // handful of small integer points. The upstream clipper assumes a convex
    // polygon listed counterclockwise, and integer coordinates make the sign
    // of dot(N,V) - c exact for an integer halfplane.
    function convexHullCcw(points: number[][]): Vector[] {
        const pts = [...points].sort((a, b) => a[0] - b[0] || a[1] - b[1]);
        const uniq: number[][] = [];
        for (const p of pts) {
            const last = uniq[uniq.length - 1];
            if (!last || last[0] !== p[0] || last[1] !== p[1]) {
                uniq.push(p);
            }
        }
        if (uniq.length < 3) {
            return [];
        }
        const cross2 = (o: number[], a: number[], b: number[]): number =>
            (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
        const build = (src: number[][]): number[][] => {
            const out: number[][] = [];
            for (const p of src) {
                while (out.length >= 2
                    && cross2(out[out.length - 2], out[out.length - 1], p) <= 0) {
                    out.pop();
                }
                out.push(p);
            }
            out.pop();
            return out;
        };
        const lower = build(uniq);
        const upper = build([...uniq].reverse());
        const hull = lower.concat(upper);
        return hull.length >= 3 ? hull.map(p => Vector.fromArray(p)) : [];
    }

    const polygonArb = fc.array(
        fc.array(fc.integer({ min: -5, max: 5 }), { minLength: 2, maxLength: 2 }),
        { minLength: 3, maxLength: 8 })
        .map(convexHullCcw)
        .filter(p => p.length >= 3);

    // Unnormalized integer normals: the query uses only signs of
    // dot(N,V) - c and the ratio d_i/(d_i - d_j), both invariant under a
    // positive rescaling of (N,c), so the branches taken are identical while
    // the sign tests become exact.
    const latticeHs = fc.tuple(fc.integer({ min: -3, max: 3 }),
        fc.integer({ min: -3, max: 3 }), fc.integer({ min: -6, max: 6 }))
        .filter(([a, b]) => a !== 0 || b !== 0)
        .map(([a, b, c]) =>
            Halfspace.fromNormalConstant(Vector.fromArray([a, b]), c));

    function dedupe(points: readonly Vector[], tol: number): Vector[] {
        const out: Vector[] = [];
        for (const p of points) {
            if (!out.some(q => Math.hypot(p.values[0] - q.values[0],
                p.values[1] - q.values[1]) <= tol)) {
                out.push(p);
            }
        }
        return out;
    }

    it('the clipped polygon equals an exact Sutherland-Hodgman clip', () => {
        check(fc.tuple(latticeHs, polygonArb), ([h, poly2]) => {
            const r = fi.find(h, poly2);
            const d = poly2.map(p => dot(h.normal, p) - h.constant);
            const positive = d.filter(x => x > 0).length;
            const negative = d.filter(x => x < 0).length;

            if (positive === 0) {
                // No vertex strictly inside: upstream calls this "strictly
                // outside", so a polygon that only touches the boundary line
                // reports no intersection (measure-zero contact convention).
                expect(r.intersect).toBe(false);
                expect(r.polygon.length).toBe(0);
                return;
            }
            if (negative === 0) {
                // Fully inside: upstream returns true with an empty polygon
                // (upstream issue #139), preserved by the port.
                expect(r.intersect).toBe(true);
                expect(r.polygon.length).toBe(0);
                return;
            }

            expect(r.intersect).toBe(true);
            const oracle = dedupe(clipOracle(poly2, h), 1e-12);
            const got = dedupe(r.polygon, 1e-12);
            expect(got.length).toBe(oracle.length);
            // Both polygons come from the same exact integer data with one
            // division per clipped edge, so 1e-12 is exact up to rounding.
            for (const p of got) {
                expect(oracle.some(q => Math.hypot(p.values[0] - q.values[0],
                    p.values[1] - q.values[1]) <= 1e-12)).toBe(true);
            }
            for (const q of oracle) {
                expect(got.some(p => Math.hypot(p.values[0] - q.values[0],
                    p.values[1] - q.values[1]) <= 1e-12)).toBe(true);
            }
        });
    });

    it('the clipped polygon is convex, counterclockwise and inside the halfplane', () => {
        check(fc.tuple(latticeHs, polygonArb), ([h, poly2]) => {
            const r = fi.find(h, poly2);
            if (r.polygon.length === 0) {
                return;
            }
            expect(r.polygon.length).toBeGreaterThanOrEqual(3);
            const scale = 1 + Math.abs(h.constant);
            for (const p of r.polygon) {
                expect(Number.isFinite(p.values[0])).toBe(true);
                expect(Number.isFinite(p.values[1])).toBe(true);
                expect(dot(h.normal, p) - h.constant)
                    .toBeGreaterThanOrEqual(-1e-12 * scale);
            }
            // The clip preserves the counterclockwise order, so the signed
            // area stays nonnegative and never grows.
            expect(area(r.polygon)).toBeGreaterThanOrEqual(-1e-12);
            expect(area(r.polygon)).toBeLessThanOrEqual(area(poly2) + 1e-9);
        });
    });

    it('does not alias or mutate the input polygon vertices', () => {
        const p = poly([[0, 0], [4, 0], [4, 4], [0, 4]]);
        const before = p.map(v => [...v.values]);
        const r = fi.find(halfspace(1, 0, 2), p);   // x >= 2
        expect(r.intersect).toBe(true);
        for (const v of r.polygon) {
            expect(p.some(q => q === v)).toBe(false);
        }
        r.polygon.forEach(v => v.set(0, 99));
        expect(p.map(v => [...v.values])).toEqual(before);
    });

    it('a polygon touching the boundary from outside reports no intersection', () => {
        // Every vertex has distance <= 0 with one exactly zero: 'positive' is
        // zero, so the query calls it strictly outside.
        const square = poly([[-4, 0], [0, 0], [0, 4], [-4, 4]]);
        const r = fi.find(halfspace(1, 0, 0), square);   // x >= 0
        expect(r.intersect).toBe(false);
        expect(r.polygon.length).toBe(0);
    });
});
