import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane';
import { Vector, dot } from '../src/Vector';
import {
    IntrConvexPolygonHyperplaneConfiguration as Cfg,
    IntrConvexPolygonHyperplaneTI,
    IntrConvexPolygonHyperplaneFI
} from '../src/IntrConvexPolygonHyperplane';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function plane2(nx: number, ny: number, c: number): Hyperplane {
    return Hyperplane.fromNormalConstant(v2(nx, ny), c);
}

function height(hp: Hyperplane, p: Vector): number {
    return dot(hp.normal, p) - hp.constant;
}

// Shoelace area of a 2D polygon (signed).
function area2(polygon: readonly Vector[]): number {
    let sum = 0;
    for (let i = 0; i < polygon.length; ++i) {
        const p = polygon[i];
        const q = polygon[(i + 1) % polygon.length];
        sum += p.values[0] * q.values[1] - q.values[0] * p.values[1];
    }
    return 0.5 * sum;
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The unit square, counterclockwise.
const square = [v2(-1, -1), v2(1, -1), v2(1, 1), v2(-1, 1)];

describe('IntrConvexPolygonHyperplaneTI', () => {
    const query = new IntrConvexPolygonHyperplaneTI();

    it('rejects polygons with fewer than three vertices', () => {
        const hp = plane2(1, 0, 0);
        for (const poly of [[], [v2(0, 0)], [v2(0, 0), v2(1, 0)]]) {
            const result = query.test(poly, hp);
            expect(result.intersect).toBe(false);
            expect(result.configuration).toBe(Cfg.INVALID_POLYGON);
        }
    });

    it('classifies a polygon strictly on the positive side', () => {
        const result = query.test(square, plane2(1, 0, -5));
        expect(result.intersect).toBe(false);
        expect(result.configuration).toBe(Cfg.POSITIVE_SIDE_STRICT);
    });

    it('classifies a polygon strictly on the negative side', () => {
        const result = query.test(square, plane2(1, 0, 5));
        expect(result.intersect).toBe(false);
        expect(result.configuration).toBe(Cfg.NEGATIVE_SIDE_STRICT);
    });

    it('classifies a polygon touching the plane in an edge', () => {
        // The plane x = -1 contains the edge from (-1,-1) to (-1,1).
        let result = query.test(square, plane2(1, 0, -1));
        expect(result.intersect).toBe(true);
        expect(result.configuration).toBe(Cfg.POSITIVE_SIDE_EDGE);

        result = query.test(square, plane2(1, 0, 1));
        expect(result.intersect).toBe(true);
        expect(result.configuration).toBe(Cfg.NEGATIVE_SIDE_EDGE);
    });

    it('classifies a polygon touching the plane in a vertex', () => {
        // The plane x + y = -2 touches the square only at (-1,-1).
        let result = query.test(square, plane2(1, 1, -2));
        expect(result.intersect).toBe(true);
        expect(result.configuration).toBe(Cfg.POSITIVE_SIDE_VERTEX);

        result = query.test(square, plane2(1, 1, 2));
        expect(result.intersect).toBe(true);
        expect(result.configuration).toBe(Cfg.NEGATIVE_SIDE_VERTEX);
    });

    it('classifies a split polygon', () => {
        const result = query.test(square, plane2(1, 0, 0));
        expect(result.intersect).toBe(true);
        expect(result.configuration).toBe(Cfg.SPLIT);
    });

    it('classifies a polygon contained by the hyperplane', () => {
        // A triangle in the plane z = 0 of 3D.
        const triangle = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)];
        const hp = Hyperplane.fromNormalConstant(v3(0, 0, 1), 0);
        const result = query.test(triangle, hp);
        expect(result.intersect).toBe(true);
        expect(result.configuration).toBe(Cfg.CONTAINED);
    });
});

describe('IntrConvexPolygonHyperplaneFI', () => {
    const query = new IntrConvexPolygonHyperplaneFI();

    it('rejects polygons with fewer than three vertices', () => {
        const result = query.find([v2(0, 0), v2(1, 1)], plane2(1, 0, 0));
        expect(result.configuration).toBe(Cfg.INVALID_POLYGON);
        expect(result.intersection).toHaveLength(0);
        expect(result.positivePolygon).toHaveLength(0);
        expect(result.negativePolygon).toHaveLength(0);
    });

    it('returns the whole polygon on the positive side', () => {
        const result = query.find(square, plane2(1, 0, -5));
        expect(result.configuration).toBe(Cfg.POSITIVE_SIDE_STRICT);
        expect(result.intersection).toHaveLength(0);
        expect(result.positivePolygon).toHaveLength(4);
        expect(result.negativePolygon).toHaveLength(0);
        expect(area2(result.positivePolygon)).toBeCloseTo(4, 12);
    });

    it('returns the whole polygon on the negative side', () => {
        const result = query.find(square, plane2(1, 0, 5));
        expect(result.configuration).toBe(Cfg.NEGATIVE_SIDE_STRICT);
        expect(result.intersection).toHaveLength(0);
        expect(result.positivePolygon).toHaveLength(0);
        expect(result.negativePolygon).toHaveLength(4);
    });

    it('reports the touching vertex', () => {
        const result = query.find(square, plane2(1, 1, -2));
        expect(result.configuration).toBe(Cfg.POSITIVE_SIDE_VERTEX);
        expect(result.intersection).toHaveLength(1);
        expect(result.intersection[0].values).toEqual([-1, -1]);
        expect(result.positivePolygon).toHaveLength(4);
    });

    it('reports the touching edge', () => {
        const result = query.find(square, plane2(1, 0, 1));
        expect(result.configuration).toBe(Cfg.NEGATIVE_SIDE_EDGE);
        expect(result.intersection).toHaveLength(2);
        const xs = result.intersection.map(p => p.values[0]);
        expect(xs).toEqual([1, 1]);
        expect(result.negativePolygon).toHaveLength(4);
    });

    it('returns the polygon as the intersection when contained', () => {
        const triangle = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0)];
        const hp = Hyperplane.fromNormalConstant(v3(0, 0, 1), 0);
        const result = query.find(triangle, hp);
        expect(result.configuration).toBe(Cfg.CONTAINED);
        expect(result.intersection).toHaveLength(3);
        expect(result.positivePolygon).toHaveLength(0);
        expect(result.negativePolygon).toHaveLength(0);
    });

    it('splits the square through its center', () => {
        const hp = plane2(1, 0, 0);
        const result = query.find(square, hp);
        expect(result.configuration).toBe(Cfg.SPLIT);

        // The two sub-polygons partition the square.
        expect(Math.abs(area2(result.positivePolygon))).toBeCloseTo(2, 12);
        expect(Math.abs(area2(result.negativePolygon))).toBeCloseTo(2, 12);

        for (const p of result.positivePolygon) {
            expect(height(hp, p)).toBeGreaterThanOrEqual(0);
        }
        for (const p of result.negativePolygon) {
            expect(height(hp, p)).toBeLessThanOrEqual(0);
        }

        // The intersection is the segment x = 0, |y| <= 1.
        expect(result.intersection).toHaveLength(2);
        const ys = result.intersection.map(p => p.values[1]).sort((a, b) => a - b);
        expect(result.intersection[0].values[0]).toBeCloseTo(0, 12);
        expect(result.intersection[1].values[0]).toBeCloseTo(0, 12);
        expect(ys[0]).toBeCloseTo(-1, 12);
        expect(ys[1]).toBeCloseTo(1, 12);
    });

    it('splits with the deeper side negative (the swap path)', () => {
        // The plane x = 0.5. The maximum positive height is 0.5 and the
        // maximum |negative| height is 1.5, so the query negates the heights,
        // clips, and swaps the sub-polygons.
        const hp = plane2(1, 0, 0.5);
        const result = query.find(square, hp);
        expect(result.configuration).toBe(Cfg.SPLIT);
        expect(Math.abs(area2(result.positivePolygon))).toBeCloseTo(1, 12);
        expect(Math.abs(area2(result.negativePolygon))).toBeCloseTo(3, 12);
        for (const p of result.positivePolygon) {
            expect(height(hp, p)).toBeGreaterThanOrEqual(-1e-15);
        }
        for (const p of result.negativePolygon) {
            expect(height(hp, p)).toBeLessThanOrEqual(1e-15);
        }
    });

    it('splits through two vertices of a hexagon', () => {
        // A regular hexagon split by the plane y = 0, which passes through
        // the vertices (1,0) and (-1,0).
        // The coordinates are written out so that the two vertices on the
        // plane have height exactly zero (Math.sin(Math.PI) is 1.2e-16).
        const h = Math.sqrt(3) / 2;
        const hexagon = [
            v2(1, 0), v2(0.5, h), v2(-0.5, h),
            v2(-1, 0), v2(-0.5, -h), v2(0.5, -h)
        ];
        const hp = plane2(0, 1, 0);
        const result = query.find(hexagon, hp);
        expect(result.configuration).toBe(Cfg.SPLIT);
        const total = Math.abs(area2(hexagon));
        expect(Math.abs(area2(result.positivePolygon))).toBeCloseTo(total / 2, 12);
        expect(Math.abs(area2(result.negativePolygon))).toBeCloseTo(total / 2, 12);
        // The vertices on the plane are shared by both sub-polygons.
        expect(result.positivePolygon).toHaveLength(4);
        expect(result.negativePolygon).toHaveLength(4);
        expect(result.intersection).toHaveLength(2);
        for (const p of result.intersection) {
            expect(Math.abs(p.values[1])).toBeLessThan(1e-15);
        }
    });

    it('splits a triangle in 3D', () => {
        const triangle = [v3(0, 0, -1), v3(2, 0, -1), v3(0, 2, 3)];
        const hp = Hyperplane.fromNormalConstant(v3(0, 0, 1), 0);
        const result = query.find(triangle, hp);
        expect(result.configuration).toBe(Cfg.SPLIT);
        for (const p of result.positivePolygon) {
            expect(height(hp, p)).toBeGreaterThanOrEqual(0);
        }
        for (const p of result.negativePolygon) {
            expect(height(hp, p)).toBeLessThanOrEqual(0);
        }
        expect(result.intersection).toHaveLength(2);
        for (const p of result.intersection) {
            expect(Math.abs(p.values[2])).toBeLessThan(1e-15);
        }
    });

    it('does not alias the input vertices', () => {
        const poly = [v2(-1, -1), v2(1, -1), v2(1, 1), v2(-1, 1)];
        const result = query.find(poly, plane2(1, 0, -5));
        result.positivePolygon[0].values[0] = 99;
        expect(poly[0].values[0]).toBe(-1);
    });

    it('partitions random convex polygons consistently (randomized)', () => {
        const rand = makeRandom(987654321);
        let numSplit = 0;
        for (let trial = 0; trial < 300; ++trial) {
            // Random convex polygon: points on a circle at increasing angles.
            const n = 3 + Math.floor(6 * rand());
            const angles: number[] = [];
            for (let i = 0; i < n; ++i) {
                angles.push(2 * Math.PI * rand());
            }
            angles.sort((a, b) => a - b);
            const radius = 0.5 + rand();
            const cx = 2 * rand() - 1, cy = 2 * rand() - 1;
            const poly = angles.map(a =>
                v2(cx + radius * Math.cos(a), cy + radius * Math.sin(a)));

            // Random line.
            const theta = 2 * Math.PI * rand();
            const hp = plane2(Math.cos(theta), Math.sin(theta),
                2 * rand() - 1);

            const ti = new IntrConvexPolygonHyperplaneTI().test(poly, hp);
            const fi = query.find(poly, hp);
            expect(fi.configuration).toBe(ti.configuration);

            if (fi.configuration === Cfg.SPLIT) {
                ++numSplit;
                const total = Math.abs(area2(poly));
                const pos = Math.abs(area2(fi.positivePolygon));
                const neg = Math.abs(area2(fi.negativePolygon));
                expect(pos + neg).toBeCloseTo(total, 9);
                for (const p of fi.positivePolygon) {
                    expect(height(hp, p)).toBeGreaterThan(-1e-12);
                }
                for (const p of fi.negativePolygon) {
                    expect(height(hp, p)).toBeLessThan(1e-12);
                }
                expect(fi.intersection).toHaveLength(2);
                for (const p of fi.intersection) {
                    expect(Math.abs(height(hp, p))).toBeLessThan(1e-12);
                }
            }
            else if (fi.configuration === Cfg.POSITIVE_SIDE_STRICT
                || fi.configuration === Cfg.POSITIVE_SIDE_VERTEX
                || fi.configuration === Cfg.POSITIVE_SIDE_EDGE) {
                expect(fi.positivePolygon).toHaveLength(poly.length);
                expect(fi.negativePolygon).toHaveLength(0);
            }
            else {
                expect(fi.positivePolygon).toHaveLength(0);
                expect(fi.negativePolygon).toHaveLength(poly.length);
            }
        }
        expect(numSplit).toBeGreaterThan(30);
    });
});
