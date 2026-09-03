import { describe, expect, it } from 'vitest';
import {
    IntrTriangle3CanonicalBox3FI,
    IntrTriangle3CanonicalBox3TI,
    defaultIntrTriangle3CanonicalBox3FIResult,
    defaultIntrTriangle3CanonicalBox3TIResult
} from '../src/IntrTriangle3CanonicalBox3.js';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, length, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { DistTriangle3CanonicalBox3 } from '../src/DistTriangle3CanonicalBox3.js';

const V3 = (x: number, y: number, z: number) => Vector.fromArray([x, y, z]);

const box = (ex: number, ey: number, ez: number) =>
    CanonicalBox.fromExtent(V3(ex, ey, ez));

function tri(a: Vector, b: Vector, c: Vector): Triangle {
    return Triangle.fromVertices(a, b, c);
}

function inBox(p: Vector, b: CanonicalBox, tol: number): boolean {
    for (let i = 0; i < 3; ++i) {
        if (Math.abs(p.values[i]) > b.extent.values[i] + tol) {
            return false;
        }
    }
    return true;
}

// Barycentric containment of a point in the plane of a triangle.
function inTrianglePlane(p: Vector, t: Triangle, tol: number): boolean {
    const e0 = sub(t.v[1], t.v[0]);
    const e1 = sub(t.v[2], t.v[0]);
    const n = cross(e0, e1);
    const nn = dot(n, n);
    if (nn === 0) {
        return true;
    }
    // Off-plane distance.
    if (Math.abs(dot(n, sub(p, t.v[0]))) > tol * Math.sqrt(nn)) {
        return false;
    }
    // Barycentric coordinates via the areas of the three subtriangles.
    const d = sub(p, t.v[0]);
    const d00 = dot(e0, e0), d01 = dot(e0, e1), d11 = dot(e1, e1);
    const d20 = dot(d, e0), d21 = dot(d, e1);
    const denom = d00 * d11 - d01 * d01;
    if (denom === 0) {
        return true;
    }
    const v = (d11 * d20 - d01 * d21) / denom;
    const w = (d00 * d21 - d01 * d20) / denom;
    const u = 1 - v - w;
    return u >= -tol && v >= -tol && w >= -tol;
}

// The area of a planar polygon in 3D.
function polygonArea3(poly: readonly Vector[]): number {
    if (poly.length < 3) {
        return 0;
    }
    let sum = V3(0, 0, 0);
    for (let i = 1; i + 1 < poly.length; ++i) {
        sum = add(sum, cross(sub(poly[i], poly[0]), sub(poly[i + 1], poly[0])));
    }
    return length(sum) / 2;
}

describe('IntrTriangle3CanonicalBox3 default results', () => {
    it('match the upstream default constructors', () => {
        expect(defaultIntrTriangle3CanonicalBox3TIResult())
            .toEqual({ intersect: false });
        const fi = defaultIntrTriangle3CanonicalBox3FIResult();
        expect(fi.insidePolygon).toEqual([]);
        expect(fi.outsidePolygons).toEqual([]);
    });
});

describe('IntrTriangle3CanonicalBox3TI known configurations', () => {
    const ti = new IntrTriangle3CanonicalBox3TI();
    const unitBox = box(1, 1, 1);

    it('reports intersection for a triangle inside the box', () => {
        const t = tri(V3(-0.5, -0.5, 0), V3(0.5, -0.5, 0), V3(0, 0.5, 0));
        expect(ti.test(t, unitBox).intersect).toBe(true);
    });

    it('reports intersection for a triangle that straddles the box', () => {
        const t = tri(V3(-5, 0, 0), V3(5, 0, 0), V3(0, 5, 0));
        expect(ti.test(t, unitBox).intersect).toBe(true);
    });

    it('rejects a triangle separated by a box face normal', () => {
        const t = tri(V3(2, 0, 0), V3(4, 0, 0), V3(3, 1, 0));
        expect(ti.test(t, unitBox).intersect).toBe(false);
    });

    it('rejects a triangle separated by its own plane', () => {
        // The plane x + y + z = 6 misses the box, whose support radius in
        // the direction (1,1,1) is 3.
        const t = tri(V3(6, 0, 0), V3(0, 6, 0), V3(0, 0, 6));
        expect(ti.test(t, unitBox).intersect).toBe(false);
        // Bringing the plane to x + y + z = 3 makes it touch the corner.
        const touching = tri(V3(3, 0, 0), V3(0, 3, 0), V3(0, 0, 3));
        expect(ti.test(touching, unitBox).intersect).toBe(true);
    });

    it('rejects a triangle separated only by an edge-cross axis', () => {
        // A triangle that clears every box face slab and whose own plane
        // meets the box, yet passes outside the corner (1,1,1).
        const t = tri(V3(2, -2, 1.4), V3(-2, 2, 1.4), V3(2, 2, 4));
        expect(ti.test(t, unitBox).intersect).toBe(false);
        // Sliding it down brings it into contact.
        const t2 = tri(V3(2, -2, 1.0), V3(-2, 2, 1.0), V3(2, 2, 4));
        expect(ti.test(t2, unitBox).intersect).toBe(true);
    });

    it('handles a degenerate box (zero extent)', () => {
        const flat = box(1, 1, 0);
        const inPlane = tri(V3(-0.5, -0.5, 0), V3(0.5, -0.5, 0),
            V3(0, 0.5, 0));
        expect(ti.test(inPlane, flat).intersect).toBe(true);
        const above = tri(V3(-0.5, -0.5, 1), V3(0.5, -0.5, 1), V3(0, 0.5, 1));
        expect(ti.test(above, flat).intersect).toBe(false);
    });

    it('handles a degenerate triangle (a segment)', () => {
        const seg = tri(V3(-3, 0, 0), V3(3, 0, 0), V3(0, 0, 0));
        expect(ti.test(seg, unitBox).intersect).toBe(true);
        const segOut = tri(V3(-3, 5, 0), V3(3, 5, 0), V3(0, 5, 0));
        expect(ti.test(segOut, unitBox).intersect).toBe(false);
    });
});

describe('IntrTriangle3CanonicalBox3TI cross-check with the distance query',
    () => {
        it('agrees with DistTriangle3CanonicalBox3 over random inputs', () => {
            const ti = new IntrTriangle3CanonicalBox3TI();
            const dq = new DistTriangle3CanonicalBox3();
            let seed = 31415926;
            const rand = () => {
                seed = (seed * 1103515245 + 12345) & 0x7fffffff;
                return seed / 0x7fffffff;
            };
            const rv = (s: number) => V3((rand() - 0.5) * s,
                (rand() - 0.5) * s, (rand() - 0.5) * s);

            let hits = 0, misses = 0;
            for (let k = 0; k < 2000; ++k) {
                const b = box(0.2 + rand(), 0.2 + rand(), 0.2 + rand());
                const t = tri(rv(4), rv(4), rv(4));
                const d = dq.compute(t, b).distance;
                if (d > 1e-9) {
                    expect(ti.test(t, b).intersect).toBe(false);
                    ++misses;
                } else if (d === 0) {
                    expect(ti.test(t, b).intersect).toBe(true);
                    ++hits;
                }
            }
            expect(hits).toBeGreaterThan(100);
            expect(misses).toBeGreaterThan(100);
        });
    });

describe('IntrTriangle3CanonicalBox3FI known configurations', () => {
    const fi = new IntrTriangle3CanonicalBox3FI();
    const unitBox = box(1, 1, 1);

    it('returns the whole triangle when it is inside the box', () => {
        const t = tri(V3(-0.5, -0.5, 0), V3(0.5, -0.5, 0), V3(0, 0.5, 0));
        const r = fi.find(t, unitBox);
        expect(r.insidePolygon.length).toBe(3);
        expect(polygonArea3(r.insidePolygon)).toBeCloseTo(0.5, 12);
        expect(r.outsidePolygons).toEqual([]);
    });

    it('clips a large triangle to the box cross-section', () => {
        // The triangle lies in the plane z = 0 and covers the whole box
        // cross-section [-1,1] x [-1,1], whose area is 4.
        const t = tri(V3(-10, -10, 0), V3(10, -10, 0), V3(0, 10, 0));
        const r = fi.find(t, unitBox);
        expect(polygonArea3(r.insidePolygon)).toBeCloseTo(4, 10);
        for (const v of r.insidePolygon) {
            expect(inBox(v, unitBox, 1e-12)).toBe(true);
            expect(inTrianglePlane(v, t, 1e-9)).toBe(true);
        }
        // Four of the box faces cut the triangle, so four outside pieces are
        // reported.
        expect(r.outsidePolygons.length).toBe(4);
    });

    it('returns an empty inside polygon for a separated triangle', () => {
        const t = tri(V3(2, 0, 0), V3(4, 0, 0), V3(3, 1, 0));
        const r = fi.find(t, unitBox);
        expect(r.insidePolygon).toEqual([]);
        expect(r.outsidePolygons.length).toBeGreaterThan(0);
        // The single outside polygon is the whole triangle.
        expect(polygonArea3(r.outsidePolygons[r.outsidePolygons.length - 1]))
            .toBeCloseTo(1, 12);
    });

    it('keeps a triangle coplanar with a box face', () => {
        // The CONTAINED configuration: the triangle lies in the plane
        // z = +1, which is a box face.
        const t = tri(V3(-0.5, -0.5, 1), V3(0.5, -0.5, 1), V3(0, 0.5, 1));
        const r = fi.find(t, unitBox);
        expect(r.insidePolygon.length).toBe(3);
        expect(polygonArea3(r.insidePolygon)).toBeCloseTo(0.5, 12);
    });

    it('can produce a clipped polygon with more than seven vertices', () => {
        // The upstream comment claims the intersection has at most 7
        // vertices. The plane of the triangle already cuts the box in a
        // hexagon, and each of the three triangle edges can add another
        // vertex, so the true bound is 9. This configuration produces an
        // octagon.
        const b = CanonicalBox.fromExtent(V3(0.8596554281933491,
            0.9707823279643348, 0.9405479910972285));
        const t = tri(
            V3(-1.996932, 1.970662, -1.980587),
            V3(1.254711, -1.558586, -0.733520),
            V3(0.658056, -0.151536, 1.434068));
        const r = fi.find(t, b);
        expect(r.insidePolygon.length).toBe(8);
        for (const v of r.insidePolygon) {
            expect(inBox(v, b, 1e-12)).toBe(true);
            expect(inTrianglePlane(v, t, 1e-9)).toBe(true);
        }
    });
});

describe('IntrTriangle3CanonicalBox3FI randomized checks', () => {
    it('clips to a polygon inside both the box and the triangle', () => {
        const ti = new IntrTriangle3CanonicalBox3TI();
        const fi = new IntrTriangle3CanonicalBox3FI();
        let seed = 606060;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };
        const rv = (s: number) => V3((rand() - 0.5) * s, (rand() - 0.5) * s,
            (rand() - 0.5) * s);

        let clipped = 0, empty = 0;
        for (let k = 0; k < 1500; ++k) {
            const b = box(0.3 + rand(), 0.3 + rand(), 0.3 + rand());
            const t = tri(rv(4), rv(4), rv(4));
            const r = fi.find(t, b);

            if (r.insidePolygon.length > 0) {
                ++clipped;
                // A nonempty clip implies the objects intersect.
                expect(ti.test(t, b).intersect).toBe(true);
                // The plane of the triangle cuts the box in a polygon of up
                // to 6 edges and the three triangle edges add up to 3 more.
                expect(r.insidePolygon.length).toBeLessThanOrEqual(9);
                for (const v of r.insidePolygon) {
                    expect(inBox(v, b, 1e-9)).toBe(true);
                    expect(inTrianglePlane(v, t, 1e-7)).toBe(true);
                }
                // The clipped polygon cannot be larger than the triangle.
                const triArea = polygonArea3(t.v);
                expect(polygonArea3(r.insidePolygon))
                    .toBeLessThanOrEqual(triArea * (1 + 1e-6) + 1e-9);
            } else {
                ++empty;
            }

            // The inside and outside pieces together account for the whole
            // triangle area.
            let total = polygonArea3(r.insidePolygon);
            for (const poly of r.outsidePolygons) {
                total += polygonArea3(poly);
            }
            expect(total).toBeCloseTo(polygonArea3(t.v), 6);
        }
        expect(clipped).toBeGreaterThan(100);
        expect(empty).toBeGreaterThan(100);
    });
});
