import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import {
    IntrTriangle3AlignedBox3TI,
    IntrTriangle3AlignedBox3FI,
    defaultIntrTriangle3AlignedBox3FIResult,
    defaultIntrTriangle3AlignedBox3TIResult,
    intrTriangle3BoxFacePlanes
} from '../src/IntrTriangle3AlignedBox3.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function tri(v0: number[], v1: number[], v2: number[]): Triangle {
    return Triangle.fromVertices(Vector.fromArray(v0), Vector.fromArray(v1),
        Vector.fromArray(v2));
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(Vector.fromArray(min), Vector.fromArray(max));
}

// True when the point is inside the box (with tolerance).
function inBox(B: AlignedBox, X: Vector, tol: number): boolean {
    for (let i = 0; i < 3; ++i) {
        if (X.values[i] < B.min.values[i] - tol
            || X.values[i] > B.max.values[i] + tol) {
            return false;
        }
    }
    return true;
}

// The signed distance from X to the plane of the triangle.
function planeDistance(T: Triangle, X: Vector): number {
    const n = cross(sub(T.v[1], T.v[0]), sub(T.v[2], T.v[0]));
    const len = Math.sqrt(dot(n, n));
    return dot(n, sub(X, T.v[0])) / len;
}

// Barycentric membership test for the triangle (with tolerance).
function inTriangle(T: Triangle, X: Vector, tol: number): boolean {
    const n = cross(sub(T.v[1], T.v[0]), sub(T.v[2], T.v[0]));
    const area2 = dot(n, n);
    for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
        const e = sub(T.v[i1], T.v[i0]);
        const d = sub(X, T.v[i0]);
        if (dot(cross(e, d), n) < -tol * area2) {
            return false;
        }
    }
    return true;
}

const ti = new IntrTriangle3AlignedBox3TI();
const fi = new IntrTriangle3AlignedBox3FI();

const unitBox = box([-1, -1, -1], [1, 1, 1]);

describe('IntrTriangle3AlignedBox3', () => {
    it('default-constructs the results as empty', () => {
        expect(defaultIntrTriangle3AlignedBox3TIResult())
            .toEqual({ intersect: false });
        const r = defaultIntrTriangle3AlignedBox3FIResult();
        expect(r.insidePolygon).toEqual([]);
        expect(r.outsidePolygons).toEqual([]);
    });

    it('keeps a triangle fully inside the box', () => {
        const T = tri([-0.5, -0.5, 0], [0.5, -0.5, 0], [0, 0.5, 0]);
        expect(ti.test(T, unitBox).intersect).toBe(true);
        const result = fi.find(T, unitBox);
        expect(result.insidePolygon.length).toBe(3);
        expect(result.outsidePolygons.length).toBe(0);
        for (let i = 0; i < 3; ++i) {
            for (let j = 0; j < 3; ++j) {
                expect(result.insidePolygon[i].values[j])
                    .toBeCloseTo(T.v[i].values[j], 12);
            }
        }
    });

    it('rejects a triangle entirely outside the box', () => {
        const T = tri([5, 5, 5], [6, 5, 5], [5, 6, 5]);
        expect(ti.test(T, unitBox).intersect).toBe(false);
        const result = fi.find(T, unitBox);
        expect(result.insidePolygon.length).toBe(0);
        expect(result.outsidePolygons.length).toBeGreaterThan(0);
    });

    it('clips a large triangle to the box cross section', () => {
        // A big triangle in the plane z = 0 that covers the whole box slice.
        const T = tri([-10, -10, 0], [10, -10, 0], [0, 10, 0]);
        expect(ti.test(T, unitBox).intersect).toBe(true);
        const result = fi.find(T, unitBox);
        // The clipped polygon is the square [-1,1]x[-1,1] at z = 0.
        expect(result.insidePolygon.length).toBe(4);
        for (const p of result.insidePolygon) {
            expect(inBox(unitBox, p, 1e-12)).toBe(true);
            expect(p.values[2]).toBeCloseTo(0, 12);
            expect(Math.max(Math.abs(p.values[0]), Math.abs(p.values[1])))
                .toBeCloseTo(1, 12);
        }
        expect(result.outsidePolygons.length).toBeGreaterThan(0);
    });

    it('handles a triangle coplanar with a box face', () => {
        // A triangle on the plane z = 1, the top face of the box. Upstream
        // treats the CONTAINED configuration as inside the box.
        const T = tri([-0.5, -0.5, 1], [0.5, -0.5, 1], [0, 0.5, 1]);
        expect(ti.test(T, unitBox).intersect).toBe(true);
        const result = fi.find(T, unitBox);
        expect(result.insidePolygon.length).toBe(3);
        for (const p of result.insidePolygon) {
            expect(p.values[2]).toBeCloseTo(1, 12);
        }
    });

    it('reports a triangle touching a face from outside', () => {
        // The triangle lies on z = 1 but is entirely outside in x, so the
        // clipping by the x face removes it.
        const T = tri([2, -0.5, 1], [3, -0.5, 1], [2.5, 0.5, 1]);
        expect(ti.test(T, unitBox).intersect).toBe(false);
        expect(fi.find(T, unitBox).insidePolygon.length).toBe(0);
    });

    it('detects a face-diagonal separating-axis configuration', () => {
        // The classic case that needs an edge-cross axis: a triangle cutting
        // a corner of the box.
        const near = tri([1.4, 0, 0], [0, 1.4, 0], [0, 0, 1.4]);
        expect(ti.test(near, unitBox).intersect).toBe(true);
        const far = tri([4, 0, 0], [0, 4, 0], [0, 0, 4]);
        expect(ti.test(far, unitBox).intersect).toBe(false);
    });

    it('builds inward-pointing face planes', () => {
        const planes = intrTriangle3BoxFacePlanes(vec(0, 0, 0),
            [Vector.unit(3, 0), Vector.unit(3, 1), Vector.unit(3, 2)],
            vec(1, 2, 3));
        expect(planes.length).toBe(6);
        // The box center has positive height above every face plane.
        for (const plane of planes) {
            expect(dot(plane.normal, vec(0, 0, 0)) - plane.constant)
                .toBeGreaterThan(0);
        }
        // The plane with the -axis[i] normal is the face at +extent[i] and
        // vice versa, as upstream constructs them.
        expect(dot(planes[0].normal, vec(1, 0, 0)) - planes[0].constant)
            .toBeCloseTo(0, 12);
        expect(dot(planes[3].normal, vec(-1, 0, 0)) - planes[3].constant)
            .toBeCloseTo(0, 12);
        expect(dot(planes[1].normal, vec(0, 2, 0)) - planes[1].constant)
            .toBeCloseTo(0, 12);
        expect(dot(planes[4].normal, vec(0, -2, 0)) - planes[4].constant)
            .toBeCloseTo(0, 12);
        expect(dot(planes[2].normal, vec(0, 0, 3)) - planes[2].constant)
            .toBeCloseTo(0, 12);
        expect(dot(planes[5].normal, vec(0, 0, -3)) - planes[5].constant)
            .toBeCloseTo(0, 12);
    });

    it('agrees with TI and produces valid geometry (randomized)', () => {
        let seed = 20240719;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (a: number, b: number): number => a + (b - a) * rand();

        const B = box([-1, -0.5, -2], [1.5, 2, 1]);
        let numClipped = 0;
        let numInside = 0;
        for (let trial = 0; trial < 1500; ++trial) {
            const T = tri(
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)],
                [rnd(-3, 3), rnd(-3, 3), rnd(-3, 3)]);
            const tiResult = ti.test(T, B).intersect;
            const fiResult = fi.find(T, B);
            const hasInside = fiResult.insidePolygon.length > 0;

            if (hasInside) {
                ++numInside;
                // The TI query must agree that the primitives intersect.
                expect(tiResult).toBe(true);
                if (fiResult.insidePolygon.length > 3) {
                    ++numClipped;
                }
                // Every vertex of the inside polygon is in the box and on the
                // triangle.
                for (const p of fiResult.insidePolygon) {
                    expect(inBox(B, p, 1e-9)).toBe(true);
                    expect(Math.abs(planeDistance(T, p))).toBeLessThan(1e-7);
                    expect(inTriangle(T, p, 1e-7)).toBe(true);
                }
                // The centroid is on the triangle and in the box too.
                let centroid = Vector.zero(3);
                for (const p of fiResult.insidePolygon) {
                    centroid = add(centroid, p);
                }
                centroid = mul(1 / fiResult.insidePolygon.length, centroid);
                expect(inBox(B, centroid, 1e-9)).toBe(true);
                expect(inTriangle(T, centroid, 1e-7)).toBe(true);
            }

            // Every outside polygon vertex is on the triangle plane.
            for (const poly of fiResult.outsidePolygons) {
                for (const p of poly) {
                    expect(Math.abs(planeDistance(T, p))).toBeLessThan(1e-7);
                }
            }
        }
        expect(numInside).toBeGreaterThan(100);
        expect(numClipped).toBeGreaterThan(50);
    });
});
