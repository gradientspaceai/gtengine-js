import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane.js';
import {
    IntrPlane3Triangle3TI,
    IntrPlane3Triangle3FI
} from '../src/IntrPlane3Triangle3.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, dot, normalize } from '../src/Vector.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function plane(normal: number[], constant: number): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalConstant(n, constant);
}

function triangle(v0: number[], v1: number[], v2: number[]): Triangle {
    return Triangle.fromVertices(Vector.fromArray(v0), Vector.fromArray(v1),
        Vector.fromArray(v2));
}

describe('IntrPlane3Triangle3', () => {
    const ti = new IntrPlane3Triangle3TI();
    const fi = new IntrPlane3Triangle3FI();

    // The plane z = 0.
    const z0 = plane([0, 0, 1], 0);

    it('reports no intersection when all vertices are strictly on one side', () => {
        const t = triangle([0, 0, 1], [1, 0, 2], [0, 1, 3]);
        expect(ti.test(z0, t).intersect).toBe(false);
        expect(ti.test(z0, t).numIntersections).toBe(0);
        const result = fi.find(z0, t);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);

        const below = triangle([0, 0, -1], [1, 0, -2], [0, 1, -3]);
        expect(ti.test(z0, below).intersect).toBe(false);
        expect(fi.find(z0, below).intersect).toBe(false);
    });

    it('clips two edges when the plane separates one vertex (n=1, p=2)', () => {
        const t = triangle([0, 0, -1], [2, 0, 1], [0, 2, 1]);
        const tiResult = ti.test(z0, t);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(2);
        expect(tiResult.isInterior).toBe(true);

        const result = fi.find(z0, t);
        expect(result.numIntersections).toBe(2);
        expect(result.isInterior).toBe(true);
        for (let i = 0; i < 2; ++i) {
            expect(result.point[i].values[2]).toBeCloseTo(0, 12);
        }
        // The clipped points are the midpoints of the two edges through
        // vertex 0.
        const xs = [result.point[0].values[0], result.point[1].values[0]]
            .sort((a, b) => a - b);
        const ys = [result.point[0].values[1], result.point[1].values[1]]
            .sort((a, b) => a - b);
        expect(xs[0]).toBeCloseTo(0, 12);
        expect(xs[1]).toBeCloseTo(1, 12);
        expect(ys[0]).toBeCloseTo(0, 12);
        expect(ys[1]).toBeCloseTo(1, 12);
    });

    it('reports a vertex touch (z=1, one side only)', () => {
        const t = triangle([0, 0, 0], [1, 0, 1], [0, 1, 2]);
        const tiResult = ti.test(z0, t);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(1);
        expect(tiResult.isInterior).toBe(false);

        const result = fi.find(z0, t);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('reports a segment through a vertex when the other two straddle', () => {
        const t = triangle([0, 0, 0], [2, 0, 1], [0, 2, -1]);
        const tiResult = ti.test(z0, t);
        expect(tiResult.numIntersections).toBe(2);
        expect(tiResult.isInterior).toBe(true);

        const result = fi.find(z0, t);
        expect(result.numIntersections).toBe(2);
        expect(result.isInterior).toBe(true);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        // The second point is the midpoint of the edge v1->v2.
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
        expect(result.point[1].values[1]).toBeCloseTo(1, 12);
        expect(result.point[1].values[2]).toBeCloseTo(0, 12);
    });

    it('reports an edge in the plane (z=2 case)', () => {
        const t = triangle([0, 0, 0], [1, 0, 0], [0, 1, 5]);
        const tiResult = ti.test(z0, t);
        expect(tiResult.numIntersections).toBe(2);
        expect(tiResult.isInterior).toBe(false);

        const result = fi.find(z0, t);
        expect(result.numIntersections).toBe(2);
        expect(result.isInterior).toBe(false);
        // Vertex 2 is off the plane, so the edge is <v0,v1>.
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values).toEqual([1, 0, 0]);
    });

    it('reports the whole triangle when it lies in the plane', () => {
        const t = triangle([0, 0, 0], [1, 0, 0], [0, 1, 0]);
        const tiResult = ti.test(z0, t);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(3);

        const result = fi.find(z0, t);
        expect(result.numIntersections).toBe(3);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values).toEqual([1, 0, 0]);
        expect(result.point[2].values).toEqual([0, 1, 0]);
    });

    it('rejects non-3D inputs', () => {
        const p2 = Hyperplane.fromNormalConstant(Vector.fromArray([1, 0]), 0);
        const t2 = Triangle.fromVertices(Vector.fromArray([0, 0]),
            Vector.fromArray([1, 0]), Vector.fromArray([0, 1]));
        expect(() => ti.test(p2, t2)).toThrow();
        expect(() => fi.find(p2, t2)).toThrow();
    });

    it('agrees with the test query and puts points on the plane and triangle', () => {
        let seed = 4242424;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        // Barycentric membership of P in the triangle, expressed by checking
        // that P is a convex combination of the vertices.
        const inTriangle = (t: Triangle, P: Vector): boolean => {
            const e0 = [
                t.v[1].values[0] - t.v[0].values[0],
                t.v[1].values[1] - t.v[0].values[1],
                t.v[1].values[2] - t.v[0].values[2]
            ];
            const e1 = [
                t.v[2].values[0] - t.v[0].values[0],
                t.v[2].values[1] - t.v[0].values[1],
                t.v[2].values[2] - t.v[0].values[2]
            ];
            const d = [
                P.values[0] - t.v[0].values[0],
                P.values[1] - t.v[0].values[1],
                P.values[2] - t.v[0].values[2]
            ];
            const a = e0[0] * e0[0] + e0[1] * e0[1] + e0[2] * e0[2];
            const b = e0[0] * e1[0] + e0[1] * e1[1] + e0[2] * e1[2];
            const c = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];
            const d0 = e0[0] * d[0] + e0[1] * d[1] + e0[2] * d[2];
            const d1 = e1[0] * d[0] + e1[1] * d[1] + e1[2] * d[2];
            const det = a * c - b * b;
            const s = (c * d0 - b * d1) / det;
            const u = (a * d1 - b * d0) / det;
            const eps = 1e-7;
            return s >= -eps && u >= -eps && s + u <= 1 + eps;
        };

        for (let trial = 0; trial < 300; ++trial) {
            const p = plane([rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                rand() * 2 - 1);
            const t = triangle(
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]);

            const tiResult = ti.test(p, t);
            const fiResult = fi.find(p, t);
            expect(fiResult.intersect).toBe(tiResult.intersect);
            expect(fiResult.numIntersections).toBe(tiResult.numIntersections);
            expect(fiResult.isInterior).toBe(tiResult.isInterior);

            for (let i = 0; i < fiResult.numIntersections; ++i) {
                const X = fiResult.point[i];
                expect(dot(p.normal, X) - p.constant).toBeCloseTo(0, 8);
                expect(inTriangle(t, X)).toBe(true);
            }
        }
    });
});
