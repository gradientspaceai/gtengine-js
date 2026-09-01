import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane';
import { Vector, dot, sub, add, mul, normalize, length } from '../src/Vector';
import { cross } from '../src/Vector3';

function makeRandom(seed: number): () => number {
    let a = seed >>> 0;
    return () => {
        a = (a + 0x6D2B79F5) >>> 0;
        let t = a;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

function randomVector(n: number, rand: () => number): Vector {
    const v = new Vector(n);
    for (let d = 0; d < n; ++d) {
        v.values[d] = 4 * rand() - 2;
    }
    return v;
}

describe('Hyperplane', () => {
    it('default-constructs the plane with normal (0,...,0,1)', () => {
        for (const n of [2, 3, 4]) {
            const plane = new Hyperplane(n);
            expect(plane.dimension).toBe(n);
            const expected = new Array<number>(n).fill(0);
            expected[n - 1] = 1;
            expect(plane.normal.values).toEqual(expected);
            expect(plane.origin.values).toEqual(new Array<number>(n).fill(0));
            expect(plane.constant).toBe(0);
        }
    });

    it('rejects dimensions less than 2', () => {
        expect(() => new Hyperplane(1)).toThrow();
    });

    it('constructs from a normal and a constant', () => {
        const normal = Vector.fromArray([0, 0, 1]);
        const plane = Hyperplane.fromNormalConstant(normal, 3);
        expect(plane.normal.values).toEqual([0, 0, 1]);
        expect(plane.origin.values).toEqual([0, 0, 3]);
        expect(plane.constant).toBe(3);

        // The vector argument is copied.
        normal.values[2] = 100;
        expect(plane.normal.values[2]).toBe(1);

        // The origin is the point on the plane closest to the origin.
        expect(dot(plane.normal, plane.origin)).toBeCloseTo(plane.constant, 12);
    });

    it('constructs from a normal and an origin', () => {
        const normal = Vector.fromArray([0, 1, 0]);
        const origin = Vector.fromArray([5, 2, -7]);
        const plane = Hyperplane.fromNormalOrigin(normal, origin);
        expect(plane.normal.values).toEqual([0, 1, 0]);
        expect(plane.origin.values).toEqual([5, 2, -7]);
        expect(plane.constant).toBe(2);

        normal.values[1] = 100;
        origin.values[0] = 100;
        expect(plane.normal.values[1]).toBe(1);
        expect(plane.origin.values[0]).toBe(5);
    });

    it('rejects mismatched sizes in fromNormalOrigin', () => {
        expect(() => Hyperplane.fromNormalOrigin(
            Vector.fromArray([1, 0]), Vector.fromArray([1, 0, 0]))).toThrow();
    });

    it('constructs a 3D plane from three points', () => {
        const p = [
            Vector.fromArray([1, 0, 0]),
            Vector.fromArray([0, 1, 0]),
            Vector.fromArray([0, 0, 1])
        ];
        const plane = Hyperplane.fromPoints(p);
        const invSqrt3 = 1 / Math.sqrt(3);
        for (let d = 0; d < 3; ++d) {
            expect(plane.normal.values[d]).toBeCloseTo(invSqrt3, 12);
        }
        expect(plane.constant).toBeCloseTo(invSqrt3, 12);

        // Each input point lies on the plane.
        for (const q of p) {
            expect(dot(plane.normal, q) - plane.constant).toBeCloseTo(0, 12);
        }

        // The normal matches UnitCross(p1-p0, p2-p0) including its sign.
        const c = cross(sub(p[1], p[0]), sub(p[2], p[0]));
        normalize(c);
        expect(plane.normal.values).toEqual(c.values);
    });

    it('constructs a 2D line from two points', () => {
        const p = [Vector.fromArray([1, 1]), Vector.fromArray([3, 1])];
        const plane = Hyperplane.fromPoints(p);
        expect(length(plane.normal)).toBeCloseTo(1, 12);
        // The edge is (2,0), so the normal is +/-(0,1).
        expect(Math.abs(plane.normal.values[0])).toBeCloseTo(0, 12);
        expect(Math.abs(plane.normal.values[1])).toBeCloseTo(1, 12);
        for (const q of p) {
            expect(dot(plane.normal, q) - plane.constant).toBeCloseTo(0, 12);
        }
        expect(dot(plane.normal, plane.origin)).toBeCloseTo(plane.constant, 12);
    });

    it('constructs a 4D hyperplane from four points', () => {
        const p = [
            Vector.fromArray([1, 0, 0, 0]),
            Vector.fromArray([0, 1, 0, 0]),
            Vector.fromArray([0, 0, 1, 0]),
            Vector.fromArray([0, 0, 0, 1])
        ];
        const plane = Hyperplane.fromPoints(p);
        expect(length(plane.normal)).toBeCloseTo(1, 10);
        const half = 0.5;
        for (let d = 0; d < 4; ++d) {
            expect(Math.abs(plane.normal.values[d])).toBeCloseTo(half, 10);
        }
        for (const q of p) {
            expect(dot(plane.normal, q) - plane.constant).toBeCloseTo(0, 10);
        }
    });

    it('produces a unit normal orthogonal to all edges (randomized)', () => {
        const rand = makeRandom(0x91A2E);
        let maxError = 0;
        for (let trial = 0; trial < 150; ++trial) {
            const n = 2 + (trial % 3);  // 2, 3, 4
            const p: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                p.push(randomVector(n, rand));
            }
            const plane = Hyperplane.fromPoints(p);

            maxError = Math.max(maxError,
                Math.abs(length(plane.normal) - 1));
            for (let i = 1; i < n; ++i) {
                maxError = Math.max(maxError,
                    Math.abs(dot(plane.normal, sub(p[i], p[0]))));
            }
            // Every input point is on the plane.
            for (let i = 0; i < n; ++i) {
                maxError = Math.max(maxError,
                    Math.abs(dot(plane.normal, p[i]) - plane.constant));
            }
            // origin = constant * normal and Dot(normal, origin) = constant.
            maxError = Math.max(maxError, Math.abs(
                dot(plane.normal, plane.origin) - plane.constant));
        }
        expect(maxError).toBeLessThan(1e-9);
    });

    it('gives signed point-plane distances with the expected sign', () => {
        // The plane z = 2 with normal (0,0,1).
        const plane = Hyperplane.fromNormalConstant(
            Vector.fromArray([0, 0, 1]), 2);
        const above = Vector.fromArray([7, -3, 5]);
        const below = Vector.fromArray([7, -3, -1]);
        const on = Vector.fromArray([7, -3, 2]);
        expect(dot(plane.normal, above) - plane.constant).toBeCloseTo(3, 12);
        expect(dot(plane.normal, below) - plane.constant).toBeCloseTo(-3, 12);
        expect(dot(plane.normal, on) - plane.constant).toBeCloseTo(0, 12);

        // The signed distance of a point measured from the plane origin is
        // the same quantity.
        expect(dot(plane.normal, sub(above, plane.origin)))
            .toBeCloseTo(3, 12);
    });

    it('reproduces the signed distance for random planes and points', () => {
        const rand = makeRandom(0x7E571);
        let maxError = 0;
        for (let trial = 0; trial < 200; ++trial) {
            const normal = randomVector(3, rand);
            if (normalize(normal) === 0) {
                continue;
            }
            const origin = randomVector(3, rand);
            const plane = Hyperplane.fromNormalOrigin(normal, origin);
            const t = 4 * rand() - 2;
            // A point at signed distance t from the plane.
            const inPlane = add(origin, mul(
                cross(normal, Vector.fromArray([1, 2, 3])), rand()));
            const X = add(inPlane, mul(plane.normal, t));
            maxError = Math.max(maxError,
                Math.abs((dot(plane.normal, X) - plane.constant) - t));
        }
        expect(maxError).toBeLessThan(1e-10);
    });

    it('yields a zero normal for degenerate point sets', () => {
        // Collinear points in 3D.
        const collinear = Hyperplane.fromPoints([
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 1, 1]),
            Vector.fromArray([2, 2, 2])
        ]);
        expect(collinear.normal.values).toEqual([0, 0, 0]);
        expect(collinear.constant).toBe(0);
        expect(collinear.origin.values).toEqual([0, 0, 0]);

        // Coincident points in 3D.
        const coincident = Hyperplane.fromPoints([
            Vector.fromArray([1, 2, 3]),
            Vector.fromArray([1, 2, 3]),
            Vector.fromArray([1, 2, 3])
        ]);
        expect(coincident.normal.values).toEqual([0, 0, 0]);

        // Coincident points in 2D.
        const coincident2 = Hyperplane.fromPoints([
            Vector.fromArray([4, 5]),
            Vector.fromArray([4, 5])
        ]);
        expect(coincident2.normal.values).toEqual([0, 0]);
        expect(coincident2.constant).toBe(0);
    });

    it('rejects malformed point arrays', () => {
        expect(() => Hyperplane.fromPoints(
            [Vector.fromArray([1, 2, 3])])).toThrow();
        expect(() => Hyperplane.fromPoints([
            Vector.fromArray([1, 2]),
            Vector.fromArray([1, 2, 3])])).toThrow();
    });

    it('clones and compares for sorted containers', () => {
        const a = Hyperplane.fromNormalConstant(
            Vector.fromArray([0, 0, 1]), 1);
        const b = Hyperplane.fromNormalConstant(
            Vector.fromArray([0, 0, 1]), 2);
        const c = Hyperplane.fromNormalConstant(
            Vector.fromArray([-1, 0, 0]), 1);

        const copy = a.clone();
        expect(copy.equals(a)).toBe(true);
        expect(copy.notEquals(a)).toBe(false);
        copy.constant = 99;
        expect(copy.equals(a)).toBe(false);
        expect(a.constant).toBe(1);

        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(false);
        expect(c.lessThan(a)).toBe(true);
        expect(a.lessThanOrEqual(a.clone())).toBe(true);
        expect(a.greaterThan(c)).toBe(true);
        expect(a.greaterThanOrEqual(a.clone())).toBe(true);
        expect(b.greaterThanOrEqual(a)).toBe(true);
    });
});
