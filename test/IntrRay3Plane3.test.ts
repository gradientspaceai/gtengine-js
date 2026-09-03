import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane';
import {
    IntrRay3Plane3TI,
    IntrRay3Plane3FI,
    defaultIntrRay3Plane3FIResult,
    defaultIntrRay3Plane3TIResult,
    intrRay3Plane3FIDoQuery
} from '../src/IntrRay3Plane3';
import { Ray } from '../src/Ray';
import { Vector, add, dot, mul, normalize } from '../src/Vector';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

const ti = new IntrRay3Plane3TI();
const fi = new IntrRay3Plane3FI();

describe('IntrRay3Plane3', () => {
    it('default-constructs results as no intersection', () => {
        expect(defaultIntrRay3Plane3TIResult()).toEqual({ intersect: false });
        const r = defaultIntrRay3Plane3FIResult();
        expect(r.intersect).toBe(false);
        expect(r.numIntersections).toBe(0);
        expect(r.parameter).toBe(0);
        expect(r.point.values).toEqual([0, 0, 0]);
    });

    it('finds a transverse intersection at a known point', () => {
        const R = ray([0, 0, 5], [0, 0, -1]);
        const P = plane([0, 0, 1], [0, 0, 2]);
        const result = fi.find(R, P);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter).toBeCloseTo(3, 12);
        expect(result.point.values).toEqual([0, 0, 2]);
        expect(ti.test(R, P).intersect).toBe(true);
    });

    it('rejects the intersection behind the ray origin', () => {
        // The line hits the plane at t = -3, which is not on the ray.
        const R = ray([0, 0, 5], [0, 0, 1]);
        const P = plane([0, 0, 1], [0, 0, 2]);
        const result = fi.find(R, P);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(ti.test(R, P).intersect).toBe(false);
    });

    it('accepts an origin exactly on the plane with t = 0', () => {
        const P = plane([0, 0, 1], [0, 0, 2]);
        for (const dir of [[0, 0, 1], [0, 0, -1], [1, 0, 0]]) {
            const R = ray([1, -1, 2], dir);
            const result = fi.find(R, P);
            expect(result.intersect).toBe(true);
            expect(ti.test(R, P).intersect).toBe(true);
            if (dir[2] !== 0) {
                expect(result.numIntersections).toBe(1);
                expect(result.parameter).toBeCloseTo(0, 12);
                expect(result.point.values[0]).toBeCloseTo(1, 12);
                expect(result.point.values[1]).toBeCloseTo(-1, 12);
                expect(result.point.values[2]).toBeCloseTo(2, 12);
            }
            else {
                // The ray is on the plane.
                expect(result.numIntersections).toBe(2147483647);
            }
        }
    });

    it('reports a parallel disjoint ray as no intersection', () => {
        const R = ray([0, 0, 5], [1, 0, 0]);
        const P = plane([0, 0, 1], [0, 0, 2]);
        expect(fi.find(R, P).intersect).toBe(false);
        expect(ti.test(R, P).intersect).toBe(false);
    });

    it('matches upstream TI logic on the sign of the ray direction', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        // Origin above the plane (signedDistance > 0).
        expect(ti.test(ray([0, 0, 3], [0, 0, -1]), P).intersect).toBe(true);
        expect(ti.test(ray([0, 0, 3], [0, 0, 1]), P).intersect).toBe(false);
        // Origin below the plane (signedDistance < 0).
        expect(ti.test(ray([0, 0, -3], [0, 0, 1]), P).intersect).toBe(true);
        expect(ti.test(ray([0, 0, -3], [0, 0, -1]), P).intersect).toBe(false);
    });

    it('exposes the DoQuery helper used by derived queries', () => {
        const P = plane([1, 1, 1], [0, 0, 0]);
        const result = defaultIntrRay3Plane3FIResult();
        intrRay3Plane3FIDoQuery(vec(1, 1, 1), vec(-1, -1, -1), P, result);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(1, 12);

        const behind = defaultIntrRay3Plane3FIResult();
        intrRay3Plane3FIDoQuery(vec(1, 1, 1), vec(1, 1, 1), P, behind);
        expect(behind.intersect).toBe(false);
        expect(behind.numIntersections).toBe(0);
    });

    it('agrees with TI and puts the point on both primitives (randomized)', () => {
        let seed = 987654321;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (): number => 6 * rand() - 3;

        let hits = 0;
        for (let trial = 0; trial < 4000; ++trial) {
            const R = ray([rnd(), rnd(), rnd()], [rnd(), rnd(), rnd()]);
            const P = plane([rnd(), rnd(), rnd()], [rnd(), rnd(), rnd()]);
            const fiResult = fi.find(R, P);
            const tiResult = ti.test(R, P);
            expect(fiResult.intersect).toBe(tiResult.intersect);
            if (fiResult.intersect && fiResult.numIntersections === 1) {
                ++hits;
                expect(fiResult.parameter).toBeGreaterThanOrEqual(0);
                // The point is on the ray.
                const onRay = add(R.origin, mul(fiResult.parameter, R.direction));
                for (let i = 0; i < 3; ++i) {
                    expect(fiResult.point.values[i]).toBeCloseTo(
                        onRay.values[i], 10);
                }
                // The point is on the plane.
                expect(dot(P.normal, fiResult.point) - P.constant)
                    .toBeCloseTo(0, 8);
            }
        }
        expect(hits).toBeGreaterThan(1000);
    });
});
