import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane';
import { Hypersphere } from '../src/Hypersphere';
import {
    IntrPlane3Sphere3TI,
    IntrPlane3Sphere3FI,
    defaultIntrPlane3Sphere3FIResult
} from '../src/IntrPlane3Sphere3';
import { Vector, dot, normalize, sub } from '../src/Vector';

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

function sphere(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(Vector.fromArray(center), radius);
}

const ti = new IntrPlane3Sphere3TI();
const fi = new IntrPlane3Sphere3FI();

describe('IntrPlane3Sphere3', () => {
    it('has an all-zero invalid circle in the default result', () => {
        const result = defaultIntrPlane3Sphere3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.isCircle).toBe(false);
        expect(result.circle.center.values).toEqual([0, 0, 0]);
        expect(result.circle.normal.values).toEqual([0, 0, 0]);
        expect(result.circle.radius).toBe(0);
    });

    it('cuts a sphere in a circle of the expected radius', () => {
        // The plane z = 3 cuts the sphere of radius 5 centered at the origin
        // in a circle of radius 4.
        const P = plane([0, 0, 1], [0, 0, 3]);
        const S = sphere([0, 0, 0], 5);
        const result = fi.find(P, S);
        expect(result.intersect).toBe(true);
        expect(result.isCircle).toBe(true);
        expect(result.circle.radius).toBeCloseTo(4, 12);
        expect(result.circle.center.values[2]).toBeCloseTo(3, 12);
        expect(result.circle.normal.values).toEqual([0, 0, 1]);
        expect(ti.test(P, S).intersect).toBe(true);
    });

    it('reports a tangent plane as a single point', () => {
        const P = plane([0, 0, 1], [0, 0, 5]);
        const S = sphere([0, 0, 0], 5);
        const result = fi.find(P, S);
        expect(result.intersect).toBe(true);
        expect(result.isCircle).toBe(false);
        expect(result.point.values[0]).toBeCloseTo(0, 12);
        expect(result.point.values[1]).toBeCloseTo(0, 12);
        expect(result.point.values[2]).toBeCloseTo(5, 12);
        expect(ti.test(P, S).intersect).toBe(true);
    });

    it('reports separation when the plane misses the sphere', () => {
        const P = plane([0, 0, 1], [0, 0, 6]);
        const S = sphere([0, 0, 0], 5);
        expect(fi.find(P, S).intersect).toBe(false);
        expect(ti.test(P, S).intersect).toBe(false);
    });

    it('agrees with the TI query and geometry on random inputs', () => {
        let state = 987654321;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numCircles = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const P = plane([rand(), rand(), rand() + 0.001],
                [rand() * 2, rand() * 2, rand() * 2]);
            const S = sphere([rand() * 3, rand() * 3, rand() * 3],
                0.5 + 2 * Math.abs(rand()));
            const result = fi.find(P, S);
            expect(ti.test(P, S).intersect).toBe(result.intersect);

            const signedDistance = dot(P.normal, S.center) - P.constant;
            expect(result.intersect).toBe(Math.abs(signedDistance) <= S.radius);

            if (result.isCircle) {
                ++numCircles;
                // The circle center is on the plane and its radius satisfies
                // r^2 + d^2 = R^2.
                expect(dot(P.normal, result.circle.center) - P.constant)
                    .toBeCloseTo(0, 10);
                expect(result.circle.radius * result.circle.radius
                    + signedDistance * signedDistance)
                    .toBeCloseTo(S.radius * S.radius, 8);
                // A point on the circle is on the sphere.
                const toCenter = sub(result.circle.center, S.center);
                expect(dot(toCenter, toCenter)
                    + result.circle.radius * result.circle.radius)
                    .toBeCloseTo(S.radius * S.radius, 8);
            }
        }
        expect(numCircles).toBeGreaterThan(50);
    });
});
