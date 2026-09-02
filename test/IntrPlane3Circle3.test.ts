import { describe, it, expect } from 'vitest';
import { Circle3 } from '../src/Circle3';
import { Hyperplane } from '../src/Hyperplane';
import {
    IntrPlane3Circle3TI,
    IntrPlane3Circle3FI,
    defaultIntrPlane3Circle3FIResult,
    intrPlane3Circle3InfinitePoints
} from '../src/IntrPlane3Circle3';
import { Vector, dot, length, normalize, sub } from '../src/Vector';

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

function circle(center: number[], normal: number[], radius: number): Circle3 {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Circle3.fromCenterNormalRadius(Vector.fromArray(center), n, radius);
}

const ti = new IntrPlane3Circle3TI();
const fi = new IntrPlane3Circle3FI();

describe('IntrPlane3Circle3', () => {
    it('has an invalid circle in the default result', () => {
        const result = defaultIntrPlane3Circle3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.circle.center.values).toEqual([0, 0, 0]);
        expect(result.circle.normal.values).toEqual([0, 0, 0]);
        expect(result.circle.radius).toBe(0);
    });

    it('finds two crossings of a circle cut by a plane', () => {
        // The unit circle in the xy-plane cut by the plane x = 0.5.
        const C = circle([0, 0, 0], [0, 0, 1], 1);
        const P = plane([1, 0, 0], [0.5, 0, 0]);
        const result = fi.find(P, C);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        const y = Math.sqrt(1 - 0.25);
        for (const point of result.point) {
            expect(point.values[0]).toBeCloseTo(0.5, 12);
            expect(Math.abs(point.values[1])).toBeCloseTo(y, 12);
            expect(point.values[2]).toBeCloseTo(0, 12);
        }
        expect(result.point[0].values[1])
            .not.toBeCloseTo(result.point[1].values[1], 6);
        expect(ti.test(P, C).intersect).toBe(true);
    });

    it('reports a tangent plane as a single point', () => {
        const C = circle([0, 0, 0], [0, 0, 1], 1);
        const P = plane([1, 0, 0], [1, 0, 0]);
        const result = fi.find(P, C);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
        expect(result.point[1].values).toEqual(result.point[0].values);
        expect(ti.test(P, C).intersect).toBe(true);
    });

    it('reports the whole circle when the planes coincide', () => {
        const C = circle([1, 2, 0], [0, 0, 1], 3);
        const P = plane([0, 0, 1], [0, 0, 0]);
        const result = fi.find(P, C);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(intrPlane3Circle3InfinitePoints);
        expect(result.circle.center.values).toEqual([1, 2, 0]);
        expect(result.circle.radius).toBe(3);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(ti.test(P, C).intersect).toBe(true);
    });

    it('reports no intersection for parallel disjoint planes', () => {
        const C = circle([0, 0, 0], [0, 0, 1], 1);
        const P = plane([0, 0, 1], [0, 0, 1]);
        expect(fi.find(P, C).intersect).toBe(false);
        expect(ti.test(P, C).intersect).toBe(false);
    });

    it('reports no intersection when the cut line misses the circle', () => {
        const C = circle([0, 0, 0], [0, 0, 1], 1);
        const P = plane([1, 0, 0], [2, 0, 0]);
        expect(fi.find(P, C).intersect).toBe(false);
        expect(ti.test(P, C).intersect).toBe(false);
    });

    it('agrees with the TI query and geometry on random inputs', () => {
        let state = 31415926;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numTwo = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const C = circle([rand() * 2, rand() * 2, rand() * 2],
                [rand(), rand(), rand() + 0.001],
                0.5 + Math.abs(rand()) * 2);
            const P = plane([rand(), rand(), rand() + 0.001],
                [rand() * 2, rand() * 2, rand() * 2]);
            const result = fi.find(P, C);
            expect(ti.test(P, C).intersect).toBe(result.intersect);

            if (result.numIntersections === 2) {
                ++numTwo;
            }
            if (result.numIntersections === 1
                || result.numIntersections === 2) {
                for (let i = 0; i < 2; ++i) {
                    const X = result.point[i];
                    // The point is on the plane, on the circle plane and at
                    // the circle radius from the circle center.
                    expect(dot(P.normal, X) - P.constant).toBeCloseTo(0, 8);
                    const d = sub(X, C.center);
                    expect(dot(C.normal, d)).toBeCloseTo(0, 8);
                    expect(length(d)).toBeCloseTo(C.radius, 7);
                }
            }
        }
        expect(numTwo).toBeGreaterThan(50);
    });

    it('accepts a circle whose normal is not axis aligned', () => {
        // A circle in the plane x + y + z = 0 cut by the plane z = 0. The
        // intersection points are on the line x + y = 0, z = 0.
        const C = circle([0, 0, 0], [1, 1, 1], 2);
        const P = plane([0, 0, 1], [0, 0, 0]);
        const result = fi.find(P, C);
        expect(result.numIntersections).toBe(2);
        for (const X of result.point) {
            expect(X.values[2]).toBeCloseTo(0, 12);
            expect(X.values[0] + X.values[1]).toBeCloseTo(0, 10);
            expect(length(sub(X, C.center))).toBeCloseTo(2, 10);
        }
    });
});
