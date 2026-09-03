import { describe, it, expect } from 'vitest';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import {
    IntrRay3Ellipsoid3TI,
    IntrRay3Ellipsoid3FI,
    defaultIntrRay3Ellipsoid3FIResult,
    intrRay3Ellipsoid3FIDoQuery
} from '../src/IntrRay3Ellipsoid3.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, length, mul, normalize, sub } from '../src/Vector.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

function ellipsoid(center: number[], extent: number[]): Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(Vector.fromArray(center),
        [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)],
        Vector.fromArray(extent));
}

// The value of (X-C)^T M (X-C) - 1: zero on the ellipsoid, negative inside.
function level(E: Hyperellipsoid, X: Vector): number {
    const d = sub(X, E.center);
    let sum = -1;
    for (let i = 0; i < 3; ++i) {
        const t = d.values[i] / E.extent.values[i];
        sum += t * t;
    }
    return sum;
}

const ti = new IntrRay3Ellipsoid3TI();
const fi = new IntrRay3Ellipsoid3FI();

describe('IntrRay3Ellipsoid3', () => {
    it('has an empty default result', () => {
        const result = defaultIntrRay3Ellipsoid3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('finds both crossings when the ray starts outside', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const result = fi.find(ray([-10, 0, 0], [1, 0, 0]), E);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(7, 12);
        expect(result.parameter[1]).toBeCloseTo(13, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-3, 12);
        expect(result.point[1].values[0]).toBeCloseTo(3, 12);
        expect(ti.test(ray([-10, 0, 0], [1, 0, 0]), E).intersect).toBe(true);
    });

    it('clips the near crossing when the ray starts inside', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const R = ray([0, 0, 0], [1, 0, 0]);
        const result = fi.find(R, E);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(3, 12);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(ti.test(R, E).intersect).toBe(true);
    });

    it('reports a grazing ray at the tangent point', () => {
        // The ray y = 2 grazes the ellipsoid at (0,2,0). Rounding decides
        // whether the discriminant is zero or a tiny positive number, so the
        // query reports one or two nearly identical crossings.
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const R = ray([-10, 2, 0], [1, 0, 0]);
        const result = fi.find(R, E);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBeGreaterThanOrEqual(1);
        for (let i = 0; i < result.numIntersections; ++i) {
            expect(result.point[i].values[0]).toBeCloseTo(0, 6);
            expect(result.point[i].values[1]).toBeCloseTo(2, 12);
        }
        expect(ti.test(R, E).intersect).toBe(true);
    });

    it('reports no intersection for a ray pointing away', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const R = ray([-10, 0, 0], [-1, 0, 0]);
        expect(fi.find(R, E).intersect).toBe(false);
        expect(ti.test(R, E).intersect).toBe(false);
        // A ray that misses the ellipsoid entirely.
        const M = ray([-10, 5, 0], [1, 0, 0]);
        expect(fi.find(M, E).intersect).toBe(false);
        expect(ti.test(M, E).intersect).toBe(false);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const result = defaultIntrRay3Ellipsoid3FIResult();
        intrRay3Ellipsoid3FIDoQuery(vec(-10, 0, 0), vec(1, 0, 0),
            ellipsoid([0, 0, 0], [3, 2, 1]), result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(7, 12);
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('agrees with the TI query and geometry on random inputs', () => {
        let state = 112233;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numHits = 0;
        let numInside = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const E = ellipsoid([rand(), rand(), rand()],
                [0.5 + Math.abs(rand()) * 2, 0.5 + Math.abs(rand()) * 2,
                    0.5 + Math.abs(rand()) * 2]);
            const R = ray([rand() * 4, rand() * 4, rand() * 4],
                [rand(), rand(), rand() + 0.001]);
            const result = fi.find(R, E);
            expect(ti.test(R, E).intersect).toBe(result.intersect);

            if (result.intersect) {
                ++numHits;
                const originInside = level(E, R.origin) <= 0;
                if (originInside) {
                    ++numInside;
                    expect(result.parameter[0]).toBe(0);
                }
                for (let i = 0; i < result.numIntersections; ++i) {
                    expect(result.parameter[i]).toBeGreaterThanOrEqual(0);
                    const onRay = add(R.origin,
                        mul(result.parameter[i], R.direction));
                    expect(length(sub(result.point[i], onRay)))
                        .toBeCloseTo(0, 10);
                    // The point is on the ellipsoid unless it is the clipped
                    // ray origin.
                    if (!(i === 0 && result.parameter[0] === 0
                        && originInside)) {
                        expect(level(E, result.point[i])).toBeCloseTo(0, 8);
                    }
                }
            }
        }
        expect(numHits).toBeGreaterThan(15);
        expect(numInside).toBeGreaterThan(0);
    });
});
