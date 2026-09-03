import { describe, it, expect } from 'vitest';
import { Hypersphere } from '../src/Hypersphere.js';
import {
    IntrRay2Circle2TI,
    IntrRay2Circle2FI,
    defaultIntrRay2Circle2FIResult,
    intrRay2Circle2FIDoQuery
} from '../src/IntrRay2Circle2.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

function circle(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(Vector.fromArray(center), radius);
}

// Squared distance from the disk center minus r^2: negative inside.
function diskDepth(c: Hypersphere, P: Vector): number {
    const d = sub(P, c.center);
    return dot(d, d) - c.radius * c.radius;
}

describe('IntrRay2Circle2', () => {
    const ti = new IntrRay2Circle2TI();
    const fi = new IntrRay2Circle2FI();

    it('finds both crossings when the ray starts outside', () => {
        const result = fi.find(ray([-5, 0], [1, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('clips the near crossing when the ray starts inside the disk', () => {
        const result = fi.find(ray([0, 0], [1, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
    });

    it('rejects a disk behind the ray origin', () => {
        const r = ray([5, 0], [1, 0]);
        const c = circle([0, 0], 1);
        expect(fi.find(r, c).intersect).toBe(false);
        expect(fi.find(r, c).numIntersections).toBe(0);
        expect(ti.test(r, c).intersect).toBe(false);
    });

    it('reports a single touching point when the ray origin is on the circle', () => {
        const result = fi.find(ray([1, 0], [1, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
    });

    it('reports the tangent point', () => {
        const result = fi.find(ray([-5, 1], [1, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(5, 9);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
    });

    it('misses a disk that the line also misses', () => {
        const r = ray([-5, 3], [1, 0]);
        const c = circle([0, 0], 1);
        expect(fi.find(r, c).intersect).toBe(false);
        expect(ti.test(r, c).intersect).toBe(false);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const result = defaultIntrRay2Circle2FIResult();
        intrRay2Circle2FIDoQuery(vec(-5, 0), vec(1, 0), circle([0, 0], 1),
            result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(4, 12);
        expect(result.point[0].values).toEqual([0, 0]);
    });

    it('agrees with brute-force sampling on random rays', () => {
        let seed = 90210;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const c = circle([0.5, -1], 1.75);
        for (let trial = 0; trial < 400; ++trial) {
            const r = ray([rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1]);
            const result = fi.find(r, c);
            expect(ti.test(r, c).intersect).toBe(result.intersect);

            // Sample the ray to detect whether it enters the disk.
            let sampledHit = false;
            for (let k = 0; k <= 4000; ++k) {
                const t = k * 0.005;
                if (diskDepth(c, add(r.origin, mul(t, r.direction))) < -1e-6) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(result.intersect).toBe(true);
            }

            for (let i = 0; i < result.numIntersections; ++i) {
                expect(result.parameter[i]).toBeGreaterThanOrEqual(-1e-12);
                const P = add(r.origin, mul(result.parameter[i], r.direction));
                expect(sub(P, result.point[i]).values[0]).toBeCloseTo(0, 9);
                expect(sub(P, result.point[i]).values[1]).toBeCloseTo(0, 9);
                expect(diskDepth(c, P)).toBeLessThan(1e-8);
            }
        }
    });
});
