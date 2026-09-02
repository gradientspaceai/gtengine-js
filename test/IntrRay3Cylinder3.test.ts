import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3';
import {
    IntrRay3Cylinder3FI,
    defaultIntrRay3Cylinder3FIResult,
    intrRay3Cylinder3FIDoQuery
} from '../src/IntrRay3Cylinder3';
import { Line } from '../src/Line';
import { Ray } from '../src/Ray';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

function cylinder(origin: number[], direction: number[], radius: number,
    height: number): Cylinder3 {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(Vector.fromArray(origin), d), radius, height);
}

// Negative strictly inside the solid cylinder, zero on the boundary.
function cylinderSignedDepth(c: Cylinder3, P: Vector): number {
    const diff = sub(P, c.axis.origin);
    const z = dot(diff, c.axis.direction);
    const radial = sub(diff, mul(z, c.axis.direction));
    const rDepth = Math.sqrt(dot(radial, radial)) - c.radius;
    const zDepth = Math.abs(z) - 0.5 * c.height;
    return Math.max(rDepth, zDepth);
}

describe('IntrRay3Cylinder3FI', () => {
    const fi = new IntrRay3Cylinder3FI();

    it('finds both crossings of the cylinder wall', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        const result = fi.find(ray([-10, 0, 0], [1, 0, 0]), c);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(9, 9);
        expect(result.parameter[1]).toBeCloseTo(11, 9);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 9);
        expect(result.point[1].values[0]).toBeCloseTo(1, 9);
    });

    it('finds the crossings of the end disks', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        const result = fi.find(ray([0, 0, -10], [0, 0, 1]), c);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(8, 9);
        expect(result.parameter[1]).toBeCloseTo(12, 9);
    });

    it('clips at the ray origin when the origin is inside', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        const result = fi.find(ray([0, 0, 0], [1, 0, 0]), c);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 9);
    });

    it('rejects a cylinder behind the ray', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        const result = fi.find(ray([10, 0, 0], [1, 0, 0]), c);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('rejects a ray whose line misses the cylinder', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        expect(fi.find(ray([-10, 5, 0], [1, 0, 0]), c).intersect).toBe(false);
        // Above the end disk.
        expect(fi.find(ray([-10, 0, 5], [1, 0, 0]), c).intersect).toBe(false);
    });

    it('rejects infinite cylinders', () => {
        const infinite = cylinder([0, 0, 0], [0, 0, 1], 1, -1);
        expect(() => fi.find(ray([-10, 0, 0], [1, 0, 0]), infinite)).toThrow();
    });

    it('exposes the DoQuery helper without computing points', () => {
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        const result = defaultIntrRay3Cylinder3FIResult();
        intrRay3Cylinder3FIDoQuery(vec(-10, 0, 0), vec(1, 0, 0), c, result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(9, 9);
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('agrees with brute-force sampling on random rays', () => {
        let seed = 2718281;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const c = cylinder([0.5, -0.5, 0.25], [0.3, 0.8, -0.5], 1.1, 2.5);
        for (let trial = 0; trial < 250; ++trial) {
            const r = ray([rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]);
            const result = fi.find(r, c);

            let sampledHit = false;
            for (let k = 0; k <= 1500; ++k) {
                const t = k * 0.01;
                if (cylinderSignedDepth(c,
                    add(r.origin, mul(t, r.direction))) < -1e-6) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(result.intersect).toBe(true);
            }

            if (result.intersect) {
                for (let i = 0; i < 2; ++i) {
                    expect(result.parameter[i]).toBeGreaterThanOrEqual(-1e-9);
                    const P = add(r.origin,
                        mul(result.parameter[i], r.direction));
                    expect(sub(P, result.point[i]).values[0])
                        .toBeCloseTo(0, 9);
                    expect(cylinderSignedDepth(c, P)).toBeLessThan(1e-7);
                }
                expect(result.parameter[0])
                    .toBeLessThanOrEqual(result.parameter[1] + 1e-12);
            }
        }
    });
});
