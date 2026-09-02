import { describe, it, expect } from 'vitest';
import { Capsule } from '../src/Capsule';
import {
    IntrRay3Capsule3TI,
    IntrRay3Capsule3FI,
    defaultIntrRay3Capsule3FIResult,
    intrRay3Capsule3FIDoQuery
} from '../src/IntrRay3Capsule3';
import { Ray } from '../src/Ray';
import { Segment } from '../src/Segment';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

function capsule(p0: number[], p1: number[], radius: number): Capsule {
    return Capsule.fromSegmentRadius(
        Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1)),
        radius);
}

// Distance from P to the capsule axis segment minus the radius: negative
// strictly inside the solid capsule.
function capsuleSignedDepth(c: Capsule, P: Vector): number {
    const p0 = c.segment.p[0], p1 = c.segment.p[1];
    const e = sub(p1, p0);
    const d = sub(P, p0);
    const ee = dot(e, e);
    let t = ee > 0 ? dot(d, e) / ee : 0;
    t = Math.max(0, Math.min(1, t));
    const closest = add(p0, mul(t, e));
    const r = sub(P, closest);
    return Math.sqrt(dot(r, r)) - c.radius;
}

describe('IntrRay3Capsule3', () => {
    const ti = new IntrRay3Capsule3TI();
    const fi = new IntrRay3Capsule3FI();

    it('finds both crossings of a capsule along the x-axis', () => {
        // A capsule with axis from (-1,0,0) to (1,0,0) and radius 1 spans
        // x in [-2,2] at y=z=0.
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const result = fi.find(ray([-10, 0, 0], [1, 0, 0]), c);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(8, 9);
        expect(result.parameter[1]).toBeCloseTo(12, 9);
        expect(result.point[0].values[0]).toBeCloseTo(-2, 9);
        expect(result.point[1].values[0]).toBeCloseTo(2, 9);
    });

    it('crosses the cylindrical wall', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const result = fi.find(ray([0, -10, 0], [0, 1, 0]), c);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(9, 9);
        expect(result.parameter[1]).toBeCloseTo(11, 9);
    });

    it('clips at the ray origin when the origin is inside', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const result = fi.find(ray([0, 0, 0], [1, 0, 0]), c);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 9);
    });

    it('rejects a capsule behind the ray', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const r = ray([10, 0, 0], [1, 0, 0]);
        const result = fi.find(r, c);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
        expect(ti.test(r, c).intersect).toBe(false);
    });

    it('reports the test query using the ray-segment distance', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        // The ray passes at distance 1 from the axis: touching counts.
        expect(ti.test(ray([-10, 1, 0], [1, 0, 0]), c).intersect).toBe(true);
        expect(ti.test(ray([-10, 1.0001, 0], [1, 0, 0]), c).intersect)
            .toBe(false);
        // A ray pointing away from the capsule.
        expect(ti.test(ray([10, 0, 0], [1, 0, 0]), c).intersect).toBe(false);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const result = defaultIntrRay3Capsule3FIResult();
        intrRay3Capsule3FIDoQuery(vec(-10, 0, 0), vec(1, 0, 0), c, result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(8, 9);
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('agrees with brute-force sampling on random rays', () => {
        let seed = 5150515;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const c = capsule([-1, 0.5, 0.25], [1.5, -0.5, 1], 0.9);
        for (let trial = 0; trial < 250; ++trial) {
            const r = ray([rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]);
            const result = fi.find(r, c);

            let sampledHit = false;
            for (let k = 0; k <= 1500; ++k) {
                const t = k * 0.01;
                if (capsuleSignedDepth(c,
                    add(r.origin, mul(t, r.direction))) < -1e-6) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(result.intersect).toBe(true);
                expect(ti.test(r, c).intersect).toBe(true);
            }

            if (result.intersect) {
                for (let i = 0; i < 2; ++i) {
                    expect(result.parameter[i]).toBeGreaterThanOrEqual(-1e-9);
                    const P = add(r.origin,
                        mul(result.parameter[i], r.direction));
                    expect(sub(P, result.point[i]).values[0])
                        .toBeCloseTo(0, 9);
                    expect(capsuleSignedDepth(c, P)).toBeLessThan(1e-7);
                }
                expect(result.parameter[0])
                    .toBeLessThanOrEqual(result.parameter[1] + 1e-12);
            }
        }
    });
});
