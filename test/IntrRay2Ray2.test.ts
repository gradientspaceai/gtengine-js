import { describe, it, expect } from 'vitest';
import { Ray } from '../src/Ray.js';
import { Vector, add, mul, normalize, sub, length } from '../src/Vector.js';
import { IntrRay2Ray2TI, IntrRay2Ray2FI } from '../src/IntrRay2Ray2.js';

const INT32_MAX = 2147483647;
const MAX_T = Number.MAX_VALUE;

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function ray(p: number[], d: number[]): Ray {
    const dir = vec(d);
    normalize(dir);
    return Ray.fromOriginDirection(vec(p), dir);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrRay2Ray2', () => {
    const ti = new IntrRay2Ray2TI();
    const fi = new IntrRay2Ray2FI();

    it('finds the single crossing of two transversal rays', () => {
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([2, -1], [0, 1]);
        const tiResult = ti.test(r0, r1);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(1);

        const result = fi.find(r0, r1);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.ray0Parameter[0]).toBeCloseTo(2, 12);
        expect(result.ray1Parameter[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(2, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('rejects a crossing that lies behind one of the rays', () => {
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([-2, -1], [0, 1]);
        expect(ti.test(r0, r1).intersect).toBe(false);
        const result = fi.find(r0, r1);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('reports no intersection for parallel distinct rays', () => {
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([0, 1], [1, 0]);
        const tiResult = ti.test(r0, r1);
        expect(tiResult.intersect).toBe(false);
        expect(tiResult.numIntersections).toBe(0);
        expect(fi.find(r0, r1).intersect).toBe(false);
    });

    it('reports a ray of intersection for collinear same-direction rays', () => {
        // ray1 starts ahead of ray0.
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([3, 0], [1, 0]);
        const tiResult = ti.test(r0, r1);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(INT32_MAX);

        const result = fi.find(r0, r1);
        expect(result.numIntersections).toBe(INT32_MAX);
        expect(result.ray0Parameter).toEqual([3, MAX_T]);
        expect(result.ray1Parameter).toEqual([0, MAX_T]);
        expect(result.point[0].values).toEqual([3, 0]);
    });

    it('handles collinear same-direction rays with ray1 behind ray0', () => {
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([-3, 0], [1, 0]);
        const result = fi.find(r0, r1);
        expect(result.numIntersections).toBe(INT32_MAX);
        expect(result.ray0Parameter).toEqual([0, MAX_T]);
        expect(result.ray1Parameter).toEqual([3, MAX_T]);
        expect(result.point[0].values).toEqual([0, 0]);
    });

    it('reports a segment for overlapping opposite collinear rays', () => {
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([4, 0], [-1, 0]);
        const tiResult = ti.test(r0, r1);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(2);

        const result = fi.find(r0, r1);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.ray0Parameter).toEqual([0, 4]);
        expect(result.ray1Parameter).toEqual([0, 4]);
        expect(result.point[0].values).toEqual([0, 0]);
        expect(result.point[1].values).toEqual([4, 0]);
    });

    it('reports no intersection for disjoint opposite collinear rays', () => {
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([-4, 0], [-1, 0]);
        const tiResult = ti.test(r0, r1);
        expect(tiResult.intersect).toBe(false);
        expect(tiResult.numIntersections).toBe(0);
        expect(fi.find(r0, r1).intersect).toBe(false);
    });

    it('handles opposite collinear rays that share an origin', () => {
        // Upstream reports numIntersections = 1 from the test query (t == 0)
        // but 2 from the find query, whose t >= 0 branch produces the
        // degenerate segment [0,0]. The port preserves both behaviors.
        const r0 = ray([0, 0], [1, 0]);
        const r1 = ray([0, 0], [-1, 0]);
        const tiResult = ti.test(r0, r1);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(1);

        const result = fi.find(r0, r1);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.ray0Parameter).toEqual([0, 0]);
        expect(result.ray1Parameter).toEqual([0, 0]);
        expect(result.point[0].values).toEqual([0, 0]);
        expect(result.point[1].values).toEqual([0, 0]);
    });

    it('agrees between the test and find queries on random configurations', () => {
        const rand = makeRandom(13579);
        let single = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const r0 = ray([6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1]);
            const r1 = ray([6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1]);
            const tiResult = ti.test(r0, r1);
            const fiResult = fi.find(r0, r1);
            expect(tiResult.intersect).toBe(fiResult.intersect);
            expect(tiResult.numIntersections).toBe(fiResult.numIntersections);

            if (fiResult.numIntersections === 1) {
                ++single;
                expect(fiResult.ray0Parameter[0]).toBeGreaterThanOrEqual(0);
                expect(fiResult.ray1Parameter[0]).toBeGreaterThanOrEqual(0);
                const p0 = add(r0.origin,
                    mul(fiResult.ray0Parameter[0], r0.direction));
                const p1 = add(r1.origin,
                    mul(fiResult.ray1Parameter[0], r1.direction));
                expect(length(sub(p0, p1))).toBeLessThan(1e-9);
                expect(length(sub(p0, fiResult.point[0]))).toBeLessThan(1e-9);
            }
        }
        expect(single).toBeGreaterThan(50);
    });
});
