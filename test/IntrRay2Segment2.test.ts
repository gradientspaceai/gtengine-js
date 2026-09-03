import { describe, it, expect } from 'vitest';
import { Ray } from '../src/Ray.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, mul, normalize, sub, length } from '../src/Vector.js';
import {
    IntrRay2Segment2TI,
    IntrRay2Segment2FI
} from '../src/IntrRay2Segment2.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function ray(p: number[], d: number[]): Ray {
    const dir = vec(d);
    normalize(dir);
    return Ray.fromOriginDirection(vec(p), dir);
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(vec(p0), vec(p1));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrRay2Segment2', () => {
    const ti = new IntrRay2Segment2TI();
    const fi = new IntrRay2Segment2FI();

    it('finds the single crossing of a transversal ray and segment', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([2, -1], [2, 1]);
        const tiResult = ti.test(r, s);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(1);

        const result = fi.find(r, s);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.rayParameter[0]).toBeCloseTo(2, 12);
        // The segment parameter is relative to the segment center.
        expect(result.segmentParameter[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[0]).toBeCloseTo(2, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('rejects a crossing beyond the segment endpoints', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([2, 2], [2, 3]);
        expect(ti.test(r, s).intersect).toBe(false);
        expect(fi.find(r, s).intersect).toBe(false);
    });

    it('rejects a crossing behind the ray origin', () => {
        const r = ray([0, 0], [-1, 0]);
        const s = seg([2, -1], [2, 1]);
        expect(ti.test(r, s).intersect).toBe(false);
        expect(fi.find(r, s).intersect).toBe(false);
    });

    it('reports no intersection for a parallel, non-collinear segment', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([1, 1], [5, 1]);
        const tiResult = ti.test(r, s);
        expect(tiResult.intersect).toBe(false);
        expect(tiResult.numIntersections).toBe(0);
        expect(fi.find(r, s).intersect).toBe(false);
    });

    it('reports the overlap of a collinear segment', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([2, 0], [5, 0]);
        const tiResult = ti.test(r, s);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(2);

        const result = fi.find(r, s);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.rayParameter[0]).toBeCloseTo(2, 12);
        expect(result.rayParameter[1]).toBeCloseTo(5, 12);
        expect(result.segmentParameter[0]).toBeCloseTo(-1.5, 12);
        expect(result.segmentParameter[1]).toBeCloseTo(1.5, 12);
        expect(result.point[0].values[0]).toBeCloseTo(2, 12);
        expect(result.point[1].values[0]).toBeCloseTo(5, 12);
    });

    it('clips a collinear segment that straddles the ray origin', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([-2, 0], [4, 0]);
        const result = fi.find(r, s);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.rayParameter[0]).toBeCloseTo(0, 12);
        expect(result.rayParameter[1]).toBeCloseTo(4, 12);
        expect(result.segmentParameter[0]).toBeCloseTo(-1, 12);
        expect(result.segmentParameter[1]).toBeCloseTo(3, 12);
    });

    it('reports no intersection for a collinear segment behind the ray', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([-5, 0], [-2, 0]);
        const tiResult = ti.test(r, s);
        expect(tiResult.intersect).toBe(false);
        expect(tiResult.numIntersections).toBe(0);
        expect(fi.find(r, s).intersect).toBe(false);
    });

    it('reports a single touching point for a collinear segment that ends at the ray origin', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([-3, 0], [0, 0]);
        const tiResult = ti.test(r, s);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(1);

        const result = fi.find(r, s);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.rayParameter[0]).toBeCloseTo(0, 12);
        expect(result.segmentParameter[0]).toBeCloseTo(1.5, 12);
        expect(result.point[0].values).toEqual([0, 0]);
    });

    it('treats a degenerate zero-length segment as its center point', () => {
        // The centered form of a zero-length segment has a zero direction
        // vector, so the line-line query classifies the two lines as "the
        // same" and only the ray-parameter interval test remains. A
        // degenerate segment on the ray is reported at its parameter.
        const r = ray([0, 0], [1, 0]);
        const onRay = seg([3, 0], [3, 0]);
        const result = fi.find(r, onRay);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.rayParameter[0]).toBeCloseTo(3, 12);
        expect(result.point[0].values[0]).toBeCloseTo(3, 12);

        // A degenerate segment behind the ray origin is rejected.
        const behind = seg([-3, 0], [-3, 0]);
        expect(fi.find(r, behind).intersect).toBe(false);
        expect(ti.test(r, behind).intersect).toBe(false);
    });

    it('agrees between the test and find queries on random configurations', () => {
        const rand = makeRandom(24680);
        let single = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const r = ray([6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1]);
            const s = seg([6 * rand() - 3, 6 * rand() - 3],
                [6 * rand() - 3, 6 * rand() - 3]);
            const tiResult = ti.test(r, s);
            const fiResult = fi.find(r, s);
            expect(tiResult.intersect).toBe(fiResult.intersect);
            expect(tiResult.numIntersections).toBe(fiResult.numIntersections);

            if (fiResult.numIntersections === 1) {
                ++single;
                expect(fiResult.rayParameter[0]).toBeGreaterThanOrEqual(0);
                const centered = s.getCenteredForm();
                expect(Math.abs(fiResult.segmentParameter[0]))
                    .toBeLessThanOrEqual(centered.extent + 1e-12);
                const p0 = add(r.origin,
                    mul(fiResult.rayParameter[0], r.direction));
                const p1 = add(centered.center,
                    mul(fiResult.segmentParameter[0], centered.direction));
                expect(length(sub(p0, p1))).toBeLessThan(1e-9);
                expect(length(sub(p0, fiResult.point[0]))).toBeLessThan(1e-9);
            }
        }
        expect(single).toBeGreaterThan(20);
    });
});
