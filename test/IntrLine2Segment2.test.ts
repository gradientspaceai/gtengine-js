import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, mul, sub, dot, length } from '../src/Vector.js';
import {
    IntrLine2Segment2TI,
    IntrLine2Segment2FI
} from '../src/IntrLine2Segment2.js';

const INT32_MAX = 2147483647;
const MAX_T = Number.MAX_VALUE;

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function line(px: number, py: number, dx: number, dy: number): Line {
    return Line.fromOriginDirection(v2(px, py), v2(dx, dy));
}

function segment(x0: number, y0: number, x1: number, y1: number): Segment {
    return Segment.fromEndpoints(v2(x0, y0), v2(x1, y1));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrLine2Segment2', () => {
    const ti = new IntrLine2Segment2TI();
    const fi = new IntrLine2Segment2FI();

    it('finds the transverse intersection point and parameters', () => {
        // The x axis crosses the segment from (2,-1) to (2,3) at (2,0).
        const l = line(0, 0, 1, 0);
        const s = segment(2, -1, 2, 3);
        const result = fi.find(l, s);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.lineParameter[0]).toBeCloseTo(2, 12);
        expect(result.lineParameter[1]).toBe(result.lineParameter[0]);
        expect(result.segmentParameter[0]).toBeCloseTo(0.25, 12);
        expect(result.segmentParameter[1]).toBe(result.segmentParameter[0]);
        expect(result.point.values[0]).toBeCloseTo(2, 12);
        expect(result.point.values[1]).toBeCloseTo(0, 12);
        expect(ti.test(l, s).numIntersections).toBe(1);
    });

    it('rejects an intersection outside the segment', () => {
        // The supporting lines cross at (2,0), which is past the endpoint.
        const l = line(0, 0, 1, 0);
        const s = segment(2, 1, 2, 3);
        expect(fi.find(l, s).intersect).toBe(false);
        expect(fi.find(l, s).numIntersections).toBe(0);
        expect(ti.test(l, s).intersect).toBe(false);
    });

    it('accepts a contact at either segment endpoint', () => {
        const l = line(0, 0, 1, 0);
        const r0 = fi.find(l, segment(2, 0, 2, 3));
        expect(r0.intersect).toBe(true);
        expect(r0.segmentParameter[0]).toBe(0);

        const r1 = fi.find(l, segment(2, -3, 2, 0));
        expect(r1.intersect).toBe(true);
        expect(r1.segmentParameter[0]).toBeCloseTo(1, 12);
    });

    it('reports no intersection for parallel but distinct components', () => {
        const l = line(0, 0, 1, 0);
        const s = segment(-1, 1, 1, 1);
        expect(ti.test(l, s).intersect).toBe(false);
        expect(ti.test(l, s).numIntersections).toBe(0);
        expect(fi.find(l, s).numIntersections).toBe(0);
    });

    it('reports the collinear case with the documented parameters', () => {
        const l = line(0, 0, 1, 0);
        const s = segment(-3, 0, 5, 0);
        const tiResult = ti.test(l, s);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(INT32_MAX);

        const result = fi.find(l, s);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(INT32_MAX);
        expect(result.lineParameter).toEqual([-MAX_T, MAX_T]);
        expect(result.segmentParameter).toEqual([0, 1]);
        expect(result.point.values).toEqual([0, 0]);
    });

    it('agrees between TI and FI on random configurations', () => {
        const rnd = makeRandom(4711);
        let numHit = 0, numMiss = 0;
        for (let k = 0; k < 500; ++k) {
            const l = line(rnd() * 6 - 3, rnd() * 6 - 3,
                rnd() * 2 - 1, rnd() * 2 - 1);
            const s = segment(rnd() * 6 - 3, rnd() * 6 - 3,
                rnd() * 6 - 3, rnd() * 6 - 3);
            if (dot(l.direction, l.direction) < 1e-8) {
                continue;
            }
            const tiResult = ti.test(l, s);
            const fiResult = fi.find(l, s);
            expect(tiResult.intersect).toBe(fiResult.intersect);
            expect(tiResult.numIntersections).toBe(fiResult.numIntersections);

            if (fiResult.numIntersections === 1) {
                const pl = add(l.origin,
                    mul(fiResult.lineParameter[0], l.direction));
                const ps = add(s.p[0], mul(fiResult.segmentParameter[0],
                    sub(s.p[1], s.p[0])));
                expect(length(sub(pl, fiResult.point))).toBeLessThan(1e-9);
                expect(length(sub(ps, fiResult.point))).toBeLessThan(1e-9);
                expect(fiResult.segmentParameter[0])
                    .toBeGreaterThanOrEqual(0);
                expect(fiResult.segmentParameter[0]).toBeLessThanOrEqual(1);
                ++numHit;
            } else {
                ++numMiss;
            }
        }
        expect(numHit).toBeGreaterThan(0);
        expect(numMiss).toBeGreaterThan(0);
    });

    it('agrees with the sign test at the segment endpoints', () => {
        const rnd = makeRandom(31415);
        for (let k = 0; k < 400; ++k) {
            const l = line(rnd() * 4 - 2, rnd() * 4 - 2,
                rnd() * 2 - 1, rnd() * 2 - 1);
            const s = segment(rnd() * 4 - 2, rnd() * 4 - 2,
                rnd() * 4 - 2, rnd() * 4 - 2);
            if (dot(l.direction, l.direction) < 1e-8) {
                continue;
            }
            const side = (p: Vector): number => {
                const q = sub(p, l.origin);
                return q.values[0] * l.direction.values[1]
                    - q.values[1] * l.direction.values[0];
            };
            const s0 = side(s.p[0]);
            const s1 = side(s.p[1]);
            if (s0 * s1 < 0) {
                // The endpoints are strictly on opposite sides of the line.
                expect(ti.test(l, s).intersect).toBe(true);
            }
        }
    });

    it('handles a zero-length segment (degenerate upstream behavior)', () => {
        // A degenerate segment is a point, so the supporting "line" of the
        // segment has a zero direction. DotPerp(D0,D1) and DotPerp(Q,D1) are
        // then both zero, and the line-line query classifies the pair as
        // "the same line". The upstream query therefore reports the
        // collinear configuration for a degenerate segment regardless of
        // whether the point is on the line. The port preserves this.
        const l = line(0, 0, 1, 0);
        const onLine = segment(3, 0, 3, 0);
        expect(ti.test(l, onLine).numIntersections).toBe(INT32_MAX);
        expect(fi.find(l, onLine).numIntersections).toBe(INT32_MAX);

        const offLine = segment(3, 1, 3, 1);
        expect(ti.test(l, offLine).intersect).toBe(true);
        expect(ti.test(l, offLine).numIntersections).toBe(INT32_MAX);
        expect(fi.find(l, offLine).numIntersections).toBe(INT32_MAX);
    });
});
