import { describe, it, expect } from 'vitest';
import { Rectangle } from '../src/Rectangle';
import { Segment } from '../src/Segment';
import { Vector, add, dot, mul, sub, normalize } from '../src/Vector';
import {
    IntrSegment3Rectangle3TI,
    IntrSegment3Rectangle3FI
} from '../src/IntrSegment3Rectangle3';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The unit square in the z = 0 plane, extents 1 along x and y.
const square = Rectangle.fromCenterAxisExtent(vec(0, 0, 0),
    [vec(1, 0, 0), vec(0, 1, 0)], Vector.fromArray([1, 1]));

describe('IntrSegment3Rectangle3', () => {
    const ti = new IntrSegment3Rectangle3TI();
    const fi = new IntrSegment3Rectangle3FI();

    it('finds a crossing at the rectangle center', () => {
        const s = segment([0, 0, -2], [0, 0, 2]);
        expect(ti.test(s, square).intersect).toBe(true);
        const result = fi.find(s, square);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(0.5, 12);
        expect(result.rectCoord[0]).toBeCloseTo(0, 12);
        expect(result.rectCoord[1]).toBeCloseTo(0, 12);
        // The third component is unused by IntrLine3Rectangle3 (issue #141).
        expect(result.rectCoord[2]).toBe(0);
        expect(result.point.values[0]).toBeCloseTo(0, 12);
        expect(result.point.values[2]).toBeCloseTo(0, 12);
    });

    it('finds an off-center crossing and reports rectangle coordinates', () => {
        const s = segment([0.5, -0.25, -1], [0.5, -0.25, 3]);
        const result = fi.find(s, square);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(0.25, 12);
        expect(result.rectCoord[0]).toBeCloseTo(0.5, 12);
        expect(result.rectCoord[1]).toBeCloseTo(-0.25, 12);
        expect(result.point.values[0]).toBeCloseTo(0.5, 12);
        expect(result.point.values[1]).toBeCloseTo(-0.25, 12);
        expect(result.point.values[2]).toBeCloseTo(0, 12);
    });

    it('accepts an intersection exactly at a rectangle corner', () => {
        const s = segment([1, 1, -1], [1, 1, 1]);
        expect(ti.test(s, square).intersect).toBe(true);
        const result = fi.find(s, square);
        expect(result.intersect).toBe(true);
        expect(result.rectCoord[0]).toBeCloseTo(1, 12);
        expect(result.rectCoord[1]).toBeCloseTo(1, 12);
    });

    it('accepts an intersection exactly at a segment endpoint', () => {
        // The endpoint (0,0,0) is on the rectangle, so t = 0.
        const s = segment([0, 0, 0], [0, 0, 4]);
        expect(ti.test(s, square).intersect).toBe(true);
        const result = fi.find(s, square);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(0, 12);
    });

    it('rejects an intersection beyond the segment endpoint', () => {
        // The supporting line hits the rectangle at t = 2.
        const s = segment([0, 0, -2], [0, 0, -1]);
        expect(ti.test(s, square).intersect).toBe(false);
        const result = fi.find(s, square);
        expect(result.intersect).toBe(false);
        expect(result.parameter).toBe(0);
        expect(result.point.values).toEqual([0, 0, 0]);
    });

    it('rejects a crossing outside the rectangle extents', () => {
        const s = segment([3, 0, -1], [3, 0, 1]);
        expect(ti.test(s, square).intersect).toBe(false);
        expect(fi.find(s, square).intersect).toBe(false);
    });

    it('reports no intersection when the segment is in the plane', () => {
        // The upstream queries deliberately report no intersection for a
        // segment lying in the plane of the rectangle.
        const s = segment([-2, 0, 0], [2, 0, 0]);
        expect(ti.test(s, square).intersect).toBe(false);
        expect(fi.find(s, square).intersect).toBe(false);
    });

    it('handles a degenerate (zero-length) segment', () => {
        // The direction is (0,0,0), so DdN = 0 and the line query reports no
        // intersection even for a point on the rectangle.
        const onIt = segment([0, 0, 0], [0, 0, 0]);
        expect(ti.test(onIt, square).intersect).toBe(false);
        expect(fi.find(onIt, square).intersect).toBe(false);
    });

    it('handles a rotated, translated rectangle', () => {
        const a0 = vec(1, 1, 0);
        normalize(a0);
        const a1 = vec(-1, 1, 1);
        normalize(a1);
        const rect = Rectangle.fromCenterAxisExtent(vec(2, -1, 3), [a0, a1],
            Vector.fromArray([1.5, 0.75]));
        // A point on the rectangle, then a segment straddling it.
        const target = add(rect.center,
            add(mul(0.6, a0), mul(-0.3, a1)));
        const n = vec(
            a0.values[1] * a1.values[2] - a0.values[2] * a1.values[1],
            a0.values[2] * a1.values[0] - a0.values[0] * a1.values[2],
            a0.values[0] * a1.values[1] - a0.values[1] * a1.values[0]);
        const s = Segment.fromEndpoints(sub(target, mul(2, n)),
            add(target, mul(2, n)));
        expect(ti.test(s, rect).intersect).toBe(true);
        const result = fi.find(s, rect);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(0.5, 10);
        expect(result.rectCoord[0]).toBeCloseTo(0.6, 10);
        expect(result.rectCoord[1]).toBeCloseTo(-0.3, 10);
        const diff = sub(result.point, target);
        expect(Math.sqrt(dot(diff, diff))).toBeLessThan(1e-10);
    });

    it('agrees with an independent plane solve on random configurations', () => {
        const rnd = makeRandom(60221408);
        let tiFiMismatch = 0;
        let referenceMismatch = 0;
        let pointMismatch = 0;
        let hits = 0;

        for (let trial = 0; trial < 400; ++trial) {
            const p0 = vec(4 * rnd() - 2, 4 * rnd() - 2, 4 * rnd() - 2);
            const p1 = vec(4 * rnd() - 2, 4 * rnd() - 2, 4 * rnd() - 2);
            const d = sub(p1, p0);
            if (dot(d, d) < 1e-4 || Math.abs(d.values[2]) < 1e-3) {
                continue;
            }
            const s = Segment.fromEndpoints(p0, p1);

            const tiResult = ti.test(s, square);
            const fiResult = fi.find(s, square);
            if (tiResult.intersect !== fiResult.intersect) {
                ++tiFiMismatch;
            }

            // Reference: solve p0.z + t * d.z = 0 and check the bounds.
            const t = -p0.values[2] / d.values[2];
            const q = add(p0, mul(t, d));
            const reference = (t >= 0 && t <= 1 &&
                Math.abs(q.values[0]) <= 1 && Math.abs(q.values[1]) <= 1);
            if (reference !== fiResult.intersect) {
                ++referenceMismatch;
            }

            if (fiResult.intersect) {
                ++hits;
                if (Math.abs(fiResult.parameter - t) > 1e-9) {
                    ++pointMismatch;
                }
                const diff = sub(fiResult.point, q);
                if (Math.sqrt(dot(diff, diff)) > 1e-9) {
                    ++pointMismatch;
                }
                if (Math.abs(fiResult.rectCoord[0] - q.values[0]) > 1e-9 ||
                    Math.abs(fiResult.rectCoord[1] - q.values[1]) > 1e-9) {
                    ++pointMismatch;
                }
            }
        }

        expect(hits).toBeGreaterThan(20);
        expect([tiFiMismatch, referenceMismatch, pointMismatch])
            .toEqual([0, 0, 0]);
    });
});
