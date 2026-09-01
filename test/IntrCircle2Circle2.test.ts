import { describe, it, expect } from 'vitest';
import { Hypersphere } from '../src/Hypersphere';
import { Vector, length, sub } from '../src/Vector';
import {
    IntrCircle2Circle2TI,
    IntrCircle2Circle2FI
} from '../src/IntrCircle2Circle2';

function circle(cx: number, cy: number, r: number): Hypersphere {
    return Hypersphere.fromCenterRadius(Vector.fromArray([cx, cy]), r);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

const INT32_MAX = 2147483647;

describe('IntrCircle2Circle2', () => {
    const ti = new IntrCircle2Circle2TI();
    const fi = new IntrCircle2Circle2FI();

    it('finds the two intersection points of transverse circles', () => {
        // Unit circles centered at (0,0) and (1,0) meet at (1/2, +-sqrt(3)/2).
        const result = fi.find(circle(0, 0, 1), circle(1, 0, 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        const half = Math.sqrt(3) / 2;
        const ys = [result.point[0].values[1], result.point[1].values[1]].sort(
            (a, b) => a - b);
        expect(result.point[0].values[0]).toBeCloseTo(0.5, 12);
        expect(result.point[1].values[0]).toBeCloseTo(0.5, 12);
        expect(ys[0]).toBeCloseTo(-half, 12);
        expect(ys[1]).toBeCloseTo(half, 12);
        expect(ti.test(circle(0, 0, 1), circle(1, 0, 1)).intersect).toBe(true);
    });

    it('reports a single point for externally tangent circles', () => {
        // |U| = R0 + R1 = 3, contact at (2,0).
        const c0 = circle(0, 0, 2), c1 = circle(3, 0, 1);
        const result = fi.find(c0, c1);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(2, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
        expect(ti.test(c0, c1).intersect).toBe(true);
    });

    it('reports a single point for internally tangent circles', () => {
        // |U| = |R0 - R1| = 1, contact at (3,0).
        const c0 = circle(0, 0, 3), c1 = circle(1, 0, 2);
        const result = fi.find(c0, c1);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(3, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('reports no intersection for separated circles', () => {
        const c0 = circle(0, 0, 1), c1 = circle(5, 0, 1);
        expect(ti.test(c0, c1).intersect).toBe(false);
        const result = fi.find(c0, c1);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('reports no intersection for nested circles', () => {
        const c0 = circle(0, 0, 5), c1 = circle(0.5, 0, 1);
        const result = fi.find(c0, c1);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        // The TI query only tests |C0-C1| <= R0+R1, so it treats a circle
        // nested in another as intersecting (this is the upstream behavior:
        // TI is a solid-disk test, FI intersects the curves).
        expect(ti.test(c0, c1).intersect).toBe(true);
    });

    it('reports infinitely many intersections for identical circles', () => {
        const c0 = circle(2, -3, 4);
        const result = fi.find(c0, circle(2, -3, 4));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(INT32_MAX);
        expect(result.circle.center.values).toEqual([2, -3]);
        expect(result.circle.radius).toBe(4);
    });

    it('handles concentric circles of different radii and zero radii', () => {
        const result = fi.find(circle(0, 0, 1), circle(0, 0, 2));
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);

        // Two coincident degenerate circles (points) are "the same circle".
        const degenerate = fi.find(circle(1, 1, 0), circle(1, 1, 0));
        expect(degenerate.numIntersections).toBe(INT32_MAX);

        // A degenerate circle on another circle is a tangency.
        const onCircle = fi.find(circle(0, 0, 1), circle(1, 0, 0));
        expect(onCircle.intersect).toBe(true);
        expect(onCircle.numIntersections).toBe(1);
        expect(onCircle.point[0].values[0]).toBeCloseTo(1, 12);
    });

    it('produces points on both circles for random transverse pairs', () => {
        const rand = makeRandom(24680);
        let numTwo = 0, numNone = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const c0 = circle(4 * rand() - 2, 4 * rand() - 2, 0.2 + 2 * rand());
            const c1 = circle(4 * rand() - 2, 4 * rand() - 2, 0.2 + 2 * rand());
            const result = fi.find(c0, c1);

            const d = length(sub(c1.center, c0.center));
            const oracle = Math.abs(c0.radius - c1.radius) <= d &&
                d <= c0.radius + c1.radius;
            expect(result.intersect).toBe(oracle);

            if (result.numIntersections === 2) {
                ++numTwo;
                for (let i = 0; i < 2; ++i) {
                    expect(length(sub(result.point[i], c0.center)))
                        .toBeCloseTo(c0.radius, 8);
                    expect(length(sub(result.point[i], c1.center)))
                        .toBeCloseTo(c1.radius, 8);
                }
            } else if (result.numIntersections === 0) {
                ++numNone;
            }

            // The TI query is the solid-disk test.
            expect(ti.test(c0, c1).intersect).toBe(d <= c0.radius + c1.radius);
        }
        expect(numTwo).toBeGreaterThan(50);
        expect(numNone).toBeGreaterThan(50);
    });
});
