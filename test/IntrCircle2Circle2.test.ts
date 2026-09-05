import { describe, it, expect } from 'vitest';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector, length, sub } from '../src/Vector.js';
import {
    IntrCircle2Circle2TI,
    IntrCircle2Circle2FI
} from '../src/IntrCircle2Circle2.js';
import { dot } from '../src/Vector.js';
import { check, fc, positive, wellScaled } from './helpers/arbitraries.js';

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

describe('IntrCircle2Circle2 verification', () => {
    const ti = new IntrCircle2Circle2TI();
    const fi = new IntrCircle2Circle2FI();

    const circleArb = fc.tuple(wellScaled(-5, 5), wellScaled(-5, 5),
        positive(4)).map(([x, y, r]) => circle(x, y, r));

    // Distance of the configuration from the two tangency surfaces
    // |U| = |R0-R1| and |U| = R0+R1, relative to the sizes involved. The
    // two-point branch computes t = sqrt(R0^2/|U|^2 - s^2), which loses half
    // its digits as that argument approaches zero, so properties about the
    // intersection points need a transversality margin.
    function transversality(c0: Hypersphere, c1: Hypersphere): number {
        const d = length(sub(c1.center, c0.center));
        const scale = c0.radius + c1.radius + d + 1;
        return Math.min(Math.abs(d - Math.abs(c0.radius - c1.radius)),
            Math.abs(c0.radius + c1.radius - d)) / scale;
    }

    it('reported points lie on both circles', () => {
        check(fc.tuple(circleArb, circleArb), ([c0, c1]) => {
            const r = fi.find(c0, c1);
            if (!r.intersect || r.numIntersections === INT32_MAX) {
                return;
            }
            if (transversality(c0, c1) < 1e-3) {
                return;
            }
            // Each point is built from one square root of a quantity whose
            // relative accuracy is bounded below by the transversality
            // margin, so 1e-7 is a generous but meaningful bound.
            for (let k = 0; k < r.numIntersections; ++k) {
                const p = r.point[k];
                expect(Number.isFinite(p.values[0])).toBe(true);
                expect(Number.isFinite(p.values[1])).toBe(true);
                expect(length(sub(p, c0.center))).toBeCloseTo(c0.radius, 7);
                expect(length(sub(p, c1.center))).toBeCloseTo(c1.radius, 7);
            }
        });
    });

    it('numIntersections matches the exact |R0-R1| <= |U| <= R0+R1 classification', () => {
        check(fc.tuple(circleArb, circleArb), ([c0, c1]) => {
            const r = fi.find(c0, c1);
            const U = sub(c1.center, c0.center);
            const uu = dot(U, U);
            const dm = c0.radius - c1.radius, dp = c0.radius + c1.radius;
            if (uu === 0 && dm === 0) {
                expect(r.numIntersections).toBe(INT32_MAX);
                expect(r.intersect).toBe(true);
                expect(r.circle.radius).toBe(c0.radius);
                return;
            }
            if (uu < dm * dm || uu > dp * dp) {
                expect(r.intersect).toBe(false);
                expect(r.numIntersections).toBe(0);
            } else {
                expect(r.intersect).toBe(true);
                expect(r.numIntersections === 1 || r.numIntersections === 2)
                    .toBe(true);
            }
        });
    });

    it('the FI point set is symmetric under argument swap', () => {
        check(fc.tuple(circleArb, circleArb), ([c0, c1]) => {
            const a = fi.find(c0, c1);
            const b = fi.find(c1, c0);
            expect(b.intersect).toBe(a.intersect);
            expect(b.numIntersections).toBe(a.numIntersections);
            if (!a.intersect || a.numIntersections === INT32_MAX) {
                return;
            }
            if (transversality(c0, c1) < 1e-3) {
                return;
            }
            // The perpendicular V = Perp(C1-C0) flips sign under the swap, so
            // the two transverse points come back in the opposite order.
            const scale = 1 + length(c0.center) + length(c1.center)
                + c0.radius + c1.radius;
            const tol = 1e-7 * scale;
            for (let k = 0; k < a.numIntersections; ++k) {
                const other = b.point[a.numIntersections - 1 - k];
                expect(length(sub(a.point[k], other)))
                    .toBeLessThanOrEqual(tol);
            }
        });
    });

    it('TI is the solid-disk test and disagrees with FI for nested circles', () => {
        // Upstream's TI query tests |C0-C1| <= R0+R1, which is overlap of the
        // two solid disks; the FI query intersects the two curves. A circle
        // strictly inside another therefore gives TI true and FI zero points.
        const outer = circle(0, 0, 5);
        const inner = circle(0.5, 0, 1);
        expect(ti.test(outer, inner).intersect).toBe(true);
        const r = fi.find(outer, inner);
        expect(r.intersect).toBe(false);
        expect(r.numIntersections).toBe(0);

        // Away from nesting the two queries agree.
        check(fc.tuple(circleArb, circleArb), ([c0, c1]) => {
            const f = fi.find(c0, c1);
            const d = length(sub(c1.center, c0.center));
            if (d < Math.abs(c0.radius - c1.radius)) {
                return;    // nested: the documented disagreement
            }
            expect(ti.test(c0, c1).intersect).toBe(f.intersect);
        });
    });

    it('tangency: touching circles report a single point on both curves', () => {
        // External tangency.
        const e = fi.find(circle(0, 0, 1), circle(3, 0, 2));
        expect(e.numIntersections).toBe(1);
        expect(e.point[0].values[0]).toBeCloseTo(1, 12);
        expect(e.point[0].values[1]).toBeCloseTo(0, 12);
        // Internal tangency (|U| = |R0 - R1|).
        const i = fi.find(circle(0, 0, 3), circle(1, 0, 2));
        expect(i.numIntersections).toBe(1);
        expect(i.point[0].values[0]).toBeCloseTo(3, 12);
        expect(i.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('degenerate circles: zero radius and coincident centers', () => {
        // A zero-radius circle is a point; it meets another circle only when
        // it lies on that circle.
        const onCurve = fi.find(circle(1, 0, 0), circle(0, 0, 1));
        expect(onCurve.intersect).toBe(true);
        expect(onCurve.numIntersections).toBe(1);
        expect(onCurve.point[0].values[0]).toBeCloseTo(1, 12);
        expect(fi.find(circle(2, 0, 0), circle(0, 0, 1)).intersect).toBe(false);

        // Identical circles report the maxInt sentinel and echo circle0.
        const same = fi.find(circle(1, 2, 3), circle(1, 2, 3));
        expect(same.numIntersections).toBe(INT32_MAX);
        expect(same.circle.center.values).toEqual([1, 2]);
        expect(same.circle.radius).toBe(3);

        // Concentric circles of different radii never meet.
        expect(fi.find(circle(0, 0, 1), circle(0, 0, 2)).intersect).toBe(false);
    });

    it('the returned circle is a copy, not an alias of the input', () => {
        const c = circle(1, 2, 3);
        const r = fi.find(c, circle(1, 2, 3));
        expect(r.circle).not.toBe(c);
        r.circle.center.set(0, 99);
        expect(c.center.values[0]).toBe(1);
    });
});
