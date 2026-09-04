import { describe, it, expect } from 'vitest';
import {
    getContainerCircle2,
    inContainerCircle2,
    mergeContainersCircle2
} from '../src/ContCircle2.js';
import { Hypersphere, type Circle2 } from '../src/Hypersphere.js';
import { Vector, add, div, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, fc, rotationFrame, seededRandom, wellScaledVector
} from './helpers/arbitraries.js';

function v(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function circle(x: number, y: number, radius: number): Circle2 {
    return Hypersphere.fromCenterRadius(v(x, y), radius);
}

describe('getContainerCircle2', () => {
    it('computes the average-center circle of a square', () => {
        const points = [v(-1, -1), v(1, -1), v(1, 1), v(-1, 1)];
        const c = getContainerCircle2(points);
        expect(c.dimension).toBe(2);
        expect(c.center.values).toEqual([0, 0]);
        expect(c.radius).toBeCloseTo(Math.SQRT2, 14);
    });

    it('degenerates to a zero-radius circle for one point', () => {
        const c = getContainerCircle2([v(3, -4)]);
        expect(c.center.values).toEqual([3, -4]);
        expect(c.radius).toBe(0);
    });

    it('uses the average of the points as the center, not the optimal center', () => {
        // Three points clustered near the origin plus one far away: the
        // average center is pulled toward the cluster, so the circle is
        // larger than the minimum-area circle.
        const points = [v(0, 0), v(0, 0), v(0, 0), v(4, 0)];
        const c = getContainerCircle2(points);
        expect(c.center.values[0]).toBeCloseTo(1, 14);
        expect(c.center.values[1]).toBe(0);
        expect(c.radius).toBeCloseTo(3, 14);
        // The minimum-area circle would have radius 2.
        expect(c.radius).toBeGreaterThan(2);
    });

    it('is invariant to the order of the points', () => {
        const points = [v(1, 2), v(-3, 4), v(0, -5), v(6, 1)];
        const a = getContainerCircle2(points);
        const b = getContainerCircle2([points[2], points[0], points[3], points[1]]);
        expect(a.center.values[0]).toBeCloseTo(b.center.values[0], 14);
        expect(a.center.values[1]).toBeCloseTo(b.center.values[1], 14);
        expect(a.radius).toBeCloseTo(b.radius, 14);
    });

    it('contains every input point, and the radius is attained (randomized)', () => {
        let seed = 987654321;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff - 0.5;
        };

        for (let trial = 0; trial < 50; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 12; ++i) {
                points.push(v(20 * rand(), 20 * rand()));
            }
            const c = getContainerCircle2(points);

            // Center is the arithmetic mean.
            const mean = [0, 0];
            for (const p of points) {
                mean[0] += p.values[0];
                mean[1] += p.values[1];
            }
            expect(c.center.values[0]).toBeCloseTo(mean[0] / points.length, 12);
            expect(c.center.values[1]).toBeCloseTo(mean[1] / points.length, 12);

            // Radius is the largest distance from the center.
            let maxDist = 0;
            for (const p of points) {
                maxDist = Math.max(maxDist, length(sub(p, c.center)));
                expect(length(sub(p, c.center))).toBeLessThanOrEqual(
                    c.radius * (1 + 1e-12));
            }
            expect(c.radius).toBeCloseTo(maxDist, 12);
        }
    });

    it('throws on an empty point set', () => {
        expect(() => getContainerCircle2([]))
            .toThrow('getContainerCircle2: no points.');
    });

    it('throws when the points are not 2D', () => {
        expect(() => getContainerCircle2([Vector.fromArray([0, 0, 0])]))
            .toThrow('getContainerCircle2: points must be 2D.');
        expect(() => getContainerCircle2([v(0, 0), Vector.fromArray([1, 1, 1])]))
            .toThrow('getContainerCircle2: points must be 2D.');
    });
});

describe('inContainerCircle2', () => {
    const c = circle(1, 2, 3);

    it('accepts the center and interior points', () => {
        expect(inContainerCircle2(v(1, 2), c)).toBe(true);
        expect(inContainerCircle2(v(3, 3), c)).toBe(true);
    });

    it('accepts boundary points (the boundary is part of the circle)', () => {
        expect(inContainerCircle2(v(4, 2), c)).toBe(true);
        expect(inContainerCircle2(v(1, -1), c)).toBe(true);
    });

    it('rejects exterior points', () => {
        expect(inContainerCircle2(v(4.0001, 2), c)).toBe(false);
        expect(inContainerCircle2(v(5, 6), c)).toBe(false);
    });

    it('handles a zero-radius circle', () => {
        const point = circle(0, 0, 0);
        expect(inContainerCircle2(v(0, 0), point)).toBe(true);
        expect(inContainerCircle2(v(1e-12, 0), point)).toBe(false);
    });

    it('throws when the inputs are not 2D', () => {
        expect(() => inContainerCircle2(Vector.fromArray([0, 0, 0]), c))
            .toThrow('inContainerCircle2: inputs must be 2D.');
    });
});

describe('mergeContainersCircle2', () => {
    it('merges two separated circles of equal radius', () => {
        const merge = mergeContainersCircle2(circle(0, 0, 1), circle(4, 0, 1));
        expect(merge.center.values[0]).toBeCloseTo(2, 14);
        expect(merge.center.values[1]).toBeCloseTo(0, 14);
        expect(merge.radius).toBeCloseTo(3, 14);
    });

    it('merges two separated circles of different radii', () => {
        // Centers 3 apart with radii 1 and 2: the merged circle spans from
        // (-1,0) to (5,0), so its center is (2,0) and its radius is 3.
        const merge = mergeContainersCircle2(circle(0, 0, 1), circle(3, 0, 2));
        expect(merge.center.values[0]).toBeCloseTo(2, 14);
        expect(merge.center.values[1]).toBeCloseTo(0, 14);
        expect(merge.radius).toBeCloseTo(3, 14);
    });

    it('returns the containing circle when one circle contains the other', () => {
        const outer = circle(0, 0, 10);
        const inner = circle(1, 1, 2);
        expect(mergeContainersCircle2(outer, inner).equals(outer)).toBe(true);
        expect(mergeContainersCircle2(inner, outer).equals(outer)).toBe(true);
    });

    it('returns the larger circle for concentric circles', () => {
        const merge = mergeContainersCircle2(circle(2, 3, 1), circle(2, 3, 5));
        expect(merge.center.values).toEqual([2, 3]);
        expect(merge.radius).toBe(5);
    });

    it('is idempotent for identical circles', () => {
        const c = circle(-1, 4, 2.5);
        expect(mergeContainersCircle2(c, c).equals(c)).toBe(true);
    });

    it('does not alias its inputs when one circle contains the other', () => {
        const outer = circle(0, 0, 10);
        const merge = mergeContainersCircle2(outer, circle(0, 0, 1));
        merge.center.values[0] = 100;
        merge.radius = 0;
        expect(outer.center.values[0]).toBe(0);
        expect(outer.radius).toBe(10);
    });

    it('is commutative and contains both inputs (randomized)', () => {
        let seed = 424242;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };

        for (let trial = 0; trial < 200; ++trial) {
            const c0 = circle(10 * (rand() - 0.5), 10 * (rand() - 0.5), 5 * rand());
            const c1 = circle(10 * (rand() - 0.5), 10 * (rand() - 0.5), 5 * rand());
            const m01 = mergeContainersCircle2(c0, c1);
            const m10 = mergeContainersCircle2(c1, c0);

            expect(m01.center.values[0]).toBeCloseTo(m10.center.values[0], 12);
            expect(m01.center.values[1]).toBeCloseTo(m10.center.values[1], 12);
            expect(m01.radius).toBeCloseTo(m10.radius, 12);

            // Each input circle is inside the merged circle.
            for (const c of [c0, c1]) {
                const d = length(sub(c.center, m01.center));
                expect(d + c.radius).toBeLessThanOrEqual(m01.radius * (1 + 1e-12) + 1e-12);
            }
        }
    });

    it('throws when the inputs are not 2D', () => {
        expect(() => mergeContainersCircle2(circle(0, 0, 1), new Hypersphere(3)))
            .toThrow('mergeContainersCircle2: inputs must be 2D.');
    });
});

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks of the port
// against the upstream ContCircle2.h semantics.
// ---------------------------------------------------------------------------

describe('ContCircle2 verification', () => {
    // Upstream builds the average-center circle: C is the mean of the points
    // and r is the largest distance from C to a point. Cross-check both
    // against an independent computation. The mean is a sum of up to 12
    // well-scaled terms, so the relative error is a few ulps.
    it('center is the mean and radius the largest distance', () => {
        check(fc.array(wellScaledVector(2), { minLength: 1, maxLength: 12 }),
            (points: Vector[]) => {
                const circle = getContainerCircle2(points);

                let mean = new Vector(2);
                for (const p of points) { mean = add(mean, p); }
                mean = div(mean, points.length);
                expectClose(circle.center.get(0), mean.get(0), 1e-12, 1e-12);
                expectClose(circle.center.get(1), mean.get(1), 1e-12, 1e-12);

                let maxDist = 0;
                for (const p of points) {
                    maxDist = Math.max(maxDist, length(sub(p, circle.center)));
                }
                expectClose(circle.radius, maxDist, 1e-12, 1e-12);

                // Containment is exact: the radius is sqrt of the largest
                // squared distance, computed the same way as inContainer.
                for (const p of points) {
                    expect(inContainerCircle2(p, circle)).toBe(true);
                }

                // Minimality for this center: any smaller radius excludes a
                // point (the maximum is attained).
                if (circle.radius > 0) {
                    const tight = Hypersphere.fromCenterRadius(
                        circle.center, circle.radius * (1 - 1e-12));
                    expect(points.some(p => !inContainerCircle2(p, tight)))
                        .toBe(true);
                }
            });
    });

    // inContainer is |P - C| <= r; cross-check against the squared form.
    it('inContainer agrees with the squared-distance test', () => {
        check(fc.tuple(wellScaledVector(2), fc.double({ min: 0.1, max: 5, noNaN: true }),
            wellScaledVector(2, -12, 12)),
            ([c, r, p]: [Vector, number, Vector]) => {
                const circle = Hypersphere.fromCenterRadius(c, r);
                const d = sub(p, c);
                const sqrLen = d.get(0) * d.get(0) + d.get(1) * d.get(1);
                // Away from the boundary the two forms cannot disagree.
                if (Math.abs(Math.sqrt(sqrLen) - r) > 1e-9) {
                    expect(inContainerCircle2(p, circle)).toBe(sqrLen < r * r);
                }
            });
    });

    // Rigid motions: rotating and translating the point set rotates and
    // translates the circle and leaves the radius unchanged.
    it('is equivariant under rigid motions', () => {
        check(fc.tuple(fc.array(wellScaledVector(2), { minLength: 1, maxLength: 10 }),
            rotationFrame(2), wellScaledVector(2)),
            ([points, frame, t]: [Vector[], Vector[], Vector]) => {
                const xform = (p: Vector): Vector =>
                    add(add(mul(p.get(0), frame[0]), mul(p.get(1), frame[1])), t);
                const c0 = getContainerCircle2(points);
                const c1 = getContainerCircle2(points.map(xform));
                const expected = xform(c0.center);
                // The rotation mixes coordinates, so a few ulps of the
                // coordinate magnitude (<= 10) are lost.
                expectClose(c1.center.get(0), expected.get(0), 1e-11, 1e-11);
                expectClose(c1.center.get(1), expected.get(1), 1e-11, 1e-11);
                expectClose(c1.radius, c0.radius, 1e-11, 1e-11);
            });
    });

    // The merged circle contains both inputs: for each input circle,
    // |Ci - Cm| + ri <= rm. This is exactly the design claim of the upstream
    // algorithm, and it holds to rounding error here (unlike the ellipse
    // merge of ContEllipse2, see upstream issue #292).
    it('merge contains both input circles', () => {
        const circleArb = fc.tuple(wellScaledVector(2),
            fc.double({ min: 0, max: 5, noNaN: true }))
            .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));
        check(fc.tuple(circleArb, circleArb),
            ([c0, c1]: [Circle2, Circle2]) => {
                const merge = mergeContainersCircle2(c0, c1);
                for (const input of [c0, c1]) {
                    const d = length(sub(input.center, merge.center));
                    expect(d + input.radius).toBeLessThanOrEqual(
                        merge.radius + 1e-12 + 1e-12 * merge.radius);
                }
                // The merge is also tight: it touches at least one input.
                const touch = [c0, c1].some(input =>
                    Math.abs(length(sub(input.center, merge.center))
                        + input.radius - merge.radius) <= 1e-9);
                expect(touch).toBe(true);
            });
    });

    // Sampling the boundary of each input confirms the analytic containment
    // test above with an independent construction.
    it('merge contains sampled boundary points of both inputs', () => {
        const rand = seededRandom(0x5eed1);
        for (let trial = 0; trial < 200; ++trial) {
            const mk = (): Circle2 => Hypersphere.fromCenterRadius(
                Vector.fromArray([10 * rand() - 5, 10 * rand() - 5]),
                4 * rand());
            const c0 = mk(), c1 = mk();
            const merge = mergeContainersCircle2(c0, c1);
            for (const input of [c0, c1]) {
                for (let k = 0; k < 16; ++k) {
                    const a = (2 * Math.PI * k) / 16;
                    const p = add(input.center, Vector.fromArray(
                        [input.radius * Math.cos(a), input.radius * Math.sin(a)]));
                    expect(length(sub(p, merge.center)))
                        .toBeLessThanOrEqual(merge.radius + 1e-9);
                }
            }
        }
    });

    // Degenerate inputs: coincident centers, and merging a circle with itself.
    it('handles coincident centers and self-merge', () => {
        check(fc.tuple(wellScaledVector(2), fc.double({ min: 0, max: 5, noNaN: true }),
            fc.double({ min: 0, max: 5, noNaN: true })),
            ([c, r0, r1]: [Vector, number, number]) => {
                const a = Hypersphere.fromCenterRadius(c, r0);
                const b = Hypersphere.fromCenterRadius(c, r1);
                // rDiffSqr >= lenSqr = 0, so the larger circle is returned.
                const merge = mergeContainersCircle2(a, b);
                expect(merge.radius).toBe(Math.max(r0, r1));
                expect(merge.center.get(0)).toBe(c.get(0));
                const self = mergeContainersCircle2(a, a);
                expect(self.radius).toBe(r0);
                // The result never aliases an input.
                expect(self.center).not.toBe(a.center);
            });
    });
});
