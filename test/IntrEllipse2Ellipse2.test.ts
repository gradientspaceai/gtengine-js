import { describe, it, expect } from 'vitest';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import { Matrix } from '../src/Matrix.js';
import { Vector, dot, sub } from '../src/Vector.js';
import {
    IntrEllipse2Ellipse2TI,
    IntrEllipse2Ellipse2FI,
    IntrEllipse2Ellipse2Classification as Cls,
    intrEllipse2Ellipse2InfinitePoints
} from '../src/IntrEllipse2Ellipse2.js';
import { add, mul } from '../src/Vector.js';
import {
    check, expectClose, fc, positive, rotationFrame, seededRandom,
    wellScaledVector
} from './helpers/arbitraries.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// An ellipse with the given center, orientation angle and extents.
function ellipse(cx: number, cy: number, angle: number, a: number,
    b: number): Hyperellipsoid {
    const c = Math.cos(angle), s = Math.sin(angle);
    return Hyperellipsoid.fromCenterAxisExtent(v2(cx, cy),
        [v2(c, s), v2(-s, c)], v2(a, b));
}

function circle(cx: number, cy: number, r: number): Hyperellipsoid {
    return ellipse(cx, cy, 0, r, r);
}

// The quadratic form (X-K)^T*M*(X-K) of the ellipse, evaluated at X. The
// value is 1 on the ellipse, less than 1 inside and greater than 1 outside.
function quadratic(e: Hyperellipsoid, x: Vector): number {
    const d = sub(x, e.center);
    let sum = 0;
    for (let i = 0; i < 2; ++i) {
        const t = dot(d, e.axis[i]) / e.extent.values[i];
        sum += t * t;
    }
    return sum;
}

// A point on the ellipse at the given parameter angle.
function pointAt(e: Hyperellipsoid, t: number): Vector {
    const u = e.axis[0], v = e.axis[1];
    const a = e.extent.values[0] * Math.cos(t);
    const b = e.extent.values[1] * Math.sin(t);
    return v2(e.center.values[0] + a * u.values[0] + b * v.values[0],
        e.center.values[1] + a * u.values[1] + b * v.values[1]);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The minimum and maximum of the quadratic form of 'outer' over a dense
// sampling of the boundary of 'inner'.
function sampleExtremes(outer: Hyperellipsoid, inner: Hyperellipsoid,
    numSamples: number = 4096): { min: number, max: number } {
    let min = Number.MAX_VALUE, max = -Number.MAX_VALUE;
    for (let i = 0; i < numSamples; ++i) {
        const q = quadratic(outer, pointAt(inner, 2 * Math.PI * i / numSamples));
        min = Math.min(min, q);
        max = Math.max(max, q);
    }
    return { min, max };
}

describe('IntrEllipse2Ellipse2TI', () => {
    const query = new IntrEllipse2Ellipse2TI();

    it('classifies identical ellipses as equal', () => {
        const e = ellipse(1, -2, 0.3, 3, 2);
        expect(query.test(e, e)).toBe(Cls.ELLIPSES_EQUAL);
    });

    it('classifies separated circles', () => {
        expect(query.test(circle(0, 0, 1), circle(5, 0, 1)))
            .toBe(Cls.ELLIPSES_SEPARATED);
    });

    it('classifies overlapping circles', () => {
        expect(query.test(circle(0, 0, 1), circle(1, 0, 1)))
            .toBe(Cls.ELLIPSES_OVERLAP);
    });

    it('classifies concentric containment both ways', () => {
        const big = circle(0, 0, 3);
        const small = circle(0, 0, 1);
        expect(query.test(big, small))
            .toBe(Cls.ELLIPSE0_STRICTLY_CONTAINS_ELLIPSE1);
        expect(query.test(small, big))
            .toBe(Cls.ELLIPSE1_STRICTLY_CONTAINS_ELLIPSE0);
    });

    it('classifies off-center containment', () => {
        const big = ellipse(0, 0, 0, 4, 3);
        const small = circle(1, 0.5, 0.5);
        expect(query.test(big, small))
            .toBe(Cls.ELLIPSE0_STRICTLY_CONTAINS_ELLIPSE1);
        expect(query.test(small, big))
            .toBe(Cls.ELLIPSE1_STRICTLY_CONTAINS_ELLIPSE0);
    });

    it('classifies overlapping ellipses of different orientation', () => {
        const e0 = ellipse(0, 0, 0, 3, 1);
        const e1 = ellipse(0, 0, Math.PI / 2, 3, 1);
        expect(query.test(e0, e1)).toBe(Cls.ELLIPSES_OVERLAP);
    });

    it('classifies external tangency of circles', () => {
        // The circles of radius 1 centered at (0,0) and (2,0) touch at (1,0).
        expect(query.test(circle(0, 0, 1), circle(2, 0, 1)))
            .toBe(Cls.ELLIPSE0_OUTSIDE_ELLIPSE1_BUT_TANGENT);
    });

    it('classifies internal tangency of circles', () => {
        // The circle of radius 1 at (1,0) touches the circle of radius 2 at
        // the origin, from the inside, at the point (2,0).
        expect(query.test(circle(0, 0, 2), circle(1, 0, 1)))
            .toBe(Cls.ELLIPSE0_CONTAINS_ELLIPSE1_BUT_TANGENT);
        expect(query.test(circle(1, 0, 1), circle(0, 0, 2)))
            .toBe(Cls.ELLIPSE1_CONTAINS_ELLIPSE0_BUT_TANGENT);
    });

    it('is unaffected by a rigid motion of both ellipses', () => {
        const rand = makeRandom(777);
        const angle = 0.7, c = Math.cos(angle), s = Math.sin(angle);
        const move = (e: Hyperellipsoid): Hyperellipsoid => {
            const k = e.center.values;
            const nc = v2(c * k[0] - s * k[1] + 3, s * k[0] + c * k[1] - 1);
            const rot = (u: Vector): Vector => v2(
                c * u.values[0] - s * u.values[1],
                s * u.values[0] + c * u.values[1]);
            return Hyperellipsoid.fromCenterAxisExtent(nc,
                [rot(e.axis[0]), rot(e.axis[1])], e.extent.clone());
        };
        for (let trial = 0; trial < 30; ++trial) {
            const e0 = ellipse(2 * rand() - 1, 2 * rand() - 1,
                Math.PI * rand(), 0.5 + rand(), 0.5 + rand());
            const e1 = ellipse(2 * rand() - 1, 2 * rand() - 1,
                Math.PI * rand(), 0.5 + rand(), 0.5 + rand());
            expect(query.test(move(e0), move(e1))).toBe(query.test(e0, e1));
        }
    });

    it('agrees with a sampled classification (randomized)', () => {
        const rand = makeRandom(20260901);
        const counts = new Map<number, number>();
        for (let trial = 0; trial < 200; ++trial) {
            const e0 = ellipse(2 * rand() - 1, 2 * rand() - 1,
                Math.PI * rand(), 0.5 + 1.5 * rand(), 0.5 + 1.5 * rand());
            const e1 = ellipse(3 * rand() - 1.5, 3 * rand() - 1.5,
                Math.PI * rand(), 0.3 + 2.2 * rand(), 0.3 + 2.2 * rand());

            // The quadratic form of e0 sampled over the boundary of e1.
            const { min, max } = sampleExtremes(e0, e1);
            // Skip configurations that are numerically close to tangency;
            // the sampled extremes are only approximate.
            if (Math.abs(min - 1) < 1e-2 || Math.abs(max - 1) < 1e-2) {
                continue;
            }

            let expected: Cls;
            if (max < 1) {
                expected = Cls.ELLIPSE0_STRICTLY_CONTAINS_ELLIPSE1;
            }
            else if (min < 1) {
                expected = Cls.ELLIPSES_OVERLAP;
            }
            else if (quadratic(e1, e0.center) < 1) {
                expected = Cls.ELLIPSE1_STRICTLY_CONTAINS_ELLIPSE0;
            }
            else {
                expected = Cls.ELLIPSES_SEPARATED;
            }

            const actual = query.test(e0, e1);
            expect(actual).toBe(expected);
            counts.set(actual, (counts.get(actual) ?? 0) + 1);
        }
        // Every non-degenerate classification is exercised.
        expect(counts.get(Cls.ELLIPSES_SEPARATED) ?? 0).toBeGreaterThan(5);
        expect(counts.get(Cls.ELLIPSES_OVERLAP) ?? 0).toBeGreaterThan(5);
        expect(counts.get(Cls.ELLIPSE0_STRICTLY_CONTAINS_ELLIPSE1) ?? 0)
            .toBeGreaterThan(0);
        expect(counts.get(Cls.ELLIPSE1_STRICTLY_CONTAINS_ELLIPSE0) ?? 0)
            .toBeGreaterThan(0);
    });

    it('throws for a non-2D hyperellipsoid', () => {
        const e3 = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([0, 0, 0]),
            [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
                Vector.fromArray([0, 0, 1])],
            Vector.fromArray([1, 1, 1]));
        expect(() => query.test(e3, e3)).toThrow();
    });
});

describe('IntrEllipse2Ellipse2FI', () => {
    const query = new IntrEllipse2Ellipse2FI();

    it('computes the standard form of an axis-aligned ellipse', () => {
        const { C, M } = query.getStandardForm(ellipse(1, 2, 0, 3, 2));
        expect(C.values).toEqual([1, 2]);
        expect(M.get(0, 0)).toBeCloseTo(1 / 9, 14);
        expect(M.get(1, 1)).toBeCloseTo(1 / 4, 14);
        expect(M.get(0, 1)).toBeCloseTo(0, 14);
        expect(M.get(1, 0)).toBeCloseTo(0, 14);
    });

    it('computes the standard form of a rotated ellipse', () => {
        const e = ellipse(-1, 0.5, 0.4, 3, 1);
        const { C, M } = query.getStandardForm(e);
        // The quadratic form built from M matches the axis/extent form.
        for (let i = 0; i < 8; ++i) {
            const x = pointAt(e, i);
            const d = sub(x, C);
            const q = d.values[0] * (M.get(0, 0) * d.values[0]
                + M.get(0, 1) * d.values[1])
                + d.values[1] * (M.get(1, 0) * d.values[0]
                    + M.get(1, 1) * d.values[1]);
            expect(q).toBeCloseTo(1, 12);
            expect(q).toBeCloseTo(quadratic(e, x), 12);
        }
    });

    it('computes the aligned box of an axis-aligned ellipse', () => {
        const box = query.computeAlignedBox(ellipse(1, 2, 0, 3, 2));
        expect(box.min.values[0]).toBeCloseTo(-2, 12);
        expect(box.max.values[0]).toBeCloseTo(4, 12);
        expect(box.min.values[1]).toBeCloseTo(0, 12);
        expect(box.max.values[1]).toBeCloseTo(4, 12);
    });

    it('computes the aligned box of a rotated ellipse', () => {
        const angle = 0.6, a = 3, b = 1;
        const e = ellipse(0, 0, angle, a, b);
        const box = query.computeAlignedBox(e);
        const c = Math.cos(angle), s = Math.sin(angle);
        const hx = Math.sqrt(a * a * c * c + b * b * s * s);
        const hy = Math.sqrt(a * a * s * s + b * b * c * c);
        expect(box.max.values[0]).toBeCloseTo(hx, 10);
        expect(box.max.values[1]).toBeCloseTo(hy, 10);
        expect(box.min.values[0]).toBeCloseTo(-hx, 10);
        expect(box.min.values[1]).toBeCloseTo(-hy, 10);
        // The extreme box points are on the ellipse.
        expect(quadratic(e, v2(hx, 0))).toBeGreaterThanOrEqual(1 - 1e-10);
    });

    it('reports infinitely many points for identical ellipses', () => {
        const e = ellipse(1, -2, 0.3, 3, 2);
        const result = query.find(e, e);
        expect(result.numPoints).toBe(intrEllipse2Ellipse2InfinitePoints);
        expect(result.intersect).toBe(false);
        for (const p of result.points) {
            expect(p.values).toEqual([0, 0]);
        }
    });

    it('reports no intersection for separated circles', () => {
        const result = query.find(circle(0, 0, 1), circle(5, 0, 1));
        expect(result.intersect).toBe(false);
        expect(result.numPoints).toBe(0);
    });

    it('reports no intersection for a circle inside a circle', () => {
        const result = query.find(circle(0, 0, 3), circle(0.5, 0, 1));
        expect(result.intersect).toBe(false);
        expect(result.numPoints).toBe(0);
    });

    it('finds the two intersections of overlapping circles', () => {
        // The unit circles at (0,0) and (1,0) meet at (1/2, +-sqrt(3)/2).
        const result = query.find(circle(0, 0, 1), circle(1, 0, 1));
        expect(result.intersect).toBe(true);
        expect(result.numPoints).toBe(2);
        const ys = [result.points[0].values[1], result.points[1].values[1]]
            .sort((p, q) => p - q);
        expect(result.points[0].values[0]).toBeCloseTo(0.5, 10);
        expect(result.points[1].values[0]).toBeCloseTo(0.5, 10);
        expect(ys[0]).toBeCloseTo(-Math.sqrt(3) / 2, 10);
        expect(ys[1]).toBeCloseTo(Math.sqrt(3) / 2, 10);
        expect(result.isTransverse[0]).toBe(true);
        expect(result.isTransverse[1]).toBe(true);
    });

    it('finds the four intersections of two crossed ellipses', () => {
        // x^2/4 + y^2 = 1 and x^2 + y^2/4 = 1 meet where x^2 = y^2 = 4/5.
        const e0 = ellipse(0, 0, 0, 2, 1);
        const e1 = ellipse(0, 0, Math.PI / 2, 2, 1);
        const result = query.find(e0, e1);
        expect(result.intersect).toBe(true);
        expect(result.numPoints).toBe(4);
        const t = Math.sqrt(4 / 5);
        for (let i = 0; i < 4; ++i) {
            expect(Math.abs(result.points[i].values[0])).toBeCloseTo(t, 7);
            expect(Math.abs(result.points[i].values[1])).toBeCloseTo(t, 7);
            expect(quadratic(e0, result.points[i])).toBeCloseTo(1, 7);
            expect(quadratic(e1, result.points[i])).toBeCloseTo(1, 7);
        }
    });

    it('finds the intersections of rotated overlapping ellipses', () => {
        const e0 = ellipse(0, 0, 0.3, 3, 1);
        const e1 = ellipse(1, 0.5, -0.7, 2, 1.5);
        const result = query.find(e0, e1);
        expect(result.intersect).toBe(true);
        expect(result.numPoints).toBeGreaterThan(0);
        for (let i = 0; i < result.numPoints; ++i) {
            expect(quadratic(e0, result.points[i])).toBeCloseTo(1, 8);
            expect(quadratic(e1, result.points[i])).toBeCloseTo(1, 8);
        }
    });

    it('gives the same answer with the early-exit test disabled', () => {
        const pairs: Array<[Hyperellipsoid, Hyperellipsoid]> = [
            [circle(0, 0, 1), circle(1, 0, 1)],
            [circle(0, 0, 1), circle(5, 0, 1)],
            [ellipse(0, 0, 0.3, 3, 1), ellipse(1, 0.5, -0.7, 2, 1.5)]
        ];
        for (const [e0, e1] of pairs) {
            const withTest = query.find(e0, e1, true);
            const without = query.find(e0, e1, false);
            expect(without.numPoints).toBe(withTest.numPoints);
            for (let i = 0; i < withTest.numPoints; ++i) {
                expect(without.points[i].values)
                    .toEqual(withTest.points[i].values);
            }
        }
    });

    it('accepts the standard forms directly', () => {
        const e0 = circle(0, 0, 1);
        const e1 = circle(1, 0, 1);
        const sf0 = query.getStandardForm(e0);
        const sf1 = query.getStandardForm(e1);
        const direct = query.findStandardForm(sf0.C, sf0.M, sf1.C, sf1.M);
        const viaEllipses = query.find(e0, e1);
        expect(direct.numPoints).toBe(viaEllipses.numPoints);
        for (let i = 0; i < direct.numPoints; ++i) {
            expect(direct.points[i].values)
                .toEqual(viaEllipses.points[i].values);
        }
    });

    it('reports a tangential contact of circles', () => {
        // The unit circles at (0,0) and (2,0) touch at (1,0).
        const result = query.find(circle(0, 0, 1), circle(2, 0, 1));
        expect(result.intersect).toBe(true);
        expect(result.numPoints).toBe(1);
        expect(result.points[0].values[0]).toBeCloseTo(1, 8);
        expect(result.points[0].values[1]).toBeCloseTo(0, 8);
        expect(result.isTransverse[0]).toBe(false);
    });

    it('detects identical standard forms even for equal-but-rebuilt ellipses',
        () => {
            const e = ellipse(0.25, -0.5, 0, 2, 1);
            const copy = ellipse(0.25, -0.5, 0, 2, 1);
            expect(query.find(e, copy).numPoints)
                .toBe(intrEllipse2Ellipse2InfinitePoints);
        });

    it('places every reported point on both ellipses (randomized)', () => {
        const rand = makeRandom(31415926);
        let totalPoints = 0;
        let numWith2 = 0, numWith4 = 0, numWith0 = 0;
        for (let trial = 0; trial < 200; ++trial) {
            const e0 = ellipse(2 * rand() - 1, 2 * rand() - 1,
                Math.PI * rand(), 0.5 + 1.5 * rand(), 0.5 + 1.5 * rand());
            const e1 = ellipse(3 * rand() - 1.5, 3 * rand() - 1.5,
                Math.PI * rand(), 0.5 + 1.5 * rand(), 0.5 + 1.5 * rand());
            const result = query.find(e0, e1);
            expect(result.numPoints).toBeLessThanOrEqual(4);
            expect(result.intersect).toBe(result.numPoints > 0);
            for (let i = 0; i < result.numPoints; ++i) {
                expect(quadratic(e0, result.points[i])).toBeCloseTo(1, 6);
                expect(quadratic(e1, result.points[i])).toBeCloseTo(1, 6);
            }
            totalPoints += result.numPoints;
            if (result.numPoints === 0) {
                ++numWith0;
            }
            else if (result.numPoints === 2) {
                ++numWith2;
            }
            else if (result.numPoints === 4) {
                ++numWith4;
            }
        }
        expect(totalPoints).toBeGreaterThan(100);
        expect(numWith0).toBeGreaterThan(10);
        expect(numWith2).toBeGreaterThan(10);
        expect(numWith4).toBeGreaterThan(0);
    });

    it('finds as many points as a boundary sign-change count (randomized)',
        () => {
            const rand = makeRandom(2718281);
            let checked = 0;
            for (let trial = 0; trial < 200; ++trial) {
                const e0 = ellipse(2 * rand() - 1, 2 * rand() - 1,
                    Math.PI * rand(), 0.5 + 1.5 * rand(), 0.5 + 1.5 * rand());
                const e1 = ellipse(3 * rand() - 1.5, 3 * rand() - 1.5,
                    Math.PI * rand(), 0.5 + 1.5 * rand(), 0.5 + 1.5 * rand());

                // Count sign changes of quadratic(e1,.) - 1 along the
                // boundary of e0; that is the number of transverse
                // intersections.
                const n = 2048;
                let numSignChanges = 0;
                let minAbs = Number.MAX_VALUE;
                let prev = quadratic(e1, pointAt(e0, 0)) - 1;
                for (let i = 1; i <= n; ++i) {
                    const cur = quadratic(e1, pointAt(e0, 2 * Math.PI * i / n))
                        - 1;
                    minAbs = Math.min(minAbs, Math.abs(cur));
                    if ((prev < 0) !== (cur < 0)) {
                        ++numSignChanges;
                    }
                    prev = cur;
                }
                // Skip near-tangential configurations, where the sampled
                // count is unreliable.
                if (minAbs < 1e-3) {
                    continue;
                }
                ++checked;
                const result = query.find(e0, e1);
                expect(result.numPoints).toBe(numSignChanges);
            }
            expect(checked).toBeGreaterThan(60);
        });

    it('is equivariant under a rigid motion of both ellipses', () => {
        const angle = 0.9, c = Math.cos(angle), s = Math.sin(angle);
        const tx = 2, ty = -3;
        const xform = (p: Vector): Vector => v2(
            c * p.values[0] - s * p.values[1] + tx,
            s * p.values[0] + c * p.values[1] + ty);
        const move = (e: Hyperellipsoid): Hyperellipsoid => {
            const rot = (u: Vector): Vector => v2(
                c * u.values[0] - s * u.values[1],
                s * u.values[0] + c * u.values[1]);
            return Hyperellipsoid.fromCenterAxisExtent(xform(e.center),
                [rot(e.axis[0]), rot(e.axis[1])], e.extent.clone());
        };
        const e0 = ellipse(0, 0, 0.3, 3, 1);
        const e1 = ellipse(1, 0.5, -0.7, 2, 1.5);
        const base = query.find(e0, e1);
        const moved = query.find(move(e0), move(e1));
        expect(moved.numPoints).toBe(base.numPoints);
        // Match the point sets (the order need not be preserved).
        for (let i = 0; i < base.numPoints; ++i) {
            const target = xform(base.points[i]);
            let best = Number.MAX_VALUE;
            for (let j = 0; j < moved.numPoints; ++j) {
                const d = sub(moved.points[j], target);
                best = Math.min(best, Math.hypot(d.values[0], d.values[1]));
            }
            expect(best).toBeLessThan(1e-8);
        }
    });

    it('does not modify the input ellipses', () => {
        const e0 = ellipse(0, 0, 0.3, 3, 1);
        const e1 = ellipse(1, 0.5, -0.7, 2, 1.5);
        const before = [e0.center.values.slice(), e1.center.values.slice()];
        query.find(e0, e1);
        expect(e0.center.values).toEqual(before[0]);
        expect(e1.center.values).toEqual(before[1]);
    });

    it('throws for a non-2D hyperellipsoid', () => {
        const e3 = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([0, 0, 0]),
            [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
                Vector.fromArray([0, 0, 1])],
            Vector.fromArray([1, 1, 1]));
        expect(() => query.find(e3, e3)).toThrow();
    });

    it('handles ellipses sharing a center but differing in shape', () => {
        // Concentric, axis-aligned, one strictly inside the other: upstream
        // takes the branch where all of e1, e2, e3, e4 are zero and reports
        // no intersection.
        const result = query.find(ellipse(0, 0, 0, 3, 2), ellipse(0, 0, 0, 1,
            0.5));
        expect(result.intersect).toBe(false);
        expect(result.numPoints).toBe(0);
    });

    it('finds the intersections of concentric crossed ellipses', () => {
        const e0 = ellipse(0, 0, 0, 3, 1);
        const e1 = ellipse(0, 0, 0, 1, 3);
        const result = query.find(e0, e1);
        expect(result.numPoints).toBe(4);
        for (let i = 0; i < 4; ++i) {
            expect(quadratic(e0, result.points[i])).toBeCloseTo(1, 8);
            expect(quadratic(e1, result.points[i])).toBeCloseTo(1, 8);
        }
    });

    it('finds the intersections when the ellipses share an axis line', () => {
        // Centers differ along x only, both axis-aligned: this exercises the
        // e4 = 0 branches (the M matrices are diagonal).
        const e0 = ellipse(0, 0, 0, 2, 1);
        const e1 = ellipse(1, 0, 0, 2, 1);
        const result = query.find(e0, e1);
        expect(result.numPoints).toBe(2);
        for (let i = 0; i < 2; ++i) {
            expect(result.points[i].values[0]).toBeCloseTo(0.5, 10);
            expect(quadratic(e0, result.points[i])).toBeCloseTo(1, 8);
            expect(quadratic(e1, result.points[i])).toBeCloseTo(1, 8);
        }
    });

    it('reports both points of a symmetric double root (upstream bug fix)',
        () => {
            // The unit circle at the origin, and the ellipse with matrix
            // [[3,1/2],[1/2,1/2]] centered at (0,1/2). Two of the four
            // intersections are (1/2, +-sqrt(3)/2), which share a first
            // coordinate. In upstream's CaseE4NotZero this is the
            // 'divisor == 0' branch, where both points are written to the
            // same slot of result.points and numPoints is incremented once,
            // so the first of the pair is lost and only 3 points are
            // reported. The port stores both, so all 4 are reported. The
            // coefficients here are dyadic rationals, so the quartic root at
            // 1/2 is exact and the branch really is taken.
            const M0 = Matrix.fromArray(2, 2, [1, 0, 0, 1]);
            const M1 = Matrix.fromArray(2, 2, [3, 0.5, 0.5, 0.5]);
            const result = query.findStandardForm(v2(0, 0), M0, v2(0, 0.5), M1);
            expect(result.intersect).toBe(true);
            expect(result.numPoints).toBe(4);

            // The symmetric pair is present with both signs.
            const pair = [];
            for (let i = 0; i < 4; ++i) {
                const p = result.points[i];
                // Every point is on the unit circle...
                expect(Math.hypot(p.values[0], p.values[1]))
                    .toBeCloseTo(1, 10);
                // ...and on the second ellipse.
                const dx = p.values[0], dy = p.values[1] - 0.5;
                expect(3 * dx * dx + 2 * 0.5 * dx * dy + 0.5 * dy * dy)
                    .toBeCloseTo(1, 10);
                if (Math.abs(p.values[0] - 0.5) < 1e-12) {
                    pair.push(p.values[1]);
                }
            }
            pair.sort((a, b) => a - b);
            expect(pair).toHaveLength(2);
            expect(pair[0]).toBeCloseTo(-Math.sqrt(3) / 2, 12);
            expect(pair[1]).toBeCloseTo(Math.sqrt(3) / 2, 12);
        });

    it('accepts a standard form built by hand', () => {
        // The unit circle at the origin and the circle of radius 2 at (2,0),
        // given directly as (center, matrix) pairs. They meet where
        // x = 1/4 and y^2 = 15/16.
        const M0 = Matrix.fromArray(2, 2, [1, 0, 0, 1]);
        const M1 = Matrix.fromArray(2, 2, [0.25, 0, 0, 0.25]);
        const result = query.findStandardForm(v2(0, 0), M0, v2(2, 0), M1);
        expect(result.numPoints).toBe(2);
        for (let i = 0; i < 2; ++i) {
            const p = result.points[i];
            expect(Math.hypot(p.values[0], p.values[1])).toBeCloseTo(1, 10);
            expect(Math.hypot(p.values[0] - 2, p.values[1])).toBeCloseTo(2, 10);
            expect(p.values[0]).toBeCloseTo(0.25, 10);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V32): properties cross-checking the port against upstream
// IntrEllipse2Ellipse2.h.
// ---------------------------------------------------------------------------

describe('IntrEllipse2Ellipse2 verification', () => {
    const tiQuery = new IntrEllipse2Ellipse2TI();
    const fiQuery = new IntrEllipse2Ellipse2FI();

    // Upstream's own header says the number type "should support exact
    // rational arithmetic in order for the polynomial root construction to be
    // robust", and the port instantiates the floating-point path. Two
    // conditioning cliffs of that path are exercised by the deterministic
    // tests at the end of this block: the TI root bracketing asserts for
    // nearly concentric ellipses, and the FI quartic loses nearly double
    // roots. The randomized properties below therefore either predict only
    // the robust direction (no intersection) or draw from a well-conditioned
    // family of configurations.

    const arbEllipse = fc.tuple(wellScaledVector(2, -3, 3), rotationFrame(2),
        fc.array(positive(3, 0.4), { minLength: 2, maxLength: 2 }))
        .map(([c, frame, e]) => Hyperellipsoid.fromCenterAxisExtent(c,
            [frame[0], frame[1]], Vector.fromArray(e)));

    // Pairs are drawn from a seeded uniform generator rather than
    // fast-check doubles: fc's boundary bias (2.999999999999994, 0.4000000000000004)
    // produces nearly circular, nearly concentric pairs, which is exactly the
    // upstream cliff described above. Axis ratio >= 1.15 and centre distance
    // >= 0.05 keep the pair out of that regime.
    function seededEllipse(rnd: () => number): Hyperellipsoid {
        const c = Vector.fromArray([6 * rnd() - 3, 6 * rnd() - 3]);
        const angle = 2 * Math.PI * rnd();
        const ca = Math.cos(angle), sa = Math.sin(angle);
        const e0 = 0.4 + 2.6 * rnd();
        const e1 = 0.4 + 2.6 * rnd();
        const big = Math.max(e0, e1), small = Math.min(e0, e1);
        const ratio = Math.max(1.15, big / small);
        return Hyperellipsoid.fromCenterAxisExtent(c,
            [Vector.fromArray([ca, sa]), Vector.fromArray([-sa, ca])],
            Vector.fromArray([small * ratio, small]));
    }
    const arbPair = fc.integer({ min: 1, max: 1000000 }).map(seed => {
        const rnd = seededRandom(seed);
        return [seededEllipse(rnd), seededEllipse(rnd)] as
            [Hyperellipsoid, Hyperellipsoid];
    }).filter(([e0, e1]) => {
        const dx = e0.center.get(0) - e1.center.get(0);
        const dy = e0.center.get(1) - e1.center.get(1);
        return dx * dx + dy * dy >= 0.05 * 0.05;
    });

    // A point on the boundary of the ellipse at the given parameter angle.
    function boundaryPoint(e: Hyperellipsoid, t: number): Vector {
        return add(e.center,
            add(mul(e.extent.get(0) * Math.cos(t), e.axis[0]),
                mul(e.extent.get(1) * Math.sin(t), e.axis[1])));
    }

    // Sample the boundary of 'e' and report the extreme values of the other
    // ellipse's quadratic form on it. The form is 1 on the boundary of the
    // other ellipse, less inside and more outside.
    function extremes(e: Hyperellipsoid, other: Hyperellipsoid, n: number):
        { min: number, max: number } {
        let lo = Number.MAX_VALUE, hi = -Number.MAX_VALUE;
        for (let k = 0; k < n; ++k) {
            const q = quadratic(other, boundaryPoint(e, (2 * Math.PI * k) / n));
            lo = Math.min(lo, q);
            hi = Math.max(hi, q);
        }
        return { min: lo, max: hi };
    }

    it('the FI points lie on both ellipses', () => {
        check(arbPair, ([e0, e1]) => {
            const res = fiQuery.find(e0, e1);
            if (res.numPoints === intrEllipse2Ellipse2InfinitePoints) {
                return;
            }
            expect(res.numPoints).toBeLessThanOrEqual(4);
            expect(res.intersect).toBe(res.numPoints > 0);
            for (let i = 0; i < res.numPoints; ++i) {
                const p = res.points[i];
                expect(Number.isNaN(p.get(0))).toBe(false);
                expect(Number.isNaN(p.get(1))).toBe(false);
                // The residual scales with the conditioning of the quartic
                // built from the standard forms, and for nearly congruent or
                // nearly circular ellipses (which fast-check shrinks towards)
                // it has no useful bound - the deterministic well-conditioned
                // family below is what checks the residual to 1e-6. Here the
                // point is only required to be a plausible intersection: on
                // both curves to within a few per cent of the quadratic form
                // and inside the bounding box of each ellipse.
                expectClose(quadratic(e0, p), 1, 0.05, 0.05);
                expectClose(quadratic(e1, p), 1, 0.05, 0.05);
                for (const e of [e0, e1]) {
                    const box = fiQuery.computeAlignedBox(e);
                    for (let d = 0; d < 2; ++d) {
                        // Same few-per-cent scale as the quadratic-form check
                        // above: a root off the curve by O(1e-3) can also sit
                        // just outside the box (observed 5e-4 for a
                        // near-tangent pair).
                        const tol = 1e-2 * (1 + Math.abs(p.get(d)));
                        expect(p.get(d))
                            .toBeGreaterThanOrEqual(box.min.get(d) - tol);
                        expect(p.get(d))
                            .toBeLessThanOrEqual(box.max.get(d) + tol);
                    }
                }
            }
        }, 100);
    });

    it('the FI query finds no points when the ellipses are robustly separated'
        + ' or nested', () => {
        check(arbPair, ([e0, e1]) => {
            const on1 = extremes(e1, e0, 512);
            const on0 = extremes(e0, e1, 512);
            const margin = 0.05;
            if ((on1.min > 1 + margin && on0.min > 1 + margin)
                || on1.max < 1 - margin || on0.max < 1 - margin) {
                const res = fiQuery.find(e0, e1);
                expect(res.numPoints).toBe(0);
                expect(res.intersect).toBe(false);
            }
        }, 100);
    });

    it('the early-exit bounding-box test only ever removes intersections',
        () => {
            check(arbPair, ([e0, e1]) => {
                const withExit = fiQuery.find(e0, e1, true);
                const without = fiQuery.find(e0, e1, false);
                const box0 = fiQuery.computeAlignedBox(e0);
                const box1 = fiQuery.computeAlignedBox(e1);
                let disjoint = false;
                for (let i = 0; i < 2; ++i) {
                    disjoint = disjoint || box0.max.get(i) < box1.min.get(i)
                        || box0.min.get(i) > box1.max.get(i);
                }
                if (disjoint) {
                    // The boxes bound the ellipses, so there is no
                    // intersection to report.
                    expect(withExit.numPoints).toBe(0);
                    expect(withExit.intersect).toBe(false);
                }
                else {
                    expect(withExit.numPoints).toBe(without.numPoints);
                    expect(withExit.intersect).toBe(without.intersect);
                    for (let i = 0; i < withExit.numPoints; ++i) {
                        expect(withExit.points[i].equals(without.points[i]))
                            .toBe(true);
                    }
                }
            });
        });

    it('the quartic solver can return spurious roots for nearly congruent'
        + ' offset ellipses (upstream conditioning limitation)', () => {
        // Two nearly identical ellipses whose bounding boxes are disjoint in
        // x, so there is no intersection. With the early-exit test disabled
        // the quartic solver returns two roots anyway; the early exit hides
        // them, which is why the query enables it by default.
        const e0 = ellipse(0, 0, 0, 0.40000000000000024, 2.999999999999992);
        const e1 = ellipse(0.9376363112298954, 0, 0, 0.40000000000000097,
            2.999999999999998);
        expect(fiQuery.find(e0, e1, true).numPoints).toBe(0);
        expect(fiQuery.find(e0, e1, false).numPoints).toBe(2);
        // The spurious points are not on the ellipses.
        const spurious = fiQuery.find(e0, e1, false);
        let offCurve = 0;
        for (let i = 0; i < spurious.numPoints; ++i) {
            if (Math.abs(quadratic(e0, spurious.points[i]) - 1) > 1e-3
                || Math.abs(quadratic(e1, spurious.points[i]) - 1) > 1e-3) {
                ++offCurve;
            }
        }
        expect(offCurve).toBeGreaterThan(0);
    });

    it('does not modify its inputs', () => {
        check(arbPair, ([e0, e1]) => {
            const c0 = e0.center.clone(), c1 = e1.center.clone();
            fiQuery.find(e0, e1);
            try {
                tiQuery.test(e0, e1);
            }
            catch {
                // Documented upstream assert; the inputs must still be intact.
            }
            expect(e0.center.equals(c0)).toBe(true);
            expect(e1.center.equals(c1)).toBe(true);
        });
    });

    it('getStandardForm uses the upstream Matrix operator/ semantics', () => {
        // Upstream evaluates M = (UUTrn / aSqr + VVTrn / bSqr) / USqrLen with
        // the Matrix operator/, which multiplies by the reciprocal of the
        // divisor. Dividing each element instead rounds differently, which
        // matters because every downstream branch of this file is selected by
        // an exact comparison against zero.
        check(arbEllipse, e => {
            const { M } = fiQuery.getStandardForm(e);
            const u = e.axis[0], v = e.axis[1];
            const aSqr = e.extent.get(0) * e.extent.get(0);
            const bSqr = e.extent.get(1) * e.extent.get(1);
            const uSqrLen = u.get(0) * u.get(0) + u.get(1) * u.get(1);
            const invA = 1 / aSqr, invB = 1 / bSqr, invU = 1 / uSqrLen;
            for (let r = 0; r < 2; ++r) {
                for (let c = 0; c < 2; ++c) {
                    const uu = u.get(r) * u.get(c);
                    const vv = v.get(r) * v.get(c);
                    expect(M.get(r, c)).toBe((uu * invA + vv * invB) * invU);
                }
            }
        });
    });

    it('the standard form is symmetric, positive definite and reproduces the'
        + ' ellipse', () => {
        check(arbEllipse, e => {
            const { C, M } = fiQuery.getStandardForm(e);
            expect(M.get(0, 1)).toBe(M.get(1, 0));
            expect(M.get(0, 0)).toBeGreaterThan(0);
            const det = M.get(0, 0) * M.get(1, 1) - M.get(0, 1) * M.get(1, 0);
            expect(det).toBeGreaterThan(0);
            expect(C.equals(e.center)).toBe(true);
            for (let k = 0; k < 8; ++k) {
                const d = sub(boundaryPoint(e, (2 * Math.PI * k) / 8), C);
                const q = d.get(0) * (M.get(0, 0) * d.get(0)
                    + M.get(0, 1) * d.get(1))
                    + d.get(1) * (M.get(1, 0) * d.get(0)
                        + M.get(1, 1) * d.get(1));
                expectClose(q, 1, 1e-9, 1e-9);
            }
        });
    });

    it('the aligned box bounds the ellipse', () => {
        check(arbEllipse, e => {
            const b = fiQuery.computeAlignedBox(e);
            for (let k = 0; k < 64; ++k) {
                const p = boundaryPoint(e, (2 * Math.PI * k) / 64);
                for (let i = 0; i < 2; ++i) {
                    expect(p.get(i))
                        .toBeGreaterThanOrEqual(b.min.get(i) - 1e-9);
                    expect(p.get(i)).toBeLessThanOrEqual(b.max.get(i) + 1e-9);
                }
            }
        });
    });

    // A deterministic family of well-conditioned configurations: the centers
    // are clearly separated, the extents are clearly distinct (no near-circle
    // and no near-congruent pair) and the orientations are unrelated. This
    // avoids the conditioning cliffs documented at the end of the block, and
    // being deterministic it cannot be shrunk onto them.
    function wellConditionedPairs(seed: number, count: number):
        Array<[Hyperellipsoid, Hyperellipsoid]> {
        const rnd = seededRandom(seed);
        const pairs: Array<[Hyperellipsoid, Hyperellipsoid]> = [];
        const make = (cx: number, cy: number, angle: number, a: number,
            b: number): Hyperellipsoid => {
            const c = Math.cos(angle), s = Math.sin(angle);
            return Hyperellipsoid.fromCenterAxisExtent(v2(cx, cy),
                [v2(c, s), v2(-s, c)], v2(a, b));
        };
        while (pairs.length < count) {
            const a0 = 0.6 + 2 * rnd(), b0 = 0.6 + 2 * rnd();
            const a1 = 0.6 + 2 * rnd(), b1 = 0.6 + 2 * rnd();
            if (Math.abs(a0 - b0) < 0.4 || Math.abs(a1 - b1) < 0.4) {
                continue;   // near-circular: e4 is nearly zero
            }
            if (Math.abs(a0 - b1) < 0.3 && Math.abs(b0 - a1) < 0.3) {
                continue;   // nearly congruent under a quarter turn
            }
            if (Math.abs(a0 - a1) < 0.3 && Math.abs(b0 - b1) < 0.3) {
                continue;   // nearly congruent
            }
            const dx = 4 * rnd() - 2, dy = 4 * rnd() - 2;
            if (Math.sqrt(dx * dx + dy * dy) < 0.5) {
                continue;   // nearly concentric
            }
            pairs.push([make(0, 0, Math.PI * rnd(), a0, b0),
                make(dx, dy, Math.PI * rnd(), a1, b1)]);
        }
        return pairs;
    }

    it('the TI classification matches a boundary sampling on well-conditioned'
        + ' configurations', () => {
        let checked = 0;
        for (const [e0, e1] of wellConditionedPairs(0xC0FFEE, 400)) {
            const on1 = extremes(e1, e0, 4096);
            const on0 = extremes(e0, e1, 4096);
            const margin = 1e-3;
            let expected: Cls | undefined;
            if (on1.min > 1 + margin && on0.min > 1 + margin) {
                expected = Cls.ELLIPSES_SEPARATED;
            }
            else if (on1.max < 1 - margin) {
                expected = Cls.ELLIPSE0_STRICTLY_CONTAINS_ELLIPSE1;
            }
            else if (on0.max < 1 - margin) {
                expected = Cls.ELLIPSE1_STRICTLY_CONTAINS_ELLIPSE0;
            }
            else if (on1.min < 1 - margin && on1.max > 1 + margin) {
                expected = Cls.ELLIPSES_OVERLAP;
            }
            if (expected !== undefined) {
                expect(tiQuery.test(e0, e1)).toBe(expected);
                ++checked;
            }
        }
        expect(checked).toBeGreaterThan(300);
    }, 30000);

    it('the number of reported points matches the boundary sign changes on'
        + ' well-conditioned configurations', () => {
        let crossingCases = 0;
        for (const [e0, e1] of wellConditionedPairs(0x5EED, 400)) {
            const n = 4096;
            const values: number[] = [];
            for (let k = 0; k < n; ++k) {
                values.push(quadratic(e0, boundaryPoint(e1,
                    (2 * Math.PI * k) / n)) - 1);
            }
            // Require every crossing to be well separated and transversal, so
            // that the count does not depend on the sampling resolution.
            const crossIndices: number[] = [];
            for (let k = 0; k < n; ++k) {
                if (values[k] * values[(k + 1) % n] < 0) {
                    crossIndices.push(k);
                }
            }
            let robust = crossIndices.length > 0;
            for (let i = 0; i < crossIndices.length && robust; ++i) {
                const j = crossIndices[(i + 1) % crossIndices.length];
                let sep = j - crossIndices[i];
                if (sep <= 0) {
                    sep += n;
                }
                robust = sep > n / 40;
                const slope = Math.abs(values[(crossIndices[i] + 1) % n]
                    - values[crossIndices[i]]) * n;
                robust = robust && slope > 0.5;
            }
            if (!robust) {
                continue;
            }
            const res = fiQuery.find(e0, e1);
            expect(res.numPoints).toBe(crossIndices.length);
            for (let i = 0; i < res.numPoints; ++i) {
                expectClose(quadratic(e0, res.points[i]), 1, 1e-6, 1e-6);
                expectClose(quadratic(e1, res.points[i]), 1, 1e-6, 1e-6);
            }

            ++crossingCases;
        }
        expect(crossingCases).toBeGreaterThan(50);
    }, 30000);

    it('is equivariant under rigid motions on well-conditioned'
        + ' configurations', () => {
        const rnd = seededRandom(0xBEEF);
        for (const [e0, e1] of wellConditionedPairs(0xD00D, 200)) {
            const angle = 2 * Math.PI * rnd();
            const c = Math.cos(angle), s = Math.sin(angle);
            const shift = v2(6 * rnd() - 3, 6 * rnd() - 3);
            const rot = (v: Vector): Vector => v2(
                c * v.get(0) - s * v.get(1), s * v.get(0) + c * v.get(1));
            const moved = (e: Hyperellipsoid): Hyperellipsoid =>
                Hyperellipsoid.fromCenterAxisExtent(add(rot(e.center), shift),
                    [rot(e.axis[0]), rot(e.axis[1])], e.extent.clone());
            const a = fiQuery.find(e0, e1);
            const b = fiQuery.find(moved(e0), moved(e1));
            if (a.numPoints === intrEllipse2Ellipse2InfinitePoints) {
                continue;
            }
            expect(b.numPoints).toBe(a.numPoints);
            for (let i = 0; i < a.numPoints; ++i) {
                const q = add(rot(a.points[i]), shift);
                let best = Number.MAX_VALUE;
                for (let j = 0; j < b.numPoints; ++j) {
                    const d = sub(b.points[j], q);
                    best = Math.min(best, Math.sqrt(dot(d, d)));
                }
                expect(best).toBeLessThan(1e-6);
            }
        }
    }, 30000);

    it('the TI query reports the overlap of an ellipse whose centre differs'
        + ' from the other along one principal direction (upstream-bug'
        + ' regression)', () => {
        // The critical-point equation p_i * (1 - s*d_i) = -s*d_i*k_i has, when
        // k_j = 0, the extra solution s = 1/d_j with p_j free. That value is a
        // pole of the function f(s) whose roots upstream computes, and
        // upstream also drops the c_j = k_j^2 = 0 term from f, so two of the
        // four critical points of the squared distance are lost and the
        // classification is taken from the remaining pair.
        //
        // Here the unit circle and the ellipse with extents (1/2, 2) centred
        // at (0, 1/2) overlap; upstream reports that the ellipse contains the
        // circle.
        const unitCircle = circle(0, 0, 1);
        for (const dy of [0.5, 0.25, 1e-3, 1e-9]) {
            const e = ellipse(0, dy, 0, 0.5, 2);
            expect(tiQuery.test(unitCircle, e)).toBe(Cls.ELLIPSES_OVERLAP);
            expect(tiQuery.test(e, unitCircle)).toBe(Cls.ELLIPSES_OVERLAP);
        }
        // The point (1/2, 1/2) is on the ellipse and strictly inside the
        // circle, and (0, 5/2) is on the ellipse and outside it.
        const e = ellipse(0, 0.5, 0, 0.5, 2);
        expectClose(quadratic(e, v2(0.5, 0.5)), 1);
        expect(quadratic(unitCircle, v2(0.5, 0.5))).toBeLessThan(1);
        expectClose(quadratic(e, v2(0, 2.5)), 1);
        expect(quadratic(unitCircle, v2(0, 2.5))).toBeGreaterThan(1);
        // The exactly concentric case takes the K = 0 path and was already
        // classified correctly.
        expect(tiQuery.test(unitCircle, ellipse(0, 0, 0, 0.5, 2)))
            .toBe(Cls.ELLIPSES_OVERLAP);
    });

    it('the TI query asserts for some nearly concentric ellipses (upstream'
        + ' limitation, preserved)', () => {
        // The bound smax = (1 + sqrt(sum)) / d1 used to bracket the last root
        // collapses onto the pole 1/d1 when sum = d0*c0 + d1*c1 is tiny, that
        // is when the two centres nearly coincide. F(smax) is then a large
        // positive value instead of a nonpositive one and the upstream
        // LogAssert fires. The port preserves the assert (it throws an Error
        // carrying the upstream message).
        const e0 = ellipse(0, -0.001, 0, 0.4000000000000001,
            0.4000000000000001);
        const e1 = Hyperellipsoid.fromCenterAxisExtent(v2(0, 0),
            [v2(-0.9999999999999725, -2.3445664311150744e-7),
                v2(2.3445664311150744e-7, -0.9999999999999725)],
            v2(0.4000000000000001, 0.800680078467124));
        expect(() => tiQuery.test(e0, e1)).toThrow('Unexpected condition.');
    });

    it('the FI query loses the crossings of a near-circular ellipse against an'
        + ' axis-aligned one (upstream conditioning limitation)', () => {
        // e1 is a circle up to one ulp, so e4 = d1 * g4 is a tiny nonzero
        // number, the quartic is nearly a quadratic and the quartic solver
        // loses the two real roots. The two boundaries do cross: the point
        // (0.4, 0) is inside e0 and (-0.4, 0) is outside it.
        const e0 = ellipse(2.674085229787872, 0, 0, 2.3341186669893186,
            0.4069937584944733);
        const e1 = Hyperellipsoid.fromCenterAxisExtent(v2(0, 0),
            [v2(0.9999994998840165, 0.0010001158516906427),
                v2(-0.0010001158516906427, 0.9999994998840165)],
            v2(0.4000000000000001, 0.40000000000000013));
        expect(quadratic(e0, v2(0.4, 0))).toBeLessThan(1);
        expect(quadratic(e0, v2(-0.4, 0))).toBeGreaterThan(1);
        // Documented failure: the query reports no intersection.
        expect(fiQuery.find(e0, e1).numPoints).toBe(0);
    });
});
