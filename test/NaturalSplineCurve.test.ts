import { describe, expect, it } from 'vitest';
import { NaturalSplineCurve } from '../src/NaturalSplineCurve';
import { Vector, sub, length as vectorLength } from '../src/Vector';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function jetOf(curve: NaturalSplineCurve, t: number, order: number): Vector[] {
    const jet = curve.createJet();
    curve.evaluate(t, order, jet);
    return jet;
}

function expectVectorClose(actual: Vector, expected: Vector,
    tolerance: number): void {
    expect(vectorLength(sub(actual, expected))).toBeLessThan(tolerance);
}

const points2D = [
    v2(0, 0), v2(1, 2), v2(3, 1), v2(4, 4), v2(6, 0), v2(7, 3)
];
const times = [0, 1, 2.5, 3, 4.25, 5];

// A closed spline requires the last point to coincide with the first so that
// the wrap-around equation describes a periodic curve.
const closedPoints = [
    v2(0, 0), v2(2, 1), v2(3, 3), v2(1, 4), v2(-1, 2), v2(0, 0)
];

describe('NaturalSplineCurve', () => {
    it('interpolates the input points at the input times', () => {
        for (const curve of [
            NaturalSplineCurve.createFree(points2D, times),
            NaturalSplineCurve.createClamped(points2D, times,
                v2(1, 2), v2(1, 3)),
            NaturalSplineCurve.createClosed(closedPoints, times)
        ]) {
            const pts = curve.getPoints();
            for (let i = 0; i < pts.length; ++i) {
                const jet = jetOf(curve, times[i], 0);
                expectVectorClose(jet[0], pts[i], 1e-10);
            }
        }
    });

    it('reports the construction state and member access', () => {
        const curve = NaturalSplineCurve.createFree(points2D, times);
        expect(curve.isConstructed()).toBe(true);
        expect(curve.getNumPoints()).toBe(6);
        expect(curve.getNumSegments()).toBe(5);
        expect(curve.getDimension()).toBe(2);
        expect(curve.getTMin()).toBeCloseTo(0, 12);
        expect(curve.getTMax()).toBeCloseTo(5, 12);
        // getPoints returns copies of the inputs, not the input objects.
        expect(curve.getPoints()[2]).not.toBe(points2D[2]);
        expectVectorClose(curve.getPoints()[2], points2D[2], 1e-15);
    });

    it('is C2 across the interior knots', () => {
        const h = 1e-6;
        for (const curve of [
            NaturalSplineCurve.createFree(points2D, times),
            NaturalSplineCurve.createClamped(points2D, times,
                v2(1, 2), v2(1, 3)),
            NaturalSplineCurve.createClosed(closedPoints, times)
        ]) {
            for (let i = 1; i + 1 < times.length; ++i) {
                const before = jetOf(curve, times[i] - h, 2);
                const after = jetOf(curve, times[i] + h, 2);
                expectVectorClose(before[0], after[0], 1e-4);
                expectVectorClose(before[1], after[1], 1e-4);
                expectVectorClose(before[2], after[2], 1e-4);
            }
        }
    });

    it('has zero second derivative at the endpoints of a free spline', () => {
        const curve = NaturalSplineCurve.createFree(points2D, times);
        const atMin = jetOf(curve, curve.getTMin(), 2);
        const atMax = jetOf(curve, curve.getTMax(), 2);
        expectVectorClose(atMin[2], v2(0, 0), 1e-10);
        expectVectorClose(atMax[2], v2(0, 0), 1e-10);
    });

    it('matches the prescribed derivatives of a clamped spline', () => {
        const d0 = v2(1.5, -0.5);
        const d1 = v2(-2, 3.25);
        const curve = NaturalSplineCurve.createClamped(points2D, times, d0, d1);
        expectVectorClose(jetOf(curve, curve.getTMin(), 1)[1], d0, 1e-10);
        expectVectorClose(jetOf(curve, curve.getTMax(), 1)[1], d1, 1e-10);
    });

    it('is periodic in the first and second derivatives when closed', () => {
        const curve = NaturalSplineCurve.createClosed(closedPoints, times);
        const atMin = jetOf(curve, curve.getTMin(), 2);
        const atMax = jetOf(curve, curve.getTMax(), 2);
        // The closed system forces C[0] == C[numSegments] and the wrap-around
        // equation enforces continuity of the first derivative.
        expectVectorClose(atMin[0], atMax[0], 1e-10);
        expectVectorClose(atMin[1], atMax[1], 1e-10);
        expectVectorClose(atMin[2], atMax[2], 1e-10);
    });

    it('is periodic for a closed spline with only two segments', () => {
        // Upstream assembles the wrap-around row of the closed-spline matrix
        // with plain assignments, so for two segments the dt[numSm1] term is
        // overwritten by the dt[0] term and the system is wrong. The port
        // accumulates the terms instead.
        const p = [v2(0, 0), v2(2, 1), v2(0, 0)];
        const t = [0, 1, 2.5];
        const curve = NaturalSplineCurve.createClosed(p, t);
        const atMin = jetOf(curve, curve.getTMin(), 2);
        const atMax = jetOf(curve, curve.getTMax(), 2);
        expectVectorClose(atMin[0], atMax[0], 1e-10);
        expectVectorClose(atMin[1], atMax[1], 1e-10);
        expectVectorClose(atMin[2], atMax[2], 1e-10);
        expectVectorClose(jetOf(curve, 1, 0)[0], v2(2, 1), 1e-10);

        // The same identities must hold at the interior knot.
        const h = 1e-6;
        const before = jetOf(curve, 1 - h, 2);
        const after = jetOf(curve, 1 + h, 2);
        expectVectorClose(before[1], after[1], 1e-4);
        expectVectorClose(before[2], after[2], 1e-4);
    });

    it('reduces to a line segment for two points', () => {
        const p = [v2(1, 1), v2(4, 5)];
        const t = [0, 2];
        const curve = NaturalSplineCurve.createFree(p, t);
        for (let i = 0; i <= 8; ++i) {
            const s = i / 8;
            const jet = jetOf(curve, 2 * s, 3);
            expectVectorClose(jet[0], v2(1 + 3 * s, 1 + 4 * s), 1e-12);
            expectVectorClose(jet[1], v2(1.5, 2), 1e-12);
            expectVectorClose(jet[2], v2(0, 0), 1e-12);
            expectVectorClose(jet[3], v2(0, 0), 1e-12);
        }
    });

    it('clamps evaluation outside the domain to the endpoints', () => {
        const curve = NaturalSplineCurve.createFree(points2D, times);
        expectVectorClose(jetOf(curve, -3, 0)[0], points2D[0], 1e-12);
        expectVectorClose(jetOf(curve, 17, 0)[0], points2D[5], 1e-12);
    });

    it('evaluates the derivatives consistently with finite differences', () => {
        // Randomized cross-check: compare the analytic jet against centered
        // differences of the position and the first derivative.
        let seed = 987654321;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 8; ++trial) {
            const p: Vector[] = [];
            const t: number[] = [0];
            for (let i = 0; i < 5; ++i) {
                p.push(v2(4 * rand() - 2, 4 * rand() - 2));
                if (i > 0) {
                    t.push(t[t.length - 1] + 0.5 + rand());
                }
            }
            const curve = (trial % 2 === 0
                ? NaturalSplineCurve.createFree(p, t)
                : NaturalSplineCurve.createClamped(p, t, v2(0.5, -0.25),
                    v2(-1, 0.75)));

            const h = 1e-5;
            for (let k = 1; k <= 6; ++k) {
                const s = curve.getTMin() +
                    (k / 7) * (curve.getTMax() - curve.getTMin());
                const jet = jetOf(curve, s, 2);
                const fdD1 = Vector.fromArray([
                    (jetOf(curve, s + h, 0)[0].values[0] -
                        jetOf(curve, s - h, 0)[0].values[0]) / (2 * h),
                    (jetOf(curve, s + h, 0)[0].values[1] -
                        jetOf(curve, s - h, 0)[0].values[1]) / (2 * h)
                ]);
                const fdD2 = Vector.fromArray([
                    (jetOf(curve, s + h, 1)[1].values[0] -
                        jetOf(curve, s - h, 1)[1].values[0]) / (2 * h),
                    (jetOf(curve, s + h, 1)[1].values[1] -
                        jetOf(curve, s - h, 1)[1].values[1]) / (2 * h)
                ]);
                expectVectorClose(jet[1], fdD1, 1e-4);
                expectVectorClose(jet[2], fdD2, 1e-4);
            }
        }
    });

    it('supports dimensions other than two', () => {
        const p = [
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 2, 3]),
            Vector.fromArray([2, -1, 4]),
            Vector.fromArray([5, 0, 1])
        ];
        const t = [0, 1, 2, 3];
        const curve = NaturalSplineCurve.createFree(p, t);
        expect(curve.getDimension()).toBe(3);
        for (let i = 0; i < 4; ++i) {
            expectVectorClose(jetOf(curve, t[i], 0)[0], p[i], 1e-10);
        }
        expectVectorClose(jetOf(curve, 0, 2)[2], Vector.fromArray([0, 0, 0]),
            1e-10);
    });

    it('rejects invalid input', () => {
        expect(() => NaturalSplineCurve.createFree([v2(0, 0)], [0])).toThrow();
        expect(() => NaturalSplineCurve.createFree(
            [v2(0, 0), v2(1, 1)], [0])).toThrow();
        expect(() => NaturalSplineCurve.createClamped(
            [v2(0, 0), v2(1, 1)], [0, 1], v2(1, 1),
            Vector.fromArray([1, 1, 1]))).toThrow();
    });
});
