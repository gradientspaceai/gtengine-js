import { describe, expect, it } from 'vitest';
import { NaturalSplineCurve } from '../src/NaturalSplineCurve.js';
import { Vector, sub, length as vectorLength } from '../src/Vector.js';

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

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). NaturalSplineCurve.h was read
// line by line against src/NaturalSplineCurve.ts; the properties below pin the
// mathematics that determines the spline uniquely (interpolation, C2 joins and
// the boundary conditions of each factory) plus the exactness results that
// follow from that uniqueness, so a transcription error in any of the three
// coefficient solvers has to show up.
import {
    check, fc, expectClose, expectVectorClose as expectVecClose,
    wellScaledVector
} from './helpers/arbitraries.js';

// Strictly increasing times with well-separated knots, and points that are
// neither subnormal nor huge: the tridiagonal solves of CreateFree and
// CreateClamped and the Gaussian elimination of CreateClosed both lose
// relative accuracy on degenerate spacing, which would make the identities
// below meaningless.
const splineData = (dim: number, minPoints = 3, maxPoints = 7):
    fc.Arbitrary<{ points: Vector[], times: number[] }> =>
    fc.integer({ min: minPoints, max: maxPoints }).chain(n =>
        fc.tuple(
            fc.array(wellScaledVector(dim, -5, 5),
                { minLength: n, maxLength: n }),
            fc.array(fc.double({ min: 0.25, max: 2, noNaN: true,
                noDefaultInfinity: true }),
            { minLength: n - 1, maxLength: n - 1 }))
            .map(([points, gaps]) => {
                const t = [0];
                for (const g of gaps) { t.push(t[t.length - 1] + g); }
                return { points, times: t };
            }));

// The one-sided limit of the jet at the right end of segment 'i'. On that
// segment the curve is a cubic, so the third-order Taylor polynomial centered
// at the segment midpoint reproduces it exactly; this avoids the finite
// difference error of sampling just to the left of the knot.
function leftLimit(curve: NaturalSplineCurve, knots: readonly number[],
    i: number): { p: Vector, d1: Vector, d2: Vector } {
    const s = 0.5 * (knots[i] + knots[i + 1]);
    const h = knots[i + 1] - s;
    const j = jetOf(curve, s, 3);
    const p = new Vector(j[0].size);
    const d1 = new Vector(j[0].size);
    const d2 = new Vector(j[0].size);
    for (let k = 0; k < p.size; ++k) {
        p.values[k] = j[0].values[k] + h * j[1].values[k] +
            0.5 * h * h * j[2].values[k] + h * h * h * j[3].values[k] / 6;
        d1.values[k] = j[1].values[k] + h * j[2].values[k] +
            0.5 * h * h * j[3].values[k];
        d2.values[k] = j[2].values[k] + h * j[3].values[k];
    }
    return { p, d1, d2 };
}

// Interpolation at every knot plus C0/C1/C2 agreement across every interior
// knot. Together with the boundary conditions these determine the spline
// uniquely, so every solver is checked against the definition rather than
// against itself.
function expectInterpolatingC2(curve: NaturalSplineCurve, points: Vector[],
    knots: number[]): void {
    for (let i = 0; i < points.length; ++i) {
        expectVecClose(jetOf(curve, knots[i], 0)[0], points[i], 1e-8, 1e-8);
    }
    for (let i = 1; i + 1 < knots.length; ++i) {
        const left = leftLimit(curve, knots, i - 1);
        const right = jetOf(curve, knots[i], 2);
        expectVecClose(left.p, right[0], 1e-8, 1e-8);
        expectVecClose(left.d1, right[1], 1e-8, 1e-8);
        expectVecClose(left.d2, right[2], 1e-8, 1e-8);
    }
}

describe('NaturalSplineCurve verification', () => {
    it('free splines interpolate, are C2 and have zero end curvature', () => {
        check(splineData(2), ({ points, times: knots }) => {
            const curve = NaturalSplineCurve.createFree(points, knots);
            expectInterpolatingC2(curve, points, knots);
            const zero = new Vector(2);
            expectVecClose(jetOf(curve, curve.getTMin(), 2)[2], zero, 1e-8, 1e-8);
            expectVecClose(jetOf(curve, curve.getTMax(), 2)[2], zero, 1e-8, 1e-8);
        });
    });

    it('clamped splines interpolate, are C2 and match the end tangents', () => {
        check(fc.tuple(splineData(3), wellScaledVector(3, -4, 4),
            wellScaledVector(3, -4, 4)),
        ([{ points, times: knots }, d0, d1]) => {
            const curve = NaturalSplineCurve.createClamped(points, knots, d0, d1);
            expectInterpolatingC2(curve, points, knots);
            expectVecClose(jetOf(curve, curve.getTMin(), 1)[1], d0, 1e-8, 1e-8);
            expectVecClose(jetOf(curve, curve.getTMax(), 1)[1], d1, 1e-8, 1e-8);
        });
    });

    it('closed splines interpolate, are C2 and are periodic at the wrap', () => {
        // A closed spline treats the wrap as a knot, so the periodicity of the
        // position, the tangent and the second derivative is what the
        // wrap-around row of the matrix buys. Two segments (three points, the
        // last repeating the first) is the case upstream gets wrong.
        check(splineData(2, 3, 7), ({ points, times: knots }) => {
            const closed = points.map(p => p.clone());
            closed[closed.length - 1] = closed[0].clone();
            const curve = NaturalSplineCurve.createClosed(closed, knots);
            expectInterpolatingC2(curve, closed, knots);
            const atMin = jetOf(curve, curve.getTMin(), 2);
            const atMax = jetOf(curve, curve.getTMax(), 2);
            expectVecClose(atMin[0], atMax[0], 1e-8, 1e-8);
            expectVecClose(atMin[1], atMax[1], 1e-8, 1e-8);
            expectVecClose(atMin[2], atMax[2], 1e-8, 1e-8);
        });
    });

    it('clamped interpolation reproduces cubic polynomials exactly', () => {
        // A clamped cubic spline is the unique C2 interpolant with the given
        // end tangents, and a cubic polynomial satisfies all of those
        // conditions, so the two must agree everywhere. This is an independent
        // closed-form cross-check of CreateClamped.
        const cubic = fc.array(wellScaledVector(2, -2, 2),
            { minLength: 4, maxLength: 4 });
        check(fc.tuple(splineData(2, 4, 6), cubic, fc.double({ min: 0, max: 1,
            noNaN: true, noDefaultInfinity: true })),
        ([{ times: knots }, c, u]) => {
            const at = (t: number): Vector => Vector.fromArray([0, 1].map(k =>
                c[0].values[k] + t * (c[1].values[k] + t * (c[2].values[k] +
                    t * c[3].values[k]))));
            const der = (t: number): Vector => Vector.fromArray([0, 1].map(k =>
                c[1].values[k] + t * (2 * c[2].values[k] +
                    3 * t * c[3].values[k])));
            const points = knots.map(at);
            const last = knots[knots.length - 1];
            const curve = NaturalSplineCurve.createClamped(points, knots,
                der(knots[0]), der(last));
            const t = knots[0] + u * (last - knots[0]);
            expectVecClose(jetOf(curve, t, 1)[0], at(t), 1e-8, 1e-8);
            expectVecClose(jetOf(curve, t, 1)[1], der(t), 1e-8, 1e-8);
        });
    });

    it('free interpolation reproduces affine functions exactly', () => {
        // An affine function has zero second derivative everywhere, so it also
        // satisfies the free (natural) boundary conditions.
        check(fc.tuple(splineData(2, 3, 7), wellScaledVector(2, -3, 3),
            wellScaledVector(2, -3, 3),
            fc.double({ min: 0, max: 1, noNaN: true, noDefaultInfinity: true })),
        ([{ times: knots }, a, b, u]) => {
            const at = (t: number): Vector => Vector.fromArray([
                a.values[0] + t * b.values[0], a.values[1] + t * b.values[1]]);
            const last = knots[knots.length - 1];
            const curve = NaturalSplineCurve.createFree(knots.map(at), knots);
            const t = knots[0] + u * (last - knots[0]);
            expectVecClose(jetOf(curve, t, 2)[0], at(t), 1e-9, 1e-9);
            expectVecClose(jetOf(curve, t, 2)[1], b, 1e-9, 1e-9);
            expectVecClose(jetOf(curve, t, 2)[2], new Vector(2), 1e-9, 1e-9);
        });
    });

    it('is equivariant under translation and scaling of the points', () => {
        // For fixed times the free spline depends linearly on the points, so
        // scaling the points scales the curve and translating them translates
        // it. A sign or index error in one coefficient of B, C or D breaks
        // this.
        check(fc.tuple(splineData(2), wellScaledVector(2, -6, 6),
            fc.double({ min: -3, max: 3, noNaN: true, noDefaultInfinity: true })
                .filter(s => Math.abs(s) > 0.1),
            fc.double({ min: 0, max: 1, noNaN: true, noDefaultInfinity: true })),
        ([{ points, times: knots }, shift, scale, u]) => {
            const base = NaturalSplineCurve.createFree(points, knots);
            const moved = NaturalSplineCurve.createFree(points.map(p =>
                Vector.fromArray([scale * p.values[0] + shift.values[0],
                    scale * p.values[1] + shift.values[1]])), knots);
            const last = knots[knots.length - 1];
            const t = knots[0] + u * (last - knots[0]);
            const p0 = jetOf(base, t, 1);
            const p1 = jetOf(moved, t, 1);
            expectClose(p1[0].values[0],
                scale * p0[0].values[0] + shift.values[0], 1e-8, 1e-8);
            expectClose(p1[0].values[1],
                scale * p0[0].values[1] + shift.values[1], 1e-8, 1e-8);
            expectClose(p1[1].values[0], scale * p0[1].values[0], 1e-8, 1e-8);
            expectClose(p1[1].values[1], scale * p0[1].values[1], 1e-8, 1e-8);
        });
    });

    it('evaluation outside the domain clamps to the endpoint values', () => {
        check(fc.tuple(splineData(2), fc.double({ min: 0.001, max: 50,
            noNaN: true, noDefaultInfinity: true })),
        ([{ points, times: knots }, d]) => {
            const curve = NaturalSplineCurve.createFree(points, knots);
            const last = knots.length - 1;
            expectVecClose(jetOf(curve, knots[0] - d, 0)[0], points[0],
                1e-9, 1e-9);
            expectVecClose(jetOf(curve, knots[last] + d, 0)[0], points[last],
                1e-8, 1e-8);
        });
    });

    it('leaves the jet slots above the requested order untouched', () => {
        // Upstream fills jet[i] only for i <= order (and zeroes i >= 4 only
        // when order >= 3), so a caller-supplied slot above the order must
        // keep its value.
        const curve = NaturalSplineCurve.createFree(points2D, times);
        const jet = curve.createJet();
        jet[2] = v2(7, 8);
        jet[3] = v2(9, 10);
        curve.evaluate(2.0, 1, jet);
        expect(jet[2].values).toEqual([7, 8]);
        expect(jet[3].values).toEqual([9, 10]);
    });
});
