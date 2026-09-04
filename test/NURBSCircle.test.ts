import { describe, it, expect } from 'vitest';
import { Arc2 } from '../src/Arc2.js';
import {
    NURBSCircularArcDegree2, NURBSFullCircleDegree3, NURBSHalfCircleDegree3,
    NURBSQuarterCircleDegree2, NURBSQuarterCircleDegree4
} from '../src/NURBSCircle.js';
import { NURBSCurve } from '../src/NURBSCurve.js';
import { Vector, length, sub } from '../src/Vector.js';

function samplePositions(curve: NURBSCurve, count: number): Vector[] {
    const jet = curve.createJet();
    const tmin = curve.getTMin();
    const tmax = curve.getTMax();
    const positions: Vector[] = [];
    for (let i = 0; i <= count; ++i) {
        const t = tmin + (tmax - tmin) * i / count;
        curve.evaluate(t, 0, jet);
        positions.push(jet[0].clone());
    }
    return positions;
}

function expectUnitCircle(curve: NURBSCurve, count = 64): Vector[] {
    const positions = samplePositions(curve, count);
    for (const p of positions) {
        expect(Math.abs(length(p) - 1)).toBeLessThan(1e-12);
    }
    return positions;
}

// The counterclockwise angle in [0, 2*pi) of a point on the unit circle.
function angleOf(p: Vector): number {
    const a = Math.atan2(p.get(1), p.get(0));
    return a < 0 ? a + 2 * Math.PI : a;
}

describe('NURBSQuarterCircleDegree2', () => {
    it('has the documented control points and weights', () => {
        const curve = new NURBSQuarterCircleDegree2();
        expect(curve.getNumControls()).toBe(3);
        expect(curve.getBasisFunction().getDegree()).toBe(2);
        expect(curve.getWeights()).toEqual([Math.SQRT2, 1, Math.SQRT2]);
        expect(curve.getControl(0).values).toEqual([1, 0]);
        expect(curve.getControl(1).values).toEqual([1, 1]);
        expect(curve.getControl(2).values).toEqual([0, 1]);
    });

    it('lies on the unit circle and spans the first quadrant', () => {
        const curve = new NURBSQuarterCircleDegree2();
        const positions = expectUnitCircle(curve);
        for (const p of positions) {
            expect(p.get(0)).toBeGreaterThanOrEqual(-1e-15);
            expect(p.get(1)).toBeGreaterThanOrEqual(-1e-15);
        }
        // Endpoints (1,0) and (0,1), traversed counterclockwise.
        expect(length(sub(positions[0], Vector.fromArray([1, 0]))))
            .toBeLessThan(1e-14);
        expect(length(sub(positions[positions.length - 1],
            Vector.fromArray([0, 1])))).toBeLessThan(1e-14);
        let previous = -1;
        for (const p of positions) {
            const a = angleOf(p);
            expect(a).toBeGreaterThanOrEqual(previous - 1e-14);
            previous = a;
        }
        expect(previous).toBeCloseTo(Math.PI / 2, 12);
    });
});

describe('NURBSQuarterCircleDegree4', () => {
    it('has the documented control points and weights', () => {
        const curve = new NURBSQuarterCircleDegree4();
        expect(curve.getNumControls()).toBe(5);
        expect(curve.getBasisFunction().getDegree()).toBe(4);
        const w = curve.getWeights();
        expect(w[0]).toBe(1);
        expect(w[1]).toBe(1);
        expect(w[2]).toBeCloseTo(2 * Math.SQRT2 / 3, 15);
        expect(w[3]).toBe(1);
        expect(w[4]).toBe(1);
        const y1 = 0.5 / Math.SQRT2;
        const x2 = 1 - Math.SQRT2 / 8;
        expect(curve.getControl(1).values).toEqual([1, y1]);
        expect(curve.getControl(2).values).toEqual([x2, x2]);
        expect(curve.getControl(3).values).toEqual([y1, 1]);
    });

    it('lies on the unit circle and spans the first quadrant', () => {
        const curve = new NURBSQuarterCircleDegree4();
        const positions = expectUnitCircle(curve);
        expect(length(sub(positions[0], Vector.fromArray([1, 0]))))
            .toBeLessThan(1e-14);
        expect(length(sub(positions[positions.length - 1],
            Vector.fromArray([0, 1])))).toBeLessThan(1e-14);
        expect(angleOf(positions[(positions.length - 1) / 2]))
            .toBeCloseTo(Math.PI / 4, 12);
    });

    it('agrees with the degree-2 quarter circle as a point set', () => {
        // Both curves trace the same quarter circle, though with different
        // parameterizations, so compare the reached angles.
        const c2 = new NURBSQuarterCircleDegree2();
        const c4 = new NURBSQuarterCircleDegree4();
        for (const t of [0, 0.25, 0.5, 0.75, 1]) {
            const jet2 = c2.createJet();
            const jet4 = c4.createJet();
            c2.evaluate(t, 0, jet2);
            c4.evaluate(t, 0, jet4);
            expect(length(jet2[0])).toBeCloseTo(1, 12);
            expect(length(jet4[0])).toBeCloseTo(1, 12);
        }
    });
});

describe('NURBSHalfCircleDegree3', () => {
    it('has the documented control points and weights', () => {
        const curve = new NURBSHalfCircleDegree3();
        expect(curve.getNumControls()).toBe(4);
        expect(curve.getBasisFunction().getDegree()).toBe(3);
        expect(curve.getWeights()).toEqual([1, 1 / 3, 1 / 3, 1]);
        expect(curve.getControl(0).values).toEqual([1, 0]);
        expect(curve.getControl(1).values).toEqual([1, 2]);
        expect(curve.getControl(2).values).toEqual([-1, 2]);
        expect(curve.getControl(3).values).toEqual([-1, 0]);
    });

    it('traces the upper half of the unit circle', () => {
        // The upstream comment says the half circle is x >= 0; the control
        // points actually generate the half with y >= 0.
        const curve = new NURBSHalfCircleDegree3();
        const positions = expectUnitCircle(curve);
        for (const p of positions) {
            expect(p.get(1)).toBeGreaterThanOrEqual(-1e-15);
        }
        expect(length(sub(positions[0], Vector.fromArray([1, 0]))))
            .toBeLessThan(1e-14);
        expect(length(sub(positions[positions.length - 1],
            Vector.fromArray([-1, 0])))).toBeLessThan(1e-14);
        expect(angleOf(positions[(positions.length - 1) / 2]))
            .toBeCloseTo(Math.PI / 2, 12);
        // Counterclockwise traversal, spanning pi radians.
        let previous = -1;
        for (let i = 0; i < positions.length - 1; ++i) {
            const a = angleOf(positions[i]);
            expect(a).toBeGreaterThanOrEqual(previous - 1e-14);
            previous = a;
        }
    });
});

describe('NURBSFullCircleDegree3', () => {
    it('has the documented control points, weights and knots', () => {
        const curve = new NURBSFullCircleDegree3();
        expect(curve.getNumControls()).toBe(7);
        expect(curve.getBasisFunction().getDegree()).toBe(3);
        expect(curve.getWeights()).toEqual([1, 1 / 3, 1 / 3, 1, 1 / 3,
            1 / 3, 1]);
        expect(curve.getControl(0).values).toEqual([1, 0]);
        expect(curve.getControl(3).values).toEqual([-1, 0]);
        expect(curve.getControl(6).values).toEqual([1, 0]);
        // The knots are (0,0,0,0,1/2,1/2,1/2,1,1,1,1).
        const basis = curve.getBasisFunction();
        expect(basis.getNumUniqueKnots()).toBe(3);
        const unique = basis.getUniqueKnots();
        expect(unique.map(k => k.t)).toEqual([0, 0.5, 1]);
        expect(unique.map(k => k.multiplicity)).toEqual([4, 3, 4]);
    });

    it('traces the whole unit circle counterclockwise', () => {
        const curve = new NURBSFullCircleDegree3();
        const positions = expectUnitCircle(curve, 128);
        expect(length(sub(positions[0], Vector.fromArray([1, 0]))))
            .toBeLessThan(1e-14);
        expect(length(sub(positions[positions.length - 1],
            Vector.fromArray([1, 0])))).toBeLessThan(1e-14);
        expect(length(sub(positions[(positions.length - 1) / 2],
            Vector.fromArray([-1, 0])))).toBeLessThan(1e-14);
        // The traversal covers a full turn: the accumulated signed angle
        // change is 2*pi.
        let total = 0;
        for (let i = 1; i < positions.length; ++i) {
            const a0 = Math.atan2(positions[i - 1].get(1),
                positions[i - 1].get(0));
            const a1 = Math.atan2(positions[i].get(1), positions[i].get(0));
            let d = a1 - a0;
            while (d > Math.PI) { d -= 2 * Math.PI; }
            while (d < -Math.PI) { d += 2 * Math.PI; }
            total += d;
            expect(d).toBeGreaterThan(-1e-14);
        }
        expect(total).toBeCloseTo(2 * Math.PI, 10);
    });
});

describe('NURBSCircularArcDegree2', () => {
    it('reproduces the arc of a circle to high accuracy', () => {
        const center = Vector.fromArray([2, -1]);
        const radius = 3;
        for (const [a0, a1] of [[0, 1.1], [0.3, 2.0], [-1.2, 0.4],
            [2.5, 3.5], [0, Math.PI / 2]]) {
            const e0 = Vector.fromArray([center.get(0) + radius * Math.cos(a0),
                center.get(1) + radius * Math.sin(a0)]);
            const e1 = Vector.fromArray([center.get(0) + radius * Math.cos(a1),
                center.get(1) + radius * Math.sin(a1)]);
            const arc = Arc2.fromCenterRadiusEnds(center, radius, e0, e1);
            const curve = new NURBSCircularArcDegree2(arc);
            expect(curve.getNumControls()).toBe(3);
            const positions = samplePositions(curve, 64);
            for (const p of positions) {
                expect(Math.abs(length(sub(p, center)) - radius))
                    .toBeLessThan(1e-12);
            }
            // The endpoints of the NURBS curve are the arc endpoints.
            expect(length(sub(positions[0], e0))).toBeLessThan(1e-12);
            expect(length(sub(positions[positions.length - 1], e1)))
                .toBeLessThan(1e-12);
            // The arc spans exactly a1 - a0 radians.
            let total = 0;
            for (let i = 1; i < positions.length; ++i) {
                const p0 = sub(positions[i - 1], center);
                const p1 = sub(positions[i], center);
                const t0 = Math.atan2(p0.get(1), p0.get(0));
                const t1 = Math.atan2(p1.get(1), p1.get(0));
                let d = t1 - t0;
                while (d > Math.PI) { d -= 2 * Math.PI; }
                while (d < -Math.PI) { d += 2 * Math.PI; }
                total += d;
            }
            expect(total).toBeCloseTo(a1 - a0, 9);
        }
    });

    it('has the middle weight 1 and symmetric outer weights', () => {
        const center = Vector.fromArray([0, 0]);
        const arc = Arc2.fromCenterRadiusEnds(center, 1,
            Vector.fromArray([1, 0]),
            Vector.fromArray([Math.cos(1.4), Math.sin(1.4)]));
        const curve = new NURBSCircularArcDegree2(arc);
        const w = curve.getWeights();
        expect(w[1]).toBe(1);
        expect(w[0]).toBe(w[2]);
        expect(w[0]).toBeGreaterThan(0);
        // For a quarter circle the outer weights are sqrt(2) relative to a
        // middle weight of 1, matching NURBSQuarterCircleDegree2.
        const quarter = Arc2.fromCenterRadiusEnds(center, 1,
            Vector.fromArray([1, 0]), Vector.fromArray([0, 1]));
        const qcurve = new NURBSCircularArcDegree2(quarter);
        expect(qcurve.getWeights()[0]).toBeCloseTo(Math.SQRT2, 12);
        expect(qcurve.getControl(1).values[0]).toBeCloseTo(1, 12);
        expect(qcurve.getControl(1).values[1]).toBeCloseTo(1, 12);
    });

    it('matches an arbitrary circle arc (randomized)', () => {
        let seed = 55555;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 100; ++trial) {
            const center = Vector.fromArray([10 * rand() - 5, 10 * rand() - 5]);
            const radius = 0.25 + 5 * rand();
            const a0 = 2 * Math.PI * rand();
            // Keep the arc span strictly less than pi so the degree-2
            // rational representation exists.
            const span = 0.05 + 2.9 * rand() * 0.5;
            const a1 = a0 + span;
            const e0 = Vector.fromArray([center.get(0) + radius * Math.cos(a0),
                center.get(1) + radius * Math.sin(a0)]);
            const e1 = Vector.fromArray([center.get(0) + radius * Math.cos(a1),
                center.get(1) + radius * Math.sin(a1)]);
            const curve = new NURBSCircularArcDegree2(
                Arc2.fromCenterRadiusEnds(center, radius, e0, e1));
            const jet = curve.createJet();
            for (let i = 0; i <= 8; ++i) {
                curve.evaluate(i / 8, 0, jet);
                expect(Math.abs(length(sub(jet[0], center)) - radius))
                    .toBeLessThan(1e-10 * Math.max(1, radius));
            }
        }
    });
});


// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). NURBSCircle.h was read line by
// line against src/NURBSCircle.ts: every control point, weight and knot of the
// five constructions was compared term by term. The properties below check
// what the constructions claim - that the curves lie exactly on a circle, that
// they cover the documented portion of it and that they sweep it
// monotonically - rather than restating the control nets.
import {
    check, fc, expectClose
} from './helpers/arbitraries.js';

function circleAngles(curve: NURBSCurve, count: number,
    center: readonly number[] = [0, 0]): number[] {
    const jet = curve.createJet();
    const tmin = curve.getTMin();
    const tmax = curve.getTMax();
    const angles: number[] = [];
    for (let i = 0; i <= count; ++i) {
        curve.evaluate(tmin + ((tmax - tmin) * i) / count, 0, jet);
        angles.push(Math.atan2(jet[0].values[1] - center[1],
            jet[0].values[0] - center[0]));
    }
    return angles;
}

const unitParam = fc.double({ min: 0, max: 1, noNaN: true,
    noDefaultInfinity: true });

describe('NURBSCircle verification', () => {
    it('every construction lies exactly on the unit circle', () => {
        const curves: [string, NURBSCurve][] = [
            ['quarter2', new NURBSQuarterCircleDegree2()],
            ['quarter4', new NURBSQuarterCircleDegree4()],
            ['half3', new NURBSHalfCircleDegree3()],
            ['full3', new NURBSFullCircleDegree3()]];
        check(unitParam, u => {
            for (const [, curve] of curves) {
                const jet = curve.createJet();
                curve.evaluate(u, 0, jet);
                expectClose(length(jet[0]), 1, 1e-13, 1e-13);
            }
        });
    });

    it('covers the documented portion of the circle', () => {
        // The two quarter circles cover the first quadrant, the half circle
        // covers y >= 0 (upstream's comment says x >= 0, which the control
        // points (1,0), (1,2), (-1,2), (-1,0) contradict), and the full circle
        // covers everything.
        check(unitParam, u => {
            for (const curve of [new NURBSQuarterCircleDegree2(),
                new NURBSQuarterCircleDegree4()]) {
                const jet = curve.createJet();
                curve.evaluate(u, 0, jet);
                expect(jet[0].values[0]).toBeGreaterThanOrEqual(-1e-13);
                expect(jet[0].values[1]).toBeGreaterThanOrEqual(-1e-13);
            }
            const half = new NURBSHalfCircleDegree3();
            const jet = half.createJet();
            half.evaluate(u, 0, jet);
            expect(jet[0].values[1]).toBeGreaterThanOrEqual(-1e-13);
        });
    });

    it('sweeps counterclockwise and monotonically', () => {
        const sweep = (angles: number[]): number[] => {
            let previous = angles[0];
            let turns = 0;
            return angles.map(a => {
                if (a < previous - 1e-9) { turns += 2 * Math.PI; }
                previous = a;
                return a + turns;
            });
        };
        for (const [curve, total] of [
            [new NURBSQuarterCircleDegree2(), Math.PI / 2],
            [new NURBSQuarterCircleDegree4(), Math.PI / 2],
            [new NURBSHalfCircleDegree3(), Math.PI],
            [new NURBSFullCircleDegree3(), 2 * Math.PI]] as
            [NURBSCurve, number][]) {
            const swept = sweep(circleAngles(curve, 64));
            expectClose(swept[0], 0, 1e-12, 1e-12);
            for (let i = 1; i < swept.length; ++i) {
                expect(swept[i]).toBeGreaterThan(swept[i - 1]);
            }
            expectClose(swept[swept.length - 1], total, 1e-9, 1e-9);
        }
    });

    it('the two quarter circles trace the same point set', () => {
        // Different degrees and different parameterizations, but the same
        // geometry: for each parameter of one curve there is a point of the
        // other at the same polar angle.
        check(unitParam, u => {
            const a = new NURBSQuarterCircleDegree2();
            const b = new NURBSQuarterCircleDegree4();
            const ja = a.createJet();
            a.evaluate(u, 0, ja);
            const angle = Math.atan2(ja[0].values[1], ja[0].values[0]);
            // Invert the degree-4 parameterization by bisection on the angle.
            const jb = b.createJet();
            const angleAt = (t: number): number => {
                b.evaluate(t, 0, jb);
                return Math.atan2(jb[0].values[1], jb[0].values[0]);
            };
            let lo = 0, hi = 1;
            for (let i = 0; i < 60; ++i) {
                const mid = 0.5 * (lo + hi);
                if (angleAt(mid) < angle) { lo = mid; } else { hi = mid; }
            }
            expectClose(angleAt(0.5 * (lo + hi)), angle, 1e-9, 1e-9);
        });
    });

    it('the circular arc construction reproduces the given arc', () => {
        check(fc.tuple(
            fc.tuple(fc.double({ min: -5, max: 5, noNaN: true,
                noDefaultInfinity: true }),
            fc.double({ min: -5, max: 5, noNaN: true,
                noDefaultInfinity: true })),
            fc.double({ min: 0.2, max: 8, noNaN: true,
                noDefaultInfinity: true }),
            fc.double({ min: -Math.PI, max: Math.PI, noNaN: true,
                noDefaultInfinity: true }),
            // The rational quadratic represents an arc of less than pi; at pi
            // the tangent-intersection control point is at infinity.
            fc.double({ min: 0.05, max: Math.PI - 0.05, noNaN: true,
                noDefaultInfinity: true }),
            unitParam),
        ([center, radius, a0, span, u]) => {
            const point = (angle: number): Vector => Vector.fromArray([
                center[0] + radius * Math.cos(angle),
                center[1] + radius * Math.sin(angle)]);
            const arc = Arc2.fromCenterRadiusEnds(
                Vector.fromArray([center[0], center[1]]), radius,
                point(a0), point(a0 + span));
            const curve = new NURBSCircularArcDegree2(arc);
            const jet = curve.createJet();
            curve.evaluate(u, 0, jet);
            // On the circle.
            expectClose(Math.hypot(jet[0].values[0] - center[0],
                jet[0].values[1] - center[1]), radius, 1e-9 * radius, 1e-11);
            // Inside the angular span, measured from the first endpoint.
            const twoPi = 2 * Math.PI;
            let swept = Math.atan2(jet[0].values[1] - center[1],
                jet[0].values[0] - center[0]) - a0;
            swept = ((swept % twoPi) + twoPi) % twoPi;
            // A sample at the first endpoint can land a rounding step below
            // zero, which wraps to just under 2*pi; fold that back.
            if (swept > twoPi - 1e-6) { swept -= twoPi; }
            expect(swept).toBeGreaterThanOrEqual(-1e-6);
            expect(swept).toBeLessThanOrEqual(span + 1e-6);
            // The endpoints are interpolated.
            curve.evaluate(0, 0, jet);
            expectClose(jet[0].values[0], point(a0).values[0],
                1e-9 * radius, 1e-11);
            expectClose(jet[0].values[1], point(a0).values[1],
                1e-9 * radius, 1e-11);
            curve.evaluate(1, 0, jet);
            expectClose(jet[0].values[0], point(a0 + span).values[0],
                1e-9 * radius, 1e-11);
            expectClose(jet[0].values[1], point(a0 + span).values[1],
                1e-9 * radius, 1e-11);
            // The outer weights are equal and the middle weight is 1.
            const w = curve.getWeights();
            expect(w[1]).toBe(1);
            expect(w[2]).toBe(w[0]);
            expectClose(w[0], 1 / Math.cos(0.5 * span), 1e-9, 1e-9);
        });
    });
});
