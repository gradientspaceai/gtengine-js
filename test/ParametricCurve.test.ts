import { describe, expect, it } from 'vitest';
import { ParametricCurve } from '../src/ParametricCurve.js';
import { Vector, length as vectorLength, sub } from '../src/Vector.js';

// A circle of radius r traversed once on [0, 2*pi]. Speed is the constant r,
// so every quantity computed by ParametricCurve has a closed form:
//   length(t0,t1) = r*(t1-t0),  totalLength = 2*pi*r,  time(s) = s/r.
class CircleCurve extends ParametricCurve {
    constructor(private radius: number) {
        super(2, 0, 2 * Math.PI);
        this.mConstructed = true;
    }

    evaluate(t: number, order: number, jet: Vector[]): void {
        const r = this.radius;
        const c = Math.cos(t);
        const s = Math.sin(t);
        jet[0].set(0, r * c);
        jet[0].set(1, r * s);
        if (order >= 1) {
            jet[1].set(0, -r * s);
            jet[1].set(1, r * c);
        }
        if (order >= 2) {
            jet[2].set(0, -r * c);
            jet[2].set(1, -r * s);
        }
        if (order >= 3) {
            jet[3].set(0, r * s);
            jet[3].set(1, -r * c);
        }
    }
}

// A piecewise-linear curve with one segment per consecutive vertex pair and
// knot times 0, 1, 2, ..., n-1. The speed is constant on each segment, so
// Romberg integration is exact there and segment lengths are known.
class PolylineCurve extends ParametricCurve {
    constructor(private points: Vector[]) {
        super(2, points.length - 1,
            Array.from({ length: points.length }, (_v, i) => i));
        this.mConstructed = true;
    }

    evaluate(t: number, order: number, jet: Vector[]): void {
        const numSegments = this.points.length - 1;
        let i = Math.floor(t);
        if (i < 0) {
            i = 0;
        }
        if (i > numSegments - 1) {
            i = numSegments - 1;
        }
        const u = t - i;
        const p = this.points[i];
        const d = sub(this.points[i + 1], p);
        jet[0].set(0, p.get(0) + u * d.get(0));
        jet[0].set(1, p.get(1) + u * d.get(1));
        if (order >= 1) {
            jet[1].set(0, d.get(0));
            jet[1].set(1, d.get(1));
        }
    }
}

// The same circle, but built with the multiple-segment constructor form
// (four quarter-circle segments). The speed is still the constant r, so the
// closed forms above still hold and the segment-crossing code paths of
// getLength/getTime can be checked against them.
class MultiSegmentCircle extends ParametricCurve {
    constructor(private radius: number,
        times: readonly number[] = [0, Math.PI / 2, Math.PI,
            (3 * Math.PI) / 2, 2 * Math.PI]) {
        super(2, times.length - 1, times);
        this.mConstructed = true;
    }

    evaluate(t: number, order: number, jet: Vector[]): void {
        const r = this.radius;
        jet[0].set(0, r * Math.cos(t));
        jet[0].set(1, r * Math.sin(t));
        if (order >= 1) {
            jet[1].set(0, -r * Math.sin(t));
            jet[1].set(1, r * Math.cos(t));
        }
    }
}

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

describe('ParametricCurve construction and member access', () => {
    it('reports the single-segment interval and dimension', () => {
        const curve = new CircleCurve(3);
        expect(curve.getDimension()).toBe(2);
        expect(curve.getTMin()).toBe(0);
        expect(curve.getTMax()).toBeCloseTo(2 * Math.PI, 15);
        expect(curve.getNumSegments()).toBe(1);
        expect(curve.isConstructed()).toBe(true);
        expect(Array.from(curve.getTimes())).toEqual([0, 2 * Math.PI]);
    });

    it('reports the multiple-segment times', () => {
        const curve = new PolylineCurve([v2(0, 0), v2(3, 0), v2(3, 4), v2(6, 8)]);
        expect(curve.getNumSegments()).toBe(3);
        expect(Array.from(curve.getTimes())).toEqual([0, 1, 2, 3]);
        expect(curve.getTMin()).toBe(0);
        expect(curve.getTMax()).toBe(3);
    });

    it('setTimeInterval applies only to two-time curves', () => {
        const single = new CircleCurve(1);
        single.setTimeInterval(-1, 5);
        expect(single.getTMin()).toBe(-1);
        expect(single.getTMax()).toBe(5);

        const multi = new PolylineCurve([v2(0, 0), v2(1, 0), v2(1, 1)]);
        multi.setTimeInterval(-10, 10);
        expect(Array.from(multi.getTimes())).toEqual([0, 1, 2]);
    });

    it('clamps the Romberg order and bisection count to at least one', () => {
        const curve = new CircleCurve(1);
        // These have no observable accessor upstream either; exercise them to
        // confirm they do not disturb the computed values.
        curve.setRombergOrder(0);
        curve.setMaxBisections(0);
        expect(curve.getLength(0, Math.PI)).toBeGreaterThan(0);
    });

    it('createJet allocates SUP_ORDER zero vectors of the curve dimension', () => {
        const curve = new CircleCurve(1);
        const jet = curve.createJet();
        expect(ParametricCurve.SUP_ORDER).toBe(4);
        expect(jet.length).toBe(4);
        for (const j of jet) {
            expect(j.size).toBe(2);
            expect(j.values).toEqual([0, 0]);
        }
    });
});

describe('ParametricCurve differential geometry (circle)', () => {
    const r = 2.5;
    const curve = new CircleCurve(r);

    it('getPosition matches the parameterization', () => {
        for (const t of [0, 0.3, 1.7, Math.PI, 5.9]) {
            const p = curve.getPosition(t);
            expect(p.get(0)).toBeCloseTo(r * Math.cos(t), 12);
            expect(p.get(1)).toBeCloseTo(r * Math.sin(t), 12);
        }
    });

    it('getTangent is unit length and orthogonal to the position', () => {
        for (const t of [0.1, 1.2, 2.4, 4.8]) {
            const tangent = curve.getTangent(t);
            expect(vectorLength(tangent)).toBeCloseTo(1, 12);
            const p = curve.getPosition(t);
            expect(p.get(0) * tangent.get(0) + p.get(1) * tangent.get(1))
                .toBeCloseTo(0, 12);
        }
    });

    it('getSpeed is the constant radius', () => {
        for (const t of [0, 1, 2, 3, 6]) {
            expect(curve.getSpeed(t)).toBeCloseTo(r, 12);
        }
    });

    it('getPosition returns a copy, not the internal jet entry', () => {
        const p0 = curve.getPosition(1);
        p0.set(0, 12345);
        const p1 = curve.getPosition(1);
        expect(p1.get(0)).toBeCloseTo(r * Math.cos(1), 12);
    });
});

describe('ParametricCurve arc length', () => {
    it('total length of a circle is 2*pi*r', () => {
        const r = 1.75;
        const curve = new CircleCurve(r);
        expect(curve.getTotalLength()).toBeCloseTo(2 * Math.PI * r, 10);
    });

    it('partial lengths of a circle are r*(t1-t0)', () => {
        const r = 3;
        const curve = new CircleCurve(r);
        expect(curve.getLength(0, Math.PI)).toBeCloseTo(r * Math.PI, 10);
        expect(curve.getLength(1, 2)).toBeCloseTo(r, 10);
        expect(curve.getLength(0.25, 5.5)).toBeCloseTo(r * 5.25, 10);
    });

    it('clamps the requested interval to [tmin, tmax]', () => {
        const r = 1;
        const curve = new CircleCurve(r);
        expect(curve.getLength(-100, 100)).toBeCloseTo(2 * Math.PI * r, 10);
    });

    it('a degenerate interval has zero length', () => {
        const curve = new CircleCurve(1);
        expect(curve.getLength(2, 2)).toBe(0);
    });

    it('sums whole and partial segments of a multi-segment curve', () => {
        // Speed 2 with knots at t = 0, 1, 3, 6, so the segment lengths are
        // 2, 4 and 6 and the total is 12.
        const curve = new MultiSegmentCircle(2, [0, 1, 3, 6]);
        expect(curve.getNumSegments()).toBe(3);
        expect(curve.getTotalLength()).toBeCloseTo(12, 10);
        expect(curve.getLength(0, 1)).toBeCloseTo(2, 10);
        expect(curve.getLength(1, 6)).toBeCloseTo(10, 10);
        // Crossing two knots: 2*(0.5) + 2*(2) + 2*(1.5) = 8.
        expect(curve.getLength(0.5, 4.5)).toBeCloseTo(8, 10);
        // Wholly inside one segment.
        expect(curve.getLength(3.25, 3.75)).toBeCloseTo(1, 10);
        // Ending exactly on a knot skips the trailing partial integration.
        expect(curve.getLength(0.5, 3)).toBeCloseTo(5, 10);
    });

    it('agrees with the single-segment result across segment boundaries', () => {
        const r = 2;
        const one = new CircleCurve(r);
        const many = new MultiSegmentCircle(r);
        expect(many.getNumSegments()).toBe(4);
        expect(many.getTotalLength()).toBeCloseTo(2 * Math.PI * r, 10);
        for (const [t0, t1] of [[0, 2], [0.5, 5.5], [1, 1.2], [2.5, 4.5],
            [0, 2 * Math.PI]]) {
            expect(many.getLength(t0, t1)).toBeCloseTo(one.getLength(t0, t1), 10);
            expect(many.getLength(t0, t1)).toBeCloseTo(r * (t1 - t0), 10);
        }
        // Exactly on a knot: the partial-segment branch is skipped.
        expect(many.getLength(Math.PI / 2, Math.PI)).toBeCloseTo(r * Math.PI / 2, 10);
    });

    it('getTime inverts the length on a multi-segment curve', () => {
        const r = 2;
        const many = new MultiSegmentCircle(r);
        for (const s of [0.5, 3, 7, 12]) {
            expect(many.getTime(s)).toBeCloseTo(s / r, 8);
        }
    });

    it('getTotalLength is cached after the segment lengths are computed', () => {
        const curve = new MultiSegmentCircle(2, [0, 1, 3]);
        const first = curve.getTotalLength();
        const second = curve.getTotalLength();
        expect(second).toBe(first);
        expect(second).toBeCloseTo(6, 10);
    });
});

describe('ParametricCurve inverse arc length', () => {
    const r = 2;
    const curve = new CircleCurve(r);

    it('maps nonpositive lengths to tmin and overlong lengths to tmax', () => {
        expect(curve.getTime(0)).toBe(curve.getTMin());
        expect(curve.getTime(-5)).toBe(curve.getTMin());
        expect(curve.getTime(1000)).toBe(curve.getTMax());
    });

    it('inverts the arc length of a circle, t = s/r', () => {
        for (const s of [0.5, 1, 3, 6, 11]) {
            expect(curve.getTime(s)).toBeCloseTo(s / r, 8);
        }
    });

    it('getTime is the inverse of getLength', () => {
        for (const t of [0.4, 1.9, 3.3, 5.7]) {
            const s = curve.getLength(curve.getTMin(), t);
            expect(curve.getTime(s)).toBeCloseTo(t, 8);
        }
    });
});

describe('ParametricCurve subdivision', () => {
    it('subdivideByTime uses uniform parameter spacing', () => {
        const curve = new CircleCurve(1);
        const points = curve.subdivideByTime(5);
        expect(points.length).toBe(5);
        for (let i = 0; i < 5; ++i) {
            const t = (2 * Math.PI * i) / 4;
            expect(points[i].get(0)).toBeCloseTo(Math.cos(t), 12);
            expect(points[i].get(1)).toBeCloseTo(Math.sin(t), 12);
        }
    });

    it('subdivideByLength is uniform in arc length for a circle', () => {
        const r = 1.5;
        const curve = new CircleCurve(r);
        const numPoints = 9;
        const points = curve.subdivideByLength(numPoints);
        expect(points.length).toBe(numPoints);

        // For a constant-speed curve, equal arc length means equal parameter
        // spacing, so consecutive chords all have the same length.
        const expected = 2 * r * Math.sin(Math.PI / (numPoints - 1));
        for (let i = 1; i < numPoints; ++i) {
            expect(vectorLength(sub(points[i], points[i - 1])))
                .toBeCloseTo(expected, 7);
        }
        // The first and last points are the endpoints of the interval.
        expect(points[0].get(0)).toBeCloseTo(r, 10);
        expect(points[numPoints - 1].get(0)).toBeCloseTo(r, 7);
    });

    it('subdivideByLength on a multi-segment curve is uniform', () => {
        // Four quarter-circle segments; the speed is the constant r on all of
        // them, so the arc-length subdivision is uniform in the parameter.
        const r = 2;
        const curve = new MultiSegmentCircle(r);
        const numPoints = 13;
        const points = curve.subdivideByLength(numPoints);
        expect(points.length).toBe(numPoints);
        const expected = 2 * r * Math.sin(Math.PI / (numPoints - 1));
        for (let i = 1; i < numPoints; ++i) {
            expect(vectorLength(sub(points[i], points[i - 1])))
                .toBeCloseTo(expected, 7);
        }
    });

    it('subdivideByTime with two points returns the endpoints', () => {
        const curve = new PolylineCurve([v2(0, 0), v2(3, 0), v2(3, 4)]);
        const points = curve.subdivideByTime(2);
        expect(points.length).toBe(2);
        expect(points[0].values).toEqual([0, 0]);
        expect(points[1].values).toEqual([3, 4]);
    });
});

describe('ParametricCurve randomized cross-check', () => {
    it('arc length agrees with a fine polygonal approximation', () => {
        // An ellipse has no closed-form arc length; compare the Romberg
        // result against a dense chord sum, which converges from below.
        class EllipseCurve extends ParametricCurve {
            constructor(private a: number, private b: number) {
                super(2, 0, 2 * Math.PI);
                this.mConstructed = true;
            }

            evaluate(t: number, order: number, jet: Vector[]): void {
                jet[0].set(0, this.a * Math.cos(t));
                jet[0].set(1, this.b * Math.sin(t));
                if (order >= 1) {
                    jet[1].set(0, -this.a * Math.sin(t));
                    jet[1].set(1, this.b * Math.cos(t));
                }
            }
        }

        let seed = 987654321;
        const rand = (): number => {
            seed = (1664525 * seed + 1013904223) >>> 0;
            return seed / 4294967296;
        };

        for (let trial = 0; trial < 3; ++trial) {
            const a = 0.5 + 2 * rand();
            const b = 0.5 + 2 * rand();
            const curve = new EllipseCurve(a, b);
            curve.setRombergOrder(12);

            const t1 = 2 * Math.PI * rand();
            const numChords = 40000;
            let chordSum = 0;
            let prev = curve.getPosition(0);
            for (let i = 1; i <= numChords; ++i) {
                const p = curve.getPosition((t1 * i) / numChords);
                chordSum += vectorLength(sub(p, prev));
                prev = p;
            }

            expect(curve.getLength(0, t1)).toBeCloseTo(chordSum, 6);
        }
    });
});
