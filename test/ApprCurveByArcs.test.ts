import { describe, expect, it } from 'vitest';
import { approximateCurveByArcs } from '../src/ApprCurveByArcs';
import { ParametricCurve } from '../src/ParametricCurve';
import { Vector, length as vectorLength, sub } from '../src/Vector';

// A circle of radius r centered at (cx,cy), traversed once on [0, 2*pi].
// The speed is the constant r, so the arc-length subdivision is exact.
class CircleCurve extends ParametricCurve {
    constructor(private radius: number, private cx: number = 0,
        private cy: number = 0) {
        super(2, 0, 2 * Math.PI);
        this.mConstructed = true;
    }

    evaluate(t: number, order: number, jet: Vector[]): void {
        const r = this.radius;
        const c = Math.cos(t);
        const s = Math.sin(t);
        jet[0].set(0, this.cx + r * c);
        jet[0].set(1, this.cy + r * s);
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

// A straight line segment from P0 to P1 on t in [0,1].
class LineCurve extends ParametricCurve {
    constructor(private p0: Vector, private p1: Vector) {
        super(2, 0, 1);
        this.mConstructed = true;
    }

    evaluate(t: number, order: number, jet: Vector[]): void {
        const d = sub(this.p1, this.p0);
        jet[0].set(0, this.p0.get(0) + t * d.get(0));
        jet[0].set(1, this.p0.get(1) + t * d.get(1));
        if (order >= 1) {
            jet[1].set(0, d.get(0));
            jet[1].set(1, d.get(1));
        }
        for (let i = 2; i <= order; ++i) {
            jet[i].set(0, 0);
            jet[i].set(1, 0);
        }
    }
}

// The Archimedean spiral X(t) = ((a + b*t)*cos(t), (a + b*t)*sin(t)), a
// curve with a genuinely varying curvature.
class SpiralCurve extends ParametricCurve {
    constructor(private a: number, private b: number, tmax: number) {
        super(2, 0, tmax);
        this.mConstructed = true;
    }

    evaluate(t: number, order: number, jet: Vector[]): void {
        const r = this.a + this.b * t;
        const c = Math.cos(t);
        const s = Math.sin(t);
        jet[0].set(0, r * c);
        jet[0].set(1, r * s);
        if (order >= 1) {
            jet[1].set(0, this.b * c - r * s);
            jet[1].set(1, this.b * s + r * c);
        }
        if (order >= 2) {
            jet[2].set(0, -2 * this.b * s - r * c);
            jet[2].set(1, 2 * this.b * c - r * s);
        }
        if (order >= 3) {
            jet[3].set(0, -3 * this.b * c + r * s);
            jet[3].set(1, -3 * this.b * s - r * c);
        }
    }
}

describe('approximateCurveByArcs', () => {
    it('rejects invalid input', () => {
        const curve = new CircleCurve(1);
        expect(() => approximateCurveByArcs(curve, 0)).toThrowError('Invalid input.');
        expect(() => approximateCurveByArcs(curve, -3)).toThrowError('Invalid input.');
    });

    it('produces the requested container sizes', () => {
        const numArcs = 5;
        const result = approximateCurveByArcs(new CircleCurve(1), numArcs);
        expect(result.times.length).toBe(2 * numArcs + 1);
        expect(result.points.length).toBe(2 * numArcs + 1);
        expect(result.arcs.length).toBe(numArcs);
    });

    it('reproduces a circle nearly exactly', () => {
        const radius = 2.5;
        const cx = -1, cy = 3;
        const numArcs = 6;
        const result = approximateCurveByArcs(
            new CircleCurve(radius, cx, cy), numArcs);

        // Every arc's circumscribed circle is the curve's own circle.
        for (const arc of result.arcs) {
            expect(arc.center.get(0)).toBeCloseTo(cx, 6);
            expect(arc.center.get(1)).toBeCloseTo(cy, 6);
            expect(arc.radius).toBeCloseTo(radius, 6);
        }

        // Each of the 2*N+1 sample points is on the circle.
        for (const p of result.points) {
            expect(vectorLength(sub(p, Vector.fromArray([cx, cy]))))
                .toBeCloseTo(radius, 6);
        }

        // The arc-length subdivision is uniform: for a circle, equal arc
        // lengths mean equal parameter increments.
        const dt = (2 * Math.PI) / (2 * numArcs);
        for (let i = 0; i <= 2 * numArcs; i += 2) {
            expect(result.times[i]).toBeCloseTo(dt * i, 5);
        }
    });

    it('places the midpoint of each arc on the bisector of its chord', () => {
        const result = approximateCurveByArcs(new SpiralCurve(1, 0.4, 6), 5);
        for (let i = 0, j0 = 0, j1 = 1, j2 = 2; i < 5;
            ++i, j0 += 2, j1 += 2, j2 += 2) {
            const P0 = result.points[j0];
            const M = result.points[j1];
            const P1 = result.points[j2];
            // |M-P0| = |M-P1| is what the bisection root guarantees.
            expect(vectorLength(sub(M, P0)))
                .toBeCloseTo(vectorLength(sub(M, P1)), 8);

            // The arc parameter is between its endpoints' parameters.
            expect(result.times[j1]).toBeGreaterThan(result.times[j0]);
            expect(result.times[j1]).toBeLessThan(result.times[j2]);
        }
    });

    it('produces a C0-continuous chain of arcs through the curve samples', () => {
        const numArcs = 7;
        const curve = new SpiralCurve(0.5, 0.3, 8);
        const result = approximateCurveByArcs(curve, numArcs);

        for (let i = 0; i < numArcs; ++i) {
            const arc = result.arcs[i];
            expect(arc.end[0].get(0)).toBe(result.points[2 * i].get(0));
            expect(arc.end[0].get(1)).toBe(result.points[2 * i].get(1));
            expect(arc.end[1].get(0)).toBe(result.points[2 * i + 2].get(0));
            expect(arc.end[1].get(1)).toBe(result.points[2 * i + 2].get(1));

            // Consecutive arcs share an endpoint.
            if (i + 1 < numArcs) {
                expect(arc.end[1].get(0)).toBe(result.arcs[i + 1].end[0].get(0));
                expect(arc.end[1].get(1)).toBe(result.arcs[i + 1].end[0].get(1));
            }

            // The three fitting points are equidistant from the arc center.
            const M = result.points[2 * i + 1];
            expect(vectorLength(sub(arc.end[0], arc.center)))
                .toBeCloseTo(arc.radius, 6);
            expect(vectorLength(sub(arc.end[1], arc.center)))
                .toBeCloseTo(arc.radius, 6);
            expect(vectorLength(sub(M, arc.center))).toBeCloseTo(arc.radius, 6);
        }

        // The times are increasing and span the curve's parameter interval.
        for (let i = 1; i < result.times.length; ++i) {
            expect(result.times[i]).toBeGreaterThan(result.times[i - 1]);
        }
        expect(result.times[0]).toBeCloseTo(curve.getTMin(), 8);
        expect(result.times[result.times.length - 1])
            .toBeCloseTo(curve.getTMax(), 8);
    });

    it('marks colinear triples as line segments', () => {
        // A straight curve makes every {P0,M,P1} colinear, so with a
        // positive epsilon threshold the arcs degenerate to segments.
        const p0 = Vector.fromArray([-2, 1]);
        const p1 = Vector.fromArray([6, 5]);
        const result = approximateCurveByArcs(new LineCurve(p0, p1), 4, 1e-12);

        for (const arc of result.arcs) {
            expect(arc.radius).toBe(Number.MAX_VALUE);
            expect(arc.center.get(0)).toBe(Number.MAX_VALUE);
            expect(arc.center.get(1)).toBe(Number.MAX_VALUE);
        }

        // The endpoints still subdivide the segment uniformly by length.
        for (let i = 0; i < 4; ++i) {
            const arc = result.arcs[i];
            expect(arc.end[0].get(0)).toBeCloseTo(-2 + 8 * (i / 4), 6);
            expect(arc.end[0].get(1)).toBeCloseTo(1 + 4 * (i / 4), 6);
        }
    });

    it('increases fidelity as the arc count grows', () => {
        // Measure the maximum distance from arc midpoints of the chords to
        // the true curve; more arcs must not make the approximation worse.
        const curve = new SpiralCurve(1, 0.25, 5);
        const errorFor = (numArcs: number): number => {
            const result = approximateCurveByArcs(curve, numArcs);
            let maxError = 0;
            for (let i = 0; i < numArcs; ++i) {
                const arc = result.arcs[i];
                // Sample the arc at its own midpoint parameter and compare
                // with the chord midpoint of the subcurve.
                const t0 = result.times[2 * i];
                const t2 = result.times[2 * i + 2];
                for (let k = 1; k < 8; ++k) {
                    const t = t0 + ((t2 - t0) * k) / 8;
                    const X = curve.getPosition(t);
                    const d = Math.abs(
                        vectorLength(sub(X, arc.center)) - arc.radius);
                    maxError = Math.max(maxError, d);
                }
            }
            return maxError;
        };

        const e4 = errorFor(4);
        const e16 = errorFor(16);
        expect(e16).toBeLessThan(e4);
        expect(e16).toBeLessThan(1e-3);
    });
});
