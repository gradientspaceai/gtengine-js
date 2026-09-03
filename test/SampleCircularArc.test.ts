import { describe, it, expect } from 'vitest';
import { SampleCircularArc } from '../src/SampleCircularArc.js';
import { Arc2 } from '../src/Arc2.js';
import { Vector } from '../src/Vector.js';

function arcFromAngle(center: readonly number[], radius: number,
    a0: number, a1: number): Arc2 {
    return Arc2.fromCenterRadiusEnds(
        Vector.fromArray([center[0], center[1]]), radius,
        Vector.fromArray([center[0] + radius * Math.cos(a0),
            center[1] + radius * Math.sin(a0)]),
        Vector.fromArray([center[0] + radius * Math.cos(a1),
            center[1] + radius * Math.sin(a1)]));
}

function angleOf(p: Vector, center: readonly number[]): number {
    return Math.atan2(p.values[1] - center[1], p.values[0] - center[0]);
}

// The CCW angles measured from a0, unwrapped into [0, 2*pi).
function sweptAngles(points: readonly Vector[], center: readonly number[],
    a0: number): number[] {
    return points.map(p => {
        let a = angleOf(p, center) - a0;
        while (a < -1e-9) { a += 2 * Math.PI; }
        return a;
    });
}

describe('SampleCircularArc', () => {
    const sampler = new SampleCircularArc();

    it('samples floor(r*angle) points for an arc of angle <= pi/2', () => {
        const r = 100;
        const arc = arcFromAngle([0, 0], r, 0, Math.PI / 2);
        const points = sampler.compute(arc);
        expect(points.length).toBe(Math.trunc(r * Math.PI / 2));
    });

    it('samples 2*floor(r*angle/2) points for an arc in (pi/2, pi]', () => {
        const r = 100;
        const angle = 0.9 * Math.PI;
        const arc = arcFromAngle([0, 0], r, 0, angle);
        const points = sampler.compute(arc);
        expect(points.length).toBe(2 * Math.trunc(r * angle / 2));
    });

    it('samples 3*floor(r*angle/3) points for an arc in (pi, 3*pi/2]', () => {
        const r = 100;
        const angle = 3 * Math.PI / 2;
        const arc = arcFromAngle([0, 0], r, 0, angle);
        const points = sampler.compute(arc);
        expect(points.length).toBe(3 * Math.trunc(r * angle / 3));
    });

    it('samples 4*floor(r*angle/4) points for an arc in (3*pi/2, 2*pi)', () => {
        const r = 100;
        const angle = 7 * Math.PI / 4;
        const arc = arcFromAngle([0, 0], r, 0, angle);
        const points = sampler.compute(arc);
        expect(points.length).toBe(4 * Math.trunc(r * angle / 4));
    });

    it('returns no points when the arc is too short to sample', () => {
        const arc = arcFromAngle([0, 0], 1, 0, 0.001);
        expect(sampler.compute(arc).length).toBe(0);
    });

    it('places every sample exactly on the circle', () => {
        for (const angle of [Math.PI / 6, Math.PI / 2, 0.75 * Math.PI,
            Math.PI, 1.25 * Math.PI, 1.9 * Math.PI]) {
            const center = [3, -7];
            const r = 50;
            const arc = arcFromAngle(center, r, 0.4, 0.4 + angle);
            const points = sampler.compute(arc);
            expect(points.length).toBeGreaterThan(0);
            for (const p of points) {
                const d = Math.hypot(p.values[0] - center[0],
                    p.values[1] - center[1]);
                expect(d).toBeCloseTo(r, 8);
            }
        }
    });

    it('starts at the first arc endpoint and excludes the second', () => {
        const center = [0, 0];
        const r = 40;
        const arc = arcFromAngle(center, r, 0, Math.PI / 3);
        const points = sampler.compute(arc);
        // u = i/numPoints for 0 <= i < numPoints, so the arc's first endpoint
        // is sample 0 and the second endpoint is not produced.
        expect(points[0].values[0]).toBeCloseTo(r, 10);
        expect(points[0].values[1]).toBeCloseTo(0, 10);
        const last = points[points.length - 1];
        const lastAngle = angleOf(last, center);
        expect(lastAngle).toBeLessThan(Math.PI / 3);
        expect(lastAngle).toBeGreaterThan(Math.PI / 3
            - 2 * (Math.PI / 3) / points.length);
    });

    it('sweeps monotonically for arcs of angle <= pi', () => {
        const center = [-2, 5];
        for (const angle of [0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]) {
            const a0 = 1.1;
            const arc = arcFromAngle(center, 80, a0, a0 + angle);
            const points = sampler.compute(arc);
            const swept = sweptAngles(points, center, a0);
            expect(swept[0]).toBeCloseTo(0, 8);
            for (let i = 1; i < swept.length; ++i) {
                expect(swept[i]).toBeGreaterThan(swept[i - 1]);
            }
            expect(swept[swept.length - 1]).toBeLessThan(angle + 1e-9);
        }
    });

    it('subdivides an arc in (pi/2, pi] at its exact bisector', () => {
        // SampleArc2 splits at the normalized P0 + P2, which is the true
        // angular bisector of the arc.
        const center = [0, 0];
        const angle = 0.8 * Math.PI;
        const arc = arcFromAngle(center, 100, 0, angle);
        const points = sampler.compute(arc);
        const half = points.length / 2;
        const swept = sweptAngles(points, center, 0);
        expect(swept[half]).toBeCloseTo(angle / 2, 8);
    });

    it('handles a half circle without degenerating', () => {
        const center = [0, 0];
        const arc = arcFromAngle(center, 100, 0, Math.PI);
        const points = sampler.compute(arc);
        expect(points.length).toBe(2 * Math.trunc(100 * Math.PI / 2));
        for (const p of points) {
            expect(Number.isFinite(p.values[0])).toBe(true);
            expect(Math.hypot(p.values[0], p.values[1])).toBeCloseTo(100, 8);
        }
    });

    it('scales with the radius and translates with the center', () => {
        const points0 = sampler.compute(arcFromAngle([0, 0], 60, 0.3, 1.2));
        const points1 = sampler.compute(arcFromAngle([10, -4], 60, 0.3, 1.2));
        expect(points1.length).toBe(points0.length);
        for (let i = 0; i < points0.length; ++i) {
            expect(points1[i].values[0]).toBeCloseTo(points0[i].values[0] + 10, 9);
            expect(points1[i].values[1]).toBeCloseTo(points0[i].values[1] - 4, 9);
        }
    });

    // Upstream defect, preserved by the port: SampleArc3 and SampleArc4 build
    // their subdivision directions as -(2*P0+P2), -(P0+2*P2) (and the
    // 3:1/1:1/1:3 analogues), which are exact only for the 1:1 midpoint. For
    // arcs longer than pi the other split directions are not the intended
    // trisectors/quadsectors, so the subarcs are not traversed in order. The
    // samples still lie on the circle, which is what this test pins down.
    it('preserves the upstream subdivision for arcs longer than pi', () => {
        const center = [0, 0];
        const arc = arcFromAngle(center, 100, 0, 3 * Math.PI / 2);
        const points = sampler.compute(arc);
        for (const p of points) {
            expect(Math.hypot(p.values[0], p.values[1])).toBeCloseTo(100, 8);
        }
        const swept = sweptAngles(points, center, 0);
        const third = points.length / 3;
        // The 1:2 split direction is at ~153.4 degrees rather than the exact
        // trisector at 90 degrees.
        expect(swept[third] * 180 / Math.PI).toBeCloseTo(153.434948, 4);
        // Consequently, the sweep is not monotone.
        let monotone = true;
        for (let i = 1; i < swept.length; ++i) {
            if (swept[i] <= swept[i - 1]) { monotone = false; break; }
        }
        expect(monotone).toBe(false);
    });

    it('preserves the upstream 1:1 split direction for arcs longer than pi',
        () => {
            // The middle quadsector, -(P0+P2), is the exact bisector even for
            // arcs longer than pi.
            const center = [0, 0];
            const angle = 7 * Math.PI / 4;
            const arc = arcFromAngle(center, 100, 0, angle);
            const points = sampler.compute(arc);
            const quarter = points.length / 4;
            const swept = sweptAngles(points, center, 0);
            expect(swept[2 * quarter]).toBeCloseTo(angle / 2, 6);
        });
});
