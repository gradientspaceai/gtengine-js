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


// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). SampleCircularArc.h was read
// line by line against src/SampleCircularArc.ts. SampleAcuteArc evaluates a
// rational quadratic Bezier with weights (w0, 1, w0), control points P0, the
// tangent intersection P1 and P2, and w0 = 1/cos(phi) for a subarc of
// half-angle phi. Solving that parameterization in closed form gives the
// sample at parameter u the polar angle
//     psi(u) = phi + 2 * atan((2u - 1) * tan(phi/2))
// measured from the start of the subarc. That identity comes from the geometry
// of the rational quadratic circle, not from the code, and is what the first
// two properties below check.
//
// The subtended angle the class works with is recovered from the endpoints by
// acos of a clamped dot product, so it is not bit-identical to the angle a
// test constructs the arc from (near a half circle acos is ill conditioned).
// The helper below reproduces that trivial recovery step so the properties
// compare against the parameterization the class actually uses.
import {
    check, fc, expectClose
} from './helpers/arbitraries.js';

function circlePoint(center: readonly number[], radius: number,
    angle: number): number[] {
    return [center[0] + radius * Math.cos(angle),
        center[1] + radius * Math.sin(angle)];
}

// The unit direction vectors of the arc endpoints, as the class computes them.
function endDirections(arc: Arc2): { P0: number[], P2: number[] } {
    const c = arc.center.values;
    const r = arc.radius;
    return {
        P0: [(arc.end[0].values[0] - c[0]) / r,
            (arc.end[0].values[1] - c[1]) / r],
        P2: [(arc.end[1].values[0] - c[0]) / r,
            (arc.end[1].values[1] - c[1]) / r]
    };
}

function subtendedAngle(arc: Arc2): number {
    const { P0, P2 } = endDirections(arc);
    const d = Math.max(-1, Math.min(P0[0] * P2[0] + P0[1] * P2[1], 1));
    return Math.acos(d);
}

// The closed-form samples of one acute subarc going from unit direction A to
// unit direction B.
function acuteSamples(center: readonly number[], radius: number,
    A: readonly number[], B: readonly number[],
    numPoints: number): number[][] {
    const d = Math.max(-1, Math.min(A[0] * B[0] + A[1] * B[1], 1));
    const phi = 0.5 * Math.acos(d);
    const halfTan = Math.tan(0.5 * phi);
    const start = Math.atan2(A[1], A[0]);
    const out: number[][] = [];
    for (let i = 0; i < numPoints; ++i) {
        const u = i / numPoints;
        const psi = phi + 2 * Math.atan((2 * u - 1) * halfTan);
        out.push(circlePoint(center, radius, start + psi));
    }
    return out;
}

function normalized(v: readonly number[]): number[] {
    const len = Math.hypot(v[0], v[1]);
    return [v[0] / len, v[1] / len];
}

const arcCenter = fc.tuple(
    fc.double({ min: -10, max: 10, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: -10, max: 10, noNaN: true, noDefaultInfinity: true }));
const arcRadius = fc.double({ min: 20, max: 200, noNaN: true,
    noDefaultInfinity: true });
const startAngle = fc.double({ min: -Math.PI, max: Math.PI, noNaN: true,
    noDefaultInfinity: true });

const arcOf = (center: readonly number[], radius: number, a0: number,
    a1: number): Arc2 =>
    Arc2.fromCenterRadiusEnds(Vector.fromArray([center[0], center[1]]),
        radius,
        Vector.fromArray(circlePoint(center, radius, a0)),
        Vector.fromArray(circlePoint(center, radius, a1)));

describe('SampleCircularArc verification', () => {
    const sampler = new SampleCircularArc();

    it('matches the closed-form rational-quadratic angles for angle <= pi/2',
        () => {
            check(fc.tuple(arcCenter, arcRadius, startAngle,
                fc.double({ min: 0.1, max: 0.5 * Math.PI - 1e-3, noNaN: true,
                    noDefaultInfinity: true })),
            ([center, radius, a0, span]) => {
                const arc = arcOf(center, radius, a0, a0 + span);
                const angle = subtendedAngle(arc);
                const { P0, P2 } = endDirections(arc);
                const points = sampler.compute(arc);
                const n = Math.trunc(radius * angle);
                expect(points.length).toBe(n);
                const expected = acuteSamples(center, radius, P0, P2, n);
                for (let i = 0; i < n; ++i) {
                    expectClose(points[i].values[0], expected[i][0],
                        1e-9 * radius, 1e-11);
                    expectClose(points[i].values[1], expected[i][1],
                        1e-9 * radius, 1e-11);
                }
            });
        });

    it('splits an arc in (pi/2, pi] at its bisector and samples both halves',
        () => {
            check(fc.tuple(arcCenter, arcRadius, startAngle,
                fc.double({ min: 0.5 * Math.PI + 1e-2, max: Math.PI - 1e-2,
                    noNaN: true, noDefaultInfinity: true })),
            ([center, radius, a0, span]) => {
                const arc = arcOf(center, radius, a0, a0 + span);
                const angle = subtendedAngle(arc);
                const { P0, P2 } = endDirections(arc);
                const points = sampler.compute(arc);
                const n = Math.trunc((radius * angle) / 2);
                expect(points.length).toBe(2 * n);
                const bisector = normalized([P0[0] + P2[0], P0[1] + P2[1]]);
                const expected = [
                    ...acuteSamples(center, radius, P0, bisector, n),
                    ...acuteSamples(center, radius, bisector, P2, n)];
                for (let i = 0; i < 2 * n; ++i) {
                    expectClose(points[i].values[0], expected[i][0],
                        1e-8 * radius, 1e-10);
                    expectClose(points[i].values[1], expected[i][1],
                        1e-8 * radius, 1e-10);
                }
            });
        });

    it('places every sample exactly on the circle for any arc', () => {
        // True even for the arcs longer than pi whose subdivision directions
        // upstream computes incorrectly (issue #183): the split directions are
        // still unit vectors, so every subarc is still a circular arc.
        check(fc.tuple(arcCenter, arcRadius, startAngle,
            fc.double({ min: 0.1, max: 2 * Math.PI - 0.1, noNaN: true,
                noDefaultInfinity: true })),
        ([center, radius, a0, span]) => {
            const points = sampler.compute(arcOf(center, radius, a0,
                a0 + span));
            for (const p of points) {
                expectClose(Math.hypot(p.values[0] - center[0],
                    p.values[1] - center[1]), radius, 1e-8 * radius, 1e-10);
            }
        });
    });

    it('sweeps monotonically from the first endpoint for angle <= pi', () => {
        check(fc.tuple(arcCenter, arcRadius, startAngle,
            fc.double({ min: 0.1, max: Math.PI - 1e-3, noNaN: true,
                noDefaultInfinity: true })),
        ([center, radius, a0, span]) => {
            const points = sampler.compute(arcOf(center, radius, a0,
                a0 + span));
            if (points.length === 0) { return; }
            const swept = points.map(p => {
                let a = Math.atan2(p.values[1] - center[1],
                    p.values[0] - center[0]) - a0;
                while (a < -1e-9) { a += 2 * Math.PI; }
                while (a > 2 * Math.PI - 1e-9) { a -= 2 * Math.PI; }
                return a;
            });
            expectClose(swept[0], 0, 1e-9, 1e-9);
            for (let i = 1; i < swept.length; ++i) {
                expect(swept[i]).toBeGreaterThan(swept[i - 1]);
            }
            expect(swept[swept.length - 1]).toBeLessThan(span + 1e-9);
        });
    });

    it('is equivariant under rigid motions of the arc', () => {
        // The sample count depends only on radius * angle, so a rotation and a
        // translation move the samples rigidly. The two arcs recover slightly
        // different subtended angles from their endpoints, so the tolerance
        // carries that difference and the counts may differ by one subarc.
        check(fc.tuple(arcCenter, arcRadius, startAngle,
            fc.double({ min: 0.1, max: 2 * Math.PI - 0.1, noNaN: true,
                noDefaultInfinity: true }),
            fc.double({ min: -Math.PI, max: Math.PI, noNaN: true,
                noDefaultInfinity: true })),
        ([center, radius, a0, span, delta]) => {
            const arcA = arcOf([0, 0], radius, a0, a0 + span);
            const arcB = arcOf(center, radius, a0 + delta,
                a0 + delta + span);
            const a = sampler.compute(arcA);
            const b = sampler.compute(arcB);
            const drift = Math.abs(subtendedAngle(arcA) -
                subtendedAngle(arcB));
            if (a.length !== b.length) {
                // Only a truncation boundary can change the count.
                expect(Math.abs(a.length - b.length)).toBeLessThanOrEqual(4);
                return;
            }
            const c = Math.cos(delta), s = Math.sin(delta);
            const tol = radius * (1e-9 + 4 * drift);
            for (let i = 0; i < a.length; ++i) {
                expectClose(b[i].values[0],
                    c * a[i].values[0] - s * a[i].values[1] + center[0],
                    tol, 1e-10);
                expectClose(b[i].values[1],
                    s * a[i].values[0] + c * a[i].values[1] + center[1],
                    tol, 1e-10);
            }
        });
    });

    it('uses the documented sample counts', () => {
        check(fc.tuple(arcRadius, fc.double({ min: 0.05,
            max: 2 * Math.PI - 0.05, noNaN: true, noDefaultInfinity: true })),
        ([radius, span]) => {
            const arc = arcOf([0, 0], radius, 0, span);
            const points = sampler.compute(arc);
            const angle = subtendedAngle(arc);
            let expected: number;
            if (span <= 0.5 * Math.PI) {
                expected = Math.trunc(radius * angle);
            }
            else if (span <= Math.PI) {
                expected = 2 * Math.trunc((radius * angle) / 2);
            }
            else if (span <= 1.5 * Math.PI) {
                expected = 3 * Math.trunc((radius * (2 * Math.PI - angle)) / 3);
            }
            else {
                expected = 4 * Math.trunc((radius * (2 * Math.PI - angle)) / 4);
            }
            expect(points.length).toBe(expected);
            for (const p of points) {
                expect(Number.isFinite(p.values[0])).toBe(true);
                expect(Number.isFinite(p.values[1])).toBe(true);
            }
        });
    });
});
