import { describe, expect, it } from 'vitest';
import { Arc2 } from '../src/Arc2.js';
import { DistLine2Arc2 } from '../src/DistLine2Arc2.js';
import { Line } from '../src/Line.js';
import { DistLine2Circle2 } from '../src/DistLine2Circle2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, positive,
    rotationFrame, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function arc(center: number[], radius: number, a0: number, a1: number): Arc2 {
    const c = v(...center);
    return Arc2.fromCenterRadiusEnds(c, radius,
        v(c.values[0] + radius * Math.cos(a0),
            c.values[1] + radius * Math.sin(a0)),
        v(c.values[0] + radius * Math.cos(a1),
            c.values[1] + radius * Math.sin(a1)));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistLine2Arc2', () => {
    const query = new DistLine2Arc2();

    it('uses a circle closest point that lies on the arc', () => {
        // The upper half of the unit circle and the line y = 3.
        const a = arc([0, 0], 1, 0, Math.PI);
        const result = query.compute(line([0, 3], [1, 0]), a);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.closest[0][1].values[0]).toBeCloseTo(0, 10);
        expect(result.closest[0][1].values[1]).toBeCloseTo(1, 10);
    });

    it('reports two closest pairs when both circle points are on the arc',
        () => {
            // The full-ish arc from angle -3 to 3 (going counterclockwise
            // through pi) and the vertical line x = 3. The line-circle query
            // returns one pair; use instead a line through the circle so that
            // the two intersection points are both on the arc.
            const a = arc([0, 0], 1, -Math.PI / 2, Math.PI / 2);
            // The vertical line x = 0.5 crosses the right half circle twice.
            const result = query.compute(line([0.5, 0], [0, 1]), a);
            expect(result.numClosestPairs).toBe(2);
            expect(result.distance).toBeCloseTo(0, 10);
            for (let j = 0; j < 2; ++j) {
                const p = result.closest[j][1];
                expect(Math.sqrt(dot(p, p))).toBeCloseTo(1, 10);
                expect(p.values[0]).toBeCloseTo(0.5, 10);
            }
        });

    it('falls back to arc endpoints when no circle point is on the arc',
        () => {
            // The upper half arc, and the line y = -3 below it.
            const a = arc([0, 0], 1, 0, Math.PI);
            const result = query.compute(line([0, -3], [1, 0]), a);
            // Both endpoints (1,0) and (-1,0) are 3 units from the line, so
            // two closest pairs are reported.
            expect(result.numClosestPairs).toBe(2);
            expect(result.distance).toBeCloseTo(3, 10);
            expect(result.closest[0][1].values[0]).toBeCloseTo(1, 10);
            expect(result.closest[1][1].values[0]).toBeCloseTo(-1, 10);
        });

    it('selects the single nearer arc endpoint', () => {
        const a = arc([0, 0], 1, 0, Math.PI / 2);
        // A line below and to the right; the endpoint (1,0) is nearer.
        const result = query.compute(line([0, -3], [1, 1]), a);
        expect(result.numClosestPairs).toBe(1);
        expect(result.closest[0][1].values[0]).toBeCloseTo(1, 10);
        expect(result.closest[0][1].values[1]).toBeCloseTo(0, 10);
    });

    it('reports the line parameter of the closest line point', () => {
        const a = arc([0, 0], 1, 0, Math.PI);
        const ln = line([0, 3], [1, 0]);
        const result = query.compute(ln, a);
        const p = add(ln.origin, mul(result.parameter[0], ln.direction));
        for (let i = 0; i < 2; ++i) {
            expect(p.values[i]).toBeCloseTo(result.closest[0][0].values[i],
                10);
        }
    });

    it('reports zero distance for a line tangent at an arc point', () => {
        const a = arc([0, 0], 1, 0, Math.PI);
        const result = query.compute(line([0, 1], [1, 0]), a);
        expect(result.distance).toBeCloseTo(0, 10);
    });

    it('agrees with a dense sampling of the arc', () => {
        const rnd = makeRandom(90210);
        const a0 = -0.3, a1 = 2.4;
        const center = [0.25, -0.5];
        const radius = 1.5;
        const a = arc(center, radius, a0, a1);

        for (let trial = 0; trial < 50; ++trial) {
            const origin = v(6 * rnd() - 3, 6 * rnd() - 3);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const ln = Line.fromOriginDirection(origin, dir);
            const result = query.compute(ln, a);
            expect(result.numClosestPairs).toBeGreaterThanOrEqual(1);

            for (let j = 0; j < result.numClosestPairs; ++j) {
                // The line point matches the parameter.
                const p = add(ln.origin, mul(result.parameter[j],
                    ln.direction));
                for (let i = 0; i < 2; ++i) {
                    expect(p.values[i]).toBeCloseTo(
                        result.closest[j][0].values[i], 7);
                }
                // The arc point is on the circle of the arc.
                const d = sub(result.closest[j][1], v(...center));
                expect(Math.sqrt(dot(d, d))).toBeCloseTo(radius, 7);
                // The pair realizes the reported distance.
                const e = sub(result.closest[j][0], result.closest[j][1]);
                expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 6);
            }

            // No sampled arc point is closer to the line.
            const n = 4000;
            const dd = dot(ln.direction, ln.direction);
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= n; ++i) {
                const t = a0 + (i / n) * (a1 - a0);
                const q = v(center[0] + radius * Math.cos(t),
                    center[1] + radius * Math.sin(t));
                const w = sub(q, ln.origin);
                const s = dot(w, ln.direction) / dd;
                const f = sub(w, mul(s, ln.direction));
                best = Math.min(best, dot(f, f));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistLine2Arc2.h.
// ---------------------------------------------------------------------------

describe('DistLine2Arc2 verification', () => {
    const query = new DistLine2Arc2();

    // An arc is generated from a center, a radius, a start angle and a sweep
    // in (0, 2*pi). wellScaledVector/positive keep the geometry away from the
    // subnormal range so that the |P-C| = r invariants hold to 1e-12.
    const arcArb = fc.tuple(wellScaledVector(2, -5, 5), positive(4, 0.5),
        finite(-Math.PI, Math.PI), finite(0.05, 2 * Math.PI - 0.05))
        .map(([c, r, a0, sweep]) => {
            const e0 = v(c.values[0] + r * Math.cos(a0),
                c.values[1] + r * Math.sin(a0));
            const e1 = v(c.values[0] + r * Math.cos(a0 + sweep),
                c.values[1] + r * Math.sin(a0 + sweep));
            return { arc: Arc2.fromCenterRadiusEnds(c, r, e0, e1), a0, sweep };
        });

    const lineArb = fc.tuple(wellScaledVector(2, -6, 6), unitVector(2))
        .map(([o, d]) => Line.fromOriginDirection(o, d));

    // Distance from a point to the line P + t*D (D is not required to be unit
    // length), computed independently of the query under test.
    function pointLineDistance(p: Vector, ln: Line): number {
        const diff = sub(p, ln.origin);
        const t = dot(diff, ln.direction) / dot(ln.direction, ln.direction);
        return length(sub(diff, mul(t, ln.direction)));
    }

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(lineArb, arcArb), ([ln, { arc }]) => {
            const r = query.compute(ln, arc);
            expect(r.numClosestPairs === 1 || r.numClosestPairs === 2)
                .toBe(true);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            for (let j = 0; j < r.numClosestPairs; ++j) {
                // closest[j][0] is the line point at parameter[j].
                expectVectorClose(r.closest[j][0],
                    add(ln.origin, mul(r.parameter[j], ln.direction)),
                    1e-9, 1e-9);
                // closest[j][1] is on the circle of the arc and on the arc.
                expectClose(length(sub(r.closest[j][1], arc.center)),
                    arc.radius, 1e-9, 1e-9);
                expect(arc.containsOnCircle(r.closest[j][1])
                    // Endpoints can fall marginally outside the halfplane
                    // test because of rounding; accept them explicitly.
                    || length(sub(r.closest[j][1], arc.end[0])) < 1e-9
                    || length(sub(r.closest[j][1], arc.end[1])) < 1e-9)
                    .toBe(true);
                // The reported distance is the distance of this pair.
                expectClose(length(sub(r.closest[j][0], r.closest[j][1])),
                    r.distance, 1e-9, 1e-9);
            }
        }, 150);
    });

    it('is not larger than the distance to any sampled arc point', () => {
        check(fc.tuple(lineArb, arcArb), ([ln, { arc, a0, sweep }]) => {
            const r = query.compute(ln, arc);
            const n = 4096;
            let best = Number.POSITIVE_INFINITY;
            for (let i = 0; i <= n; ++i) {
                const a = a0 + (i / n) * sweep;
                const p = v(arc.center.values[0] + arc.radius * Math.cos(a),
                    arc.center.values[1] + arc.radius * Math.sin(a));
                best = Math.min(best, pointLineDistance(p, ln));
            }
            expect(r.distance).toBeLessThanOrEqual(best + 1e-9);
            // The distance to a line is 1-Lipschitz in arc length, so the
            // sampled minimum overestimates the true minimum by at most half
            // an arc-length step; a wrong branch shows up as a larger gap.
            const step = (arc.radius * sweep) / n;
            expect(r.distance).toBeGreaterThanOrEqual(best - step);
        }, 100);
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(lineArb, arcArb, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([ln, { arc }, frame, shift]) => {
            const xf = (p: Vector): Vector => add(shift,
                v(frame[0].values[0] * p.values[0]
                    + frame[1].values[0] * p.values[1],
                    frame[0].values[1] * p.values[0]
                    + frame[1].values[1] * p.values[1]));
            const xfDir = (d: Vector): Vector =>
                v(frame[0].values[0] * d.values[0]
                    + frame[1].values[0] * d.values[1],
                    frame[0].values[1] * d.values[0]
                    + frame[1].values[1] * d.values[1]);
            const movedLine = Line.fromOriginDirection(xf(ln.origin),
                xfDir(ln.direction));
            const movedArc = Arc2.fromCenterRadiusEnds(xf(arc.center),
                arc.radius, xf(arc.end[0]), xf(arc.end[1]));
            const r0 = query.compute(ln, arc);
            const r1 = query.compute(movedLine, movedArc);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
        }, 100);
    });

    it('matches the line-circle query when the arc is nearly a full circle',
        () => {
            // An arc whose sweep is close to 2*pi contains almost every
            // circle point, so the arc query must agree with the circle
            // query for lines that do not pass near the omitted sliver.
            const c = v(1, -2);
            const r = 3;
            const a0 = 0.3;
            const a1 = a0 + 2 * Math.PI - 1e-6;
            const wide = Arc2.fromCenterRadiusEnds(c, r,
                v(c.values[0] + r * Math.cos(a0),
                    c.values[1] + r * Math.sin(a0)),
                v(c.values[0] + r * Math.cos(a1),
                    c.values[1] + r * Math.sin(a1)));
            const circle = Hypersphere.fromCenterRadius(c, r);
            const lc = new DistLine2Circle2();
            check(lineArb, ln => {
                const ra = query.compute(ln, wide);
                const rc = lc.compute(ln, circle);
                expectClose(ra.distance, rc.distance, 1e-9, 1e-9);
            }, 100);
        });
});
