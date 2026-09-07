import { describe, it, expect } from 'vitest';
import { Arc2 } from '../src/Arc2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { IntrLine2Circle2FI } from '../src/IntrLine2Circle2.js';
import {
    IntrSegment2Arc2TI,
    IntrSegment2Arc2FI,
    defaultIntrSegment2Arc2FIResult
} from '../src/IntrSegment2Arc2.js';
import { Line } from '../src/Line.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { dotPerp } from '../src/Vector2.js';
import {
    check, expectClose, expectVectorClose, fc, positive,
    rotationFrame, segment as arbSegment, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

// The arc of the unit circle in the first quadrant, from (1,0) to (0,1).
function firstQuadrantArc(): Arc2 {
    return Arc2.fromCenterRadiusEnds(vec(0, 0), 1, vec(1, 0), vec(0, 1));
}

const ti = new IntrSegment2Arc2TI();
const fi = new IntrSegment2Arc2FI();

describe('IntrSegment2Arc2', () => {
    it('has an empty default result', () => {
        const result = defaultIntrSegment2Arc2FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
        expect(result.point[0].values).toEqual([0, 0]);
        expect(result.point[1].values).toEqual([0, 0]);
    });

    it('finds a crossing at a known point', () => {
        // The segment from (0,0) to (2,2) crosses the arc at
        // (sqrt(2)/2, sqrt(2)/2). The centered form has origin (1,1), unit
        // direction (1,1)/sqrt(2) and extent sqrt(2), so the parameter is
        // 1 - sqrt(2).
        const result = fi.find(seg([0, 0], [2, 2]), firstQuadrantArc());
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(Math.SQRT1_2, 12);
        expect(result.point[0].values[1]).toBeCloseTo(Math.SQRT1_2, 12);
        expect(result.parameter[0]).toBeCloseTo(1 - Math.SQRT2, 12);
        expect(ti.test(seg([0, 0], [2, 2]), firstQuadrantArc()).intersect)
            .toBe(true);
    });

    it('discards crossings outside the segment', () => {
        // The full line through (2,0) and (-2,0) crosses the circle at
        // (1,0) and (-1,0), but the segment from (2,0) to (1.5,0) reaches
        // neither.
        const S = seg([2, 0], [1.5, 0]);
        const A = firstQuadrantArc();
        expect(fi.find(S, A).intersect).toBe(false);
        expect(ti.test(S, A).intersect).toBe(false);
        // Extending the segment past (1,0) produces one crossing.
        expect(fi.find(seg([2, 0], [0.5, 0]), A).numIntersections).toBe(1);
    });

    it('discards crossings that are off the arc', () => {
        // The segment from (-2,0) to (2,0) crosses the circle at (-1,0) and
        // (1,0); only (1,0) is on the first-quadrant arc.
        const result = fi.find(seg([-2, 0], [2, 0]), firstQuadrantArc());
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('finds two crossings when both are on the arc and the segment', () => {
        const A = Arc2.fromCenterRadiusEnds(vec(0, 0), 1, vec(0, -1),
            vec(0, 1));
        const result = fi.find(seg([0.5, -2], [0.5, 2]), A);
        expect(result.numIntersections).toBe(2);
        const y = Math.sqrt(1 - 0.25);
        expect(result.point[0].values[1]).toBeCloseTo(-y, 12);
        expect(result.point[1].values[1]).toBeCloseTo(y, 12);
    });

    it('agrees with the line-circle curve query on random inputs', () => {
        let state = 5150;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        const lcQuery = new IntrLine2Circle2FI();
        let numOne = 0;
        let numTwo = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const center = vec(rand(), rand());
            const radius = 0.5 + Math.abs(rand()) * 2;
            const a0 = rand() * Math.PI;
            const a1 = a0 + 0.3 + Math.abs(rand()) * 2;
            const A = Arc2.fromCenterRadiusEnds(center, radius,
                add(center, vec(radius * Math.cos(a0), radius * Math.sin(a0))),
                add(center, vec(radius * Math.cos(a1), radius * Math.sin(a1))));
            const S = seg([rand() * 4, rand() * 4], [rand() * 4, rand() * 4]);

            const result = fi.find(S, A);
            expect(ti.test(S, A).intersect).toBe(result.intersect);

            const { center: segOrigin, direction: segDirection,
                extent: segExtent } = S.getCenteredForm();
            const circleResult = lcQuery.find(
                Line.fromOriginDirection(segOrigin, segDirection),
                Hypersphere.fromCenterRadius(center, radius));
            let expected = 0;
            for (let i = 0; i < circleResult.numIntersections; ++i) {
                if (Math.abs(circleResult.parameter[i]) <= segExtent
                    && A.containsOnCircle(circleResult.point[i])) {
                    ++expected;
                }
            }
            expect(result.numIntersections).toBe(expected);
            if (expected === 1) {
                ++numOne;
            }
            else if (expected === 2) {
                ++numTwo;
            }

            for (let i = 0; i < result.numIntersections; ++i) {
                // The point is on the circle, on the segment and on the arc.
                expect(length(sub(result.point[i], center)))
                    .toBeCloseTo(radius, 8);
                expect(Math.abs(result.parameter[i]))
                    .toBeLessThanOrEqual(segExtent + 1e-12);
                const onSegment = add(segOrigin,
                    mul(result.parameter[i], segDirection));
                expect(length(sub(result.point[i], onSegment)))
                    .toBeCloseTo(0, 10);
                expect(A.containsOnCircle(result.point[i])).toBe(true);
            }
        }
        expect(numOne).toBeGreaterThan(20);
        expect(numTwo).toBeGreaterThanOrEqual(1);
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrSegment2Arc2.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrSegment2Arc2 verification', () => {
    const tiQ = new IntrSegment2Arc2TI();
    const fiQ = new IntrSegment2Arc2FI();

    const arbArc = fc.tuple(wellScaledVector(2), positive(4),
        fc.double({ min: -Math.PI, max: Math.PI, noNaN: true }),
        fc.double({ min: 0.05, max: 2 * Math.PI - 0.05, noNaN: true }))
        .map(([c, r, a0, span]) => {
            const e0 = Vector.fromArray([c.values[0] + r * Math.cos(a0),
                c.values[1] + r * Math.sin(a0)]);
            const a1 = a0 + span;
            const e1 = Vector.fromArray([c.values[0] + r * Math.cos(a1),
                c.values[1] + r * Math.sin(a1)]);
            return Arc2.fromCenterRadiusEnds(c, r, e0, e1);
        });

    it('TI and FI agree on intersect', () => {
        check(fc.tuple(arbSegment(2), arbArc), ([S, A]) => {
            expect(tiQ.test(S, A).intersect).toBe(fiQ.find(S, A).intersect);
        });
    });

    it('the reported points use the centered segment form and are on the arc', () => {
        check(fc.tuple(arbSegment(2), arbArc), ([S, A]) => {
            const r = fiQ.find(S, A);
            const cf = S.getCenteredForm();
            const scale = 1 + A.radius + length(A.center) + length(cf.center);
            for (let i = 0; i < r.numIntersections; ++i) {
                const t = r.parameter[i];
                // The documented convention: P = C + t*D with |t| <= extent.
                expect(Math.abs(t)).toBeLessThanOrEqual(cf.extent);
                const X = r.point[i];
                expectVectorClose(X, add(cf.center, mul(t, cf.direction)),
                    1e-12 * scale, 1e-12);
                expectClose(length(sub(X, A.center)), A.radius,
                    1e-9 * scale, 1e-9);
                expect(A.containsOnCircle(X)).toBe(true);
                // The point is also between the segment endpoints.
                const d0 = length(sub(X, S.p[0]));
                const d1 = length(sub(X, S.p[1]));
                expectClose(d0 + d1, length(sub(S.p[1], S.p[0])),
                    1e-9 * scale, 1e-9);
            }
        });
    });

    it('is the line-circle query filtered to |t| <= extent and to the arc', () => {
        const lcQuery = new IntrLine2Circle2FI();
        check(fc.tuple(arbSegment(2), arbArc), ([S, A]) => {
            const cf = S.getCenteredForm();
            const lc = lcQuery.find(
                Line.fromOriginDirection(cf.center, cf.direction),
                Hypersphere.fromCenterRadius(A.center, A.radius));
            const kept: number[] = [];
            for (let i = 0; i < lc.numIntersections; ++i) {
                if (Math.abs(lc.parameter[i]) <= cf.extent
                    && A.containsOnCircle(lc.point[i])) {
                    kept.push(lc.parameter[i]);
                }
            }
            const r = fiQ.find(S, A);
            expect(r.numIntersections).toBe(kept.length);
            expect(r.intersect).toBe(kept.length > 0);
            for (let i = 0; i < kept.length; ++i) {
                expect(r.parameter[i]).toBe(kept[i]);
            }
        });
    });

    it('never reports an off-arc hit when an endpoint is inside the disk', () => {
        // The upstream defect (#304): the solid-disk segment-circle query
        // clips the t-interval to [-e,e], so a segment endpoint inside the
        // disk is reported as an "intersection point" and handed to
        // Arc2::Contains, which assumes its argument is on the circle.
        check(fc.tuple(arbArc, unitVector(2),
            fc.double({ min: 0, max: 0.9, noNaN: true }),
            fc.double({ min: -Math.PI, max: Math.PI, noNaN: true }),
            positive(6, 0.1)),
            ([A, d, frac, phi, len]) => {
                const inside = add(A.center, Vector.fromArray([
                    frac * A.radius * Math.cos(phi),
                    frac * A.radius * Math.sin(phi)]));
                const S = Segment.fromEndpoints(inside,
                    add(inside, mul(len, d)));
                const r = fiQ.find(S, A);
                const scale = 1 + A.radius + length(A.center);
                for (let i = 0; i < r.numIntersections; ++i) {
                    expectClose(length(sub(r.point[i], A.center)), A.radius,
                        1e-9 * scale, 1e-9);
                    expect(A.containsOnCircle(r.point[i])).toBe(true);
                }
                expect(r.numIntersections).toBeLessThanOrEqual(1);
            });
    });

    it('pins the endpoint-inside-the-disk regression', () => {
        // The unit circle centred at the origin; the arc is the upper half.
        // The segment from (0,-0.5) to (0,-0.9) is entirely inside the disk
        // and meets neither the circle nor the arc. Upstream's solid-disk
        // clipping reported the endpoint (0,-0.5) and Arc2::Contains accepted
        // it (dotPerp((0,-0.5)-(1,0), (-2,0)) = -1*-0.5... > 0).
        const A = Arc2.fromCenterRadiusEnds(Vector.zero(2), 1,
            Vector.fromArray([1, 0]), Vector.fromArray([-1, 0]));
        const S = Segment.fromEndpoints(Vector.fromArray([0, -0.5]),
            Vector.fromArray([0, -0.9]));
        const r = fiQ.find(S, A);
        expect(r.intersect).toBe(false);
        expect(r.numIntersections).toBe(0);
        expect(tiQ.test(S, A).intersect).toBe(false);

        // A segment reaching up through the arc does intersect it.
        const S2 = Segment.fromEndpoints(Vector.fromArray([0, -0.5]),
            Vector.fromArray([0, 2]));
        const r2 = fiQ.find(S2, A);
        expect(r2.intersect).toBe(true);
        expect(r2.numIntersections).toBe(1);
        expectVectorClose(r2.point[0], Vector.fromArray([0, 1]), 1e-12, 1e-12);
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(arbSegment(2), arbArc, rotationFrame(2),
            wellScaledVector(2)),
            ([S, A, Rot, t]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(Rot[0], v), dot(Rot[1], v)]);
                const map = (v: Vector) => add(rot(v), t);
                const cf = S.getCenteredForm();
                const lc = new IntrLine2Circle2FI().find(
                    Line.fromOriginDirection(cf.center, cf.direction),
                    Hypersphere.fromCenterRadius(A.center, A.radius));
                if (lc.numIntersections === 1) { return; }  // tangency
                const sc = 1 + A.radius + length(A.center) + length(cf.center);
                for (let i = 0; i < lc.numIntersections; ++i) {
                    if (Math.abs(Math.abs(lc.parameter[i]) - cf.extent)
                        < 1e-6 * sc) {
                        return;
                    }
                    const dp0 = dotPerp(sub(lc.point[i], A.end[0]),
                        sub(A.end[1], A.end[0]));
                    if (Math.abs(dp0) < 1e-6 * A.radius * A.radius) {
                        return;
                    }
                }
                const S2 = Segment.fromEndpoints(map(S.p[0]), map(S.p[1]));
                const A2 = Arc2.fromCenterRadiusEnds(map(A.center), A.radius,
                    map(A.end[0]), map(A.end[1]));
                expect(fiQ.find(S2, A2).numIntersections)
                    .toBe(fiQ.find(S, A).numIntersections);
            });
    });

    it('reports a degenerate (point) segment outside the circle as empty', () => {
        check(fc.tuple(arbArc, wellScaledVector(2)), ([A, p]) => {
            if (Math.abs(length(sub(p, A.center)) - A.radius) < 1e-6) {
                return;
            }
            const S = Segment.fromEndpoints(p, p.clone());
            const r = fiQ.find(S, A);
            expect(r.intersect).toBe(false);
            expect(r.numIntersections).toBe(0);
        });
    });
});
