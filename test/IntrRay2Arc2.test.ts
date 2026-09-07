import { describe, it, expect } from 'vitest';
import { Arc2 } from '../src/Arc2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { IntrLine2Circle2FI } from '../src/IntrLine2Circle2.js';
import { Line } from '../src/Line.js';
import {
    IntrRay2Arc2TI,
    IntrRay2Arc2FI,
    defaultIntrRay2Arc2FIResult
} from '../src/IntrRay2Arc2.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { dotPerp } from '../src/Vector2.js';
import {
    check, expectClose, expectVectorClose, fc, positive, ray as arbRay,
    rotationFrame, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

// The arc of the unit circle in the first quadrant, from (1,0) to (0,1).
function firstQuadrantArc(): Arc2 {
    return Arc2.fromCenterRadiusEnds(vec(0, 0), 1, vec(1, 0), vec(0, 1));
}

const ti = new IntrRay2Arc2TI();
const fi = new IntrRay2Arc2FI();

describe('IntrRay2Arc2', () => {
    it('has an empty default result', () => {
        const result = defaultIntrRay2Arc2FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('finds the single crossing of an arc at a known point', () => {
        // A ray from the origin along (1,1) meets the arc at
        // (sqrt(2)/2, sqrt(2)/2) at parameter 1.
        const result = fi.find(ray([0, 0], [1, 1]), firstQuadrantArc());
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(Math.SQRT1_2, 12);
        expect(result.point[0].values[1]).toBeCloseTo(Math.SQRT1_2, 12);
        expect(ti.test(ray([0, 0], [1, 1]), firstQuadrantArc()).intersect)
            .toBe(true);
    });

    it('discards circle crossings that are off the arc', () => {
        // A ray along -x from (2,0) crosses the circle at (1,0) and (-1,0).
        // Only (1,0) is on the first-quadrant arc.
        const R = ray([2, 0], [-1, 0]);
        const A = firstQuadrantArc();
        const result = fi.find(R, A);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('finds two crossings when both are on the arc', () => {
        // The arc from (0,-1) to (0,1) through (1,0) is the right half of the
        // circle. A ray along +x from (0.5,-2) does not hit it; a ray along
        // +y from (0.5,-2) hits it twice.
        const A = Arc2.fromCenterRadiusEnds(vec(0, 0), 1, vec(0, -1),
            vec(0, 1));
        const result = fi.find(ray([0.5, -2], [0, 1]), A);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        const y = Math.sqrt(1 - 0.25);
        expect(result.point[0].values[1]).toBeCloseTo(-y, 12);
        expect(result.point[1].values[1]).toBeCloseTo(y, 12);
    });

    it('reports no intersection for a ray pointing away from the arc', () => {
        const R = ray([2, 0], [1, 0]);
        const A = firstQuadrantArc();
        expect(fi.find(R, A).intersect).toBe(false);
        expect(ti.test(R, A).intersect).toBe(false);
        // A ray that misses the circle entirely.
        const M = ray([0, 5], [1, 0]);
        expect(fi.find(M, A).intersect).toBe(false);
        expect(ti.test(M, A).intersect).toBe(false);
    });

    it('agrees with the circle query filtered by the arc on random inputs', () => {
        let state = 8080808;
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
            const R = ray([rand() * 4, rand() * 4], [rand(), rand() + 0.001]);

            const result = fi.find(R, A);
            expect(ti.test(R, A).intersect).toBe(result.intersect);

            // The reference computes the ray-vs-circular-curve intersections
            // directly: the line-circle crossings with nonnegative parameter.
            const circleResult = lcQuery.find(
                Line.fromOriginDirection(R.origin, R.direction),
                Hypersphere.fromCenterRadius(center, radius));
            let expected = 0;
            for (let i = 0; i < circleResult.numIntersections; ++i) {
                if (circleResult.parameter[i] >= 0
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
                // The point is on the circle, on the ray and on the arc.
                expect(length(sub(result.point[i], center)))
                    .toBeCloseTo(radius, 8);
                expect(result.parameter[i]).toBeGreaterThanOrEqual(-1e-12);
                const onRay = add(R.origin,
                    mul(result.parameter[i], R.direction));
                expect(length(sub(result.point[i], onRay)))
                    .toBeCloseTo(0, 10);
                expect(A.containsOnCircle(result.point[i])).toBe(true);
            }
        }
        expect(numOne).toBeGreaterThan(20);
        expect(numTwo).toBeGreaterThan(1);
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrRay2Arc2.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrRay2Arc2 verification', () => {
    const tiQ = new IntrRay2Arc2TI();
    const fiQ = new IntrRay2Arc2FI();

    // An arc on a circle with a positive angular span in (0, 2*pi).
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
        check(fc.tuple(arbRay(2), arbArc), ([R, A]) => {
            expect(tiQ.test(R, A).intersect).toBe(fiQ.find(R, A).intersect);
        });
    });

    it('the reported points are on the ray, on the circle and on the arc', () => {
        check(fc.tuple(arbRay(2), arbArc), ([R, A]) => {
            const r = fiQ.find(R, A);
            expect(r.numIntersections).toBe(r.intersect
                ? r.numIntersections : 0);
            const scale = 1 + A.radius + length(A.center) + length(R.origin);
            for (let i = 0; i < r.numIntersections; ++i) {
                const t = r.parameter[i];
                expect(t).toBeGreaterThanOrEqual(0);
                const X = r.point[i];
                expectVectorClose(X, add(R.origin, mul(t, R.direction)),
                    1e-12 * scale, 1e-12);
                // On the circle carrying the arc.
                expectClose(length(sub(X, A.center)), A.radius,
                    1e-9 * scale, 1e-9);
                // On the arc.
                expect(A.containsOnCircle(X)).toBe(true);
            }
        });
    });

    it('is the line-circle query filtered to t >= 0 and to the arc', () => {
        const lcQuery = new IntrLine2Circle2FI();
        check(fc.tuple(arbRay(2), arbArc), ([R, A]) => {
            const circle = Hypersphere.fromCenterRadius(A.center, A.radius);
            const line = Line.fromOriginDirection(R.origin, R.direction);
            const lc = lcQuery.find(line, circle);
            const kept: number[] = [];
            for (let i = 0; i < lc.numIntersections; ++i) {
                if (lc.parameter[i] >= 0 && A.containsOnCircle(lc.point[i])) {
                    kept.push(lc.parameter[i]);
                }
            }
            const r = fiQ.find(R, A);
            expect(r.numIntersections).toBe(kept.length);
            expect(r.intersect).toBe(kept.length > 0);
            for (let i = 0; i < kept.length; ++i) {
                expect(r.parameter[i]).toBe(kept[i]);
            }
        });
    });

    it('never reports an off-arc hit when the ray origin is inside the disk', () => {
        // The upstream defect (#304): reusing the solid-disk ray-circle query
        // hands the ray ORIGIN to Arc2::Contains when the origin is inside
        // the disk, so an intersection is reported at a point that is not on
        // the circle at all. The port intersects the circular curve instead.
        check(fc.tuple(arbArc, unitVector(2), fc.double({ min: 0, max: 0.95,
            noNaN: true }), fc.double({ min: -Math.PI, max: Math.PI,
            noNaN: true })),
            ([A, d, frac, phi]) => {
                const origin = add(A.center, Vector.fromArray([
                    frac * A.radius * Math.cos(phi),
                    frac * A.radius * Math.sin(phi)]));
                const R = Ray.fromOriginDirection(origin, d);
                const r = fiQ.find(R, A);
                const scale = 1 + A.radius + length(A.center);
                for (let i = 0; i < r.numIntersections; ++i) {
                    expectClose(length(sub(r.point[i], A.center)), A.radius,
                        1e-9 * scale, 1e-9);
                    expect(A.containsOnCircle(r.point[i])).toBe(true);
                }
                // A ray whose origin is strictly inside the disk leaves the
                // circle exactly once.
                expect(r.numIntersections).toBeLessThanOrEqual(1);
            });
    });

    it('pins the origin-inside-the-disk regression', () => {
        // The unit circle centred at the origin; the arc is the upper half
        // (from (1,0) counterclockwise to (-1,0)). A ray from (0,-0.5) going
        // in -y direction leaves the circle at (0,-1), which is NOT on the
        // arc, so there is no intersection. Upstream's solid-disk clipping
        // reported the ray origin (0,-0.5) as an "intersection" and
        // Arc2::Contains accepted it.
        const A = Arc2.fromCenterRadiusEnds(Vector.zero(2), 1,
            Vector.fromArray([1, 0]), Vector.fromArray([-1, 0]));
        const R = Ray.fromOriginDirection(Vector.fromArray([0, -0.5]),
            Vector.fromArray([0, -1]));
        const r = fiQ.find(R, A);
        expect(r.intersect).toBe(false);
        expect(r.numIntersections).toBe(0);
        expect(tiQ.test(R, A).intersect).toBe(false);

        // Reversing the ray direction exits through the arc at (0,1).
        const R2 = Ray.fromOriginDirection(Vector.fromArray([0, -0.5]),
            Vector.fromArray([0, 1]));
        const r2 = fiQ.find(R2, A);
        expect(r2.intersect).toBe(true);
        expect(r2.numIntersections).toBe(1);
        expectVectorClose(r2.point[0], Vector.fromArray([0, 1]), 1e-12, 1e-12);
        expect(r2.parameter[0]).toBeCloseTo(1.5, 12);
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(arbRay(2), arbArc, rotationFrame(2),
            wellScaledVector(2)),
            ([R, A, Rot, t]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(Rot[0], v), dot(Rot[1], v)]);
                const map = (v: Vector) => add(rot(v), t);
                const R2 = Ray.fromOriginDirection(map(R.origin),
                    rot(R.direction));
                const A2 = Arc2.fromCenterRadiusEnds(map(A.center), A.radius,
                    map(A.end[0]), map(A.end[1]));
                // Skip the configurations that are within rounding of a
                // tangency, of the ray origin or of an arc endpoint, where
                // the strict tests can flip under a rigid motion. The guard
                // looks at the unfiltered line-circle roots, not only at the
                // ones the query kept.
                const lc = new IntrLine2Circle2FI().find(
                    Line.fromOriginDirection(R.origin, R.direction),
                    Hypersphere.fromCenterRadius(A.center, A.radius));
                if (lc.numIntersections === 1) { return; }  // tangency
                const sc = 1 + A.radius + length(A.center) + length(R.origin);
                for (let i = 0; i < lc.numIntersections; ++i) {
                    if (Math.abs(lc.parameter[i]) < 1e-6 * sc) { return; }
                    const dp0 = dotPerp(sub(lc.point[i], A.end[0]),
                        sub(A.end[1], A.end[0]));
                    if (Math.abs(dp0) < 1e-6 * A.radius * A.radius) {
                        return;
                    }
                }
                const r0 = fiQ.find(R, A);
                const r1 = fiQ.find(R2, A2);
                expect(r1.numIntersections).toBe(r0.numIntersections);
                for (let i = 0; i < r0.numIntersections; ++i) {
                    expectClose(r1.parameter[i], r0.parameter[i], 1e-7, 1e-7);
                }
            });
    });

    it('reports a ray that misses the circle as empty', () => {
        check(fc.tuple(arbArc, unitVector(2), positive(5, 1)),
            ([A, d, extra]) => {
                // Start far away and aim away from the arc.
                const origin = add(A.center,
                    mul(A.radius + extra + 1, d));
                const R = Ray.fromOriginDirection(origin, d);
                const r = fiQ.find(R, A);
                expect(r.intersect).toBe(false);
                expect(r.numIntersections).toBe(0);
                expect(r.point[0].values).toEqual([0, 0]);
                expect(r.point[1].values).toEqual([0, 0]);
            });
    });

    it('reports the default result without aliasing the query internals', () => {
        const d0 = defaultIntrRay2Arc2FIResult();
        const d1 = defaultIntrRay2Arc2FIResult();
        expect(d0.point[0]).not.toBe(d1.point[0]);
        const A = Arc2.fromCenterRadiusEnds(Vector.zero(2), 1,
            Vector.fromArray([1, 0]), Vector.fromArray([-1, 0]));
        const R = Ray.fromOriginDirection(Vector.fromArray([-5, 0.5]),
            Vector.fromArray([1, 0]));
        const r0 = fiQ.find(R, A);
        const saved = r0.point[0].clone();
        fiQ.find(R, A);
        expect(r0.point[0].values).toEqual(saved.values);
    });
});
