import { describe, it, expect } from 'vitest';
import { Cone } from '../src/Cone.js';
import {
    IntrLine3Cone3FI,
    IntrLine3Cone3FIResultType,
    defaultIntrLine3Cone3FIResult,
    intrLine3Cone3Convert,
    intrLine3Cone3ConvertPoint
} from '../src/IntrLine3Cone3.js';
import { Line } from '../src/Line.js';
import { QFNumber } from '../src/QFNumber.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, negate, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, line as arbLine, positive,
    rotationFrame, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function line(origin: number[], direction: number[]): Line {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Line.fromOriginDirection(Vector.fromArray(origin), d);
}

// A cone with apex at 'origin', axis 'direction' and the given half-angle.
// A negative maxHeight means the cone is infinite.
function cone(origin: number[], direction: number[], angle: number,
    minHeight: number, maxHeight: number): Cone {
    const C = new Cone(3);
    const d = Vector.fromArray(direction);
    normalize(d);
    C.ray = Ray.fromOriginDirection(Vector.fromArray(origin), d);
    C.setAngle(angle);
    if (maxHeight < 0) {
        if (minHeight > 0) {
            C.makeInfiniteTruncatedCone(minHeight);
        }
        else {
            C.makeInfiniteCone();
        }
    }
    else if (minHeight > 0) {
        C.makeConeFrustum(minHeight, maxHeight);
    }
    else {
        C.makeFiniteCone(maxHeight);
    }
    return C;
}

// True when X is in the solid cone (with a small tolerance).
function inSolidCone(C: Cone, X: Vector, tolerance: number): boolean {
    const diff = sub(X, C.ray.origin);
    const h = dot(C.ray.direction, diff);
    if (h < C.getMinHeight() - tolerance) {
        return false;
    }
    if (C.isFinite() && h > C.getMaxHeight() + tolerance) {
        return false;
    }
    return h >= length(diff) * C.cosAngle - tolerance;
}

const fi = new IntrLine3Cone3FI();

describe('IntrLine3Cone3', () => {
    it('has an empty default result', () => {
        const result = defaultIntrLine3Cone3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.type).toBe(IntrLine3Cone3FIResultType.isEmpty);
        expect(intrLine3Cone3Convert(result.t[0])).toBe(0);
        expect(intrLine3Cone3ConvertPoint(result.P[0]).values)
            .toEqual([0, 0, 0]);
    });

    it('converts quadratic-field values to numbers', () => {
        // 1 + 2*sqrt(9) = 7.
        expect(intrLine3Cone3Convert(new QFNumber(1, 2, 9))).toBeCloseTo(7, 12);
    });

    it('cuts an infinite cone in a segment', () => {
        // The 45-degree cone about +z has radius 1 at height 1, so the line
        // x -> (x,0,1) meets it for x in [-1,1].
        const C = cone([0, 0, 0], [0, 0, 1], Math.PI / 4, 0, -1);
        const result = fi.find(line([0, 0, 1], [1, 0, 0]), C);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(IntrLine3Cone3FIResultType.isSegment);
        expect(intrLine3Cone3Convert(result.t[0])).toBeCloseTo(-1, 12);
        expect(intrLine3Cone3Convert(result.t[1])).toBeCloseTo(1, 12);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
        expect(P0.values[0]).toBeCloseTo(-1, 12);
        expect(P1.values[0]).toBeCloseTo(1, 12);
        expect(P0.values[2]).toBeCloseTo(1, 12);
        expect(P1.values[2]).toBeCloseTo(1, 12);
    });

    it('reports a positive ray for a line along the infinite cone axis', () => {
        const C = cone([0, 0, 0], [0, 0, 1], Math.PI / 4, 0, -1);
        const result = fi.find(line([0, 0, -1], [0, 0, 1]), C);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(IntrLine3Cone3FIResultType.isRayPositive);
        expect(intrLine3Cone3Convert(result.t[0])).toBeCloseTo(1, 12);
        // P[0] is the ray origin (the cone apex) and P[1] is the direction.
        expect(intrLine3Cone3ConvertPoint(result.P[0]).values)
            .toEqual([0, 0, 0]);
        expect(intrLine3Cone3ConvertPoint(result.P[1]).values)
            .toEqual([0, 0, 1]);
    });

    it('reports a negative ray when the line direction opposes the axis', () => {
        const C = cone([0, 0, 0], [0, 0, 1], Math.PI / 4, 0, -1);
        const result = fi.find(line([0, 0, 1], [0, 0, -1]), C);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(IntrLine3Cone3FIResultType.isRayNegative);
        expect(intrLine3Cone3Convert(result.t[1])).toBeCloseTo(1, 12);
        // P[0] is the ray endpoint (the apex) and P[1] is the direction.
        expect(intrLine3Cone3ConvertPoint(result.P[0]).values)
            .toEqual([0, 0, 0]);
        expect(intrLine3Cone3ConvertPoint(result.P[1]).values)
            .toEqual([0, 0, -1]);
    });

    it('misses the positive cone when the line is behind the apex', () => {
        const C = cone([0, 0, 0], [0, 0, 1], Math.PI / 4, 0, -1);
        const result = fi.find(line([0, 0, -1], [1, 0, 0]), C);
        expect(result.intersect).toBe(false);
        expect(result.type).toBe(IntrLine3Cone3FIResultType.isEmpty);
        expect(intrLine3Cone3ConvertPoint(result.P[0]).values)
            .toEqual([0, 0, 0]);
    });

    it('misses a cone the line does not reach', () => {
        // A narrow cone about +z: the line at height 1 with |x| >= 1 misses.
        const C = cone([0, 0, 0], [0, 0, 1], Math.PI / 8, 0, -1);
        const result = fi.find(line([0, 5, 1], [1, 0, 0]), C);
        expect(result.intersect).toBe(false);
        expect(result.type).toBe(IntrLine3Cone3FIResultType.isEmpty);
    });

    it('clamps to a finite cone height range', () => {
        // The cone is truncated at height 2, so a line at height 3 misses.
        const finiteCone = cone([0, 0, 0], [0, 0, 1], Math.PI / 4, 0, 2);
        expect(fi.find(line([0, 0, 3], [1, 0, 0]), finiteCone).type)
            .toBe(IntrLine3Cone3FIResultType.isEmpty);
        // A line at height 1 still cuts a segment.
        const inside = fi.find(line([0, 0, 1], [1, 0, 0]), finiteCone);
        expect(inside.type).toBe(IntrLine3Cone3FIResultType.isSegment);
    });

    it('clamps a cone frustum axis line to the height range', () => {
        const frustum = cone([0, 0, 0], [0, 0, 1], Math.PI / 4, 1, 3);
        const result = fi.find(line([0, 0, 0], [0, 0, 1]), frustum);
        expect(result.type).toBe(IntrLine3Cone3FIResultType.isSegment);
        expect(intrLine3Cone3Convert(result.t[0])).toBeCloseTo(1, 12);
        expect(intrLine3Cone3Convert(result.t[1])).toBeCloseTo(3, 12);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
        expect(P0.values[2]).toBeCloseTo(1, 12);
        expect(P1.values[2]).toBeCloseTo(3, 12);
    });

    it('reports the apex when a line meets only the cone vertex', () => {
        // A line through the apex perpendicular to the axis is outside the
        // double-sided cone except at the vertex.
        const C = cone([0, 0, 0], [0, 0, 1], Math.PI / 4, 0, -1);
        const result = fi.find(line([0, 0, 0], [1, 0, 0]), C);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(IntrLine3Cone3FIResultType.isPoint);
        expect(intrLine3Cone3Convert(result.t[0])).toBe(0);
        expect(intrLine3Cone3ConvertPoint(result.P[0]).values)
            .toEqual([0, 0, 0]);
        expect(intrLine3Cone3ConvertPoint(result.P[1]).values)
            .toEqual([0, 0, 0]);

        // For a frustum whose heights start at 1, the vertex is out of range.
        const frustum = cone([0, 0, 0], [0, 0, 1], Math.PI / 4, 1, 3);
        expect(fi.find(line([0, 0, 0], [1, 0, 0]), frustum).type)
            .toBe(IntrLine3Cone3FIResultType.isEmpty);
    });

    it('clamps a tangent point by its own height, not the vertex height', () => {
        // The line x = 1, z = 1 with direction (0,1,0) is tangent to the
        // 45-degree cone about +z at (1,0,1), whose height is 1. The
        // discriminant is exactly zero only when cos^2(angle) is exactly 1/2,
        // which cos(pi/4)^2 is not in double precision, so the test sets the
        // squared cosine directly. This exercises the port's fix for the
        // upstream vertex test (see IntrLine3Cone3.ts).
        const makeCone = (minHeight: number, maxHeight: number) => {
            const C = cone([0, 0, 0], [0, 0, 1], Math.PI / 4, minHeight,
                maxHeight);
            C.cosAngleSqr = 0.5;
            return C;
        };

        const infinite = makeCone(0, -1);
        const result = fi.find(line([1, 0, 1], [0, 1, 0]), infinite);
        expect(result.type).toBe(IntrLine3Cone3FIResultType.isPoint);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        expect(P0.values[0]).toBeCloseTo(1, 12);
        expect(P0.values[1]).toBeCloseTo(0, 12);
        expect(P0.values[2]).toBeCloseTo(1, 12);

        // The tangent point is at height 1, which is inside the [0.5, 3]
        // range of this frustum. Upstream would clamp with the vertex height
        // 0 and report no intersection.
        const inRange = fi.find(line([1, 0, 1], [0, 1, 0]), makeCone(0.5, 3));
        expect(inRange.type).toBe(IntrLine3Cone3FIResultType.isPoint);

        // The same tangent point is outside the [2, 3] range.
        const outOfRange = fi.find(line([1, 0, 1], [0, 1, 0]), makeCone(2, 3));
        expect(outOfRange.type).toBe(IntrLine3Cone3FIResultType.isEmpty);
    });

    it('validates reported segments and rays against the cone on random inputs', () => {
        let state = 19283746;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numSegments = 0;
        let numRays = 0;
        let numEmpty = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const finite = (trial % 3 === 0);
            const C = cone([rand(), rand(), rand()],
                [rand(), rand(), rand() + 0.001],
                0.2 + Math.abs(rand()) * 1.1,
                0, finite ? 1 + Math.abs(rand()) * 3 : -1);
            const L = line([rand() * 3, rand() * 3, rand() * 3],
                [rand(), rand(), rand() + 0.001]);
            const result = fi.find(L, C);

            const pointAt = (t: number) => add(L.origin, mul(t, L.direction));

            if (result.type === IntrLine3Cone3FIResultType.isSegment) {
                ++numSegments;
                const t0 = intrLine3Cone3Convert(result.t[0]);
                const t1 = intrLine3Cone3Convert(result.t[1]);
                expect(t0).toBeLessThanOrEqual(t1);
                // The reported endpoints match the parameters.
                const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
                expect(length(sub(P0, pointAt(t0)))).toBeCloseTo(0, 8);
                // The interior of the segment is in the solid cone and points
                // just outside are not.
                for (let k = 1; k < 5; ++k) {
                    const t = t0 + (t1 - t0) * k / 5;
                    expect(inSolidCone(C, pointAt(t), 1e-8)).toBe(true);
                }
                const span = Math.max(t1 - t0, 1e-3);
                expect(inSolidCone(C, pointAt(t0 - 0.05 * span), -1e-8))
                    .toBe(false);
                expect(inSolidCone(C, pointAt(t1 + 0.05 * span), -1e-8))
                    .toBe(false);
            }
            else if (result.type
                === IntrLine3Cone3FIResultType.isRayPositive) {
                ++numRays;
                const t0 = intrLine3Cone3Convert(result.t[0]);
                for (const t of [t0 + 0.5, t0 + 5, t0 + 50]) {
                    expect(inSolidCone(C, pointAt(t), 1e-7)).toBe(true);
                }
                expect(inSolidCone(C, pointAt(t0 - 0.5), -1e-8)).toBe(false);
            }
            else if (result.type
                === IntrLine3Cone3FIResultType.isRayNegative) {
                ++numRays;
                const t1 = intrLine3Cone3Convert(result.t[1]);
                for (const t of [t1 - 0.5, t1 - 5, t1 - 50]) {
                    expect(inSolidCone(C, pointAt(t), 1e-7)).toBe(true);
                }
                expect(inSolidCone(C, pointAt(t1 + 0.5), -1e-8)).toBe(false);
            }
            else if (result.type === IntrLine3Cone3FIResultType.isEmpty) {
                ++numEmpty;
                expect(result.intersect).toBe(false);
                // A dense sampling of the line finds no point in the cone.
                for (let k = -60; k <= 60; ++k) {
                    expect(inSolidCone(C, pointAt(k * 0.5), -1e-9)).toBe(false);
                }
            }
        }
        expect(numSegments).toBeGreaterThan(20);
        expect(numRays).toBeGreaterThan(10);
        expect(numEmpty).toBeGreaterThan(20);
    });

    it('is consistent under reversal of the line direction', () => {
        let state = 5555;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        for (let trial = 0; trial < 200; ++trial) {
            const C = cone([rand(), rand(), rand()],
                [rand(), rand(), rand() + 0.001],
                0.2 + Math.abs(rand()) * 1.1, 0, -1);
            const origin = [rand() * 3, rand() * 3, rand() * 3];
            const direction = [rand(), rand(), rand() + 0.001];
            const forward = fi.find(line(origin, direction), C);
            const backward = fi.find(line(origin,
                [-direction[0], -direction[1], -direction[2]]), C);
            expect(backward.intersect).toBe(forward.intersect);

            if (forward.type === IntrLine3Cone3FIResultType.isSegment) {
                expect(backward.type)
                    .toBe(IntrLine3Cone3FIResultType.isSegment);
                expect(intrLine3Cone3Convert(backward.t[0]))
                    .toBeCloseTo(-intrLine3Cone3Convert(forward.t[1]), 8);
                expect(intrLine3Cone3Convert(backward.t[1]))
                    .toBeCloseTo(-intrLine3Cone3Convert(forward.t[0]), 8);
            }
            else if (forward.type
                === IntrLine3Cone3FIResultType.isRayPositive) {
                expect(backward.type)
                    .toBe(IntrLine3Cone3FIResultType.isRayNegative);
                expect(intrLine3Cone3Convert(backward.t[1]))
                    .toBeCloseTo(-intrLine3Cone3Convert(forward.t[0]), 8);
            }
            else if (forward.type
                === IntrLine3Cone3FIResultType.isRayNegative) {
                expect(backward.type)
                    .toBe(IntrLine3Cone3FIResultType.isRayPositive);
                expect(intrLine3Cone3Convert(backward.t[0]))
                    .toBeCloseTo(-intrLine3Cone3Convert(forward.t[1]), 8);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrLine3Cone3.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrLine3Cone3 verification', () => {
    const fiQ = new IntrLine3Cone3FI();

    // The line generators use wellScaledVector: the shared 'vector' generator
    // emits subnormal components, and a line origin 5e-324 off the cone axis
    // underflows every coefficient of the quadratic (c0 = c1 = 0), which puts
    // the query in the exactly-through-the-vertex configuration.
    const arbWellLine = fc.tuple(wellScaledVector(3), unitVector(3))
        .map(([o, d]) => Line.fromOriginDirection(o, d));

    const arbAngle = fc.double({ min: 0.15, max: 1.35, noNaN: true });

    const arbInfiniteCone = fc.tuple(wellScaledVector(3), unitVector(3),
        arbAngle).map(([v, d, a]) =>
            Cone.fromRayAngle(Ray.fromOriginDirection(v, d), a));

    const arbFrustum = fc.tuple(wellScaledVector(3), unitVector(3), arbAngle,
        positive(2), positive(4, 0.2))
        .map(([v, d, a, hmin, extra]) => Cone.fromRayAngleMinMaxHeight(
            Ray.fromOriginDirection(v, d), a, hmin, hmin + extra));

    const arbFiniteCone = fc.tuple(wellScaledVector(3), unitVector(3),
        arbAngle, positive(5, 0.2))
        .map(([v, d, a, hmax]) => Cone.fromRayAngleMinMaxHeight(
            Ray.fromOriginDirection(v, d), a, 0, hmax));

    // The signed "insideness" of X relative to the (double-sided-free) solid
    // cone: G(X) = Dot(D, X-V) - cosAngle*|X-V|. It is nonnegative exactly on
    // the positive solid cone. The height clamp is applied separately.
    function coneG(cone: Cone, X: Vector): number {
        const d = sub(X, cone.ray.origin);
        return dot(cone.ray.direction, d) - cone.cosAngle * length(d);
    }

    function coneHeight(cone: Cone, X: Vector): number {
        return dot(cone.ray.direction, sub(X, cone.ray.origin));
    }

    function pointOnLine(L: Line, t: number): Vector {
        return add(L.origin, mul(t, L.direction));
    }

    // Does the line point at parameter t lie in the solid cone (with a
    // relative tolerance to absorb rounding near the boundary)?
    function inCone(cone: Cone, L: Line, t: number, tol: number): boolean {
        const X = pointOnLine(L, t);
        const scale = 1 + length(sub(X, cone.ray.origin));
        const h = coneHeight(cone, X);
        if (h < cone.getMinHeight() - tol * scale) { return false; }
        if (cone.isFinite() && h > cone.getMaxHeight() + tol * scale) {
            return false;
        }
        return coneG(cone, X) >= -tol * scale;
    }

    function tOf(q: QFNumber): number {
        return intrLine3Cone3Convert(q);
    }

    // The distance from the cone vertex to the line, relative to the size of
    // the configuration. A line through (or numerically through) the vertex
    // makes the discriminant of the quadratic exactly zero in exact
    // arithmetic, and rounding then routes the query into the discr > 0 or
    // discr < 0 branches, whose analysis assumes a nondegenerate
    // configuration. Those cases are excluded from the properties below and
    // pinned separately.
    function relVertexDistance(cone: Cone, L: Line): number {
        const PmV = sub(L.origin, cone.ray.origin);
        const perp = sub(PmV, mul(dot(PmV, L.direction), L.direction));
        return length(perp) / (1 + length(PmV));
    }

    it('the reported set is contained in the solid cone', () => {
        check(fc.tuple(arbWellLine,
            fc.oneof(arbInfiniteCone, arbFiniteCone, arbFrustum)),
            ([L, cone]) => {
                if (relVertexDistance(cone, L) < 1e-6) { return; }
                const r = fiQ.find(L, cone);
                expect(r.intersect).toBe(
                    r.type !== IntrLine3Cone3FIResultType.isEmpty);
                if (!r.intersect) { return; }
                const t0 = tOf(r.t[0]);
                const t1 = tOf(r.t[1]);
                const tol = 1e-7;
                switch (r.type) {
                    case IntrLine3Cone3FIResultType.isPoint:
                        expect(t1).toBe(t0);
                        expect(inCone(cone, L, t0, tol)).toBe(true);
                        break;
                    case IntrLine3Cone3FIResultType.isSegment:
                        expect(t0).toBeLessThanOrEqual(t1);
                        for (let i = 0; i <= 20; ++i) {
                            const t = t0 + ((t1 - t0) * i) / 20;
                            expect(inCone(cone, L, t, tol)).toBe(true);
                        }
                        break;
                    case IntrLine3Cone3FIResultType.isRayPositive:
                        for (let i = 0; i <= 20; ++i) {
                            expect(inCone(cone, L, t0 + i * 0.5, tol))
                                .toBe(true);
                        }
                        break;
                    case IntrLine3Cone3FIResultType.isRayNegative:
                        for (let i = 0; i <= 20; ++i) {
                            expect(inCone(cone, L, t1 - i * 0.5, tol))
                                .toBe(true);
                        }
                        break;
                    default:
                        break;
                }
            }, 100);
    });

    it('agrees with a dense sampling of the line against the solid cone', () => {
        check(fc.tuple(arbWellLine,
            fc.oneof(arbFiniteCone, arbFrustum)),
            ([L, cone]) => {
                if (relVertexDistance(cone, L) < 1e-6) { return; }
                // The finite cone is bounded, so a sampling of the line over
                // a window that contains it is conclusive.
                const V = cone.ray.origin;
                const reach = length(sub(V, L.origin))
                    + cone.getMaxHeight() / Math.max(cone.cosAngle, 1e-3) + 2;
                const N = 4000;
                let lo = Number.POSITIVE_INFINITY;
                let hi = Number.NEGATIVE_INFINITY;
                let count = 0;
                for (let i = 0; i <= N; ++i) {
                    const t = -reach + (2 * reach * i) / N;
                    if (inCone(cone, L, t, 0)) {
                        lo = Math.min(lo, t);
                        hi = Math.max(hi, t);
                        ++count;
                    }
                }
                const r = fiQ.find(L, cone);
                if (count === 0) { return; }  // sampling proves nothing
                expect(r.intersect).toBe(true);
                expect(r.type).not.toBe(
                    IntrLine3Cone3FIResultType.isRayPositive);
                expect(r.type).not.toBe(
                    IntrLine3Cone3FIResultType.isRayNegative);
                const t0 = tOf(r.t[0]);
                const t1 = tOf(r.t[1]);
                const step = (2 * reach) / N;
                // The reported interval contains the sampled hits.
                expect(t0).toBeLessThanOrEqual(lo + step);
                expect(t1).toBeGreaterThanOrEqual(hi - step);
            }, 60);
    });

    it('reports a ray for a line inside an infinite cone', () => {
        check(fc.tuple(arbInfiniteCone, fc.double({ min: 0.05, max: 0.9,
            noNaN: true }), positive(4, 0.1)),
            ([cone, frac, dist]) => {
                // A line through a point strictly inside the cone, parallel
                // to the axis: it enters the cone and stays inside.
                const V = cone.ray.origin;
                const D = cone.ray.direction;
                const helper = Math.abs(D.values[0]) < 0.9
                    ? Vector.fromArray([1, 0, 0]) : Vector.fromArray([0, 1, 0]);
                const u = cross(D, helper);
                normalize(u);
                const radius = frac * dist * cone.tanAngle;
                const P = add(add(V, mul(dist, D)), mul(radius, u));
                expect(coneG(cone, P)).toBeGreaterThan(0);
                const L = Line.fromOriginDirection(P, D);
                const r = fiQ.find(L, cone);
                expect(r.intersect).toBe(true);
                expect(r.type).toBe(IntrLine3Cone3FIResultType.isRayPositive);
                const t0 = tOf(r.t[0]);
                for (let i = 0; i <= 10; ++i) {
                    expect(inCone(cone, L, t0 + i, 1e-7)).toBe(true);
                }
                // Just before the ray origin the point is outside.
                expect(coneG(cone, pointOnLine(L, t0 - 1)))
                    .toBeLessThan(1e-9);
            });
    });

    it('reverses the reported set when the line direction is negated', () => {
        // DoQuery negates the direction when Dot(U,D) < 0 and mirrors the
        // t-values; the geometric set must be the same.
        check(fc.tuple(arbWellLine,
            fc.oneof(arbInfiniteCone, arbFiniteCone, arbFrustum)),
            ([L, cone]) => {
                const L2 = Line.fromOriginDirection(L.origin,
                    negate(L.direction));
                const r0 = fiQ.find(L, cone);
                const r1 = fiQ.find(L2, cone);
                expect(r1.intersect).toBe(r0.intersect);
                if (!r0.intersect) { return; }
                switch (r0.type) {
                    case IntrLine3Cone3FIResultType.isPoint:
                        expect(r1.type).toBe(IntrLine3Cone3FIResultType.isPoint);
                        expectClose(tOf(r1.t[0]), -tOf(r0.t[0]), 1e-9, 1e-9);
                        break;
                    case IntrLine3Cone3FIResultType.isSegment:
                        expect(r1.type)
                            .toBe(IntrLine3Cone3FIResultType.isSegment);
                        expectClose(tOf(r1.t[0]), -tOf(r0.t[1]), 1e-9, 1e-9);
                        expectClose(tOf(r1.t[1]), -tOf(r0.t[0]), 1e-9, 1e-9);
                        break;
                    case IntrLine3Cone3FIResultType.isRayPositive:
                        expect(r1.type)
                            .toBe(IntrLine3Cone3FIResultType.isRayNegative);
                        expectClose(tOf(r1.t[1]), -tOf(r0.t[0]), 1e-9, 1e-9);
                        break;
                    case IntrLine3Cone3FIResultType.isRayNegative:
                        expect(r1.type)
                            .toBe(IntrLine3Cone3FIResultType.isRayPositive);
                        expectClose(tOf(r1.t[0]), -tOf(r0.t[1]), 1e-9, 1e-9);
                        break;
                    default:
                        break;
                }
            }, 100);
    });

    it('the computed points match origin + t * direction', () => {
        check(fc.tuple(arbWellLine,
            fc.oneof(arbInfiniteCone, arbFiniteCone, arbFrustum)),
            ([L, cone]) => {
                const r = fiQ.find(L, cone);
                if (!r.intersect) {
                    for (let i = 0; i < 3; ++i) {
                        expect(intrLine3Cone3Convert(r.P[0][i])).toBe(0);
                        expect(intrLine3Cone3Convert(r.P[1][i])).toBe(0);
                    }
                    return;
                }
                const P0 = intrLine3Cone3ConvertPoint(r.P[0]);
                const scale = 1 + length(L.origin);
                if (r.type === IntrLine3Cone3FIResultType.isRayNegative) {
                    expectVectorClose(P0, pointOnLine(L, tOf(r.t[1])),
                        1e-9 * scale, 1e-9);
                    // P[1] is the line direction for the ray types.
                    expectVectorClose(intrLine3Cone3ConvertPoint(r.P[1]),
                        L.direction, 1e-12, 1e-12);
                } else if (r.type
                    === IntrLine3Cone3FIResultType.isRayPositive) {
                    expectVectorClose(P0, pointOnLine(L, tOf(r.t[0])),
                        1e-9 * scale, 1e-9);
                    expectVectorClose(intrLine3Cone3ConvertPoint(r.P[1]),
                        L.direction, 1e-12, 1e-12);
                } else {
                    expectVectorClose(P0, pointOnLine(L, tOf(r.t[0])),
                        1e-9 * scale, 1e-9);
                    expectVectorClose(intrLine3Cone3ConvertPoint(r.P[1]),
                        pointOnLine(L, tOf(r.t[1])), 1e-9 * scale, 1e-9);
                }
            }, 100);
    });

    it('takes the vertex branch only for a genuine vertex hit', () => {
        // The tangency at X != V is covered above; here the complementary
        // case is checked. A line THROUGH the vertex that is inside the
        // double cone (c2 > 0) must still take the vertex branch, which
        // reports a ray with origin V (Block 5 of the upstream algorithm).
        const axis = Ray.fromOriginDirection(Vector.zero(3),
            Vector.fromArray([0, 0, 1]));
        const infinite = Cone.fromRayAngle(axis, Math.PI / 4);
        infinite.cosAngleSqr = 0.5;
        const through = Line.fromOriginDirection(Vector.zero(3),
            Vector.fromArray([0, 0, 1]));
        const rv = fiQ.find(through, infinite);
        expect(rv.type).toBe(IntrLine3Cone3FIResultType.isRayPositive);
        expectVectorClose(intrLine3Cone3ConvertPoint(rv.P[0]), Vector.zero(3),
            1e-12, 1e-12);

        // The same line against a frustum is clamped to the height range.
        const frustum = Cone.fromRayAngleMinMaxHeight(axis, Math.PI / 4, 1, 3);
        frustum.cosAngleSqr = 0.5;
        const rf = fiQ.find(through, frustum);
        expect(rf.type).toBe(IntrLine3Cone3FIResultType.isSegment);
        expectClose(intrLine3Cone3Convert(rf.t[0]), 1, 1e-12, 1e-12);
        expectClose(intrLine3Cone3Convert(rf.t[1]), 3, 1e-12, 1e-12);

        // A line through the vertex that is OUTSIDE the double cone (c2 < 0)
        // meets it only at V (Block 4), which the height range can reject.
        const across = Line.fromOriginDirection(Vector.zero(3),
            Vector.fromArray([1, 0, 0]));
        expect(fiQ.find(across, infinite).type)
            .toBe(IntrLine3Cone3FIResultType.isPoint);
        expect(fiQ.find(across, frustum).type)
            .toBe(IntrLine3Cone3FIResultType.isEmpty);
    });


    it('pins the through-the-vertex discriminant fragility (upstream)', () => {
        // When the line passes through the cone vertex the quadratic
        // Q(t) = c2*t^2 + 2*c1*t + c0 has a double root, so discr = 0 in
        // exact arithmetic. discr = c1*c1 - c0*c2 is computed as the
        // difference of two nearly equal products, so rounding routes the
        // query into the discr > 0 or discr < 0 branches instead, whose case
        // analysis assumes a nondegenerate configuration:
        //
        // (a) c2 < 0 (the line direction is OUTSIDE the cone angle). The set
        //     Q >= 0 is the segment between the roots, which degenerates to
        //     the vertex. With discr rounded above zero the heights of the
        //     two roots straddle zero, so upstream takes Block 3 and calls
        //     SetRayClamp, whose analysis is valid only for c2 > 0. It
        //     reports a whole segment up to the cone's maximum height where
        //     the true intersection is the single vertex point.
        const coneA = Cone.fromRayAngleMinMaxHeight(Ray.fromOriginDirection(
            Vector.zero(3),
            Vector.fromArray([-0.9980401154430762, 0.06257737583480703, 0])),
            0.19085073139022773, 0, 0.20000000000000004);
        const lineA = Line.fromOriginDirection(
            Vector.fromArray([0, 3.7938468539782884, 0]),
            Vector.fromArray([0, 1, 0]));
        // The line contains the cone vertex exactly.
        expectVectorClose(add(lineA.origin,
            mul(-lineA.origin.values[1], lineA.direction)), Vector.zero(3),
            0, 0);
        const rA = fiQ.find(lineA, coneA);
        expect(rA.type).toBe(IntrLine3Cone3FIResultType.isSegment);
        // The far endpoint is 3.19 away from the vertex and well outside the
        // cone (G = -2.94); the correct answer is the vertex point alone.
        const far = pointOnLine(lineA, tOf(rA.t[1]));
        expect(length(sub(far, coneA.ray.origin))).toBeGreaterThan(3);
        expect(coneG(coneA, far)).toBeLessThan(-2);
        // Forcing the exactly-degenerate branch (discr == 0) gives the right
        // answer: the line meets the double cone only at the vertex, which is
        // in the height range [0, 0.2].
        expect(inCone(coneA, lineA, -3.7938468539782884, 1e-12)).toBe(true);

        // (b) c2 > 0 (the line is the cone axis, direction inside the cone).
        //     Rounding the same difference below zero takes the
        //     no-real-roots branch and reports no intersection at all, while
        //     the line covers the whole cone.
        const coneB = Cone.fromRayAngleMinMaxHeight(Ray.fromOriginDirection(
            Vector.fromArray([0.21849652746281664, 0, 0]),
            Vector.fromArray([1, 0, 0])), 0.878707887102512, 0,
            0.20000000000000004);
        const lineB = Line.fromOriginDirection(Vector.zero(3),
            Vector.fromArray([1, 0, 0]));
        const rB = fiQ.find(lineB, coneB);
        expect(rB.type).toBe(IntrLine3Cone3FIResultType.isEmpty);
        // The line really does run along the cone axis through the solid.
        expect(inCone(coneB, lineB, 0.3, 0)).toBe(true);
    });

    it('reports the empty result with zeroed members', () => {
        const d = defaultIntrLine3Cone3FIResult();
        expect(d.intersect).toBe(false);
        expect(d.type).toBe(IntrLine3Cone3FIResultType.isEmpty);
        // A line that misses the cone entirely.
        const cone = Cone.fromRayAngleMinMaxHeight(Ray.fromOriginDirection(
            Vector.zero(3), Vector.fromArray([0, 0, 1])), Math.PI / 6, 0, 1);
        const L = Line.fromOriginDirection(Vector.fromArray([0, 0, 10]),
            Vector.fromArray([1, 0, 0]));
        const r = fiQ.find(L, cone);
        expect(r.intersect).toBe(false);
        expect(r.type).toBe(IntrLine3Cone3FIResultType.isEmpty);
        expect(intrLine3Cone3Convert(r.t[0])).toBe(0);
        expect(intrLine3Cone3Convert(r.t[1])).toBe(0);
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(arbWellLine, arbFiniteCone, rotationFrame(3),
            wellScaledVector(3)),
            ([L, cone, Rot, t]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(Rot[0], v), dot(Rot[1], v), dot(Rot[2], v)]);
                const map = (v: Vector) => add(rot(v), t);
                const L2 = Line.fromOriginDirection(map(L.origin),
                    rot(L.direction));
                const cone2 = Cone.fromRayAngleMinMaxHeight(
                    Ray.fromOriginDirection(map(cone.ray.origin),
                        rot(cone.ray.direction)), cone.angle,
                    cone.getMinHeight(), cone.getMaxHeight());
                if (relVertexDistance(cone, L) < 1e-6) { return; }
                const r0 = fiQ.find(L, cone);
                const r1 = fiQ.find(L2, cone2);
                if (r0.type === IntrLine3Cone3FIResultType.isPoint
                    || r1.type === IntrLine3Cone3FIResultType.isPoint) {
                    return;  // tangency is a discontinuity
                }
                expect(r1.intersect).toBe(r0.intersect);
                if (r0.type === IntrLine3Cone3FIResultType.isSegment
                    && r1.type === IntrLine3Cone3FIResultType.isSegment) {
                    const scale = 1 + Math.abs(tOf(r0.t[0]))
                        + Math.abs(tOf(r0.t[1]));
                    expectClose(tOf(r1.t[0]), tOf(r0.t[0]), 1e-6 * scale,
                        1e-6);
                    expectClose(tOf(r1.t[1]), tOf(r0.t[1]), 1e-6 * scale,
                        1e-6);
                }
            }, 100);
    });
});
