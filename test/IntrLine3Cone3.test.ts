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
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';

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
