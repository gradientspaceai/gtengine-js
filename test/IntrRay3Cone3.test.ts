import { describe, it, expect } from 'vitest';
import { Cone } from '../src/Cone';
import {
    IntrLine3Cone3FI,
    IntrLine3Cone3FIResultType,
    intrLine3Cone3Convert,
    intrLine3Cone3ConvertPoint
} from '../src/IntrLine3Cone3';
import {
    IntrRay3Cone3FI,
    defaultIntrRay3Cone3FIResult
} from '../src/IntrRay3Cone3';
import { Line } from '../src/Line';
import { Ray } from '../src/Ray';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector';

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
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

const fi = new IntrRay3Cone3FI();
const lineFI = new IntrLine3Cone3FI();

const T = IntrLine3Cone3FIResultType;
const quarterPi = Math.PI / 4;

describe('IntrRay3Cone3', () => {
    it('default-constructs the result as empty', () => {
        const r = defaultIntrRay3Cone3FIResult();
        expect(r.intersect).toBe(false);
        expect(r.type).toBe(T.isEmpty);
    });

    it('finds the known segment for a ray crossing a finite cone', () => {
        // Cone with apex at the origin, axis +z, half-angle pi/4, heights
        // [0,4]. The ray starts at (-5,0,2) heading in +x; the cone cross
        // section at z = 2 is the disk of radius 2, so the segment endpoints
        // are (-2,0,2) and (2,0,2).
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const R = ray([-5, 0, 2], [1, 0, 0]);
        const result = fi.find(R, C);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(T.isSegment);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
        expect(P0.values[0]).toBeCloseTo(-2, 10);
        expect(P1.values[0]).toBeCloseTo(2, 10);
        expect(intrLine3Cone3Convert(result.t[0])).toBeCloseTo(3, 10);
        expect(intrLine3Cone3Convert(result.t[1])).toBeCloseTo(7, 10);
    });

    it('clips the line result at the ray origin (block 14)', () => {
        // The same configuration but the ray origin is inside the cone, so
        // the line segment [3,7] is clipped to [0,7] in ray parameters.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const R = ray([0, 0, 2], [1, 0, 0]);
        const result = fi.find(R, C);
        expect(result.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(result.t[0])).toBeCloseTo(0, 12);
        expect(intrLine3Cone3Convert(result.t[1])).toBeCloseTo(2, 10);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        expect(P0.values[0]).toBeCloseTo(0, 12);
    });

    it('reports empty when the intersection is entirely behind the origin', () => {
        // Block 15: the whole line-cone segment has t < 0.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const R = ray([5, 0, 2], [1, 0, 0]);
        const result = fi.find(R, C);
        expect(result.intersect).toBe(false);
        expect(result.type).toBe(T.isEmpty);
    });

    it('clips to a sliver when the origin is just inside the boundary', () => {
        // The cone cross section at z = 2 is the disk of radius 2, so a ray
        // starting just inside the boundary and heading out exits almost
        // immediately. An exactly-on-the-boundary origin (upstream block 16)
        // is not representable in floating point for this cone.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const inside = fi.find(ray([2 - 1e-6, 0, 2], [1, 0, 0]), C);
        expect(inside.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(inside.t[0])).toBeCloseTo(0, 12);
        expect(intrLine3Cone3Convert(inside.t[1])).toBeCloseTo(1e-6, 12);

        const outside = fi.find(ray([2 + 1e-6, 0, 2], [1, 0, 0]), C);
        expect(outside.intersect).toBe(false);
        expect(outside.type).toBe(T.isEmpty);
    });

    it('keeps a positive ray for an infinite cone (block 17)', () => {
        // The ray starts on the axis inside an infinite cone, so the whole
        // ray is inside.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1);
        const R = ray([0, 0, 1], [0, 0, 1]);
        const result = fi.find(R, C);
        expect(result.type).toBe(T.isRayPositive);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        expect(P0.values[2]).toBeCloseTo(1, 12);
        expect(intrLine3Cone3Convert(result.t[0])).toBeCloseTo(0, 12);
    });

    it('turns a negative ray into a segment or a point (blocks 18-20)', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1);

        // Block 18: the ray travels down the axis from inside the cone, so
        // it exits at the apex.
        const R18 = ray([0, 0, 1], [0, 0, -1]);
        const r18 = fi.find(R18, C);
        expect(r18.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(r18.t[0])).toBeCloseTo(0, 12);
        expect(intrLine3Cone3Convert(r18.t[1])).toBeCloseTo(1, 10);
        const Q1 = intrLine3Cone3ConvertPoint(r18.P[1]);
        expect(length(Q1)).toBeCloseTo(0, 10);

        // Block 19: below the apex heading further down, nothing is hit.
        const R19 = ray([0, 0, -1], [0, 0, -1]);
        expect(fi.find(R19, C).type).toBe(T.isEmpty);

        // Block 20: starting exactly at the apex heading down touches only
        // the apex.
        const R20 = ray([0, 0, 0], [0, 0, -1]);
        const r20 = fi.find(R20, C);
        expect(r20.type).toBe(T.isPoint);
        const P0 = intrLine3Cone3ConvertPoint(r20.P[0]);
        expect(length(P0)).toBeCloseTo(0, 12);
    });

    it('keeps a ray that starts off-axis inside an infinite cone', () => {
        // The line x = 0.5, y = 0 enters the 45-degree cone at z = 0.5, so a
        // downward ray from (0.5,0,3) is clipped to the segment ending there.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1);
        const down = fi.find(ray([0.5, 0, 3], [0, 0, -1]), C);
        expect(down.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(down.t[0])).toBeCloseTo(0, 12);
        expect(intrLine3Cone3Convert(down.t[1])).toBeCloseTo(2.5, 8);

        // The same line traversed upward from inside is an unbounded ray.
        const up = fi.find(ray([0.5, 0, 3], [0, 0, 1]), C);
        expect(up.type).toBe(T.isRayPositive);
        expect(intrLine3Cone3Convert(up.t[0])).toBeCloseTo(0, 12);
    });

    it('clamps to the cone height range for a cone frustum', () => {
        // A frustum with heights [1,3]. A ray down the axis from above must
        // produce the segment between the caps.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 1, 3);
        const R = ray([0, 0, 10], [0, 0, -1]);
        const result = fi.find(R, C);
        expect(result.type).toBe(T.isSegment);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
        const zs = [P0.values[2], P1.values[2]].sort((a, b) => a - b);
        expect(zs[0]).toBeCloseTo(1, 10);
        expect(zs[1]).toBeCloseTo(3, 10);
    });

    it('agrees with the line query on the ray portion (randomized)', () => {
        let seed = 13572468;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (): number => 6 * rand() - 3;

        const cones = [
            cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1),
            cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 3),
            cone([0.5, -0.25, 0], [1, 1, 1], 0.6, 0.5, 2.5),
            cone([0, 0, 0], [0, 1, 0], 0.4, 0.25, -1)
        ];

        let numSegments = 0;
        for (let trial = 0; trial < 2000; ++trial) {
            const C = cones[trial % cones.length];
            const R = ray([rnd(), rnd(), rnd()], [rnd(), rnd(), rnd()]);
            const result = fi.find(R, C);
            const lineResult = lineFI.find(
                Line.fromOriginDirection(R.origin, R.direction), C);

            if (result.intersect) {
                // The ray result is never larger than the line result.
                expect(lineResult.intersect).toBe(true);
            }

            if (result.type === T.isSegment) {
                ++numSegments;
                const t0 = intrLine3Cone3Convert(result.t[0]);
                const t1 = intrLine3Cone3Convert(result.t[1]);
                expect(t0).toBeGreaterThanOrEqual(-1e-12);
                expect(t1).toBeGreaterThanOrEqual(t0 - 1e-12);
                const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
                const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
                expect(inSolidCone(C, P0, 1e-8)).toBe(true);
                expect(inSolidCone(C, P1, 1e-8)).toBe(true);
                // The midpoint is inside the solid cone too (convexity).
                const mid = mul(0.5, add(P0, P1));
                expect(inSolidCone(C, mid, 1e-8)).toBe(true);
                // The endpoints are on the ray.
                for (const [t, P] of [[t0, P0], [t1, P1]] as [number, Vector][]) {
                    const X = add(R.origin, mul(t, R.direction));
                    for (let i = 0; i < 3; ++i) {
                        expect(P.values[i]).toBeCloseTo(X.values[i], 8);
                    }
                }
            }
            else if (result.type === T.isPoint) {
                const t0 = intrLine3Cone3Convert(result.t[0]);
                expect(t0).toBeGreaterThanOrEqual(-1e-12);
                const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
                expect(inSolidCone(C, P0, 1e-8)).toBe(true);
            }
            else if (result.type === T.isRayPositive) {
                const t0 = intrLine3Cone3Convert(result.t[0]);
                expect(t0).toBeGreaterThanOrEqual(-1e-12);
                const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
                expect(inSolidCone(C, P0, 1e-8)).toBe(true);
                // Walking along the ray direction stays inside.
                const far = add(P0, mul(100, R.direction));
                expect(inSolidCone(C, far, 1e-6)).toBe(true);
            }
            else {
                expect(result.type).toBe(T.isEmpty);
            }
        }
        expect(numSegments).toBeGreaterThan(50);
    });
});
