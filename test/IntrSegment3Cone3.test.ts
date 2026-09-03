import { describe, it, expect } from 'vitest';
import { Cone } from '../src/Cone.js';
import {
    IntrLine3Cone3FI,
    IntrLine3Cone3FIResultType,
    intrLine3Cone3Convert,
    intrLine3Cone3ConvertPoint
} from '../src/IntrLine3Cone3.js';
import {
    IntrSegment3Cone3FI,
    defaultIntrSegment3Cone3FIResult
} from '../src/IntrSegment3Cone3.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
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

const fi = new IntrSegment3Cone3FI();
const lineFI = new IntrLine3Cone3FI();

const T = IntrLine3Cone3FIResultType;
const quarterPi = Math.PI / 4;

describe('IntrSegment3Cone3', () => {
    it('default-constructs the result as empty', () => {
        const r = defaultIntrSegment3Cone3FIResult();
        expect(r.intersect).toBe(false);
        expect(r.type).toBe(T.isEmpty);
    });

    it('finds the known chord of a finite cone', () => {
        // The cone cross section at z = 2 is the disk of radius 2. The
        // segment from (-5,0,2) to (5,0,2) is parameterized over [0,1], so
        // the chord is t in [0.3,0.7].
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const S = segment([-5, 0, 2], [5, 0, 2]);
        const result = fi.find(S, C);
        expect(result.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(result.t[0])).toBeCloseTo(0.3, 10);
        expect(intrLine3Cone3Convert(result.t[1])).toBeCloseTo(0.7, 10);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
        expect(P0.values[0]).toBeCloseTo(-2, 10);
        expect(P1.values[0]).toBeCloseTo(2, 10);
        expect(P0.values[2]).toBeCloseTo(2, 12);
    });

    it('clips both ends to the segment interval (block 24)', () => {
        // The segment is strictly inside the chord, so the result is the
        // whole segment [0,1].
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const S = segment([-1, 0, 2], [1, 0, 2]);
        const result = fi.find(S, C);
        expect(result.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(result.t[0])).toBeCloseTo(0, 12);
        expect(intrLine3Cone3Convert(result.t[1])).toBeCloseTo(1, 12);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
        expect(P0.values[0]).toBeCloseTo(-1, 12);
        expect(P1.values[0]).toBeCloseTo(1, 12);
    });

    it('reports empty for segments on either side of the chord (block 23)', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        expect(fi.find(segment([-9, 0, 2], [-5, 0, 2]), C).type).toBe(T.isEmpty);
        expect(fi.find(segment([5, 0, 2], [9, 0, 2]), C).type).toBe(T.isEmpty);
    });

    it('reports empty for a segment entirely below the cone', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, 4);
        const result = fi.find(segment([-3, -3, -2], [3, 3, -2]), C);
        expect(result.intersect).toBe(false);
    });

    it('clips an unbounded line result to the segment (blocks 27 and 30)', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1);

        // Block 30: the line result is a negative ray ending at the apex. The
        // segment runs downward from (0.5,0,3) past the entry point z = 0.5.
        const down = fi.find(segment([0.5, 0, 3], [0.5, 0, -3]), C);
        expect(down.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(down.t[0])).toBeCloseTo(0, 12);
        // The segment direction has length 6 and the entry is 2.5 below the
        // start, so t = 2.5/6.
        expect(intrLine3Cone3Convert(down.t[1])).toBeCloseTo(2.5 / 6, 8);

        // Block 27: the line result is a positive ray, and the segment runs
        // upward from inside, so the whole segment is the intersection.
        const up = fi.find(segment([0.5, 0, 3], [0.5, 0, 9]), C);
        expect(up.type).toBe(T.isSegment);
        expect(intrLine3Cone3Convert(up.t[0])).toBeCloseTo(0, 12);
        expect(intrLine3Cone3Convert(up.t[1])).toBeCloseTo(1, 12);
    });

    it('reports empty when an unbounded line result misses the segment', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1);
        // Block 29: the negative-ray line result ends at the apex, well
        // before the segment starts.
        const result = fi.find(segment([0.5, 0, -1], [0.5, 0, -5]), C);
        expect(result.type).toBe(T.isEmpty);
    });

    it('touches the apex in a single point (block 31)', () => {
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 0, -1);
        const result = fi.find(segment([0, 0, 0], [0, 0, -4]), C);
        expect(result.type).toBe(T.isPoint);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        expect(length(P0)).toBeCloseTo(0, 12);
    });

    it('clamps to the cone frustum height range', () => {
        // The frustum has heights [1,3]; a segment down the axis from above
        // must be clipped to the two cap planes.
        const C = cone([0, 0, 0], [0, 0, 1], quarterPi, 1, 3);
        const result = fi.find(segment([0, 0, 10], [0, 0, -10]), C);
        expect(result.type).toBe(T.isSegment);
        const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
        const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
        const zs = [P0.values[2], P1.values[2]].sort((a, b) => a - b);
        expect(zs[0]).toBeCloseTo(1, 10);
        expect(zs[1]).toBeCloseTo(3, 10);
    });

    it('agrees with the line query and stays in the segment (randomized)', () => {
        let seed = 55501234;
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
            const p0 = [rnd(), rnd(), rnd()];
            const p1 = [rnd(), rnd(), rnd()];
            const S = segment(p0, p1);
            const result = fi.find(S, C);

            // The result is never a ray: the segment is bounded.
            expect(result.type).not.toBe(T.isRayPositive);
            expect(result.type).not.toBe(T.isRayNegative);

            if (result.intersect) {
                // The line query must find at least as much.
                const dir = sub(S.p[1], S.p[0]);
                const lineResult = lineFI.find(
                    Line.fromOriginDirection(S.p[0], dir), C);
                expect(lineResult.intersect).toBe(true);
            }

            if (result.type === T.isSegment) {
                ++numSegments;
                const t0 = intrLine3Cone3Convert(result.t[0]);
                const t1 = intrLine3Cone3Convert(result.t[1]);
                expect(t0).toBeGreaterThanOrEqual(-1e-12);
                expect(t1).toBeLessThanOrEqual(1 + 1e-12);
                expect(t1).toBeGreaterThanOrEqual(t0 - 1e-12);
                const P0 = intrLine3Cone3ConvertPoint(result.P[0]);
                const P1 = intrLine3Cone3ConvertPoint(result.P[1]);
                expect(inSolidCone(C, P0, 1e-8)).toBe(true);
                expect(inSolidCone(C, P1, 1e-8)).toBe(true);
                expect(inSolidCone(C, mul(0.5, add(P0, P1)), 1e-8)).toBe(true);
                // The endpoints are on the segment.
                const dir = sub(S.p[1], S.p[0]);
                for (const [t, P] of [[t0, P0], [t1, P1]] as [number, Vector][]) {
                    const X = add(S.p[0], mul(t, dir));
                    for (let i = 0; i < 3; ++i) {
                        expect(P.values[i]).toBeCloseTo(X.values[i], 8);
                    }
                }
            }
            else if (result.type === T.isPoint) {
                const t0 = intrLine3Cone3Convert(result.t[0]);
                expect(t0).toBeGreaterThanOrEqual(-1e-12);
                expect(t0).toBeLessThanOrEqual(1 + 1e-12);
                expect(inSolidCone(C,
                    intrLine3Cone3ConvertPoint(result.P[0]), 1e-8)).toBe(true);
            }
            else {
                expect(result.type).toBe(T.isEmpty);
            }
        }
        expect(numSegments).toBeGreaterThan(30);
    });
});
