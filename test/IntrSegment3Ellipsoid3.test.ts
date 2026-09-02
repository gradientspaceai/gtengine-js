import { describe, it, expect } from 'vitest';
import { Hyperellipsoid } from '../src/Hyperellipsoid';
import {
    IntrSegment3Ellipsoid3TI,
    IntrSegment3Ellipsoid3FI,
    defaultIntrSegment3Ellipsoid3FIResult,
    intrSegment3Ellipsoid3FIDoQuery
} from '../src/IntrSegment3Ellipsoid3';
import { Segment } from '../src/Segment';
import { Vector, add, length, mul, sub } from '../src/Vector';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function ellipsoid(center: number[], extent: number[]): Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(Vector.fromArray(center),
        [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)],
        Vector.fromArray(extent));
}

// The value of (X-C)^T M (X-C) - 1: zero on the ellipsoid, negative inside.
function level(E: Hyperellipsoid, X: Vector): number {
    const d = sub(X, E.center);
    let sum = -1;
    for (let i = 0; i < 3; ++i) {
        const t = d.values[i] / E.extent.values[i];
        sum += t * t;
    }
    return sum;
}

// An independent test: sample the segment densely and look for a point in the
// solid ellipsoid.
function bruteForceIntersect(S: Segment, E: Hyperellipsoid): boolean {
    const delta = sub(S.p[1], S.p[0]);
    const n = 4000;
    for (let i = 0; i <= n; ++i) {
        if (level(E, add(S.p[0], mul(i / n, delta))) <= 0) {
            return true;
        }
    }
    return false;
}

const ti = new IntrSegment3Ellipsoid3TI();
const fi = new IntrSegment3Ellipsoid3FI();

describe('IntrSegment3Ellipsoid3', () => {
    it('has an empty default result', () => {
        const result = defaultIntrSegment3Ellipsoid3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('finds both crossings of a straddling segment', () => {
        // The segment spans x in [-10,10] at y = z = 0; its centered form has
        // origin (0,0,0), direction (1,0,0) and extent 10, so the crossings
        // are at parameters -3 and 3.
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const S = seg([-10, 0, 0], [10, 0, 0]);
        const result = fi.find(S, E);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-3, 12);
        expect(result.parameter[1]).toBeCloseTo(3, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-3, 12);
        expect(result.point[1].values[0]).toBeCloseTo(3, 12);
        expect(ti.test(S, E).intersect).toBe(true);
    });

    it('clips a crossing when an endpoint is inside the ellipsoid', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        // The segment from the center to (10,0,0): centered form origin
        // (5,0,0), direction (1,0,0), extent 5, so the ellipsoid crossing at
        // x = 3 is at parameter -2 and the clipped near end is at -5.
        const S = seg([0, 0, 0], [10, 0, 0]);
        const result = fi.find(S, E);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-5, 12);
        expect(result.parameter[1]).toBeCloseTo(-2, 12);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(3, 12);
        expect(ti.test(S, E).intersect).toBe(true);
    });

    it('reports no intersection for a segment that stops short', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const S = seg([-10, 0, 0], [-4, 0, 0]);
        expect(fi.find(S, E).intersect).toBe(false);
        expect(ti.test(S, E).intersect).toBe(false);
        // A segment whose supporting line misses the ellipsoid.
        const M = seg([-10, 5, 0], [10, 5, 0]);
        expect(fi.find(M, E).intersect).toBe(false);
        expect(ti.test(M, E).intersect).toBe(false);
    });

    it('reports a segment fully inside the ellipsoid', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const S = seg([-1, 0, 0], [1, 0, 0]);
        const result = fi.find(S, E);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(ti.test(S, E).intersect).toBe(true);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const result = defaultIntrSegment3Ellipsoid3FIResult();
        intrSegment3Ellipsoid3FIDoQuery(vec(0, 0, 0), vec(1, 0, 0), 10,
            ellipsoid([0, 0, 0], [3, 2, 1]), result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-3, 12);
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('agrees with the TI query and a dense sampling on random inputs', () => {
        let state = 24680;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numHits = 0;
        for (let trial = 0; trial < 150; ++trial) {
            const E = ellipsoid([rand(), rand(), rand()],
                [0.5 + Math.abs(rand()) * 2, 0.5 + Math.abs(rand()) * 2,
                    0.5 + Math.abs(rand()) * 2]);
            const S = seg([rand() * 4, rand() * 4, rand() * 4],
                [rand() * 4, rand() * 4, rand() * 4]);
            const result = fi.find(S, E);
            expect(ti.test(S, E).intersect).toBe(result.intersect);
            expect(result.intersect).toBe(bruteForceIntersect(S, E));

            if (result.intersect) {
                ++numHits;
                const { center: segOrigin, direction: segDirection,
                    extent: segExtent } = S.getCenteredForm();
                for (let i = 0; i < result.numIntersections; ++i) {
                    expect(Math.abs(result.parameter[i]))
                        .toBeLessThanOrEqual(segExtent + 1e-12);
                    const onSegment = add(segOrigin,
                        mul(result.parameter[i], segDirection));
                    expect(length(sub(result.point[i], onSegment)))
                        .toBeCloseTo(0, 10);
                    // The point is on or inside the ellipsoid; it is on the
                    // boundary unless it is a clipped segment endpoint.
                    expect(level(E, result.point[i])).toBeLessThan(1e-8);
                    const clipped =
                        Math.abs(Math.abs(result.parameter[i]) - segExtent)
                        < 1e-14;
                    if (!clipped) {
                        expect(level(E, result.point[i])).toBeCloseTo(0, 8);
                    }
                }
            }
        }
        expect(numHits).toBeGreaterThan(10);
    });
});
