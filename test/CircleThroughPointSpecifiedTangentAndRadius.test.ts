import { describe, it, expect } from 'vitest';
import { circleThroughPointSpecifiedTangentAndRadius } from '../src/CircleThroughPointSpecifiedTangentAndRadius';
import { Vector, dot, length, normalize, sub } from '../src/Vector';

const v2 = (x: number, y: number): Vector => Vector.fromArray([x, y]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The signed distance from X to the line Dot(N,Y-A) = 0 for unit-length N.
function signedDistance(X: Vector, A: Vector, N: Vector): number {
    return dot(N, sub(X, A));
}

describe('circleThroughPointSpecifiedTangentAndRadius', () => {
    const A = v2(0, 0);
    const N = v2(0, 1);

    it('case (1): P on the line gives two circles tangent at P', () => {
        const P = v2(3, 0);
        const result = circleThroughPointSpecifiedTangentAndRadius(P, A, N, 2);
        expect(result.numCircles).toBe(2);
        expect(result.circle[0].center.values).toEqual([3, -2]);
        expect(result.circle[1].center.values).toEqual([3, 2]);
        expect(result.circle[0].radius).toBe(2);
        expect(result.circle[1].radius).toBe(2);
    });

    it('case (2): s = r gives two circles whose centers are P -/+ r*Perp(N)', () => {
        // Perp((0,1)) = (1,-0), so the centers are (-1,1) and (1,1).
        const P = v2(0, 1);
        const result = circleThroughPointSpecifiedTangentAndRadius(P, A, N, 1);
        expect(result.numCircles).toBe(2);
        expect(result.circle[0].center.values[0]).toBeCloseTo(-1, 12);
        expect(result.circle[0].center.values[1]).toBeCloseTo(1, 12);
        expect(result.circle[1].center.values[0]).toBeCloseTo(1, 12);
        expect(result.circle[1].center.values[1]).toBeCloseTo(1, 12);
    });

    it('case (3): s = 2*r gives a single circle with center P - r*N', () => {
        const P = v2(0, 2);
        const result = circleThroughPointSpecifiedTangentAndRadius(P, A, N, 1);
        expect(result.numCircles).toBe(1);
        expect(result.circle[0].center.values).toEqual([0, 1]);
        expect(result.circle[0].radius).toBe(1);
        expect(result.circle[1].center.values).toEqual([0, 0]);
        expect(result.circle[1].radius).toBe(0);
    });

    it('case (4): s > 2*r gives no circles', () => {
        const result = circleThroughPointSpecifiedTangentAndRadius(v2(0, 5), A, N, 1);
        expect(result.numCircles).toBe(0);
        for (const circle of result.circle) {
            expect(circle.center.values).toEqual([0, 0]);
            expect(circle.radius).toBe(0);
        }
    });

    it('case (5a): 0 < s < r gives two circles (hand-worked example)', () => {
        // P = (0,1), r = 2, line y = 0. The bisector origin is (0,2) with
        // direction Perp(N) = (1,0) and h = sqrt(r^2 - (r-s)^2) = sqrt(3).
        const P = v2(0, 1);
        const result = circleThroughPointSpecifiedTangentAndRadius(P, A, N, 2);
        expect(result.numCircles).toBe(2);
        expect(result.circle[0].center.values[0]).toBeCloseTo(-Math.sqrt(3), 12);
        expect(result.circle[0].center.values[1]).toBeCloseTo(2, 12);
        expect(result.circle[1].center.values[0]).toBeCloseTo(Math.sqrt(3), 12);
        expect(result.circle[1].center.values[1]).toBeCloseTo(2, 12);
    });

    it('case (5b): r < s < 2*r gives two circles', () => {
        // P = (0,3), r = 2, line y = 0. B = (0,2), (r-s) = -1, h = sqrt(3).
        const P = v2(0, 3);
        const result = circleThroughPointSpecifiedTangentAndRadius(P, A, N, 2);
        expect(result.numCircles).toBe(2);
        expect(result.circle[0].center.values[0]).toBeCloseTo(-Math.sqrt(3), 12);
        expect(result.circle[0].center.values[1]).toBeCloseTo(2, 12);
        expect(result.circle[1].center.values[0]).toBeCloseTo(Math.sqrt(3), 12);
        expect(result.circle[1].center.values[1]).toBeCloseTo(2, 12);
    });

    it('mirrors the result when P is on the negative side of the line', () => {
        const above = circleThroughPointSpecifiedTangentAndRadius(v2(0, 1), A, N, 2);
        const below = circleThroughPointSpecifiedTangentAndRadius(v2(0, -1), A, N, 2);
        expect(below.numCircles).toBe(2);
        // P = (0,-1) is the reflection of (0,1) through the origin, which is
        // on the line, so the whole configuration is reflected through the
        // origin and the centers negate. The circles do not trade places:
        // both the normal and the bisector direction Perp(N) negate.
        for (let i = 0; i < 2; ++i) {
            expect(below.circle[i].center.values[0])
                .toBeCloseTo(-above.circle[i].center.values[0], 12);
            expect(below.circle[i].center.values[1])
                .toBeCloseTo(-above.circle[i].center.values[1], 12);
        }
    });

    it('does not modify the caller normal when P is below the line', () => {
        const normal = v2(0, 1);
        circleThroughPointSpecifiedTangentAndRadius(v2(0, -1), A, normal, 2);
        expect(normal.values).toEqual([0, 1]);
    });

    it('gives centers at distance r from both P and the tangent line', () => {
        const random = makeRandom(19981998);
        let numTwo = 0;
        let numZero = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const origin = v2(6 * random() - 3, 6 * random() - 3);
            const normal = v2(2 * random() - 1, 2 * random() - 1);
            if (length(normal) < 1e-3) {
                continue;
            }
            normalize(normal);
            const P = v2(10 * random() - 5, 10 * random() - 5);
            const r = 0.25 + 4 * random();
            const result = circleThroughPointSpecifiedTangentAndRadius(P, origin, normal, r);

            const s = Math.abs(signedDistance(P, origin, normal));
            if (s > 2 * r) {
                expect(result.numCircles).toBe(0);
                ++numZero;
                continue;
            }

            expect(result.numCircles).toBeGreaterThan(0);
            if (result.numCircles === 2) {
                ++numTwo;
            }
            for (let i = 0; i < result.numCircles; ++i) {
                const C = result.circle[i].center;
                // The circle passes through P.
                expect(length(sub(C, P))).toBeCloseTo(r, 7);
                // The line is tangent: the distance from C to the line is r.
                expect(Math.abs(signedDistance(C, origin, normal))).toBeCloseTo(r, 7);
                expect(result.circle[i].radius).toBe(r);
            }
        }
        // The random sample exercises both the solution and no-solution paths.
        expect(numTwo).toBeGreaterThan(0);
        expect(numZero).toBeGreaterThan(0);
    });

    it('is unaffected by the sign of the input normal', () => {
        const P = v2(1, 1.5);
        const r = 2;
        const positive = circleThroughPointSpecifiedTangentAndRadius(P, A, v2(0, 1), r);
        const negative = circleThroughPointSpecifiedTangentAndRadius(P, A, v2(0, -1), r);
        expect(positive.numCircles).toBe(negative.numCircles);
        for (let i = 0; i < 2; ++i) {
            expect(positive.circle[i].center.values[0])
                .toBeCloseTo(negative.circle[i].center.values[0], 10);
            expect(positive.circle[i].center.values[1])
                .toBeCloseTo(negative.circle[i].center.values[1], 10);
        }
    });

    it('throws when an input does not have size 2', () => {
        expect(() => circleThroughPointSpecifiedTangentAndRadius(
            Vector.fromArray([0, 0, 0]), A, N, 1)).toThrow();
    });
});
