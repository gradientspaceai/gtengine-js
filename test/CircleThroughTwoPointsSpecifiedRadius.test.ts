import { describe, it, expect } from 'vitest';
import { circleThroughTwoPointsSpecifiedRadius } from '../src/CircleThroughTwoPointsSpecifiedRadius.js';
import { Vector, length, sub } from '../src/Vector.js';

const v2 = (x: number, y: number): Vector => Vector.fromArray([x, y]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('circleThroughTwoPointsSpecifiedRadius', () => {
    it('returns two circles when r > |P-Q|/2 (hand-worked example)', () => {
        // P = (-1,0), Q = (1,0), r = sqrt(2). The bisector is the y-axis and
        // the centers are (0,-1) and (0,1); each is at distance sqrt(2) from
        // both P and Q.
        const P = v2(-1, 0);
        const Q = v2(1, 0);
        const r = Math.SQRT2;
        const result = circleThroughTwoPointsSpecifiedRadius(P, Q, r);
        expect(result.numCircles).toBe(2);
        expect(result.circle[0].center.values[0]).toBeCloseTo(0, 12);
        expect(result.circle[0].center.values[1]).toBeCloseTo(-1, 12);
        expect(result.circle[1].center.values[0]).toBeCloseTo(0, 12);
        expect(result.circle[1].center.values[1]).toBeCloseTo(1, 12);
        expect(result.circle[0].radius).toBe(r);
        expect(result.circle[1].radius).toBe(r);
    });

    it('returns one circle when r = |P-Q|/2, centered at the midpoint', () => {
        const P = v2(-1, 2);
        const Q = v2(1, 2);
        const result = circleThroughTwoPointsSpecifiedRadius(P, Q, 1);
        expect(result.numCircles).toBe(1);
        expect(result.circle[0].center.values[0]).toBeCloseTo(0, 12);
        expect(result.circle[0].center.values[1]).toBeCloseTo(2, 12);
        expect(result.circle[0].radius).toBe(1);

        // The unused circle has zero members.
        expect(result.circle[1].center.values).toEqual([0, 0]);
        expect(result.circle[1].radius).toBe(0);
    });

    it('returns no circles when r < |P-Q|/2', () => {
        const result = circleThroughTwoPointsSpecifiedRadius(v2(-1, 0), v2(1, 0), 0.5);
        expect(result.numCircles).toBe(0);
        for (const circle of result.circle) {
            expect(circle.center.values).toEqual([0, 0]);
            expect(circle.radius).toBe(0);
        }
    });

    it('returns no circles for a zero radius with distinct points', () => {
        const result = circleThroughTwoPointsSpecifiedRadius(v2(0, 0), v2(3, 4), 0);
        expect(result.numCircles).toBe(0);
    });

    it('returns no circles when P = Q (degenerate input)', () => {
        for (const r of [0, 1, 100]) {
            const result = circleThroughTwoPointsSpecifiedRadius(v2(2, 3), v2(2, 3), r);
            expect(result.numCircles).toBe(0);
            expect(result.circle[0].radius).toBe(0);
            expect(result.circle[1].radius).toBe(0);
        }
    });

    it('is symmetric in P and Q up to swapping the two centers', () => {
        const P = v2(1, 2);
        const Q = v2(4, -1);
        const r = 5;
        const forward = circleThroughTwoPointsSpecifiedRadius(P, Q, r);
        const reverse = circleThroughTwoPointsSpecifiedRadius(Q, P, r);
        expect(forward.numCircles).toBe(2);
        expect(reverse.numCircles).toBe(2);
        // Perp(P-Q) = -Perp(Q-P), so the two centers trade places.
        for (let i = 0; i < 2; ++i) {
            expect(forward.circle[i].center.values[0])
                .toBeCloseTo(reverse.circle[1 - i].center.values[0], 12);
            expect(forward.circle[i].center.values[1])
                .toBeCloseTo(reverse.circle[1 - i].center.values[1], 12);
        }
    });

    it('produces centers equidistant from P and Q at the given radius', () => {
        const random = makeRandom(20260830);
        for (let trial = 0; trial < 200; ++trial) {
            const P = v2(20 * random() - 10, 20 * random() - 10);
            const Q = v2(20 * random() - 10, 20 * random() - 10);
            const halfLength = 0.5 * length(sub(P, Q));
            // Choose r > |P-Q|/2 so that two circles exist.
            const r = halfLength + 0.1 + 5 * random();
            const result = circleThroughTwoPointsSpecifiedRadius(P, Q, r);
            expect(result.numCircles).toBe(2);
            for (let i = 0; i < 2; ++i) {
                const C = result.circle[i].center;
                expect(length(sub(C, P))).toBeCloseTo(r, 8);
                expect(length(sub(C, Q))).toBeCloseTo(r, 8);
                expect(result.circle[i].radius).toBe(r);
            }
            // The two centers are distinct and symmetric about the midpoint.
            const c0 = result.circle[0].center;
            const c1 = result.circle[1].center;
            expect(0.5 * (c0.values[0] + c1.values[0]))
                .toBeCloseTo(0.5 * (P.values[0] + Q.values[0]), 8);
            expect(0.5 * (c0.values[1] + c1.values[1]))
                .toBeCloseTo(0.5 * (P.values[1] + Q.values[1]), 8);
        }
    });

    it('returns no circles when the radius is too small (randomized)', () => {
        const random = makeRandom(777);
        for (let trial = 0; trial < 100; ++trial) {
            const P = v2(10 * random() - 5, 10 * random() - 5);
            const Q = v2(10 * random() - 5, 10 * random() - 5);
            const halfLength = 0.5 * length(sub(P, Q));
            const r = halfLength * (0.1 + 0.8 * random());
            expect(circleThroughTwoPointsSpecifiedRadius(P, Q, r).numCircles).toBe(0);
        }
    });

    it('throws when a point does not have size 2', () => {
        expect(() => circleThroughTwoPointsSpecifiedRadius(
            Vector.fromArray([0, 0, 0]), v2(1, 1), 2)).toThrow();
    });
});
