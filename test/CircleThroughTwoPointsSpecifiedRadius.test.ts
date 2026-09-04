import { describe, it, expect } from 'vitest';
import { circleThroughTwoPointsSpecifiedRadius } from '../src/CircleThroughTwoPointsSpecifiedRadius.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { perp } from '../src/Vector2.js';
import {
    check, expectClose, expectVectorClose, fc, scaled, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md).
//
// The defining constraints are |C - P| = |C - Q| = r, so the properties check
// them directly and, separately, that the solution set is complete: the
// centers must be the only two points of the perpendicular bisector at
// distance r from P. The construction is a bisector parameterization plus one
// square root, so the relative error of the reported distances is a few ulps;
// the tolerances below are relative to the scale of the configuration.
// ---------------------------------------------------------------------------

const twoPointConfiguration = fc.record({
    P: wellScaledVector(2, -10, 10),
    Q: wellScaledVector(2, -10, 10),
    // The radius is expressed as a multiple of |P-Q|/2 so that all three
    // branches (two circles, one circle, none) are exercised.
    factor: scaled(0.1, 3, 64)
}).filter(({ P, Q }) => length(sub(P, Q)) > 1e-2);

describe('circleThroughTwoPointsSpecifiedRadius verification', () => {
    it('the reported circles satisfy the defining constraints', () => {
        check(twoPointConfiguration, ({ P, Q, factor }) => {
            const half = 0.5 * length(sub(P, Q));
            const r = factor * half;
            const { numCircles, circle } = circleThroughTwoPointsSpecifiedRadius(P, Q, r);
            const scale = Math.max(r, half, 1);
            for (let i = 0; i < numCircles; ++i) {
                expect(circle[i].radius).toBe(r);
                expectClose(length(sub(circle[i].center, P)), r, 1e-12 * scale, 1e-12);
                expectClose(length(sub(circle[i].center, Q)), r, 1e-12 * scale, 1e-12);
            }
            // Entries past numCircles are zeroed, as documented.
            for (let i = numCircles; i < 2; ++i) {
                expect(circle[i].radius).toBe(0);
                expect(circle[i].center.values).toEqual([0, 0]);
            }
        });
    });

    it('the count follows the r versus |P-Q|/2 comparison', () => {
        check(twoPointConfiguration, ({ P, Q, factor }) => {
            const half = 0.5 * length(sub(P, Q));
            const r = factor * half;
            const { numCircles } = circleThroughTwoPointsSpecifiedRadius(P, Q, r);
            // The branch is taken on the sign of r^2/|P-Q|^2 - 1/4, which for
            // factor away from 1 has the same sign as r - |P-Q|/2.
            if (Math.abs(factor - 1) > 1e-6) {
                expect(numCircles).toBe(r > half ? 2 : 0);
            } else {
                expect(numCircles).toBeGreaterThanOrEqual(1);
            }
            // Coincident points never produce a circle.
            expect(circleThroughTwoPointsSpecifiedRadius(P, P.clone(), r).numCircles)
                .toBe(0);
        });
    });

    it('the two centers are the only solutions on the bisector', () => {
        check(twoPointConfiguration, ({ P, Q, factor }) => {
            const half = 0.5 * length(sub(P, Q));
            const r = (1 + factor) * half;   // r > |P-Q|/2, so two circles
            const { numCircles, circle } = circleThroughTwoPointsSpecifiedRadius(P, Q, r);
            expect(numCircles).toBe(2);
            // Independent solution: the centers are M +- d*U where M is the
            // midpoint, U is the unit bisector direction and
            // d = sqrt(r^2 - |P-Q|^2/4).
            const M = mul(add(P, Q), 0.5);
            const U = perp(sub(P, Q));
            normalize(U);
            const d = Math.sqrt(r * r - half * half);
            const expected = [sub(M, mul(U, d)), add(M, mul(U, d))];
            const scale = Math.max(r, 1);
            // Upstream orders the centers as origin - root*Perp(P-Q) first.
            expectVectorClose(circle[0].center, expected[0], 1e-10 * scale, 1e-12);
            expectVectorClose(circle[1].center, expected[1], 1e-10 * scale, 1e-12);
            // The centers straddle the segment's line, so they are distinct.
            expect(dot(U, sub(circle[0].center, M)))
                .toBeLessThan(dot(U, sub(circle[1].center, M)));
        });
    });

    it('is symmetric under swapping P and Q (the centers swap)', () => {
        check(twoPointConfiguration, ({ P, Q, factor }) => {
            const r = (1 + factor) * 0.5 * length(sub(P, Q));
            const forward = circleThroughTwoPointsSpecifiedRadius(P, Q, r);
            const reverse = circleThroughTwoPointsSpecifiedRadius(Q, P, r);
            expect(reverse.numCircles).toBe(forward.numCircles);
            // Perp(Q-P) = -Perp(P-Q), so the roles of circle[0] and circle[1]
            // are exchanged. The midpoint is computed as 0.5*(P+Q) in both
            // calls, but floating-point addition is commutative, so the
            // centers agree exactly.
            expect(reverse.circle[0].center.values)
                .toEqual(forward.circle[1].center.values);
            expect(reverse.circle[1].center.values)
                .toEqual(forward.circle[0].center.values);
        });
    });

    it('is equivariant under rigid motions of the plane', () => {
        check(fc.tuple(twoPointConfiguration, scaled(-Math.PI, Math.PI, 32),
            wellScaledVector(2, -5, 5)), ([{ P, Q, factor }, angle, t]) => {
            const r = (1 + factor) * 0.5 * length(sub(P, Q));
            const c = Math.cos(angle), s = Math.sin(angle);
            const move = (v: Vector): Vector => Vector.fromArray([
                c * v.values[0] - s * v.values[1] + t.values[0],
                s * v.values[0] + c * v.values[1] + t.values[1]]);
            const direct = circleThroughTwoPointsSpecifiedRadius(P, Q, r);
            const moved = circleThroughTwoPointsSpecifiedRadius(move(P), move(Q), r);
            expect(moved.numCircles).toBe(direct.numCircles);
            const scale = Math.max(r, 10);
            for (let i = 0; i < direct.numCircles; ++i) {
                expectVectorClose(moved.circle[i].center,
                    move(direct.circle[i].center), 1e-10 * scale, 1e-12);
            }
        });
    });

    it('degenerate inputs return no circles with zeroed output', () => {
        check(fc.tuple(wellScaledVector(2, -5, 5), scaled(0, 5, 32)), ([P, r]) => {
            for (const result of [
                circleThroughTwoPointsSpecifiedRadius(P, P.clone(), r),
                circleThroughTwoPointsSpecifiedRadius(P, P.clone(), 0)]) {
                expect(result.numCircles).toBe(0);
                expect(result.circle[0].center.values).toEqual([0, 0]);
                expect(result.circle[0].radius).toBe(0);
                expect(result.circle[1].center.values).toEqual([0, 0]);
                expect(result.circle[1].radius).toBe(0);
            }
        });
        // r exactly |P-Q|/2 gives the single midpoint circle.
        const single = circleThroughTwoPointsSpecifiedRadius(v2(-3, 4), v2(1, 4), 2);
        expect(single.numCircles).toBe(1);
        expect(single.circle[0].center.values).toEqual([-1, 4]);
        expect(single.circle[0].radius).toBe(2);
        expect(single.circle[1].center.values).toEqual([0, 0]);
        expect(single.circle[1].radius).toBe(0);
    });
});
