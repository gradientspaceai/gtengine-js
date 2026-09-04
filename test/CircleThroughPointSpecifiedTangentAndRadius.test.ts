import { describe, it, expect } from 'vitest';
import { circleThroughPointSpecifiedTangentAndRadius } from '../src/CircleThroughPointSpecifiedTangentAndRadius.js';
import { Vector, add, dot, length, mul, negate, normalize, sub } from '../src/Vector.js';
import { perp } from '../src/Vector2.js';
import {
    check, expectClose, fc, scaled, wellScaledVector
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

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md).
//
// The reference below is derived from the constraint set rather than from the
// implementation: write a candidate center as C = A + x*D + y*N with
// D = Perp(N). Tangency to the line Dot(N,X-A) = 0 forces |y| = r, and
// |C - P| = r with P = A + u*D + s*N forces (x-u)^2 = r^2 - (y-s)^2. For
// s >= 0 the y = -r branch is solvable only when s = 0, so the complete
// solution set is
//   s = 0        : {P - r*N, P + r*N}
//   0 < s < 2*r  : {B - h*D, B + h*D}, B = A + u*D + r*N, h = sqrt(r^2-(r-s)^2)
//   s = 2*r      : {B}
//   s > 2*r      : {}
// which is exactly what the function must report, including its case (2)
// (s = r, where B = P).
//
// The generator keeps s away from 2*r because h = sqrt(r^2 - (r-s)^2) loses
// half its significant digits as s approaches 2*r; the exact boundary cases
// are pinned separately with axis-aligned configurations where the branch
// tests are exact.
// ---------------------------------------------------------------------------

interface TangentSolution { centers: Vector[] }

function referenceCenters(P: Vector, A: Vector, N: Vector, r: number): TangentSolution {
    let n = N.clone();
    let s = dot(n, sub(P, A));
    if (s === 0) {
        return { centers: [sub(P, mul(n, r)), add(P, mul(n, r))] };
    }
    if (s < 0) {
        n = negate(n);
        s = -s;
    }
    if (s > 2 * r) {
        return { centers: [] };
    }
    const D = perp(n);
    const u = dot(D, sub(P, A));
    const B = add(add(A, mul(D, u)), mul(n, r));
    const argument = r * r - (r - s) * (r - s);
    if (argument > 0) {
        const h = Math.sqrt(argument);
        return { centers: [sub(B, mul(D, h)), add(B, mul(D, h))] };
    }
    return { centers: [B] };
}

// Compare two unordered center sets with a scale-relative tolerance.
function expectSameCenters(actual: Vector[], expected: Vector[], tol: number): void {
    expect(actual.length).toBe(expected.length);
    const used = new Array<boolean>(expected.length).fill(false);
    for (const a of actual) {
        let matched = -1;
        for (let i = 0; i < expected.length; ++i) {
            if (!used[i] && length(sub(a, expected[i])) <= tol) { matched = i; break; }
        }
        expect(matched, `no match for center (${a.values}) in `
            + `[${expected.map(e => e.values.join(','))}]`).toBeGreaterThanOrEqual(0);
        used[matched] = true;
    }
}

// A configuration built in the (D, N) frame of the line so that the signed
// distance s of P is controlled directly.
const tangentConfiguration = fc.record({
    origin: wellScaledVector(2, -6, 6),
    angle: scaled(-Math.PI, Math.PI, 64),
    r: scaled(0.25, 5, 32),
    u: scaled(-6, 6, 64),
    // s = sign * factor * 2r with factor in [0.05, 0.95], so P is strictly
    // between the line and the far tangent line, on either side.
    factor: scaled(0.05, 0.95, 32),
    sign: fc.constantFrom(-1, 1)
}).map(({ origin, angle, r, u, factor, sign }) => {
    const N = v2(Math.cos(angle), Math.sin(angle));
    const D = perp(N);
    const s = sign * factor * 2 * r;
    const P = add(add(origin, mul(D, u)), mul(N, s));
    return { P, A: origin, N, r };
});

describe('circleThroughPointSpecifiedTangentAndRadius verification', () => {
    it('reports exactly the solutions of the constraint set', () => {
        check(tangentConfiguration, ({ P, A, N, r }) => {
            const { numCircles, circle } = circleThroughPointSpecifiedTangentAndRadius(P, A, N, r);
            const reference = referenceCenters(P, A, N, r);
            expect(numCircles).toBe(reference.centers.length);
            const scale = Math.max(r, length(sub(P, A)), 1);
            expectSameCenters(circle.slice(0, numCircles).map(c => c.center),
                reference.centers, 1e-9 * scale);
            for (let i = numCircles; i < 2; ++i) {
                expect(circle[i].radius).toBe(0);
                expect(circle[i].center.values).toEqual([0, 0]);
            }
        });
    });

    it('every reported circle passes through P and is tangent to the line', () => {
        check(tangentConfiguration, ({ P, A, N, r }) => {
            const { numCircles, circle } = circleThroughPointSpecifiedTangentAndRadius(P, A, N, r);
            expect(numCircles).toBe(2);   // 0 < s < 2r for this generator
            const scale = Math.max(r, length(sub(P, A)), 1);
            for (let i = 0; i < numCircles; ++i) {
                expect(circle[i].radius).toBe(r);
                expectClose(length(sub(circle[i].center, P)), r, 1e-9 * scale, 1e-12);
                expectClose(Math.abs(dot(N, sub(circle[i].center, A))), r,
                    1e-9 * scale, 1e-12);
                // The center is on the same side of the line as P.
                expect(dot(N, sub(circle[i].center, A)) * dot(N, sub(P, A)))
                    .toBeGreaterThan(0);
            }
            // The two centers are distinct and symmetric about the foot of P.
            const midpoint = mul(add(circle[0].center, circle[1].center), 0.5);
            expectClose(dot(perp(N), sub(midpoint, P)), 0, 1e-9 * scale, 1e-12);
        });
    });

    it('is insensitive to the sign of the normal and to rigid motions', () => {
        check(fc.tuple(tangentConfiguration, scaled(-Math.PI, Math.PI, 32),
            wellScaledVector(2, -4, 4)), ([{ P, A, N, r }, angle, t]) => {
            const flipped = circleThroughPointSpecifiedTangentAndRadius(
                P, A, negate(N), r);
            const direct = circleThroughPointSpecifiedTangentAndRadius(P, A, N, r);
            expect(flipped.numCircles).toBe(direct.numCircles);
            const scale = Math.max(r, length(sub(P, A)), 1);
            expectSameCenters(flipped.circle.slice(0, flipped.numCircles).map(c => c.center),
                direct.circle.slice(0, direct.numCircles).map(c => c.center),
                1e-9 * scale);

            const c = Math.cos(angle), s = Math.sin(angle);
            const rotate = (v: Vector): Vector => Vector.fromArray([
                c * v.values[0] - s * v.values[1], s * v.values[0] + c * v.values[1]]);
            const move = (v: Vector): Vector => add(rotate(v), t);
            const moved = circleThroughPointSpecifiedTangentAndRadius(
                move(P), move(A), rotate(N), r);
            expect(moved.numCircles).toBe(direct.numCircles);
            expectSameCenters(moved.circle.slice(0, moved.numCircles).map(c => c.center),
                direct.circle.slice(0, direct.numCircles).map(c => c.center)
                    .map(move), 1e-9 * Math.max(scale, 5));
        });
    });

    it('does not modify the caller normal', () => {
        check(tangentConfiguration, ({ P, A, N, r }) => {
            const before = [...N.values];
            circleThroughPointSpecifiedTangentAndRadius(P, A, N, r);
            expect(N.values).toEqual(before);
            // Upstream takes N by value, so its caller is likewise unaffected.
            circleThroughPointSpecifiedTangentAndRadius(
                sub(A, sub(P, A)), A, N, r);
            expect(N.values).toEqual(before);
        });
    });

    it('the exact branch values follow the case analysis', () => {
        // Axis-aligned configurations make s = Dot(N, P-A) exact.
        const A0 = v2(2, -1);
        const N0 = v2(0, 1);
        const r = 2;
        // Case (1): s = 0.
        const onLine = circleThroughPointSpecifiedTangentAndRadius(v2(5, -1), A0, N0, r);
        expect(onLine.numCircles).toBe(2);
        expect(onLine.circle[0].center.values).toEqual([5, -3]);
        expect(onLine.circle[1].center.values).toEqual([5, 1]);
        // Case (2): s = r; the centers are P -/+ r*Perp(N) = P -/+ (2,0)... .
        const atRadius = circleThroughPointSpecifiedTangentAndRadius(v2(5, 1), A0, N0, r);
        expect(atRadius.numCircles).toBe(2);
        expect(atRadius.circle[0].center.values).toEqual([3, 1]);
        expect(atRadius.circle[1].center.values).toEqual([7, 1]);
        // Case (3): s = 2*r, one circle centered at P - r*N.
        const atDiameter = circleThroughPointSpecifiedTangentAndRadius(v2(5, 3), A0, N0, r);
        expect(atDiameter.numCircles).toBe(1);
        expect(atDiameter.circle[0].center.values).toEqual([5, 1]);
        expect(atDiameter.circle[1].center.values).toEqual([0, 0]);
        expect(atDiameter.circle[1].radius).toBe(0);
        // Case (4): s > 2*r, no circles.
        const tooFar = circleThroughPointSpecifiedTangentAndRadius(v2(5, 4), A0, N0, r);
        expect(tooFar.numCircles).toBe(0);
        expect(tooFar.circle[0].center.values).toEqual([0, 0]);
        expect(tooFar.circle[0].radius).toBe(0);
        // The negative side mirrors the positive side.
        const below = circleThroughPointSpecifiedTangentAndRadius(v2(5, -5), A0, N0, r);
        expect(below.numCircles).toBe(1);
        expect(below.circle[0].center.values).toEqual([5, -3]);
    });
});
