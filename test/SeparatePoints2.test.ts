import { describe, expect, it } from 'vitest';
import { SeparatePoints2 } from '../src/SeparatePoints2';
import type { Line2 } from '../src/Line';
import { Vector, dot } from '../src/Vector';
import { perp } from '../src/Vector2';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// The signed distance of 'p' from the line, using the perpendicular of the
// (unit-length) line direction as the normal.
function signedDistance(line: Line2, p: Vector): number {
    const n = perp(line.direction);
    return dot(n, p) - dot(n, line.origin);
}

// Verify that 'line' separates the two point sets: one set has all signed
// distances <= 0 and the other has all signed distances >= 0, and at least
// one point of each set is strictly off the line on its own side.
function expectSeparates(line: Line2, points0: readonly Vector[],
    points1: readonly Vector[]): void {
    const eps = 1e-12;
    const d0 = points0.map(p => signedDistance(line, p));
    const d1 = points1.map(p => signedDistance(line, p));

    const zeroSide0Neg = d0.every(d => d <= eps);
    const zeroSide0Pos = d0.every(d => d >= -eps);
    expect(zeroSide0Neg || zeroSide0Pos).toBe(true);

    if (zeroSide0Neg && !zeroSide0Pos) {
        // Set 0 is on the negative side, so set 1 must be on the positive.
        expect(d1.every(d => d >= -eps)).toBe(true);
        expect(d1.some(d => d > eps)).toBe(true);
    }
    else if (zeroSide0Pos && !zeroSide0Neg) {
        expect(d1.every(d => d <= eps)).toBe(true);
        expect(d1.some(d => d < -eps)).toBe(true);
    }
    else {
        // Set 0 lies entirely on the line, which cannot happen for a
        // 2-dimensional hull.
        expect(false).toBe(true);
    }

    // The direction must be unit length.
    expect(dot(line.direction, line.direction)).toBeCloseTo(1, 12);
}

const unitSquare: Vector[] = [
    v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1), v2(0.5, 0.5)
];

function translate(points: readonly Vector[], dx: number, dy: number): Vector[] {
    return points.map(p => v2(p.get(0) + dx, p.get(1) + dy));
}

describe('SeparatePoints2', () => {
    it('separates two disjoint axis-aligned squares', () => {
        const query = new SeparatePoints2();
        const points1 = translate(unitSquare, 3, 0);
        const result = query.compute(unitSquare, points1);
        expect(result.separated).toBe(true);
        expectSeparates(result.separatingLine, unitSquare, points1);
    });

    it('separates squares that are disjoint along several directions', () => {
        const query = new SeparatePoints2();
        const offsets: Array<[number, number]> = [
            [3, 0], [-3, 0], [0, 3], [0, -3], [2.5, 2.5], [-2.5, 3.5]
        ];
        for (const [dx, dy] of offsets) {
            const points1 = translate(unitSquare, dx, dy);
            const result = query.compute(unitSquare, points1);
            expect(result.separated).toBe(true);
            expectSeparates(result.separatingLine, unitSquare, points1);
        }
    });

    it('separates squares that touch along an edge', () => {
        // Touching is allowed: the separating line contains the shared edge,
        // and the sets have no interior overlap.
        const query = new SeparatePoints2();
        const points1 = translate(unitSquare, 1, 0);
        const result = query.compute(unitSquare, points1);
        expect(result.separated).toBe(true);
        expectSeparates(result.separatingLine, unitSquare, points1);
    });

    it('reports non-separable for overlapping squares', () => {
        const query = new SeparatePoints2();
        expect(query.compute(unitSquare, translate(unitSquare, 0.5, 0.5))
            .separated).toBe(false);
        expect(query.compute(unitSquare, translate(unitSquare, 0.5, 0))
            .separated).toBe(false);
        expect(query.compute(unitSquare, unitSquare).separated).toBe(false);
    });

    it('reports non-separable when one set contains the other', () => {
        const query = new SeparatePoints2();
        const big = [v2(-5, -5), v2(5, -5), v2(5, 5), v2(-5, 5)];
        expect(query.compute(big, unitSquare).separated).toBe(false);
        expect(query.compute(unitSquare, big).separated).toBe(false);
    });

    it('reports non-separable when a hull is 0- or 1-dimensional', () => {
        const query = new SeparatePoints2();
        const singlePoint = [v2(7, 7), v2(7, 7), v2(7, 7)];
        const collinear = [v2(5, 5), v2(6, 6), v2(7, 7), v2(8, 8)];

        // A degenerate hull cannot support the separating-axis test, so the
        // query reports 'not separated' even though the sets are disjoint.
        expect(query.compute(unitSquare, singlePoint).separated).toBe(false);
        expect(query.compute(singlePoint, unitSquare).separated).toBe(false);
        expect(query.compute(unitSquare, collinear).separated).toBe(false);
        expect(query.compute(collinear, unitSquare).separated).toBe(false);
    });

    it('separates rotated hexagons that are disjoint', () => {
        const query = new SeparatePoints2();
        const hexagon = (cx: number, cy: number, r: number, phase: number) => {
            const pts: Vector[] = [];
            for (let i = 0; i < 6; ++i) {
                const a = phase + (2 * Math.PI * i) / 6;
                pts.push(v2(cx + r * Math.cos(a), cy + r * Math.sin(a)));
            }
            return pts;
        };
        const h0 = hexagon(0, 0, 1, 0.3);
        const h1 = hexagon(4, 1, 1.5, 1.1);
        const result = query.compute(h0, h1);
        expect(result.separated).toBe(true);
        expectSeparates(result.separatingLine, h0, h1);

        // Overlapping hexagons are not separable.
        const h2 = hexagon(1, 0.5, 1.5, 1.1);
        expect(query.compute(h0, h2).separated).toBe(false);
    });

    it('agrees with a brute-force separating-axis test on random sets', () => {
        // Independent check: two convex sets are disjoint if and only if
        // some edge normal of one of the hulls separates them. The brute
        // force uses all edge directions of both point sets' hulls, computed
        // by projecting all raw points (projections of the hull bound the
        // projections of all points).
        let seed = 987654321;
        const rnd = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const makeSet = (cx: number, cy: number, n: number): Vector[] => {
            const pts: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                pts.push(v2(cx + 2 * rnd() - 1, cy + 2 * rnd() - 1));
            }
            return pts;
        };
        const bruteDisjoint = (a: readonly Vector[],
            b: readonly Vector[]): boolean => {
            const all = [...a, ...b];
            for (let i = 0; i < all.length; ++i) {
                for (let j = 0; j < all.length; ++j) {
                    if (i === j) {
                        continue;
                    }
                    const nx = -(all[j].get(1) - all[i].get(1));
                    const ny = all[j].get(0) - all[i].get(0);
                    const len = Math.hypot(nx, ny);
                    if (len < 1e-12) {
                        continue;
                    }
                    const proj = (p: Vector) =>
                        (nx * p.get(0) + ny * p.get(1)) / len;
                    const a0 = Math.min(...a.map(proj));
                    const a1 = Math.max(...a.map(proj));
                    const b0 = Math.min(...b.map(proj));
                    const b1 = Math.max(...b.map(proj));
                    if (a1 < b0 - 1e-9 || b1 < a0 - 1e-9) {
                        return true;
                    }
                }
            }
            return false;
        };

        const query = new SeparatePoints2();
        let numSeparated = 0;
        let numOverlapping = 0;
        for (let trial = 0; trial < 40; ++trial) {
            const a = makeSet(0, 0, 8);
            const b = makeSet(2.5 * rnd(), 2.5 * rnd(), 8);
            const result = query.compute(a, b);
            expect(result.separated).toBe(bruteDisjoint(a, b));
            if (result.separated) {
                ++numSeparated;
                expectSeparates(result.separatingLine, a, b);
            }
            else {
                ++numOverlapping;
            }
        }
        // Both outcomes must actually be exercised by the trials.
        expect(numSeparated).toBeGreaterThan(0);
        expect(numOverlapping).toBeGreaterThan(0);
    });
});
