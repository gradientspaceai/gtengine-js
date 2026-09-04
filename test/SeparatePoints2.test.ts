import { describe, expect, it } from 'vitest';
import { SeparatePoints2 } from '../src/SeparatePoints2.js';
import type { Line2 } from '../src/Line.js';
import { Vector, dot } from '../src/Vector.js';
import { perp } from '../src/Vector2.js';
import { fc, check, latticeVector, expectClose } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V11): property-based cross-checks against an exact bigint
// separating-axis oracle. The generators are integer lattices, so every
// projection below is exact and the decision is not a tolerance.
// ---------------------------------------------------------------------------

type BigPoint2 = [bigint, bigint];

const toBig2 = (p: Vector): BigPoint2 =>
    [BigInt(p.values[0]), BigInt(p.values[1])];

/**
 * The convex hull of integer points as a counterclockwise vertex list,
 * computed by an independent monotone chain in bigint. Returns null when the
 * points are 0- or 1-dimensional, the cases SeparatePoints2 declines.
 */
function exactHullCCW(points: readonly Vector[]): BigPoint2[] | null {
    const pts = points.map(toBig2);
    pts.sort((p, q) => (p[0] < q[0] ? -1 : p[0] > q[0] ? 1
        : p[1] < q[1] ? -1 : p[1] > q[1] ? 1 : 0));
    const uniq: BigPoint2[] = [];
    for (const p of pts) {
        const last = uniq[uniq.length - 1];
        if (!last || last[0] !== p[0] || last[1] !== p[1]) { uniq.push(p); }
    }
    if (uniq.length < 3) { return null; }
    const cross = (o: BigPoint2, a: BigPoint2, b: BigPoint2): bigint =>
        (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    const build = (src: BigPoint2[]): BigPoint2[] => {
        const out: BigPoint2[] = [];
        for (const p of src) {
            while (out.length >= 2
                && cross(out[out.length - 2], out[out.length - 1], p) <= 0n) {
                out.pop();
            }
            out.push(p);
        }
        out.pop();
        return out;
    };
    const hull = build(uniq).concat(build(uniq.slice().reverse()));
    return hull.length >= 3 ? hull : null;
}

/**
 * The exact separating-axis decision, formulated by projection extremes
 * rather than by upstream's side counting: the sets are separated when some
 * directed hull edge <a,b> has the whole other hull in the closed halfplane
 * its outward normal Perp(b - a) points into, with at least one vertex
 * strictly outside. That is exactly the acceptance condition of upstream's
 * OnSameSide/WhichSide pair for a convex counterclockwise hull.
 */
function exactSeparated(hull0: readonly BigPoint2[],
    hull1: readonly BigPoint2[]): boolean {
    const testEdges = (owner: readonly BigPoint2[],
        other: readonly BigPoint2[]): boolean => {
        for (let j = 0; j < owner.length; ++j) {
            const a = owner[j];
            const b = owner[(j + 1) % owner.length];
            const nx = b[1] - a[1];
            const ny = a[0] - b[0];          // Perp(b - a) = (dy, -dx)
            let lo: bigint | null = null;
            let hi: bigint | null = null;
            for (const q of other) {
                const s = nx * (q[0] - a[0]) + ny * (q[1] - a[1]);
                if (lo === null || s < lo) { lo = s; }
                if (hi === null || s > hi) { hi = s; }
            }
            if (lo !== null && hi !== null && lo >= 0n && hi > 0n) {
                return true;
            }
        }
        return false;
    };
    return testEdges(hull0, hull1) || testEdges(hull1, hull0);
}

const latticeSet = (count: number, range: number, cx = 0, cy = 0):
    fc.Arbitrary<Vector[]> =>
    fc.array(latticeVector(2, -range, range), { minLength: count, maxLength: count })
        .map(vs => vs.map(v =>
            Vector.fromArray([v.values[0] + cx, v.values[1] + cy])));

describe('SeparatePoints2 verification', () => {
    it('agrees with an exact separating-axis oracle', () => {
        // The offsets make separated and overlapping draws roughly equally
        // likely.
        const pair = fc.oneof(
            fc.tuple(latticeSet(6, 4), latticeSet(6, 4, 5, 2)),
            fc.tuple(latticeSet(7, 5), latticeSet(7, 5, 9, 0)),
            fc.tuple(latticeSet(5, 3), latticeSet(8, 6, 2, 2)),
            fc.tuple(latticeSet(6, 4), latticeSet(6, 4)));
        check(pair, ([points0, points1]) => {
            const query = new SeparatePoints2();
            const result = query.compute(points0, points1);
            const hull0 = exactHullCCW(points0);
            const hull1 = exactHullCCW(points1);
            if (hull0 === null || hull1 === null) {
                // A 0- or 1-dimensional hull: upstream declines.
                expect(result.separated).toBe(false);
                return true;
            }
            expect(result.separated).toBe(exactSeparated(hull0, hull1));
            return true;
        }, 200);
    }, 30000);

    it('returns a line that really separates the point sets', () => {
        const pair = fc.oneof(
            fc.tuple(latticeSet(6, 4), latticeSet(6, 4, 6, 3)),
            fc.tuple(latticeSet(7, 5), latticeSet(7, 5, 11, 0)),
            fc.tuple(latticeSet(5, 3), latticeSet(8, 6, 4, 4)));
        let separations = 0;
        check(pair, ([points0, points1]) => {
            const query = new SeparatePoints2();
            const result = query.compute(points0, points1);
            if (!result.separated) { return true; }
            ++separations;

            const line = result.separatingLine;
            // The direction is unit length and the origin is one of the
            // input points of the owning set.
            expectClose(dot(line.direction, line.direction), 1, 1e-12, 1e-12);
            const normal = perp(line.direction);
            const signed = (p: Vector): number =>
                dot(normal, Vector.fromArray([
                    p.values[0] - line.origin.values[0],
                    p.values[1] - line.origin.values[1]]));

            // Every point of one set is on the closed nonpositive side and
            // every point of the other on the closed nonnegative side. The
            // coordinates are small integers and the direction is a rounded
            // unit vector, so the round-off in the signed distance is a few
            // ulps of the coordinate range.
            const tol = 1e-12;
            const max0 = Math.max(...points0.map(signed));
            const min0 = Math.min(...points0.map(signed));
            const max1 = Math.max(...points1.map(signed));
            const min1 = Math.min(...points1.map(signed));
            const zeroBelow = max0 <= tol && min1 >= -tol;
            const zeroAbove = min0 >= -tol && max1 <= tol;
            expect(zeroBelow || zeroAbove).toBe(true);
            // At least one point is strictly off the line on its side, which
            // is what OnSameSide's "+1" requires.
            expect(zeroBelow ? max1 > tol : min1 < -tol).toBe(true);
            return true;
        }, 150);
        expect(separations).toBeGreaterThan(20);
    }, 30000);

    it('is symmetric in its arguments', () => {
        const pair = fc.oneof(
            fc.tuple(latticeSet(6, 4), latticeSet(6, 4, 5, 2)),
            fc.tuple(latticeSet(7, 5), latticeSet(7, 5)));
        check(pair, ([points0, points1]) => {
            const query = new SeparatePoints2();
            expect(query.compute(points0, points1).separated)
                .toBe(query.compute(points1, points0).separated);
            return true;
        }, 150);
    }, 30000);

    it('is invariant under integer translation and lattice symmetry', () => {
        check(fc.tuple(latticeSet(6, 4), latticeSet(6, 4, 5, 2),
            latticeVector(2, -20, 20), fc.integer({ min: 0, max: 7 })),
            ([points0, points1, shift, g]) => {
                const query = new SeparatePoints2();
                const base = query.compute(points0, points1).separated;
                const transform = (p: Vector): Vector => {
                    const x = p.values[0] + shift.values[0];
                    const y = p.values[1] + shift.values[1];
                    const flipped = (g & 4) !== 0 ? [y, x] : [x, y];
                    switch (g & 3) {
                        case 0: return Vector.fromArray([flipped[0], flipped[1]]);
                        case 1: return Vector.fromArray([-flipped[1], flipped[0]]);
                        case 2: return Vector.fromArray([-flipped[0], -flipped[1]]);
                        default: return Vector.fromArray([flipped[1], -flipped[0]]);
                    }
                };
                expect(query.compute(points0.map(transform),
                    points1.map(transform)).separated).toBe(base);
                return true;
            }, 150);
    }, 30000);

    it('never separates a set from a superset that contains it', () => {
        // The hull of points0 is contained in the hull of points0 + extra,
        // so their interiors overlap and no separating line can exist.
        check(fc.tuple(latticeSet(6, 5), latticeSet(4, 5)),
            ([points0, extra]) => {
                const query = new SeparatePoints2();
                const both = points0.concat(extra);
                const hull0 = exactHullCCW(points0);
                if (hull0 === null) { return true; }
                expect(query.compute(points0, both).separated).toBe(false);
                expect(query.compute(both, points0).separated).toBe(false);
                return true;
            }, 120);
    }, 30000);

    it('separates half-plane-disjoint sets by construction', () => {
        // Every point of set 0 has x <= 0 and every point of set 1 has
        // x >= 1, so the sets are strictly separated by the y-axis.
        check(fc.tuple(latticeSet(6, 5), latticeSet(6, 5)),
            ([a, b]) => {
                const points0 = a.map(p => Vector.fromArray(
                    [-Math.abs(p.values[0]), p.values[1]]));
                const points1 = b.map(p => Vector.fromArray(
                    [Math.abs(p.values[0]) + 1, p.values[1]]));
                const hull0 = exactHullCCW(points0);
                const hull1 = exactHullCCW(points1);
                if (hull0 === null || hull1 === null) { return true; }
                expect(new SeparatePoints2().compute(points0, points1)
                    .separated).toBe(true);
                return true;
            }, 120);
    }, 30000);
});
