import { describe, expect, it } from 'vitest';
import { MinimumWidthPoints2 } from '../src/MinimumWidthPoints2.js';
import type { OrientedBox2 } from '../src/OrientedBox.js';
import { Vector, dot, sub } from '../src/Vector.js';
import { fc, check, latticeVector, expectClose } from './helpers/arbitraries.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function makeRandom(seed: number): () => number {
    let s = seed;
    return () => {
        s = (1103515245 * s + 12345) % 2147483648;
        return s / 2147483648;
    };
}

// Andrew's monotone chain, an independent convex hull implementation.
function convexHull(points: readonly Vector[]): Vector[] {
    const pts = [...points].sort((a, b) =>
        a.get(0) !== b.get(0) ? a.get(0) - b.get(0) : a.get(1) - b.get(1));
    const unique: Vector[] = [];
    for (const p of pts) {
        const last = unique[unique.length - 1];
        if (last === undefined || last.get(0) !== p.get(0) ||
            last.get(1) !== p.get(1)) {
            unique.push(p);
        }
    }
    if (unique.length < 3) {
        return unique;
    }
    const cross = (o: Vector, a: Vector, b: Vector): number =>
        (a.get(0) - o.get(0)) * (b.get(1) - o.get(1)) -
        (a.get(1) - o.get(1)) * (b.get(0) - o.get(0));
    const build = (src: readonly Vector[]): Vector[] => {
        const chain: Vector[] = [];
        for (const p of src) {
            while (chain.length >= 2 &&
                cross(chain[chain.length - 2], chain[chain.length - 1], p) <= 0) {
                chain.pop();
            }
            chain.push(p);
        }
        chain.pop();
        return chain;
    };
    return [...build(unique), ...build([...unique].reverse())];
}

// The width of the point set: the minimum over the hull edge directions of
// the maximum perpendicular distance from the edge line.
function bruteForceWidth(points: readonly Vector[]): number {
    const hull = convexHull(points);
    let minWidth = Number.MAX_VALUE;
    const n = hull.length;
    for (let i0 = n - 1, i1 = 0; i1 < n; i0 = i1++) {
        const dx = hull[i1].get(0) - hull[i0].get(0);
        const dy = hull[i1].get(1) - hull[i0].get(1);
        const len = Math.hypot(dx, dy);
        if (len === 0) {
            continue;
        }
        const ux = dx / len;
        const uy = dy / len;
        let maxWidth = 0;
        for (const p of points) {
            const w = Math.abs(ux * (p.get(1) - hull[i0].get(1)) -
                uy * (p.get(0) - hull[i0].get(0)));
            maxWidth = Math.max(maxWidth, w);
        }
        minWidth = Math.min(minWidth, maxWidth);
    }
    return minWidth;
}

function expectContainsAll(box: OrientedBox2, points: readonly Vector[],
    tol = 1e-9): void {
    for (const p of points) {
        const diff = sub(p, box.center);
        for (let i = 0; i < 2; ++i) {
            expect(Math.abs(dot(diff, box.axis[i])))
                .toBeLessThanOrEqual(box.extent.get(i) + tol);
        }
    }
}

// The axes must be orthonormal, and axis[1] = -Perp(axis[0]) as documented.
function expectOrthonormalFrame(box: OrientedBox2): void {
    expect(dot(box.axis[0], box.axis[0])).toBeCloseTo(1, 12);
    expect(dot(box.axis[1], box.axis[1])).toBeCloseTo(1, 12);
    expect(dot(box.axis[0], box.axis[1])).toBeCloseTo(0, 12);
    expect(box.axis[1].get(0)).toBeCloseTo(box.axis[0].get(1), 12);
    expect(box.axis[1].get(1)).toBeCloseTo(-box.axis[0].get(0), 12);
}

function rotate(points: readonly Vector[], angle: number, cx = 0,
    cy = 0): Vector[] {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return points.map(p => v2(
        cx + c * p.get(0) - s * p.get(1),
        cy + s * p.get(0) + c * p.get(1)));
}

describe('MinimumWidthPoints2', () => {
    it('computes the width of an axis-aligned rectangle point set', () => {
        const points = [
            v2(0, 0), v2(4, 0), v2(4, 1), v2(0, 1), v2(2, 0.5), v2(3, 0.25)
        ];
        for (const useCalipers of [true, false]) {
            const box = new MinimumWidthPoints2().compute(points, useCalipers);
            // The width is the short side, 1, so extent[0] = 0.5.
            expect(box.extent.get(0)).toBeCloseTo(0.5, 12);
            expect(box.extent.get(1)).toBeCloseTo(2, 12);
            expect(box.center.get(0)).toBeCloseTo(2, 12);
            expect(box.center.get(1)).toBeCloseTo(0.5, 12);
            expectOrthonormalFrame(box);
            expectContainsAll(box, points);
        }
    });

    it('computes the width of a rotated rectangle point set', () => {
        const base = [v2(-3, -0.5), v2(3, -0.5), v2(3, 0.5), v2(-3, 0.5),
            v2(0, 0)];
        for (const angle of [0.2, 0.9, 1.8, 2.7, 3.9, 5.5]) {
            const points = rotate(base, angle, -4, 6);
            for (const useCalipers of [true, false]) {
                const box = new MinimumWidthPoints2()
                    .compute(points, useCalipers);
                expect(2 * box.extent.get(0)).toBeCloseTo(1, 9);
                expect(2 * box.extent.get(1)).toBeCloseTo(6, 9);
                expect(box.center.get(0)).toBeCloseTo(-4, 9);
                expect(box.center.get(1)).toBeCloseTo(6, 9);
                expectOrthonormalFrame(box);
                expectContainsAll(box, points, 1e-9);
            }
        }
    });

    it('computes the width of an equilateral triangle', () => {
        // The width of an equilateral triangle of side s is its altitude,
        // s*sqrt(3)/2.
        const s = 2;
        const h = (s * Math.sqrt(3)) / 2;
        const points = [v2(0, 0), v2(s, 0), v2(s / 2, h)];
        for (const useCalipers of [true, false]) {
            const box = new MinimumWidthPoints2().compute(points, useCalipers);
            expect(2 * box.extent.get(0)).toBeCloseTo(h, 12);
            expectOrthonormalFrame(box);
            expectContainsAll(box, points);
        }
    });

    it('computes the width of a right triangle', () => {
        // The width of the triangle with legs 3 and 4 is the altitude to the
        // hypotenuse, 3*4/5 = 2.4.
        const points = [v2(0, 0), v2(4, 0), v2(0, 3)];
        for (const useCalipers of [true, false]) {
            const box = new MinimumWidthPoints2().compute(points, useCalipers);
            expect(2 * box.extent.get(0)).toBeCloseTo(2.4, 12);
            expectContainsAll(box, points);
        }
    });

    it('computes the width of a regular polygon sampled on a circle', () => {
        // For a regular n-gon of circumradius R the width is 2*R*cos(pi/(2n))
        // when n is odd and 2*R*cos(pi/n) when n is even.
        const R = 3;
        for (const n of [5, 6, 7, 8, 9, 12]) {
            const points: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                const a = (2 * Math.PI * i) / n;
                points.push(v2(R * Math.cos(a), R * Math.sin(a)));
            }
            const expected = n % 2 === 1
                ? R * (1 + Math.cos(Math.PI / n))
                : 2 * R * Math.cos(Math.PI / n);
            for (const useCalipers of [true, false]) {
                const box = new MinimumWidthPoints2()
                    .compute(points, useCalipers);
                expect(2 * box.extent.get(0)).toBeCloseTo(expected, 9);
                expectContainsAll(box, points, 1e-9);
            }
        }
    });

    it('handles a 0-dimensional point set', () => {
        const points = [v2(-2, 5), v2(-2, 5), v2(-2, 5), v2(-2, 5)];
        const box = new MinimumWidthPoints2().compute(points);
        expect(box.center.values).toEqual([-2, 5]);
        expect(box.extent.values).toEqual([0, 0]);
        expect(box.axis[0].values).toEqual([1, 0]);
        expect(box.axis[1].values).toEqual([0, 1]);
    });

    it('handles a 1-dimensional point set', () => {
        const points = [v2(1, 0), v2(3, 0), v2(-1, 0), v2(0, 0)];
        const box = new MinimumWidthPoints2().compute(points);
        expect(box.extent.get(0)).toBe(0);
        expect(box.extent.get(1)).toBeCloseTo(2, 12);
        expect(box.center.get(0)).toBeCloseTo(1, 12);
        expect(box.center.get(1)).toBeCloseTo(0, 12);
        expect(dot(box.axis[0], box.axis[1])).toBeCloseTo(0, 12);
        expectContainsAll(box, points);
    });

    it('rejects inputs with fewer than three points', () => {
        const query = new MinimumWidthPoints2();
        expect(() => query.compute([v2(0, 0), v2(1, 1)])).toThrow(/Invalid input/);
        expect(() => query.computeIndexed(
            [v2(0, 0), v2(1, 0), v2(0, 1)], [0, 1])).toThrow(/Invalid input/);
    });

    it('supports an indexed subset of the input points', () => {
        const polygon = [v2(0, 0), v2(6, 0), v2(6, 2), v2(0, 2)];
        const points = [
            v2(-100, -100), v2(0, 0), v2(6, 0), v2(50, 50), v2(6, 2), v2(0, 2)
        ];
        const boxDirect = new MinimumWidthPoints2().compute(polygon);
        const boxIndexed = new MinimumWidthPoints2()
            .computeIndexed(points, [1, 2, 4, 5]);
        expect(boxIndexed.extent.values).toEqual(boxDirect.extent.values);
        expect(boxIndexed.center.values).toEqual(boxDirect.center.values);

        // An empty index array uses the points directly.
        const boxEmpty = new MinimumWidthPoints2().computeIndexed(polygon, []);
        expect(boxEmpty.extent.values).toEqual(boxDirect.extent.values);
    });

    it('matches a brute-force hull-edge search on random point sets', () => {
        const rnd = makeRandom(20250119);
        for (let trial = 0; trial < 60; ++trial) {
            const numPoints = 3 + Math.floor(rnd() * 40);
            const points: Vector[] = [];
            for (let i = 0; i < numPoints; ++i) {
                points.push(v2(20 * rnd() - 10, 20 * rnd() - 10));
            }

            const box = new MinimumWidthPoints2().compute(points, true);
            const expected = bruteForceWidth(points);
            expect(2 * box.extent.get(0)).toBeCloseTo(expected, 8);
            expectOrthonormalFrame(box);
            expectContainsAll(box, points, 1e-8);
        }
    });

    it('gives the same result for the calipers and brute-force searches', () => {
        const rnd = makeRandom(555321);
        for (let trial = 0; trial < 40; ++trial) {
            const numPoints = 6 + Math.floor(rnd() * 30);
            const points: Vector[] = [];
            for (let i = 0; i < numPoints; ++i) {
                if (trial % 2 === 0) {
                    const a = 2 * Math.PI * rnd();
                    const r = 4 * Math.sqrt(rnd());
                    points.push(v2(r * Math.cos(a), r * Math.sin(a)));
                }
                else {
                    points.push(v2(12 * rnd(), 2 * rnd()));
                }
            }

            const fast = new MinimumWidthPoints2().compute(points, true);
            const slow = new MinimumWidthPoints2().compute(points, false);
            expect(fast.extent.get(0)).toBeCloseTo(slow.extent.get(0), 9);
            expect(fast.extent.get(1)).toBeCloseTo(slow.extent.get(1), 9);
            expect(fast.center.get(0)).toBeCloseTo(slow.center.get(0), 9);
            expect(fast.center.get(1)).toBeCloseTo(slow.center.get(1), 9);
        }
    });

    it('is invariant under rigid motions on random point sets', () => {
        const rnd = makeRandom(9090);
        for (let trial = 0; trial < 25; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 15; ++i) {
                points.push(v2(8 * rnd() - 4, 8 * rnd() - 4));
            }
            const moved = rotate(points, 2 * Math.PI * rnd(), -13, 29);
            const w0 = new MinimumWidthPoints2().compute(points)
                .extent.get(0);
            const w1 = new MinimumWidthPoints2().compute(moved).extent.get(0);
            expect(w1).toBeCloseTo(w0, 8);
        }
    });

    it('is never wider than the minimum-area box is thin', () => {
        // The minimum width is a lower bound for the thinner extent of any
        // bounding box, in particular of the axis-aligned bounding box.
        const rnd = makeRandom(112233);
        for (let trial = 0; trial < 30; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 12; ++i) {
                points.push(v2(10 * rnd(), 4 * rnd()));
            }
            const width = 2 * new MinimumWidthPoints2().compute(points)
                .extent.get(0);
            const ys = points.map(p => p.get(1));
            const xs = points.map(p => p.get(0));
            const aabbThin = Math.min(
                Math.max(...xs) - Math.min(...xs),
                Math.max(...ys) - Math.min(...ys));
            expect(width).toBeLessThanOrEqual(aabbThin + 1e-9);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V11): property-based cross-checks against exact bigint
// oracles. The generators are integer lattices, so the candidate widths below
// are exact rationals (under one square root) rather than tolerances.
// ---------------------------------------------------------------------------

type ExactRatio = { num: bigint; den: bigint };

function ratioLess(a: ExactRatio, b: ExactRatio): boolean {
    return a.num * b.den < b.num * a.den;
}

/**
 * The exact minimum width of an integer point set: the minimum over every
 * direction spanned by a pair of distinct points of the spread of the points
 * perpendicular to that direction. The minimum-width strip is flush with a
 * convex-hull edge and every hull edge is one of those pairs, so the brute
 * force is the true minimum; every other pair yields a valid bounding strip,
 * which can only be wider. The returned ratio is the squared width.
 */
function bruteForceSqrWidthExact(points: readonly Vector[]): ExactRatio | null {
    const pts = points.map(p =>
        [BigInt(p.values[0]), BigInt(p.values[1])] as [bigint, bigint]);
    let best: ExactRatio | null = null;
    for (let i = 0; i < pts.length; ++i) {
        for (let j = i + 1; j < pts.length; ++j) {
            const dx = pts[j][0] - pts[i][0];
            const dy = pts[j][1] - pts[i][1];
            if (dx === 0n && dy === 0n) { continue; }
            let lo = -dy * pts[0][0] + dx * pts[0][1];
            let hi = lo;
            for (const p of pts) {
                const w = -dy * p[0] + dx * p[1];
                if (w < lo) { lo = w; }
                if (w > hi) { hi = w; }
            }
            const spread = hi - lo;
            const candidate: ExactRatio =
                { num: spread * spread, den: dx * dx + dy * dy };
            if (best === null || ratioLess(candidate, best)) { best = candidate; }
        }
    }
    return best;
}

function expectBoxEncloses(box: OrientedBox2, points: readonly Vector[],
    tolerance: number): void {
    for (const p of points) {
        const d = sub(p, box.center);
        for (let i = 0; i < 2; ++i) {
            expect(Math.abs(dot(d, box.axis[i])))
                .toBeLessThanOrEqual(box.extent.get(i) + tolerance);
        }
    }
}

const latticeCloud = (count: number, range: number): fc.Arbitrary<Vector[]> =>
    fc.array(latticeVector(2, -range, range),
        { minLength: count, maxLength: count });

describe('MinimumWidthPoints2 verification', () => {
    it('attains the exact minimum width over all strip directions', () => {
        check(fc.oneof(latticeCloud(6, 4), latticeCloud(9, 7),
            latticeCloud(13, 9)), points => {
            const query = new MinimumWidthPoints2();
            const box = query.compute(points);
            const expected = bruteForceSqrWidthExact(points);
            if (expected === null) {
                // All points coincide: the 0-dimensional path.
                expect(box.extent.get(0)).toBe(0);
                expect(box.extent.get(1)).toBe(0);
                return true;
            }
            const expectedWidth =
                Math.sqrt(Number(expected.num) / Number(expected.den));
            // The width is the exact rational squared width rounded once to
            // double and then square-rooted, so only a few ulps separate the
            // two evaluations.
            expectClose(2 * box.extent.get(0), expectedWidth, 1e-11, 1e-11);
            expectBoxEncloses(box, points, 1e-9);

            // The frame is orthonormal and right-handed:
            // axis[0] = -Perp(axis[1]).
            expectClose(dot(box.axis[0], box.axis[1]), 0, 1e-12, 1e-12);
            expectClose(dot(box.axis[0], box.axis[0]), 1, 1e-12, 1e-12);
            expectClose(dot(box.axis[1], box.axis[1]), 1, 1e-12, 1e-12);
            if (box.extent.get(0) > 0) {
                expectClose(box.axis[0].values[0], -box.axis[1].values[1],
                    1e-12, 1e-12);
                expectClose(box.axis[0].values[1], box.axis[1].values[0],
                    1e-12, 1e-12);
            }
            return true;
        }, 120);
    }, 30000);

    it('reports the tightest height along the width direction', () => {
        // 2*extent[1] must be the exact spread of the points along axis[1],
        // which is the length of the minimum-width strip's supporting edge
        // span.
        check(fc.oneof(latticeCloud(7, 5), latticeCloud(11, 8)), points => {
            const query = new MinimumWidthPoints2();
            const box = query.compute(points);
            if (box.extent.get(0) === 0) { return true; }
            let lo = Infinity, hi = -Infinity;
            for (const p of points) {
                const h = dot(sub(p, box.center), box.axis[1]);
                if (h < lo) { lo = h; }
                if (h > hi) { hi = h; }
            }
            expectClose(hi - lo, 2 * box.extent.get(1), 1e-9, 1e-9);
            expectClose(hi + lo, 0, 1e-9, 1e-9);
            return true;
        }, 100);
    }, 30000);

    it('agrees between the calipers and brute-force searches', () => {
        check(fc.oneof(latticeCloud(7, 5), latticeCloud(12, 9)), points => {
            const query = new MinimumWidthPoints2();
            const fast = query.compute(points, true);
            const slow = query.compute(points, false);
            // The width must agree; the height need not, because the two
            // searches can settle on different hull edges when the minimum
            // width is attained more than once.
            expectClose(fast.extent.get(0), slow.extent.get(0), 1e-11, 1e-11);
            expectBoxEncloses(fast, points, 1e-9);
            expectBoxEncloses(slow, points, 1e-9);
            return true;
        }, 100);
    }, 30000);

    it('is invariant under the lattice symmetry group', () => {
        check(fc.tuple(latticeCloud(9, 6), fc.integer({ min: 0, max: 7 })),
            ([points, g]) => {
                const query = new MinimumWidthPoints2();
                const base = query.compute(points);
                const transform = (p: Vector): Vector => {
                    const x = p.values[0], y = p.values[1];
                    const flipped = (g & 4) !== 0 ? [y, x] : [x, y];
                    switch (g & 3) {
                        case 0: return Vector.fromArray([flipped[0], flipped[1]]);
                        case 1: return Vector.fromArray([-flipped[1], flipped[0]]);
                        case 2: return Vector.fromArray([-flipped[0], -flipped[1]]);
                        default: return Vector.fromArray([flipped[1], -flipped[0]]);
                    }
                };
                const moved = query.compute(points.map(transform));
                // Only the width is invariant. When several hull edges
                // attain the minimum width, the search picks one of them and
                // the reported height (extent[1]) is that strip's length,
                // which differs between the tied choices.
                expectClose(moved.extent.get(0), base.extent.get(0),
                    1e-12, 1e-12);
                expectBoxEncloses(moved, points.map(transform), 1e-9);
                return true;
            }, 120);
    }, 30000);

    it('handles collinear input with the documented degenerate frame', () => {
        check(fc.tuple(latticeVector(2, -8, 8), latticeVector(2, -5, 5),
            fc.array(fc.integer({ min: -6, max: 6 }),
                { minLength: 3, maxLength: 8 })),
            ([origin, direction, ts]) => {
                if (direction.values[0] === 0 && direction.values[1] === 0) {
                    return true;
                }
                const points = ts.map(t => Vector.fromArray([
                    origin.values[0] + t * direction.values[0],
                    origin.values[1] + t * direction.values[1]]));
                const query = new MinimumWidthPoints2();
                const box = query.compute(points);
                expect(box.extent.get(0)).toBe(0);
                expectBoxEncloses(box, points, 1e-9);
                // The degenerate branch uses Perp(direction) where the
                // general branch uses -Perp(U); the handedness therefore
                // differs, an upstream quirk the port preserves. It is
                // harmless because extent[0] is zero.
                expectClose(dot(box.axis[0], box.axis[1]), 0, 1e-12, 1e-12);
                return true;
            }, 150);
    }, 30000);

    it('computeIndexed on a subset matches compute on that subset', () => {
        check(fc.tuple(latticeCloud(12, 7), fc.array(fc.nat(11),
            { minLength: 3, maxLength: 8 })), ([points, picks]) => {
            const indices = picks.map(i => i % points.length);
            const query = new MinimumWidthPoints2();
            const indexed = query.computeIndexed(points, indices);
            const direct = query.compute(indices.map(i => points[i]));
            // computeIndexed forwards to compute on the compacted points, so
            // both extents must agree exactly here.
            expect(indexed.extent.values).toEqual(direct.extent.values);
            expect(indexed.center.values).toEqual(direct.center.values);
            return true;
        }, 100);
    }, 30000);
});
