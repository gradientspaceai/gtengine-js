import { describe, expect, it } from 'vitest';
import { MinimumAreaBox2 } from '../src/MinimumAreaBox2';
import type { OrientedBox2 } from '../src/OrientedBox';
import { Vector, dot, sub } from '../src/Vector';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// A deterministic pseudorandom generator so the tests are reproducible.
function makeRandom(seed: number): () => number {
    let s = seed;
    return () => {
        s = (1103515245 * s + 12345) % 2147483648;
        return s / 2147483648;
    };
}

// Andrew's monotone chain, an independent convex hull implementation for the
// brute-force cross-checks. The result is counterclockwise ordered.
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

// The minimum-area box over all hull-edge-aligned candidate boxes. This is
// the textbook O(n^2) statement of the problem: the minimum-area enclosing
// box of a convex polygon has a side flush with one of the polygon edges.
function bruteForceMinArea(points: readonly Vector[]):
    { area: number, extent0: number, extent1: number } {
    const hull = convexHull(points);
    let best = Number.MAX_VALUE;
    let bestE0 = 0;
    let bestE1 = 0;
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
        let umin = Number.MAX_VALUE;
        let umax = -Number.MAX_VALUE;
        let vmin = Number.MAX_VALUE;
        let vmax = -Number.MAX_VALUE;
        for (const p of points) {
            const u = ux * p.get(0) + uy * p.get(1);
            const v = -uy * p.get(0) + ux * p.get(1);
            umin = Math.min(umin, u);
            umax = Math.max(umax, u);
            vmin = Math.min(vmin, v);
            vmax = Math.max(vmax, v);
        }
        const area = (umax - umin) * (vmax - vmin);
        if (area < best) {
            best = area;
            bestE0 = 0.5 * (umax - umin);
            bestE1 = 0.5 * (vmax - vmin);
        }
    }
    return { area: best, extent0: bestE0, extent1: bestE1 };
}

// Every input point must lie inside the box.
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

// The axes must be orthonormal and right-handed.
function expectOrthonormalFrame(box: OrientedBox2): void {
    expect(dot(box.axis[0], box.axis[0])).toBeCloseTo(1, 12);
    expect(dot(box.axis[1], box.axis[1])).toBeCloseTo(1, 12);
    expect(dot(box.axis[0], box.axis[1])).toBeCloseTo(0, 12);
    const det = box.axis[0].get(0) * box.axis[1].get(1) -
        box.axis[0].get(1) * box.axis[1].get(0);
    expect(det).toBeCloseTo(1, 12);
}

// The unordered pair of extents, so a box can be compared regardless of
// which axis is the "bottom" direction.
function sortedExtents(box: OrientedBox2): [number, number] {
    const a = box.extent.get(0);
    const b = box.extent.get(1);
    return a <= b ? [a, b] : [b, a];
}

function rotate(points: readonly Vector[], angle: number, cx = 0,
    cy = 0): Vector[] {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return points.map(p => v2(
        cx + c * p.get(0) - s * p.get(1),
        cy + s * p.get(0) + c * p.get(1)));
}

describe('MinimumAreaBox2', () => {
    it('computes the box of an axis-aligned rectangle point set', () => {
        const points = [
            v2(0, 0), v2(4, 0), v2(4, 2), v2(0, 2),
            v2(1, 1), v2(2, 0.5), v2(3, 1.5)
        ];
        const query = new MinimumAreaBox2();
        const box = query.compute(points);

        expect(query.getArea()).toBeCloseTo(8, 12);
        expect(sortedExtents(box)).toEqual([
            expect.closeTo(1, 12), expect.closeTo(2, 12)
        ]);
        expect(box.center.get(0)).toBeCloseTo(2, 12);
        expect(box.center.get(1)).toBeCloseTo(1, 12);
        expectOrthonormalFrame(box);
        expectContainsAll(box, points);
        expect(4 * box.extent.get(0) * box.extent.get(1))
            .toBeCloseTo(query.getArea(), 12);
    });

    it('computes the box of a rotated rectangle point set', () => {
        const base = [
            v2(-3, -1), v2(3, -1), v2(3, 1), v2(-3, 1),
            v2(0, 0), v2(1.5, 0.5)
        ];
        for (const angle of [0.1, 0.7, 1.3, 2.2, 3.0, 4.5]) {
            const points = rotate(base, angle, 5, -2);
            const query = new MinimumAreaBox2();
            const box = query.compute(points);

            expect(query.getArea()).toBeCloseTo(12, 9);
            const [e0, e1] = sortedExtents(box);
            expect(e0).toBeCloseTo(1, 9);
            expect(e1).toBeCloseTo(3, 9);
            expect(box.center.get(0)).toBeCloseTo(5, 9);
            expect(box.center.get(1)).toBeCloseTo(-2, 9);
            expectOrthonormalFrame(box);
            expectContainsAll(box, points, 1e-9);
        }
    });

    it('computes the box of a triangle (area is twice the triangle area)', () => {
        // A well-known result: the minimum-area enclosing rectangle of a
        // triangle is flush with one of its edges and has exactly twice the
        // triangle's area.
        const triangles: Vector[][] = [
            [v2(0, 0), v2(4, 0), v2(1, 3)],
            [v2(0, 0), v2(5, 1), v2(2, 4)],
            [v2(-2, -1), v2(3, -2), v2(0, 5)],
            [v2(0, 0), v2(1, 0), v2(0.5, 10)]
        ];
        for (const t of triangles) {
            const triArea = 0.5 * Math.abs(
                (t[1].get(0) - t[0].get(0)) * (t[2].get(1) - t[0].get(1)) -
                (t[1].get(1) - t[0].get(1)) * (t[2].get(0) - t[0].get(0)));
            const query = new MinimumAreaBox2();
            const box = query.compute(t);
            expect(query.getArea()).toBeCloseTo(2 * triArea, 9);
            expectOrthonormalFrame(box);
            expectContainsAll(box, t);
        }
    });

    it('computes the box of a square rotated 45 degrees', () => {
        // The minimum-area box of a diamond is the diamond's own bounding
        // box in the rotated frame, of area 2 for unit "radius".
        const points = [v2(1, 0), v2(0, 1), v2(-1, 0), v2(0, -1)];
        const query = new MinimumAreaBox2();
        const box = query.compute(points);
        expect(query.getArea()).toBeCloseTo(2, 12);
        expect(sortedExtents(box)).toEqual([
            expect.closeTo(Math.SQRT1_2, 12), expect.closeTo(Math.SQRT1_2, 12)
        ]);
        expect(box.center.get(0)).toBeCloseTo(0, 12);
        expect(box.center.get(1)).toBeCloseTo(0, 12);
        expectContainsAll(box, points);
    });

    it('handles a 0-dimensional point set', () => {
        const points = [v2(3, -7), v2(3, -7), v2(3, -7)];
        const query = new MinimumAreaBox2();
        const box = query.compute(points);
        expect(box.center.values).toEqual([3, -7]);
        expect(box.extent.values).toEqual([0, 0]);
        expect(box.axis[0].values).toEqual([1, 0]);
        expect(box.axis[1].values).toEqual([0, 1]);
        expect([...query.getHull()]).toEqual([0]);
        expect(query.getArea()).toBe(0);
    });

    it('handles a 1-dimensional point set and reports the extreme indices', () => {
        // The hull origin is points[1], not points[0]: upstream leaves the
        // extreme indices at 0, which is wrong. The port reports the actual
        // extreme points.
        const points = [v2(2, 2), v2(0, 0), v2(1, 1), v2(4, 4), v2(3, 3)];
        const query = new MinimumAreaBox2();
        const box = query.compute(points);

        expect(box.center.get(0)).toBeCloseTo(2, 12);
        expect(box.center.get(1)).toBeCloseTo(2, 12);
        expect(box.extent.get(0)).toBeCloseTo(Math.sqrt(8), 12);
        expect(box.extent.get(1)).toBe(0);
        expectOrthonormalFrame(box);

        const hull = [...query.getHull()];
        expect(hull).toHaveLength(2);
        // The extremes are points[1] = (0,0) and points[3] = (4,4).
        expect(new Set(hull)).toEqual(new Set([1, 3]));
        expect(query.getArea()).toBe(0);
        expectContainsAll(box, points);
    });

    it('resets the area and support indices for a degenerate query', () => {
        const query = new MinimumAreaBox2();
        query.compute([v2(0, 0), v2(4, 0), v2(4, 2), v2(0, 2)]);
        expect(query.getArea()).toBeCloseTo(8, 12);

        // Upstream would leave the stale area of 8 in place here.
        query.compute([v2(9, 9), v2(9, 9)]);
        expect(query.getArea()).toBe(0);
        expect([...query.getSupportIndices()]).toEqual([0, 0, 0, 0]);
    });

    it('supports a caller-supplied convex polygon, with and without indices',
        () => {
            const polygon = [v2(0, 0), v2(4, 0), v2(4, 3), v2(0, 3)];
            const direct = new MinimumAreaBox2();
            const boxDirect = direct.computeConvexPolygon(polygon, []);
            expect(direct.getArea()).toBeCloseTo(12, 12);
            expect([...direct.getHull()]).toEqual([0, 1, 2, 3]);
            expectContainsAll(boxDirect, polygon);

            // The same polygon as a subset of a larger point array.
            const points = [
                v2(100, 100), v2(0, 0), v2(-50, 7), v2(4, 0), v2(4, 3), v2(0, 3)
            ];
            const indexed = new MinimumAreaBox2();
            const boxIndexed = indexed.computeConvexPolygon(points, [1, 3, 4, 5]);
            expect(indexed.getArea()).toBeCloseTo(12, 12);
            expectContainsAll(boxIndexed, polygon);
            expect(boxIndexed.center.values)
                .toEqual(boxDirect.center.values);
        });

    it('returns a degenerate box for invalid caller-supplied polygons', () => {
        const query = new MinimumAreaBox2();
        const box = query.computeConvexPolygon([v2(0, 0), v2(1, 1)], []);
        expect(box.center.values).toEqual([0, 0]);
        expect(box.extent.values).toEqual([0, 0]);
        expect(box.axis[0].values).toEqual([1, 0]);
        expect(box.axis[1].values).toEqual([0, 1]);

        const box2 = query.computeConvexPolygon(
            [v2(0, 0), v2(1, 0), v2(1, 1)], [0, 1]);
        expect(box2.extent.values).toEqual([0, 0]);
    });

    it('keeps the polygon corners when the input has duplicate points', () => {
        // The upstream collinear-point removal drops a genuine corner when
        // the preceding edge is the zero-length edge of a duplicate point.
        // The port compares against the most recent nonzero edge instead.
        const polygon = [
            v2(0, 0), v2(0, 0), v2(1, 0), v2(1, 1), v2(1, 1), v2(0, 1)
        ];
        const query = new MinimumAreaBox2();
        const box = query.computeConvexPolygon(polygon, []);
        expect(query.getArea()).toBeCloseTo(1, 12);
        expect(sortedExtents(box)).toEqual([
            expect.closeTo(0.5, 12), expect.closeTo(0.5, 12)
        ]);
        expectContainsAll(box, polygon);
    });

    it('drops interior collinear points of a caller-supplied polygon', () => {
        const polygon = [
            v2(0, 0), v2(2, 0), v2(4, 0), v2(4, 2), v2(2, 2), v2(0, 2)
        ];
        const query = new MinimumAreaBox2();
        const box = query.computeConvexPolygon(polygon, []);
        expect(query.getArea()).toBeCloseTo(8, 12);
        expectContainsAll(box, polygon);
    });

    it('matches a brute-force hull-edge search on random point sets', () => {
        const rnd = makeRandom(20240119);
        for (let trial = 0; trial < 60; ++trial) {
            const numPoints = 5 + Math.floor(rnd() * 40);
            const points: Vector[] = [];
            for (let i = 0; i < numPoints; ++i) {
                points.push(v2(20 * rnd() - 10, 20 * rnd() - 10));
            }

            const query = new MinimumAreaBox2();
            const box = query.compute(points);
            const brute = bruteForceMinArea(points);

            // The minimizing direction is compared through the area, not
            // through the extents: the area-versus-angle curve can be flat
            // enough near its minimum that two very different aspect ratios
            // tie to within round-off, and either is a correct answer.
            expect(query.getArea()).toBeCloseTo(brute.area, 8);
            expect(query.getArea()).toBeLessThanOrEqual(brute.area + 1e-8);
            const [e0, e1] = sortedExtents(box);
            expect(4 * e0 * e1).toBeCloseTo(query.getArea(), 8);
            expectOrthonormalFrame(box);
            expectContainsAll(box, points, 1e-8);
        }
    });

    it('gives the same result for the O(n) and O(n^2) searches', () => {
        const rnd = makeRandom(777001);
        for (let trial = 0; trial < 40; ++trial) {
            const numPoints = 6 + Math.floor(rnd() * 30);
            const points: Vector[] = [];
            for (let i = 0; i < numPoints; ++i) {
                // A mix of a disk and a rectangle so both regular and
                // irregular hulls are exercised.
                if (trial % 2 === 0) {
                    const a = 2 * Math.PI * rnd();
                    const r = 5 * Math.sqrt(rnd());
                    points.push(v2(r * Math.cos(a), r * Math.sin(a)));
                }
                else {
                    points.push(v2(10 * rnd(), 3 * rnd()));
                }
            }

            const fast = new MinimumAreaBox2();
            const boxFast = fast.compute(points, true);
            const slow = new MinimumAreaBox2();
            const boxSlow = slow.compute(points, false);

            expect(fast.getArea()).toBeCloseTo(slow.getArea(), 9);
            expect(sortedExtents(boxFast)[0])
                .toBeCloseTo(sortedExtents(boxSlow)[0], 9);
            expect(sortedExtents(boxFast)[1])
                .toBeCloseTo(sortedExtents(boxSlow)[1], 9);
            expect(boxFast.center.get(0)).toBeCloseTo(boxSlow.center.get(0), 9);
            expect(boxFast.center.get(1)).toBeCloseTo(boxSlow.center.get(1), 9);
        }
    });

    it('is invariant under rigid motions on random point sets', () => {
        const rnd = makeRandom(31337);
        for (let trial = 0; trial < 25; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 14; ++i) {
                points.push(v2(8 * rnd() - 4, 8 * rnd() - 4));
            }
            const angle = 2 * Math.PI * rnd();
            const moved = rotate(points, angle, 17, -23);

            const q0 = new MinimumAreaBox2();
            q0.compute(points);
            const q1 = new MinimumAreaBox2();
            q1.compute(moved);
            expect(q1.getArea()).toBeCloseTo(q0.getArea(), 8);
        }
    });

    it('reports support indices that index the hull vertices', () => {
        const rnd = makeRandom(4242);
        for (let trial = 0; trial < 20; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 20; ++i) {
                points.push(v2(6 * rnd() - 3, 6 * rnd() - 3));
            }
            const query = new MinimumAreaBox2();
            const box = query.compute(points);
            const hull = query.getHull();
            const support = query.getSupportIndices();

            // Each support index selects a hull vertex, and that vertex must
            // lie on the corresponding box side.
            for (let k = 0; k < 4; ++k) {
                expect(support[k]).toBeGreaterThanOrEqual(0);
                expect(support[k]).toBeLessThan(hull.length);
                const p = query.getPoints()[hull[support[k]]];
                const diff = sub(p, box.center);
                // Sides are bottom, right, top, left of the (axis0, axis1)
                // frame, i.e. axis1 = -extent1, axis0 = +extent0,
                // axis1 = +extent1, axis0 = -extent0.
                const projections = [
                    -dot(diff, box.axis[1]), dot(diff, box.axis[0]),
                    dot(diff, box.axis[1]), -dot(diff, box.axis[0])
                ];
                const extents = [
                    box.extent.get(1), box.extent.get(0),
                    box.extent.get(1), box.extent.get(0)
                ];
                expect(projections[k]).toBeCloseTo(extents[k], 8);
            }
        }
    });
});
