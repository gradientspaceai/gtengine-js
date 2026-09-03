import { describe, expect, it } from 'vitest';
import { MinimumWidthPoints2 } from '../src/MinimumWidthPoints2.js';
import type { OrientedBox2 } from '../src/OrientedBox.js';
import { Vector, dot, sub } from '../src/Vector.js';

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
