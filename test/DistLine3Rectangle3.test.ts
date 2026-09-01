import { describe, expect, it } from 'vitest';
import { DistLine3Rectangle3 } from '../src/DistLine3Rectangle3';
import { Line } from '../src/Line';
import { Rectangle } from '../src/Rectangle';
import { Vector, add, dot, mul, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistLine3Rectangle3', () => {
    const query = new DistLine3Rectangle3();
    // The unit-ish rectangle in the z = 0 plane, center at the origin, with
    // extents 2 and 1.
    const rect = Rectangle.fromCenterAxisExtent(v(0, 0, 0),
        [v(1, 0, 0), v(0, 1, 0)], v(2, 1));

    it('reports zero distance for a line piercing the rectangle', () => {
        const result = query.compute(line([0.5, -0.25, -3], [0, 0, 1]), rect);
        expect(result.distance).toBe(0);
        expect(result.parameter).toBeCloseTo(3, 12);
        expect(result.cartesian[0]).toBeCloseTo(0.5, 12);
        expect(result.cartesian[1]).toBeCloseTo(-0.25, 12);
    });

    it('measures a line parallel to the plane of the rectangle', () => {
        const result = query.compute(line([0, 0, 5], [1, 0, 0]), rect);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(0, 12);
    });

    it('measures a line piercing the plane outside the rectangle', () => {
        // The line pierces the plane at (5,0,0); the closest rectangle point
        // is the edge point (2,0,0).
        const result = query.compute(line([5, 0, -1], [0, 0, 1]), rect);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(2, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 10);
        expect(result.cartesian[0]).toBeCloseTo(2, 10);
        expect(result.cartesian[1]).toBeCloseTo(0, 10);
    });

    it('reports a corner when the plane intersection is past a corner', () => {
        const result = query.compute(line([6, 5, -1], [0, 0, 1]), rect);
        expect(result.closest[1].values[0]).toBeCloseTo(2, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 10);
        expect(result.distance).toBeCloseTo(Math.sqrt(16 + 16), 10);
    });

    it('reports cartesian coordinates consistent with the closest point',
        () => {
            const result = query.compute(line([6, 5, -1], [0, 0, 1]), rect);
            const q = add(rect.center,
                add(mul(result.cartesian[0], rect.axis[0]),
                    mul(result.cartesian[1], rect.axis[1])));
            for (let i = 0; i < 3; ++i) {
                expect(q.values[i]).toBeCloseTo(result.closest[1].values[i],
                    9);
            }
        });

    it('handles a rectangle with a zero extent (a segment)', () => {
        const degenerate = Rectangle.fromCenterAxisExtent(v(0, 0, 0),
            [v(1, 0, 0), v(0, 1, 0)], v(2, 0));
        const result = query.compute(line([0, 3, 0], [1, 0, 0]), degenerate);
        expect(result.distance).toBeCloseTo(3, 10);
    });

    it('agrees with a dense sampling of the rectangle', () => {
        const rnd = makeRandom(31337);
        const c = Math.SQRT1_2;
        const axis = [v(c, c, 0), v(-c / Math.SQRT2, c / Math.SQRT2,
            Math.SQRT1_2)];
        // Re-orthonormalize the second axis against the first.
        const a0 = axis[0];
        let a1 = sub(axis[1], mul(dot(axis[1], a0), a0));
        a1 = mul(a1, 1 / Math.sqrt(dot(a1, a1)));
        const extent = v(1.5, 0.75);
        const r = Rectangle.fromCenterAxisExtent(v(0.25, -0.5, 0.75),
            [a0, a1], extent);

        for (let trial = 0; trial < 30; ++trial) {
            const origin = v(6 * rnd() - 3, 6 * rnd() - 3, 6 * rnd() - 3);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const ln = Line.fromOriginDirection(origin, dir);
            const result = query.compute(ln, r);

            // The cartesian coordinates are in range and reproduce closest[1].
            expect(Math.abs(result.cartesian[0])).toBeLessThanOrEqual(
                extent.values[0] + 1e-8);
            expect(Math.abs(result.cartesian[1])).toBeLessThanOrEqual(
                extent.values[1] + 1e-8);
            const q = add(r.center, add(mul(result.cartesian[0], r.axis[0]),
                mul(result.cartesian[1], r.axis[1])));
            for (let i = 0; i < 3; ++i) {
                expect(q.values[i]).toBeCloseTo(result.closest[1].values[i],
                    7);
            }

            const onLine = add(ln.origin, mul(result.parameter, ln.direction));
            for (let i = 0; i < 3; ++i) {
                expect(onLine.values[i]).toBeCloseTo(
                    result.closest[0].values[i], 7);
            }

            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 7);

            // No sampled rectangle point is closer to the line.
            const n = 60;
            const dd = dot(ln.direction, ln.direction);
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= n; ++i) {
                for (let j = 0; j <= n; ++j) {
                    const p = add(r.center, add(
                        mul((2 * i / n - 1) * extent.values[0], r.axis[0]),
                        mul((2 * j / n - 1) * extent.values[1], r.axis[1])));
                    const w = sub(p, ln.origin);
                    const s = dot(w, ln.direction) / dd;
                    const f = sub(w, mul(s, ln.direction));
                    best = Math.min(best, dot(f, f));
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});
