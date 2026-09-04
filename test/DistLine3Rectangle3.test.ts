import { describe, expect, it } from 'vitest';
import { DistLine3Rectangle3 } from '../src/DistLine3Rectangle3.js';
import { Line } from '../src/Line.js';
import { Rectangle } from '../src/Rectangle.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistLine3Rectangle3.h.
// ---------------------------------------------------------------------------

describe('DistLine3Rectangle3 verification', () => {
    const query = new DistLine3Rectangle3();

    const rectArb = fc.tuple(wellScaledVector(3, -5, 5), rotationFrame(3),
        fc.tuple(finite(0, 4), finite(0, 4)))
        .map(([c, frame, e]) => Rectangle.fromCenterAxisExtent(c,
            [frame[0], frame[1]], v(e[0], e[1])));

    const lineArb = fc.tuple(wellScaledVector(3, -8, 8), unitVector(3))
        .map(([o, d]) => Line.fromOriginDirection(o, d));

    // Independent closed-form distance from a point to a solid rectangle
    // whose axes are orthonormal: clamp the in-plane coordinates.
    function pointRectDistance(p: Vector, rect: Rectangle): number {
        const delta = sub(p, rect.center);
        let inPlane = new Vector(3);
        for (let i = 0; i < 2; ++i) {
            const y = dot(rect.axis[i], delta);
            const e = rect.extent.values[i];
            const s = Math.min(Math.max(y, -e), e);
            inPlane = add(inPlane, mul(s, rect.axis[i]));
        }
        return length(sub(delta, inPlane));
    }

    function ternaryMin(f: (t: number) => number, lo: number,
        hi: number): number {
        let a = lo, b = hi;
        for (let i = 0; i < 200; ++i) {
            const m0 = a + (b - a) / 3;
            const m1 = b - (b - a) / 3;
            if (f(m0) <= f(m1)) { b = m1; } else { a = m0; }
        }
        return f(0.5 * (a + b));
    }

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(lineArb, rectArb), ([ln, rect]) => {
            const r = query.compute(ln, rect);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            expectClose(length(sub(r.closest[0], r.closest[1])), r.distance,
                1e-9, 1e-9);
            expectVectorClose(r.closest[0],
                add(ln.origin, mul(r.parameter, ln.direction)), 1e-9, 1e-9);
            // closest[1] is the rectangle point with coordinates cartesian[].
            expect(Math.abs(r.cartesian[0]))
                .toBeLessThanOrEqual(rect.extent.values[0] + 1e-9);
            expect(Math.abs(r.cartesian[1]))
                .toBeLessThanOrEqual(rect.extent.values[1] + 1e-9);
            expectVectorClose(r.closest[1], add(rect.center,
                add(mul(r.cartesian[0], rect.axis[0]),
                    mul(r.cartesian[1], rect.axis[1]))), 1e-8, 1e-8);
        });
    });

    it('matches an independent convex minimization along the line', () => {
        check(fc.tuple(lineArb, rectArb), ([ln, rect]) => {
            const r = query.compute(ln, rect);
            const best = ternaryMin(
                t => pointRectDistance(add(ln.origin, mul(t, ln.direction)),
                    rect), -100, 100);
            expectClose(r.distance, best, 1e-7, 1e-7);
        }, 100);
    });

    it('reports zero distance when the line crosses the rectangle', () => {
        check(fc.tuple(rectArb, finite(-0.95, 0.95), finite(-0.95, 0.95),
            unitVector(3)), ([rect, u0, u1, dir]) => {
            const target = add(rect.center,
                add(mul(u0 * rect.extent.values[0], rect.axis[0]),
                    mul(u1 * rect.extent.values[1], rect.axis[1])));
            const ln = Line.fromOriginDirection(add(target, mul(3, dir)), dir);
            const r = query.compute(ln, rect);
            expectClose(r.distance, 0, 1e-9, 1e-9);
        });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(lineArb, rectArb, rotationFrame(3),
            wellScaledVector(3, -6, 6)), ([ln, rect, frame, shift]) => {
            const rot = (p: Vector): Vector =>
                add(add(mul(p.values[0], frame[0]), mul(p.values[1], frame[1])),
                    mul(p.values[2], frame[2]));
            const movedLine = Line.fromOriginDirection(
                add(shift, rot(ln.origin)), rot(ln.direction));
            const movedRect = Rectangle.fromCenterAxisExtent(
                add(shift, rot(rect.center)),
                [rot(rect.axis[0]), rot(rect.axis[1])], rect.extent);
            const r0 = query.compute(ln, rect);
            const r1 = query.compute(movedLine, movedRect);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
        });
    });

    it('reduces to the point-line distance for a degenerate rectangle', () => {
        check(fc.tuple(lineArb, wellScaledVector(3, -5, 5), rotationFrame(3)),
            ([ln, c, frame]) => {
                const rect = Rectangle.fromCenterAxisExtent(c,
                    [frame[0], frame[1]], v(0, 0));
                const r = query.compute(ln, rect);
                const diff = sub(c, ln.origin);
                const t = dot(diff, ln.direction)
                    / dot(ln.direction, ln.direction);
                expectClose(r.distance,
                    length(sub(diff, mul(t, ln.direction))), 1e-9, 1e-9);
                expectVectorClose(r.closest[1], c, 1e-9, 1e-9);
            });
    });

    it('handles a line parallel to the plane of the rectangle', () => {
        check(fc.tuple(rectArb, finite(-3, 3), finite(0.1, 5),
            finite(-Math.PI, Math.PI)), ([rect, offset, height, angle]) => {
            const normal = cross(rect.axis[0], rect.axis[1]);
            const inPlane = add(mul(Math.cos(angle), rect.axis[0]),
                mul(Math.sin(angle), rect.axis[1]));
            const perp = add(mul(-Math.sin(angle), rect.axis[0]),
                mul(Math.cos(angle), rect.axis[1]));
            const origin = add(rect.center,
                add(mul(height, normal), mul(offset, perp)));
            const ln = Line.fromOriginDirection(origin, inPlane);
            const r = query.compute(ln, rect);
            const best = ternaryMin(
                t => pointRectDistance(add(ln.origin, mul(t, ln.direction)),
                    rect), -100, 100);
            expectClose(r.distance, best, 1e-7, 1e-7);
        }, 100);
    });
});
