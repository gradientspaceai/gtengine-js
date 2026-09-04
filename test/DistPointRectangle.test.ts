import { describe, expect, it } from 'vitest';
import { DistPointRectangle } from '../src/DistPointRectangle.js';
import { Rectangle } from '../src/Rectangle.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { check, expectClose, expectVectorClose, fc, positive, rotationFrame, seededRandom, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

// A rectangle in 3D whose plane is z = 0, centered at the origin, with the
// standard axes and the given extents.
function unitRect(e0: number, e1: number): Rectangle {
    return Rectangle.fromCenterAxisExtent(v(0, 0, 0),
        [v(1, 0, 0), v(0, 1, 0)], v(e0, e1));
}

describe('DistPointRectangle', () => {
    const query = new DistPointRectangle();

    it('measures the plane offset for a point over the rectangle', () => {
        const result = query.compute(v(0.5, -0.25, 3), unitRect(1, 2));
        expect(result.cartesian).toEqual([0.5, -0.25]);
        expect(result.closest[1].values).toEqual([0.5, -0.25, 0]);
        expect(result.distance).toBeCloseTo(3, 12);
    });

    it('reports zero distance for a point on the rectangle', () => {
        const result = query.compute(v(-1, 2, 0), unitRect(1, 2));
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.cartesian).toEqual([-1, 2]);
    });

    it('clamps the coordinates outside the extents', () => {
        const result = query.compute(v(5, -7, 0), unitRect(1, 2));
        expect(result.cartesian).toEqual([1, -2]);
        expect(result.closest[1].values).toEqual([1, -2, 0]);
        expect(result.distance).toBeCloseTo(Math.sqrt(16 + 25), 12);
    });

    it('clamps a single coordinate (edge region)', () => {
        const result = query.compute(v(0.5, 9, 4), unitRect(1, 2));
        expect(result.cartesian).toEqual([0.5, 2]);
        expect(result.distance).toBeCloseTo(Math.sqrt(49 + 16), 12);
    });

    it('handles a rotated rectangle in 2D', () => {
        const a0 = v(1, 1);
        const a1 = v(-1, 1);
        normalize(a0);
        normalize(a1);
        const rect = Rectangle.fromCenterAxisExtent(v(2, 3), [a0, a1],
            v(1, 1));
        // The point is 5 units along a0 from the center, so it is clamped to
        // the extent 1 and the distance is 4.
        const point = add(rect.center, mul(5, a0));
        const result = query.compute(point, rect);
        expect(result.cartesian[0]).toBeCloseTo(1, 12);
        expect(result.cartesian[1]).toBeCloseTo(0, 12);
        expect(result.distance).toBeCloseTo(4, 12);
    });

    it('agrees with a sampled minimum over the rectangle', () => {
        let seed = 777;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 8 - 4;
        };
        const rect = unitRect(1.5, 0.75);
        for (let trial = 0; trial < 40; ++trial) {
            const point = v(rand(), rand(), rand());
            const result = query.compute(point, rect);

            expect(Math.abs(result.cartesian[0])).toBeLessThanOrEqual(1.5);
            expect(Math.abs(result.cartesian[1])).toBeLessThanOrEqual(0.75);

            let best = Number.MAX_VALUE;
            for (let i = 0; i <= 120; ++i) {
                const s0 = -1.5 + 3 * i / 120;
                for (let j = 0; j <= 120; ++j) {
                    const s1 = -0.75 + 1.5 * j / 120;
                    const q = add(rect.center,
                        add(mul(s0, rect.axis[0]), mul(s1, rect.axis[1])));
                    const diff = sub(point, q);
                    best = Math.min(best, dot(diff, diff));
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of
// DistPointRectangle.ts against the upstream header DistPointRectangle.h.
// ---------------------------------------------------------------------------

function rot(R: readonly Vector[], p: Vector): Vector {
    let q = mul(p.values[0], R[0]);
    for (let i = 1; i < R.length; ++i) {
        q = add(q, mul(p.values[i], R[i]));
    }
    return q;
}

// A 3D rectangle whose two axes are the first two vectors of a rotation
// frame, so the axes are exactly unit length and orthogonal (upstream
// requires unit-length W[]).
const rectangle3 = fc.tuple(wellScaledVector(3, -8, 8), rotationFrame(3),
    positive(5), positive(5))
    .map(([c, R, e0, e1]) => Rectangle.fromCenterAxisExtent(c, [R[0], R[1]],
        Vector.fromArray([e0, e1])));

describe('DistPointRectangle verification', () => {
    const query = new DistPointRectangle();

    it('result is self consistent and the cartesian coordinates are clamped',
        () => {
            check(fc.tuple(wellScaledVector(3, -8, 8), rectangle3),
                ([p, rect]) => {
                    const r = query.compute(p, rect);
                    expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12,
                        1e-12);
                    const diff = sub(r.closest[0], r.closest[1]);
                    expectClose(r.sqrDistance, dot(diff, diff), 1e-12, 1e-12);
                    expectVectorClose(r.closest[0], p, 0, 0);
                    expect(r.closest[0]).not.toBe(p);
                    for (let i = 0; i < 2; ++i) {
                        expect(Math.abs(r.cartesian[i]))
                            .toBeLessThanOrEqual(rect.extent.values[i] + 1e-12);
                    }
                    // closest[1] is reconstructed from the cartesian
                    // coordinates.
                    const q = add(rect.center,
                        add(mul(r.cartesian[0], rect.axis[0]),
                            mul(r.cartesian[1], rect.axis[1])));
                    expectVectorClose(r.closest[1], q, 1e-12, 1e-12);
                });
        });

    it('matches the independent clamp-in-frame formula', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), rectangle3), ([p, rect]) => {
            const diff = sub(p, rect.center);
            let expected = dot(diff, diff);
            for (let i = 0; i < 2; ++i) {
                const t = dot(rect.axis[i], diff);
                const s = Math.min(Math.max(t, -rect.extent.values[i]),
                    rect.extent.values[i]);
                // Remove the in-plane component and add back the clamped
                // residual: the axes are orthonormal, so the contributions
                // are independent.
                expected += s * s - 2 * s * t;
            }
            expectClose(query.compute(p, rect).sqrDistance, expected, 1e-9,
                1e-9);
        });
    });

    it('is minimal over sampled rectangle points', () => {
        const rand = seededRandom(0x51d4);
        check(fc.tuple(wellScaledVector(3, -8, 8), rectangle3), ([p, rect]) => {
            const r = query.compute(p, rect);
            for (let k = 0; k < 24; ++k) {
                const s0 = rect.extent.values[0] * (2 * rand() - 1);
                const s1 = rect.extent.values[1] * (2 * rand() - 1);
                const q = add(rect.center, add(mul(s0, rect.axis[0]),
                    mul(s1, rect.axis[1])));
                const diff = sub(p, q);
                const sqr = dot(diff, diff);
                expect(r.sqrDistance).toBeLessThanOrEqual(sqr + 1e-9 * (1 + sqr));
            }
        }, 60);
    });

    it('is invariant under rigid motions', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), rectangle3, rotationFrame(3),
            wellScaledVector(3, -5, 5)), ([p, rect, R, tr]) => {
            const moved = Rectangle.fromCenterAxisExtent(
                add(rot(R, rect.center), tr),
                [rot(R, rect.axis[0]), rot(R, rect.axis[1])], rect.extent);
            const r0 = query.compute(p, rect);
            const r1 = query.compute(add(rot(R, p), tr), moved);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expectClose(r0.cartesian[0], r1.cartesian[0], 1e-9, 1e-8);
            expectClose(r0.cartesian[1], r1.cartesian[1], 1e-9, 1e-8);
        });
    });

    it('reports zero distance for points of the rectangle', () => {
        check(fc.tuple(rectangle3, fc.double({ min: -1, max: 1, noNaN: true }),
            fc.double({ min: -1, max: 1, noNaN: true })),
        ([rect, f0, f1]) => {
            const s0 = f0 * rect.extent.values[0];
            const s1 = f1 * rect.extent.values[1];
            const p = add(rect.center, add(mul(s0, rect.axis[0]),
                mul(s1, rect.axis[1])));
            const r = query.compute(p, rect);
            expect(r.distance).toBeLessThanOrEqual(1e-12);
            expectClose(r.cartesian[0], s0, 1e-9, 1e-9);
            expectClose(r.cartesian[1], s1, 1e-9, 1e-9);
        });
    });
});
