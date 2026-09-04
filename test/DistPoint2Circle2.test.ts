import { describe, expect, it } from 'vitest';
import { DistPoint2Circle2 } from '../src/DistPoint2Circle2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { check, expectClose, expectVectorClose, fc, positive, rotationFrame, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function circle(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(v(...center), radius);
}

describe('DistPoint2Circle2', () => {
    const query = new DistPoint2Circle2();

    it('measures the distance from a point outside the circle', () => {
        const result = query.compute(v(5, 0), circle([0, 0], 2));
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.sqrDistance).toBeCloseTo(9, 12);
        expect(result.equidistant).toBe(false);
        expect(result.closest[1].values[0]).toBeCloseTo(2, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 12);
    });

    it('measures the distance from a point inside the circle', () => {
        const result = query.compute(v(0, 1), circle([0, 0], 4));
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(0, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(4, 12);
    });

    it('reports zero distance for a point on the circle', () => {
        const result = query.compute(v(3, 4), circle([0, 0], 5));
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.equidistant).toBe(false);
    });

    it('flags the circle center as equidistant', () => {
        const result = query.compute(v(2, -3), circle([2, -3], 7));
        expect(result.equidistant).toBe(true);
        expect(result.distance).toBe(7);
        expect(result.sqrDistance).toBe(49);
        // The reported point is C + r*(1,0).
        expect(result.closest[1].values).toEqual([9, -3]);
    });

    it('handles a zero-radius (degenerate) circle', () => {
        const result = query.compute(v(3, 4), circle([0, 0], 0));
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(0, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 12);
    });

    it('places the closest point on the circle at the reported distance',
        () => {
            let seed = 20250;
            const rand = () => {
                seed = (seed * 1103515245 + 12345) % 2147483648;
                return seed / 2147483648 * 10 - 5;
            };
            const c = circle([1, -2], 2.5);
            for (let trial = 0; trial < 100; ++trial) {
                const point = v(rand(), rand());
                const result = query.compute(point, c);

                // The closest point is on the circle.
                const radial = sub(result.closest[1], c.center);
                expect(Math.sqrt(dot(radial, radial))).toBeCloseTo(c.radius,
                    9);

                // It realizes the reported distance.
                const diff = sub(result.closest[0], result.closest[1]);
                expect(Math.sqrt(dot(diff, diff))).toBeCloseTo(
                    result.distance, 9);

                // No sampled circle point is closer.
                for (let k = 0; k < 720; ++k) {
                    const s = k * Math.PI / 360;
                    const q = v(c.center.values[0] + c.radius * Math.cos(s),
                        c.center.values[1] + c.radius * Math.sin(s));
                    const d = sub(point, q);
                    expect(dot(d, d)).toBeGreaterThanOrEqual(
                        result.sqrDistance - 1e-9);
                }
            }
        });
});

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of
// DistPoint2Circle2.ts against the upstream header DistPoint2Circle2.h.
// ---------------------------------------------------------------------------

function rot(R: readonly Vector[], p: Vector): Vector {
    return add(mul(p.values[0], R[0]), mul(p.values[1], R[1]));
}

const circle2 = fc.tuple(wellScaledVector(2, -8, 8), positive(6))
    .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));

/**
 * A (point, circle) pair with the point strictly off the center, so the
 * single-closest-point branch is taken. wellScaledVector snaps tiny
 * components to exactly zero, so an unfiltered pair can hit the equidistant
 * branch.
 */
const offCenter = fc.tuple(wellScaledVector(2, -8, 8), circle2)
    .filter(([p, circle]) => p.notEquals(circle.center));

describe('DistPoint2Circle2 verification', () => {
    const query = new DistPoint2Circle2();

    it('matches | |P-C| - r | and puts closest[1] on the circle', () => {
        check(offCenter, ([p, circle]) => {
            const r = query.compute(p, circle);
            const diff = sub(p, circle.center);
            const expected = Math.abs(length(diff) - circle.radius);
            expectClose(r.distance, expected, 1e-12, 1e-12);
            expectClose(r.sqrDistance, r.distance * r.distance, 1e-12, 1e-12);
            expectVectorClose(r.closest[0], p, 0, 0);
            expect(r.closest[0]).not.toBe(p);
            // The closest circle point is at distance r from the center.
            expectClose(length(sub(r.closest[1], circle.center)), circle.radius,
                1e-9, 1e-9);
            // ... and the reported distance is |closest[0] - closest[1]|.
            expectClose(r.distance, length(sub(r.closest[0], r.closest[1])),
                1e-9, 1e-9);
            expect(r.equidistant).toBe(false);
        });
    });

    it('is minimal over sampled circle points', () => {
        check(offCenter, ([p, circle]) => {
            const r = query.compute(p, circle);
            for (let k = 0; k < 32; ++k) {
                const a = (2 * Math.PI * k) / 32;
                const q = add(circle.center, Vector.fromArray(
                    [circle.radius * Math.cos(a), circle.radius * Math.sin(a)]));
                const diff = sub(p, q);
                expect(r.distance).toBeLessThanOrEqual(
                    length(diff) + 1e-9 * (1 + length(diff)));
            }
        }, 60);
    });

    it('flags the circle center as equidistant and reports C + r*(1,0)', () => {
        check(circle2, circle => {
            const r = query.compute(circle.center.clone(), circle);
            expect(r.equidistant).toBe(true);
            expect(r.distance).toBe(circle.radius);
            expect(r.sqrDistance).toBe(circle.radius * circle.radius);
            expectVectorClose(r.closest[1],
                add(circle.center, Vector.fromArray([circle.radius, 0])), 0, 0);
        });
    });

    it('is equivariant under rotations about the circle center', () => {
        check(fc.tuple(offCenter, rotationFrame(2)),
            ([[p, circle], R]) => {
                const r0 = query.compute(p, circle);
                const moved = Hypersphere.fromCenterRadius(
                    rot(R, circle.center), circle.radius);
                const r1 = query.compute(rot(R, p), moved);
                expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
                expectVectorClose(rot(R, r0.closest[1]), r1.closest[1], 1e-8,
                    1e-8);
            });
    });
});
