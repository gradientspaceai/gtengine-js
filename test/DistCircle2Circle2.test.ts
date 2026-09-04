import { describe, expect, it } from 'vitest';
import { DistCircle2Circle2 } from '../src/DistCircle2Circle2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { check, expectClose, expectVectorClose, fc, positive, rotationFrame, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function circle(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(v(...center), radius);
}

function distanceToCenter(p: Vector, c: Hypersphere): number {
    const d = sub(p, c.center);
    return Math.sqrt(dot(d, d));
}

describe('DistCircle2Circle2', () => {
    const query = new DistCircle2Circle2();

    it('measures strictly separated circles', () => {
        const c0 = circle([0, 0], 1);
        const c1 = circle([10, 0], 2);
        const result = query.compute(c0, c1);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(7, 12);
        expect(result.concentric).toBe(false);
        // closest[0][0] is on the larger circle (c1), closest[0][1] on c0.
        expect(distanceToCenter(result.closest[0][0], c1)).toBeCloseTo(2, 10);
        expect(distanceToCenter(result.closest[0][1], c0)).toBeCloseTo(1, 10);
    });

    it('reports a single coincident pair for externally tangent circles',
        () => {
            const c0 = circle([0, 0], 2);
            const c1 = circle([5, 0], 3);
            const result = query.compute(c0, c1);
            expect(result.numClosestPairs).toBe(1);
            expect(result.distance).toBeCloseTo(0, 12);
            expect(result.closest[0][0].values[0]).toBeCloseTo(2, 10);
            expect(result.closest[0][1].values[0]).toBeCloseTo(2, 10);
        });

    it('measures a strictly nested circle', () => {
        const c0 = circle([0, 0], 10);
        const c1 = circle([1, 0], 2);
        const result = query.compute(c0, c1);
        expect(result.numClosestPairs).toBe(1);
        // r0 - r1 - |C1-C0| = 10 - 2 - 1 = 7.
        expect(result.distance).toBeCloseTo(7, 12);
        expect(distanceToCenter(result.closest[0][0], c0)).toBeCloseTo(10, 10);
        expect(distanceToCenter(result.closest[0][1], c1)).toBeCloseTo(2, 10);
    });

    it('reports a single coincident pair for internally tangent circles',
        () => {
            const c0 = circle([0, 0], 5);
            const c1 = circle([2, 0], 3);
            const result = query.compute(c0, c1);
            expect(result.numClosestPairs).toBe(1);
            expect(result.distance).toBeCloseTo(0, 12);
            expect(result.closest[0][0].values[0]).toBeCloseTo(5, 10);
            expect(result.closest[0][1].values[0]).toBeCloseTo(5, 10);
        });

    it('reports the two intersection points for crossing circles', () => {
        // Unit circles at (0,0) and (1,0) meet at (0.5, +/- sqrt(3)/2).
        const c0 = circle([0, 0], 1);
        const c1 = circle([1, 0], 1);
        const result = query.compute(c0, c1);
        expect(result.numClosestPairs).toBe(2);
        expect(result.distance).toBe(0);
        const ys = [result.closest[0][0].values[1],
            result.closest[1][0].values[1]].sort((a, b) => a - b);
        expect(result.closest[0][0].values[0]).toBeCloseTo(0.5, 10);
        expect(result.closest[1][0].values[0]).toBeCloseTo(0.5, 10);
        expect(ys[0]).toBeCloseTo(-Math.sqrt(3) / 2, 10);
        expect(ys[1]).toBeCloseTo(Math.sqrt(3) / 2, 10);
        // Both intersection points lie on both circles.
        for (const pair of result.closest) {
            expect(distanceToCenter(pair[0], c0)).toBeCloseTo(1, 9);
            expect(distanceToCenter(pair[0], c1)).toBeCloseTo(1, 9);
        }
    });

    it('flags concentric circles', () => {
        const c0 = circle([3, 4], 5);
        const c1 = circle([3, 4], 2);
        const result = query.compute(c0, c1);
        expect(result.concentric).toBe(true);
        expect(result.cocircular).toBe(false);
        expect(result.numClosestPairs).toBe(2);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[0][0].values).toEqual([-2, 4]);
        expect(result.closest[0][1].values).toEqual([1, 4]);
        expect(result.closest[1][0].values).toEqual([8, 4]);
        expect(result.closest[1][1].values).toEqual([5, 4]);
    });

    it('flags cocircular circles', () => {
        const result = query.compute(circle([1, 1], 3), circle([1, 1], 3));
        expect(result.concentric).toBe(true);
        expect(result.cocircular).toBe(true);
        expect(result.distance).toBe(0);
        expect(result.numClosestPairs).toBe(2);
    });

    it('gives the same distance whichever circle is passed first', () => {
        const c0 = circle([0, 0], 1);
        const c1 = circle([4, 3], 2);
        const a = query.compute(c0, c1);
        const b = query.compute(c1, c0);
        expect(b.distance).toBeCloseTo(a.distance, 12);
        expect(b.numClosestPairs).toBe(a.numClosestPairs);
    });

    it('always reports points that lie on their circles', () => {
        let seed = 1357;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 8 - 4;
        };
        for (let trial = 0; trial < 120; ++trial) {
            const c0 = circle([rand(), rand()], Math.abs(rand()) + 0.2);
            const c1 = circle([rand(), rand()], Math.abs(rand()) + 0.2);
            const result = query.compute(c0, c1);
            // The larger circle is used as 'circle0' internally, so identify
            // the roles by radius.
            const big = c0.radius >= c1.radius ? c0 : c1;
            const small = c0.radius >= c1.radius ? c1 : c0;
            for (let j = 0; j < result.numClosestPairs; ++j) {
                expect(distanceToCenter(result.closest[j][0], big))
                    .toBeCloseTo(big.radius, 8);
                expect(distanceToCenter(result.closest[j][1], small))
                    .toBeCloseTo(small.radius, 8);
            }
            const diff = sub(result.closest[0][0], result.closest[0][1]);
            expect(Math.sqrt(dot(diff, diff))).toBeCloseTo(result.distance, 8);
            expect(result.sqrDistance).toBeCloseTo(
                result.distance * result.distance, 8);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of
// DistCircle2Circle2.ts against the upstream header DistCircle2Circle2.h.
// ---------------------------------------------------------------------------

const circle2 = fc.tuple(wellScaledVector(2, -8, 8), positive(6))
    .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));

// Upstream's operator() passes the larger-radius circle to DoQuery first, and
// DoQuery writes the point on its own first argument into closest[j][0].
// Roles are therefore identified by radius, not by argument position; see the
// "argument order" test below.
function byRadius(c0: Hypersphere, c1: Hypersphere):
    [Hypersphere, Hypersphere] {
    return c0.radius >= c1.radius ? [c0, c1] : [c1, c0];
}

const circlePair = fc.tuple(circle2, circle2);

describe('DistCircle2Circle2 verification', () => {
    const query = new DistCircle2Circle2();

    it('matches the closed-form separated/contained/intersecting cases', () => {
        check(circlePair, ([c0, c1]) => {
            const r = query.compute(c0, c1);
            const pair = byRadius(c0, c1);
            const d = length(sub(pair[1].center, pair[0].center));
            const rSum = pair[0].radius + pair[1].radius;
            const rDif = pair[0].radius - pair[1].radius;
            if (d === 0) {
                expect(r.concentric).toBe(true);
                expect(r.numClosestPairs).toBe(2);
                expectClose(r.distance, Math.abs(rDif), 1e-12, 1e-12);
                expect(r.cocircular).toBe(c0.radius === c1.radius);
            }
            else if (d >= rSum) {
                expect(r.numClosestPairs).toBe(1);
                expectClose(r.distance, d - rSum, 1e-12, 1e-12);
            }
            else if (d <= rDif) {
                expect(r.numClosestPairs).toBe(1);
                expectClose(r.distance, rDif - d, 1e-12, 1e-12);
            }
            else {
                expect(r.numClosestPairs).toBe(2);
                expect(r.distance).toBe(0);
            }
            expectClose(r.sqrDistance, r.distance * r.distance, 1e-12, 1e-12);
            if (d !== 0) {
                expect(r.concentric).toBe(false);
                expect(r.cocircular).toBe(false);
            }
        });
    });

    it('places every reported point on its own circle', () => {
        check(circlePair, ([c0, c1]) => {
            const r = query.compute(c0, c1);
            const pair = byRadius(c0, c1);
            for (let j = 0; j < r.numClosestPairs; ++j) {
                expectClose(length(sub(r.closest[j][0], pair[0].center)),
                    pair[0].radius, 1e-8, 1e-8);
                expectClose(length(sub(r.closest[j][1], pair[1].center)),
                    pair[1].radius, 1e-8, 1e-8);
                // The reported distance is realized by the pair.
                expectClose(length(sub(r.closest[j][0], r.closest[j][1])),
                    r.distance, 1e-8, 1e-8);
            }
        });
    });

    it('is minimal over sampled circle point pairs', () => {
        check(circlePair, ([c0, c1]) => {
            const r = query.compute(c0, c1);
            const onCircle = (c: Hypersphere, a: number): Vector =>
                add(c.center, Vector.fromArray(
                    [c.radius * Math.cos(a), c.radius * Math.sin(a)]));
            for (let i = 0; i < 24; ++i) {
                for (let k = 0; k < 24; ++k) {
                    const p = onCircle(c0, (2 * Math.PI * i) / 24);
                    const q = onCircle(c1, (2 * Math.PI * k) / 24);
                    const gap = length(sub(p, q));
                    expect(r.distance)
                        .toBeLessThanOrEqual(gap + 1e-9 * (1 + gap));
                }
            }
        }, 40);
    });

    it('reports the same distance for either argument order', () => {
        check(circlePair, ([c0, c1]) => {
            const a = query.compute(c0, c1);
            const b = query.compute(c1, c0);
            expectClose(a.distance, b.distance, 0, 0);
            expect(a.numClosestPairs).toBe(b.numClosestPairs);
            expect(a.concentric).toBe(b.concentric);
            expect(a.cocircular).toBe(b.cocircular);
        });
    });

    it('keys closest[j][0] to the larger circle, not to argument 0', () => {
        // Upstream issue: operator() swaps its arguments when
        // circle1.radius > circle0.radius but DoQuery writes the point on its
        // own first argument into closest[j][0], so passing the smaller
        // circle first puts the *second* circle's point in closest[j][0].
        // This contradicts the upstream file comment and is preserved here.
        const small = Hypersphere.fromCenterRadius(v(0, 0), 1);
        const big = Hypersphere.fromCenterRadius(v(10, 0), 3);
        const r = query.compute(small, big);
        expect(r.numClosestPairs).toBe(1);
        // closest[0][0] is on the big circle even though it was argument 1.
        expectClose(length(sub(r.closest[0][0], big.center)), big.radius,
            1e-12, 1e-12);
        expectClose(length(sub(r.closest[0][1], small.center)), small.radius,
            1e-12, 1e-12);
        expectClose(r.distance, 10 - 4, 1e-12, 1e-12);
    });

    it('reports both intersection points when the circles cross', () => {
        check(circlePair.filter(([c0, c1]) => {
            const pair = byRadius(c0, c1);
            const d = length(sub(pair[1].center, pair[0].center));
            return d > 0 && d < pair[0].radius + pair[1].radius
                && d > pair[0].radius - pair[1].radius;
        }), ([c0, c1]) => {
            const r = query.compute(c0, c1);
            expect(r.numClosestPairs).toBe(2);
            expect(r.distance).toBe(0);
            for (let j = 0; j < 2; ++j) {
                // Each pair is a single point lying on both circles.
                expectVectorClose(r.closest[j][0], r.closest[j][1], 0, 0);
                expectClose(length(sub(r.closest[j][0], c0.center)),
                    c0.radius, 1e-7, 1e-7);
                expectClose(length(sub(r.closest[j][0], c1.center)),
                    c1.radius, 1e-7, 1e-7);
            }
        });
    });

    it('reports the documented pairs for concentric circles', () => {
        check(fc.tuple(wellScaledVector(2, -8, 8), positive(6), positive(6)),
            ([c, ra, rb]) => {
                const c0 = Hypersphere.fromCenterRadius(c, ra);
                const c1 = Hypersphere.fromCenterRadius(c.clone(), rb);
                const r = query.compute(c0, c1);
                const pair = byRadius(c0, c1);
                expect(r.concentric).toBe(true);
                expect(r.cocircular).toBe(ra === rb);
                expect(r.numClosestPairs).toBe(2);
                expectClose(r.distance, Math.abs(ra - rb), 1e-12, 1e-12);
                expectVectorClose(r.closest[0][0],
                    sub(c, Vector.fromArray([pair[0].radius, 0])), 0, 0);
                expectVectorClose(r.closest[1][1],
                    add(c, Vector.fromArray([pair[1].radius, 0])), 0, 0);
            });
    });

    it('is equivariant under rigid motions of the plane', () => {
        check(fc.tuple(circlePair, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([pair, R, tr]) => {
            const move = (c: Hypersphere): Hypersphere =>
                Hypersphere.fromCenterRadius(
                    add(add(mul(c.center.values[0], R[0]),
                        mul(c.center.values[1], R[1])), tr), c.radius);
            expectClose(query.compute(pair[0], pair[1]).distance,
                query.compute(move(pair[0]), move(pair[1])).distance, 1e-8,
                1e-8);
        });
    });
});
