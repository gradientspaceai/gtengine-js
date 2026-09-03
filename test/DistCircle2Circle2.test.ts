import { describe, expect, it } from 'vitest';
import { DistCircle2Circle2 } from '../src/DistCircle2Circle2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector, dot, sub } from '../src/Vector.js';

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
