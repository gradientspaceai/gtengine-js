import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistPointCanonicalBox } from '../src/DistPointCanonicalBox.js';
import { Vector, dot, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(...extent: number[]): CanonicalBox {
    return CanonicalBox.fromExtent(v(...extent));
}

describe('DistPointCanonicalBox', () => {
    const query = new DistPointCanonicalBox();

    it('reports zero distance for a point inside the box', () => {
        const result = query.compute(v(0.5, -0.5, 0.25), box(1, 1, 1));
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
        expect(result.closest[1].values).toEqual([0.5, -0.5, 0.25]);
    });

    it('reports zero distance for a point on the boundary', () => {
        const result = query.compute(v(1, 0, 0), box(1, 2, 3));
        expect(result.distance).toBe(0);
    });

    it('measures a face-region distance', () => {
        const result = query.compute(v(4, 0.5, -0.5), box(1, 1, 1));
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[1].values).toEqual([1, 0.5, -0.5]);
    });

    it('measures an edge-region distance', () => {
        const result = query.compute(v(4, 5, 0), box(1, 1, 1));
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values).toEqual([1, 1, 0]);
    });

    it('measures a corner-region distance', () => {
        const result = query.compute(v(-3, -5, -13), box(1, 1, 1));
        expect(result.closest[1].values).toEqual([-1, -1, -1]);
        expect(result.sqrDistance).toBeCloseTo(4 + 16 + 144, 12);
    });

    it('works in 2D and 5D', () => {
        const r2 = query.compute(v(3, 4), box(0, 0));
        expect(r2.distance).toBeCloseTo(5, 12);
        const r5 = query.compute(v(2, 2, 2, 2, 2), box(1, 1, 1, 1, 1));
        expect(r5.sqrDistance).toBeCloseTo(5, 12);
    });

    it('leaves the input point untouched in closest[0]', () => {
        const point = v(9, 9, 9);
        const result = query.compute(point, box(1, 1, 1));
        expect(result.closest[0].values).toEqual([9, 9, 9]);
        expect(point.values).toEqual([9, 9, 9]);
    });

    it('agrees with a clamped-coordinate reference', () => {
        let seed = 31337;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 8 - 4;
        };
        const b = box(1.5, 0.5, 2);
        for (let trial = 0; trial < 200; ++trial) {
            const point = v(rand(), rand(), rand());
            const result = query.compute(point, b);
            const clamped = v(
                Math.min(Math.max(point.values[0], -1.5), 1.5),
                Math.min(Math.max(point.values[1], -0.5), 0.5),
                Math.min(Math.max(point.values[2], -2), 2));
            const diff = sub(point, clamped);
            expect(result.sqrDistance).toBeCloseTo(dot(diff, diff), 10);
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[1].values[i]).toBeCloseTo(
                    clamped.values[i], 12);
            }
        }
    });
});
