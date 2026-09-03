import { describe, expect, it } from 'vitest';
import { DistPointRay } from '../src/DistPointRay.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

describe('DistPointRay', () => {
    const query = new DistPointRay();

    it('uses an interior ray point when the projection is positive', () => {
        const result = query.compute(v(3, 4), ray([0, 0], [1, 0]));
        expect(result.parameter).toBeCloseTo(3, 12);
        expect(result.closest[1].values).toEqual([3, 0]);
        expect(result.distance).toBeCloseTo(4, 12);
    });

    it('clamps to the ray origin when the projection is nonpositive', () => {
        const result = query.compute(v(-3, 4), ray([0, 0], [1, 0]));
        expect(result.parameter).toBe(0);
        expect(result.closest[1].values).toEqual([0, 0]);
        expect(result.distance).toBeCloseTo(5, 12);
    });

    it('clamps at exactly the origin projection', () => {
        const result = query.compute(v(0, 2), ray([0, 0], [1, 0]));
        expect(result.parameter).toBe(0);
        expect(result.distance).toBeCloseTo(2, 12);
    });

    it('handles a non-unit direction', () => {
        // Ray (1,1) + t*(0,3). The point (4,7) projects to t = 2.
        const result = query.compute(v(4, 7), ray([1, 1], [0, 3]));
        expect(result.parameter).toBeCloseTo(2, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(7, 12);
        expect(result.distance).toBeCloseTo(3, 12);
    });

    it('reports zero distance for a point on the ray', () => {
        const r = ray([1, 2, 3], [1, -1, 2]);
        const point = add(r.origin, mul(1.5, r.direction));
        const result = query.compute(point, r);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter).toBeCloseTo(1.5, 12);
    });

    it('matches a sampled minimum over the ray', () => {
        let seed = 987;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const point = v(rand(), rand(), rand());
            const r = ray([rand(), rand(), rand()],
                [rand(), rand() + 2, rand()]);
            const result = query.compute(point, r);

            expect(result.parameter).toBeGreaterThanOrEqual(0);
            const onRay = add(r.origin, mul(result.parameter, r.direction));
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[1].values[i]).toBeCloseTo(
                    onRay.values[i], 9);
            }

            let best = Number.MAX_VALUE;
            for (let k = 0; k <= 4000; ++k) {
                const t = k * 0.005;
                const q = add(r.origin, mul(t, r.direction));
                const d = sub(point, q);
                best = Math.min(best, dot(d, d));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
        }
    });
});
