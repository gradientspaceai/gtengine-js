import { describe, expect, it } from 'vitest';
import { DistPointSegment } from '../src/DistPointSegment.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

describe('DistPointSegment', () => {
    const query = new DistPointSegment();

    it('uses an interior point when the projection is inside [0,1]', () => {
        const result = query.compute(v(2, 3), seg([0, 0], [4, 0]));
        expect(result.parameter).toBeCloseTo(0.5, 12);
        expect(result.closest[1].values).toEqual([2, 0]);
        expect(result.distance).toBeCloseTo(3, 12);
    });

    it('clamps to the second endpoint', () => {
        const result = query.compute(v(9, 3), seg([0, 0], [4, 0]));
        expect(result.parameter).toBe(1);
        expect(result.closest[1].values).toEqual([4, 0]);
        expect(result.distance).toBeCloseTo(Math.sqrt(25 + 9), 12);
    });

    it('clamps to the first endpoint', () => {
        const result = query.compute(v(-3, 4), seg([0, 0], [4, 0]));
        expect(result.parameter).toBe(0);
        expect(result.closest[1].values).toEqual([0, 0]);
        expect(result.distance).toBeCloseTo(5, 12);
    });

    it('handles a degenerate (zero-length) segment', () => {
        const result = query.compute(v(3, 4), seg([1, 1], [1, 1]));
        // Dot(direction, point - p1) = 0 >= 0, so the query reports t = 1.
        expect(result.parameter).toBe(1);
        expect(result.closest[1].values).toEqual([1, 1]);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 9), 12);
    });

    it('reports zero distance for a point on the segment', () => {
        const s = seg([1, 2, 3], [4, 0, -1]);
        const d = sub(s.p[1], s.p[0]);
        const point = add(s.p[0], mul(0.25, d));
        const result = query.compute(point, s);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter).toBeCloseTo(0.25, 12);
    });

    it('agrees with a sampled minimum over the segment', () => {
        let seed = 4242;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const point = v(rand(), rand(), rand());
            const s = seg([rand(), rand(), rand()], [rand(), rand(), rand()]);
            const result = query.compute(point, s);

            expect(result.parameter).toBeGreaterThanOrEqual(0);
            expect(result.parameter).toBeLessThanOrEqual(1);

            const dir = sub(s.p[1], s.p[0]);
            const onSeg = add(s.p[0], mul(result.parameter, dir));
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[1].values[i]).toBeCloseTo(
                    onSeg.values[i], 9);
            }

            let best = Number.MAX_VALUE;
            for (let k = 0; k <= 2000; ++k) {
                const t = k / 2000;
                const q = add(s.p[0], mul(t, dir));
                const diff = sub(point, q);
                best = Math.min(best, dot(diff, diff));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
        }
    });
});
