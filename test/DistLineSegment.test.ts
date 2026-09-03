import { describe, expect, it } from 'vitest';
import { DistLineSegment } from '../src/DistLineSegment.js';
import { Line } from '../src/Line.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

describe('DistLineSegment', () => {
    const query = new DistLineSegment();

    it('uses interior points when the segment parameter is in [0,1]', () => {
        // The x-axis and the segment from (2,3,4) to (2,3,-4). The closest
        // segment point is (2,3,0) at s1 = 0.5, and the distance is 3.
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            seg([2, 3, 4], [2, 3, -4]));
        expect(result.parameter[1]).toBeCloseTo(0.5, 10);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[0].values[0]).toBeCloseTo(2, 10);
    });

    it('clamps to the endpoint Q1', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            seg([0, 3, 5], [0, 3, 1]));
        expect(result.parameter[1]).toBe(1);
        expect(result.closest[1].values).toEqual([0, 3, 1]);
        expect(result.distance).toBeCloseTo(Math.sqrt(9 + 1), 10);
    });

    it('clamps to the endpoint Q0', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            seg([0, 3, 1], [0, 3, 5]));
        expect(result.parameter[1]).toBe(0);
        expect(result.closest[1].values).toEqual([0, 3, 1]);
        expect(result.distance).toBeCloseTo(Math.sqrt(9 + 1), 10);
    });

    it('handles a line parallel to the segment', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            seg([2, 4, 0], [7, 4, 0]));
        expect(result.parameter[1]).toBe(0);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values).toEqual([2, 4, 0]);
    });

    it('handles a degenerate (zero-length) segment', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            seg([2, 5, 0], [2, 5, 0]));
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values).toEqual([2, 5, 0]);
    });

    it('reports zero distance when the segment crosses the line', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            seg([3, -1, 0], [3, 1, 0]));
        expect(result.distance).toBeCloseTo(0, 10);
        expect(result.parameter[1]).toBeCloseTo(0.5, 10);
    });

    it('agrees with a sampled minimum over the segment', () => {
        let seed = 8080;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const l = line([rand(), rand(), rand()],
                [rand() + 4, rand(), rand()]);
            const s = seg([rand(), rand(), rand()],
                [rand(), rand(), rand()]);
            const result = query.compute(l, s);

            expect(result.parameter[1]).toBeGreaterThanOrEqual(0);
            expect(result.parameter[1]).toBeLessThanOrEqual(1);

            const segDir = sub(s.p[1], s.p[0]);
            const a00 = dot(l.direction, l.direction);
            let best = Number.MAX_VALUE;
            for (let k = 0; k <= 2000; ++k) {
                const s1 = k / 2000;
                const q = add(s.p[0], mul(s1, segDir));
                const w = sub(q, l.origin);
                const t = dot(l.direction, w) / a00;
                const d = sub(w, mul(t, l.direction));
                best = Math.min(best, dot(d, d));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
        }
    });
});
