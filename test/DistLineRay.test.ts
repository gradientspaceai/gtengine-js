import { describe, expect, it } from 'vitest';
import { DistLineRay } from '../src/DistLineRay.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

describe('DistLineRay', () => {
    const query = new DistLineRay();

    it('uses interior points when the ray parameter is nonnegative', () => {
        // The x-axis line and the ray (0,3,0)+s*(0,-1,1). A ray point is
        // (0,3-s,s), whose distance to the x-axis is sqrt((3-s)^2+s^2). The
        // minimum is at s = 1.5, giving the point (0,1.5,1.5) and the
        // distance sqrt(4.5).
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            ray([0, 3, 0], [0, -1, 1]));
        expect(result.parameter[1]).toBeCloseTo(1.5, 10);
        expect(result.distance).toBeCloseTo(Math.sqrt(4.5), 10);
        expect(result.closest[1].values[1]).toBeCloseTo(1.5, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(1.5, 10);
    });

    it('clamps to the ray origin when the unconstrained s1 is negative',
        () => {
            // The ray points away from the line.
            const result = query.compute(line([0, 0, 0], [1, 0, 0]),
                ray([0, 3, 0], [0, 1, 0]));
            expect(result.parameter[1]).toBe(0);
            expect(result.distance).toBeCloseTo(3, 12);
            expect(result.closest[1].values).toEqual([0, 3, 0]);
        });

    it('handles a line and ray that are parallel', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            ray([5, 4, 0], [2, 0, 0]));
        expect(result.parameter[1]).toBe(0);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values).toEqual([5, 4, 0]);
    });

    it('reports zero distance when the ray meets the line', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            ray([2, 3, 0], [0, -1, 0]));
        expect(result.distance).toBeCloseTo(0, 10);
        expect(result.parameter[1]).toBeCloseTo(3, 10);
    });

    it('agrees with a sampled minimum', () => {
        let seed = 606;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 40; ++trial) {
            const l = line([rand(), rand(), rand()],
                [rand() + 4, rand(), rand()]);
            const r = ray([rand(), rand(), rand()],
                [rand(), rand() + 4, rand()]);
            const result = query.compute(l, r);

            expect(result.parameter[1]).toBeGreaterThanOrEqual(0);

            // For each sampled ray point, the exact line distance is at least
            // the reported distance.
            const a00 = dot(l.direction, l.direction);
            let best = Number.MAX_VALUE;
            for (let k = 0; k <= 1500; ++k) {
                const s1 = k * 0.01;
                const q = add(r.origin, mul(s1, r.direction));
                const w = sub(q, l.origin);
                const t = dot(l.direction, w) / a00;
                const d = sub(w, mul(t, l.direction));
                best = Math.min(best, dot(d, d));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
        }
    });
});
