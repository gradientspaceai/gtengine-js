import { describe, expect, it } from 'vitest';
import { DistLineLine } from '../src/DistLineLine.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

describe('DistLineLine', () => {
    const query = new DistLineLine();

    it('computes the distance between perpendicular skew lines', () => {
        // The x-axis and the line (0,0,5)+s*(0,1,0) are skew with distance 5.
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            line([0, 0, 5], [0, 1, 0]));
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(0, 12);
    });

    it('computes zero distance for intersecting lines', () => {
        const result = query.compute(line([0, 0, 0], [1, 1, 0]),
            line([2, 0, 0], [-1, 1, 0]));
        expect(result.distance).toBeCloseTo(0, 10);
        expect(result.closest[0].values[0]).toBeCloseTo(1, 10);
        expect(result.closest[0].values[1]).toBeCloseTo(1, 10);
    });

    it('handles parallel lines by choosing s1 = 0', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            line([3, 4, 0], [2, 0, 0]));
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBe(0);
        expect(result.closest[1].values).toEqual([3, 4, 0]);
        expect(result.closest[0].values[0]).toBeCloseTo(3, 12);
    });

    it('handles coincident lines (zero distance, arbitrary pair)', () => {
        const result = query.compute(line([1, 2, 3], [1, 1, 1]),
            line([2, 3, 4], [2, 2, 2]));
        expect(result.distance).toBeCloseTo(0, 10);
    });

    it('is symmetric in the distance', () => {
        const l0 = line([1, -2, 3], [0.5, 1.5, -2]);
        const l1 = line([-4, 0, 1], [2, -1, 0.25]);
        const a = query.compute(l0, l1);
        const b = query.compute(l1, l0);
        expect(b.distance).toBeCloseTo(a.distance, 10);
        expect(b.parameter[0]).toBeCloseTo(a.parameter[1], 10);
        expect(b.parameter[1]).toBeCloseTo(a.parameter[0], 10);
    });

    it('produces a connecting segment orthogonal to both directions', () => {
        let seed = 5150;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const l0 = line([rand(), rand(), rand()],
                [rand() + 4, rand(), rand()]);
            const l1 = line([rand(), rand(), rand()],
                [rand(), rand() + 4, rand()]);
            const result = query.compute(l0, l1);

            const diff = sub(result.closest[0], result.closest[1]);
            expect(dot(diff, l0.direction)).toBeCloseTo(0, 8);
            expect(dot(diff, l1.direction)).toBeCloseTo(0, 8);
            expect(result.sqrDistance).toBeCloseTo(dot(diff, diff), 9);

            // The closest points lie on their lines.
            const c0 = add(l0.origin, mul(result.parameter[0], l0.direction));
            const c1 = add(l1.origin, mul(result.parameter[1], l1.direction));
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[0].values[i]).toBeCloseTo(c0.values[i],
                    9);
                expect(result.closest[1].values[i]).toBeCloseTo(c1.values[i],
                    9);
            }

            // Perturbations do not reduce the squared distance.
            for (const ds of [-0.05, 0.05]) {
                const a = add(l0.origin,
                    mul(result.parameter[0] + ds, l0.direction));
                const d = sub(a, result.closest[1]);
                expect(dot(d, d)).toBeGreaterThanOrEqual(
                    result.sqrDistance - 1e-9);
            }
        }
    });
});
