import { describe, expect, it } from 'vitest';
import { DistPointLine } from '../src/DistPointLine.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

describe('DistPointLine', () => {
    const query = new DistPointLine();

    it('computes the perpendicular distance to an axis line', () => {
        const result = query.compute(v(0, 2), line([0, 0], [1, 0]));
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.sqrDistance).toBeCloseTo(4, 12);
        expect(result.parameter).toBeCloseTo(0, 12);
        expect(result.closest[1].values).toEqual([0, 0]);
    });

    it('reports the input point in closest[0]', () => {
        const point = v(3, -1, 7);
        const result = query.compute(point, line([0, 0, 0], [0, 0, 1]));
        expect(result.closest[0].values).toEqual([3, -1, 7]);
        expect(result.parameter).toBeCloseTo(7, 12);
        expect(result.distance).toBeCloseTo(Math.sqrt(10), 12);
    });

    it('handles a non-unit-length direction (parameter scales)', () => {
        // The line is (1,1) + t*(2,0). The closest point to (5,4) is (5,1),
        // reached at t = 2.
        const result = query.compute(v(5, 4), line([1, 1], [2, 0]));
        expect(result.parameter).toBeCloseTo(2, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(5, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 12);
        expect(result.distance).toBeCloseTo(3, 12);
    });

    it('reports zero distance for a point on the line', () => {
        const l = line([1, 2, 3], [1, -1, 2]);
        const point = add(l.origin, mul(-0.75, l.direction));
        const result = query.compute(point, l);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter).toBeCloseTo(-0.75, 12);
    });

    it('is translation invariant', () => {
        const point = v(2, -3, 5);
        const l = line([1, 0, -1], [0.3, 0.4, -0.5]);
        const shift = v(10, -20, 30);
        const shifted = Line.fromOriginDirection(add(l.origin, shift),
            l.direction);
        const r0 = query.compute(point, l);
        const r1 = query.compute(add(point, shift), shifted);
        expect(r1.distance).toBeCloseTo(r0.distance, 12);
        expect(r1.parameter).toBeCloseTo(r0.parameter, 12);
    });

    it('produces a closest point whose offset is orthogonal to the line',
        () => {
            let seed = 12345;
            const rand = () => {
                seed = (seed * 1103515245 + 12345) % 2147483648;
                return seed / 2147483648 * 4 - 2;
            };
            for (let trial = 0; trial < 50; ++trial) {
                const point = v(rand(), rand(), rand());
                const l = line([rand(), rand(), rand()],
                    [rand() + 3, rand(), rand()]);
                const result = query.compute(point, l);
                const diff = sub(result.closest[0], result.closest[1]);
                expect(dot(diff, l.direction)).toBeCloseTo(0, 9);
                expect(result.sqrDistance).toBeCloseTo(dot(diff, diff), 9);
                expect(result.distance).toBeCloseTo(
                    Math.sqrt(result.sqrDistance), 12);

                // The closest point lies on the line.
                const onLine = add(l.origin, mul(result.parameter,
                    l.direction));
                expect(result.closest[1].values[0]).toBeCloseTo(
                    onLine.values[0], 9);
                expect(result.closest[1].values[1]).toBeCloseTo(
                    onLine.values[1], 9);
                expect(result.closest[1].values[2]).toBeCloseTo(
                    onLine.values[2], 9);

                // No nearby line point is closer.
                for (const dt of [-0.1, -0.01, 0.01, 0.1]) {
                    const other = add(l.origin,
                        mul(result.parameter + dt, l.direction));
                    const od = sub(point, other);
                    expect(dot(od, od)).toBeGreaterThanOrEqual(
                        result.sqrDistance - 1e-12);
                }
            }
        });
});
