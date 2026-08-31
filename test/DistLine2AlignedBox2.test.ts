import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { DistLine2AlignedBox2 } from '../src/DistLine2AlignedBox2';
import { Line } from '../src/Line';
import { Vector, add, dot, mul, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(v(...min), v(...max));
}

// The exact squared distance from a point to a solid aligned box.
function pointBoxSqrDistance(p: Vector, b: AlignedBox): number {
    let sum = 0;
    for (let i = 0; i < 2; ++i) {
        const c = Math.min(Math.max(p.values[i], b.min.values[i]),
            b.max.values[i]);
        sum += (p.values[i] - c) * (p.values[i] - c);
    }
    return sum;
}

// A ternary search for the minimum of the convex function
// t -> pointBoxSqrDistance(origin + t*direction, box).
function sampledMinimum(l: Line, b: AlignedBox): number {
    let lo = -1e6;
    let hi = 1e6;
    const f = (t: number) =>
        pointBoxSqrDistance(add(l.origin, mul(t, l.direction)), b);
    for (let i = 0; i < 200; ++i) {
        const m0 = lo + (hi - lo) / 3;
        const m1 = hi - (hi - lo) / 3;
        if (f(m0) < f(m1)) {
            hi = m1;
        }
        else {
            lo = m0;
        }
    }
    return f(0.5 * (lo + hi));
}

describe('DistLine2AlignedBox2', () => {
    const query = new DistLine2AlignedBox2();

    it('measures a horizontal line above the box', () => {
        const result = query.compute(line([0, 5], [1, 0]),
            box([-1, -1], [1, 1]));
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 10);
    });

    it('measures a vertical line to the right of the box', () => {
        const result = query.compute(line([7, 0], [0, 1]),
            box([-1, -1], [1, 1]));
        expect(result.distance).toBeCloseTo(6, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 10);
        expect(result.closest[0].values[0]).toBeCloseTo(7, 10);
    });

    it('reports zero distance for a line through the box', () => {
        const result = query.compute(line([0, 0], [1, 1]),
            box([-1, -1], [1, 1]));
        expect(result.distance).toBeCloseTo(0, 10);
    });

    it('measures a diagonal line missing a corner', () => {
        // The line x + y = 4, i.e. (4,0) + t*(-1,1), and the unit box. The
        // closest box point is the corner (1,1), at distance
        // (4-2)/sqrt(2) = sqrt(2).
        const result = query.compute(line([4, 0], [-1, 1]),
            box([-1, -1], [1, 1]));
        expect(result.distance).toBeCloseTo(Math.SQRT2, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 10);
    });

    it('handles a degenerate zero direction as a point query', () => {
        const result = query.compute(line([5, 9], [0, 0]),
            box([-1, -1], [1, 1]));
        expect(result.parameter).toBe(0);
        expect(result.closest[0].values).toEqual([5, 9]);
        expect(result.closest[1].values).toEqual([1, 1]);
        expect(result.sqrDistance).toBeCloseTo(16 + 64, 10);
    });

    it('handles an off-center box', () => {
        const b = box([10, 20], [12, 26]);
        // The infinite line y = 23 crosses the box, so the distance is 0.
        expect(query.compute(line([0, 23], [1, 0]), b).distance)
            .toBeCloseTo(0, 10);
        // The line y = 30 is above the box.
        const result = query.compute(line([0, 30], [1, 0]), b);
        expect(result.distance).toBeCloseTo(4, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(26, 10);
    });

    it('gives the same answer for a line and its reversed direction', () => {
        // The line x + y = 12 misses the box, whose far corner is (3,4). The
        // distance is (12-7)/sqrt(2).
        const b = box([-1, -2], [3, 4]);
        const base = line([12, 0], [-1, 1]);
        const flipped = line([12, 0], [1, -1]);
        const a = query.compute(base, b);
        const c = query.compute(flipped, b);
        expect(a.distance).toBeCloseTo(5 / Math.SQRT2, 10);
        expect(c.distance).toBeCloseTo(a.distance, 10);
        expect(a.closest[1].values).toEqual([3, 4]);
        expect(c.closest[1].values[0]).toBeCloseTo(3, 8);
        expect(c.closest[1].values[1]).toBeCloseTo(4, 8);
    });

    it('matches a numeric minimization and reports consistent points', () => {
        let seed = 90210;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 10 - 5;
        };
        for (let trial = 0; trial < 150; ++trial) {
            const lo = [rand(), rand()];
            const b = box(lo, [lo[0] + Math.abs(rand()) + 0.1,
                lo[1] + Math.abs(rand()) + 0.1]);
            const dir = [rand(), rand()];
            if (Math.hypot(dir[0], dir[1]) < 0.25) {
                continue;
            }
            const l = line([rand(), rand()], dir);
            const result = query.compute(l, b);

            expect(result.sqrDistance).toBeCloseTo(sampledMinimum(l, b), 7);

            // The line point matches the reported parameter.
            const onLine = add(l.origin, mul(result.parameter, l.direction));
            expect(result.closest[0].values[0]).toBeCloseTo(onLine.values[0],
                7);
            expect(result.closest[0].values[1]).toBeCloseTo(onLine.values[1],
                7);

            // The box point lies in the box.
            for (let i = 0; i < 2; ++i) {
                expect(result.closest[1].values[i]).toBeGreaterThanOrEqual(
                    b.min.values[i] - 1e-8);
                expect(result.closest[1].values[i]).toBeLessThanOrEqual(
                    b.max.values[i] + 1e-8);
            }

            const diff = sub(result.closest[0], result.closest[1]);
            expect(dot(diff, diff)).toBeCloseTo(result.sqrDistance, 7);
        }
    });
});
