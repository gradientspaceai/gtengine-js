import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistLine3CanonicalBox3 } from '../src/DistLine3CanonicalBox3.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function box(...extent: number[]): CanonicalBox {
    return CanonicalBox.fromExtent(v(...extent));
}

// The exact squared distance from a point to the solid canonical box.
function pointBoxSqrDistance(p: Vector, b: CanonicalBox): number {
    let sum = 0;
    for (let i = 0; i < 3; ++i) {
        const e = b.extent.values[i];
        const c = Math.min(Math.max(p.values[i], -e), e);
        sum += (p.values[i] - c) * (p.values[i] - c);
    }
    return sum;
}

// A ternary search for the minimum of the convex function
// t -> pointBoxSqrDistance(origin + t*direction, box).
function sampledMinimum(l: Line, b: CanonicalBox): number {
    let lo = -1e6;
    let hi = 1e6;
    const f = (t: number) =>
        pointBoxSqrDistance(add(l.origin, mul(t, l.direction)), b);
    for (let i = 0; i < 250; ++i) {
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

describe('DistLine3CanonicalBox3', () => {
    const query = new DistLine3CanonicalBox3();
    const unitBox = box(1, 1, 1);

    it('measures an axis-parallel line above the box', () => {
        const result = query.compute(line([0, 0, 4], [1, 0, 0]), unitBox);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 10);
    });

    it('reports zero distance for a line through the box', () => {
        const result = query.compute(line([0, 0, 0], [1, 2, 3]), unitBox);
        expect(result.distance).toBeCloseTo(0, 10);
        expect(result.sqrDistance).toBeCloseTo(0, 10);
    });

    it('handles a degenerate zero direction as a point query', () => {
        const result = query.compute(line([5, -9, 0.5], [0, 0, 0]), unitBox);
        expect(result.parameter).toBe(0);
        expect(result.closest[0].values).toEqual([5, -9, 0.5]);
        expect(result.closest[1].values).toEqual([1, -1, 0.5]);
        expect(result.sqrDistance).toBeCloseTo(16 + 64, 10);
    });

    it('handles a line parallel to an axis but off the box in two axes',
        () => {
            const result = query.compute(line([4, 5, 0], [0, 0, 1]), unitBox);
            expect(result.sqrDistance).toBeCloseTo(9 + 16, 10);
            expect(result.distance).toBeCloseTo(5, 10);
            expect(result.closest[1].values[0]).toBeCloseTo(1, 10);
            expect(result.closest[1].values[1]).toBeCloseTo(1, 10);
        });

    it('measures a diagonal line missing a corner', () => {
        // The line (4,0,0)+t*(-1,1,0) lies in z = 0 and its distance to the
        // unit box is the same as in the 2D case: sqrt(2).
        const result = query.compute(line([4, 0, 0], [-1, 1, 0]), unitBox);
        expect(result.distance).toBeCloseTo(Math.SQRT2, 8);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 8);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 8);
    });

    it('gives the same distance for a line and its reversed direction', () => {
        const b = box(1, 2, 0.5);
        const a = query.compute(line([-3, 4, 5], [1, -2, 3]), b);
        const c = query.compute(line([-3, 4, 5], [-1, 2, -3]), b);
        expect(c.distance).toBeCloseTo(a.distance, 9);
        for (let i = 0; i < 3; ++i) {
            expect(c.closest[1].values[i]).toBeCloseTo(a.closest[1].values[i],
                7);
        }
    });

    it('matches a numeric minimization over the line', () => {
        let seed = 13579;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 8 - 4;
        };
        for (let trial = 0; trial < 200; ++trial) {
            const b = box(Math.abs(rand()) + 0.1, Math.abs(rand()) + 0.1,
                Math.abs(rand()) + 0.1);
            const dir = [rand(), rand(), rand()];
            if (Math.hypot(dir[0], dir[1], dir[2]) < 0.5) {
                continue;
            }
            const l = line([rand(), rand(), rand()], dir);
            const result = query.compute(l, b);

            expect(result.sqrDistance).toBeCloseTo(sampledMinimum(l, b), 6);
            expect(result.distance).toBeCloseTo(
                Math.sqrt(result.sqrDistance), 10);

            // The line point matches the reported parameter.
            const onLine = add(l.origin, mul(result.parameter, l.direction));
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[0].values[i]).toBeCloseTo(
                    onLine.values[i], 7);
            }

            // The box point lies in the box.
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(result.closest[1].values[i]))
                    .toBeLessThanOrEqual(b.extent.values[i] + 1e-7);
            }

            // The reported pair realizes the reported squared distance.
            const diff = sub(result.closest[0], result.closest[1]);
            expect(dot(diff, diff)).toBeCloseTo(result.sqrDistance, 6);
        }
    });

    it('exercises the axis-aligned direction cases', () => {
        let seed = 24601;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 8 - 4;
        };
        const dirs = [
            [1, 0, 0], [0, 1, 0], [0, 0, 1],
            [1, 1, 0], [1, 0, 1], [0, 1, 1],
            [-1, 0, 0], [0, -1, 0], [0, 0, -1],
            [-1, 1, 0], [1, 0, -1], [0, -1, -1]
        ];
        const b = box(1.25, 0.5, 2);
        for (const dir of dirs) {
            for (let trial = 0; trial < 10; ++trial) {
                const l = line([rand(), rand(), rand()], dir);
                const result = query.compute(l, b);
                expect(result.sqrDistance).toBeCloseTo(sampledMinimum(l, b),
                    6);
                const diff = sub(result.closest[0], result.closest[1]);
                expect(dot(diff, diff)).toBeCloseTo(result.sqrDistance, 6);
            }
        }
    });
});
