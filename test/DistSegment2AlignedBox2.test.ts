import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { DistSegment2AlignedBox2 } from '../src/DistSegment2AlignedBox2';
import { Segment } from '../src/Segment';
import { Vector, add, length, mul, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(v(...min), v(...max));
}

// The exact squared distance from a point to the solid aligned box, computed
// independently of the library.
function pointBoxSqrDistance(p: Vector, b: AlignedBox): number {
    let sqrDistance = 0;
    for (let i = 0; i < p.size; ++i) {
        const delta = Math.max(0, b.min.values[i] - p.values[i],
            p.values[i] - b.max.values[i]);
        sqrDistance += delta * delta;
    }
    return sqrDistance;
}

// Verify that the reported box point is inside the box.
function verifyBoxPoint(b: AlignedBox, q: Vector): void {
    for (let i = 0; i < q.size; ++i) {
        expect(q.values[i]).toBeGreaterThanOrEqual(b.min.values[i] - 1e-9);
        expect(q.values[i]).toBeLessThanOrEqual(b.max.values[i] + 1e-9);
    }
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

// The distance from a segment point to the solid box is a convex function of
// the segment parameter, so a ternary search over [0,1] converges to the true
// minimum independently of the query under test.
function bruteForce(s: Segment, b: AlignedBox): number {
    const direction = sub(s.p[1], s.p[0]);
    const f = (t: number): number =>
        pointBoxSqrDistance(add(s.p[0], mul(t, direction)), b);
    let lo = 0;
    let hi = 1;
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
    return Math.sqrt(f(0.5 * (lo + hi)));
}

// Verify that the reported closest points are consistent with the reported
// distance and lie on their primitives.
function verifyClosest(s: Segment, b: AlignedBox,
    result: { distance: number, parameter: number, closest: [Vector, Vector] }
): void {
    expect(result.parameter).toBeGreaterThanOrEqual(0);
    expect(result.parameter).toBeLessThanOrEqual(1);
    const onSeg = add(s.p[0], mul(result.parameter, sub(s.p[1], s.p[0])));
    for (let i = 0; i < onSeg.size; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onSeg.values[i], 9);
    }
    verifyBoxPoint(b, result.closest[1]);
    expect(length(sub(result.closest[0], result.closest[1])))
        .toBeCloseTo(result.distance, 9);
}

describe('DistSegment2AlignedBox2', () => {
    const query = new DistSegment2AlignedBox2();
    const unitBox = box([-1, -1], [1, 1]);

    it('clamps to the first endpoint', () => {
        const s = segment([3, 0], [5, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.parameter).toBe(0);
        expect(result.closest[0].values).toEqual([3, 0]);
        expect(result.closest[1].values).toEqual([1, 0]);
    });

    it('clamps to the second endpoint', () => {
        const s = segment([5, 0], [3, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.parameter).toBe(1);
        expect(result.closest[0].values).toEqual([3, 0]);
    });

    it('uses an interior parameter when the projection is inside', () => {
        const s = segment([-3, 3], [3, 3]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.parameter).toBeGreaterThan(0);
        expect(result.parameter).toBeLessThan(1);
        verifyClosest(s, unitBox, result);
    });

    it('reports zero distance when the segment crosses the box', () => {
        const s = segment([-3, 0], [3, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('handles a degenerate zero-length segment', () => {
        const s = segment([4, 0], [4, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[0].values).toEqual([4, 0]);
        expect(result.closest[1].values).toEqual([1, 0]);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 19283746;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 200; ++trial) {
            const b = box([-1 - rand(), -1 - rand()], [1 + rand(), 1 + rand()]);
            const s = segment([8 * rand() - 4, 8 * rand() - 4],
                [8 * rand() - 4, 8 * rand() - 4]);
            const result = query.compute(s, b);
            const brute = bruteForce(s, b);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-6);
            expect(result.distance).toBeCloseTo(brute, 6);
            verifyClosest(s, b, result);
        }
    });
});
