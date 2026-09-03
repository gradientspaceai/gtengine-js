import { describe, expect, it } from 'vitest';
import { DistSegment2OrientedBox2 } from '../src/DistSegment2OrientedBox2.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(center: number[], axis: number[][],
    extent: number[]): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(v(...center),
        axis.map(a => v(...a)), v(...extent));
}

// The exact squared distance from a point to the solid oriented box, computed
// independently of the library.
function pointBoxSqrDistance(p: Vector, b: OrientedBox): number {
    const delta = sub(p, b.center);
    let sqrDistance = 0;
    for (let i = 0; i < p.size; ++i) {
        const y = dot(delta, b.axis[i]);
        const excess = Math.max(0, Math.abs(y) - b.extent.values[i]);
        sqrDistance += excess * excess;
    }
    return sqrDistance;
}

// Verify that the reported box point is inside the box.
function verifyBoxPoint(b: OrientedBox, q: Vector): void {
    const delta = sub(q, b.center);
    for (let i = 0; i < q.size; ++i) {
        expect(Math.abs(dot(delta, b.axis[i])))
            .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
    }
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

// The distance from a segment point to the solid box is a convex function of
// the segment parameter, so a ternary search over [0,1] converges to the true
// minimum independently of the query under test.
function bruteForce(s: Segment, b: OrientedBox): number {
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
function verifyClosest(s: Segment, b: OrientedBox,
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

describe('DistSegment2OrientedBox2', () => {
    const query = new DistSegment2OrientedBox2();
    const axisAligned = box([0, 0], [[1, 0], [0, 1]], [1, 1]);
    const c = Math.SQRT1_2;
    const rotated = box([0, 0], [[c, c], [-c, c]], [1, 1]);

    it('clamps to the first endpoint', () => {
        const s = segment([3, 0], [5, 0]);
        const result = query.compute(s, axisAligned);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.parameter).toBe(0);
    });

    it('clamps to the second endpoint', () => {
        const s = segment([5, 0], [3, 0]);
        const result = query.compute(s, axisAligned);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.parameter).toBe(1);
    });

    it('handles a box rotated by 45 degrees', () => {
        const s = segment([4, 0], [3, 0]);
        const result = query.compute(s, rotated);
        expect(result.distance).toBeCloseTo(3 - Math.SQRT2, 10);
        verifyClosest(s, rotated, result);
    });

    it('handles a degenerate zero-length segment', () => {
        const s = segment([0, 4], [0, 4]);
        const result = query.compute(s, rotated);
        expect(result.distance).toBeCloseTo(4 - Math.SQRT2, 10);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 66554433;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 200; ++trial) {
            const angle = 2 * Math.PI * rand();
            const ca = Math.cos(angle);
            const sa = Math.sin(angle);
            const b = box([2 * rand() - 1, 2 * rand() - 1],
                [[ca, sa], [-sa, ca]],
                [0.2 + rand(), 0.2 + rand()]);
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
