import { describe, expect, it } from 'vitest';
import { DistPointRectangle } from '../src/DistPointRectangle';
import { DistSegment3Rectangle3 } from '../src/DistSegment3Rectangle3';
import type { DistSegment3Rectangle3Result }
    from '../src/DistSegment3Rectangle3';
import { Rectangle } from '../src/Rectangle';
import { Segment } from '../src/Segment';
import { Vector, add, dot, mul, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

// The rectangle in the z = 0 plane with the standard x and y axes.
function xyRectangle(e0: number, e1: number, center: number[] = [0, 0, 0]):
    Rectangle {
    return Rectangle.fromCenterAxisExtent(v(...center),
        [v(1, 0, 0), v(0, 1, 0)], v(e0, e1));
}

// An orthonormal frame; the first two vectors are the rectangle axes.
function frame(a: number, b: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    return [
        v(ca, sa, 0),
        v(-sa * cb, ca * cb, sb),
        v(sa * sb, -ca * sb, cb)
    ];
}

// The exact squared distance from a point to a solid rectangle: clamp the
// rectangle coordinates of the point to the extents.
function pointRectangleSqrDistance(p: Vector, r: Rectangle): number {
    const delta = sub(p, r.center);
    let closest = r.center.clone();
    for (let i = 0; i < 2; ++i) {
        const e = r.extent.values[i];
        const s = Math.min(Math.max(dot(r.axis[i], delta), -e), e);
        closest = add(closest, mul(s, r.axis[i]));
    }
    const d = sub(p, closest);
    return dot(d, d);
}

// The squared distance from segment(t) to the solid rectangle is a convex
// function of t on [0,1], so a ternary search finds its minimum. A dense
// sampling is used as a second opinion.
function bruteForceSqrDistance(s: Segment, r: Rectangle): number {
    const f = (u: number) => pointRectangleSqrDistance(
        add(s.p[0], mul(u, sub(s.p[1], s.p[0]))), r);

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
    let best = f(0.5 * (lo + hi));
    for (let i = 0; i <= 2000; ++i) {
        const value = f(i / 2000);
        if (value < best) {
            best = value;
        }
    }
    return best;
}

// Verify the internal consistency of a result: the closest points lie on
// their primitives, the reported W-coordinates describe closest[1], and the
// pair realizes the reported distance.
function expectConsistent(result: DistSegment3Rectangle3Result,
    s: Segment, r: Rectangle): void {
    expect(result.parameter).toBeGreaterThanOrEqual(-1e-12);
    expect(result.parameter).toBeLessThanOrEqual(1 + 1e-12);

    const onSegment = add(s.p[0],
        mul(result.parameter, sub(s.p[1], s.p[0])));
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onSegment.values[i], 6);
    }

    let onRectangle = r.center.clone();
    for (let i = 0; i < 2; ++i) {
        expect(Math.abs(result.cartesian[i]))
            .toBeLessThanOrEqual(r.extent.values[i] + 1e-9);
        onRectangle = add(onRectangle,
            mul(result.cartesian[i], r.axis[i]));
    }
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[1].values[i])
            .toBeCloseTo(onRectangle.values[i], 6);
    }

    const diff = sub(result.closest[0], result.closest[1]);
    expect(Math.sqrt(dot(diff, diff))).toBeCloseTo(result.distance, 6);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 8);
}

// A small deterministic linear congruential generator.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistSegment3Rectangle3', () => {
    const query = new DistSegment3Rectangle3();
    const rect = xyRectangle(2, 1);

    it('measures a segment parallel to the rectangle plane', () => {
        const s = segment([-1, 0, 3], [1, 0, 3]);
        const result = query.compute(s, rect);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(0, 12);
        expectConsistent(result, s, rect);
    });

    it('reports zero distance when the segment crosses the rectangle', () => {
        const s = segment([0.5, 0.25, -2], [0.5, 0.25, 2]);
        const result = query.compute(s, rect);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter).toBeCloseTo(0.5, 10);
        expect(result.cartesian[0]).toBeCloseTo(0.5, 10);
        expect(result.cartesian[1]).toBeCloseTo(0.25, 10);
        expectConsistent(result, s, rect);
    });

    it('reports zero distance when an endpoint lies on the rectangle', () => {
        const s = segment([1, 0.5, 0], [1, 0.5, 4]);
        const result = query.compute(s, rect);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, rect);
    });

    it('misses the rectangle when the crossing point is outside it', () => {
        // The line crosses z = 0 at (5,0,0), well beyond the x extent 2, so
        // the closest rectangle point is on the edge x = 2.
        const s = segment([5, 0, -1], [5, 0, 1]);
        const result = query.compute(s, rect);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.cartesian[0]).toBeCloseTo(2, 10);
        expect(result.cartesian[1]).toBeCloseTo(0, 10);
        expectConsistent(result, s, rect);
    });

    it('clamps to the first endpoint when the line minimum is behind it',
        () => {
            const s = segment([6, 0, 0], [9, 0, 0]);
            const result = query.compute(s, rect);
            expect(result.parameter).toBe(0);
            expect(result.distance).toBeCloseTo(4, 10);
            expect(result.cartesian[0]).toBeCloseTo(2, 10);
            expectConsistent(result, s, rect);
        });

    it('clamps to the second endpoint when the line minimum is beyond it',
        () => {
            const s = segment([9, 0, 0], [6, 0, 0]);
            const result = query.compute(s, rect);
            expect(result.parameter).toBe(1);
            expect(result.distance).toBeCloseTo(4, 10);
            expectConsistent(result, s, rect);
        });

    it('measures the distance to a rectangle corner', () => {
        // The nearest rectangle point to the segment is the corner (2,1,0).
        const s = segment([5, 4, 0], [7, 6, 0]);
        const result = query.compute(s, rect);
        expect(result.parameter).toBe(0);
        expect(result.distance).toBeCloseTo(Math.hypot(3, 3), 10);
        expect(result.cartesian[0]).toBeCloseTo(2, 10);
        expect(result.cartesian[1]).toBeCloseTo(1, 10);
        expectConsistent(result, s, rect);
    });

    it('gives the same distance for both segment orientations', () => {
        const axes = frame(0.6, 1.3);
        const r = Rectangle.fromCenterAxisExtent(v(1, -1, 2),
            [axes[0], axes[1]], v(1.5, 0.75));
        const cases: Array<[number[], number[]]> = [
            [[5, 0, 0], [8, 0, 0]],
            [[-1, 0, 3], [1, 0, 3]],
            [[-3, -3, -3], [3, 3, 3]],
            [[2, 2, -5], [2, 2, 5]],
            [[0.5, 0.5, 0.5], [7, 1, -3]]
        ];
        for (const [p0, p1] of cases) {
            const forward = query.compute(segment(p0, p1), r);
            const backward = query.compute(segment(p1, p0), r);
            expect(forward.distance).toBeCloseTo(backward.distance, 9);
        }
    });

    it('handles a degenerate zero-length segment off the rectangle', () => {
        const s = segment([3, 4, 5], [3, 4, 5]);
        const result = query.compute(s, rect);
        const expected = new DistPointRectangle().compute(v(3, 4, 5), rect);
        expect(result.distance).toBeCloseTo(expected.distance, 10);
        expect(result.closest[0].values[0]).toBeCloseTo(3, 12);
        expectConsistent(result, s, rect);
    });

    it('handles a degenerate zero-length segment on the rectangle', () => {
        const s = segment([0.5, -0.25, 0], [0.5, -0.25, 0]);
        const result = query.compute(s, rect);
        expect(result.distance).toBeCloseTo(0, 10);
        expectConsistent(result, s, rect);
    });

    it('agrees with a brute-force sampling for random configurations', () => {
        const random = makeRandom(87084);
        for (let trial = 0; trial < 250; ++trial) {
            const axes = frame(2 * Math.PI * random(), Math.PI * random());
            const r = Rectangle.fromCenterAxisExtent(
                v(3 * random() - 1.5, 3 * random() - 1.5, 3 * random() - 1.5),
                [axes[0], axes[1]],
                v(0.2 + 2 * random(), 0.2 + 2 * random()));
            const s = segment(
                [6 * random() - 3, 6 * random() - 3, 6 * random() - 3],
                [6 * random() - 3, 6 * random() - 3, 6 * random() - 3]);
            const result = query.compute(s, r);
            const expected = Math.sqrt(bruteForceSqrDistance(s, r));
            expect(Math.abs(result.distance - expected)).toBeLessThan(1e-6);
            expectConsistent(result, s, r);
        }
    });
});
