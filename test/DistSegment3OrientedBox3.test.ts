import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistPointOrientedBox } from '../src/DistPointOrientedBox.js';
import { DistSegment3CanonicalBox3 } from '../src/DistSegment3CanonicalBox3.js';
import { DistSegment3OrientedBox3 } from '../src/DistSegment3OrientedBox3.js';
import type { DistSegment3OrientedBox3Result }
    from '../src/DistSegment3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

// An axis-aligned oriented box (identity axes) with the given center and
// extents.
function alignedObb(center: number[], extent: number[]): OrientedBox {
    return OrientedBox.fromCenterAxisExtent(v(...center),
        [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], v(...extent));
}

// A rotation about the z axis by angle a, then about the x axis by angle b,
// applied to the standard basis; the result is an orthonormal frame.
function frame(a: number, b: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    // Rz(a) columns, then Rx(b) applied to them.
    const u0 = v(ca, sa, 0);
    const u1 = v(-sa * cb, ca * cb, sb);
    const u2 = v(sa * sb, -ca * sb, cb);
    return [u0, u1, u2];
}

// The exact squared distance from a point to a solid oriented box: clamp the
// box coordinates of the point to the extents.
function pointBoxSqrDistance(p: Vector, b: OrientedBox): number {
    const delta = sub(p, b.center);
    let closest = b.center.clone();
    for (let i = 0; i < 3; ++i) {
        const e = b.extent.values[i];
        const y = Math.min(Math.max(dot(b.axis[i], delta), -e), e);
        closest = add(closest, mul(y, b.axis[i]));
    }
    const d = sub(p, closest);
    return dot(d, d);
}

// The squared distance from segment(t) to the solid box is a convex function
// of t on [0,1], so a ternary search finds its minimum. A dense sampling is
// used as a second opinion.
function bruteForceSqrDistance(s: Segment, b: OrientedBox): number {
    const f = (u: number) => pointBoxSqrDistance(
        add(s.p[0], mul(u, sub(s.p[1], s.p[0]))), b);

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
// their primitives and realize the reported distance.
function expectConsistent(result: DistSegment3OrientedBox3Result,
    s: Segment, b: OrientedBox): void {
    expect(result.parameter).toBeGreaterThanOrEqual(-1e-12);
    expect(result.parameter).toBeLessThanOrEqual(1 + 1e-12);

    const onSegment = add(s.p[0],
        mul(result.parameter, sub(s.p[1], s.p[0])));
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onSegment.values[i], 6);
    }

    // The closest box point must have |y[i]| <= e[i] in the box frame.
    const delta = sub(result.closest[1], b.center);
    for (let i = 0; i < 3; ++i) {
        expect(Math.abs(dot(b.axis[i], delta)))
            .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
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

describe('DistSegment3OrientedBox3', () => {
    const query = new DistSegment3OrientedBox3();
    const unitBox = alignedObb([0, 0, 0], [1, 1, 1]);

    it('measures a segment parallel to a face', () => {
        const s = segment([-2, 0, 4], [2, 0, 4]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 10);
        expectConsistent(result, s, unitBox);
    });

    it('reports zero distance for a segment passing through the box', () => {
        const s = segment([-3, 0, 0], [3, 0, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitBox);
    });

    it('reports zero distance for a segment strictly inside the box', () => {
        const s = segment([-0.5, 0.25, 0.1], [0.5, -0.25, -0.1]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitBox);
    });

    it('reports zero distance when an endpoint touches a face', () => {
        const s = segment([1, 0, 0], [4, 0, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitBox);
    });

    it('clamps to the first endpoint when the line minimum is behind it',
        () => {
            const s = segment([3, 0, 0], [6, 0, 0]);
            const result = query.compute(s, unitBox);
            expect(result.parameter).toBe(0);
            expect(result.distance).toBeCloseTo(2, 10);
            expect(result.closest[1].values[0]).toBeCloseTo(1, 10);
            expectConsistent(result, s, unitBox);
        });

    it('clamps to the second endpoint when the line minimum is beyond it',
        () => {
            const s = segment([6, 0, 0], [3, 0, 0]);
            const result = query.compute(s, unitBox);
            expect(result.parameter).toBe(1);
            expect(result.distance).toBeCloseTo(2, 10);
            expectConsistent(result, s, unitBox);
        });

    it('measures the distance to a box corner', () => {
        const s = segment([3, 3, 3], [4, 5, 6]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(Math.sqrt(12), 10);
        expect(result.parameter).toBe(0);
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i]).toBeCloseTo(1, 10);
        }
        expectConsistent(result, s, unitBox);
    });

    it('measures a rotated box analytically', () => {
        // A box rotated 45 degrees about z, with extents (1,1,1). Its
        // "corner" nearest to +x lies at distance sqrt(2) from the center
        // along x. A segment parallel to the y axis at x = 5 is therefore at
        // distance 5 - sqrt(2).
        const c = Math.SQRT1_2;
        const box45 = OrientedBox.fromCenterAxisExtent(v(0, 0, 0),
            [v(c, c, 0), v(-c, c, 0), v(0, 0, 1)], v(1, 1, 1));
        const s = segment([5, -3, 0], [5, 3, 0]);
        const result = query.compute(s, box45);
        expect(result.distance).toBeCloseTo(5 - Math.SQRT2, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(Math.SQRT2, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 10);
        expectConsistent(result, s, box45);
    });

    it('is translation and rotation invariant', () => {
        const axes = frame(0.7, -1.1);
        const rotated = OrientedBox.fromCenterAxisExtent(v(3, -2, 5), axes,
            v(1, 2, 0.5));
        const cbox = CanonicalBox.fromExtent(v(1, 2, 0.5));
        const cQuery = new DistSegment3CanonicalBox3();
        const random = makeRandom(20260901);
        for (let trial = 0; trial < 60; ++trial) {
            // Build a segment in the box frame, then map it to world space.
            const local: Vector[] = [];
            for (let k = 0; k < 2; ++k) {
                local.push(v(6 * random() - 3, 6 * random() - 3,
                    6 * random() - 3));
            }
            const world = local.map((q) => {
                let w = rotated.center.clone();
                for (let i = 0; i < 3; ++i) {
                    w = add(w, mul(q.values[i], axes[i]));
                }
                return w;
            });
            const expected = cQuery.compute(
                Segment.fromEndpoints(local[0], local[1]), cbox);
            const actual = query.compute(
                Segment.fromEndpoints(world[0], world[1]), rotated);
            expect(actual.distance).toBeCloseTo(expected.distance, 9);
            expect(actual.parameter).toBeCloseTo(expected.parameter, 9);
        }
    });

    it('gives the same distance for both segment orientations', () => {
        const axes = frame(0.35, 0.9);
        const b = OrientedBox.fromCenterAxisExtent(v(1, 1, -1), axes,
            v(2, 0.5, 1));
        const cases: Array<[number[], number[]]> = [
            [[3, 0, 0], [6, 0, 0]],
            [[-2, 0, 4], [2, 0, 4]],
            [[-3, -3, -3], [3, 3, 3]],
            [[2, 2, -5], [2, 2, 5]],
            [[0.5, 0.5, 0.5], [7, 1, -3]]
        ];
        for (const [p0, p1] of cases) {
            const forward = query.compute(segment(p0, p1), b);
            const backward = query.compute(segment(p1, p0), b);
            // Only the distance is orientation independent; when several
            // segment points realize the minimum, the reported parameter
            // depends on the traversal direction.
            expect(forward.distance).toBeCloseTo(backward.distance, 10);
        }
    });

    it('handles a degenerate zero-length segment outside the box', () => {
        const axes = frame(1.2, 0.4);
        const b = OrientedBox.fromCenterAxisExtent(v(0, 1, 0), axes,
            v(1, 2, 3));
        const s = segment([3, 4, 5], [3, 4, 5]);
        const result = query.compute(s, b);
        const expected = new DistPointOrientedBox().compute(v(3, 4, 5), b);
        expect(result.distance).toBeCloseTo(expected.distance, 12);
        expect(result.closest[0].values[0]).toBeCloseTo(3, 12);
        expectConsistent(result, s, b);
    });

    it('handles a degenerate zero-length segment inside the box', () => {
        const s = segment([0.25, -0.5, 0.75], [0.25, -0.5, 0.75]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitBox);
    });

    it('handles a degenerate (flat) box', () => {
        const flat = alignedObb([0, 0, 0], [2, 3, 0]);
        const s = segment([0, 0, 1], [0, 0, 5]);
        const result = query.compute(s, flat);
        expect(result.distance).toBeCloseTo(1, 10);
        expect(result.parameter).toBe(0);
        expectConsistent(result, s, flat);
    });

    it('agrees with a brute-force sampling for random configurations', () => {
        const random = makeRandom(87042);
        for (let trial = 0; trial < 250; ++trial) {
            const axes = frame(2 * Math.PI * random(), Math.PI * random());
            const b = OrientedBox.fromCenterAxisExtent(
                v(4 * random() - 2, 4 * random() - 2, 4 * random() - 2),
                axes,
                v(0.2 + 2 * random(), 0.2 + 2 * random(),
                    0.2 + 2 * random()));
            const s = segment(
                [8 * random() - 4, 8 * random() - 4, 8 * random() - 4],
                [8 * random() - 4, 8 * random() - 4, 8 * random() - 4]);
            const result = query.compute(s, b);
            const expected = Math.sqrt(bruteForceSqrDistance(s, b));
            expect(Math.abs(result.distance - expected)).toBeLessThan(1e-6);
            expectConsistent(result, s, b);
        }
    });
});
