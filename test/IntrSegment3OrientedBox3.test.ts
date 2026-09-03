import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import {
    IntrSegment3AlignedBox3TI,
    IntrSegment3AlignedBox3FI
} from '../src/IntrSegment3AlignedBox3.js';
import {
    IntrSegment3OrientedBox3TI,
    IntrSegment3OrientedBox3FI,
    defaultIntrSegment3OrientedBox3TIResult,
    defaultIntrSegment3OrientedBox3FIResult
} from '../src/IntrSegment3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

// An orthonormal frame: a rotation about z by 'a' followed by a rotation
// about the rotated x axis by 'b'.
function frame(a: number, b: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    return [
        vec(ca, sa, 0),
        vec(-sa * cb, ca * cb, sb),
        vec(sa * sb, -ca * sb, cb)
    ];
}

function obb(center: number[], axes: Vector[], extent: number[]):
    OrientedBox {
    return OrientedBox.fromCenterAxisExtent(Vector.fromArray(center), axes,
        Vector.fromArray(extent));
}

// The signed "outside" amount of a point relative to the box: negative when
// strictly inside, zero on the boundary, positive when outside.
function boxSignedDepth(p: Vector, box: OrientedBox): number {
    const delta = sub(p, box.center);
    let worst = -Number.MAX_VALUE;
    for (let i = 0; i < 3; ++i) {
        const value = Math.abs(dot(box.axis[i], delta)) -
            box.extent.values[i];
        if (value > worst) {
            worst = value;
        }
    }
    return worst;
}

// The most deeply contained sampled point of the segment, expressed as the
// signed depth above; negative means some sample is strictly inside.
function bestSampledDepth(s: Segment, box: OrientedBox): number {
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= 4000; ++i) {
        const p = add(s.p[0], mul(i / 4000, sub(s.p[1], s.p[0])));
        const value = boxSignedDepth(p, box);
        if (value < best) {
            best = value;
        }
    }
    return best;
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrSegment3OrientedBox3', () => {
    const ti = new IntrSegment3OrientedBox3TI();
    const fi = new IntrSegment3OrientedBox3FI();
    const identityAxes = [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)];
    const unit = obb([0, 0, 0], identityAxes, [1, 1, 1]);

    it('default results report no intersection', () => {
        const tiResult = defaultIntrSegment3OrientedBox3TIResult();
        expect(tiResult.intersect).toBe(false);
        const fiResult = defaultIntrSegment3OrientedBox3FIResult();
        expect(fiResult.intersect).toBe(false);
        expect(fiResult.numIntersections).toBe(0);
        expect(fiResult.parameter).toEqual([0, 0]);
    });

    it('clips a segment crossing an axis-aligned box', () => {
        const s = segment([-3, 0, 0], [3, 0, 0]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The parameters are relative to the centered form, extent 3.
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports the whole segment when it lies inside the box', () => {
        const s = segment([-0.5, -0.25, 0.1], [0.5, 0.25, -0.1]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(-0.5, 12);
        expect(result.point[1].values[0]).toBeCloseTo(0.5, 12);
    });

    it('clips a segment with only one endpoint inside the box', () => {
        // The segment straddles the face x = 1.
        const s = segment([0, 0, 0], [4, 0, 0]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports a tangent (single-point) contact on a face', () => {
        // The segment lies in the plane x = 1 and grazes the box face.
        const s = segment([1, -3, 0], [1, 3, 0]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        for (const p of result.point) {
            expect(p.values[0]).toBeCloseTo(1, 12);
        }
    });

    it('reports a single-point contact at a corner', () => {
        // The segment touches the corner (1,1,1) and immediately leaves.
        const s = segment([1, 1, 1], [3, 4, 5]);
        expect(ti.test(s, unit).intersect).toBe(true);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        for (let i = 0; i < 3; ++i) {
            expect(result.point[0].values[i]).toBeCloseTo(1, 10);
        }
    });

    it('reports no intersection for a segment that stops short of the box',
        () => {
            const s = segment([2, 0, 0], [5, 0, 0]);
            expect(ti.test(s, unit).intersect).toBe(false);
            expect(fi.find(s, unit).intersect).toBe(false);
        });

    it('reports no intersection when the line hits but the segment does not',
        () => {
            // The line through these points passes through the box, but the
            // overlap with the segment t-interval is empty.
            const s = segment([3, 0, 0], [6, 0, 0]);
            expect(ti.test(s, unit).intersect).toBe(false);
            expect(fi.find(s, unit).intersect).toBe(false);
        });

    it('is sensitive to the box orientation', () => {
        // A box rotated 45 degrees about z reaches out to x = sqrt(2), so a
        // point at x = 1.2 on the x axis is inside it but outside the
        // axis-aligned unit box.
        const c = Math.SQRT1_2;
        const box45 = obb([0, 0, 0], [vec(c, c, 0), vec(-c, c, 0),
            vec(0, 0, 1)], [1, 1, 1]);
        const s = segment([1.2, 0, 0], [3, 0, 0]);
        expect(ti.test(s, box45).intersect).toBe(true);
        expect(ti.test(s, unit).intersect).toBe(false);
        const result = fi.find(s, box45);
        expect(result.intersect).toBe(true);
        expect(result.point[0].values[0]).toBeCloseTo(1.2, 10);
        expect(result.point[1].values[0]).toBeCloseTo(Math.SQRT2, 10);
    });

    it('handles a degenerate zero-length segment', () => {
        const inside = segment([0.25, -0.5, 0.75], [0.25, -0.5, 0.75]);
        expect(ti.test(inside, unit).intersect).toBe(true);
        const insideResult = fi.find(inside, unit);
        expect(insideResult.intersect).toBe(true);
        for (let i = 0; i < 2; ++i) {
            expect(insideResult.point[i].values[0]).toBeCloseTo(0.25, 12);
            expect(insideResult.point[i].values[1]).toBeCloseTo(-0.5, 12);
            expect(insideResult.point[i].values[2]).toBeCloseTo(0.75, 12);
        }

        const outside = segment([3, 4, 5], [3, 4, 5]);
        expect(ti.test(outside, unit).intersect).toBe(false);
        expect(fi.find(outside, unit).intersect).toBe(false);
    });

    it('matches the aligned-box query when the axes are the identity', () => {
        const abTI = new IntrSegment3AlignedBox3TI();
        const abFI = new IntrSegment3AlignedBox3FI();
        const random = makeRandom(90134);
        for (let trial = 0; trial < 300; ++trial) {
            const min = [3 * random() - 3, 3 * random() - 3,
                3 * random() - 3];
            const max = [min[0] + 0.2 + 3 * random(),
                min[1] + 0.2 + 3 * random(), min[2] + 0.2 + 3 * random()];
            const aligned = AlignedBox.fromMinMax(Vector.fromArray(min),
                Vector.fromArray(max));
            const { center, extent } = aligned.getCenteredForm();
            const oriented = OrientedBox.fromCenterAxisExtent(center,
                identityAxes, extent);
            const s = segment(
                [6 * random() - 3, 6 * random() - 3, 6 * random() - 3],
                [6 * random() - 3, 6 * random() - 3, 6 * random() - 3]);

            expect(ti.test(s, oriented).intersect)
                .toBe(abTI.test(s, aligned).intersect);
            const expected = abFI.find(s, aligned);
            const actual = fi.find(s, oriented);
            expect(actual.intersect).toBe(expected.intersect);
            if (actual.intersect) {
                expect(actual.numIntersections)
                    .toBe(expected.numIntersections);
                for (let i = 0; i < 2; ++i) {
                    expect(actual.parameter[i])
                        .toBeCloseTo(expected.parameter[i], 9);
                    for (let j = 0; j < 3; ++j) {
                        expect(actual.point[i].values[j])
                            .toBeCloseTo(expected.point[i].values[j], 9);
                    }
                }
            }
        }
    });

    it('agrees with sampling, and TI agrees with FI, on random configurations',
        () => {
            const random = makeRandom(90135);
            let intersections = 0;
            for (let trial = 0; trial < 400; ++trial) {
                const axes = frame(2 * Math.PI * random(), Math.PI * random());
                const box = obb(
                    [3 * random() - 1.5, 3 * random() - 1.5,
                        3 * random() - 1.5],
                    axes,
                    [0.2 + 1.5 * random(), 0.2 + 1.5 * random(),
                        0.2 + 1.5 * random()]);
                const s = segment(
                    [6 * random() - 3, 6 * random() - 3, 6 * random() - 3],
                    [6 * random() - 3, 6 * random() - 3, 6 * random() - 3]);

                const tiResult = ti.test(s, box);
                const fiResult = fi.find(s, box);
                expect(tiResult.intersect).toBe(fiResult.intersect);

                const depth = bestSampledDepth(s, box);
                if (depth < -1e-6) {
                    // A sampled segment point is strictly inside the box.
                    expect(tiResult.intersect).toBe(true);
                }
                else if (depth > 1e-2) {
                    // Every segment point is well outside the box; the dense
                    // sampling cannot have missed a thin crossing.
                    expect(tiResult.intersect).toBe(false);
                }

                if (fiResult.intersect) {
                    ++intersections;
                    // The reported points are on the segment and on or in
                    // the box, and the parameters are ordered.
                    expect(fiResult.parameter[0])
                        .toBeLessThanOrEqual(fiResult.parameter[1] + 1e-12);
                    const { center, direction, extent } =
                        s.getCenteredForm();
                    for (let i = 0; i < 2; ++i) {
                        expect(Math.abs(fiResult.parameter[i]))
                            .toBeLessThanOrEqual(extent + 1e-9);
                        const onSegment = add(center,
                            mul(fiResult.parameter[i], direction));
                        for (let j = 0; j < 3; ++j) {
                            expect(fiResult.point[i].values[j])
                                .toBeCloseTo(onSegment.values[j], 9);
                        }
                        expect(boxSignedDepth(fiResult.point[i], box))
                            .toBeLessThan(1e-9);
                    }
                }
            }
            // Guard against a degenerate test that never intersects.
            expect(intersections).toBeGreaterThan(20);
        });
});
