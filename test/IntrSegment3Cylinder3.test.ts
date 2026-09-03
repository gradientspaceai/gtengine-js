import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3.js';
import {
    IntrSegment3Cylinder3FI,
    defaultIntrSegment3Cylinder3FIResult,
    intrSegment3Cylinder3FIDoQuery
} from '../src/IntrSegment3Cylinder3.js';
import { IntrLine3Cylinder3FI } from '../src/IntrLine3Cylinder3.js';
import { Line } from '../src/Line.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function cylinder(origin: number[], direction: number[], radius: number,
    height: number): Cylinder3 {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(Vector.fromArray(origin), d), radius, height);
}

// The signed "outside" amount of a point relative to the solid cylinder:
// negative strictly inside, zero on the boundary, positive outside.
function cylinderSignedDepth(p: Vector, c: Cylinder3): number {
    const diff = sub(p, c.axis.origin);
    const z = dot(diff, c.axis.direction);
    const radial = sub(diff, mul(z, c.axis.direction));
    const rDepth = Math.sqrt(dot(radial, radial)) - c.radius;
    const zDepth = Math.abs(z) - 0.5 * c.height;
    return Math.max(rDepth, zDepth);
}

// The most deeply contained sampled point of the segment.
function bestSampledDepth(s: Segment, c: Cylinder3): number {
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= 4000; ++i) {
        const p = add(s.p[0], mul(i / 4000, sub(s.p[1], s.p[0])));
        const value = cylinderSignedDepth(p, c);
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

describe('IntrSegment3Cylinder3', () => {
    const fi = new IntrSegment3Cylinder3FI();
    // A cylinder about the z axis: radius 1, height 2, so |z| <= 1.
    const unit = cylinder([0, 0, 0], [0, 0, 1], 1, 2);

    it('default results report no intersection', () => {
        const result = defaultIntrSegment3Cylinder3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('clips a segment crossing the cylinder wall', () => {
        const s = segment([-3, 0, 0], [3, 0, 0]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The parameters are relative to the centered form, extent 3.
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('clips a segment crossing both end caps', () => {
        const s = segment([0, 0, -3], [0, 0, 3]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[2]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[2]).toBeCloseTo(1, 12);
    });

    it('reports the whole segment when it lies inside the cylinder', () => {
        const s = segment([-0.5, 0, 0], [0.5, 0, 0]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-0.5, 12);
        expect(result.parameter[1]).toBeCloseTo(0.5, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-0.5, 12);
        expect(result.point[1].values[0]).toBeCloseTo(0.5, 12);
    });

    it('clips a segment whose endpoints straddle an end cap', () => {
        const s = segment([0, 0, 0], [0, 0, 4]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The centered form has center (0,0,2), extent 2, so the overlap of
        // the line interval with [-2,2] is [-2,-1].
        expect(result.parameter[0]).toBeCloseTo(-2, 12);
        expect(result.parameter[1]).toBeCloseTo(-1, 12);
        expect(result.point[0].values[2]).toBeCloseTo(0, 12);
        expect(result.point[1].values[2]).toBeCloseTo(1, 12);
    });

    it('reports a tangent contact on the cylinder wall', () => {
        const s = segment([1, -3, 0], [1, 3, 0]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports a tangent contact on the cap rim', () => {
        const s = segment([1, -3, 1], [1, 3, 1]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[2]).toBeCloseTo(1, 12);
    });

    it('reports no intersection for a segment that misses the cylinder',
        () => {
            const s = segment([5, 0, 0], [5, 0, 3]);
            expect(fi.find(s, unit).intersect).toBe(false);
        });

    it('reports no intersection when the line hits but the segment does not',
        () => {
            // The line x = 0, y = 0 passes through the cylinder, but the
            // segment lies entirely above the top cap.
            const s = segment([0, 0, 3], [0, 0, 5]);
            expect(fi.find(s, unit).intersect).toBe(false);
        });

    it('clips a diagonal segment symmetrically', () => {
        const s = segment([-3, -3, -3], [3, 3, 3]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The direction is (1,1,1)/sqrt(3); the wall x^2+y^2 = 1 is met at
        // x = y = +-1/sqrt(2), which is inside |z| <= 1.
        const h = Math.SQRT1_2;
        for (let i = 0; i < 3; ++i) {
            expect(result.point[0].values[i]).toBeCloseTo(-h, 10);
            expect(result.point[1].values[i]).toBeCloseTo(h, 10);
        }
        expect(result.parameter[0]).toBeCloseTo(-h * Math.sqrt(3), 10);
        expect(result.parameter[1]).toBeCloseTo(h * Math.sqrt(3), 10);
    });

    it('handles a tilted, translated cylinder', () => {
        const c = cylinder([1, 2, 3], [1, 1, 0], 0.5, 4);
        // A segment along the cylinder axis, longer than the cylinder.
        const h = Math.SQRT1_2;
        const s = segment([1 - 5 * h, 2 - 5 * h, 3], [1 + 5 * h, 2 + 5 * h, 3]);
        const result = fi.find(s, c);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-2, 10);
        expect(result.parameter[1]).toBeCloseTo(2, 10);
        expect(cylinderSignedDepth(result.point[0], c)).toBeCloseTo(0, 9);
        expect(cylinderSignedDepth(result.point[1], c)).toBeCloseTo(0, 9);
    });

    it('rejects a zero-length segment outside the cylinder', () => {
        const s = segment([5, 5, 5], [5, 5, 5]);
        expect(fi.find(s, unit).intersect).toBe(false);
    });

    it('reports a degenerate zero-length segment inside the cylinder', () => {
        // A zero-length segment has no well-defined centered-form direction
        // (Normalize of the zero vector leaves it zero), so the quadratic in
        // the line-cylinder DoQuery is identically zero and the reported
        // parameters are NaN. This matches the upstream C++ behavior; see
        // the port notes in IntrSegment3Cylinder3.ts.
        const s = segment([0.2, 0.1, 0.3], [0.2, 0.1, 0.3]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(Number.isNaN(result.parameter[0])).toBe(true);
        expect(Number.isNaN(result.parameter[1])).toBe(true);
    });

    it('rejects an infinite cylinder', () => {
        const infinite = cylinder([0, 0, 0], [0, 0, 1], 1, 2);
        infinite.makeInfiniteCylinder();
        expect(() => fi.find(segment([-3, 0, 0], [3, 0, 0]), infinite))
            .toThrow('Infinite cylinders are not yet supported.');
    });

    it('exposes the DoQuery helper, which matches the public query', () => {
        const s = segment([-3, 0, 0.5], [3, 0, 0.5]);
        const { center, direction, extent } = s.getCenteredForm();
        const result = defaultIntrSegment3Cylinder3FIResult();
        intrSegment3Cylinder3FIDoQuery(center, direction, extent, unit,
            result);
        const expected = fi.find(s, unit);
        expect(result.intersect).toBe(expected.intersect);
        expect(result.numIntersections).toBe(expected.numIntersections);
        expect(result.parameter).toEqual(expected.parameter);
        // DoQuery does not compute the points; the public query does.
        expect(result.point[0].values[0]).toBe(0);
        expect(expected.point[0].values[0]).toBeCloseTo(-1, 12);
    });

    it('restricts the line query to the segment t-interval', () => {
        // Whenever the segment query reports an intersection, the line query
        // must report one too, and the segment parameters must lie inside
        // both the line interval and [-extent, extent].
        const lineQuery = new IntrLine3Cylinder3FI();
        const random = makeRandom(90941);
        let intersections = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const c = cylinder(
                [2 * random() - 1, 2 * random() - 1, 2 * random() - 1],
                [2 * random() - 1, 2 * random() - 1, 2 * random() - 1],
                0.2 + random(), 0.4 + 2 * random());
            const s = segment(
                [4 * random() - 2, 4 * random() - 2, 4 * random() - 2],
                [4 * random() - 2, 4 * random() - 2, 4 * random() - 2]);
            const { center, direction, extent } = s.getCenteredForm();
            const segResult = fi.find(s, c);
            const lineResult = lineQuery.find(
                Line.fromOriginDirection(center, direction), c);
            if (segResult.intersect) {
                ++intersections;
                expect(lineResult.intersect).toBe(true);
                for (let i = 0; i < 2; ++i) {
                    expect(segResult.parameter[i])
                        .toBeGreaterThanOrEqual(
                            Math.max(lineResult.parameter[0], -extent) - 1e-9);
                    expect(segResult.parameter[i])
                        .toBeLessThanOrEqual(
                            Math.min(lineResult.parameter[1], extent) + 1e-9);
                }
            }
        }
        expect(intersections).toBeGreaterThan(20);
    });

    it('agrees with sampling on random configurations', () => {
        const random = makeRandom(90942);
        let intersections = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const c = cylinder(
                [2 * random() - 1, 2 * random() - 1, 2 * random() - 1],
                [2 * random() - 1, 2 * random() - 1, 2 * random() - 1],
                0.2 + random(), 0.4 + 2 * random());
            const s = segment(
                [5 * random() - 2.5, 5 * random() - 2.5, 5 * random() - 2.5],
                [5 * random() - 2.5, 5 * random() - 2.5, 5 * random() - 2.5]);

            const result = fi.find(s, c);
            const depth = bestSampledDepth(s, c);
            if (depth < -1e-6) {
                expect(result.intersect).toBe(true);
            }
            else if (depth > 1e-2) {
                expect(result.intersect).toBe(false);
            }

            if (result.intersect) {
                ++intersections;
                expect(result.parameter[0])
                    .toBeLessThanOrEqual(result.parameter[1] + 1e-12);
                const { center, direction, extent } = s.getCenteredForm();
                for (let i = 0; i < 2; ++i) {
                    expect(Math.abs(result.parameter[i]))
                        .toBeLessThanOrEqual(extent + 1e-9);
                    const onSegment = add(center,
                        mul(result.parameter[i], direction));
                    for (let j = 0; j < 3; ++j) {
                        expect(result.point[i].values[j])
                            .toBeCloseTo(onSegment.values[j], 9);
                    }
                    expect(cylinderSignedDepth(result.point[i], c))
                        .toBeLessThan(1e-9);
                }
            }
        }
        expect(intersections).toBeGreaterThan(20);
    });
});
