import { describe, it, expect } from 'vitest';
import { Capsule } from '../src/Capsule';
import {
    IntrSegment3Capsule3TI,
    IntrSegment3Capsule3FI,
    defaultIntrSegment3Capsule3FIResult,
    intrSegment3Capsule3FIDoQuery
} from '../src/IntrSegment3Capsule3';
import { Segment } from '../src/Segment';
import { Vector, add, dot, mul, sub } from '../src/Vector';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function capsule(p0: number[], p1: number[], radius: number): Capsule {
    return Capsule.fromSegmentRadius(segment(p0, p1), radius);
}

function capsuleSignedDepth(c: Capsule, P: Vector): number {
    const p0 = c.segment.p[0], p1 = c.segment.p[1];
    const e = sub(p1, p0);
    const d = sub(P, p0);
    const ee = dot(e, e);
    let t = ee > 0 ? dot(d, e) / ee : 0;
    t = Math.max(0, Math.min(1, t));
    const closest = add(p0, mul(t, e));
    const r = sub(P, closest);
    return Math.sqrt(dot(r, r)) - c.radius;
}

describe('IntrSegment3Capsule3', () => {
    const ti = new IntrSegment3Capsule3TI();
    const fi = new IntrSegment3Capsule3FI();

    it('finds both crossings of a capsule along the x-axis', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        // Centered form: origin (0,0,0), direction (1,0,0), extent 10.
        const result = fi.find(segment([-10, 0, 0], [10, 0, 0]), c);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-2, 9);
        expect(result.parameter[1]).toBeCloseTo(2, 9);
        expect(result.point[0].values[0]).toBeCloseTo(-2, 9);
        expect(result.point[1].values[0]).toBeCloseTo(2, 9);
    });

    it('clips against the segment endpoints', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        // Centered form: origin (-5,0,0), direction (1,0,0), extent 5.
        const result = fi.find(segment([-10, 0, 0], [0, 0, 0]), c);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(3, 9);
        expect(result.parameter[1]).toBeCloseTo(5, 9);
        expect(result.point[0].values[0]).toBeCloseTo(-2, 9);
        expect(result.point[1].values[0]).toBeCloseTo(0, 9);
    });

    it('rejects a segment that stops short of the capsule', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const s = segment([-10, 0, 0], [-5, 0, 0]);
        const result = fi.find(s, c);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
        expect(ti.test(s, c).intersect).toBe(false);
    });

    it('rejects a segment whose line misses the capsule', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const s = segment([-10, 5, 0], [10, 5, 0]);
        expect(fi.find(s, c).intersect).toBe(false);
        expect(ti.test(s, c).intersect).toBe(false);
    });

    it('reports the test query using the segment-segment distance', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        expect(ti.test(segment([-10, 1, 0], [10, 1, 0]), c).intersect)
            .toBe(true);
        expect(ti.test(segment([-10, 1.0001, 0], [10, 1.0001, 0]), c)
            .intersect).toBe(false);
        // A degenerate segment inside the capsule.
        expect(ti.test(segment([0, 0, 0], [0, 0, 0]), c).intersect).toBe(true);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const c = capsule([-1, 0, 0], [1, 0, 0], 1);
        const result = defaultIntrSegment3Capsule3FIResult();
        intrSegment3Capsule3FIDoQuery(vec(0, 0, 0), vec(1, 0, 0), 10, c,
            result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-2, 9);
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('agrees with brute-force sampling on random segments', () => {
        let seed = 31415926;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const c = capsule([-1, 0.5, 0.25], [1.5, -0.5, 1], 0.9);
        for (let trial = 0; trial < 250; ++trial) {
            const s = segment([rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4]);
            const result = fi.find(s, c);

            let sampledHit = false;
            for (let k = 0; k <= 400; ++k) {
                const t = k / 400;
                const P = add(s.p[0], mul(t, sub(s.p[1], s.p[0])));
                if (capsuleSignedDepth(c, P) < -1e-6) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(result.intersect).toBe(true);
                expect(ti.test(s, c).intersect).toBe(true);
            }

            if (result.intersect) {
                const { center, direction, extent } = s.getCenteredForm();
                for (let i = 0; i < 2; ++i) {
                    expect(Math.abs(result.parameter[i]))
                        .toBeLessThanOrEqual(extent + 1e-9);
                    const P = add(center,
                        mul(result.parameter[i], direction));
                    expect(sub(P, result.point[i]).values[0])
                        .toBeCloseTo(0, 9);
                    expect(capsuleSignedDepth(c, P)).toBeLessThan(1e-7);
                }
            }
        }
    });
});
