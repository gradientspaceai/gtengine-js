import { describe, it, expect } from 'vitest';
import { Hypersphere } from '../src/Hypersphere';
import {
    IntrSegment2Circle2TI,
    IntrSegment2Circle2FI,
    defaultIntrSegment2Circle2FIResult,
    intrSegment2Circle2FIDoQuery
} from '../src/IntrSegment2Circle2';
import { Segment } from '../src/Segment';
import { Vector, add, dot, mul, sub } from '../src/Vector';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function circle(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(Vector.fromArray(center), radius);
}

function diskDepth(c: Hypersphere, P: Vector): number {
    const d = sub(P, c.center);
    return dot(d, d) - c.radius * c.radius;
}

describe('IntrSegment2Circle2', () => {
    const ti = new IntrSegment2Circle2TI();
    const fi = new IntrSegment2Circle2FI();

    it('finds both crossings of a long chord through a disk', () => {
        // The centered form has origin (0,0), direction (1,0), extent 5.
        const result = fi.find(segment([-5, 0], [5, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('clips against the segment endpoints', () => {
        // The segment stops halfway across the disk.
        const result = fi.find(segment([-5, 0], [0, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // Centered form: origin (-2.5,0), direction (1,0), extent 2.5.
        expect(result.parameter[0]).toBeCloseTo(1.5, 12);
        expect(result.parameter[1]).toBeCloseTo(2.5, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(0, 12);
    });

    it('reports a segment fully inside the disk as no boundary crossing', () => {
        const s = segment([-0.2, 0], [0.2, 0]);
        const c = circle([0, 0], 1);
        const result = fi.find(s, c);
        // The line-disk t-interval is [-1,1] and the segment interval is
        // [-0.2,0.2]; the overlap is the whole segment.
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-0.2, 12);
        expect(result.parameter[1]).toBeCloseTo(0.2, 12);
        expect(ti.test(s, c).intersect).toBe(true);
    });

    it('rejects a segment that stops short of the disk', () => {
        const s = segment([-5, 0], [-2, 0]);
        const c = circle([0, 0], 1);
        expect(fi.find(s, c).intersect).toBe(false);
        expect(fi.find(s, c).numIntersections).toBe(0);
        expect(ti.test(s, c).intersect).toBe(false);
    });

    it('rejects a segment whose line misses the circle', () => {
        const s = segment([-5, 3], [5, 3]);
        const c = circle([0, 0], 1);
        expect(fi.find(s, c).intersect).toBe(false);
        expect(ti.test(s, c).intersect).toBe(false);
    });

    it('reports a single touching point', () => {
        // The segment endpoint just reaches the circle.
        const result = fi.find(segment([-5, 0], [-1, 0]), circle([0, 0], 1));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 9);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const result = defaultIntrSegment2Circle2FIResult();
        intrSegment2Circle2FIDoQuery(vec(0, 0), vec(1, 0), 5,
            circle([0, 0], 1), result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.point[0].values).toEqual([0, 0]);
    });

    it('agrees with brute-force sampling on random segments', () => {
        let seed = 1122334;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const c = circle([-0.75, 0.5], 1.25);
        for (let trial = 0; trial < 400; ++trial) {
            const s = segment([rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 8 - 4, rand() * 8 - 4]);
            const result = fi.find(s, c);
            expect(ti.test(s, c).intersect).toBe(result.intersect);

            let sampledHit = false;
            for (let k = 0; k <= 400; ++k) {
                const t = k / 400;
                const P = add(s.p[0], mul(t, sub(s.p[1], s.p[0])));
                if (diskDepth(c, P) < -1e-6) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(result.intersect).toBe(true);
            }

            const { center, direction, extent } = s.getCenteredForm();
            for (let i = 0; i < result.numIntersections; ++i) {
                expect(Math.abs(result.parameter[i]))
                    .toBeLessThanOrEqual(extent + 1e-9);
                const P = add(center, mul(result.parameter[i], direction));
                expect(sub(P, result.point[i]).values[0]).toBeCloseTo(0, 9);
                expect(sub(P, result.point[i]).values[1]).toBeCloseTo(0, 9);
                expect(diskDepth(c, P)).toBeLessThan(1e-8);
            }
        }
    });
});
