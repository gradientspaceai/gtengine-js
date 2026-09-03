import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane';
import {
    IntrSegment3Plane3TI,
    IntrSegment3Plane3FI,
    defaultIntrSegment3Plane3FIResult,
    defaultIntrSegment3Plane3TIResult,
    intrSegment3Plane3FIDoQuery
} from '../src/IntrSegment3Plane3';
import { Segment } from '../src/Segment';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

const ti = new IntrSegment3Plane3TI();
const fi = new IntrSegment3Plane3FI();

describe('IntrSegment3Plane3', () => {
    it('default-constructs results as no intersection', () => {
        expect(defaultIntrSegment3Plane3TIResult()).toEqual({ intersect: false });
        const r = defaultIntrSegment3Plane3FIResult();
        expect(r.intersect).toBe(false);
        expect(r.numIntersections).toBe(0);
    });

    it('finds a transverse intersection at a known point', () => {
        const S = segment([0, 0, -1], [0, 0, 3]);
        const P = plane([0, 0, 1], [0, 0, 1]);
        const result = fi.find(S, P);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        // The centered form has center (0,0,1), direction (0,0,1), extent 2,
        // so the intersection parameter is 0.
        expect(result.parameter).toBeCloseTo(0, 12);
        expect(result.point.values[2]).toBeCloseTo(1, 12);
        expect(ti.test(S, P).intersect).toBe(true);
    });

    it('rejects a segment entirely on one side of the plane', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        for (const S of [segment([0, 0, 1], [0, 0, 5]),
            segment([0, 0, -5], [0, 0, -1])]) {
            expect(fi.find(S, P).intersect).toBe(false);
            expect(fi.find(S, P).numIntersections).toBe(0);
            expect(ti.test(S, P).intersect).toBe(false);
        }
    });

    it('accepts an endpoint touching the plane', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        const S0 = segment([1, 2, 0], [1, 2, 4]);
        expect(ti.test(S0, P).intersect).toBe(true);
        const r0 = fi.find(S0, P);
        expect(r0.intersect).toBe(true);
        expect(r0.point.values[2]).toBeCloseTo(0, 12);

        const S1 = segment([1, 2, -4], [1, 2, 0]);
        expect(ti.test(S1, P).intersect).toBe(true);
        expect(fi.find(S1, P).intersect).toBe(true);
    });

    it('reports a segment lying in the plane as infinitely many hits', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        const S = segment([-1, -1, 0], [2, 3, 0]);
        expect(ti.test(S, P).intersect).toBe(true);
        const r = fi.find(S, P);
        expect(r.intersect).toBe(true);
        expect(r.numIntersections).toBe(2147483647);
    });

    it('reports a parallel disjoint segment as no intersection', () => {
        const P = plane([0, 0, 1], [0, 0, 2]);
        const S = segment([-1, -1, 0], [2, 3, 0]);
        expect(ti.test(S, P).intersect).toBe(false);
        expect(fi.find(S, P).intersect).toBe(false);
    });

    it('exposes the DoQuery helper used by derived queries', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        // The line hits the plane at parameter 5, outside the extent 2.
        const outside = defaultIntrSegment3Plane3FIResult();
        intrSegment3Plane3FIDoQuery(vec(0, 0, -5), vec(0, 0, 1), 2, P, outside);
        expect(outside.intersect).toBe(false);
        expect(outside.numIntersections).toBe(0);

        const inside = defaultIntrSegment3Plane3FIResult();
        intrSegment3Plane3FIDoQuery(vec(0, 0, -1), vec(0, 0, 1), 2, P, inside);
        expect(inside.intersect).toBe(true);
        expect(inside.parameter).toBeCloseTo(1, 12);
    });

    it('agrees with TI and puts the point on both primitives (randomized)', () => {
        let seed = 24680135;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (): number => 4 * rand() - 2;

        let hits = 0;
        for (let trial = 0; trial < 4000; ++trial) {
            const p0 = [rnd(), rnd(), rnd()];
            const p1 = [rnd(), rnd(), rnd()];
            const S = segment(p0, p1);
            const P = plane([rnd(), rnd(), rnd()], [rnd(), rnd(), rnd()]);
            const fiResult = fi.find(S, P);
            const tiResult = ti.test(S, P);
            expect(fiResult.intersect).toBe(tiResult.intersect);
            if (fiResult.intersect && fiResult.numIntersections === 1) {
                ++hits;
                // The point is on the plane.
                expect(dot(P.normal, fiResult.point) - P.constant)
                    .toBeCloseTo(0, 8);
                // The point is on the segment: X = p0 + s * (p1 - p0) with
                // s in [0,1].
                const d = sub(S.p[1], S.p[0]);
                const s = dot(sub(fiResult.point, S.p[0]), d) / dot(d, d);
                expect(s).toBeGreaterThanOrEqual(-1e-9);
                expect(s).toBeLessThanOrEqual(1 + 1e-9);
                const onSeg = add(S.p[0], mul(s, d));
                for (let i = 0; i < 3; ++i) {
                    expect(fiResult.point.values[i]).toBeCloseTo(
                        onSeg.values[i], 8);
                }
            }
        }
        expect(hits).toBeGreaterThan(500);
    });
});
