import { describe, it, expect } from 'vitest';
import { Capsule } from '../src/Capsule.js';
import { Line } from '../src/Line.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, normalize, sub, length } from '../src/Vector.js';
import {
    IntrLine3Capsule3TI,
    IntrLine3Capsule3FI
} from '../src/IntrLine3Capsule3.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function line(p: number[], d: number[]): Line {
    const dir = vec(d);
    normalize(dir);
    return Line.fromOriginDirection(vec(p), dir);
}

function capsule(p0: number[], p1: number[], radius: number): Capsule {
    return Capsule.fromSegmentRadius(
        Segment.fromEndpoints(vec(p0), vec(p1)), radius);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Distance from a point to a segment, computed independently of the library.
function distanceToSegment(x: Vector, p0: Vector, p1: Vector): number {
    const d = sub(p1, p0);
    const dd = dot(d, d);
    let t = dd > 0 ? dot(sub(x, p0), d) / dd : 0;
    t = Math.min(1, Math.max(0, t));
    return length(sub(x, add(p0, mul(t, d))));
}

describe('IntrLine3Capsule3', () => {
    const ti = new IntrLine3Capsule3TI();
    const fi = new IntrLine3Capsule3FI();

    // Capsule with axis segment from (0,0,-2) to (0,0,2) and radius 1.
    const cap = capsule([0, 0, -2], [0, 0, 2], 1);

    it('finds the chord of a line crossing the capsule wall', () => {
        const l = line([0, 0, 0], [1, 0, 0]);
        expect(ti.test(l, cap).intersect).toBe(true);
        const result = fi.find(l, cap);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('handles a line parallel to the capsule axis', () => {
        const l = line([0.5, 0, -10], [0, 0, 1]);
        const result = fi.find(l, cap);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        const zOffset = Math.sqrt(1 - 0.25) + 2;
        expect(result.parameter[0]).toBeCloseTo(10 - zOffset, 12);
        expect(result.parameter[1]).toBeCloseTo(10 + zOffset, 12);
        // Both points are exactly on the capsule boundary.
        for (const p of result.point) {
            expect(distanceToSegment(p, vec([0, 0, -2]), vec([0, 0, 2])))
                .toBeCloseTo(1, 10);
        }
    });

    it('reports no intersection for a parallel line outside the radius', () => {
        const l = line([2, 0, -10], [0, 0, 1]);
        expect(ti.test(l, cap).intersect).toBe(false);
        const result = fi.find(l, cap);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('finds the chord through a hemispherical cap only', () => {
        // The plane z = 2.5 cuts the top cap in a circle of radius
        // sqrt(1 - 0.25).
        const l = line([0, 0, 2.5], [1, 0, 0]);
        const result = fi.find(l, cap);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        const half = Math.sqrt(0.75);
        expect(result.parameter[0]).toBeCloseTo(-half, 10);
        expect(result.parameter[1]).toBeCloseTo(half, 10);
    });

    it('reports a single point for a line tangent to the capsule wall', () => {
        const l = line([1, -10, 0], [0, 1, 0]);
        const result = fi.find(l, cap);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(10, 10);
        expect(result.parameter[1]).toBeCloseTo(result.parameter[0], 12);
        expect(result.point[0].values[0]).toBeCloseTo(1, 10);
        expect(result.point[0].values[1]).toBeCloseTo(0, 10);
    });

    it('misses a capsule entirely', () => {
        const l = line([5, 5, 0], [1, 0, 0]);
        expect(ti.test(l, cap).intersect).toBe(false);
        expect(fi.find(l, cap).intersect).toBe(false);
    });

    it('degrades on a zero-length capsule segment (upstream limitation)', () => {
        // A zero-length segment has no centered-form direction, so the
        // capsule coordinate frame built by the find query is degenerate. The
        // distance-based test query is still correct, but the find query
        // returns parameters that do not correspond to the sphere. The port
        // preserves upstream behavior; callers must supply a nondegenerate
        // capsule segment.
        const sphere = capsule([1, 2, 3], [1, 2, 3], 2);
        const l = line([1, 2, -10], [0, 0, 1]);
        expect(ti.test(l, sphere).intersect).toBe(true);

        const result = fi.find(l, sphere);
        // The true intersections would be at t = 11 and t = 15.
        expect(result.parameter[0]).not.toBeCloseTo(11, 6);
        expect(result.parameter[1]).not.toBeCloseTo(15, 6);
        for (const p of result.point) {
            expect(length(sub(p, vec([1, 2, 3])))).not.toBeCloseTo(2, 6);
        }
    });

    it('agrees with the distance-based test query and dense sampling', () => {
        const rand = makeRandom(8675309);
        const p0 = vec([-1, 0.5, 0]);
        const p1 = vec([2, -1, 1.5]);
        const cap2 = Capsule.fromSegmentRadius(
            Segment.fromEndpoints(p0, p1), 0.8);
        for (let trial = 0; trial < 120; ++trial) {
            const l = line(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            const result = fi.find(l, cap2);
            expect(ti.test(l, cap2).intersect).toBe(result.intersect);

            let tLo = Number.POSITIVE_INFINITY;
            let tHi = Number.NEGATIVE_INFINITY;
            const n = 20000;
            for (let k = 0; k <= n; ++k) {
                const t = -10 + (20 * k) / n;
                const x = add(l.origin, mul(t, l.direction));
                if (distanceToSegment(x, p0, p1) <= 0.8) {
                    if (t < tLo) { tLo = t; }
                    if (t > tHi) { tHi = t; }
                }
            }

            if (tLo <= tHi) {
                expect(result.intersect).toBe(true);
                expect(result.parameter[0]).toBeLessThanOrEqual(tLo + 1e-9);
                expect(result.parameter[1]).toBeGreaterThanOrEqual(tHi - 1e-9);
                expect(tLo - result.parameter[0]).toBeLessThan(3e-3);
                expect(result.parameter[1] - tHi).toBeLessThan(3e-3);
            }

            if (result.intersect && result.numIntersections === 2) {
                // The endpoints of the chord lie on the capsule boundary.
                for (const p of result.point) {
                    expect(distanceToSegment(p, p0, p1)).toBeCloseTo(0.8, 8);
                }
            }
        }
    });
});
