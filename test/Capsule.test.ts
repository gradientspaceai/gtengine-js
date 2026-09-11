import { describe, it, expect } from 'vitest';
import { Capsule } from '../src/Capsule.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, sub, mul, length, normalize } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { DistPointSegment } from '../src/DistPointSegment.js';
import { check, compareKeys, expectClose, expectStrictWeakOrder, fc, finite,
    positive, segment } from './helpers/arbitraries.js';

// Component equality that treats -0 and +0 as equal, as the C++ comparisons
// do. (Segment's default p0 is the negation of a unit vector, so its zero
// components are negative zeros.)
function expectVector(v: Vector, expected: readonly number[]): void {
    expect(v.size).toBe(expected.length);
    expect(v.equals(Vector.fromArray(expected))).toBe(true);
}

describe('Capsule construction', () => {
    it('the default constructor is the unit-radius capsule about [-e0, e0]', () => {
        const capsule = new Capsule(3);
        expect(capsule.dimension).toBe(3);
        expect(capsule.radius).toBe(1);
        expectVector(capsule.segment.p[0], [-1, 0, 0]);
        expectVector(capsule.segment.p[1], [1, 0, 0]);
    });

    it('works for any dimension N', () => {
        expectVector(new Capsule(2).segment.p[0], [-1, 0]);
        expectVector(new Capsule(4).segment.p[1], [1, 0, 0, 0]);
        expect(new Capsule(4).dimension).toBe(4);
    });

    it('fromSegmentRadius copies the segment', () => {
        const segment = Segment.fromEndpoints(Vector.fromArray([0, 0, 0]),
            Vector.fromArray([0, 0, 4]));
        const capsule = Capsule.fromSegmentRadius(segment, 2);
        segment.p[1].set(2, 99);
        expect(capsule.segment.p[1].values).toEqual([0, 0, 4]);
        expect(capsule.radius).toBe(2);
    });

    it('clone is a deep copy', () => {
        const capsule = new Capsule(3);
        const copy = capsule.clone();
        copy.radius = 7;
        copy.segment.p[0].set(0, -5);
        expect(capsule.radius).toBe(1);
        expectVector(capsule.segment.p[0], [-1, 0, 0]);
    });
});

describe('Capsule comparisons', () => {
    it('equals compares segment and radius', () => {
        const a = new Capsule(3);
        const b = new Capsule(3);
        expect(a.equals(b)).toBe(true);
        expect(a.notEquals(b)).toBe(false);

        b.radius = 2;
        expect(a.equals(b)).toBe(false);
        expect(a.notEquals(b)).toBe(true);

        const c = new Capsule(3);
        c.segment = Segment.fromEndpoints(Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 0, 0]));
        expect(a.equals(c)).toBe(false);
    });

    it('lessThan orders by segment first, then radius', () => {
        const a = new Capsule(3);
        const bigRadius = new Capsule(3);
        bigRadius.radius = 5;
        expect(a.lessThan(bigRadius)).toBe(true);
        expect(bigRadius.lessThan(a)).toBe(false);

        // A segment ordered before the default one dominates the radius.
        const earlier = new Capsule(3);
        earlier.segment = Segment.fromEndpoints(Vector.fromArray([-2, 0, 0]),
            Vector.fromArray([1, 0, 0]));
        earlier.radius = 100;
        expect(earlier.lessThan(a)).toBe(true);
        expect(a.greaterThan(earlier)).toBe(true);
    });

    it('the derived comparisons are consistent', () => {
        const a = new Capsule(3);
        const b = new Capsule(3);
        b.radius = 3;
        expect(a.lessThanOrEqual(b)).toBe(true);
        expect(a.lessThanOrEqual(a.clone())).toBe(true);
        expect(a.greaterThanOrEqual(a.clone())).toBe(true);
        expect(b.greaterThanOrEqual(a)).toBe(true);
        expect(b.greaterThan(a)).toBe(true);
        expect(a.greaterThan(a.clone())).toBe(false);
    });
});

describe('Capsule verification', () => {
    const capsule = (n: number) => fc.tuple(segment(n), positive(5))
        .map(([s, r]) => Capsule.fromSegmentRadius(s, r));
    const key = (c: Capsule) =>
        [...c.segment.p[0].values, ...c.segment.p[1].values, c.radius];

    it('the comparisons follow the (segment, radius) member order', () => {
        check(fc.tuple(capsule(3), capsule(3)), ([a, b]) => {
            const cmp = compareKeys(key(a), key(b));
            expect(a.lessThan(b)).toBe(cmp < 0);
            expect(a.greaterThan(b)).toBe(cmp > 0);
            expect(a.lessThanOrEqual(b)).toBe(cmp <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(cmp >= 0);
            expect(a.equals(b)).toBe(cmp === 0);
        });
    });

    it('a differing radius alone orders by the radius', () => {
        check(fc.tuple(segment(3), positive(5), positive(5)),
            ([s, r0, r1]) => {
                const a = Capsule.fromSegmentRadius(s, r0);
                const b = Capsule.fromSegmentRadius(s, r1);
                expect(a.lessThan(b)).toBe(r0 < r1);
            });
    });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(capsule(2), { minLength: 4, maxLength: 6 }), cs => {
            expectStrictWeakOrder(cs, (x, y) => x.lessThan(y));
        }, 50);
    });

    it('the surface is at distance radius from the segment', () => {
        // Cross-check of the primitive's definition against the library's
        // point-segment distance query: C + r*U for a unit U perpendicular to
        // the segment direction is at distance r from the segment.
        const query = new DistPointSegment();
        check(fc.tuple(capsule(3), finite(0, 1)), ([cap, t]) => {
            const d = sub(cap.segment.p[1], cap.segment.p[0]);
            const dir = d.clone();
            normalize(dir);
            const perp = cross(dir, Vector.fromArray([
                dir.get(1) + 1, dir.get(2) - 1, dir.get(0) + 2]));
            if (length(perp) < 1e-6) { return; }
            normalize(perp);
            const base = add(cap.segment.p[0], mul(t, d));
            const x = add(base, mul(cap.radius, perp));
            const result = query.compute(x, cap.segment);
            expectClose(result.distance, cap.radius, 1e-9, 1e-9);
        }, 100);
    });

    it('equals is element equality: a NaN endpoint breaks self-equality',
        () => {
            const cap = Capsule.fromSegmentRadius(
                Segment.fromEndpoints(Vector.fromArray([NaN, 0, 0]),
                    Vector.fromArray([1, 0, 0])), 1);
            expect(cap.equals(cap)).toBe(false);
            expect(cap.lessThan(cap)).toBe(false);
            expect(cap.lessThanOrEqual(cap)).toBe(true);
        });

    it('the factory and clone copy the segment', () => {
        check(segment(3), s => {
            const cap = Capsule.fromSegmentRadius(s, 2);
            const cloned = cap.clone();
            s.p[0].set(0, 999);
            cap.segment.p[1].set(1, 888);
            expect(cap.segment.p[0].get(0)).not.toBe(999);
            expect(cloned.segment.p[1].get(1)).not.toBe(888);
        });
    });
});
