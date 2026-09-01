import { describe, it, expect } from 'vitest';
import { Capsule } from '../src/Capsule';
import { Segment } from '../src/Segment';
import { Vector } from '../src/Vector';

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
