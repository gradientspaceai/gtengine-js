import { describe, it, expect } from 'vitest';
import { Segment } from '../src/Segment.js';
import { Vector, add, sub, mul, length } from '../src/Vector.js';

// Numeric equality that treats -0 and +0 as equal, as the C++ comparisons do.
function expectVector(v: Vector, expected: readonly number[]): void {
    expect(v.size).toBe(expected.length);
    expect(v.equals(Vector.fromArray(expected))).toBe(true);
}

describe('Segment construction', () => {
    it('the default constructor is the segment from (-1,0,..) to (1,0,..)', () => {
        const segment = new Segment(3);
        expect(segment.dimension).toBe(3);
        expect(segment.p.length).toBe(2);
        expectVector(segment.p[0], [-1, 0, 0]);
        expectVector(segment.p[1], [1, 0, 0]);
    });

    it('fromEndpoints copies the input vectors', () => {
        const p0 = Vector.fromArray([1, 2]);
        const p1 = Vector.fromArray([4, 6]);
        const segment = Segment.fromEndpoints(p0, p1);
        p0.set(0, 99);
        p1.set(1, 99);
        expect(segment.p[0].values).toEqual([1, 2]);
        expect(segment.p[1].values).toEqual([4, 6]);
    });

    it('fromPointArray takes the two endpoints as an array', () => {
        const segment = Segment.fromPointArray([Vector.fromArray([0, 0]),
            Vector.fromArray([1, 1])]);
        expect(segment.p[0].values).toEqual([0, 0]);
        expect(segment.p[1].values).toEqual([1, 1]);
        expect(() => Segment.fromPointArray([new Vector(2)]))
            .toThrow('Segment: invalid number of endpoints.');
    });

    it('fromEndpoints throws on mismatched sizes', () => {
        expect(() => Segment.fromEndpoints(new Vector(2), new Vector(3)))
            .toThrow('Segment: mismatched sizes.');
    });

    it('clone is a deep copy', () => {
        const segment = Segment.fromEndpoints(Vector.fromArray([0, 0]),
            Vector.fromArray([2, 2]));
        const copy = segment.clone();
        copy.p[1].set(0, 9);
        expect(segment.p[1].values).toEqual([2, 2]);
        expect(copy.p[1].values).toEqual([9, 2]);
    });
});

describe('Segment centered form', () => {
    it('getCenteredForm matches hand-computed values for a 3-4-5 segment', () => {
        const segment = Segment.fromEndpoints(Vector.fromArray([1, 2]),
            Vector.fromArray([4, 6]));
        const { center, direction, extent } = segment.getCenteredForm();
        expect(center.values).toEqual([2.5, 4]);
        expect(direction.get(0)).toBeCloseTo(0.6, 15);
        expect(direction.get(1)).toBeCloseTo(0.8, 15);
        expect(length(direction)).toBeCloseTo(1, 15);
        expect(extent).toBeCloseTo(2.5, 15);
    });

    it('setCenteredForm produces C -+ e*D', () => {
        const segment = new Segment(2);
        segment.setCenteredForm(Vector.fromArray([1, 1]),
            Vector.fromArray([0, 1]), 3);
        expectVector(segment.p[0], [1, -2]);
        expectVector(segment.p[1], [1, 4]);
    });

    it('fromCenteredForm agrees with setCenteredForm', () => {
        const center = Vector.fromArray([-1, 2, 3]);
        const direction = Vector.fromArray([0, 0, 1]);
        const built = Segment.fromCenteredForm(center, direction, 4);
        const set = new Segment(3);
        set.setCenteredForm(center, direction, 4);
        expect(built.equals(set)).toBe(true);
        expectVector(built.p[0], [-1, 2, -1]);
        expectVector(built.p[1], [-1, 2, 7]);
    });

    it('endpoints round-trip through the centered form', () => {
        const p0 = Vector.fromArray([-3, 5, 1]);
        const p1 = Vector.fromArray([7, -2, 4]);
        const segment = Segment.fromEndpoints(p0, p1);
        const { center, direction, extent } = segment.getCenteredForm();
        const q0 = sub(center, mul(extent, direction));
        const q1 = add(center, mul(extent, direction));
        // Upstream notes that round-off can make q0 != p0 exactly, so the
        // comparison is approximate.
        for (let i = 0; i < 3; ++i) {
            expect(q0.get(i)).toBeCloseTo(p0.get(i), 12);
            expect(q1.get(i)).toBeCloseTo(p1.get(i), 12);
        }
        expect(extent).toBeCloseTo(0.5 * length(sub(p1, p0)), 12);
    });

    it('a degenerate segment has zero extent and zero direction', () => {
        const p = Vector.fromArray([2, 2]);
        const segment = Segment.fromEndpoints(p, p);
        const { center, direction, extent } = segment.getCenteredForm();
        expect(center.values).toEqual([2, 2]);
        expect(extent).toBe(0);
        // Vector.normalize sets a zero-length vector to zero.
        expect(direction.values).toEqual([0, 0]);
    });
});

describe('Segment comparisons', () => {
    const a = Segment.fromEndpoints(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 0]));
    const sameAsA = Segment.fromEndpoints(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 0]));
    // Larger first endpoint.
    const b = Segment.fromEndpoints(Vector.fromArray([0, 1]),
        Vector.fromArray([1, 0]));
    // Same first endpoint, larger second endpoint.
    const c = Segment.fromEndpoints(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 1]));

    it('equals/notEquals compare the endpoint arrays', () => {
        expect(a.equals(sameAsA)).toBe(true);
        expect(a.notEquals(sameAsA)).toBe(false);
        expect(a.equals(c)).toBe(false);
        expect(a.notEquals(b)).toBe(true);
    });

    it('orders lexicographically by the endpoint array', () => {
        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(false);
        expect(a.lessThan(c)).toBe(true);
        expect(c.lessThan(a)).toBe(false);
        expect(a.lessThanOrEqual(sameAsA)).toBe(true);
        expect(a.greaterThanOrEqual(sameAsA)).toBe(true);
        expect(a.greaterThan(sameAsA)).toBe(false);
        expect(c.greaterThan(a)).toBe(true);
        expect(c.greaterThanOrEqual(a)).toBe(true);
    });

    it('reversing the endpoints yields a different segment', () => {
        const reversed = Segment.fromEndpoints(a.p[1], a.p[0]);
        expect(a.equals(reversed)).toBe(false);
        expect(a.lessThan(reversed)).toBe(true);
    });
});
