import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Vector } from '../src/Vector.js';

describe('AlignedBox construction', () => {
    it('the default constructor sets min to -1 and max to +1', () => {
        const box = new AlignedBox(3);
        expect(box.dimension).toBe(3);
        expect(box.min.values).toEqual([-1, -1, -1]);
        expect(box.max.values).toEqual([1, 1, 1]);
    });

    it('fromMinMax copies the input vectors', () => {
        const min = Vector.fromArray([-2, -3]);
        const max = Vector.fromArray([4, 5]);
        const box = AlignedBox.fromMinMax(min, max);
        min.set(0, 99);
        max.set(1, 99);
        expect(box.min.values).toEqual([-2, -3]);
        expect(box.max.values).toEqual([4, 5]);
        expect(box.dimension).toBe(2);
    });

    it('fromMinMax throws on mismatched sizes', () => {
        expect(() => AlignedBox.fromMinMax(new Vector(2), new Vector(3)))
            .toThrow('AlignedBox: mismatched sizes.');
    });

    it('clone is a deep copy', () => {
        const box = AlignedBox.fromMinMax(Vector.fromArray([0, 0]),
            Vector.fromArray([1, 2]));
        const copy = box.clone();
        copy.min.set(0, -7);
        expect(box.min.values).toEqual([0, 0]);
        expect(copy.min.values).toEqual([-7, 0]);
    });
});

describe('AlignedBox centered form', () => {
    it('matches hand-computed center and extent', () => {
        const box = AlignedBox.fromMinMax(Vector.fromArray([-1, 2, 3]),
            Vector.fromArray([3, 6, 4]));
        const { center, extent } = box.getCenteredForm();
        expect(center.values).toEqual([1, 4, 3.5]);
        expect(extent.values).toEqual([2, 2, 0.5]);
    });

    it('round-trips a box built from the centered form', () => {
        const box = AlignedBox.fromMinMax(Vector.fromArray([-4, 0]),
            Vector.fromArray([2, 10]));
        const { center, extent } = box.getCenteredForm();
        const min = Vector.fromArray([center.get(0) - extent.get(0),
            center.get(1) - extent.get(1)]);
        const max = Vector.fromArray([center.get(0) + extent.get(0),
            center.get(1) + extent.get(1)]);
        expect(AlignedBox.fromMinMax(min, max).equals(box)).toBe(true);
    });
});

describe('AlignedBox vertices', () => {
    it('produces the 4 corners of a 2D box in bit-pattern order', () => {
        const box = AlignedBox.fromMinMax(Vector.fromArray([-1, -2]),
            Vector.fromArray([3, 5]));
        const vertex = box.getVertices();
        expect(vertex.length).toBe(4);
        expect(vertex[0].values).toEqual([-1, -2]);
        expect(vertex[1].values).toEqual([3, -2]);
        expect(vertex[2].values).toEqual([-1, 5]);
        expect(vertex[3].values).toEqual([3, 5]);
    });

    it('produces the 8 corners of a 3D box, vertex[i][d] chosen by bit d', () => {
        const min = Vector.fromArray([0, 10, 20]);
        const max = Vector.fromArray([1, 11, 21]);
        const box = AlignedBox.fromMinMax(min, max);
        const vertex = box.getVertices();
        expect(vertex.length).toBe(8);
        for (let i = 0; i < 8; ++i) {
            for (let d = 0; d < 3; ++d) {
                const expected = (i & (1 << d)) > 0
                    ? max.get(d) : min.get(d);
                expect(vertex[i].get(d)).toBe(expected);
            }
        }
    });
});

describe('AlignedBox comparisons', () => {
    const a = AlignedBox.fromMinMax(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 1]));
    const sameAsA = AlignedBox.fromMinMax(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 1]));
    // Larger min.
    const b = AlignedBox.fromMinMax(Vector.fromArray([0, 1]),
        Vector.fromArray([1, 1]));
    // Same min as a, larger max.
    const c = AlignedBox.fromMinMax(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 2]));

    it('equals/notEquals compare min and max', () => {
        expect(a.equals(sameAsA)).toBe(true);
        expect(a.notEquals(sameAsA)).toBe(false);
        expect(a.equals(c)).toBe(false);
        expect(a.notEquals(c)).toBe(true);
    });

    it('orders by min first, then by max', () => {
        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(false);
        expect(a.lessThan(c)).toBe(true);
        expect(c.lessThan(a)).toBe(false);
        expect(a.lessThan(sameAsA)).toBe(false);
    });

    it('the non-strict comparisons are the negations of the strict ones', () => {
        expect(a.lessThanOrEqual(sameAsA)).toBe(true);
        expect(a.greaterThanOrEqual(sameAsA)).toBe(true);
        expect(a.greaterThan(sameAsA)).toBe(false);
        expect(c.greaterThan(a)).toBe(true);
        expect(c.greaterThanOrEqual(a)).toBe(true);
        expect(a.lessThanOrEqual(c)).toBe(true);
    });
});
