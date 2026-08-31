import { describe, it, expect } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox';
import { Vector } from '../src/Vector';

describe('CanonicalBox construction', () => {
    it('the default constructor sets all members to zero', () => {
        const box = new CanonicalBox(3);
        expect(box.dimension).toBe(3);
        expect(box.extent.values).toEqual([0, 0, 0]);
    });

    it('fromExtent copies the input vector', () => {
        const extent = Vector.fromArray([1, 2, 3]);
        const box = CanonicalBox.fromExtent(extent);
        extent.set(0, 99);
        expect(box.extent.values).toEqual([1, 2, 3]);
    });

    it('clone is a deep copy', () => {
        const box = CanonicalBox.fromExtent(Vector.fromArray([4, 5]));
        const copy = box.clone();
        copy.extent.set(1, 0);
        expect(box.extent.values).toEqual([4, 5]);
        expect(copy.extent.values).toEqual([4, 0]);
    });
});

describe('CanonicalBox vertices', () => {
    it('produces the 4 corners of a 2D box in bit-pattern order', () => {
        const box = CanonicalBox.fromExtent(Vector.fromArray([2, 3]));
        const vertex = box.getVertices();
        expect(vertex.length).toBe(4);
        expect(vertex[0].values).toEqual([-2, -3]);
        expect(vertex[1].values).toEqual([2, -3]);
        expect(vertex[2].values).toEqual([-2, 3]);
        expect(vertex[3].values).toEqual([2, 3]);
    });

    it('vertex[i][d] has sign 2*b[d]-1 in 3D', () => {
        const extent = Vector.fromArray([1, 2, 4]);
        const box = CanonicalBox.fromExtent(extent);
        const vertex = box.getVertices();
        expect(vertex.length).toBe(8);
        for (let i = 0; i < 8; ++i) {
            for (let d = 0; d < 3; ++d) {
                const sign = (i & (1 << d)) > 0 ? 1 : -1;
                expect(vertex[i].get(d)).toBe(sign * extent.get(d));
            }
        }
    });

    it('a degenerate extent collapses the corresponding coordinate', () => {
        const box = CanonicalBox.fromExtent(Vector.fromArray([5, 0]));
        const vertex = box.getVertices();
        for (const v of vertex) {
            expect(Math.abs(v.get(0))).toBe(5);
            expect(v.get(1)).toBe(0);
        }
    });
});

describe('CanonicalBox comparisons', () => {
    const a = CanonicalBox.fromExtent(Vector.fromArray([1, 2]));
    const sameAsA = CanonicalBox.fromExtent(Vector.fromArray([1, 2]));
    const b = CanonicalBox.fromExtent(Vector.fromArray([1, 3]));

    it('equals/notEquals compare the extents', () => {
        expect(a.equals(sameAsA)).toBe(true);
        expect(a.notEquals(sameAsA)).toBe(false);
        expect(a.equals(b)).toBe(false);
        expect(a.notEquals(b)).toBe(true);
    });

    it('orders lexicographically by extent', () => {
        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(false);
        expect(a.lessThanOrEqual(sameAsA)).toBe(true);
        expect(a.greaterThanOrEqual(sameAsA)).toBe(true);
        expect(b.greaterThan(a)).toBe(true);
        expect(a.greaterThan(b)).toBe(false);
    });
});
