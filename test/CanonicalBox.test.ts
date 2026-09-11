import { describe, it, expect } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { Vector } from '../src/Vector.js';
import { check, compareKeys, expectStrictWeakOrder, fc, vector }
    from './helpers/arbitraries.js';

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

describe('CanonicalBox verification', () => {
    const canonicalBox = (n: number) =>
        vector(n, 0, 5).map(e => CanonicalBox.fromExtent(e));
    const key = (b: CanonicalBox) => [...b.extent.values];

    it('getVertices emits 2^N corners with sign[d] = 2*b[d] - 1', () => {
        for (const n of [1, 2, 3, 4]) {
            check(canonicalBox(n), box => {
                const vertex = box.getVertices();
                expect(vertex.length).toBe(1 << n);
                for (let i = 0; i < vertex.length; ++i) {
                    for (let d = 0; d < n; ++d) {
                        const sign = ((i >> d) & 1) === 1 ? 1 : -1;
                        // '+ 0' normalizes -0, which toBe would distinguish.
                        expect(vertex[i].get(d) + 0)
                            .toBe(sign * box.extent.get(d) + 0);
                        // Every corner is on the boundary: |x[d]| = e[d].
                        expect(Math.abs(vertex[i].get(d)))
                            .toBe(box.extent.get(d));
                    }
                }
            }, 50);
        }
    });

    it('the comparisons follow the extent member order', () => {
        check(fc.tuple(canonicalBox(3), canonicalBox(3)), ([a, b]) => {
            const c = compareKeys(key(a), key(b));
            expect(a.lessThan(b)).toBe(c < 0);
            expect(a.greaterThan(b)).toBe(c > 0);
            expect(a.lessThanOrEqual(b)).toBe(c <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(c >= 0);
            expect(a.equals(b)).toBe(c === 0);
        });
    });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(canonicalBox(2), { minLength: 4, maxLength: 6 }),
            boxes => {
                expectStrictWeakOrder(boxes, (x, y) => x.lessThan(y));
            }, 50);
    });

    it('equals is element equality, so a NaN box does not equal itself', () => {
        const box = CanonicalBox.fromExtent(Vector.fromArray([NaN, 1]));
        expect(box.equals(box)).toBe(false);
        expect(box.lessThan(box)).toBe(false);
        expect(box.lessThanOrEqual(box)).toBe(true);
    });

    it('fromExtent and clone are independent of their inputs', () => {
        check(vector(3, 0, 5), e => {
            const box = CanonicalBox.fromExtent(e);
            const copy = box.clone();
            e.set(0, 999);
            box.extent.set(1, 777);
            expect(box.extent.get(0)).not.toBe(999);
            expect(copy.extent.get(1)).not.toBe(777);
        });
    });

    it('a zero extent produces 2^N coincident corners at the origin', () => {
        const box = new CanonicalBox(3);
        const vertex = box.getVertices();
        expect(vertex.length).toBe(8);
        for (const v of vertex) {
            expect(v.values.map(x => x + 0)).toEqual([0, 0, 0]);
        }
    });
});
