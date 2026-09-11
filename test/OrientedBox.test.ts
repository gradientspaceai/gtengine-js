import { describe, it, expect } from 'vitest';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, dot, sub } from '../src/Vector.js';
import { check, compareKeys, expectClose, expectStrictWeakOrder, fc,
    orientedBox } from './helpers/arbitraries.js';

// Numeric equality that treats -0 and +0 as equal, as the C++ comparisons do.
function expectVector(v: Vector, expected: readonly number[]): void {
    expect(v.size).toBe(expected.length);
    expect(v.equals(Vector.fromArray(expected))).toBe(true);
}

describe('OrientedBox construction', () => {
    it('the default constructor is the unit-extent axis-aligned box', () => {
        const box = new OrientedBox(3);
        expect(box.dimension).toBe(3);
        expect(box.center.values).toEqual([0, 0, 0]);
        expect(box.axis.length).toBe(3);
        expect(box.axis[0].values).toEqual([1, 0, 0]);
        expect(box.axis[1].values).toEqual([0, 1, 0]);
        expect(box.axis[2].values).toEqual([0, 0, 1]);
        expect(box.extent.values).toEqual([1, 1, 1]);
    });

    it('fromCenterAxisExtent copies the inputs', () => {
        const center = Vector.fromArray([1, 2]);
        const axis = [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])];
        const extent = Vector.fromArray([3, 4]);
        const box = OrientedBox.fromCenterAxisExtent(center, axis, extent);
        center.set(0, 99);
        axis[0].set(0, 99);
        axis[1] = Vector.fromArray([9, 9]);
        extent.set(0, 99);
        expect(box.center.values).toEqual([1, 2]);
        expect(box.axis[0].values).toEqual([1, 0]);
        expect(box.axis[1].values).toEqual([0, 1]);
        expect(box.extent.values).toEqual([3, 4]);
    });

    it('fromCenterAxisExtent throws on mismatched sizes', () => {
        expect(() => OrientedBox.fromCenterAxisExtent(new Vector(3),
            [new Vector(3), new Vector(3)], new Vector(3)))
            .toThrow('OrientedBox: mismatched sizes.');
        expect(() => OrientedBox.fromCenterAxisExtent(new Vector(2),
            [new Vector(2), new Vector(3)], new Vector(2)))
            .toThrow('OrientedBox: mismatched sizes.');
    });

    it('clone is a deep copy', () => {
        const box = new OrientedBox(2);
        const copy = box.clone();
        copy.axis[0].set(0, -1);
        copy.extent.set(1, 5);
        expect(box.axis[0].values).toEqual([1, 0]);
        expect(box.extent.values).toEqual([1, 1]);
    });
});

describe('OrientedBox vertices', () => {
    it('matches hand-computed corners for an axis-aligned 2D box', () => {
        const box = OrientedBox.fromCenterAxisExtent(
            Vector.fromArray([1, 2]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([3, 4]));
        const vertex = box.getVertices();
        expect(vertex.length).toBe(4);
        expectVector(vertex[0], [-2, -2]);
        expectVector(vertex[1], [4, -2]);
        expectVector(vertex[2], [-2, 6]);
        expectVector(vertex[3], [4, 6]);
    });

    it('matches hand-computed corners for a rotated 2D box', () => {
        // Axes rotated by 90 degrees: U0 = (0,1), U1 = (-1,0).
        const box = OrientedBox.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([0, 1]), Vector.fromArray([-1, 0])],
            Vector.fromArray([1, 2]));
        const vertex = box.getVertices();
        expectVector(vertex[0], [2, -1]);
        expectVector(vertex[1], [2, 1]);
        expectVector(vertex[2], [-2, -1]);
        expectVector(vertex[3], [-2, 1]);
    });

    it('every 3D vertex has box coordinates (+-e0, +-e1, +-e2)', () => {
        const invSqrt2 = 1 / Math.sqrt(2);
        const axis = [
            Vector.fromArray([invSqrt2, invSqrt2, 0]),
            Vector.fromArray([-invSqrt2, invSqrt2, 0]),
            Vector.fromArray([0, 0, 1])
        ];
        const extent = Vector.fromArray([1, 2, 3]);
        const box = OrientedBox.fromCenterAxisExtent(
            Vector.fromArray([5, -1, 2]), axis, extent);
        const vertex = box.getVertices();
        expect(vertex.length).toBe(8);
        for (let i = 0; i < 8; ++i) {
            const delta = sub(vertex[i], box.center);
            for (let d = 0; d < 3; ++d) {
                const sign = (i & (1 << d)) > 0 ? 1 : -1;
                expect(dot(delta, axis[d]))
                    .toBeCloseTo(sign * extent.get(d), 12);
            }
        }
    });
});

describe('OrientedBox comparisons', () => {
    function box(center: number[], axis0: number[], axis1: number[],
        extent: number[]): OrientedBox {
        return OrientedBox.fromCenterAxisExtent(Vector.fromArray(center),
            [Vector.fromArray(axis0), Vector.fromArray(axis1)],
            Vector.fromArray(extent));
    }

    const a = box([0, 0], [1, 0], [0, 1], [1, 1]);
    const sameAsA = box([0, 0], [1, 0], [0, 1], [1, 1]);
    // Larger center.
    const b = box([0, 1], [1, 0], [0, 1], [1, 1]);
    // Same center, larger axis array (second axis is larger).
    const c = box([0, 0], [1, 0], [0, 2], [1, 1]);
    // Same center and axes, larger extent.
    const d = box([0, 0], [1, 0], [0, 1], [1, 2]);

    it('equals/notEquals compare center, axes and extent', () => {
        expect(a.equals(sameAsA)).toBe(true);
        expect(a.notEquals(sameAsA)).toBe(false);
        expect(a.equals(c)).toBe(false);
        expect(a.equals(d)).toBe(false);
        expect(a.notEquals(b)).toBe(true);
    });

    it('orders by center, then axes, then extent', () => {
        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(false);
        expect(a.lessThan(c)).toBe(true);
        expect(c.lessThan(a)).toBe(false);
        expect(a.lessThan(d)).toBe(true);
        expect(d.lessThan(a)).toBe(false);
        expect(a.lessThan(sameAsA)).toBe(false);
    });

    it('the non-strict comparisons are the negations of the strict ones', () => {
        expect(a.lessThanOrEqual(sameAsA)).toBe(true);
        expect(a.greaterThanOrEqual(sameAsA)).toBe(true);
        expect(a.greaterThan(sameAsA)).toBe(false);
        expect(d.greaterThan(a)).toBe(true);
        expect(d.greaterThanOrEqual(a)).toBe(true);
    });
});

describe('OrientedBox verification', () => {
    const key = (b: OrientedBox) => [...b.center.values,
        ...b.axis.flatMap(a => [...a.values]), ...b.extent.values];

    it('getVertices emits 2^N corners with sign[d] = 2*b[d] - 1', () => {
        for (const n of [2, 3]) {
            check(orientedBox(n), box => {
                const vertex = box.getVertices();
                expect(vertex.length).toBe(1 << n);
                for (let i = 0; i < vertex.length; ++i) {
                    // Recover the box coordinates of the corner; for an
                    // orthonormal frame y[d] = Dot(U[d], V - C).
                    const delta = sub(vertex[i], box.center);
                    for (let d = 0; d < n; ++d) {
                        const sign = ((i >> d) & 1) === 1 ? 1 : -1;
                        expectClose(dot(box.axis[d], delta),
                            sign * box.extent.get(d), 1e-12, 1e-12);
                    }
                }
            }, 50);
        }
    });

    it('opposite corners average to the center', () => {
        check(orientedBox(3), box => {
            const vertex = box.getVertices();
            for (let i = 0; i < 8; ++i) {
                const opposite = vertex[7 - i];
                for (let d = 0; d < 3; ++d) {
                    expectClose(0.5 * (vertex[i].get(d) + opposite.get(d)),
                        box.center.get(d), 1e-12, 1e-12);
                }
            }
        }, 50);
    });

    it('the comparisons follow the (center, axis, extent) member order',
        () => {
            check(fc.tuple(orientedBox(2), orientedBox(2)), ([a, b]) => {
                const c = compareKeys(key(a), key(b));
                expect(a.lessThan(b)).toBe(c < 0);
                expect(a.greaterThan(b)).toBe(c > 0);
                expect(a.lessThanOrEqual(b)).toBe(c <= 0);
                expect(a.greaterThanOrEqual(b)).toBe(c >= 0);
                expect(a.equals(b)).toBe(c === 0);
            });
        });

    it('the axis array is compared before the extent', () => {
        // Same center, axes differing in the first component of axis[0]: the
        // axis comparison must decide the order, whatever the extents are.
        const c = Vector.fromArray([0, 0]);
        const lo = OrientedBox.fromCenterAxisExtent(c,
            [Vector.fromArray([0, 1]), Vector.fromArray([-1, 0])],
            Vector.fromArray([9, 9]));
        const hi = OrientedBox.fromCenterAxisExtent(c,
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        expect(lo.lessThan(hi)).toBe(true);
        expect(hi.lessThan(lo)).toBe(false);
    });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(orientedBox(2), { minLength: 4, maxLength: 5 }),
            boxes => {
                expectStrictWeakOrder(boxes, (x, y) => x.lessThan(y));
            }, 30);
    });

    it('equals is element equality: a NaN axis breaks self-equality', () => {
        // Regression: upstream operator== is 'center == box.center &&
        // axis == box.axis && extent == box.extent' where std::array's
        // operator== applies Vector's element-by-element ==. Comparing the
        // axes with the lexicographic ordering instead reports a NaN axis as
        // equal to itself.
        const box = OrientedBox.fromCenterAxisExtent(Vector.fromArray([0, 0]),
            [Vector.fromArray([NaN, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        expect(box.equals(box)).toBe(false);
        expect(box.notEquals(box)).toBe(true);
        // The ordering operators are lexicographic and unaffected.
        expect(box.lessThan(box)).toBe(false);
        expect(box.lessThanOrEqual(box)).toBe(true);
    });

    it('fromCenterAxisExtent and clone copy every vector', () => {
        check(orientedBox(3), box => {
            const axis = box.axis.map(a => a.clone());
            const copy = OrientedBox.fromCenterAxisExtent(box.center, axis,
                box.extent);
            const cloned = copy.clone();
            axis[0].set(0, 999);
            copy.axis[1].set(1, 888);
            expect(copy.axis[0].get(0)).not.toBe(999);
            expect(cloned.axis[1].get(1)).not.toBe(888);
        }, 50);
    });
});
