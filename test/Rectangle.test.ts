import { describe, it, expect } from 'vitest';
import { Rectangle } from '../src/Rectangle.js';
import { Vector, dot, sub } from '../src/Vector.js';
import { check, compareKeys, expectClose, expectStrictWeakOrder, fc, positive,
    rotationFrame, vector } from './helpers/arbitraries.js';

// Numeric equality that treats -0 and +0 as equal, as the C++ comparisons do.
function expectVector(v: Vector, expected: readonly number[]): void {
    expect(v.size).toBe(expected.length);
    expect(v.equals(Vector.fromArray(expected))).toBe(true);
}

describe('Rectangle construction', () => {
    it('the default constructor is the unit-extent rectangle in the x-y plane', () => {
        const rectangle = new Rectangle(3);
        expect(rectangle.dimension).toBe(3);
        expect(rectangle.center.values).toEqual([0, 0, 0]);
        expect(rectangle.axis.length).toBe(2);
        expect(rectangle.axis[0].values).toEqual([1, 0, 0]);
        expect(rectangle.axis[1].values).toEqual([0, 1, 0]);
        expect(rectangle.extent.values).toEqual([1, 1]);
    });

    it('the extent always has two components, whatever N is', () => {
        expect(new Rectangle(2).extent.size).toBe(2);
        expect(new Rectangle(4).extent.size).toBe(2);
        expect(new Rectangle(4).center.size).toBe(4);
        expect(new Rectangle(4).axis[1].values).toEqual([0, 1, 0, 0]);
    });

    it('fromCenterAxisExtent copies the inputs', () => {
        const center = Vector.fromArray([1, 1, 0]);
        const axis = [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0])];
        const extent = Vector.fromArray([2, 3]);
        const rectangle = Rectangle.fromCenterAxisExtent(center, axis, extent);
        center.set(0, 99);
        axis[0].set(0, 99);
        extent.set(0, 99);
        expect(rectangle.center.values).toEqual([1, 1, 0]);
        expect(rectangle.axis[0].values).toEqual([1, 0, 0]);
        expect(rectangle.extent.values).toEqual([2, 3]);
    });

    it('fromCenterAxisExtent throws on mismatched sizes', () => {
        expect(() => Rectangle.fromCenterAxisExtent(new Vector(3),
            [new Vector(3), new Vector(3)], new Vector(3)))
            .toThrow('Rectangle: mismatched sizes.');
        expect(() => Rectangle.fromCenterAxisExtent(new Vector(3),
            [new Vector(3), new Vector(2)], new Vector(2)))
            .toThrow('Rectangle: mismatched sizes.');
    });

    it('clone is a deep copy', () => {
        const rectangle = new Rectangle(2);
        const copy = rectangle.clone();
        copy.center.set(0, 5);
        copy.extent.set(0, 9);
        expect(rectangle.center.values).toEqual([0, 0]);
        expect(rectangle.extent.values).toEqual([1, 1]);
    });
});

describe('Rectangle vertices', () => {
    it('matches hand-computed corners of an axis-aligned rectangle in 3D', () => {
        const rectangle = Rectangle.fromCenterAxisExtent(
            Vector.fromArray([1, 1, 0]),
            [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0])],
            Vector.fromArray([2, 3]));
        const vertex = rectangle.getVertices();
        expect(vertex.length).toBe(4);
        expectVector(vertex[0], [-1, -2, 0]);
        expectVector(vertex[1], [3, -2, 0]);
        expectVector(vertex[2], [-1, 4, 0]);
        expectVector(vertex[3], [3, 4, 0]);
    });

    it('vertex[i] has rectangle coordinates (+-e0, +-e1) per bit pattern', () => {
        const invSqrt2 = 1 / Math.sqrt(2);
        const axis = [
            Vector.fromArray([invSqrt2, invSqrt2, 0]),
            Vector.fromArray([0, 0, 1])
        ];
        const extent = Vector.fromArray([2, 5]);
        const rectangle = Rectangle.fromCenterAxisExtent(
            Vector.fromArray([1, -2, 3]), axis, extent);
        const vertex = rectangle.getVertices();
        for (let i = 0; i < 4; ++i) {
            const delta = sub(vertex[i], rectangle.center);
            for (let d = 0; d < 2; ++d) {
                const sign = (i & (1 << d)) > 0 ? 1 : -1;
                expect(dot(delta, axis[d]))
                    .toBeCloseTo(sign * extent.get(d), 12);
            }
        }
    });
});

describe('Rectangle comparisons', () => {
    function rect(center: number[], axis0: number[], axis1: number[],
        extent: number[]): Rectangle {
        return Rectangle.fromCenterAxisExtent(Vector.fromArray(center),
            [Vector.fromArray(axis0), Vector.fromArray(axis1)],
            Vector.fromArray(extent));
    }

    const a = rect([0, 0], [1, 0], [0, 1], [1, 1]);
    const sameAsA = rect([0, 0], [1, 0], [0, 1], [1, 1]);
    // Larger center.
    const b = rect([1, 0], [1, 0], [0, 1], [1, 1]);
    // Same center, larger second axis.
    const c = rect([0, 0], [1, 0], [0, 2], [1, 1]);
    // Same center and axes, larger extent.
    const d = rect([0, 0], [1, 0], [0, 1], [1, 3]);

    it('equals/notEquals compare center, axes and extent', () => {
        expect(a.equals(sameAsA)).toBe(true);
        expect(a.notEquals(sameAsA)).toBe(false);
        expect(a.equals(b)).toBe(false);
        expect(a.equals(c)).toBe(false);
        expect(a.equals(d)).toBe(false);
        expect(a.notEquals(d)).toBe(true);
    });

    it('orders by center, then axes, then extent', () => {
        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(false);
        expect(a.lessThan(c)).toBe(true);
        expect(c.lessThan(a)).toBe(false);
        expect(a.lessThan(d)).toBe(true);
        expect(d.lessThan(a)).toBe(false);
        expect(a.lessThanOrEqual(sameAsA)).toBe(true);
        expect(a.greaterThanOrEqual(sameAsA)).toBe(true);
        expect(a.greaterThan(sameAsA)).toBe(false);
        expect(d.greaterThan(a)).toBe(true);
    });
});

describe('Rectangle verification', () => {
    const rectangle = (n: number) =>
        fc.tuple(vector(n), rotationFrame(n === 2 ? 2 : 3),
            fc.array(positive(5), { minLength: 2, maxLength: 2 }))
            .map(([c, frame, ext]) => Rectangle.fromCenterAxisExtent(c,
                [frame[0], frame[1]], Vector.fromArray(ext)));
    const key = (r: Rectangle) => [...r.center.values,
        ...r.axis.flatMap(a => [...a.values]), ...r.extent.values];

    it('getVertices emits the four corners in bit-pattern order', () => {
        // vertex[i] = C + sign[0]*e0*A0 + sign[1]*e1*A1 with sign[d] = 2b[d]-1
        // and i = b[1]b[0]. NOTE: the boundary traversal order is 0,1,3,2;
        // several consumers (V22, V30) rely on this bit ordering.
        for (const n of [2, 3]) {
            check(rectangle(n), rect => {
                const vertex = rect.getVertices();
                expect(vertex.length).toBe(4);
                for (let i = 0; i < 4; ++i) {
                    const delta = sub(vertex[i], rect.center);
                    for (let d = 0; d < 2; ++d) {
                        const sign = ((i >> d) & 1) === 1 ? 1 : -1;
                        expectClose(dot(rect.axis[d], delta),
                            sign * rect.extent.get(d), 1e-12, 1e-12);
                    }
                }
            }, 50);
        }
    });

    it('the diagonals of the vertex quad share the center', () => {
        check(rectangle(3), rect => {
            for (const [i, j] of [[0, 3], [1, 2]]) {
                for (let d = 0; d < 3; ++d) {
                    expectClose(0.5 * (rect.getVertices()[i].get(d)
                        + rect.getVertices()[j].get(d)),
                        rect.center.get(d), 1e-12, 1e-12);
                }
            }
        }, 50);
    });

    it('the comparisons follow the (center, axis, extent) member order',
        () => {
            check(fc.tuple(rectangle(3), rectangle(3)), ([a, b]) => {
                const c = compareKeys(key(a), key(b));
                expect(a.lessThan(b)).toBe(c < 0);
                expect(a.greaterThan(b)).toBe(c > 0);
                expect(a.lessThanOrEqual(b)).toBe(c <= 0);
                expect(a.greaterThanOrEqual(b)).toBe(c >= 0);
                expect(a.equals(b)).toBe(c === 0);
            });
        });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(rectangle(2), { minLength: 4, maxLength: 5 }), rs => {
            expectStrictWeakOrder(rs, (x, y) => x.lessThan(y));
        }, 30);
    });

    it('equals is element equality, so a NaN axis breaks self-equality', () => {
        const rect = Rectangle.fromCenterAxisExtent(Vector.fromArray([0, 0]),
            [Vector.fromArray([NaN, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([1, 1]));
        expect(rect.equals(rect)).toBe(false);
        expect(rect.lessThan(rect)).toBe(false);
        expect(rect.lessThanOrEqual(rect)).toBe(true);
    });

    it('fromCenterAxisExtent and clone copy every vector', () => {
        check(rectangle(3), rect => {
            const axis = rect.axis.map(a => a.clone());
            const copy = Rectangle.fromCenterAxisExtent(rect.center, axis,
                rect.extent);
            const cloned = copy.clone();
            axis[0].set(0, 999);
            copy.axis[1].set(1, 888);
            copy.extent.set(0, 777);
            expect(copy.axis[0].get(0)).not.toBe(999);
            expect(cloned.axis[1].get(1)).not.toBe(888);
            expect(cloned.extent.get(0)).not.toBe(777);
        }, 50);
    });
});
