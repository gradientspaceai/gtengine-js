import { describe, it, expect } from 'vitest';
import { Rectangle } from '../src/Rectangle';
import { Vector, dot, sub } from '../src/Vector';

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
