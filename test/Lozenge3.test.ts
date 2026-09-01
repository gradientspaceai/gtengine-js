import { describe, it, expect } from 'vitest';
import { Lozenge3 } from '../src/Lozenge3';
import { Rectangle } from '../src/Rectangle';
import { Vector } from '../src/Vector';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

describe('Lozenge3 construction', () => {
    it('the default constructor is the unit rectangle with radius 1', () => {
        const lozenge = new Lozenge3();
        expect(lozenge.radius).toBe(1);
        expect(lozenge.rectangle.dimension).toBe(3);
        expect(lozenge.rectangle.center.values).toEqual([0, 0, 0]);
        expect(lozenge.rectangle.axis[0].values).toEqual([1, 0, 0]);
        expect(lozenge.rectangle.axis[1].values).toEqual([0, 1, 0]);
        expect(lozenge.rectangle.extent.values).toEqual([1, 1]);
    });

    it('fromRectangleRadius copies the rectangle', () => {
        const rectangle = Rectangle.fromCenterAxisExtent(v3(1, 2, 3),
            [v3(1, 0, 0), v3(0, 0, 1)], Vector.fromArray([2, 5]));
        const lozenge = Lozenge3.fromRectangleRadius(rectangle, 3);
        rectangle.center.set(0, 99);
        rectangle.extent.set(0, 99);
        expect(lozenge.rectangle.center.values).toEqual([1, 2, 3]);
        expect(lozenge.rectangle.extent.values).toEqual([2, 5]);
        expect(lozenge.radius).toBe(3);
    });

    it('rejects a rectangle that is not in 3D', () => {
        expect(() => Lozenge3.fromRectangleRadius(new Rectangle(2),
            1)).toThrow();
    });

    it('clone is a deep copy', () => {
        const lozenge = new Lozenge3();
        const copy = lozenge.clone();
        copy.radius = 4;
        copy.rectangle.center.set(1, 9);
        expect(lozenge.radius).toBe(1);
        expect(lozenge.rectangle.center.values).toEqual([0, 0, 0]);
    });
});

describe('Lozenge3 comparisons', () => {
    it('equals compares the rectangle and the radius', () => {
        const a = new Lozenge3();
        const b = new Lozenge3();
        expect(a.equals(b)).toBe(true);
        expect(a.notEquals(b)).toBe(false);

        b.radius = 2;
        expect(a.equals(b)).toBe(false);
        expect(a.notEquals(b)).toBe(true);

        const c = new Lozenge3();
        c.rectangle.extent.set(0, 3);
        expect(a.equals(c)).toBe(false);
    });

    it('lessThan orders by rectangle first, then radius', () => {
        const a = new Lozenge3();

        const bigRadius = new Lozenge3();
        bigRadius.radius = 5;
        expect(a.lessThan(bigRadius)).toBe(true);
        expect(bigRadius.lessThan(a)).toBe(false);

        // An earlier rectangle dominates a larger radius.
        const earlier = new Lozenge3();
        earlier.rectangle.center = v3(-1, 0, 0);
        earlier.radius = 100;
        expect(earlier.lessThan(a)).toBe(true);
        expect(a.greaterThan(earlier)).toBe(true);
    });

    it('the derived comparisons are consistent', () => {
        const a = new Lozenge3();
        const b = new Lozenge3();
        b.radius = 3;
        expect(a.lessThanOrEqual(b)).toBe(true);
        expect(a.lessThanOrEqual(a.clone())).toBe(true);
        expect(b.greaterThanOrEqual(a)).toBe(true);
        expect(b.greaterThan(a)).toBe(true);
        expect(a.greaterThan(a.clone())).toBe(false);
        expect(a.greaterThanOrEqual(a.clone())).toBe(true);
    });
});
