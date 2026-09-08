import { describe, it, expect } from 'vitest';
import { Lozenge3 } from '../src/Lozenge3.js';
import { Rectangle } from '../src/Rectangle.js';
import { Vector, add, mul } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { DistPointRectangle } from '../src/DistPointRectangle.js';
import { check, compareKeys, expectClose, expectStrictWeakOrder, fc, finite,
    positive, rotationFrame, vector } from './helpers/arbitraries.js';

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

describe('Lozenge3 verification', () => {
    const rectangle3 = () => fc.tuple(vector(3, -5, 5), rotationFrame(3),
        fc.array(positive(5), { minLength: 2, maxLength: 2 }))
        .map(([c, frame, ext]) => Rectangle.fromCenterAxisExtent(c,
            [frame[0], frame[1]], Vector.fromArray(ext)));
    const lozenge = () => fc.tuple(rectangle3(), positive(5))
        .map(([r, radius]) => Lozenge3.fromRectangleRadius(r, radius));
    const key = (l: Lozenge3) => [...l.rectangle.center.values,
        ...l.rectangle.axis[0].values, ...l.rectangle.axis[1].values,
        ...l.rectangle.extent.values, l.radius];

    it('the comparisons follow the (rectangle, radius) member order', () => {
        check(fc.tuple(lozenge(), lozenge()), ([a, b]) => {
            const cmp = compareKeys(key(a), key(b));
            expect(a.lessThan(b)).toBe(cmp < 0);
            expect(a.greaterThan(b)).toBe(cmp > 0);
            expect(a.lessThanOrEqual(b)).toBe(cmp <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(cmp >= 0);
            expect(a.equals(b)).toBe(cmp === 0);
        });
    });

    it('a differing radius alone orders by the radius', () => {
        check(fc.tuple(rectangle3(), positive(5), positive(5)),
            ([r, r0, r1]) => {
                const a = Lozenge3.fromRectangleRadius(r, r0);
                const b = Lozenge3.fromRectangleRadius(r, r1);
                expect(a.lessThan(b)).toBe(r0 < r1);
            });
    });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(lozenge(), { minLength: 4, maxLength: 5 }), ls => {
            expectStrictWeakOrder(ls, (x, y) => x.lessThan(y));
        }, 30);
    });

    it('the surface is at distance radius from the rectangle', () => {
        // The definition: points equidistant from the rectangle. A point on
        // the rectangle displaced by radius along the rectangle normal is on
        // the lozenge boundary.
        const query = new DistPointRectangle();
        check(fc.tuple(lozenge(), finite(-1, 1), finite(-1, 1)),
            ([l, s0, s1]) => {
                const rect = l.rectangle;
                const normal = cross(rect.axis[0], rect.axis[1]);
                const p = add(rect.center,
                    add(mul(s0 * rect.extent.get(0), rect.axis[0]),
                        mul(s1 * rect.extent.get(1), rect.axis[1])));
                const x = add(p, mul(l.radius, normal));
                expectClose(query.compute(x, rect).distance, l.radius,
                    1e-9, 1e-9);
            }, 100);
    });

    it('equals is element equality, so a NaN axis breaks self-equality', () => {
        const rect = Rectangle.fromCenterAxisExtent(
            Vector.fromArray([0, 0, 0]),
            [Vector.fromArray([NaN, 0, 0]), Vector.fromArray([0, 1, 0])],
            Vector.fromArray([1, 1]));
        const l = Lozenge3.fromRectangleRadius(rect, 1);
        expect(l.equals(l)).toBe(false);
        expect(l.lessThan(l)).toBe(false);
        expect(l.lessThanOrEqual(l)).toBe(true);
    });

    it('the factory and clone copy the rectangle', () => {
        check(rectangle3(), rect => {
            const l = Lozenge3.fromRectangleRadius(rect, 1);
            const cloned = l.clone();
            rect.center.set(0, 999);
            l.rectangle.axis[0].set(1, 888);
            expect(l.rectangle.center.get(0)).not.toBe(999);
            expect(cloned.rectangle.axis[0].get(1)).not.toBe(888);
        }, 50);
    });
});
