import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';
import { DistPointLine } from '../src/DistPointLine.js';
import { check, compareKeys, expectClose, expectStrictWeakOrder, fc, finite,
    line, positive } from './helpers/arbitraries.js';
import { Line } from '../src/Line.js';
import { Vector, add, mul } from '../src/Vector.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

describe('Cylinder3 construction', () => {
    it('the default constructor uses the default 3D line, radius 1, height 1', () => {
        const cylinder = new Cylinder3();
        expect(cylinder.axis.origin.values).toEqual([0, 0, 0]);
        // Line3<T>() has direction (1,0,0), despite the upstream comment.
        expect(cylinder.axis.direction.values).toEqual([1, 0, 0]);
        expect(cylinder.radius).toBe(1);
        expect(cylinder.height).toBe(1);
    });

    it('fromAxisRadiusHeight copies the line', () => {
        const axis = Line.fromOriginDirection(v3(1, 2, 3), v3(0, 0, 1));
        const cylinder = Cylinder3.fromAxisRadiusHeight(axis, 2, 6);
        axis.origin.set(0, 99);
        axis.direction.set(2, 99);
        expect(cylinder.axis.origin.values).toEqual([1, 2, 3]);
        expect(cylinder.axis.direction.values).toEqual([0, 0, 1]);
        expect(cylinder.radius).toBe(2);
        expect(cylinder.height).toBe(6);
    });

    it('rejects an axis line that is not 3D', () => {
        expect(() => Cylinder3.fromAxisRadiusHeight(new Line(2), 1,
            1)).toThrow();
    });

    it('clone is a deep copy', () => {
        const cylinder = new Cylinder3();
        const copy = cylinder.clone();
        copy.axis.origin.set(0, 5);
        copy.height = 8;
        expect(cylinder.axis.origin.values).toEqual([0, 0, 0]);
        expect(cylinder.height).toBe(1);
    });
});

describe('Cylinder3 finite/infinite state', () => {
    it('the default cylinder is finite', () => {
        const cylinder = new Cylinder3();
        expect(cylinder.isFinite()).toBe(true);
        expect(cylinder.isInfinite()).toBe(false);
    });

    it('makeInfiniteCylinder sets the height sentinel to -1', () => {
        const cylinder = new Cylinder3();
        cylinder.makeInfiniteCylinder();
        expect(cylinder.height).toBe(-1);
        expect(cylinder.isFinite()).toBe(false);
        expect(cylinder.isInfinite()).toBe(true);
    });

    it('makeFiniteCylinder ignores negative heights', () => {
        const cylinder = new Cylinder3();
        cylinder.makeInfiniteCylinder();
        cylinder.makeFiniteCylinder(-3);
        expect(cylinder.height).toBe(-1);
        expect(cylinder.isInfinite()).toBe(true);

        cylinder.makeFiniteCylinder(4);
        expect(cylinder.height).toBe(4);
        expect(cylinder.isFinite()).toBe(true);
    });

    it('zero height is finite (a degenerate disk)', () => {
        const cylinder = new Cylinder3();
        cylinder.makeFiniteCylinder(0);
        expect(cylinder.height).toBe(0);
        expect(cylinder.isFinite()).toBe(true);
        expect(cylinder.isInfinite()).toBe(false);
    });
});

describe('Cylinder3 comparisons', () => {
    const base = new Cylinder3();

    it('equals compares axis, radius and height', () => {
        expect(base.equals(new Cylinder3())).toBe(true);
        expect(base.notEquals(new Cylinder3())).toBe(false);

        const other = base.clone();
        other.height = 2;
        expect(base.equals(other)).toBe(false);
        expect(base.notEquals(other)).toBe(true);
    });

    it('lessThan orders by axis, then radius, then height', () => {
        const earlierAxis = base.clone();
        earlierAxis.axis = Line.fromOriginDirection(v3(-1, 0, 0),
            v3(1, 0, 0));
        earlierAxis.radius = 100;
        earlierAxis.height = 100;
        expect(earlierAxis.lessThan(base)).toBe(true);

        const smallRadius = base.clone();
        smallRadius.radius = 0.5;
        smallRadius.height = 100;
        expect(smallRadius.lessThan(base)).toBe(true);

        const smallHeight = base.clone();
        smallHeight.height = 0.5;
        expect(smallHeight.lessThan(base)).toBe(true);
        expect(base.lessThan(smallHeight)).toBe(false);
    });

    it('the derived comparisons are consistent', () => {
        const taller = base.clone();
        taller.height = 5;
        expect(base.lessThanOrEqual(taller)).toBe(true);
        expect(base.lessThanOrEqual(base.clone())).toBe(true);
        expect(taller.greaterThan(base)).toBe(true);
        expect(taller.greaterThanOrEqual(base)).toBe(true);
        expect(base.greaterThan(base.clone())).toBe(false);
        expect(base.greaterThanOrEqual(base.clone())).toBe(true);
    });
});

describe('Cylinder3 verification', () => {
    const cylinder = () => fc.tuple(line(3), positive(5), positive(5))
        .map(([a, r, h]) => Cylinder3.fromAxisRadiusHeight(a, r, h));
    const key = (c: Cylinder3) => [...c.axis.origin.values,
        ...c.axis.direction.values, c.radius, c.height];

    it('the height sentinel partitions finite from infinite cylinders', () => {
        // Upstream uses height = -1 (not infinity()) for infinite cylinders,
        // so IsFinite is height >= 0 and IsInfinite is height < 0.
        check(finite(-5, 5), h => {
            const c = Cylinder3.fromAxisRadiusHeight(new Line(3), 1, h);
            expect(c.isFinite()).toBe(h >= 0);
            expect(c.isInfinite()).toBe(h < 0);
            expect(c.isFinite()).toBe(!c.isInfinite());
        });
    });

    it('makeFiniteCylinder ignores negative heights and makeInfinite wins',
        () => {
            check(finite(-5, 5), h => {
                const c = new Cylinder3();
                c.makeInfiniteCylinder();
                expect(c.height).toBe(-1);
                expect(c.isInfinite()).toBe(true);
                c.makeFiniteCylinder(h);
                if (h >= 0) {
                    expect(c.height).toBe(h);
                    expect(c.isFinite()).toBe(true);
                } else {
                    // The negative height is rejected; the sentinel stands.
                    expect(c.height).toBe(-1);
                    expect(c.isInfinite()).toBe(true);
                }
            });
        });

    it('the finite cylinder wall is at distance radius from the axis', () => {
        // The definition: points at distance R from the axis line. Checked
        // against the library's point-line distance query.
        const query = new DistPointLine();
        check(fc.tuple(cylinder(), finite(-3, 3), finite(-Math.PI, Math.PI)),
            ([c, t, angle]) => {
                const basis = [c.axis.direction.clone(), new Vector(3),
                    new Vector(3)];
                computeOrthogonalComplement3(1, basis);
                const x = add(add(c.axis.origin, mul(t, c.axis.direction)),
                    add(mul(c.radius * Math.cos(angle), basis[1]),
                        mul(c.radius * Math.sin(angle), basis[2])));
                expectClose(query.compute(x, c.axis).distance, c.radius,
                    1e-9, 1e-9);
            }, 100);
    });

    it('the comparisons follow the (axis, radius, height) member order', () => {
        check(fc.tuple(cylinder(), cylinder()), ([a, b]) => {
            const cmp = compareKeys(key(a), key(b));
            expect(a.lessThan(b)).toBe(cmp < 0);
            expect(a.greaterThan(b)).toBe(cmp > 0);
            expect(a.lessThanOrEqual(b)).toBe(cmp <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(cmp >= 0);
            expect(a.equals(b)).toBe(cmp === 0);
        });
    });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(cylinder(), { minLength: 4, maxLength: 5 }), cs => {
            expectStrictWeakOrder(cs, (x, y) => x.lessThan(y));
        }, 30);
    });

    it('equals is element equality, so a NaN axis breaks self-equality', () => {
        const c = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(Vector.fromArray([0, 0, 0]),
                Vector.fromArray([NaN, 0, 1])), 1, 1);
        expect(c.equals(c)).toBe(false);
        expect(c.lessThan(c)).toBe(false);
        expect(c.lessThanOrEqual(c)).toBe(true);
    });

    it('the factory and clone copy the axis line', () => {
        check(line(3), axis => {
            const c = Cylinder3.fromAxisRadiusHeight(axis, 1, 1);
            const cloned = c.clone();
            axis.origin.set(0, 999);
            c.axis.direction.set(1, 888);
            expect(c.axis.origin.get(0)).not.toBe(999);
            expect(cloned.axis.direction.get(1)).not.toBe(888);
        });
    });
});
