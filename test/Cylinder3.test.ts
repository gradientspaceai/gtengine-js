import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3';
import { Line } from '../src/Line';
import { Vector } from '../src/Vector';

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
