import { describe, it, expect } from 'vitest';
import { AxisAngle } from '../src/AxisAngle';
import { Vector } from '../src/Vector';

describe('AxisAngle construction', () => {
    it('default constructor produces the 3D zero axis and zero angle', () => {
        const aa = new AxisAngle();
        expect(aa.axis.size).toBe(3);
        expect(aa.axis.values).toEqual([0, 0, 0]);
        expect(aa.angle).toBe(0);
    });

    it('stores a 3D axis and angle', () => {
        const axis = Vector.fromArray([0, 0, 1]);
        const aa = new AxisAngle(axis, Math.PI / 2);
        expect(aa.axis.values).toEqual([0, 0, 1]);
        expect(aa.angle).toBe(Math.PI / 2);
    });

    it('stores a 4D affine axis (x,y,z,0)', () => {
        const axis = Vector.fromArray([1, 0, 0, 0]);
        const aa = new AxisAngle(axis, 0.25);
        expect(aa.axis.size).toBe(4);
        expect(aa.axis.values).toEqual([1, 0, 0, 0]);
        expect(aa.angle).toBe(0.25);
    });

    it('copies the axis (C++ value semantics)', () => {
        const axis = Vector.fromArray([0, 1, 0]);
        const aa = new AxisAngle(axis, 1);
        axis.set(0, 99);
        expect(aa.axis.values).toEqual([0, 1, 0]);
    });

    it('throws for dimensions other than 3 or 4', () => {
        expect(() => new AxisAngle(Vector.fromArray([1, 0]), 1))
            .toThrow('Dimension must be 3 or 4.');
        expect(() => new AxisAngle(new Vector(5), 1))
            .toThrow('Dimension must be 3 or 4.');
    });

    it('fields are mutable, as the upstream public members', () => {
        const aa = new AxisAngle();
        aa.axis = Vector.fromArray([1, 0, 0, 0]);
        aa.angle = Math.PI;
        expect(aa.axis.values).toEqual([1, 0, 0, 0]);
        expect(aa.angle).toBe(Math.PI);
    });
});
