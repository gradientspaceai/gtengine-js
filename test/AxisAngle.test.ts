import { describe, it, expect } from 'vitest';
import { check, finite, fc } from './helpers/arbitraries.js';
import { AxisAngle } from '../src/AxisAngle.js';
import { Vector } from '../src/Vector.js';

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

// ---------------------------------------------------------------------------
// Verification wave: property-based checks against upstream AxisAngle.h.
// ---------------------------------------------------------------------------

describe('AxisAngle verification', () => {
    it('stores a copy of the axis and the angle for N = 3 and N = 4', () => {
        check(fc.tuple(fc.constantFrom(3, 4), finite()), ([n, angle]) => {
            const axis = Vector.fromArray(
                Array.from({ length: n }, (_, i) => i + 1));
            const aa = new AxisAngle(axis, angle);
            expect(aa.axis.size).toBe(n);
            expect(aa.axis.values).toEqual(axis.values);
            expect(aa.angle).toBe(angle);
            // Upstream copies the axis into the member (C++ value semantics),
            // so later writes to the input must not be visible.
            axis.set(0, 12345);
            expect(aa.axis.get(0)).toBe(1);
            // ... nor must writes to the member change the input.
            aa.axis.set(1, -7);
            expect(axis.get(1)).toBe(2);
        });
    });

    it('the default constructor is the 3D zero axis with zero angle', () => {
        const aa = new AxisAngle();
        expect(aa.axis.size).toBe(3);
        expect(aa.axis.values).toEqual([0, 0, 0]);
        expect(aa.angle).toBe(0);
        // Distinct instances must not share the default axis.
        const bb = new AxisAngle();
        aa.axis.set(0, 1);
        expect(bb.axis.get(0)).toBe(0);
    });

    it('rejects every dimension other than 3 and 4', () => {
        // Upstream: static_assert(N == 3 || N == 4).
        check(fc.integer({ min: 0, max: 8 }).filter(n => n !== 3 && n !== 4),
            n => {
                expect(() => new AxisAngle(new Vector(n), 1))
                    .toThrow('Dimension must be 3 or 4.');
            });
    });
});
