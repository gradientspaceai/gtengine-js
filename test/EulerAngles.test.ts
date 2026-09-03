import { describe, it, expect } from 'vitest';
import { check, finite, fc } from './helpers/arbitraries.js';
import { EulerAngles, EulerResult } from '../src/EulerAngles.js';

describe('EulerResult', () => {
    it('has the four upstream enumerators in declaration order', () => {
        expect(EulerResult.INVALID).toBe(0);
        expect(EulerResult.UNIQUE).toBe(1);
        expect(EulerResult.NOT_UNIQUE_SUM).toBe(2);
        expect(EulerResult.NOT_UNIQUE_DIF).toBe(3);
    });
});

describe('EulerAngles construction', () => {
    it('default constructor zeroes members and is INVALID', () => {
        const ea = new EulerAngles();
        expect(ea.axis).toEqual([0, 0, 0]);
        expect(ea.angle).toEqual([0, 0, 0]);
        expect(ea.result).toBe(EulerResult.INVALID);
    });

    it('6-argument constructor stores axes/angles and is UNIQUE', () => {
        const ea = new EulerAngles(2, 0, 1, 0.1, -0.2, 0.3);
        expect(ea.axis).toEqual([2, 0, 1]);
        expect(ea.angle).toEqual([0.1, -0.2, 0.3]);
        expect(ea.result).toBe(EulerResult.UNIQUE);
    });

    it('members are mutable, matching the upstream public fields', () => {
        const ea = new EulerAngles();
        ea.axis = [0, 1, 0];
        ea.angle = [Math.PI, 0, -Math.PI];
        ea.result = EulerResult.NOT_UNIQUE_SUM;
        expect(ea.axis).toEqual([0, 1, 0]);
        expect(ea.angle[0]).toBe(Math.PI);
        expect(ea.result).toBe(EulerResult.NOT_UNIQUE_SUM);
    });
});

// ---------------------------------------------------------------------------
// Verification wave: property-based checks against upstream EulerAngles.h.
// ---------------------------------------------------------------------------

describe('EulerAngles verification', () => {
    it('the six-argument constructor stores the axes/angles and UNIQUE', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 2 }),
            fc.integer({ min: 0, max: 2 }), fc.integer({ min: 0, max: 2 }),
            finite(), finite(), finite()), ([i0, i1, i2, a0, a1, a2]) => {
            const e = new EulerAngles(i0, i1, i2, a0, a1, a2);
            expect(e.axis).toEqual([i0, i1, i2]);
            expect(e.angle).toEqual([a0, a1, a2]);
            // Upstream sets result(EulerResult::UNIQUE) unconditionally in
            // this constructor; it does not validate the axis triple.
            expect(e.result).toBe(EulerResult.UNIQUE);
        });
    });

    it('the default constructor is the INVALID zero state', () => {
        const e = new EulerAngles();
        expect(e.axis).toEqual([0, 0, 0]);
        expect(e.angle).toEqual([0, 0, 0]);
        expect(e.result).toBe(EulerResult.INVALID);
    });

    it('instances do not share their axis/angle arrays', () => {
        // C++ std::array members are per-object values; the port must not
        // hand out a shared array.
        const a = new EulerAngles();
        const b = new EulerAngles();
        a.axis[0] = 2;
        a.angle[1] = 1.5;
        expect(b.axis).toEqual([0, 0, 0]);
        expect(b.angle).toEqual([0, 0, 0]);

        const c = new EulerAngles(0, 1, 2, 0.1, 0.2, 0.3);
        const d = new EulerAngles(0, 1, 2, 0.1, 0.2, 0.3);
        c.angle[2] = 9;
        expect(d.angle[2]).toBe(0.3);
    });

    it('the enumerators keep their upstream declaration order', () => {
        // Order matters: conversion code compares against these values.
        expect([EulerResult.INVALID, EulerResult.UNIQUE,
            EulerResult.NOT_UNIQUE_SUM, EulerResult.NOT_UNIQUE_DIF])
            .toEqual([0, 1, 2, 3]);
    });
});
