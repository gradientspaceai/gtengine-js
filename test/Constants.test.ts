import { describe, it, expect } from 'vitest';
import {
    GTE_C_PI, GTE_C_HALF_PI, GTE_C_QUARTER_PI, GTE_C_TWO_PI,
    GTE_C_INV_PI, GTE_C_INV_TWO_PI, GTE_C_INV_HALF_PI,
    GTE_C_DEG_TO_RAD, GTE_C_RAD_TO_DEG,
    GTE_C_SQRT_2, GTE_C_INV_SQRT_2,
    GTE_C_LN_2, GTE_C_INV_LN_2, GTE_C_LN_10, GTE_C_INV_LN_10
} from '../src/Constants.js';

describe('Constants', () => {
    it('pi constants match the exact double values', () => {
        expect(GTE_C_PI).toBe(Math.PI);
        expect(GTE_C_HALF_PI).toBe(Math.PI / 2);
        expect(GTE_C_QUARTER_PI).toBe(Math.PI / 4);
        expect(GTE_C_TWO_PI).toBe(2 * Math.PI);
        expect(GTE_C_INV_PI).toBe(1 / Math.PI);
        expect(GTE_C_INV_TWO_PI).toBe(1 / (2 * Math.PI));
        expect(GTE_C_INV_HALF_PI).toBe(2 / Math.PI);
    });

    it('degree/radian conversions match the exact double values', () => {
        expect(GTE_C_DEG_TO_RAD).toBe(Math.PI / 180);
        expect(GTE_C_RAD_TO_DEG).toBe(180 / Math.PI);
        // Round trips are exact for these particular constants at 90 deg.
        expect(90 * GTE_C_DEG_TO_RAD).toBeCloseTo(GTE_C_HALF_PI, 15);
    });

    it('common constants match the exact double values', () => {
        expect(GTE_C_SQRT_2).toBe(Math.sqrt(2));
        expect(GTE_C_LN_2).toBe(Math.LN2);
        expect(GTE_C_INV_LN_2).toBe(Math.LOG2E);
        expect(GTE_C_LN_10).toBe(Math.LN10);
    });

    it('keeps the upstream literals where they differ from the Math constants', () => {
        // GTE's 1/sqrt(2) is the correctly rounded 1/Math.sqrt(2), which is
        // one ulp below Math.SQRT1_2.
        expect(GTE_C_INV_SQRT_2).toBe(1 / Math.sqrt(2));
        expect(GTE_C_INV_SQRT_2).not.toBe(Math.SQRT1_2);
        expect(GTE_C_INV_SQRT_2).toBe(0.7071067811865475);

        expect(GTE_C_INV_LN_10).toBe(0.43429448190325176);
        expect(Math.abs(GTE_C_INV_LN_10 - Math.LOG10E)).toBeLessThanOrEqual(Number.EPSILON);
    });

    it('reciprocal identities hold to double precision', () => {
        expect(GTE_C_PI * GTE_C_INV_PI).toBeCloseTo(1, 15);
        expect(GTE_C_TWO_PI * GTE_C_INV_TWO_PI).toBeCloseTo(1, 15);
        expect(GTE_C_HALF_PI * GTE_C_INV_HALF_PI).toBeCloseTo(1, 15);
        expect(GTE_C_SQRT_2 * GTE_C_INV_SQRT_2).toBeCloseTo(1, 15);
        expect(GTE_C_LN_2 * GTE_C_INV_LN_2).toBeCloseTo(1, 15);
        expect(GTE_C_LN_10 * GTE_C_INV_LN_10).toBeCloseTo(1, 15);
        expect(GTE_C_DEG_TO_RAD * GTE_C_RAD_TO_DEG).toBeCloseTo(1, 15);
    });
});
