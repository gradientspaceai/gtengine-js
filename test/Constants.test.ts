import { describe, it, expect } from 'vitest';
import { check, finite, expectClose, wellScaled } from './helpers/arbitraries.js';
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

describe('Constants verification', () => {
    // Every literal in the port must be byte-identical to the upstream literal.
    // Upstream Constants.h (File Version 8.0.2026.08.01) declares these values.
    const upstream: ReadonlyArray<readonly [number, number]> = [
        [GTE_C_PI, 3.1415926535897931],
        [GTE_C_HALF_PI, 1.5707963267948966],
        [GTE_C_QUARTER_PI, 0.7853981633974483],
        [GTE_C_TWO_PI, 6.2831853071795862],
        [GTE_C_INV_PI, 0.3183098861837907],
        [GTE_C_INV_TWO_PI, 0.15915494309189535],
        [GTE_C_INV_HALF_PI, 0.63661977236758138],
        [GTE_C_DEG_TO_RAD, 0.017453292519943295],
        [GTE_C_RAD_TO_DEG, 57.295779513082321],
        [GTE_C_SQRT_2, 1.4142135623730951],
        [GTE_C_INV_SQRT_2, 0.7071067811865475],
        [GTE_C_LN_2, 0.6931471805599453],
        [GTE_C_INV_LN_2, 1.4426950408889634],
        [GTE_C_LN_10, 2.3025850929940459],
        [GTE_C_INV_LN_10, 0.43429448190325176]
    ];

    it('every constant is bit-identical to the upstream literal', () => {
        const bits = (x: number) => {
            const dv = new DataView(new ArrayBuffer(8));
            dv.setFloat64(0, x);
            return dv.getBigUint64(0);
        };
        for (const [ported, literal] of upstream) {
            expect(bits(ported)).toBe(bits(literal));
        }
    });

    it('every constant is within one ulp of the mathematical value', () => {
        // Independent computation of each constant; the upstream literals are
        // correctly rounded, so they agree to within one ulp of these.
        const ulp = (x: number) => Math.abs(x) * Number.EPSILON;
        const pairs: ReadonlyArray<readonly [number, number]> = [
            [GTE_C_PI, Math.PI],
            [GTE_C_HALF_PI, Math.PI / 2],
            [GTE_C_QUARTER_PI, Math.PI / 4],
            [GTE_C_TWO_PI, 2 * Math.PI],
            [GTE_C_INV_PI, 1 / Math.PI],
            [GTE_C_INV_TWO_PI, 1 / (2 * Math.PI)],
            [GTE_C_INV_HALF_PI, 2 / Math.PI],
            [GTE_C_DEG_TO_RAD, Math.PI / 180],
            [GTE_C_RAD_TO_DEG, 180 / Math.PI],
            [GTE_C_SQRT_2, Math.SQRT2],
            [GTE_C_INV_SQRT_2, Math.SQRT1_2],
            [GTE_C_LN_2, Math.LN2],
            [GTE_C_INV_LN_2, Math.LOG2E],
            [GTE_C_LN_10, Math.LN10],
            [GTE_C_INV_LN_10, Math.LOG10E]
        ];
        for (const [ported, reference] of pairs) {
            expect(Math.abs(ported - reference)).toBeLessThanOrEqual(ulp(reference));
        }
    });

    it('degree-radian conversion agrees with pi/180 for random angles', () => {
        // wellScaled: a subnormal deg makes an ulp-relative bound meaningless.
        check(wellScaled(-1e4, 1e4), deg => {
            expectClose(deg * GTE_C_DEG_TO_RAD, deg * Math.PI / 180, 0, 4 * Number.EPSILON);
        });
    });

    it('degree-radian conversion round trips', () => {
        // Two roundings of a product; a few ulps of relative error is expected.
        check(wellScaled(-1e4, 1e4), deg => {
            expectClose(deg * GTE_C_DEG_TO_RAD * GTE_C_RAD_TO_DEG, deg, 1e-12,
                8 * Number.EPSILON);
        });
    });

    it('reciprocal constants agree with division for random values', () => {
        check(wellScaled(-1e3, 1e3), x => {
            expectClose(x * GTE_C_INV_PI, x / Math.PI, 0, 4 * Number.EPSILON);
            expectClose(x * GTE_C_INV_TWO_PI, x / (2 * Math.PI), 0, 4 * Number.EPSILON);
            expectClose(x * GTE_C_INV_HALF_PI, x / (Math.PI / 2), 0, 4 * Number.EPSILON);
            expectClose(x * GTE_C_INV_SQRT_2, x / Math.sqrt(2), 0, 4 * Number.EPSILON);
        });
    });

    it('log constants convert bases for random positive values', () => {
        check(finite(1e-6, 1e6).filter(x => x > 0), x => {
            expectClose(Math.log(x) * GTE_C_INV_LN_2, Math.log2(x), 1e-12, 1e-14);
            expectClose(Math.log(x) * GTE_C_INV_LN_10, Math.log10(x), 1e-12, 1e-14);
            expectClose(Math.log2(x) * GTE_C_LN_2, Math.log(x), 1e-12, 1e-14);
            expectClose(Math.log10(x) * GTE_C_LN_10, Math.log(x), 1e-12, 1e-14);
        });
    });
});
