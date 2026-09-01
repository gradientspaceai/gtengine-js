import { describe, it, expect } from 'vitest';
import {
    BSNumber,
    BSNumberRoundingMode,
    convertBSNumber,
    BSRational,
    convertBSRational,
    convertBSRationalToBSNumber,
    convertBSRationalToNumber,
    convertBSRationalToFloat32,
    BSPrecision,
    BSPrecisionParameters,
    BSPrecisionType,
} from '../src/ArbitraryPrecision';

// The umbrella module must expose the same working set as the direct modules.
import * as Umbrella from '../src/ArbitraryPrecision';
import { BSNumber as DirectBSNumber } from '../src/BSNumber';
import { BSRational as DirectBSRational } from '../src/BSRational';
import { BSPrecision as DirectBSPrecision } from '../src/BSPrecision';

describe('ArbitraryPrecision (umbrella module)', () => {
    it('re-exports exactly the symbols of the ported includes', () => {
        expect(Object.keys(Umbrella).sort()).toEqual([
            'BSNumber',
            'BSNumberRoundingMode',
            'BSPrecision',
            'BSPrecisionParameters',
            'BSPrecisionType',
            'BSRational',
            'convertBSNumber',
            'convertBSRational',
            'convertBSRationalToBSNumber',
            'convertBSRationalToFloat32',
            'convertBSRationalToNumber',
        ]);
    });

    it('re-exports the identical bindings, not copies', () => {
        expect(BSNumber).toBe(DirectBSNumber);
        expect(BSRational).toBe(DirectBSRational);
        expect(BSPrecision).toBe(DirectBSPrecision);
    });

    it('exports the enums with their upstream values', () => {
        expect(BSNumberRoundingMode.FE_TONEAREST).toBe(0);
        expect(BSNumberRoundingMode.FE_UPWARD).toBe(3);
        expect(BSPrecisionType.IS_DOUBLE).toBe(1);
    });

    it('does not re-export the omitted UInteger*32 backends', () => {
        const names = Object.keys(Umbrella);
        expect(names.some(n => n.startsWith('UInteger'))).toBe(false);
    });
});

describe('ArbitraryPrecision: BSNumber through the umbrella', () => {
    it('performs exact arithmetic on values a double cannot hold', () => {
        // 2^60 + 1 is not representable in binary64; BSNumber is exact.
        const big = BSNumber.fromBigInt(1n << 60n);
        const one = BSNumber.fromNumber(1);
        const sum = big.add(one);
        expect(sum.getUInteger()).toBe((1n << 60n) + 1n);
        expect(sum.getSign()).toBe(1);
        expect(sum.sub(one).equals(big)).toBe(true);
    });

    it('multiplies exactly and compares', () => {
        const a = BSNumber.fromNumber(3);
        const b = BSNumber.fromNumber(-1.5);
        const p = a.mul(b);
        expect(p.toNumber()).toBe(-4.5);
        expect(p.lessThan(BSNumber.fromNumber(0))).toBe(true);
        expect(p.greaterThan(BSNumber.fromNumber(-5))).toBe(true);
    });

    it('handles the zero and negation degeneracies', () => {
        const zero = new BSNumber();
        expect(zero.getSign()).toBe(0);
        expect(zero.toNumber()).toBe(0);
        expect(zero.add(BSNumber.fromNumber(7)).toNumber()).toBe(7);
        expect(BSNumber.fromNumber(7).negated().toNumber()).toBe(-7);
    });

    it('rounds with convertBSNumber', () => {
        // 3 bits of precision cannot hold 2^4 + 1 = 10001b; round to nearest.
        const x = BSNumber.fromBigInt(17n);
        const rounded = convertBSNumber(x, 3, BSNumberRoundingMode.FE_TONEAREST);
        expect(rounded.toNumber()).toBe(16);
        const up = convertBSNumber(x, 3, BSNumberRoundingMode.FE_UPWARD);
        expect(up.toNumber()).toBe(20);
    });
});

describe('ArbitraryPrecision: BSRational through the umbrella', () => {
    it('adds fractions exactly (1/3 + 1/6 = 1/2)', () => {
        const third = BSRational.fromNumber(1, 3);
        const sixth = BSRational.fromNumber(1, 6);
        const sum = third.add(sixth);
        expect(sum.equals(BSRational.fromNumber(1, 2))).toBe(true);
        expect(sum.toNumber()).toBe(0.5);
    });

    it('divides and multiplies as inverse operations', () => {
        const a = BSRational.fromNumber(7, 11);
        const b = BSRational.fromNumber(-3, 5);
        expect(a.div(b).mul(b).equals(a)).toBe(true);
        expect(a.mul(b).toNumber()).toBeCloseTo(-21 / 55, 15);
    });

    it('orders fractions that round to the same double', () => {
        // 1/3 and (1/3 + 1/2^80) are the same binary64 value but distinct
        // rationals.
        const third = BSRational.fromNumber(1, 3);
        const eps = BSRational.fromBSNumber(BSNumber.fromNumber(1),
            BSNumber.fromBigInt(1n << 80n));
        const bigger = third.add(eps);
        expect(third.lessThan(bigger)).toBe(true);
        expect(bigger.toNumber()).toBe(third.toNumber());
    });

    it('rejects division by zero', () => {
        expect(() => BSRational.fromNumber(1, 0)).toThrow();
        expect(() => BSRational.fromNumber(1, 2).div(new BSRational())).toThrow();
    });

    it('converts through the umbrella conversion helpers', () => {
        const r = BSRational.fromNumber(1, 4);
        expect(convertBSRationalToNumber(r, BSNumberRoundingMode.FE_TONEAREST)).toBe(0.25);
        expect(convertBSRationalToFloat32(r, BSNumberRoundingMode.FE_TONEAREST)).toBe(0.25);
        const asNumber = convertBSRationalToBSNumber(r, 8,
            BSNumberRoundingMode.FE_TONEAREST);
        expect(asNumber.toNumber()).toBe(0.25);
        const asRational = convertBSRational(r, 8, BSNumberRoundingMode.FE_TONEAREST);
        expect(asRational.toNumber()).toBe(0.25);
    });
});

describe('ArbitraryPrecision: BSPrecision through the umbrella', () => {
    it('reports the double parameters', () => {
        const p = new BSPrecision(BSPrecisionType.IS_DOUBLE);
        expect(p.bsn.minExponent).toBe(-1074);
        expect(p.bsn.maxExponent).toBe(1023);
        expect(p.bsn.maxBits).toBe(53);
    });

    it('propagates precision through a sum', () => {
        const a = new BSPrecision(BSPrecisionType.IS_FLOAT);
        const sum = a.add(a);
        expect(sum.bsn.maxBits).toBeGreaterThanOrEqual(a.bsn.maxBits);
        expect(sum.bsn.minExponent).toBe(a.bsn.minExponent);
    });

    it('computes maxWords from maxBits', () => {
        expect(new BSPrecisionParameters(0, 0, 0).maxWords).toBe(0);
        expect(new BSPrecisionParameters(0, 0, 32).maxWords).toBe(1);
        expect(new BSPrecisionParameters(0, 0, 33).maxWords).toBe(2);
    });
});
