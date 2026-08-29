import { describe, it, expect } from 'vitest';
import {
    isArbitraryPrecision, hasDivisionOperator,
    type ArbitraryPrecisionNumber
} from '../src/TypeTraits';

// Stand-ins for the future BSNumber (no division) and BSRational (division)
// arbitrary-precision types.
class MockBSNumber implements ArbitraryPrecisionNumber {
    readonly isArbitraryPrecision = true as const;
    readonly hasDivisionOperator = false;
}

class MockBSRational implements ArbitraryPrecisionNumber {
    readonly isArbitraryPrecision = true as const;
    readonly hasDivisionOperator = true;
}

describe('TypeTraits', () => {
    it('numbers are not arbitrary precision (is_arbitrary_precision<double> == false)', () => {
        expect(isArbitraryPrecision(0)).toBe(false);
        expect(isArbitraryPrecision(1.5)).toBe(false);
        expect(isArbitraryPrecision(Number.MAX_VALUE)).toBe(false);
        expect(isArbitraryPrecision(NaN)).toBe(false);
    });

    it('numbers support division (has_division_operator<double> == true)', () => {
        expect(hasDivisionOperator(0)).toBe(true);
        expect(hasDivisionOperator(-2.25)).toBe(true);
    });

    it('marker-implementing values are arbitrary precision', () => {
        expect(isArbitraryPrecision(new MockBSNumber())).toBe(true);
        expect(isArbitraryPrecision(new MockBSRational())).toBe(true);
    });

    it('arbitrary-precision values report their own division support', () => {
        expect(hasDivisionOperator(new MockBSNumber())).toBe(false);
        expect(hasDivisionOperator(new MockBSRational())).toBe(true);
    });

    it('unrelated values are neither arbitrary precision nor divisible', () => {
        expect(isArbitraryPrecision(null)).toBe(false);
        expect(isArbitraryPrecision(undefined)).toBe(false);
        expect(isArbitraryPrecision('1.5')).toBe(false);
        expect(isArbitraryPrecision({})).toBe(false);
        expect(isArbitraryPrecision({ isArbitraryPrecision: false })).toBe(false);
        expect(hasDivisionOperator(null)).toBe(false);
        expect(hasDivisionOperator({})).toBe(false);
        expect(hasDivisionOperator('3')).toBe(false);
    });
});
