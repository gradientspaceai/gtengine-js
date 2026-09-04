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
} from '../src/ArbitraryPrecision.js';

// The umbrella module must expose the same working set as the direct modules.
import * as Umbrella from '../src/ArbitraryPrecision.js';
import { BSNumber as DirectBSNumber } from '../src/BSNumber.js';
import { BSRational as DirectBSRational } from '../src/BSRational.js';
import { BSPrecision as DirectBSPrecision } from '../src/BSPrecision.js';
import * as BSNumberModule from '../src/BSNumber.js';
import * as BSRationalModule from '../src/BSRational.js';
import * as BSPrecisionModule from '../src/BSPrecision.js';
import { check, fc } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). Upstream ArbitraryPrecision.h
// is a pure umbrella header: it includes UIntegerALU32.h, UIntegerAP32.h,
// UIntegerFP32.h, BSNumber.h, BSRational.h and BSPrecision.h and declares
// nothing of its own. The port omits the three UInteger*32 storage backends
// (the significands are bigint), so the module must re-export the complete
// surface of exactly BSNumber.ts, BSRational.ts and BSPrecision.ts.
// ---------------------------------------------------------------------------

describe('ArbitraryPrecision verification', () => {
    it('re-exports every binding of the three included modules and no more',
        () => {
            const expected = new Set([
                ...Object.keys(BSNumberModule),
                ...Object.keys(BSRationalModule),
                ...Object.keys(BSPrecisionModule)
            ]);
            const actual = new Set(Object.keys(Umbrella));
            expect([...actual].sort()).toEqual([...expected].sort());
            const umbrella = Umbrella as Record<string, unknown>;
            for (const name of actual) {
                const source =
                    (BSNumberModule as Record<string, unknown>)[name]
                    ?? (BSRationalModule as Record<string, unknown>)[name]
                    ?? (BSPrecisionModule as Record<string, unknown>)[name];
                expect(umbrella[name]).toBe(source);
            }
        });

    it('is transparent: the umbrella bindings compute the same values', () => {
        const smallInt = fc.integer({ min: -200, max: 200 });
        const denominator = fc.integer({ min: 1, max: 40 });
        check(fc.tuple(smallInt, denominator, smallInt, denominator),
            ([n0, d0, n1, d1]) => {
                const a = BSRational.fromNumber(n0, d0);
                const b = BSRational.fromNumber(n1, d1);
                const direct = DirectBSRational.fromNumber(n0, d0)
                    .add(DirectBSRational.fromNumber(n1, d1));
                expect(a.add(b).equals(direct)).toBe(true);

                // The exact sum, cross-checked against bigint arithmetic.
                const num = BigInt(n0) * BigInt(d1) + BigInt(n1) * BigInt(d0);
                const den = BigInt(d0) * BigInt(d1);
                expect(a.add(b).mul(BSRational.fromBigInt(den))
                    .equals(BSRational.fromBigInt(num))).toBe(true);

                // BSNumber multiplication of the integer numerators is exact
                // and both products fit a double, so toNumber is lossless.
                // (BSNumber has a single unsigned zero, so -0 becomes 0.)
                expect(BSNumber.fromBigInt(BigInt(n0))
                    .mul(BSNumber.fromBigInt(BigInt(n1))).toNumber())
                    .toBe(n0 * n1 === 0 ? 0 : n0 * n1);
            });
    });

    it('routes the conversion helpers to the same rounding behavior', () => {
        check(fc.tuple(fc.integer({ min: -500, max: 500 }),
            fc.integer({ min: 1, max: 60 }),
            fc.integer({ min: 1, max: 20 })), ([n, d, precision]) => {
                const x = BSRational.fromNumber(n, d);
                const down = convertBSRationalToBSNumber(x, precision,
                    BSNumberRoundingMode.FE_DOWNWARD);
                const up = convertBSRationalToBSNumber(x, precision,
                    BSNumberRoundingMode.FE_UPWARD);
                // The rational-output overload wraps the BSNumber overload.
                expect(convertBSRational(x, precision,
                    BSNumberRoundingMode.FE_TONEAREST)
                    .equals(BSRational.fromBSNumber(
                        convertBSRationalToBSNumber(x, precision,
                            BSNumberRoundingMode.FE_TONEAREST)))).toBe(true);
                // The directed results bracket the exact value.
                expect(BSRational.fromBSNumber(down).lessThanOrEqual(x))
                    .toBe(true);
                expect(x.lessThanOrEqual(BSRational.fromBSNumber(up)))
                    .toBe(true);
                // The double and float conversions are the 53- and 24-bit
                // instances of the same routine.
                expect(convertBSRationalToNumber(x,
                    BSNumberRoundingMode.FE_TONEAREST))
                    .toBe(convertBSRationalToBSNumber(x, 53,
                        BSNumberRoundingMode.FE_TONEAREST).toNumber());
                expect(convertBSRationalToFloat32(x,
                    BSNumberRoundingMode.FE_TONEAREST))
                    .toBe(convertBSRationalToBSNumber(x, 24,
                        BSNumberRoundingMode.FE_TONEAREST).toFloat32());
                // convertBSNumber at full double precision is the identity
                // on a number that came from a double.
                const asNumber = BSNumber.fromNumber(n / d);
                expect(convertBSNumber(asNumber, 53,
                    BSNumberRoundingMode.FE_TONEAREST)
                    .equals(asNumber)).toBe(true);
            });
    });

    it('propagates BSPrecision through the umbrella as the direct module does',
        () => {
            const parameters = fc.tuple(fc.integer({ min: -60, max: 0 }),
                fc.integer({ min: 0, max: 60 }), fc.integer({ min: 1, max: 60 }));
            check(fc.tuple(parameters, parameters), ([p0, p1]) => {
                const build = (Ctor: typeof BSPrecision,
                    p: [number, number, number]) => new Ctor(p[0], p[1], p[2]);
                for (const op of ['add', 'sub', 'mul'] as const) {
                    const viaUmbrella =
                        build(BSPrecision, p0)[op](build(BSPrecision, p1));
                    const direct = build(DirectBSPrecision, p0)[op](
                        build(DirectBSPrecision, p1));
                    expect(viaUmbrella.bsn.minExponent)
                        .toBe(direct.bsn.minExponent);
                    expect(viaUmbrella.bsn.maxExponent)
                        .toBe(direct.bsn.maxExponent);
                    expect(viaUmbrella.bsn.maxBits).toBe(direct.bsn.maxBits);
                    // maxWords is derived from maxBits by the same formula.
                    expect(viaUmbrella.bsn.maxWords)
                        .toBe(new BSPrecisionParameters(0, 0,
                            viaUmbrella.bsn.maxBits).maxWords);
                }
            });
        });
});
