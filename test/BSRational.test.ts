import { describe, it, expect } from 'vitest';
import {
    BSRational, convertBSRational, convertBSRationalToBSNumber,
    convertBSRationalToFloat32, convertBSRationalToNumber
} from '../src/BSRational.js';
import { BSNumber, BSNumberRoundingMode } from '../src/BSNumber.js';
import { isArbitraryPrecision, hasDivisionOperator } from '../src/TypeTraits.js';
import { check, fc } from './helpers/arbitraries.js';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// ---------------------------------------------------------------------------
// An independent exact-rational arithmetic implemented with plain bigint
// numerator/denominator pairs. Every BSRational result is cross-checked
// against this.
// ---------------------------------------------------------------------------

type Rat = { n: bigint, d: bigint };  // d > 0, not necessarily reduced

function ratMake(n: bigint, d: bigint): Rat {
    if (d < 0n) {
        return { n: -n, d: -d };
    }
    return { n, d };
}

function ratAdd(a: Rat, b: Rat): Rat {
    return ratMake(a.n * b.d + a.d * b.n, a.d * b.d);
}

function ratSub(a: Rat, b: Rat): Rat {
    return ratMake(a.n * b.d - a.d * b.n, a.d * b.d);
}

function ratMul(a: Rat, b: Rat): Rat {
    return ratMake(a.n * b.n, a.d * b.d);
}

function ratDiv(a: Rat, b: Rat): Rat {
    return ratMake(a.n * b.d, a.d * b.n);
}

// -1, 0 or +1 for a < b, a == b, a > b.
function ratCompare(a: Rat, b: Rat): number {
    const left = a.n * b.d;
    const right = a.d * b.n;
    return left < right ? -1 : (left > right ? 1 : 0);
}

// The exact value of a BSNumber as num * 2^exp.
function bsnExact(x: BSNumber): { num: bigint, exp: number } {
    const magnitude = x.getUInteger();
    return { num: (x.getSign() < 0 ? -magnitude : magnitude), exp: x.getBiasedExponent() };
}

// The exact value of a BSRational as a bigint fraction.
function exact(x: BSRational): Rat {
    const n = bsnExact(x.getNumerator());
    const d = bsnExact(x.getDenominator());
    // (n.num * 2^n.exp) / (d.num * 2^d.exp). Move the power of two into
    // whichever of the numerator or denominator keeps both integral.
    const e = n.exp - d.exp;
    if (e >= 0) {
        return ratMake(n.num << BigInt(e), d.num);
    }
    return ratMake(n.num, d.num << BigInt(-e));
}

// Every nonzero BSNumber must have an odd positive unsigned integer; zero
// must have all members zero.
function expectValidBSNumber(x: BSNumber): void {
    expect(x.isValid()).toBe(true);
}

function expectValid(x: BSRational): void {
    expectValidBSNumber(x.getNumerator());
    expectValidBSNumber(x.getDenominator());
    // The denominator is always positive.
    expect(x.getDenominator().getSign()).toBe(1);
    // The denominator is normalized to have exponent zero, i.e. it lies in
    // the interval [1,2).
    expect(x.getDenominator().getExponent()).toBe(0);
}

function expectExactlyEqual(x: BSRational, r: Rat): void {
    const e = exact(x);
    expect(ratCompare(e, r)).toBe(0);
}

// ---------------------------------------------------------------------------

describe('BSRational: construction', () => {
    it('default constructs zero', () => {
        const zero = new BSRational();
        expect(zero.getSign()).toBe(0);
        expect(zero.toNumber()).toBe(0);
        expectValid(zero);
        expect(zero.getNumerator().getSign()).toBe(0);
    });

    it('fromNumber with a single argument sets the denominator to one', () => {
        const x = BSRational.fromNumber(-3.5);
        expect(x.toNumber()).toBe(-3.5);
        expect(x.getDenominator().toNumber()).toBe(1);
        expect(x.getSign()).toBe(-1);
        expectValid(x);
    });

    it('fromNumber with a denominator produces the exact quotient', () => {
        const third = BSRational.fromNumber(1, 3);
        expectExactlyEqual(third, { n: 1n, d: 3n });
        expect(third.toNumber()).toBe(1 / 3);
        expectValid(third);
    });

    it('normalizes a negative denominator onto the numerator', () => {
        const x = BSRational.fromNumber(3, -4);
        expect(x.getSign()).toBe(-1);
        expect(x.getDenominator().getSign()).toBe(1);
        expect(x.toNumber()).toBe(-0.75);
        expectExactlyEqual(x, { n: -3n, d: 4n });
        expectValid(x);
    });

    it('rejects a zero denominator', () => {
        expect(() => BSRational.fromNumber(1, 0)).toThrow(/Division by zero/);
    });

    it('keeps a zero numerator canonical (upstream would corrupt it)', () => {
        // Upstream unconditionally subtracts the denominator exponent from
        // the numerator biased exponent, which for a zero numerator makes an
        // invalid BSNumber; the port guards it.
        const x = BSRational.fromNumber(0, 1024);
        expect(x.getSign()).toBe(0);
        expect(x.toNumber()).toBe(0);
        expectValid(x);
        expect(x.equals(new BSRational())).toBe(true);
    });

    it('fromBigInt handles integers wider than a double', () => {
        const big = (1n << 100n) + 1n;
        const x = BSRational.fromBigInt(big, 7n);
        expectExactlyEqual(x, { n: big, d: 7n });
        expectValid(x);
    });

    it('fromBSNumber copies its inputs', () => {
        const n = BSNumber.fromNumber(5);
        const d = BSNumber.fromNumber(2);
        const x = BSRational.fromBSNumber(n, d);
        n.setSign(-1);
        expect(x.toNumber()).toBe(2.5);
    });

    it('fromFloat32 rounds the inputs to binary32', () => {
        const x = BSRational.fromFloat32(0.1);
        expect(x.toNumber()).toBe(Math.fround(0.1));
        expect(x.toNumber()).not.toBe(0.1);
    });

    it('normalizes the denominator exponent to zero', () => {
        // 3 * 2^10 / (5 * 2^20). The denominator exponent must be zero and
        // the numerator absorbs the difference.
        const x = BSRational.fromNumber(3 * 1024, 5 * 1048576);
        expect(x.getDenominator().getExponent()).toBe(0);
        expectExactlyEqual(x, { n: 3n * 1024n, d: 5n * 1048576n });
    });
});

describe('BSRational: fromString', () => {
    it('parses integers', () => {
        expect(BSRational.fromString('123').toNumber()).toBe(123);
        expect(BSRational.fromString('+123').toNumber()).toBe(123);
        expect(BSRational.fromString('-123').toNumber()).toBe(-123);
    });

    it('parses "x.y" exactly', () => {
        const x = BSRational.fromString('1.5');
        expectExactlyEqual(x, { n: 3n, d: 2n });
        expect(x.toNumber()).toBe(1.5);
    });

    it('parses "x." as the integer x', () => {
        const x = BSRational.fromString('12.');
        expect(x.toNumber()).toBe(12);
        expectExactlyEqual(x, { n: 12n, d: 1n });
    });

    it('parses ".y"', () => {
        const x = BSRational.fromString('.25');
        expectExactlyEqual(x, { n: 1n, d: 4n });
        expect(x.toNumber()).toBe(0.25);
    });

    it('represents 0.1 exactly, unlike the double 0.1', () => {
        const x = BSRational.fromString('0.1');
        expectExactlyEqual(x, { n: 1n, d: 10n });
        // The nearest double to the exact 1/10 is the double literal 0.1.
        expect(x.toNumber()).toBe(0.1);
        // But the exact values differ.
        expect(x.equals(BSRational.fromNumber(0.1))).toBe(false);
    });

    it('parses a signed fraction', () => {
        const x = BSRational.fromString('-3.75');
        expectExactlyEqual(x, { n: -15n, d: 4n });
    });

    it('keeps a signed zero canonical', () => {
        const x = BSRational.fromString('-0.0');
        expect(x.getSign()).toBe(0);
        expectValid(x);
    });

    it('rejects malformed input', () => {
        expect(() => BSRational.fromString('')).toThrow(/A number must be specified/);
        expect(() => BSRational.fromString('+')).toThrow(/Invalid number format/);
        expect(() => BSRational.fromString('1.2a')).toThrow(/Invalid number format/);
    });

    it('agrees with exact rational arithmetic for many decimals', () => {
        const cases: Array<[string, bigint, bigint]> = [
            ['0.5', 1n, 2n],
            ['2.125', 17n, 8n],
            ['0.001', 1n, 1000n],
            ['12.34', 1234n, 100n],
            ['-7.007', -7007n, 1000n],
            ['100.0', 100n, 1n]
        ];
        for (const [s, n, d] of cases) {
            expectExactlyEqual(BSRational.fromString(s), { n, d });
        }
    });
});

describe('BSRational: comparisons', () => {
    it('compares numbers with unlike denominators', () => {
        const a = BSRational.fromNumber(1, 3);
        const b = BSRational.fromNumber(2, 6);
        expect(a.equals(b)).toBe(true);
        expect(a.notEquals(b)).toBe(false);
        expect(a.lessThan(b)).toBe(false);
        expect(a.lessThanOrEqual(b)).toBe(true);
        expect(a.greaterThanOrEqual(b)).toBe(true);
    });

    it('orders across signs', () => {
        const neg = BSRational.fromNumber(-1, 7);
        const zero = new BSRational();
        const pos = BSRational.fromNumber(1, 1000000);
        expect(neg.lessThan(zero)).toBe(true);
        expect(zero.lessThan(pos)).toBe(true);
        expect(neg.lessThan(pos)).toBe(true);
        expect(pos.lessThan(neg)).toBe(false);
        expect(zero.lessThan(zero)).toBe(false);
        expect(pos.greaterThan(neg)).toBe(true);
    });

    it('handles the left-aligned magnitude comparison of BSNumber', () => {
        // 3/4 and 11/16 have the same exponent for the cross-multiplied
        // numerators but different bit counts; the comparison must be the
        // left-aligned one that BSNumber implements.
        const a = BSRational.fromNumber(3, 4);
        const b = BSRational.fromNumber(11, 16);
        expect(a.greaterThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(true);
        expect(ratCompare(exact(a), exact(b))).toBe(1);
    });

    it('cross-checks against bigint rational comparison', () => {
        const rand = makeRandom(20260901);
        for (let trial = 0; trial < 300; ++trial) {
            const an = Math.trunc(rand() * 200) - 100;
            const ad = Math.trunc(rand() * 99) + 1;
            const bn = Math.trunc(rand() * 200) - 100;
            const bd = Math.trunc(rand() * 99) + 1;
            const a = BSRational.fromNumber(an, ad);
            const b = BSRational.fromNumber(bn, bd);
            const c = ratCompare(ratMake(BigInt(an), BigInt(ad)),
                ratMake(BigInt(bn), BigInt(bd)));
            expect(a.lessThan(b)).toBe(c < 0);
            expect(a.equals(b)).toBe(c === 0);
            expect(a.greaterThan(b)).toBe(c > 0);
            expect(a.lessThanOrEqual(b)).toBe(c <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(c >= 0);
        }
    });
});

describe('BSRational: arithmetic', () => {
    it('adds, subtracts, multiplies and divides exactly', () => {
        const a = BSRational.fromNumber(1, 3);
        const b = BSRational.fromNumber(1, 7);
        expectExactlyEqual(a.add(b), { n: 10n, d: 21n });
        expectExactlyEqual(a.sub(b), { n: 4n, d: 21n });
        expectExactlyEqual(a.mul(b), { n: 1n, d: 21n });
        expectExactlyEqual(a.div(b), { n: 7n, d: 3n });
    });

    it('returns a canonical zero when the result is zero', () => {
        const a = BSRational.fromNumber(1, 3);
        const diff = a.sub(a);
        expect(diff.getSign()).toBe(0);
        expect(diff.getDenominator().toNumber()).toBe(1);
        expectValid(diff);

        const prod = a.mul(new BSRational());
        expect(prod.getSign()).toBe(0);
        expect(prod.getDenominator().toNumber()).toBe(1);

        const quot = new BSRational().div(a);
        expect(quot.getSign()).toBe(0);
        expectValid(quot);
    });

    it('rejects division by zero', () => {
        expect(() => BSRational.fromNumber(1).div(new BSRational()))
            .toThrow(/Division by zero/);
    });

    it('normalizes the sign after division by a negative', () => {
        const q = BSRational.fromNumber(1, 2).div(BSRational.fromNumber(-1, 4));
        expect(q.getDenominator().getSign()).toBe(1);
        expect(q.toNumber()).toBe(-2);
        expectExactlyEqual(q, { n: -2n, d: 1n });
    });

    it('satisfies the field identities', () => {
        const rand = makeRandom(777);
        for (let trial = 0; trial < 60; ++trial) {
            const mk = (): BSRational => BSRational.fromNumber(
                Math.trunc(rand() * 40) - 20, Math.trunc(rand() * 19) + 1);
            const a = mk(), b = mk(), c = mk();

            // Commutativity and associativity.
            expect(a.add(b).equals(b.add(a))).toBe(true);
            expect(a.mul(b).equals(b.mul(a))).toBe(true);
            expect(a.add(b).add(c).equals(a.add(b.add(c)))).toBe(true);
            expect(a.mul(b).mul(c).equals(a.mul(b.mul(c)))).toBe(true);

            // Distributivity.
            expect(a.mul(b.add(c)).equals(a.mul(b).add(a.mul(c)))).toBe(true);

            // Subtraction is addition of the negation.
            expect(a.sub(b).equals(a.add(b.negated()))).toBe(true);

            // Division inverts multiplication.
            if (b.getSign() !== 0) {
                expect(a.div(b).mul(b).equals(a)).toBe(true);
            }
        }
    });

    it('cross-checks arithmetic against bigint rational arithmetic', () => {
        const rand = makeRandom(31415);
        for (let trial = 0; trial < 400; ++trial) {
            const an = Math.trunc(rand() * 400) - 200;
            const ad = Math.trunc(rand() * 200) + 1;
            const bn = Math.trunc(rand() * 400) - 200;
            const bd = Math.trunc(rand() * 200) + 1;
            const a = BSRational.fromNumber(an, ad);
            const b = BSRational.fromNumber(bn, bd);
            const ra = ratMake(BigInt(an), BigInt(ad));
            const rb = ratMake(BigInt(bn), BigInt(bd));

            expectExactlyEqual(a.add(b), ratAdd(ra, rb));
            expectExactlyEqual(a.sub(b), ratSub(ra, rb));
            expectExactlyEqual(a.mul(b), ratMul(ra, rb));
            if (bn !== 0) {
                expectExactlyEqual(a.div(b), ratDiv(ra, rb));
            }
            expectValid(a.add(b));
            expectValid(a.mul(b));
        }
    });

    it('accumulates a long exact sum that doubles cannot represent', () => {
        // sum_{k=1}^{20} 1/k as an exact rational.
        let sum = new BSRational();
        let ref: Rat = { n: 0n, d: 1n };
        for (let k = 1; k <= 20; ++k) {
            sum = sum.add(BSRational.fromNumber(1, k));
            ref = ratAdd(ref, ratMake(1n, BigInt(k)));
        }
        expectExactlyEqual(sum, ref);
        expect(sum.toNumber()).toBeCloseTo(Number(ref.n) / Number(ref.d), 12);
    });

    it('negates, in place and out of place', () => {
        const x = BSRational.fromNumber(5, 8);
        const y = x.negated();
        expect(y.toNumber()).toBe(-0.625);
        expect(x.toNumber()).toBe(0.625);
        x.negate();
        expect(x.toNumber()).toBe(-0.625);
        expect(x.getDenominator().getSign()).toBe(1);
    });

    it('setSign forces the numerator sign and a positive denominator', () => {
        const x = BSRational.fromNumber(3, 4);
        x.setSign(-1);
        expect(x.getSign()).toBe(-1);
        expect(x.toNumber()).toBe(-0.75);
    });

    it('clone is a deep copy', () => {
        const x = BSRational.fromNumber(1, 3);
        const y = x.clone();
        y.negate();
        expect(x.toNumber()).toBe(1 / 3);
        expect(y.toNumber()).toBe(-1 / 3);
    });

    it('fma, robustSOP and robustDOP are exact', () => {
        const u = BSRational.fromNumber(1, 3);
        const v = BSRational.fromNumber(3, 5);
        const w = BSRational.fromNumber(1, 7);
        const z = BSRational.fromNumber(7, 11);
        expectExactlyEqual(BSRational.fma(u, v, w),
            ratAdd(ratMul({ n: 1n, d: 3n }, { n: 3n, d: 5n }), { n: 1n, d: 7n }));
        expectExactlyEqual(BSRational.robustSOP(u, v, w, z),
            ratAdd(ratMul({ n: 1n, d: 3n }, { n: 3n, d: 5n }),
                ratMul({ n: 1n, d: 7n }, { n: 7n, d: 11n })));
        expectExactlyEqual(BSRational.robustDOP(u, v, w, z),
            ratSub(ratMul({ n: 1n, d: 3n }, { n: 3n, d: 5n }),
                ratMul({ n: 1n, d: 7n }, { n: 7n, d: 11n })));
    });

    it('computes an exact 2x2 determinant sign that doubles get wrong', () => {
        // The classic robust-predicate example: the double evaluation of
        // a*d - b*c cancels catastrophically.
        const a = BSRational.fromNumber(1e15 + 1);
        const b = BSRational.fromNumber(1e15);
        const c = BSRational.fromNumber(1e15);
        const d = BSRational.fromNumber(1e15 - 1);
        const det = BSRational.robustDOP(a, d, b, c);
        expect(det.getSign()).toBe(-1);
        expect(det.toNumber()).toBe(-1);
        // The double evaluation of the same expression loses the sign.
        expect((1e15 + 1) * (1e15 - 1) - 1e15 * 1e15).toBe(0);
    });
});

describe('BSRational: conversions to floating point', () => {
    it('round-trips exact doubles', () => {
        const values = [0, 1, -1, 0.5, -0.5, 3.25, 1e-300, 1e300,
            Number.MAX_SAFE_INTEGER, -1234.5678];
        for (const v of values) {
            expect(BSRational.fromNumber(v).toNumber()).toBe(v);
        }
    });

    it('rounds 1/3 and 1/10 to the nearest double', () => {
        expect(BSRational.fromNumber(1, 3).toNumber()).toBe(1 / 3);
        expect(BSRational.fromString('0.1').toNumber()).toBe(0.1);
        expect(BSRational.fromNumber(2, 3).toNumber()).toBe(2 / 3);
    });

    it('cross-checks random quotients against the double quotient', () => {
        const rand = makeRandom(98765);
        for (let trial = 0; trial < 300; ++trial) {
            // Use small integers so that the double quotient is the
            // correctly rounded value of the exact quotient.
            const n = Math.trunc(rand() * 2000) - 1000;
            const d = Math.trunc(rand() * 1000) + 1;
            expect(BSRational.fromNumber(n, d).toNumber()).toBe(n / d);
        }
    });

    it('toFloat32 matches Math.fround of the exact quotient', () => {
        const rand = makeRandom(2468);
        for (let trial = 0; trial < 200; ++trial) {
            const n = Math.trunc(rand() * 2000) - 1000;
            const d = Math.trunc(rand() * 1000) + 1;
            expect(BSRational.fromNumber(n, d).toFloat32()).toBe(Math.fround(n / d));
        }
    });
});

describe('BSRational: convert with an explicit precision and rounding mode', () => {
    it('reproduces the exact value when the precision suffices', () => {
        const x = BSRational.fromNumber(3, 4);  // 0.75 = 1.1 * 2^-1
        for (const mode of [BSNumberRoundingMode.FE_TONEAREST,
            BSNumberRoundingMode.FE_DOWNWARD, BSNumberRoundingMode.FE_UPWARD,
            BSNumberRoundingMode.FE_TOWARDZERO]) {
            const n = convertBSRationalToBSNumber(x, 8, mode);
            expect(n.toNumber()).toBe(0.75);
            expect(n.isValid()).toBe(true);
        }
    });

    it('converts zero to zero', () => {
        const n = convertBSRationalToBSNumber(new BSRational(), 10,
            BSNumberRoundingMode.FE_TONEAREST);
        expect(n.getSign()).toBe(0);
    });

    it('rejects a nonpositive precision', () => {
        expect(() => convertBSRationalToBSNumber(BSRational.fromNumber(1), 0,
            BSNumberRoundingMode.FE_TONEAREST)).toThrow(/Precision must be positive/);
    });

    it('rejects an unsupported rounding mode', () => {
        expect(() => convertBSRationalToBSNumber(BSRational.fromNumber(1, 3), 4,
            99 as BSNumberRoundingMode))
            .toThrow(/Implementation-dependent rounding mode not supported/);
    });

    it('rounds 1/3 to 4 bits in each mode', () => {
        // 1/3 = 0.01010101..._2 = 1.0101..._2 * 2^-2. Four bits of the
        // significand: 1.010|1... so the tie-free bits round up to 1.011.
        const x = BSRational.fromNumber(1, 3);
        const nearest = convertBSRationalToBSNumber(x, 4,
            BSNumberRoundingMode.FE_TONEAREST);
        expect(nearest.toNumber()).toBe(11 / 32);  // 1.011 * 2^-2
        const toward = convertBSRationalToBSNumber(x, 4,
            BSNumberRoundingMode.FE_TOWARDZERO);
        expect(toward.toNumber()).toBe(10 / 32);   // 1.010 * 2^-2
        const down = convertBSRationalToBSNumber(x, 4,
            BSNumberRoundingMode.FE_DOWNWARD);
        expect(down.toNumber()).toBe(10 / 32);
        const up = convertBSRationalToBSNumber(x, 4,
            BSNumberRoundingMode.FE_UPWARD);
        expect(up.toNumber()).toBe(11 / 32);
    });

    it('rounds -1/3 to 4 bits in each mode', () => {
        const x = BSRational.fromNumber(-1, 3);
        expect(convertBSRationalToBSNumber(x, 4,
            BSNumberRoundingMode.FE_TOWARDZERO).toNumber()).toBe(-10 / 32);
        // Downward on a negative rounds the magnitude up.
        expect(convertBSRationalToBSNumber(x, 4,
            BSNumberRoundingMode.FE_DOWNWARD).toNumber()).toBe(-11 / 32);
        // Upward on a negative truncates the magnitude.
        expect(convertBSRationalToBSNumber(x, 4,
            BSNumberRoundingMode.FE_UPWARD).toNumber()).toBe(-10 / 32);
        expect(convertBSRationalToBSNumber(x, 4,
            BSNumberRoundingMode.FE_TONEAREST).toNumber()).toBe(-11 / 32);
    });

    it('brackets the exact value between the downward and upward results', () => {
        const rand = makeRandom(13579);
        for (let trial = 0; trial < 150; ++trial) {
            const n = Math.trunc(rand() * 2000) - 1000;
            const d = Math.trunc(rand() * 1000) + 1;
            if (n === 0) {
                continue;
            }
            const x = BSRational.fromNumber(n, d);
            const precision = 3 + Math.trunc(rand() * 20);
            const down = BSRational.fromBSNumber(convertBSRationalToBSNumber(
                x, precision, BSNumberRoundingMode.FE_DOWNWARD));
            const up = BSRational.fromBSNumber(convertBSRationalToBSNumber(
                x, precision, BSNumberRoundingMode.FE_UPWARD));
            const near = BSRational.fromBSNumber(convertBSRationalToBSNumber(
                x, precision, BSNumberRoundingMode.FE_TONEAREST));
            const zero = BSRational.fromBSNumber(convertBSRationalToBSNumber(
                x, precision, BSNumberRoundingMode.FE_TOWARDZERO));
            expect(down.lessThanOrEqual(x)).toBe(true);
            expect(x.lessThanOrEqual(up)).toBe(true);
            expect(down.lessThanOrEqual(near)).toBe(true);
            expect(near.lessThanOrEqual(up)).toBe(true);
            // Truncation never increases the magnitude.
            expect(BSRational.fabs(zero).lessThanOrEqual(BSRational.fabs(x)))
                .toBe(true);
            // The nearest result is at least as close as either directed one.
            const eNear = BSRational.fabs(near.sub(x));
            expect(eNear.lessThanOrEqual(BSRational.fabs(down.sub(x)))).toBe(true);
            expect(eNear.lessThanOrEqual(BSRational.fabs(up.sub(x)))).toBe(true);
        }
    });

    it('convertBSRational wraps the BSNumber result', () => {
        const r = convertBSRational(BSRational.fromNumber(1, 3), 4,
            BSNumberRoundingMode.FE_TONEAREST);
        expect(r.toNumber()).toBe(11 / 32);
        expect(r.getDenominator().toNumber()).toBe(1);
    });

    it('convertBSRationalToNumber and ...ToFloat32 use 53 and 24 bits', () => {
        const x = BSRational.fromNumber(1, 3);
        expect(convertBSRationalToNumber(x, BSNumberRoundingMode.FE_TONEAREST))
            .toBe(1 / 3);
        expect(convertBSRationalToFloat32(x, BSNumberRoundingMode.FE_TONEAREST))
            .toBe(Math.fround(1 / 3));
        // 1/3 has an infinite binary expansion, so rounding upward gives a
        // value strictly larger than the nearest double and truncation gives
        // one no larger than it.
        expect(convertBSRationalToNumber(x, BSNumberRoundingMode.FE_UPWARD))
            .toBeGreaterThan(1 / 3);
        expect(convertBSRationalToNumber(x, BSNumberRoundingMode.FE_TOWARDZERO))
            .toBeLessThanOrEqual(1 / 3);
        expect(convertBSRationalToFloat32(x, BSNumberRoundingMode.FE_UPWARD))
            .toBeGreaterThan(1 / 3);
    });

    it('handles an exactly representable rational (loop terminates early)', () => {
        // 5/8 has a finite binary expansion, exercising the early break of
        // the bit loop.
        const x = BSRational.fromNumber(5, 8);
        expect(convertBSRationalToBSNumber(x, 40,
            BSNumberRoundingMode.FE_TONEAREST).toNumber()).toBe(0.625);
        expect(convertBSRationalToBSNumber(x, 3,
            BSNumberRoundingMode.FE_TONEAREST).toNumber()).toBe(0.625);
    });
});

describe('BSRational: math functions', () => {
    it('fabs is exact', () => {
        const x = BSRational.fromNumber(-1, 3);
        const a = BSRational.fabs(x);
        expectExactlyEqual(a, { n: 1n, d: 3n });
        expect(BSRational.fabs(BSRational.fromNumber(1, 3)).equals(a)).toBe(true);
        // fabs does not alias its input.
        a.negate();
        expect(x.toNumber()).toBe(-1 / 3);
    });

    it('frexp splits x into a fraction in [1/2,1) and an exponent', () => {
        const rand = makeRandom(555);
        for (let trial = 0; trial < 100; ++trial) {
            const n = Math.trunc(rand() * 2000) - 1000;
            const d = Math.trunc(rand() * 1000) + 1;
            if (n === 0) {
                continue;
            }
            const x = BSRational.fromNumber(n, d);
            const { result, exponent } = BSRational.frexp(x);
            // x = result * 2^exponent, exactly.
            expect(BSRational.ldexp(result, exponent).equals(x)).toBe(true);
            const mag = BSRational.fabs(result);
            expect(mag.greaterThanOrEqual(BSRational.fromNumber(1, 2))).toBe(true);
            expect(mag.lessThan(BSRational.fromNumber(1))).toBe(true);
            // Matches the double frexp exponent.
            const e = Math.floor(Math.log2(Math.abs(n / d))) + 1;
            expect(exponent).toBe(e);
        }
    });

    it('frexp of zero is zero', () => {
        const { result, exponent } = BSRational.frexp(new BSRational());
        expect(result.getSign()).toBe(0);
        expect(exponent).toBe(0);
    });

    it('frexp does not modify its input', () => {
        const x = BSRational.fromNumber(3, 7);
        BSRational.frexp(x);
        expect(x.toNumber()).toBe(3 / 7);
    });

    it('ldexp scales by a power of two exactly', () => {
        const x = BSRational.fromNumber(1, 3);
        expectExactlyEqual(BSRational.ldexp(x, 5), { n: 32n, d: 3n });
        expectExactlyEqual(BSRational.ldexp(x, -5), { n: 1n, d: 96n });
        // The input is unchanged.
        expect(x.toNumber()).toBe(1 / 3);
    });

    it('routes the transcendental functions through double', () => {
        const half = BSRational.fromNumber(1, 2);
        expect(BSRational.sqrt(BSRational.fromNumber(2)).toNumber())
            .toBe(Math.sqrt(2));
        expect(BSRational.sin(half).toNumber()).toBe(Math.sin(0.5));
        expect(BSRational.cos(half).toNumber()).toBe(Math.cos(0.5));
        expect(BSRational.tan(half).toNumber()).toBe(Math.tan(0.5));
        expect(BSRational.asin(half).toNumber()).toBe(Math.asin(0.5));
        expect(BSRational.acos(half).toNumber()).toBe(Math.acos(0.5));
        expect(BSRational.atan(half).toNumber()).toBe(Math.atan(0.5));
        expect(BSRational.atan2(half, BSRational.fromNumber(1)).toNumber())
            .toBe(Math.atan2(0.5, 1));
        expect(BSRational.sinh(half).toNumber()).toBe(Math.sinh(0.5));
        expect(BSRational.cosh(half).toNumber()).toBe(Math.cosh(0.5));
        expect(BSRational.tanh(half).toNumber()).toBe(Math.tanh(0.5));
        expect(BSRational.asinh(half).toNumber()).toBe(Math.asinh(0.5));
        expect(BSRational.acosh(BSRational.fromNumber(2)).toNumber())
            .toBe(Math.acosh(2));
        expect(BSRational.atanh(half).toNumber()).toBe(Math.atanh(0.5));
        expect(BSRational.exp(half).toNumber()).toBe(Math.exp(0.5));
        expect(BSRational.exp2(half).toNumber()).toBe(Math.pow(2, 0.5));
        expect(BSRational.log(half).toNumber()).toBe(Math.log(0.5));
        expect(BSRational.log2(half).toNumber()).toBe(Math.log2(0.5));
        expect(BSRational.log10(half).toNumber()).toBe(Math.log10(0.5));
        expect(BSRational.pow(BSRational.fromNumber(3), BSRational.fromNumber(4))
            .toNumber()).toBe(81);
    });

    it('floor, ceil and fmod', () => {
        const x = BSRational.fromNumber(7, 2);
        expect(BSRational.floor(x).toNumber()).toBe(3);
        expect(BSRational.ceil(x).toNumber()).toBe(4);
        expect(BSRational.fmod(x, BSRational.fromNumber(1)).toNumber()).toBe(0.5);
        const y = BSRational.fromNumber(-7, 2);
        expect(BSRational.floor(y).toNumber()).toBe(-4);
        expect(BSRational.ceil(y).toNumber()).toBe(-3);
    });

    it('remainder rounds the quotient to the nearest even integer', () => {
        // 3.5 / 1 = 3.5, ties to even gives n = 4, so the remainder is -0.5.
        expect(BSRational.remainder(BSRational.fromNumber(7, 2),
            BSRational.fromNumber(1)).toNumber()).toBe(-0.5);
        // 2.5 / 1 = 2.5, ties to even gives n = 2, so the remainder is 0.5.
        expect(BSRational.remainder(BSRational.fromNumber(5, 2),
            BSRational.fromNumber(1)).toNumber()).toBe(0.5);
        expect(BSRational.remainder(BSRational.fromNumber(-5, 2),
            BSRational.fromNumber(1)).toNumber()).toBe(-0.5);
        expect(BSRational.remainder(BSRational.fromNumber(-7, 2),
            BSRational.fromNumber(1)).toNumber()).toBe(0.5);
        expect(BSRational.remainder(BSRational.fromNumber(9), BSRational.fromNumber(4))
            .toNumber()).toBe(1);
    });

    it('supplies the gte helper functions', () => {
        expect(BSRational.isign(BSRational.fromNumber(-3))).toBe(-1);
        expect(BSRational.isign(new BSRational())).toBe(0);
        expect(BSRational.isign(BSRational.fromNumber(3))).toBe(1);
        expect(BSRational.sign(BSRational.fromNumber(-3)).toNumber()).toBe(-1);
        expect(BSRational.sqr(BSRational.fromNumber(3)).toNumber()).toBe(9);
        expect(BSRational.invsqrt(BSRational.fromNumber(4)).toNumber()).toBe(0.5);
        expect(BSRational.saturate(BSRational.fromNumber(3, 2)).toNumber()).toBe(1);
        expect(BSRational.saturate(BSRational.fromNumber(-1)).toNumber()).toBe(0);
        expect(BSRational.clamp(BSRational.fromNumber(5), BSRational.fromNumber(0),
            BSRational.fromNumber(2)).toNumber()).toBe(2);
        expect(BSRational.cospi(BSRational.fromNumber(1)).toNumber()).toBe(-1);
        expect(BSRational.sinpi(new BSRational()).toNumber()).toBe(0);
        expect(BSRational.exp10(BSRational.fromNumber(2)).toNumber())
            .toBeCloseTo(100, 10);
        expect(BSRational.atandivpi(BSRational.fromNumber(1)).toNumber())
            .toBeCloseTo(0.25, 12);
        expect(BSRational.atan2divpi(BSRational.fromNumber(1),
            BSRational.fromNumber(1)).toNumber()).toBeCloseTo(0.25, 12);
    });
});

describe('BSRational: type traits', () => {
    it('is an arbitrary-precision type with division', () => {
        const x = new BSRational();
        expect(isArbitraryPrecision(x)).toBe(true);
        expect(hasDivisionOperator(x)).toBe(true);
        // BSNumber has no division operator, BSRational does.
        expect(hasDivisionOperator(new BSNumber())).toBe(false);
        expect(hasDivisionOperator(1.5)).toBe(true);
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md): property-based cross-checks
// of every BSRational.h function against exact bigint rational arithmetic.
// ---------------------------------------------------------------------------

describe('BSRational verification', () => {
    // Rationals built from small integers, which exercise the shared-factor
    // and zero-numerator paths often.
    const smallRational: fc.Arbitrary<BSRational> =
        fc.tuple(fc.integer({ min: -64, max: 64 }),
            fc.integer({ min: 1, max: 64 }))
            .map(([n, d]) => BSRational.fromNumber(n, d));

    // Rationals far wider than a double, which exercise the bigint paths.
    const wideRational: fc.Arbitrary<BSRational> =
        fc.tuple(fc.bigInt({ min: -(2n ** 90n), max: 2n ** 90n }),
            fc.bigInt({ min: 1n, max: 2n ** 90n }))
            .map(([n, d]) => BSRational.fromBigInt(n, d));

    const anyRational: fc.Arbitrary<BSRational> =
        fc.oneof(smallRational, wideRational);

    const nonzeroRational: fc.Arbitrary<BSRational> =
        anyRational.filter(x => x.getSign() !== 0);

    function ratAbs(a: Rat): Rat {
        return ratMake(a.n < 0n ? -a.n : a.n, a.d);
    }

    it('keeps the class invariant and the exact value through every operation',
        () => {
            check(fc.tuple(anyRational, nonzeroRational), ([a, b]) => {
                const ea = exact(a);
                const eb = exact(b);
                const cases: [BSRational, Rat][] = [
                    [a.add(b), ratAdd(ea, eb)],
                    [a.sub(b), ratSub(ea, eb)],
                    [a.mul(b), ratMul(ea, eb)],
                    [a.div(b), ratDiv(ea, eb)],
                    [a.negated(), ratMake(-ea.n, ea.d)],
                    [BSRational.fabs(a), ratAbs(ea)]
                ];
                for (const [result, expected] of cases) {
                    expectValid(result);
                    expectExactlyEqual(result, expected);
                }
            });
        });

    it('orders exactly as the bigint comparison does', () => {
        check(fc.tuple(anyRational, anyRational), ([a, b]) => {
            const c = ratCompare(exact(a), exact(b));
            expect(a.lessThan(b)).toBe(c < 0);
            expect(a.greaterThan(b)).toBe(c > 0);
            expect(a.equals(b)).toBe(c === 0);
            expect(a.notEquals(b)).toBe(c !== 0);
            expect(a.lessThanOrEqual(b)).toBe(c <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(c >= 0);
            // getSign is the comparison against zero.
            expect(Math.sign(a.getSign())).toBe(
                ratCompare(exact(a), { n: 0n, d: 1n }));
        });
    });

    it('satisfies the field axioms exactly', () => {
        check(fc.tuple(smallRational, smallRational, smallRational),
            ([a, b, c]) => {
                // Associativity and commutativity of + and *.
                expect(a.add(b).add(c).equals(a.add(b.add(c)))).toBe(true);
                expect(a.mul(b).mul(c).equals(a.mul(b.mul(c)))).toBe(true);
                expect(a.add(b).equals(b.add(a))).toBe(true);
                expect(a.mul(b).equals(b.mul(a))).toBe(true);
                // Distributivity.
                expect(a.mul(b.add(c)).equals(a.mul(b).add(a.mul(c)))).toBe(true);
                // Additive and multiplicative inverses.
                expect(a.sub(a).getSign()).toBe(0);
                if (a.getSign() !== 0) {
                    expect(a.div(a).equals(BSRational.fromNumber(1))).toBe(true);
                }
                // fma / robustSOP / robustDOP are the exact expressions.
                expect(BSRational.fma(a, b, c).equals(a.mul(b).add(c))).toBe(true);
                expect(BSRational.robustSOP(a, b, a, c)
                    .equals(a.mul(b).add(a.mul(c)))).toBe(true);
                expect(BSRational.robustDOP(a, b, a, c)
                    .equals(a.mul(b).sub(a.mul(c)))).toBe(true);
            });
    });

    it('frexp and ldexp are exact and invert each other', () => {
        check(anyRational, a => {
            const { result, exponent } = BSRational.frexp(a);
            expectValid(result);
            // x = result * 2^exponent.
            const back = BSRational.ldexp(result, exponent);
            expect(back.equals(a)).toBe(true);
            if (a.getSign() === 0) {
                expect(exponent).toBe(0);
                return;
            }
            // |result| is in [1/2, 1).
            const magnitude = BSRational.fabs(result);
            expect(magnitude.greaterThanOrEqual(
                BSRational.fromNumber(1, 2))).toBe(true);
            expect(magnitude.lessThan(BSRational.fromNumber(1))).toBe(true);
        });
    });

    it('brackets the exact value with the directed rounding modes', () => {
        check(fc.tuple(anyRational, fc.integer({ min: 1, max: 24 })),
            ([a, precision]) => {
                const down = convertBSRationalToBSNumber(a, precision,
                    BSNumberRoundingMode.FE_DOWNWARD);
                const up = convertBSRationalToBSNumber(a, precision,
                    BSNumberRoundingMode.FE_UPWARD);
                const near = convertBSRationalToBSNumber(a, precision,
                    BSNumberRoundingMode.FE_TONEAREST);
                const zero = convertBSRationalToBSNumber(a, precision,
                    BSNumberRoundingMode.FE_TOWARDZERO);
                for (const x of [down, up, near, zero]) {
                    expectValidBSNumber(x);
                    expect(x.getNumBits()).toBeLessThanOrEqual(precision);
                }

                const ea = exact(a);
                const rd = exact(BSRational.fromBSNumber(down));
                const ru = exact(BSRational.fromBSNumber(up));
                const rn = exact(BSRational.fromBSNumber(near));
                const rz = exact(BSRational.fromBSNumber(zero));
                expect(ratCompare(rd, ea)).toBeLessThanOrEqual(0);
                expect(ratCompare(ru, ea)).toBeGreaterThanOrEqual(0);
                // Round to nearest lands on one of the two directed results
                // and is no farther from the exact value than the other.
                expect(ratCompare(rn, rd) === 0 || ratCompare(rn, ru) === 0)
                    .toBe(true);
                const other = ratCompare(rn, rd) === 0 ? ru : rd;
                expect(ratCompare(ratAbs(ratSub(rn, ea)),
                    ratAbs(ratSub(other, ea)))).toBeLessThanOrEqual(0);
                // Round toward zero truncates the magnitude.
                expect(ratCompare(rz, a.getSign() < 0 ? ru : rd)).toBe(0);
                // An exactly representable value is fixed by every mode.
                if (ratCompare(rd, ru) === 0) {
                    expect(ratCompare(rd, ea)).toBe(0);
                }
            });
    });

    it('converts a double-valued rational back to that double exactly', () => {
        check(fc.tuple(fc.double({ min: -1e12, max: 1e12, noNaN: true }),
            fc.integer({ min: -40, max: 40 })), ([value, shift]) => {
                const x = BSRational.fromNumber(value);
                // A signed zero converts to the unsigned zero, as it does in
                // C++ (BSNumber has a single zero with sign 0).
                expect(x.toNumber()).toBe(value === 0 ? 0 : value);
                // Scaling by a power of two is exact when it neither
                // overflows nor underflows the double range.
                const expected = value * Math.pow(2, shift);
                if (Number.isFinite(expected) && expected !== 0
                    && Math.abs(expected) > 1e-280) {
                    expect(BSRational.ldexp(x, shift).toNumber()).toBe(expected);
                }
            });
    });

    it('fromString agrees with the exact decimal fraction', () => {
        check(fc.tuple(fc.integer({ min: -999, max: 999 }),
            fc.stringMatching(/^[0-9]{0,6}$/)), ([intPart, frcPart]) => {
                const sign = intPart < 0 ? '-' : '';
                const text = sign + String(Math.abs(intPart)) + '.' + frcPart;
                const x = BSRational.fromString(text);
                expectValid(x);
                const scale = 10n ** BigInt(frcPart.length);
                const magnitude = BigInt(Math.abs(intPart)) * scale
                    + (frcPart.length > 0 ? BigInt(frcPart) : 0n);
                expectExactlyEqual(x, ratMake(
                    intPart < 0 ? -magnitude : magnitude, scale));
            });
    });

    // Regression for a port defect: BSRational.remainder formed the quotient
    // as a double, so above 2^53 the nearest integer n was off by one and
    // 'dx - n*dy' rounded on top of that. Upstream converts to double and
    // calls std::remainder, whose result is exact.
    it('remainder is the exact IEEE remainder of the two doubles', () => {
        // std::remainder(1e17, 3) == 1; the double computation returned 0.
        expect(BSRational.remainder(BSRational.fromNumber(1e17),
            BSRational.fromNumber(3)).toNumber()).toBe(1);

        const exactRemainder = (dx: number, dy: number): Rat => {
            const ex = exact(BSRational.fromNumber(dx));
            const ey = exact(BSRational.fromNumber(dy));
            // n = the integer nearest to x/y, ties to even, computed on
            // exact bigint fractions.
            const num = ex.n * ey.d;
            const den = ex.d * ey.n;
            const sign = den < 0n ? -1n : 1n;
            const a = num * sign;
            const b = den * sign;
            let q = a / b;
            let r = a - q * b;
            if (r < 0n) { q -= 1n; r += b; }   // floor division
            if (2n * r > b || (2n * r === b && (q & 1n) === 1n)) { q += 1n; }
            return ratSub(ex, ratMul({ n: q, d: 1n }, ey));
        };

        check(fc.tuple(fc.double({ min: -1e18, max: 1e18, noNaN: true }),
            fc.double({ min: -1e6, max: 1e6, noNaN: true })
                .filter(y => Math.abs(y) > 1e-3)),
            ([dx, dy]) => {
                const got = BSRational.remainder(BSRational.fromNumber(dx),
                    BSRational.fromNumber(dy));
                expectExactlyEqual(got, exactRemainder(dx, dy));
            });
    });

    it('setSign, negate and the zero canonicalization agree with the value',
        () => {
            check(fc.tuple(anyRational, fc.constantFrom(-1, 1)),
                ([a, s]) => {
                    const b = a.clone();
                    b.setSign(s);
                    expect(b.getDenominator().getSign()).toBe(1);
                    if (a.getSign() !== 0) {
                        expectValid(b);
                        expect(b.getSign()).toBe(s);
                        const magnitude = ratAbs(exact(a));
                        expectExactlyEqual(b, s < 0
                            ? ratMake(-magnitude.n, magnitude.d) : magnitude);
                    }

                    const c = a.clone();
                    c.negate();
                    expectValid(c);
                    // Negating zero stays a canonical zero (int32 -0 is 0).
                    expect(Object.is(c.getSign(), -0)).toBe(false);
                    expectExactlyEqual(c, ratMake(-exact(a).n, exact(a).d));
                });
        });

    it('clone keeps C++ value semantics', () => {
        check(anyRational, a => {
            const before = exact(a);
            const b = a.clone();
            b.negate();
            b.getNumerator().setBiasedExponent(
                b.getNumerator().getBiasedExponent() + 3);
            // Mutating the clone must not disturb the original.
            expect(ratCompare(exact(a), before)).toBe(0);
            expectValid(a);
        });
    });
});
