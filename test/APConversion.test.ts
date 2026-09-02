import { describe, it, expect } from 'vitest';
import { APConversion } from '../src/APConversion';
import type { APConversionQFN1 } from '../src/APConversion';
import { BSRational } from '../src/BSRational';

function r(numerator: number, denominator: number = 1): BSRational {
    return BSRational.fromNumber(numerator, denominator);
}

const zero = r(0);

// The sign of t - sqrt(aSqr) for aSqr >= 0, computed exactly.
function compareToSqrt(t: BSRational, aSqr: BSRational): number {
    if (t.getSign() < 0) {
        return -1;
    }
    const difference = t.mul(t).sub(aSqr);
    return difference.getSign();
}

// The sign of t - (a + b), where a = sqrt(aSqr) and b = sqrt(bSqr), computed
// exactly. With L = t^2 - (a^2 + b^2), the quantity t^2 - (a+b)^2 = L - 2*a*b
// is negative when L <= 0 and otherwise has the sign of L^2 - 4*a^2*b^2.
function compareToApB(t: BSRational, aSqr: BSRational, bSqr: BSRational): number {
    if (t.getSign() < 0) {
        return -1;
    }
    const el = t.mul(t).sub(aSqr.add(bSqr));
    if (el.getSign() <= 0) {
        return -1;
    }
    return el.mul(el).sub(BSRational.ldexp(aSqr.mul(bSqr), 2)).getSign();
}

// The sign of t - (a - b), where a = sqrt(aSqr) >= b = sqrt(bSqr) >= 0,
// computed exactly. Here t^2 - (a-b)^2 = L + 2*a*b is positive when L >= 0
// and otherwise has the sign of 4*a^2*b^2 - L^2.
function compareToAmB(t: BSRational, aSqr: BSRational, bSqr: BSRational): number {
    if (t.getSign() < 0) {
        return -1;
    }
    const el = t.mul(t).sub(aSqr.add(bSqr));
    if (el.getSign() >= 0) {
        return +1;
    }
    return BSRational.ldexp(aSqr.mul(bSqr), 2).sub(el.mul(el)).getSign();
}

const PRECISION = 60;
const MAX_ITERATIONS = 24;

function makeConverter(): APConversion {
    return new APConversion(PRECISION, MAX_ITERATIONS);
}

// The convergence threshold 2^-precision that the estimators use.
function threshold(precision: number = PRECISION): BSRational {
    return BSRational.ldexp(r(1), -precision);
}

describe('APConversion construction', () => {
    it('validates the constructor arguments and exposes them', () => {
        expect(() => new APConversion(0, 8)).toThrow(/Invalid precision/);
        expect(() => new APConversion(-1, 8)).toThrow(/Invalid precision/);
        expect(() => new APConversion(32, 0)).toThrow(/Invalid maximum iterations/);

        const converter = new APConversion(32, 8);
        expect(converter.getPrecision()).toBe(32);
        expect(converter.getMaxIterations()).toBe(8);

        converter.setPrecision(16);
        converter.setMaxIterations(4);
        expect(converter.getPrecision()).toBe(16);
        expect(converter.getMaxIterations()).toBe(4);
        expect(() => converter.setPrecision(0)).toThrow(/Invalid precision/);
        expect(() => converter.setMaxIterations(0)).toThrow(/Invalid maximum iterations/);

        // The new precision is honored: the bounding interval for sqrt(2) is
        // wider at 16 bits than at 60 bits.
        const coarse = converter.estimateSqrt(r(2));
        const fine = makeConverter().estimateSqrt(r(2));
        expect(coarse.aMax.sub(coarse.aMin).lessThan(threshold(16))).toBe(true);
        expect(fine.aMax.sub(fine.aMin).lessThan(threshold())).toBe(true);
    });
});

describe('APConversion.estimateSqrt', () => {
    it('brackets sqrt(a^2) to the requested precision', () => {
        const converter = makeConverter();
        // Both even and odd frexp exponents of a^2 are exercised: 2 = 0.5*2^2
        // (even), 5 = 0.625*2^3 (odd), 1/3 = 0.666*2^-1 (odd),
        // 12345/64 (even), 4 and 9/16 are perfect squares.
        const inputs = [r(2), r(5), r(1, 3), r(12345, 64), r(4), r(9, 16),
            r(1), r(1, 1024), r(1000000)];
        for (const aSqr of inputs) {
            const { numIterates, aMin, aMax } = converter.estimateSqrt(aSqr);
            expect(numIterates).toBeGreaterThanOrEqual(1);
            expect(numIterates).toBeLessThanOrEqual(MAX_ITERATIONS);
            // Rigorous bracketing: aMin <= sqrt(aSqr) <= aMax.
            expect(compareToSqrt(aMin, aSqr)).toBeLessThanOrEqual(0);
            expect(compareToSqrt(aMax, aSqr)).toBeGreaterThanOrEqual(0);
            // The endpoints agree to the requested precision.
            expect(aMax.sub(aMin).lessThan(threshold())).toBe(true);
            // Cross-check against the FPU.
            expect(aMin.toNumber()).toBeCloseTo(Math.sqrt(aSqr.toNumber()), 12);
        }
    });

    it('is exact for perfect squares of dyadic rationals', () => {
        const converter = makeConverter();
        const { aMin, aMax } = converter.estimateSqrt(r(9, 16));
        expect(aMin.lessThanOrEqual(r(3, 4))).toBe(true);
        expect(r(3, 4).lessThanOrEqual(aMax)).toBe(true);
    });

    it('handles the zero input', () => {
        const converter = makeConverter();
        const { aMin, aMax } = converter.estimateSqrt(zero);
        expect(aMin.getSign()).toBe(0);
        expect(aMax.getSign()).toBeGreaterThanOrEqual(0);
        expect(aMax.lessThan(threshold())).toBe(true);
    });

    it('averages the bounds for the single-value estimate', () => {
        const converter = makeConverter();
        for (const aSqr of [r(2), r(3), r(7, 5)]) {
            const bounds = converter.estimateSqrt(aSqr);
            const { numIterates, a } = converter.estimateSqrtValue(aSqr);
            expect(numIterates).toBe(bounds.numIterates);
            expect(bounds.aMin.lessThanOrEqual(a)).toBe(true);
            expect(a.lessThanOrEqual(bounds.aMax)).toBe(true);
            expect(a.toNumber()).toBeCloseTo(Math.sqrt(aSqr.toNumber()), 14);
        }
    });
});

describe('APConversion.estimateApB', () => {
    it('brackets sqrt(a^2) + sqrt(b^2) to the requested precision', () => {
        const converter = makeConverter();
        const inputs: Array<[BSRational, BSRational]> = [
            [r(2), r(3)], [r(4), r(9)], [r(1, 3), r(1, 7)],
            [r(100), r(1)], [r(2), r(2)], [r(12345), r(1, 64)]
        ];
        for (const [aSqr, bSqr] of inputs) {
            const { numIterates, tMin, tMax } = converter.estimateApB(aSqr, bSqr);
            expect(numIterates).toBeGreaterThanOrEqual(1);
            expect(numIterates).toBeLessThanOrEqual(MAX_ITERATIONS);
            expect(compareToApB(tMin, aSqr, bSqr)).toBeLessThanOrEqual(0);
            expect(compareToApB(tMax, aSqr, bSqr)).toBeGreaterThanOrEqual(0);
            expect(tMax.sub(tMin).lessThan(threshold())).toBe(true);
            expect(tMin.toNumber()).toBeCloseTo(
                Math.sqrt(aSqr.toNumber()) + Math.sqrt(bSqr.toNumber()), 12);
        }
    });

    it('brackets the exact value 2 + 3 for a^2 = 4 and b^2 = 9', () => {
        const { tMin, tMax } = makeConverter().estimateApB(r(4), r(9));
        expect(tMin.lessThanOrEqual(r(5))).toBe(true);
        expect(r(5).lessThanOrEqual(tMax)).toBe(true);
    });
});

describe('APConversion.estimateAmB', () => {
    it('brackets sqrt(a^2) - sqrt(b^2) to the requested precision', () => {
        const converter = makeConverter();
        // The first three inputs have a positive second derivative at a-b
        // (a^2 is much larger than b^2); the rest have a negative one. The
        // last pair is chosen near the boundary a^2/b^2 = (7+sqrt(45))/2,
        // where the initial guess can fall outside the expected basin and
        // the bisection code runs.
        const inputs: Array<[BSRational, BSRational]> = [
            [r(100), r(1)], [r(7), r(1)], [r(1000), r(3)],
            [r(3), r(2)], [r(2), r(2)], [r(1, 3), r(1, 7)],
            [BSRational.fromNumber((7 + Math.sqrt(45)) / 2), r(1)]
        ];
        for (const [aSqr, bSqr] of inputs) {
            const { numIterates, tMin, tMax } = converter.estimateAmB(aSqr, bSqr);
            expect(numIterates).toBeGreaterThanOrEqual(1);
            expect(numIterates).toBeLessThanOrEqual(MAX_ITERATIONS);
            expect(compareToAmB(tMin, aSqr, bSqr)).toBeLessThanOrEqual(0);
            expect(compareToAmB(tMax, aSqr, bSqr)).toBeGreaterThanOrEqual(0);
            expect(tMax.sub(tMin).lessThan(threshold())).toBe(true);
            expect(tMax.toNumber()).toBeCloseTo(
                Math.sqrt(aSqr.toNumber()) - Math.sqrt(bSqr.toNumber()), 12);
        }
    });

    it('brackets the exact value 10 - 1 for a^2 = 100 and b^2 = 1', () => {
        const { tMin, tMax } = makeConverter().estimateAmB(r(100), r(1));
        expect(tMin.lessThanOrEqual(r(9))).toBe(true);
        expect(r(9).lessThanOrEqual(tMax)).toBe(true);
    });

    it('brackets zero when a^2 equals b^2', () => {
        const { tMin, tMax } = makeConverter().estimateAmB(r(5), r(5));
        expect(tMin.lessThanOrEqual(zero)).toBe(true);
        expect(zero.lessThanOrEqual(tMax)).toBe(true);
        expect(tMax.sub(tMin).lessThan(threshold())).toBe(true);
    });
});

describe('APConversion.estimate for quadratic field numbers', () => {
    function qfn(x0: BSRational, x1: BSRational, d: BSRational): APConversionQFN1 {
        return { x: [x0, x1], d };
    }

    it('brackets x + y*sqrt(d) for positive and negative y', () => {
        const converter = makeConverter();
        const cases: Array<[BSRational, BSRational, BSRational]> = [
            [r(1, 2), r(3), r(2)],        // 0.5 + 3*sqrt(2)
            [r(1, 2), r(-3), r(2)],       // 0.5 - 3*sqrt(2)
            [r(-7, 3), r(1, 5), r(11)],   // -7/3 + sqrt(11)/5
            [r(0), r(1), r(5)]            // sqrt(5)
        ];
        for (const [x, y, d] of cases) {
            const { numIterates, qMin, qMax } = converter.estimate(qfn(x, y, d));
            expect(numIterates).toBeGreaterThanOrEqual(1);
            expect(numIterates).toBeLessThanOrEqual(MAX_ITERATIONS);

            // The exact test: y*sqrt(d) = +-sqrt(y^2*d), so compare the
            // shifted endpoints with the square root of y^2*d.
            const aSqr = y.mul(y).mul(d);
            if (y.getSign() > 0) {
                expect(compareToSqrt(qMin.sub(x), aSqr)).toBeLessThanOrEqual(0);
                expect(compareToSqrt(qMax.sub(x), aSqr)).toBeGreaterThanOrEqual(0);
            } else {
                expect(compareToSqrt(x.sub(qMax), aSqr)).toBeLessThanOrEqual(0);
                expect(compareToSqrt(x.sub(qMin), aSqr)).toBeGreaterThanOrEqual(0);
            }
            expect(qMax.sub(qMin).lessThan(threshold())).toBe(true);

            const expected = x.toNumber()
                + y.toNumber() * Math.sqrt(d.toNumber());
            expect(qMin.toNumber()).toBeCloseTo(expected, 12);
            expect(qMax.toNumber()).toBeCloseTo(expected, 12);
        }
    });

    it('is degenerate when d or y is zero', () => {
        const converter = makeConverter();
        for (const q of [qfn(r(7, 2), r(3), zero), qfn(r(7, 2), zero, r(2))]) {
            const { numIterates, qMin, qMax } = converter.estimate(q);
            expect(numIterates).toBe(0);
            expect(qMin.toNumber()).toBe(3.5);
            expect(qMax.toNumber()).toBe(3.5);
            const value = converter.estimateValue(q);
            expect(value.numIterates).toBe(0);
            expect(value.qEstimate.toNumber()).toBe(3.5);
        }
    });

    it('averages the bounds for the single-value estimate', () => {
        const converter = makeConverter();
        const q = qfn(r(1, 2), r(3), r(2));
        const bounds = converter.estimate(q);
        const { numIterates, qEstimate } = converter.estimateValue(q);
        expect(numIterates).toBe(bounds.numIterates);
        expect(bounds.qMin.lessThanOrEqual(qEstimate)).toBe(true);
        expect(qEstimate.lessThanOrEqual(bounds.qMax)).toBe(true);
        expect(qEstimate.toNumber()).toBeCloseTo(0.5 + 3 * Math.SQRT2, 14);
    });
});

describe('APConversion randomized cross-checks', () => {
    it('brackets the exact roots for random rational inputs', () => {
        const converter = new APConversion(40, 32);
        const width = threshold(40);
        // A deterministic linear congruential generator.
        let seed = 24680135;
        const next = (bound: number): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed % bound;
        };

        for (let trial = 0; trial < 12; ++trial) {
            let aSqr = r(next(50) + 1, next(20) + 1);
            let bSqr = r(next(50) + 1, next(20) + 1);
            if (aSqr.lessThan(bSqr)) {
                const swap = aSqr; aSqr = bSqr; bSqr = swap;
            }

            const root = converter.estimateSqrt(aSqr);
            expect(compareToSqrt(root.aMin, aSqr)).toBeLessThanOrEqual(0);
            expect(compareToSqrt(root.aMax, aSqr)).toBeGreaterThanOrEqual(0);
            expect(root.aMax.sub(root.aMin).lessThan(width)).toBe(true);
            expect(root.aMin.toNumber()).toBeCloseTo(Math.sqrt(aSqr.toNumber()), 10);

            const sum = converter.estimateApB(aSqr, bSqr);
            expect(compareToApB(sum.tMin, aSqr, bSqr)).toBeLessThanOrEqual(0);
            expect(compareToApB(sum.tMax, aSqr, bSqr)).toBeGreaterThanOrEqual(0);
            expect(sum.tMax.sub(sum.tMin).lessThan(width)).toBe(true);

            const difference = converter.estimateAmB(aSqr, bSqr);
            expect(compareToAmB(difference.tMin, aSqr, bSqr)).toBeLessThanOrEqual(0);
            expect(compareToAmB(difference.tMax, aSqr, bSqr)).toBeGreaterThanOrEqual(0);
            expect(difference.tMax.sub(difference.tMin).lessThan(width)).toBe(true);

            // The quadratic field number x + y*sqrt(d) with y = 1, d = aSqr.
            const x = r(next(11) - 5, next(4) + 1);
            const q = converter.estimate({ x: [x, r(1)], d: aSqr });
            expect(compareToSqrt(q.qMin.sub(x), aSqr)).toBeLessThanOrEqual(0);
            expect(compareToSqrt(q.qMax.sub(x), aSqr)).toBeGreaterThanOrEqual(0);
        }
    });
});
