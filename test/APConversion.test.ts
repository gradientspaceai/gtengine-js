import { describe, it, expect } from 'vitest';
import { APConversion } from '../src/APConversion.js';
import type { APConversionQFN1 } from '../src/APConversion.js';
import { BSRational } from '../src/BSRational.js';
import { check, fc } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). Every bound is checked with
// exact rational predicates rather than floating-point comparisons: for
// t >= 0 the sign of f(t) = (t^2 - (a^2+b^2))^2 - 4*a^2*b^2 together with the
// sign of t^2 - (a^2+b^2) decides t against a+b and against |a-b|, because
// the roots of f on the nonnegative axis are exactly |a-b| and a+b.
// ---------------------------------------------------------------------------

describe('APConversion verification', () => {
    const PRECISION = 24;
    const MAX_ITERATIONS = 64;

    function threshold(precision: number): BSRational {
        return BSRational.ldexp(BSRational.fromNumber(1), -precision);
    }

    // Positive rationals with small numerators and denominators.
    const positiveSqr: fc.Arbitrary<BSRational> =
        fc.tuple(fc.integer({ min: 1, max: 400 }),
            fc.integer({ min: 1, max: 32 }))
            .map(([n, d]) => BSRational.fromNumber(n, d));

    const zero = new BSRational();

    function sqr(t: BSRational): BSRational {
        return t.mul(t);
    }

    // t <= sqrt(A) for A >= 0.
    function leSqrt(t: BSRational, A: BSRational): boolean {
        return t.getSign() <= 0 || sqr(t).lessThanOrEqual(A);
    }

    // t >= sqrt(A) for A >= 0.
    function geSqrt(t: BSRational, A: BSRational): boolean {
        return t.getSign() >= 0 && sqr(t).greaterThanOrEqual(A);
    }

    function quartic(t: BSRational, aSqr: BSRational, bSqr: BSRational): BSRational {
        const s = aSqr.add(bSqr);
        const inner = sqr(t).sub(s);
        return sqr(inner).sub(BSRational.fromNumber(4).mul(aSqr).mul(bSqr));
    }

    // t <= sqrt(aSqr) + sqrt(bSqr).
    function leSum(t: BSRational, aSqr: BSRational, bSqr: BSRational): boolean {
        if (t.getSign() <= 0) { return true; }
        const beyond = quartic(t, aSqr, bSqr).getSign() > 0
            && sqr(t).greaterThan(aSqr.add(bSqr));
        return !beyond;
    }

    // t >= sqrt(aSqr) + sqrt(bSqr).
    function geSum(t: BSRational, aSqr: BSRational, bSqr: BSRational): boolean {
        return t.getSign() >= 0
            && quartic(t, aSqr, bSqr).getSign() >= 0
            && sqr(t).greaterThanOrEqual(aSqr.add(bSqr));
    }

    // t <= |sqrt(aSqr) - sqrt(bSqr)|.
    function leDiff(t: BSRational, aSqr: BSRational, bSqr: BSRational): boolean {
        if (t.getSign() <= 0) { return true; }
        return quartic(t, aSqr, bSqr).getSign() >= 0
            && sqr(t).lessThanOrEqual(aSqr.add(bSqr));
    }

    // t >= |sqrt(aSqr) - sqrt(bSqr)|.
    function geDiff(t: BSRational, aSqr: BSRational, bSqr: BSRational): boolean {
        if (t.getSign() < 0) { return false; }
        const below = quartic(t, aSqr, bSqr).getSign() > 0
            && sqr(t).lessThan(aSqr.add(bSqr));
        return !below;
    }

    it('brackets sqrt(aSqr) to the requested precision', () => {
        const conv = new APConversion(PRECISION, MAX_ITERATIONS);
        const eps = threshold(PRECISION);
        check(positiveSqr, aSqr => {
            const { numIterates, aMin, aMax } = conv.estimateSqrt(aSqr);
            expect(numIterates).toBeLessThanOrEqual(MAX_ITERATIONS);
            expect(aMin.getSign()).toBeGreaterThan(0);
            // aMin <= sqrt(aSqr) <= aMax, exactly.
            expect(leSqrt(aMin, aSqr)).toBe(true);
            expect(geSqrt(aMax, aSqr)).toBe(true);
            expect(aMax.sub(aMin).lessThan(eps)).toBe(true);

            // The single-value overload is the midpoint of the bracket.
            const single = conv.estimateSqrtValue(aSqr);
            expect(single.numIterates).toBe(numIterates);
            expect(single.a.equals(
                BSRational.ldexp(aMin.add(aMax), -1))).toBe(true);
            expect(aMin.lessThanOrEqual(single.a)).toBe(true);
            expect(single.a.lessThanOrEqual(aMax)).toBe(true);
        }, 40);
    });

    it('brackets sqrt(aSqr) + sqrt(bSqr) to the requested precision', () => {
        const conv = new APConversion(PRECISION, MAX_ITERATIONS);
        const eps = threshold(PRECISION);
        check(fc.tuple(positiveSqr, positiveSqr), ([aSqr, bSqr]) => {
            const { numIterates, tMin, tMax } = conv.estimateApB(aSqr, bSqr);
            expect(numIterates).toBeLessThanOrEqual(MAX_ITERATIONS);
            expect(leSum(tMin, aSqr, bSqr)).toBe(true);
            expect(geSum(tMax, aSqr, bSqr)).toBe(true);
            expect(tMax.sub(tMin).lessThan(eps)).toBe(true);
        }, 30);
    });

    it('brackets sqrt(aSqr) - sqrt(bSqr) to the requested precision', () => {
        const conv = new APConversion(PRECISION, MAX_ITERATIONS);
        const eps = threshold(PRECISION);
        // Upstream requires aSqr >= bSqr so that a - b >= 0.
        const ordered = fc.tuple(positiveSqr, positiveSqr)
            .map(([p, q]) => (p.lessThan(q) ? [q, p] : [p, q]) as
                [BSRational, BSRational]);
        check(ordered, ([aSqr, bSqr]) => {
            const { numIterates, tMin, tMax } = conv.estimateAmB(aSqr, bSqr);
            expect(numIterates).toBeLessThanOrEqual(MAX_ITERATIONS);
            expect(leDiff(tMin, aSqr, bSqr)).toBe(true);
            expect(geDiff(tMax, aSqr, bSqr)).toBe(true);
            expect(tMax.sub(tMin).lessThan(eps)).toBe(true);
        }, 30);
    });

    it('brackets a quadratic field number x + y*sqrt(d)', () => {
        const conv = new APConversion(PRECISION, MAX_ITERATIONS);
        const eps = threshold(PRECISION);
        const coefficient: fc.Arbitrary<BSRational> =
            fc.tuple(fc.integer({ min: -20, max: 20 }),
                fc.integer({ min: 1, max: 8 }))
                .map(([n, d]) => BSRational.fromNumber(n, d));

        check(fc.tuple(coefficient, coefficient, positiveSqr),
            ([x, y, d]) => {
                const q: APConversionQFN1 = { x: [x, y], d };
                const { numIterates, qMin, qMax } = conv.estimate(q);
                expect(qMin.lessThanOrEqual(qMax)).toBe(true);

                if (y.getSign() === 0) {
                    expect(numIterates).toBe(0);
                    expect(qMin.equals(x)).toBe(true);
                    expect(qMax.equals(x)).toBe(true);
                    return;
                }

                // The exact value is x + sign(y)*sqrt(y^2*d).
                const aSqr = y.mul(y).mul(d);
                if (y.getSign() > 0) {
                    expect(leSqrt(qMin.sub(x), aSqr)).toBe(true);
                    expect(geSqrt(qMax.sub(x), aSqr)).toBe(true);
                } else {
                    // qMin <= x - sqrt(aSqr) <= qMax.
                    expect(geSqrt(x.sub(qMin), aSqr)).toBe(true);
                    expect(leSqrt(x.sub(qMax), aSqr)).toBe(true);
                }
                expect(qMax.sub(qMin).lessThan(eps)).toBe(true);

                const single = conv.estimateValue(q);
                expect(single.qEstimate.equals(
                    BSRational.ldexp(qMin.add(qMax), -1))).toBe(true);
            }, 30);
    });

    it('tightens the bracket as the precision grows', () => {
        const aSqr = BSRational.fromNumber(2);
        let previous: BSRational | null = null;
        for (const precision of [8, 16, 32, 64]) {
            const conv = new APConversion(precision, MAX_ITERATIONS);
            expect(conv.getPrecision()).toBe(precision);
            expect(conv.getMaxIterations()).toBe(MAX_ITERATIONS);
            const { aMin, aMax } = conv.estimateSqrt(aSqr);
            const width = aMax.sub(aMin);
            expect(width.lessThan(threshold(precision))).toBe(true);
            expect(leSqrt(aMin, aSqr)).toBe(true);
            expect(geSqrt(aMax, aSqr)).toBe(true);
            if (previous !== null) {
                expect(width.lessThanOrEqual(previous)).toBe(true);
            }
            previous = width;
        }

        // setPrecision replaces the threshold.
        const conv = new APConversion(4, MAX_ITERATIONS);
        conv.setPrecision(50);
        expect(conv.getPrecision()).toBe(50);
        const tight = conv.estimateSqrt(aSqr);
        expect(tight.aMax.sub(tight.aMin).lessThan(threshold(50))).toBe(true);
        conv.setMaxIterations(3);
        expect(conv.getMaxIterations()).toBe(3);
    });

    it('is exact for squares of dyadic rationals', () => {
        const conv = new APConversion(PRECISION, MAX_ITERATIONS);
        check(fc.tuple(fc.integer({ min: 1, max: 64 }),
            fc.integer({ min: -6, max: 6 })), ([n, e]) => {
                const a = BSRational.ldexp(BSRational.fromNumber(n), e);
                const { aMin, aMax } = conv.estimateSqrt(a.mul(a));
                // The exact root lies in the bracket and the bracket is tight.
                expect(aMin.lessThanOrEqual(a)).toBe(true);
                expect(a.lessThanOrEqual(aMax)).toBe(true);
            }, 40);
    });

    // Regression for the port fix of the upstream stale-square defect
    // (issue #280). When a^2/b^2 is close to (7 + 3*sqrt(5))/2 the quantity
    // f''(a-b)/8 = a^2 - 3ab + b^2 is tiny, so the 53-bit initial bounds fall
    // on the wrong side of the inflection point and the bisection loop of
    // estimateAmB runs. With few iterations the loop exits by exhaustion
    // rather than by its break, and upstream then feeds the Newton bound a
    // tMinSqr (respectively tMaxSqr) that no longer corresponds to tMin
    // (tMax). The resulting "bound" is not a bound: the exact a-b falls
    // outside [tMin, tMax]. The convergents of (7 + 3*sqrt(5))/2 below all
    // fail the bracket assertion on the pre-fix source.
    it('keeps the bracket when the bisection loop exhausts its iterations',
        () => {
            const convergents: [bigint, bigint][] = [
                [10983760033n, 1602508992n],
                [64300051206n, 9381251041n],
                [516002918640n, 75283811239n],
                [3020733700601n, 440719107401n],
                [141910183877041n, 20704416796806n],
                [6666757908520326n, 972666870342481n],
                [53500214605455696n, 7805576116155895n]
            ];
            for (let maxIterations = 1; maxIterations <= 8; ++maxIterations) {
                const conv = new APConversion(40, maxIterations);
                for (const [p, q] of convergents) {
                    const aSqr = BSRational.fromBigInt(p);
                    const bSqr = BSRational.fromBigInt(q);
                    const { tMin, tMax } = conv.estimateAmB(aSqr, bSqr);
                    expect(leDiff(tMin, aSqr, bSqr)).toBe(true);
                    expect(geDiff(tMax, aSqr, bSqr)).toBe(true);
                }
            }
        });

    it('rejects invalid construction arguments', () => {
        expect(() => new APConversion(0, 1)).toThrow(/Invalid precision/);
        expect(() => new APConversion(-1, 1)).toThrow(/Invalid precision/);
        expect(() => new APConversion(1, 0)).toThrow(/Invalid maximum/);
        const conv = new APConversion(8, 8);
        expect(() => conv.setPrecision(0)).toThrow(/Invalid precision/);
        expect(() => conv.setMaxIterations(0)).toThrow(/Invalid maximum/);
    });
});
