import { describe, it, expect } from 'vitest';
import {
    expEstimate, expEstimateRR, getExpEstimateMaxError
} from '../src/ExpEstimate.js';
import { exp2Estimate, exp2EstimateRR } from '../src/Exp2Estimate.js';
import { GTE_C_INV_LN_2 } from '../src/Constants.js';
import { check, fc } from './helpers/arbitraries.js';

const DEGREES = [1, 2, 3, 4, 5, 6, 7] as const;

// Upstream-documented maximum errors of the estimates on [0,log(2)]. They
// differ slightly from the Exp2Estimate bounds at degrees 5-7 because the
// error is measured against exp rather than 2^x.
const MAX_ERROR = [
    8.6071332055935e-2,  // degree 1
    3.8132476831059e-3,  // degree 2
    1.4694877755229e-4,  // degree 3
    4.7617792662269e-6,  // degree 4
    1.3162098766451e-7,  // degree 5
    3.1590550175765e-9,  // degree 6
    6.7157168714971e-11  // degree 7
];

const LN_2 = Math.LN2;

describe('expEstimate', () => {
    it('stays within the documented max error on a dense grid of [0,log(2)]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            let maxObserved = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = (LN_2 * i) / samples;
                const error = Math.abs(expEstimate(x, degree) - Math.exp(x));
                if (error > maxObserved) {
                    maxObserved = error;
                }
            }
            expect(maxObserved).toBeLessThanOrEqual(bound * (1 + 1e-6) + 1e-15);
            // The minimax bound is tight, so a dense grid must come close.
            expect(maxObserved).toBeGreaterThan(0.5 * bound);
        }
    });

    it('is exact at x = 0 for every degree', () => {
        for (const degree of DEGREES) {
            expect(expEstimate(0, degree)).toBe(1);
        }
    });

    it('delegates to exp2Estimate via the identity exp(x) = 2^(x/ln 2)', () => {
        for (const degree of DEGREES) {
            for (const x of [0, 0.1, 0.35, LN_2]) {
                expect(expEstimate(x, degree))
                    .toBe(exp2Estimate(x * GTE_C_INV_LN_2, degree));
            }
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => expEstimate(0.5, 0)).toThrow('Invalid degree.');
        expect(() => expEstimate(0.5, 8)).toThrow('Invalid degree.');
        expect(() => expEstimate(0.5, 2.5)).toThrow('Invalid degree.');
    });
});

describe('expEstimateRR', () => {
    it('delegates to exp2EstimateRR via the identity exp(x) = 2^(x/ln 2)', () => {
        for (const degree of DEGREES) {
            for (const x of [-3.7, -1, 0, 0.25, 1, 5.5]) {
                expect(expEstimateRR(x, degree))
                    .toBe(exp2EstimateRR(x * GTE_C_INV_LN_2, degree));
            }
        }
    });

    it('stays within the scaled max error over [-10,10]', () => {
        const samples = 2048;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            for (let i = 0; i <= samples; ++i) {
                const x = -10 + (20 * i) / samples;
                const u = x * GTE_C_INV_LN_2;
                const error = Math.abs(expEstimateRR(x, degree) - Math.exp(x));
                // |ldexp(p(y), n) - 2^u| <= 2^n * bound for n = floor(u); a
                // relative cushion absorbs the rounding of u = x/ln(2).
                const scaledBound = Math.pow(2, Math.floor(u)) * bound;
                expect(error).toBeLessThanOrEqual(
                    scaledBound * (1 + 1e-6) + 4e-15 * Math.exp(x) + 1e-300);
            }
        }
    });

    it('approximates powers of two at multiples of log(2)', () => {
        // x = n*ln(2) maps to u near the integer n. Because ln(2) is
        // rounded, u may land just below n, in which case the polynomial is
        // evaluated near the endpoint y = 1; the relative error is still
        // bounded by the documented max error.
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            for (const n of [-8, -1, 1, 8]) {
                const x = n * LN_2;
                const relError = Math.abs(expEstimateRR(x, degree) / Math.pow(2, n) - 1);
                expect(relError).toBeLessThanOrEqual(bound * (1 + 1e-6) + 1e-13);
            }
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => expEstimateRR(1, 0)).toThrow('Invalid degree.');
        expect(() => expEstimateRR(1, 8)).toThrow('Invalid degree.');
    });
});

describe('getExpEstimateMaxError', () => {
    it('returns the documented bounds', () => {
        for (const degree of DEGREES) {
            expect(getExpEstimateMaxError(degree)).toBe(MAX_ERROR[degree - 1]);
        }
    });

    it('is strictly decreasing in the degree', () => {
        for (let degree = 2; degree <= 7; ++degree) {
            expect(getExpEstimateMaxError(degree))
                .toBeLessThan(getExpEstimateMaxError(degree - 1));
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => getExpEstimateMaxError(0)).toThrow('Invalid degree.');
        expect(() => getExpEstimateMaxError(8)).toThrow('Invalid degree.');
    });
});

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream ExpEstimate.h.
// ---------------------------------------------------------------------------

describe('ExpEstimate verification', () => {
    it('is exactly the exp2 estimate of x/log(2), for every input', () => {
        // Upstream is a one-line forward; a stray extra rounding would show.
        const wide = fc.double({ min: -600, max: 600, noNaN: true });
        for (const degree of DEGREES) {
            check(wide, x =>
                expEstimate(x, degree) === exp2Estimate(x * GTE_C_INV_LN_2, degree)
                && expEstimateRR(x, degree)
                    === exp2EstimateRR(x * GTE_C_INV_LN_2, degree));
        }
    });

    it('is within the documented bound on fast-check samples of [0,log(2)]', () => {
        const inDomain = fc.double({ min: 0, max: LN_2, noNaN: true });
        for (const degree of DEGREES) {
            check(inDomain, x =>
                Math.abs(expEstimate(x, degree) - Math.exp(x))
                    <= MAX_ERROR[degree - 1]);
        }
    });

    it('has a max error over [0,log(2)] that decreases with the degree', () => {
        // Measured independently of the published table.
        const samples = 40000;
        let previous = Number.POSITIVE_INFINITY;
        for (const degree of DEGREES) {
            let observed = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = (i * LN_2) / samples;
                observed = Math.max(observed,
                    Math.abs(expEstimate(x, degree) - Math.exp(x)));
            }
            expect(observed).toBeLessThanOrEqual(MAX_ERROR[degree - 1]);
            expect(observed).toBeGreaterThan(0.9 * MAX_ERROR[degree - 1]);
            expect(observed).toBeLessThan(previous);
            previous = observed;
        }
    }, 30000);

    it('bounds the relative error of the range-reduced form', () => {
        // exp(x) = 2^(x/ln 2) and the exp2 reduction keeps the relative error
        // at the tabulated absolute bound of the polynomial on [0,1]. The
        // extra rounding of x*(1/ln 2) contributes |x|*eps to the exponent,
        // hence a relative term of |x|*ln(2)*eps.
        const wide = fc.double({ min: -200, max: 200, noNaN: true });
        for (const degree of DEGREES) {
            check(wide, x => {
                const exact = Math.exp(x);
                const relative = Math.abs(expEstimateRR(x, degree) - exact)
                    / exact;
                const cushion = Math.abs(x) * Number.EPSILON;
                return relative <= MAX_ERROR[degree - 1] + cushion + 1e-15;
            });
        }
    });

    it('is strictly increasing over a wide range for every degree', () => {
        const samples = 4000;
        for (const degree of DEGREES) {
            let previous = Number.NEGATIVE_INFINITY;
            for (let i = 0; i <= samples; ++i) {
                const v = expEstimateRR(-50 + (100 * i) / samples, degree);
                expect(v).toBeGreaterThan(previous);
                previous = v;
            }
        }
    }, 30000);

    it('keeps exp(a)*exp(b) close to exp(a+b) at the top degree', () => {
        // A cross-check that does not restate the implementation: the
        // multiplicative law holds to the relative accuracy of the estimate.
        const small = fc.double({ min: -5, max: 5, noNaN: true });
        check(fc.tuple(small, small), ([a, b]) => {
            const lhs = expEstimateRR(a, 7) * expEstimateRR(b, 7);
            const rhs = expEstimateRR(a + b, 7);
            return Math.abs(lhs - rhs) <= 3 * MAX_ERROR[6] * Math.abs(rhs)
                + 1e-13 * Math.abs(rhs);
        });
    });
});
