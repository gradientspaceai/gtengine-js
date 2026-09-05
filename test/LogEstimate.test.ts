import { describe, it, expect } from 'vitest';
import { logEstimate, logEstimateRR, getLogEstimateMaxError } from '../src/LogEstimate.js';
import { log2Estimate, log2EstimateRR, getLog2EstimateMaxError } from '../src/Log2Estimate.js';
import { GTE_C_LN_2 } from '../src/Constants.js';
import { expEstimateRR } from '../src/ExpEstimate.js';
import { check, fc } from './helpers/arbitraries.js';

const DEGREES = [1, 2, 3, 4, 5, 6, 7, 8] as const;

// Upstream forwards GetLogEstimateMaxError to GetLog2EstimateMaxError, so
// these log2 bounds are what the API reports.
const LOG2_MAX_ERROR = [
    8.6071332055935e-2,  // degree 1
    7.6362868906659e-3,  // degree 2
    8.7902902652948e-4,  // degree 3
    1.1318551356388e-4,  // degree 4
    1.5521274483455e-5,  // degree 5
    2.2162052037978e-6,  // degree 6
    3.2546558681457e-7,  // degree 7
    4.8798286744756e-8   // degree 8
];

describe('logEstimate', () => {
    it('stays within the tight log(2)-scaled bound on a dense grid of [1,2]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            // log(x) - logEstimate(x) = log(2)*(log2(x) - log2Estimate(x)),
            // so the true bound is log(2) times the reported log2 bound.
            const bound = LOG2_MAX_ERROR[degree - 1] * GTE_C_LN_2;
            let maxObserved = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = 1 + i / samples;
                const error = Math.abs(logEstimate(x, degree) - Math.log(x));
                if (error > maxObserved) {
                    maxObserved = error;
                }
            }
            expect(maxObserved).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-15);
            // The minimax bound is tight, so a dense grid must come close.
            expect(maxObserved).toBeGreaterThan(0.5 * bound);
        }
    });

    it('also satisfies the (looser) reported max error', () => {
        for (const degree of DEGREES) {
            const bound = getLogEstimateMaxError(degree);
            for (let i = 0; i <= 512; ++i) {
                const x = 1 + i / 512;
                expect(Math.abs(logEstimate(x, degree) - Math.log(x)))
                    .toBeLessThanOrEqual(bound);
            }
        }
    });

    it('is exactly log2Estimate(x)*log(2)', () => {
        for (const degree of DEGREES) {
            for (const x of [1, 1.125, 1.5, 1.75, 2]) {
                expect(logEstimate(x, degree)).toBe(log2Estimate(x, degree) * GTE_C_LN_2);
            }
        }
    });

    it('is exactly 0 at x = 1 for every degree', () => {
        // The polynomial is evaluated at t = x - 1 = 0 and multiplied by t.
        for (const degree of DEGREES) {
            expect(logEstimate(1, degree)).toBe(0);
        }
    });

    it('is monotone increasing on [1,2] for every degree', () => {
        for (const degree of DEGREES) {
            let previous = logEstimate(1, degree);
            for (let i = 1; i <= 1024; ++i) {
                const value = logEstimate(1 + i / 1024, degree);
                expect(value).toBeGreaterThan(previous);
                previous = value;
            }
        }
    });

    it('increases in accuracy with degree', () => {
        const worst = DEGREES.map((degree) => {
            let m = 0;
            for (let i = 0; i <= 512; ++i) {
                const x = 1 + i / 512;
                m = Math.max(m, Math.abs(logEstimate(x, degree) - Math.log(x)));
            }
            return m;
        });
        for (let i = 1; i < worst.length; ++i) {
            expect(worst[i]).toBeLessThan(worst[i - 1]);
        }
    });
});

describe('logEstimateRR', () => {
    it('agrees with logEstimate on the reduction interval [1,2)', () => {
        for (const degree of DEGREES) {
            for (const x of [1, 1.25, 1.5, 1.9999]) {
                expect(logEstimateRR(x, degree)).toBeCloseTo(logEstimate(x, degree), 15);
            }
        }
    });

    it('is exactly log2EstimateRR(x)*log(2)', () => {
        for (const degree of DEGREES) {
            for (const x of [1e-8, 0.3, 1, 7, 1024, 1e12]) {
                expect(logEstimateRR(x, degree)).toBe(log2EstimateRR(x, degree) * GTE_C_LN_2);
            }
        }
    });

    it('matches log(x) within the scaled bound over many decades', () => {
        for (const degree of DEGREES) {
            // The log2 range reduction adds an exact integer exponent, so the
            // error is that of the [1,2] estimate scaled by log(2).
            const bound = LOG2_MAX_ERROR[degree - 1] * GTE_C_LN_2;
            for (let e = -40; e <= 40; ++e) {
                for (const m of [1, 1.3, 1.7, 1.95]) {
                    const x = m * Math.pow(2, e);
                    const error = Math.abs(logEstimateRR(x, degree) - Math.log(x));
                    expect(error).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-12);
                }
            }
        }
    });

    it('reproduces powers of two exactly up to the estimate at y = 1', () => {
        // For x = 2^p, the reduction yields y = 1 with logEstimate(1) = 0, so
        // the result is exactly p*log(2).
        for (const degree of DEGREES) {
            for (const p of [-10, -1, 0, 1, 5, 30]) {
                expect(logEstimateRR(Math.pow(2, p), degree)).toBe(p * GTE_C_LN_2);
            }
        }
    });
});

describe('getLogEstimateMaxError', () => {
    it('forwards the log2 bounds, as upstream does', () => {
        for (const degree of DEGREES) {
            expect(getLogEstimateMaxError(degree)).toBe(LOG2_MAX_ERROR[degree - 1]);
            expect(getLogEstimateMaxError(degree)).toBe(getLog2EstimateMaxError(degree));
        }
    });

    it('rejects degrees outside [1,8] and non-integers', () => {
        expect(() => getLogEstimateMaxError(0)).toThrow(/Invalid degree/);
        expect(() => getLogEstimateMaxError(9)).toThrow(/Invalid degree/);
        expect(() => getLogEstimateMaxError(3.5)).toThrow(/Invalid degree/);
        expect(() => logEstimate(1.5, 0)).toThrow(/Invalid degree/);
        expect(() => logEstimateRR(1.5, 9)).toThrow(/Invalid degree/);
    });
});

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream LogEstimate.h.
// ---------------------------------------------------------------------------

describe('LogEstimate verification', () => {
    const anyPositive = fc.tuple(
        fc.double({ min: 1, max: 2, noNaN: true, noDefaultInfinity: true }),
        fc.integer({ min: -1080, max: 1020 })
    ).map(([m, e]) => m * Math.pow(2, e)).filter(x => x > 0 && Number.isFinite(x));

    it('is exactly the log2 estimate scaled by log(2), for every input', () => {
        // Upstream is a one-line forward in both variants.
        for (const degree of DEGREES) {
            check(anyPositive, x =>
                logEstimate(x, degree) === log2Estimate(x, degree) * GTE_C_LN_2
                && logEstimateRR(x, degree)
                    === log2EstimateRR(x, degree) * GTE_C_LN_2);
        }
    });

    it('has an error of exactly log(2) times the log2 error', () => {
        // log(x) - logEstimate(x) = log(2) * (log2(x) - log2Estimate(x)), so
        // the reported (log2) bound is loose by the factor log(2). This pins
        // the upstream choice of forwarding GetLogEstimateMaxError.
        const inDomain = fc.double({ min: 1, max: 2, noNaN: true });
        for (const degree of DEGREES) {
            check(inDomain, x => {
                const bound = LOG2_MAX_ERROR[degree - 1];
                const error = Math.abs(logEstimate(x, degree) - Math.log(x));
                return error <= GTE_C_LN_2 * bound * (1 + 1e-9) + 1e-16
                    && error <= bound;
            });
        }
    });

    it('has a max error over [1,2] that decreases with the degree', () => {
        const samples = 40000;
        let previous = Number.POSITIVE_INFINITY;
        for (const degree of DEGREES) {
            let observed = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = 1 + i / samples;
                observed = Math.max(observed,
                    Math.abs(logEstimate(x, degree) - Math.log(x)));
            }
            const tight = GTE_C_LN_2 * LOG2_MAX_ERROR[degree - 1];
            expect(observed).toBeLessThanOrEqual(tight * (1 + 1e-9));
            expect(observed).toBeGreaterThan(0.9 * tight);
            expect(observed).toBeLessThan(previous);
            previous = observed;
        }
    }, 30000);

    it('bounds the absolute error of the range-reduced form everywhere', () => {
        for (const degree of DEGREES) {
            check(anyPositive, x =>
                Math.abs(logEstimateRR(x, degree) - Math.log(x))
                    <= LOG2_MAX_ERROR[degree - 1] * (1 + 1e-9));
        }
    });

    it('inverts expEstimateRR to the accuracy of the two estimates', () => {
        // A cross-check against the other direction of the same identity:
        // log(exp(x)) = x, computed through two independent estimate files.
        const moderate = fc.double({ min: -20, max: 20, noNaN: true });
        check(moderate, x => {
            const e = expEstimateRR(x, 7);
            if (!(e > 0) || !Number.isFinite(e)) { return true; }
            const back = logEstimateRR(e, 8);
            // Relative error r of exp maps to an absolute error log(1+r) in
            // the logarithm, plus the logarithm's own absolute bound.
            return Math.abs(back - x)
                <= LOG2_MAX_ERROR[7] + 6.8e-11 + 1e-12;
        });
    });

    it('turns products into sums to the accuracy of the estimate', () => {
        const moderate = fc.double({ min: 1e-30, max: 1e30, noNaN: true })
            .filter(x => x > 0);
        for (const degree of [3, 8]) {
            check(fc.tuple(moderate, moderate), ([a, b]) => {
                const product = a * b;
                if (!Number.isFinite(product) || product === 0) { return true; }
                const lhs = logEstimateRR(product, degree);
                const rhs = logEstimateRR(a, degree) + logEstimateRR(b, degree);
                return Math.abs(lhs - rhs)
                    <= 3 * LOG2_MAX_ERROR[degree - 1] + 1e-9;
            });
        }
    });
});
