import { describe, it, expect } from 'vitest';
import { asinEstimate, getASinEstimateMaxError } from '../src/ASinEstimate.js';
import { acosEstimate, getACosEstimateMaxError } from '../src/ACosEstimate.js';
import { GTE_C_HALF_PI } from '../src/Constants.js';
import { sinEstimate } from '../src/SinEstimate.js';
import { check, fc } from './helpers/arbitraries.js';

const DEGREES = [1, 2, 3, 4, 5, 6, 7, 8] as const;

// Upstream-documented maximum errors of the estimates on [0,1].
const MAX_ERROR = [
    9.0128265558586e-3,  // degree 1
    8.1851275863202e-4,  // degree 2
    8.8200141836567e-5,  // degree 3
    1.0563052499871e-5,  // degree 4
    1.3535063235066e-6,  // degree 5
    1.8169471743823e-7,  // degree 6
    2.5231622315797e-8,  // degree 7
    3.5952707963527e-9   // degree 8
];

describe('asinEstimate', () => {
    it('stays within the documented max error on a dense grid of [0,1]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            let maxObserved = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = i / samples;
                const error = Math.abs(asinEstimate(x, degree) - Math.asin(x));
                if (error > maxObserved) {
                    maxObserved = error;
                }
            }
            expect(maxObserved).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-15);
            // The minimax bound is tight, so a dense grid must come close.
            expect(maxObserved).toBeGreaterThan(0.5 * bound);
        }
    });

    it('is exactly pi/2 - acosEstimate(x) for every degree', () => {
        for (const degree of DEGREES) {
            for (const x of [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1]) {
                expect(asinEstimate(x, degree))
                    .toBe(GTE_C_HALF_PI - acosEstimate(x, degree));
            }
        }
    });

    it('returns exactly pi/2 at x = 1, where asin(1) = pi/2', () => {
        // acosEstimate(1) = sqrt(1-1)*p(1) = 0 exactly, so the estimate is
        // the double closest to pi/2.
        for (const degree of DEGREES) {
            expect(asinEstimate(1, degree)).toBe(GTE_C_HALF_PI);
        }
    });

    it('returns 0 at x = 0 up to the constant-coefficient rounding', () => {
        // The degree-D polynomial's constant coefficient is the double
        // closest to pi/2, so the estimate at x = 0 is pi/2 - that value.
        for (const degree of DEGREES) {
            expect(asinEstimate(0, degree)).toBe(0);
        }
    });

    it('has the same error magnitude as the acos estimate (identity check)', () => {
        // asin(x) - asinEstimate(x) = acosEstimate(x) - acos(x).
        for (const degree of DEGREES) {
            for (let i = 0; i <= 64; ++i) {
                const x = i / 64;
                const asinErr = Math.asin(x) - asinEstimate(x, degree);
                const acosErr = acosEstimate(x, degree) - Math.acos(x);
                expect(Math.abs(asinErr - acosErr)).toBeLessThan(1e-15);
            }
        }
    });

    it('increases in accuracy with degree', () => {
        const worst = DEGREES.map((degree) => {
            let m = 0;
            for (let i = 0; i <= 512; ++i) {
                const x = i / 512;
                m = Math.max(m, Math.abs(asinEstimate(x, degree) - Math.asin(x)));
            }
            return m;
        });
        for (let i = 1; i < worst.length; ++i) {
            expect(worst[i]).toBeLessThan(worst[i - 1]);
        }
    });
});

describe('getASinEstimateMaxError', () => {
    it('returns the documented bounds', () => {
        for (const degree of DEGREES) {
            expect(getASinEstimateMaxError(degree)).toBe(MAX_ERROR[degree - 1]);
        }
    });

    it('agrees with the acos bounds to 7 significant digits', () => {
        // The two upstream tables are the same quantity; they differ only in
        // the last digits of the published decimals.
        for (const degree of DEGREES) {
            const a = getASinEstimateMaxError(degree);
            const b = getACosEstimateMaxError(degree);
            expect(Math.abs(a - b) / b).toBeLessThan(1e-7);
        }
    });

    it('rejects degrees outside [1,8] and non-integers', () => {
        expect(() => getASinEstimateMaxError(0)).toThrow(/Invalid degree/);
        expect(() => getASinEstimateMaxError(9)).toThrow(/Invalid degree/);
        expect(() => getASinEstimateMaxError(2.5)).toThrow(/Invalid degree/);
        expect(() => asinEstimate(0.5, 0)).toThrow(/Invalid degree/);
        expect(() => asinEstimate(0.5, 9)).toThrow(/Invalid degree/);
    });
});

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream ASinEstimate.h.
// ---------------------------------------------------------------------------

describe('ASinEstimate verification', () => {
    const inDomain = fc.double({ min: 0, max: 1, noNaN: true });

    it('is exactly pi/2 minus the acos estimate, for every input', () => {
        // Upstream is a one-line forward; the identity must be bit exact.
        for (const degree of DEGREES) {
            check(inDomain, x =>
                asinEstimate(x, degree)
                    === GTE_C_HALF_PI - acosEstimate(x, degree));
        }
    });

    it('is within the documented bound on fast-check samples of [0,1]', () => {
        for (const degree of DEGREES) {
            check(inDomain, x =>
                Math.abs(asinEstimate(x, degree) - Math.asin(x))
                    <= MAX_ERROR[degree - 1]);
        }
    });

    it('has a max error over [0,1] that decreases with the degree', () => {
        // Measured independently of the published table. asin and acos share
        // the same polynomial, so the two tables must agree to within the
        // rounding of the pi/2 constant, which this also confirms.
        const samples = 40000;
        let previous = Number.POSITIVE_INFINITY;
        for (const degree of DEGREES) {
            let observed = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = i / samples;
                observed = Math.max(observed,
                    Math.abs(asinEstimate(x, degree) - Math.asin(x)));
            }
            expect(observed).toBeLessThanOrEqual(MAX_ERROR[degree - 1]);
            expect(observed).toBeGreaterThan(0.9 * MAX_ERROR[degree - 1]);
            expect(observed).toBeLessThan(previous);
            previous = observed;
        }
    }, 30000);

    it('is strictly increasing on [0,1] for every degree', () => {
        const samples = 5000;
        for (const degree of DEGREES) {
            let previous = Number.NEGATIVE_INFINITY;
            for (let i = 0; i <= samples; ++i) {
                const v = asinEstimate(i / samples, degree);
                expect(v).toBeGreaterThan(previous);
                previous = v;
            }
        }
    }, 30000);

    it('sums with the acos estimate to exactly pi/2', () => {
        // The complementary identity holds bit exactly because the port
        // subtracts the same acos value it adds back.
        for (const degree of DEGREES) {
            check(inDomain, x =>
                asinEstimate(x, degree) + acosEstimate(x, degree)
                    === GTE_C_HALF_PI);
        }
    });

    it('matches sinEstimate as a left inverse to the combined accuracy', () => {
        // sin(asin(x)) = x through two independent files; the derivative of
        // sin at asin(x) is sqrt(1-x^2) <= 1, so the errors simply add.
        for (const degree of [6, 8]) {
            check(inDomain, x => {
                const a = asinEstimate(x, degree);
                const back = sinEstimate(a, 11);
                return Math.abs(back - x)
                    <= MAX_ERROR[degree - 1] + 2e-11 + 1e-12;
            });
        }
    });
});
