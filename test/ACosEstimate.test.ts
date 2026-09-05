import { describe, it, expect } from 'vitest';
import { acosEstimate, getACosEstimateMaxError } from '../src/ACosEstimate.js';
import { cosEstimate } from '../src/CosEstimate.js';
import { check, fc } from './helpers/arbitraries.js';

const DEGREES = [1, 2, 3, 4, 5, 6, 7, 8] as const;

// Upstream-documented maximum errors of the estimates on [0,1].
const MAX_ERROR = [
    9.0128265558585e-3,  // degree 1
    8.1851275863199e-4,  // degree 2
    8.8200141836526e-5,  // degree 3
    1.0563052499802e-5,  // degree 4
    1.3535063234649e-6,  // degree 5
    1.8169471727170e-7,  // degree 6
    2.5231622347022e-8,  // degree 7
    3.5952707477805e-9   // degree 8
];

describe('acosEstimate', () => {
    it('stays within the documented max error on a dense grid of [0,1]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            let maxObserved = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = i / samples;
                const error = Math.abs(acosEstimate(x, degree) - Math.acos(x));
                if (error > maxObserved) {
                    maxObserved = error;
                }
            }
            expect(maxObserved).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-15);
            // The minimax bound is tight, so a dense grid must come close.
            expect(maxObserved).toBeGreaterThan(0.5 * bound);
        }
    });

    it('is exact at x = 1 for every degree', () => {
        // f(1) = sqrt(1-1)*p(1) = 0 = acos(1) exactly.
        for (const degree of DEGREES) {
            expect(acosEstimate(1, degree)).toBe(0);
        }
    });

    it('returns exactly pi/2 at x = 0 for every degree', () => {
        // The constant coefficient is the double closest to pi/2 and
        // sqrt(1-0) = 1, so f(0) = coeff[0] exactly.
        for (const degree of DEGREES) {
            expect(acosEstimate(0, degree)).toBe(1.5707963267948966);
        }
    });

    it('matches an explicit evaluation of the degree-2 coefficients', () => {
        const c = [
            +1.5707963267948966,
            -2.0347053865798365e-1,
            +4.6887774236182234e-2
        ];
        for (const x of [0, 0.25, 0.5, 0.75, 1]) {
            const expected = (c[0] + x * (c[1] + x * c[2])) * Math.sqrt(1 - x);
            expect(acosEstimate(x, 2)).toBe(expected);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => acosEstimate(0.5, 0)).toThrow('Invalid degree.');
        expect(() => acosEstimate(0.5, 9)).toThrow('Invalid degree.');
        expect(() => acosEstimate(0.5, 2.5)).toThrow('Invalid degree.');
    });
});

describe('getACosEstimateMaxError', () => {
    it('returns the documented bounds', () => {
        for (const degree of DEGREES) {
            expect(getACosEstimateMaxError(degree)).toBe(MAX_ERROR[degree - 1]);
        }
    });

    it('is strictly decreasing in the degree', () => {
        for (let degree = 2; degree <= 8; ++degree) {
            expect(getACosEstimateMaxError(degree))
                .toBeLessThan(getACosEstimateMaxError(degree - 1));
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => getACosEstimateMaxError(0)).toThrow('Invalid degree.');
        expect(() => getACosEstimateMaxError(9)).toThrow('Invalid degree.');
    });
});

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream ACosEstimate.h.
// ---------------------------------------------------------------------------

describe('ACosEstimate verification', () => {
    const inDomain = fc.double({ min: 0, max: 1, noNaN: true });

    it('is within the documented bound on fast-check samples of [0,1]', () => {
        for (const degree of DEGREES) {
            check(inDomain, x =>
                Math.abs(acosEstimate(x, degree) - Math.acos(x))
                    <= MAX_ERROR[degree - 1]);
        }
    });

    it('has a max error over [0,1] that decreases with the degree', () => {
        // Measured here rather than read from the table: the RotationEstimate
        // tables were found to understate the true error at high degrees, so
        // every table in this group is re-derived.
        const samples = 40000;
        let previous = Number.POSITIVE_INFINITY;
        for (const degree of DEGREES) {
            let observed = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = i / samples;
                observed = Math.max(observed,
                    Math.abs(acosEstimate(x, degree) - Math.acos(x)));
            }
            expect(observed).toBeLessThanOrEqual(MAX_ERROR[degree - 1]);
            expect(observed).toBeGreaterThan(0.9 * MAX_ERROR[degree - 1]);
            expect(observed).toBeLessThan(previous);
            previous = observed;
        }
    }, 30000);

    it('is strictly decreasing on [0,1] for every degree', () => {
        // acos decreases from pi/2 to 0; the sqrt(1-x) factor makes this a
        // property of the whole expression, not just the polynomial.
        const samples = 5000;
        for (const degree of DEGREES) {
            let previous = Number.POSITIVE_INFINITY;
            for (let i = 0; i <= samples; ++i) {
                const v = acosEstimate(i / samples, degree);
                expect(v).toBeLessThan(previous);
                previous = v;
            }
        }
    }, 30000);

    it('carries the sqrt(1-x) factor, so the estimate vanishes only at x = 1', () => {
        // The upstream form is sqrt(1-x)*p(x), which is what makes the
        // estimate exact at the endpoint where acos has a square-root
        // singularity; dropping the factor would leave a nonzero value there.
        for (const degree of DEGREES) {
            expect(acosEstimate(1, degree)).toBe(0);
            check(fc.double({ min: 0, max: 0.999, noNaN: true }), x =>
                acosEstimate(x, degree) > 0);
        }
    });

    it('has a ratio to sqrt(1-x) that stays bounded near the endpoint', () => {
        // acos(x)/sqrt(1-x) -> sqrt(2) as x -> 1, which is what the minimax
        // polynomial approximates; the ratio must not blow up.
        for (const degree of DEGREES) {
            check(inDomain.filter(x => x < 1), x => {
                const ratio = acosEstimate(x, degree) / Math.sqrt(1 - x);
                return ratio > 1.3 && ratio < 1.6;
            });
        }
    });

    it('matches cosEstimate as a left inverse to the combined accuracy', () => {
        // cos(acos(x)) = x, computed through two independent files. cos has
        // slope -sin(a) at a = acos(x), so the error transfers with the
        // factor sqrt(1-x^2) <= 1.
        for (const degree of [6, 8]) {
            check(inDomain, x => {
                const a = acosEstimate(x, degree);
                const back = cosEstimate(a, 10);
                return Math.abs(back - x)
                    <= MAX_ERROR[degree - 1] + 2.8e-10 + 1e-12;
            });
        }
    });
});
