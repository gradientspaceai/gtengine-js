import { describe, it, expect } from 'vitest';
import { acosEstimate, getACosEstimateMaxError } from '../src/ACosEstimate.js';

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
