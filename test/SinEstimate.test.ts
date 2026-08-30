import { describe, it, expect } from 'vitest';
import {
    sinEstimate, sinEstimateRR, getSinEstimateMaxError
} from '../src/SinEstimate';
import { GTE_C_HALF_PI } from '../src/Constants';

const DEGREES = [3, 5, 7, 9, 11] as const;

// Upstream-documented maximum errors of the estimates on [-pi/2,pi/2].
const MAX_ERROR: Record<number, number> = {
    3: 1.3481903639146e-2,
    5: 1.4001209384651e-4,
    7: 1.0205878939740e-6,
    9: 5.2010783457846e-9,
    11: 1.9323431743601e-11
};

describe('sinEstimate', () => {
    it('stays within the documented max error on a dense grid of [-pi/2,pi/2]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree];
            let maxObserved = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = -GTE_C_HALF_PI + (2 * GTE_C_HALF_PI * i) / samples;
                const error = Math.abs(sinEstimate(x, degree) - Math.sin(x));
                if (error > maxObserved) {
                    maxObserved = error;
                }
            }
            expect(maxObserved).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-15);
            // The minimax bound is tight, so a dense grid must come close.
            expect(maxObserved).toBeGreaterThan(0.5 * bound);
        }
    });

    it('is exact at x = 0 for every degree', () => {
        for (const degree of DEGREES) {
            expect(sinEstimate(0, degree)).toBe(0);
        }
    });

    it('is an odd function (exact sign symmetry)', () => {
        // Only odd powers appear: p(x) = x*q(x^2), so p(-x) = -p(x) exactly.
        for (const degree of DEGREES) {
            for (const x of [0.1, 0.5, 1, 1.5]) {
                expect(sinEstimate(-x, degree)).toBe(-sinEstimate(x, degree));
            }
        }
    });

    it('satisfies the constraint p(pi/2) = 1 up to roundoff', () => {
        for (const degree of DEGREES) {
            expect(Math.abs(sinEstimate(GTE_C_HALF_PI, degree) - 1))
                .toBeLessThanOrEqual(1e-13);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => sinEstimate(0.5, 1)).toThrow('Invalid degree.');
        expect(() => sinEstimate(0.5, 4)).toThrow('Invalid degree.');
        expect(() => sinEstimate(0.5, 13)).toThrow('Invalid degree.');
        expect(() => sinEstimate(0.5, 3.5)).toThrow('Invalid degree.');
    });
});

describe('sinEstimateRR', () => {
    it('agrees with sinEstimate on [-pi/2,pi/2]', () => {
        // remainder(x, 2*pi) = x exactly there, so the same branch is taken.
        for (const degree of DEGREES) {
            for (const x of [-GTE_C_HALF_PI, -1, 0, 0.5, GTE_C_HALF_PI]) {
                expect(sinEstimateRR(x, degree)).toBe(sinEstimate(x, degree));
            }
        }
    });

    it('stays within the max error over a wide range of inputs', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree];
            for (let i = 0; i <= samples; ++i) {
                const x = -20 + (40 * i) / samples;
                const error = Math.abs(sinEstimateRR(x, degree) - Math.sin(x));
                // Small additive cushion for the rounded pi constants used
                // in the range reduction.
                expect(error).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-13);
            }
        }
    });

    it('respects periodicity: matches values shifted by 2*pi within tolerance', () => {
        for (const degree of DEGREES) {
            for (const x of [0.3, 1.2, 2.9]) {
                const base = sinEstimateRR(x, degree);
                const shifted = sinEstimateRR(x + 2 * Math.PI, degree);
                expect(Math.abs(shifted - base)).toBeLessThanOrEqual(1e-13);
            }
        }
    });

    it('is nearly zero at multiples of pi', () => {
        for (const degree of DEGREES) {
            expect(Math.abs(sinEstimateRR(Math.PI, degree))).toBeLessThanOrEqual(1e-13);
            expect(Math.abs(sinEstimateRR(-3 * Math.PI, degree))).toBeLessThanOrEqual(1e-13);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => sinEstimateRR(1, 2)).toThrow('Invalid degree.');
        expect(() => sinEstimateRR(1, 13)).toThrow('Invalid degree.');
    });
});

describe('getSinEstimateMaxError', () => {
    it('returns the documented bounds', () => {
        for (const degree of DEGREES) {
            expect(getSinEstimateMaxError(degree)).toBe(MAX_ERROR[degree]);
        }
    });

    it('is strictly decreasing in the degree', () => {
        for (let i = 1; i < DEGREES.length; ++i) {
            expect(getSinEstimateMaxError(DEGREES[i]))
                .toBeLessThan(getSinEstimateMaxError(DEGREES[i - 1]));
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => getSinEstimateMaxError(1)).toThrow('Invalid degree.');
        expect(() => getSinEstimateMaxError(13)).toThrow('Invalid degree.');
    });
});
