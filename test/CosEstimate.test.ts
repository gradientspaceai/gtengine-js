import { describe, it, expect } from 'vitest';
import {
    cosEstimate, cosEstimateRR, getCosEstimateMaxError
} from '../src/CosEstimate.js';
import { GTE_C_HALF_PI } from '../src/Constants.js';

const DEGREES = [2, 4, 6, 8, 10] as const;

// Upstream-documented maximum errors of the estimates on [-pi/2,pi/2].
const MAX_ERROR: Record<number, number> = {
    2: 5.6009595954128e-2,
    4: 9.1879932449727e-4,
    6: 9.2028470144446e-6,
    8: 5.9804535233743e-8,
    10: 2.7008567604626e-10
};

describe('cosEstimate', () => {
    it('stays within the documented max error on a dense grid of [-pi/2,pi/2]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree];
            let maxObserved = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = -GTE_C_HALF_PI + (2 * GTE_C_HALF_PI * i) / samples;
                const error = Math.abs(cosEstimate(x, degree) - Math.cos(x));
                if (error > maxObserved) {
                    maxObserved = error;
                }
            }
            expect(maxObserved).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-15);
            // The minimax bound is tight, so a dense grid must come close.
            expect(maxObserved).toBeGreaterThan(0.5 * bound);
        }
    });

    it('is exact at x = 0 for every degree (constant term 1)', () => {
        for (const degree of DEGREES) {
            expect(cosEstimate(0, degree)).toBe(1);
        }
    });

    it('is an even function (exact sign symmetry)', () => {
        // Only even powers appear: p(x) = q(x^2), so p(-x) = p(x) exactly.
        for (const degree of DEGREES) {
            for (const x of [0.1, 0.5, 1, 1.5]) {
                expect(cosEstimate(-x, degree)).toBe(cosEstimate(x, degree));
            }
        }
    });

    it('satisfies the constraint p(pi/2) = 0 up to roundoff', () => {
        for (const degree of DEGREES) {
            expect(Math.abs(cosEstimate(GTE_C_HALF_PI, degree)))
                .toBeLessThanOrEqual(1e-13);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => cosEstimate(0.5, 0)).toThrow('Invalid degree.');
        expect(() => cosEstimate(0.5, 3)).toThrow('Invalid degree.');
        expect(() => cosEstimate(0.5, 12)).toThrow('Invalid degree.');
        expect(() => cosEstimate(0.5, 4.5)).toThrow('Invalid degree.');
    });
});

describe('cosEstimateRR', () => {
    it('agrees with cosEstimate on [-pi/2,pi/2]', () => {
        // remainder(x, 2*pi) = x exactly there, so the same branch is taken.
        for (const degree of DEGREES) {
            for (const x of [-GTE_C_HALF_PI, -1, 0, 0.5, GTE_C_HALF_PI]) {
                expect(cosEstimateRR(x, degree)).toBe(cosEstimate(x, degree));
            }
        }
    });

    it('stays within the max error over a wide range of inputs', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree];
            for (let i = 0; i <= samples; ++i) {
                const x = -20 + (40 * i) / samples;
                const error = Math.abs(cosEstimateRR(x, degree) - Math.cos(x));
                // Small additive cushion for the rounded pi constants used
                // in the range reduction.
                expect(error).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-13);
            }
        }
    });

    it('respects periodicity: matches values shifted by 2*pi within tolerance', () => {
        for (const degree of DEGREES) {
            for (const x of [0.3, 1.2, 2.9]) {
                const base = cosEstimateRR(x, degree);
                const shifted = cosEstimateRR(x + 2 * Math.PI, degree);
                expect(Math.abs(shifted - base)).toBeLessThanOrEqual(1e-13);
            }
        }
    });

    it('handles the sign flip near pi correctly', () => {
        // cos(pi) = -1; the reduction maps r near pi through -cos(pi - r).
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree];
            expect(Math.abs(cosEstimateRR(Math.PI, degree) + 1))
                .toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-13);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => cosEstimateRR(1, 1)).toThrow('Invalid degree.');
        expect(() => cosEstimateRR(1, 12)).toThrow('Invalid degree.');
    });
});

describe('getCosEstimateMaxError', () => {
    it('returns the documented bounds', () => {
        for (const degree of DEGREES) {
            expect(getCosEstimateMaxError(degree)).toBe(MAX_ERROR[degree]);
        }
    });

    it('is strictly decreasing in the degree', () => {
        for (let i = 1; i < DEGREES.length; ++i) {
            expect(getCosEstimateMaxError(DEGREES[i]))
                .toBeLessThan(getCosEstimateMaxError(DEGREES[i - 1]));
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => getCosEstimateMaxError(0)).toThrow('Invalid degree.');
        expect(() => getCosEstimateMaxError(11)).toThrow('Invalid degree.');
    });
});
