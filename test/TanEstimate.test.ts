import { describe, it, expect } from 'vitest';
import {
    tanEstimate, tanEstimateRR, getTanEstimateMaxError
} from '../src/TanEstimate.js';
import { GTE_C_QUARTER_PI } from '../src/Constants.js';

const DEGREES = [3, 5, 7, 9, 11, 13] as const;

// Upstream-documented maximum errors of the estimates on [-pi/4,pi/4].
const MAX_ERROR: Record<number, number> = {
    3: 1.1661892256205e-2,
    5: 5.8431854390146e-4,
    7: 3.5418688397793e-5,
    9: 2.2988173248307e-6,
    11: 1.5426258070939e-7,
    13: 1.0550265105991e-8
};

describe('tanEstimate', () => {
    it('stays within the documented max error on a dense grid of [-pi/4,pi/4]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree];
            let maxObserved = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = -GTE_C_QUARTER_PI + (2 * GTE_C_QUARTER_PI * i) / samples;
                const error = Math.abs(tanEstimate(x, degree) - Math.tan(x));
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
            expect(tanEstimate(0, degree)).toBe(0);
        }
    });

    it('is an odd function (exact sign symmetry)', () => {
        // Only odd powers appear: p(x) = x*q(x^2), so p(-x) = -p(x) exactly.
        for (const degree of DEGREES) {
            for (const x of [0.1, 0.3, 0.5, GTE_C_QUARTER_PI]) {
                expect(tanEstimate(-x, degree)).toBe(-tanEstimate(x, degree));
            }
        }
    });

    it('satisfies the constraint p(pi/4) = 1 up to roundoff', () => {
        for (const degree of DEGREES) {
            expect(Math.abs(tanEstimate(GTE_C_QUARTER_PI, degree) - 1))
                .toBeLessThanOrEqual(1e-14);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => tanEstimate(0.5, 1)).toThrow('Invalid degree.');
        expect(() => tanEstimate(0.5, 4)).toThrow('Invalid degree.');
        expect(() => tanEstimate(0.5, 15)).toThrow('Invalid degree.');
        expect(() => tanEstimate(0.5, 3.5)).toThrow('Invalid degree.');
    });
});

describe('tanEstimateRR', () => {
    it('agrees with tanEstimate on [-pi/4,pi/4]', () => {
        // remainder(x, pi) = x exactly there, so the same branch is taken.
        for (const degree of DEGREES) {
            for (const x of [-GTE_C_QUARTER_PI, -0.5, 0, 0.5, GTE_C_QUARTER_PI]) {
                expect(tanEstimateRR(x, degree)).toBe(tanEstimate(x, degree));
            }
        }
    });

    it('stays within a derivative-scaled error bound over a wide range', () => {
        // For reduced |y| <= pi/4 the polynomial bound applies directly.
        // For |y| in (pi/4,pi/2) the identity tan(y) = (1+p)/(1-p) with
        // p = tan(y - pi/4) + e amplifies the polynomial error by
        // 2/(1-p)^2 = (1+tan(y)^2)/(1+p^2) <= 1 + tan(y)^2 = sec(y)^2, so
        // |est - tan(x)| <= bound * sec(x)^2 up to roundoff. Points too
        // close to the poles are skipped.
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree];
            for (let i = 0; i <= samples; ++i) {
                const x = -10 + (20 * i) / samples;
                if (Math.abs(Math.cos(x)) < 0.01) {
                    continue;
                }
                const tanx = Math.tan(x);
                const tolerance = bound * (1 + tanx * tanx) * (1 + 1e-6) + 1e-12;
                expect(Math.abs(tanEstimateRR(x, degree) - tanx))
                    .toBeLessThanOrEqual(tolerance);
            }
        }
    });

    it('uses the shift identities accurately at sample angles beyond pi/4', () => {
        // tan(pi/3) = sqrt(3), tan(-pi/3) = -sqrt(3).
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree];
            expect(Math.abs(tanEstimateRR(Math.PI / 3, degree) - Math.sqrt(3)))
                .toBeLessThanOrEqual(4 * bound * (1 + 1e-6) + 1e-13);
            expect(Math.abs(tanEstimateRR(-Math.PI / 3, degree) + Math.sqrt(3)))
                .toBeLessThanOrEqual(4 * bound * (1 + 1e-6) + 1e-13);
        }
    });

    it('respects periodicity: matches values shifted by pi within tolerance', () => {
        for (const degree of DEGREES) {
            for (const x of [0.3, 0.7, -0.5]) {
                const base = tanEstimateRR(x, degree);
                const shifted = tanEstimateRR(x + Math.PI, degree);
                expect(Math.abs(shifted - base)).toBeLessThanOrEqual(1e-12);
            }
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => tanEstimateRR(1, 2)).toThrow('Invalid degree.');
        expect(() => tanEstimateRR(1, 15)).toThrow('Invalid degree.');
    });
});

describe('getTanEstimateMaxError', () => {
    it('returns the documented bounds', () => {
        for (const degree of DEGREES) {
            expect(getTanEstimateMaxError(degree)).toBe(MAX_ERROR[degree]);
        }
    });

    it('is strictly decreasing in the degree', () => {
        for (let i = 1; i < DEGREES.length; ++i) {
            expect(getTanEstimateMaxError(DEGREES[i]))
                .toBeLessThan(getTanEstimateMaxError(DEGREES[i - 1]));
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => getTanEstimateMaxError(1)).toThrow('Invalid degree.');
        expect(() => getTanEstimateMaxError(15)).toThrow('Invalid degree.');
    });
});
