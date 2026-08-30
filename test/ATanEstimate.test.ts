import { describe, it, expect } from 'vitest';
import {
    atanEstimate, atanEstimateRR, getATanEstimateMaxError
} from '../src/ATanEstimate';

const DEGREES = [3, 5, 7, 9, 11, 13] as const;

// Upstream-documented maximum errors of the estimates on [-1,1].
const MAX_ERROR: Record<number, number> = {
    3: 1.5970326392625e-2,
    5: 1.3509832247375e-3,
    7: 1.5051227215525e-4,
    9: 1.8921598624725e-5,
    11: 2.5477725020825e-6,
    13: 3.5859106295450e-7
};

describe('atanEstimate', () => {
    it('stays within the documented max error on a dense grid of [-1,1]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree];
            let maxObserved = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = -1 + (2 * i) / samples;
                const error = Math.abs(atanEstimate(x, degree) - Math.atan(x));
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
            expect(atanEstimate(0, degree)).toBe(0);
        }
    });

    it('is an odd function (exact sign symmetry)', () => {
        // Only odd powers appear: p(x) = x*q(x^2), so p(-x) = -p(x) exactly.
        for (const degree of DEGREES) {
            for (const x of [0.1, 0.25, 0.5, 0.99, 1]) {
                expect(atanEstimate(-x, degree)).toBe(-atanEstimate(x, degree));
            }
        }
    });

    it('satisfies the constraint p(1) = pi/4 up to roundoff', () => {
        for (const degree of DEGREES) {
            expect(Math.abs(atanEstimate(1, degree) - Math.PI / 4))
                .toBeLessThanOrEqual(1e-13);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => atanEstimate(0.5, 1)).toThrow('Invalid degree.');
        expect(() => atanEstimate(0.5, 4)).toThrow('Invalid degree.');
        expect(() => atanEstimate(0.5, 15)).toThrow('Invalid degree.');
        expect(() => atanEstimate(0.5, 3.5)).toThrow('Invalid degree.');
    });
});

describe('atanEstimateRR', () => {
    it('agrees with atanEstimate on [-1,1]', () => {
        for (const degree of DEGREES) {
            for (const x of [-1, -0.5, 0, 0.5, 1]) {
                expect(atanEstimateRR(x, degree)).toBe(atanEstimate(x, degree));
            }
        }
    });

    it('stays within the max error over a wide range of inputs', () => {
        // atan(x) = +-pi/2 - atan(1/x) maps |x| > 1 into the polynomial
        // domain, so the documented bound applies (plus a tiny cushion for
        // the rounded pi/2 constant and the division).
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree];
            for (let i = 0; i <= samples; ++i) {
                const x = -100 + (200 * i) / samples;
                const error = Math.abs(atanEstimateRR(x, degree) - Math.atan(x));
                expect(error).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-15);
            }
        }
    });

    it('approaches +-pi/2 for large magnitude inputs', () => {
        for (const degree of DEGREES) {
            expect(Math.abs(atanEstimateRR(1e12, degree) - Math.PI / 2))
                .toBeLessThanOrEqual(1e-11);
            expect(Math.abs(atanEstimateRR(-1e12, degree) + Math.PI / 2))
                .toBeLessThanOrEqual(1e-11);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => atanEstimateRR(2, 2)).toThrow('Invalid degree.');
        expect(() => atanEstimateRR(2, 15)).toThrow('Invalid degree.');
    });
});

describe('getATanEstimateMaxError', () => {
    it('returns the documented bounds', () => {
        for (const degree of DEGREES) {
            expect(getATanEstimateMaxError(degree)).toBe(MAX_ERROR[degree]);
        }
    });

    it('is strictly decreasing in the degree', () => {
        for (let i = 1; i < DEGREES.length; ++i) {
            expect(getATanEstimateMaxError(DEGREES[i]))
                .toBeLessThan(getATanEstimateMaxError(DEGREES[i - 1]));
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => getATanEstimateMaxError(1)).toThrow('Invalid degree.');
        expect(() => getATanEstimateMaxError(15)).toThrow('Invalid degree.');
    });
});
