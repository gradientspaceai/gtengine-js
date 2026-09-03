import { describe, it, expect } from 'vitest';
import {
    exp2Estimate, exp2EstimateRR, getExp2EstimateMaxError
} from '../src/Exp2Estimate.js';

const DEGREES = [1, 2, 3, 4, 5, 6, 7] as const;

// Upstream-documented maximum errors of the minimax polynomials on [0,1].
const MAX_ERROR = [
    8.6071332055935e-2,   // degree 1
    3.8132476831059e-3,   // degree 2
    1.4694877755229e-4,   // degree 3
    4.7617792662269e-6,   // degree 4
    1.3162098788655e-7,   // degree 5
    3.1590552396211e-9,   // degree 6
    6.7157390759576e-11   // degree 7
];

describe('exp2Estimate', () => {
    it('stays within the documented max error on a dense grid of [0,1]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            let maxObserved = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = i / samples;
                const error = Math.abs(exp2Estimate(x, degree) - Math.pow(2, x));
                if (error > maxObserved) {
                    maxObserved = error;
                }
            }
            // A tiny cushion accounts for rounding in the polynomial
            // evaluation itself; the documented bound is for exact
            // arithmetic.
            expect(maxObserved).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-15);
            // The minimax bound is tight (equioscillation), so a dense grid
            // must come close to it.
            expect(maxObserved).toBeGreaterThan(0.5 * bound);
        }
    });

    it('is exact at x = 0 for every degree', () => {
        // The constant coefficient is exactly 1, and Horner evaluation at
        // x = 0 returns it unchanged, so 2^0 = 1 is reproduced exactly.
        for (const degree of DEGREES) {
            expect(exp2Estimate(0, degree)).toBe(1);
        }
    });

    it('approximates 2^1 = 2 at x = 1 within the max error', () => {
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            expect(Math.abs(exp2Estimate(1, degree) - 2)).toBeLessThanOrEqual(
                bound * (1 + 1e-8) + 1e-15);
        }
        // Degree 1 interpolates the endpoints: p(x) = 1 + x, so p(1) = 2
        // exactly.
        expect(exp2Estimate(1, 1)).toBe(2);
    });

    it('matches an explicit Horner evaluation of the degree-3 coefficients', () => {
        // Independent spot check with the upstream degree-3 coefficients.
        const c = [
            1.0,
            6.9589012084456225e-1,
            2.2486494900110188e-1,
            7.9244930154334980e-2
        ];
        for (const x of [0, 0.125, 0.5, 0.75, 1]) {
            const expected = c[0] + x * (c[1] + x * (c[2] + x * c[3]));
            expect(exp2Estimate(x, 3)).toBe(expected);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => exp2Estimate(0.5, 0)).toThrow('Invalid degree.');
        expect(() => exp2Estimate(0.5, 8)).toThrow('Invalid degree.');
        expect(() => exp2Estimate(0.5, 2.5)).toThrow('Invalid degree.');
    });
});

describe('exp2EstimateRR', () => {
    it('stays within the scaled max error over [-10, 10]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            for (let i = 0; i <= samples; ++i) {
                const x = -10 + (20 * i) / samples;
                const error = Math.abs(exp2EstimateRR(x, degree) - Math.pow(2, x));
                // |ldexp(p(y), n) - 2^x| = 2^n |p(y) - 2^y| <= 2^n * bound
                // for n = floor(x).
                const scaledBound = Math.pow(2, Math.floor(x)) * bound;
                expect(error).toBeLessThanOrEqual(scaledBound * (1 + 1e-8) + 1e-300);
            }
        }
    });

    it('is exact at integer inputs for every degree', () => {
        // At integer x the reduced argument is y = 0, the polynomial value
        // is exactly 1, and ldexp(1, x) is an exact power of two.
        for (const degree of DEGREES) {
            for (let n = -100; n <= 100; n += 10) {
                expect(exp2EstimateRR(n, degree)).toBe(Math.pow(2, n));
            }
        }
    });

    it('agrees with exp2Estimate on [0,1)', () => {
        for (const degree of DEGREES) {
            for (const x of [0, 0.25, 0.5, 0.999]) {
                expect(exp2EstimateRR(x, degree)).toBe(exp2Estimate(x, degree));
            }
        }
    });

    it('handles extreme exponents like ldexp', () => {
        for (const degree of DEGREES) {
            // Overflow to infinity.
            expect(exp2EstimateRR(1100, degree)).toBe(Infinity);
            // Underflow: subnormal results and eventually zero.
            expect(exp2EstimateRR(-1074, degree)).toBe(Math.pow(2, -1074));
            expect(exp2EstimateRR(-1200, degree)).toBe(0);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => exp2EstimateRR(1.5, 0)).toThrow('Invalid degree.');
        expect(() => exp2EstimateRR(1.5, 8)).toThrow('Invalid degree.');
    });
});

describe('getExp2EstimateMaxError', () => {
    it('returns the documented bounds', () => {
        for (const degree of DEGREES) {
            expect(getExp2EstimateMaxError(degree)).toBe(MAX_ERROR[degree - 1]);
        }
    });

    it('is strictly decreasing in the degree', () => {
        for (let degree = 2; degree <= 7; ++degree) {
            expect(getExp2EstimateMaxError(degree))
                .toBeLessThan(getExp2EstimateMaxError(degree - 1));
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => getExp2EstimateMaxError(0)).toThrow('Invalid degree.');
        expect(() => getExp2EstimateMaxError(8)).toThrow('Invalid degree.');
    });
});
