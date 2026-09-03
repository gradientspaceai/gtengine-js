import { describe, it, expect } from 'vitest';
import {
    log2Estimate, log2EstimateRR, getLog2EstimateMaxError
} from '../src/Log2Estimate.js';

const DEGREES = [1, 2, 3, 4, 5, 6, 7, 8] as const;

// Upstream-documented maximum errors of the estimates on [1,2].
const MAX_ERROR = [
    8.6071332055935e-2,  // degree 1
    7.6362868906659e-3,  // degree 2
    8.7902902652948e-4,  // degree 3
    1.1318551356388e-4,  // degree 4
    1.5521274483455e-5,  // degree 5
    2.2162052037978e-6,  // degree 6
    3.2546558681457e-7,  // degree 7
    4.8798286744756e-8   // degree 8
];

describe('log2Estimate', () => {
    it('stays within the documented max error on a dense grid of [1,2]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            let maxObserved = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = 1 + i / samples;
                const error = Math.abs(log2Estimate(x, degree) - Math.log2(x));
                if (error > maxObserved) {
                    maxObserved = error;
                }
            }
            expect(maxObserved).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-15);
            // The minimax bound is tight, so a dense grid must come close.
            expect(maxObserved).toBeGreaterThan(0.5 * bound);
        }
    });

    it('is exact at x = 1 for every degree (p has no constant term)', () => {
        for (const degree of DEGREES) {
            expect(log2Estimate(1, degree)).toBe(0);
        }
    });

    it('is exact at x = 2 for degree 1 (p(t) = t)', () => {
        expect(log2Estimate(2, 1)).toBe(1);
    });

    it('matches an explicit Horner evaluation of the degree-3 coefficients', () => {
        const c = [
            +1.4228653756681227,
            -5.8208556916449616e-1,
            +1.5922019349637218e-1
        ];
        for (const x of [1, 1.25, 1.5, 1.75, 2]) {
            const t = x - 1;
            const expected = (c[0] + t * (c[1] + t * c[2])) * t;
            expect(log2Estimate(x, 3)).toBe(expected);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => log2Estimate(1.5, 0)).toThrow('Invalid degree.');
        expect(() => log2Estimate(1.5, 9)).toThrow('Invalid degree.');
        expect(() => log2Estimate(1.5, 2.5)).toThrow('Invalid degree.');
    });
});

describe('log2EstimateRR', () => {
    it('agrees with log2Estimate on [1,2)', () => {
        for (const degree of DEGREES) {
            for (const x of [1, 1.25, 1.5, 1.999]) {
                expect(log2EstimateRR(x, degree)).toBe(log2Estimate(x, degree));
            }
        }
    });

    it('stays within the absolute max error over a wide range of x > 0', () => {
        // log2(x) = log2(y) + p with the reduction exact, so the polynomial
        // bound applies unchanged; a tiny cushion covers the final addition
        // of the (possibly large) integer exponent.
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            for (let e = -60; e <= 60; e += 3) {
                for (const frac of [1, 1.2, 1.7, 1.95]) {
                    const x = frac * Math.pow(2, e);
                    const error = Math.abs(log2EstimateRR(x, degree) - Math.log2(x));
                    expect(error).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-13);
                }
            }
        }
    });

    it('is exact at powers of two for every degree', () => {
        // x = 2^k reduces to y = 1 where the polynomial vanishes exactly,
        // leaving the integer exponent k.
        for (const degree of DEGREES) {
            for (const k of [-100, -7, -1, 0, 1, 7, 100]) {
                expect(log2EstimateRR(Math.pow(2, k), degree)).toBe(k);
            }
        }
    });

    it('is exact for subnormal powers of two', () => {
        for (const degree of DEGREES) {
            expect(log2EstimateRR(Math.pow(2, -1060), degree)).toBe(-1060);
            expect(log2EstimateRR(Math.pow(2, -1074), degree)).toBe(-1074);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => log2EstimateRR(3, 0)).toThrow('Invalid degree.');
        expect(() => log2EstimateRR(3, 9)).toThrow('Invalid degree.');
    });
});

describe('getLog2EstimateMaxError', () => {
    it('returns the documented bounds', () => {
        for (const degree of DEGREES) {
            expect(getLog2EstimateMaxError(degree)).toBe(MAX_ERROR[degree - 1]);
        }
    });

    it('is strictly decreasing in the degree', () => {
        for (let degree = 2; degree <= 8; ++degree) {
            expect(getLog2EstimateMaxError(degree))
                .toBeLessThan(getLog2EstimateMaxError(degree - 1));
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => getLog2EstimateMaxError(0)).toThrow('Invalid degree.');
        expect(() => getLog2EstimateMaxError(9)).toThrow('Invalid degree.');
    });
});
