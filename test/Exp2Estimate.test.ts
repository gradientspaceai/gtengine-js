import { describe, it, expect } from 'vitest';
import {
    exp2Estimate, exp2EstimateRR, getExp2EstimateMaxError
} from '../src/Exp2Estimate.js';
import { check, fc } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream Exp2Estimate.h.
// ---------------------------------------------------------------------------

/**
 * std::frexp by repeated halving: x = f * 2^p with f in [1/2,1). Every step
 * is exact (scaling by two never rounds, subnormals included), so this is an
 * independent reference for the bit-twiddling reduction inside the port.
 */
function frexpByLoop(x: number): { f: number; p: number } {
    if (x === 0 || !Number.isFinite(x)) { return { f: x, p: 0 }; }
    let f = x;
    let p = 0;
    while (Math.abs(f) >= 1) { f /= 2; ++p; }
    while (Math.abs(f) < 0.5) { f *= 2; --p; }
    return { f, p };
}

describe('Exp2Estimate verification', () => {
    const inDomain = fc.double({ min: 0, max: 1, noNaN: true });

    it('is within the documented bound on fast-check samples of [0,1]', () => {
        for (const degree of DEGREES) {
            check(inDomain, x =>
                Math.abs(exp2Estimate(x, degree) - Math.pow(2, x))
                    <= MAX_ERROR[degree - 1]);
        }
    });

    it('has a max error over [0,1] that decreases with the degree', () => {
        // Measured here rather than read from the table: the RotationEstimate
        // tables were found to understate the real error for high degrees.
        const samples = 40000;
        let previous = Number.POSITIVE_INFINITY;
        for (const degree of DEGREES) {
            let observed = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = i / samples;
                observed = Math.max(observed,
                    Math.abs(exp2Estimate(x, degree) - Math.pow(2, x)));
            }
            expect(observed).toBeLessThanOrEqual(MAX_ERROR[degree - 1]);
            expect(observed).toBeGreaterThan(0.9 * MAX_ERROR[degree - 1]);
            expect(observed).toBeLessThan(previous);
            previous = observed;
        }
    }, 30000);

    it('is strictly increasing on [0,1] for every degree', () => {
        // 2^x is increasing and the minimax error is far below the increment
        // over one grid step, so a coefficient transposition would show up.
        const samples = 5000;
        for (const degree of DEGREES) {
            let previous = Number.NEGATIVE_INFINITY;
            for (let i = 0; i <= samples; ++i) {
                const v = exp2Estimate(i / samples, degree);
                expect(v).toBeGreaterThan(previous);
                previous = v;
            }
        }
    }, 30000);

    it('splits the range reduction into floor and fraction, as upstream does', () => {
        // Upstream: p = floor(x), y = x - p, result = ldexp(poly(y), p). The
        // reference multiplies by an exactly computed power of two instead of
        // the port's two-step scaling.
        const wide = fc.double({ min: -900, max: 900, noNaN: true });
        for (const degree of [1, 4, 7]) {
            check(wide, x => {
                const p = Math.floor(x);
                const poly = exp2Estimate(x - p, degree);
                const expected = poly * Math.pow(2, p);
                return exp2EstimateRR(x, degree) === expected;
            });
        }
    });

    it('has a relative error bounded by the absolute error of the polynomial', () => {
        // 2^x = 2^p * 2^y with 2^y >= 1, so the reduction cannot amplify the
        // relative error beyond the tabulated absolute bound.
        const wide = fc.double({ min: -300, max: 300, noNaN: true });
        for (const degree of DEGREES) {
            check(wide, x => {
                const exact = Math.pow(2, x);
                const relative = Math.abs(exp2EstimateRR(x, degree) - exact)
                    / exact;
                return relative <= MAX_ERROR[degree - 1] * (1 + 1e-12);
            });
        }
    });

    it('reproduces powers of two exactly, including subnormal ones', () => {
        // y = 0 makes the polynomial exactly its constant term 1.
        for (const degree of DEGREES) {
            for (let p = -1074; p <= 1023; p += 7) {
                expect(exp2EstimateRR(p, degree)).toBe(Math.pow(2, p));
            }
        }
    }, 30000);

    it('agrees with an independent frexp on the exponent it produces', () => {
        // The estimate of 2^x lies in the binade selected by floor(x), which
        // is the property the ldexp step has to preserve.
        const wide = fc.double({ min: -500, max: 500, noNaN: true });
        for (const degree of [3, 7]) {
            check(wide, x => {
                const value = exp2EstimateRR(x, degree);
                const { p } = frexpByLoop(value);
                const expected = frexpByLoop(Math.pow(2, x)).p;
                return Math.abs(p - expected) <= 1;
            });
        }
    });
});
