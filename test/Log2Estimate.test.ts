import { describe, it, expect } from 'vitest';
import {
    log2Estimate, log2EstimateRR, getLog2EstimateMaxError
} from '../src/Log2Estimate.js';
import { check, fc } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream Log2Estimate.h.
// ---------------------------------------------------------------------------

/**
 * std::frexp by repeated halving: x = f * 2^p with f in [1/2,1). Scaling by
 * two is always exact, subnormals included, so this is an independent
 * reference for the exponent-field extraction inside the port.
 */
function frexpByLoop(x: number): { f: number; p: number } {
    if (x === 0 || !Number.isFinite(x)) { return { f: x, p: 0 }; }
    let f = x;
    let p = 0;
    while (Math.abs(f) >= 1) { f /= 2; ++p; }
    while (Math.abs(f) < 0.5) { f *= 2; --p; }
    return { f, p };
}

describe('Log2Estimate verification', () => {
    // Positive doubles spread over every binade, subnormals included.
    const anyPositive = fc.tuple(
        fc.double({ min: 1, max: 2, noNaN: true, noDefaultInfinity: true }),
        fc.integer({ min: -1080, max: 1020 })
    ).map(([m, e]) => m * Math.pow(2, e)).filter(x => x > 0 && Number.isFinite(x));

    it('is within the documented bound on fast-check samples of [1,2]', () => {
        const inDomain = fc.double({ min: 1, max: 2, noNaN: true });
        for (const degree of DEGREES) {
            check(inDomain, x =>
                Math.abs(log2Estimate(x, degree) - Math.log2(x))
                    <= MAX_ERROR[degree - 1]);
        }
    });

    it('has a max error over [1,2] that decreases with the degree', () => {
        // Measured independently of the published table.
        const samples = 40000;
        let previous = Number.POSITIVE_INFINITY;
        for (const degree of DEGREES) {
            let observed = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = 1 + i / samples;
                observed = Math.max(observed,
                    Math.abs(log2Estimate(x, degree) - Math.log2(x)));
            }
            expect(observed).toBeLessThanOrEqual(MAX_ERROR[degree - 1]);
            expect(observed).toBeGreaterThan(0.9 * MAX_ERROR[degree - 1]);
            expect(observed).toBeLessThan(previous);
            previous = observed;
        }
    }, 30000);

    it('is strictly increasing on [1,2] for every degree', () => {
        const samples = 5000;
        for (const degree of DEGREES) {
            let previous = Number.NEGATIVE_INFINITY;
            for (let i = 0; i <= samples; ++i) {
                const v = log2Estimate(1 + i / samples, degree);
                expect(v).toBeGreaterThan(previous);
                previous = v;
            }
        }
    }, 30000);

    it('reduces with frexp exactly, as upstream does', () => {
        // Upstream: y = 2*frexp(x, &p), --p, result = poly(y) + p. The
        // reference frexp here halves in a loop instead of reading the
        // exponent field.
        for (const degree of DEGREES) {
            check(anyPositive, x => {
                const { f, p } = frexpByLoop(x);
                const y = 2 * f;
                return log2EstimateRR(x, degree)
                    === log2Estimate(y, degree) + (p - 1);
            });
        }
    });

    it('bounds the absolute error of the range-reduced form everywhere', () => {
        // The reduction contributes an exact integer, so the absolute error
        // of log2EstimateRR is the absolute error of the polynomial on [1,2].
        for (const degree of DEGREES) {
            check(anyPositive, x =>
                Math.abs(log2EstimateRR(x, degree) - Math.log2(x))
                    <= MAX_ERROR[degree - 1] * (1 + 1e-12));
        }
    });

    it('turns products into sums to the accuracy of the estimate', () => {
        // log2(a*b) = log2(a) + log2(b): an identity the implementation does
        // not encode anywhere, so it cross-checks the reduction and the
        // polynomial together.
        const moderate = fc.double({ min: 1e-30, max: 1e30, noNaN: true })
            .filter(x => x > 0);
        for (const degree of [3, 8]) {
            check(fc.tuple(moderate, moderate), ([a, b]) => {
                const product = a * b;
                if (!Number.isFinite(product) || product === 0) { return true; }
                const lhs = log2EstimateRR(product, degree);
                const rhs = log2EstimateRR(a, degree) + log2EstimateRR(b, degree);
                return Math.abs(lhs - rhs) <= 3 * MAX_ERROR[degree - 1] + 1e-9;
            });
        }
    });

    it('is exact at every power of two, subnormals included', () => {
        for (const degree of DEGREES) {
            for (let p = -1074; p <= 1023; p += 3) {
                expect(log2EstimateRR(Math.pow(2, p), degree)).toBe(p);
            }
        }
    }, 30000);
});
