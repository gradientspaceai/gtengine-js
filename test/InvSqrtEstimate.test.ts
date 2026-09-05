import { describe, it, expect } from 'vitest';
import {
    invSqrtEstimate, invSqrtEstimateRR, getInvSqrtEstimateMaxError
} from '../src/InvSqrtEstimate.js';
import { GTE_C_INV_SQRT_2, GTE_C_SQRT_2 } from '../src/Constants.js';
import { sqrtEstimateRR } from '../src/SqrtEstimate.js';
import { check, fc } from './helpers/arbitraries.js';

const DEGREES = [1, 2, 3, 4, 5, 6, 7, 8] as const;

// Upstream-documented maximum errors of the estimates on [1,2].
const MAX_ERROR = [
    3.7814314552702e-2,  // degree 1
    4.1953446330581e-3,  // degree 2
    5.6307702007275e-4,  // degree 3
    8.1513919990229e-5,  // degree 4
    1.2289367490981e-5,  // degree 5
    1.9001451476708e-6,  // degree 6
    2.9887737629242e-7,  // degree 7
    4.7597402907940e-8   // degree 8
];

describe('invSqrtEstimate', () => {
    it('stays within the documented max error on a dense grid of [1,2]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            let maxObserved = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = 1 + i / samples;
                const error = Math.abs(invSqrtEstimate(x, degree) - 1 / Math.sqrt(x));
                if (error > maxObserved) {
                    maxObserved = error;
                }
            }
            expect(maxObserved).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-15);
            // The minimax bound is tight, so a dense grid must come close.
            expect(maxObserved).toBeGreaterThan(0.5 * bound);
        }
    });

    it('is exact at x = 1 for every degree (constant term 1)', () => {
        for (const degree of DEGREES) {
            expect(invSqrtEstimate(1, degree)).toBe(1);
        }
    });

    it('matches an explicit Horner evaluation of the degree-2 coefficients', () => {
        const c = [
            +1.0,
            -4.4539812104566801e-1,
            +1.5250490223221547e-1
        ];
        for (const x of [1, 1.25, 1.5, 1.75, 2]) {
            const t = x - 1;
            const expected = c[0] + t * (c[1] + t * c[2]);
            expect(invSqrtEstimate(x, 2)).toBe(expected);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => invSqrtEstimate(1.5, 0)).toThrow('Invalid degree.');
        expect(() => invSqrtEstimate(1.5, 9)).toThrow('Invalid degree.');
        expect(() => invSqrtEstimate(1.5, 2.5)).toThrow('Invalid degree.');
    });
});

describe('invSqrtEstimateRR', () => {
    it('agrees with invSqrtEstimate on [1,2)', () => {
        for (const degree of DEGREES) {
            for (const x of [1, 1.25, 1.5, 1.999]) {
                expect(invSqrtEstimateRR(x, degree)).toBe(invSqrtEstimate(x, degree));
            }
        }
    });

    it('stays within the relative max error over a wide range of x > 0', () => {
        // 1/sqrt(x) = adj * 2^p / sqrt(y) with y in [1,2), 1/sqrt(y) in
        // (1/sqrt(2), 1], so |est - exact|/exact <= bound * sqrt(2).
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            for (let e = -40; e <= 40; ++e) {
                for (const frac of [1, 1.1, 1.5, 1.9]) {
                    const x = frac * Math.pow(2, e);
                    const exact = 1 / Math.sqrt(x);
                    const relError = Math.abs(invSqrtEstimateRR(x, degree) / exact - 1);
                    expect(relError).toBeLessThanOrEqual(
                        Math.SQRT2 * bound * (1 + 1e-8) + 1e-15);
                }
            }
        }
    });

    it('is exact at powers of four for every degree', () => {
        // x = 4^k reduces to y = 1 with an even exponent, so the exact
        // power of two 2^-k is produced.
        for (const degree of DEGREES) {
            for (const k of [-10, -1, 0, 1, 2, 10]) {
                expect(invSqrtEstimateRR(Math.pow(4, k), degree))
                    .toBe(Math.pow(2, -k));
            }
        }
    });

    it('returns the rounded 1/sqrt(2) constant at x = 2', () => {
        // x = 2 reduces to y = 1 with an odd exponent, so the result is
        // exactly the adjustment constant.
        for (const degree of DEGREES) {
            expect(invSqrtEstimateRR(2, degree)).toBe(GTE_C_INV_SQRT_2);
        }
    });

    it('handles subnormal inputs', () => {
        const degree = 8;
        const bound = MAX_ERROR[degree - 1];
        for (const x of [Math.pow(2, -1040), Math.pow(2, -1060)]) {
            const exact = 1 / Math.sqrt(x);
            const relError = Math.abs(invSqrtEstimateRR(x, degree) / exact - 1);
            expect(relError).toBeLessThanOrEqual(Math.SQRT2 * bound * (1 + 1e-8) + 1e-15);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => invSqrtEstimateRR(2, 0)).toThrow('Invalid degree.');
        expect(() => invSqrtEstimateRR(2, 9)).toThrow('Invalid degree.');
    });
});

describe('getInvSqrtEstimateMaxError', () => {
    it('returns the documented bounds', () => {
        for (const degree of DEGREES) {
            expect(getInvSqrtEstimateMaxError(degree)).toBe(MAX_ERROR[degree - 1]);
        }
    });

    it('is strictly decreasing in the degree', () => {
        for (let degree = 2; degree <= 8; ++degree) {
            expect(getInvSqrtEstimateMaxError(degree))
                .toBeLessThan(getInvSqrtEstimateMaxError(degree - 1));
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => getInvSqrtEstimateMaxError(0)).toThrow('Invalid degree.');
        expect(() => getInvSqrtEstimateMaxError(9)).toThrow('Invalid degree.');
    });
});

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream InvSqrtEstimate.h.
// ---------------------------------------------------------------------------

/** std::frexp by repeated halving; every step is exact. */
function frexpByLoop(x: number): { f: number; p: number } {
    if (x === 0 || !Number.isFinite(x)) { return { f: x, p: 0 }; }
    let f = x;
    let p = 0;
    while (Math.abs(f) >= 1) { f /= 2; ++p; }
    while (Math.abs(f) < 0.5) { f *= 2; --p; }
    return { f, p };
}

describe('InvSqrtEstimate verification', () => {
    const anyPositive = fc.tuple(
        fc.double({ min: 1, max: 2, noNaN: true, noDefaultInfinity: true }),
        fc.integer({ min: -1020, max: 1020 })
    ).map(([m, e]) => m * Math.pow(2, e)).filter(x => x > 0 && Number.isFinite(x));

    it('is within the documented bound on fast-check samples of [1,2]', () => {
        const inDomain = fc.double({ min: 1, max: 2, noNaN: true });
        for (const degree of DEGREES) {
            check(inDomain, x =>
                Math.abs(invSqrtEstimate(x, degree) - 1 / Math.sqrt(x))
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
                    Math.abs(invSqrtEstimate(x, degree) - 1 / Math.sqrt(x)));
            }
            expect(observed).toBeLessThanOrEqual(MAX_ERROR[degree - 1]);
            expect(observed).toBeGreaterThan(0.9 * MAX_ERROR[degree - 1]);
            expect(observed).toBeLessThan(previous);
            previous = observed;
        }
    }, 30000);

    it('is strictly decreasing on [1,2] for every degree', () => {
        const samples = 5000;
        for (const degree of DEGREES) {
            let previous = Number.POSITIVE_INFINITY;
            for (let i = 0; i <= samples; ++i) {
                const v = invSqrtEstimate(1 + i / samples, degree);
                expect(v).toBeLessThan(previous);
                previous = v;
            }
        }
    }, 30000);

    it('splits the exponent with the odd/even adjustment upstream uses', () => {
        // Upstream: y = 2*frexp(x,&p), --p, adj = (1&p)/sqrt(2) + (1&~p),
        // p = -(p >> 1), result = adj*ldexp(poly(y), p). Note the negation of
        // the *shifted* exponent: `p >> 1` floors, so odd negative exponents
        // do not behave like a truncating division.
        for (const degree of DEGREES) {
            check(anyPositive, x => {
                const { f, p: pf } = frexpByLoop(x);
                const y = 2 * f;
                const p = pf - 1;
                const adj = Math.abs(p % 2) === 1 ? GTE_C_INV_SQRT_2 : 1;
                const half = -Math.floor(p / 2);
                const expected = adj
                    * (invSqrtEstimate(y, degree) * Math.pow(2, half));
                return invSqrtEstimateRR(x, degree) === expected;
            });
        }
    });

    it('bounds the relative error of the range-reduced form everywhere', () => {
        // 1/sqrt(x) = adj*2^(-p/2)/sqrt(y) with 1/sqrt(y) >= 1/sqrt(2), so
        // the relative error is at most sqrt(2) times the tabulated bound.
        for (const degree of DEGREES) {
            check(anyPositive, x => {
                const exact = 1 / Math.sqrt(x);
                return Math.abs(invSqrtEstimateRR(x, degree) - exact) / exact
                    <= GTE_C_SQRT_2 * MAX_ERROR[degree - 1] * (1 + 1e-12);
            });
        }
    });

    it('is the reciprocal of sqrtEstimateRR to the combined accuracy', () => {
        // Cross-check between two files of this group: the product of the two
        // estimates approximates 1.
        for (const degree of [4, 8]) {
            check(anyPositive.filter(x => x > 1e-150 && x < 1e150), x => {
                const product = sqrtEstimateRR(x, degree)
                    * invSqrtEstimateRR(x, degree);
                return Math.abs(product - 1)
                    <= 3 * MAX_ERROR[degree - 1] + 1e-15;
            });
        }
    });

    it('preserves the scaling law invsqrt(4x) = invsqrt(x)/2 exactly', () => {
        for (const degree of DEGREES) {
            check(anyPositive.filter(x => x > 1e-280 && x < 1e280), x =>
                invSqrtEstimateRR(4 * x, degree)
                    === 0.5 * invSqrtEstimateRR(x, degree));
        }
    });
});
