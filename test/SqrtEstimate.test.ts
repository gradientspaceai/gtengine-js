import { describe, it, expect } from 'vitest';
import {
    sqrtEstimate, sqrtEstimateRR, getSqrtEstimateMaxError
} from '../src/SqrtEstimate.js';
import { GTE_C_SQRT_2 } from '../src/Constants.js';
import { check, fc } from './helpers/arbitraries.js';

const DEGREES = [1, 2, 3, 4, 5, 6, 7, 8] as const;

// Upstream-documented maximum errors of the estimates on [1,2].
const MAX_ERROR = [
    1.7766952966369e-2,  // degree 1
    1.1795695163111e-3,  // degree 2
    1.1309620116485e-4,  // degree 3
    1.2741170151820e-5,  // degree 4
    1.5725569051384e-6,  // degree 5
    2.0584162152560e-7,  // degree 6
    2.8072338675856e-8,  // degree 7
    3.9468401880072e-9   // degree 8
];

describe('sqrtEstimate', () => {
    it('stays within the documented max error on a dense grid of [1,2]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            let maxObserved = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = 1 + i / samples;
                const error = Math.abs(sqrtEstimate(x, degree) - Math.sqrt(x));
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
            expect(sqrtEstimate(1, degree)).toBe(1);
        }
    });

    it('matches an explicit Horner evaluation of the degree-1 coefficients', () => {
        // p(t) = 1 + t*(sqrt(2)-1), interpolating sqrt at both endpoints.
        const c1 = 4.1421356237309505e-1;
        for (const x of [1, 1.25, 1.5, 2]) {
            const t = x - 1;
            expect(sqrtEstimate(x, 1)).toBe(1 + c1 * t);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => sqrtEstimate(1.5, 0)).toThrow('Invalid degree.');
        expect(() => sqrtEstimate(1.5, 9)).toThrow('Invalid degree.');
        expect(() => sqrtEstimate(1.5, 2.5)).toThrow('Invalid degree.');
    });
});

describe('sqrtEstimateRR', () => {
    it('agrees with sqrtEstimate on [1,2)', () => {
        for (const degree of DEGREES) {
            for (const x of [1, 1.25, 1.5, 1.999]) {
                expect(sqrtEstimateRR(x, degree)).toBe(sqrtEstimate(x, degree));
            }
        }
    });

    it('stays within the relative max error over a wide range of x > 0', () => {
        // sqrt(x) = adj * 2^p * sqrt(y) with y in [1,2), sqrt(y) >= 1, so
        // |est - exact|/exact <= bound.
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree - 1];
            for (let e = -40; e <= 40; ++e) {
                for (const frac of [1, 1.1, 1.5, 1.9]) {
                    const x = frac * Math.pow(2, e);
                    const exact = Math.sqrt(x);
                    const relError = Math.abs(sqrtEstimateRR(x, degree) / exact - 1);
                    expect(relError).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-15);
                }
            }
        }
    });

    it('is exact at powers of four for every degree', () => {
        // x = 4^k reduces to y = 1 with an even exponent, so the exact
        // power of two 2^k is produced.
        for (const degree of DEGREES) {
            for (const k of [-10, -1, 0, 1, 2, 10]) {
                expect(sqrtEstimateRR(Math.pow(4, k), degree)).toBe(Math.pow(2, k));
            }
        }
    });

    it('returns the rounded sqrt(2) constant at x = 2', () => {
        // x = 2 reduces to y = 1 with an odd exponent, so the result is
        // exactly the adjustment constant.
        for (const degree of DEGREES) {
            expect(sqrtEstimateRR(2, degree)).toBe(GTE_C_SQRT_2);
        }
    });

    it('handles subnormal inputs', () => {
        const degree = 8;
        const bound = MAX_ERROR[degree - 1];
        for (const x of [Math.pow(2, -1040), Math.pow(2, -1060)]) {
            const exact = Math.sqrt(x);
            const relError = Math.abs(sqrtEstimateRR(x, degree) / exact - 1);
            expect(relError).toBeLessThanOrEqual(bound * (1 + 1e-8) + 1e-15);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => sqrtEstimateRR(2, 0)).toThrow('Invalid degree.');
        expect(() => sqrtEstimateRR(2, 9)).toThrow('Invalid degree.');
    });
});

describe('getSqrtEstimateMaxError', () => {
    it('returns the documented bounds', () => {
        for (const degree of DEGREES) {
            expect(getSqrtEstimateMaxError(degree)).toBe(MAX_ERROR[degree - 1]);
        }
    });

    it('is strictly decreasing in the degree', () => {
        for (let degree = 2; degree <= 8; ++degree) {
            expect(getSqrtEstimateMaxError(degree))
                .toBeLessThan(getSqrtEstimateMaxError(degree - 1));
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => getSqrtEstimateMaxError(0)).toThrow('Invalid degree.');
        expect(() => getSqrtEstimateMaxError(9)).toThrow('Invalid degree.');
    });
});

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream SqrtEstimate.h.
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

describe('SqrtEstimate verification', () => {
    // Positive doubles across the binades the reference below can rescale
    // without overflowing (the port itself has no such restriction).
    const anyPositive = fc.tuple(
        fc.double({ min: 1, max: 2, noNaN: true, noDefaultInfinity: true }),
        fc.integer({ min: -1020, max: 1020 })
    ).map(([m, e]) => m * Math.pow(2, e)).filter(x => x > 0 && Number.isFinite(x));

    it('is within the documented bound on fast-check samples of [1,2]', () => {
        const inDomain = fc.double({ min: 1, max: 2, noNaN: true });
        for (const degree of DEGREES) {
            check(inDomain, x =>
                Math.abs(sqrtEstimate(x, degree) - Math.sqrt(x))
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
                    Math.abs(sqrtEstimate(x, degree) - Math.sqrt(x)));
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
                const v = sqrtEstimate(1 + i / samples, degree);
                expect(v).toBeGreaterThan(previous);
                previous = v;
            }
        }
    }, 30000);

    it('splits the exponent with the odd/even adjustment upstream uses', () => {
        // Upstream: y = 2*frexp(x,&p), --p, adj = (1&p)*sqrt(2) + (1&~p),
        // p >>= 1, result = adj*ldexp(poly(y), p). The reference restates the
        // integer bit tricks as arithmetic on the exponent: `1 & p` is the
        // parity of p (including for negative p) and `p >> 1` is the floor of
        // p/2, not truncation toward zero.
        for (const degree of DEGREES) {
            check(anyPositive, x => {
                const { f, p: pf } = frexpByLoop(x);
                const y = 2 * f;
                const p = pf - 1;
                const adj = Math.abs(p % 2) === 1 ? GTE_C_SQRT_2 : 1;
                const half = Math.floor(p / 2);
                const expected = adj
                    * (sqrtEstimate(y, degree) * Math.pow(2, half));
                return sqrtEstimateRR(x, degree) === expected;
            });
        }
    });

    it('bounds the relative error of the range-reduced form everywhere', () => {
        // sqrt(x) = adj*2^(p/2)*sqrt(y) with sqrt(y) >= 1, so the reduction
        // carries the absolute bound on [1,2] over as a relative bound.
        for (const degree of DEGREES) {
            check(anyPositive, x => {
                const exact = Math.sqrt(x);
                return Math.abs(sqrtEstimateRR(x, degree) - exact) / exact
                    <= MAX_ERROR[degree - 1] * (1 + 1e-12);
            });
        }
    });

    it('squares back to the input to twice the relative accuracy', () => {
        // An identity the implementation never uses: r^2 = x to within
        // 2*maxError relative, since (1+e)^2 = 1 + 2e + O(e^2).
        for (const degree of [4, 8]) {
            check(anyPositive, x => {
                const r = sqrtEstimateRR(x, degree);
                const sqr = r * r;
                if (!Number.isFinite(sqr) || sqr === 0) { return true; }
                return Math.abs(sqr - x) / x
                    <= 2.1 * MAX_ERROR[degree - 1] + 1e-15;
            });
        }
    });

    it('preserves the scaling law sqrt(4x) = 2*sqrt(x) exactly', () => {
        // Multiplying by four shifts p by two, which leaves the polynomial
        // argument and the parity adjustment untouched.
        for (const degree of DEGREES) {
            check(anyPositive.filter(x => x > 1e-280 && x < 1e280), x =>
                sqrtEstimateRR(4 * x, degree) === 2 * sqrtEstimateRR(x, degree));
        }
    });

    it('preserves the documented x = 0 quirk (upstream issue)', () => {
        // Upstream's frexp(0) yields p = 0, so the reduction evaluates the
        // polynomial at t = -1 and scales by sqrt(2)/2 - a meaningless
        // nonzero value, despite the documented x >= 0 constraint. Preserved.
        for (const degree of DEGREES) {
            const v = sqrtEstimateRR(0, degree);
            expect(Number.isFinite(v)).toBe(true);
            expect(v).not.toBe(0);
        }
    });
});
