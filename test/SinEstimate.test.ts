import { describe, it, expect } from 'vitest';
import {
    sinEstimate, sinEstimateRR, getSinEstimateMaxError
} from '../src/SinEstimate.js';
import { GTE_C_HALF_PI, GTE_C_PI, GTE_C_TWO_PI } from '../src/Constants.js';
import { check, fc } from './helpers/arbitraries.js';
import { exactRemainder } from './helpers/exact.js';

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

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream SinEstimate.h.
// ---------------------------------------------------------------------------

describe('SinEstimate verification', () => {
    // A generator that covers every binade rather than a bounded interval:
    // the range reduction is where the port's arithmetic differs most from
    // the C++ std::remainder it stands in for, and the difference only shows
    // up once |x| is large enough for the quotient x/(2*pi) to lose bits.
    const anyMagnitude = fc.tuple(
        fc.double({ min: 1, max: 2, noNaN: true, noDefaultInfinity: true }),
        fc.integer({ min: -40, max: 62 }),
        fc.boolean()
    ).map(([m, e, neg]) => (neg ? -m : m) * Math.pow(2, e));

    it('reduces the argument exactly, as std::remainder does', () => {
        // sinEstimateRR(x) must equal sinEstimate applied to the reduction of
        // x by the *double* 2*pi computed in infinite precision. Computing
        // the quotient in binary64 (x - Math.round(x/y)*y) is not the same
        // function: it is off by 1e-6 near |x| = 1e10 and leaves [-pi,pi]
        // entirely near |x| = 1e14, which pushes the polynomial outside the
        // interval it approximates.
        const reduce = (x: number): number => {
            const r = exactRemainder(x, GTE_C_TWO_PI);
            if (r > GTE_C_HALF_PI) { return GTE_C_PI - r; }
            if (r < -GTE_C_HALF_PI) { return -GTE_C_PI - r; }
            return r;
        };
        for (const degree of DEGREES) {
            check(anyMagnitude, x =>
                sinEstimateRR(x, degree) === sinEstimate(reduce(x), degree));
        }
        // Values that break the rounded-quotient reduction outright.
        for (const x of [1e10, 1e12, 1e14, 1e16, -9769346117865922,
            -26789498295071.3, 123456789012345, 4.5e15]) {
            for (const degree of DEGREES) {
                expect(sinEstimateRR(x, degree))
                    .toBe(sinEstimate(reduce(x), degree));
            }
        }
    });

    it('never leaves the range of a sine, whatever the magnitude of x', () => {
        // The reduced argument stays in [-pi/2,pi/2], where every estimate is
        // bounded by 1 + maxError. A reduction that overshoots the period
        // evaluates the polynomial outside its interval and overshoots 1.
        for (const degree of DEGREES) {
            check(anyMagnitude, x => {
                const v = sinEstimateRR(x, degree);
                return Math.abs(v) <= 1 + MAX_ERROR[degree];
            });
        }
    });

    it('is within the documented bound on fast-check samples of [-pi/2,pi/2]', () => {
        const inDomain = fc.double(
            { min: -GTE_C_HALF_PI, max: GTE_C_HALF_PI, noNaN: true });
        for (const degree of DEGREES) {
            check(inDomain, x =>
                Math.abs(sinEstimate(x, degree) - Math.sin(x))
                    <= MAX_ERROR[degree]);
        }
    });

    it('has a max error over the domain that decreases with the degree', () => {
        // Independently measured, not read from the table: the RotationEstimate
        // tables were found to understate the true error, so this checks the
        // published numbers rather than trusting them.
        const samples = 20000;
        let previous = Number.POSITIVE_INFINITY;
        for (const degree of DEGREES) {
            let observed = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = -GTE_C_HALF_PI + (2 * GTE_C_HALF_PI * i) / samples;
                observed = Math.max(observed,
                    Math.abs(sinEstimate(x, degree) - Math.sin(x)));
            }
            expect(observed).toBeLessThanOrEqual(MAX_ERROR[degree]);
            expect(observed).toBeGreaterThan(0.9 * MAX_ERROR[degree]);
            expect(observed).toBeLessThan(previous);
            previous = observed;
        }
    }, 30000);

    it('is exactly odd, for every input and degree', () => {
        for (const degree of DEGREES) {
            check(fc.double({ min: -1e3, max: 1e3, noNaN: true }), x =>
                Object.is(sinEstimate(-x, degree) + 0,
                    -sinEstimate(x, degree) + 0));
        }
    });

    it('matches the shape of the upstream Horner loop', () => {
        // p(x) = x * (c0 + c1*x^2 + ... + cn*x^{2n}) with c0 = 1, so the
        // derivative at 0 is exactly 1 and p(x)/x -> 1.
        for (const degree of DEGREES) {
            expect(sinEstimate(1e-8, degree) / 1e-8).toBeCloseTo(1, 12);
        }
    });

    it('agrees with the exact-reduction reference on the RR periodicity', () => {
        // sinEstimateRR is periodic with period equal to the double 2*pi, up
        // to the rounding of the shift itself.
        for (const degree of DEGREES) {
            check(fc.double({ min: -50, max: 50, noNaN: true }), x => {
                const a = sinEstimateRR(x, degree);
                const b = sinEstimateRR(x + 10 * GTE_C_TWO_PI, degree);
                return Math.abs(a - b) <= 1e-13;
            });
        }
    });
});
