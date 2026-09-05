import { describe, it, expect } from 'vitest';
import {
    cosEstimate, cosEstimateRR, getCosEstimateMaxError
} from '../src/CosEstimate.js';
import { GTE_C_HALF_PI, GTE_C_PI, GTE_C_TWO_PI } from '../src/Constants.js';
import { sinEstimateRR } from '../src/SinEstimate.js';
import { check, fc } from './helpers/arbitraries.js';
import { exactRemainder } from './helpers/exact.js';

// Upstream-documented maximum errors of SinEstimate, used by the
// cofunction cross-check below.
const SIN_MAX_ERROR: Record<number, number> = {
    5: 1.4001209384651e-4,
    9: 5.2010783457846e-9,
    11: 1.9323431743601e-11
};

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

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream CosEstimate.h.
// ---------------------------------------------------------------------------

describe('CosEstimate verification', () => {
    // Covers every binade rather than a bounded interval: the range reduction
    // is where the port stands in for C++ std::remainder, and the two only
    // part company once the quotient x/(2*pi) needs more than 53 bits.
    const anyMagnitude = fc.tuple(
        fc.double({ min: 1, max: 2, noNaN: true, noDefaultInfinity: true }),
        fc.integer({ min: -40, max: 62 }),
        fc.boolean()
    ).map(([m, e, neg]) => (neg ? -m : m) * Math.pow(2, e));

    const reduce = (x: number): number => {
        const r = exactRemainder(x, GTE_C_TWO_PI);
        if (r > GTE_C_HALF_PI) { return GTE_C_PI - r; }
        if (r < -GTE_C_HALF_PI) { return -GTE_C_PI - r; }
        return r;
    };
    const signOf = (x: number): number => {
        const r = exactRemainder(x, GTE_C_TWO_PI);
        return (r > GTE_C_HALF_PI || r < -GTE_C_HALF_PI) ? -1 : 1;
    };

    it('reduces the argument exactly, as std::remainder does', () => {
        // cosEstimateRR(x) must equal +-cosEstimate of the reduction of x by
        // the *double* 2*pi computed in infinite precision. The rounded
        // quotient x - Math.round(x/y)*y is a different function: it is off
        // by 1e-6 near |x| = 1e10 and leaves [-pi,pi] near |x| = 1e14.
        for (const degree of DEGREES) {
            check(anyMagnitude, x => {
                const expected = signOf(x) * cosEstimate(reduce(x), degree);
                return cosEstimateRR(x, degree) + 0 === expected + 0;
            });
        }
        for (const x of [1e10, 1e12, 1e14, 1e16, -9769346117865922,
            -26789498295071.3, 123456789012345, 4.5e15]) {
            for (const degree of DEGREES) {
                expect(cosEstimateRR(x, degree) + 0)
                    .toBe(signOf(x) * cosEstimate(reduce(x), degree) + 0);
            }
        }
    });

    it('never leaves the range of a cosine, whatever the magnitude of x', () => {
        for (const degree of DEGREES) {
            check(anyMagnitude, x =>
                Math.abs(cosEstimateRR(x, degree)) <= 1 + MAX_ERROR[degree]);
        }
    });

    it('is within the documented bound on fast-check samples of [-pi/2,pi/2]', () => {
        const inDomain = fc.double(
            { min: -GTE_C_HALF_PI, max: GTE_C_HALF_PI, noNaN: true });
        for (const degree of DEGREES) {
            check(inDomain, x =>
                Math.abs(cosEstimate(x, degree) - Math.cos(x))
                    <= MAX_ERROR[degree]);
        }
    });

    it('has a max error over the domain that decreases with the degree', () => {
        // Measured independently of the published table (the RotationEstimate
        // tables were found to understate the true error).
        const samples = 20000;
        let previous = Number.POSITIVE_INFINITY;
        for (const degree of DEGREES) {
            let observed = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = -GTE_C_HALF_PI + (2 * GTE_C_HALF_PI * i) / samples;
                observed = Math.max(observed,
                    Math.abs(cosEstimate(x, degree) - Math.cos(x)));
            }
            expect(observed).toBeLessThanOrEqual(MAX_ERROR[degree]);
            expect(observed).toBeGreaterThan(0.9 * MAX_ERROR[degree]);
            expect(observed).toBeLessThan(previous);
            previous = observed;
        }
    }, 30000);

    it('is exactly even, for every input and degree', () => {
        // Only even powers of x appear, and the polynomial sees x*x, so the
        // equality is bit exact rather than approximate.
        for (const degree of DEGREES) {
            check(fc.double({ min: -1e3, max: 1e3, noNaN: true }), x =>
                cosEstimate(-x, degree) === cosEstimate(x, degree));
        }
    });

    it('is exactly 1 at x = 0 (the constant term is 1)', () => {
        for (const degree of DEGREES) {
            expect(cosEstimate(0, degree)).toBe(1);
        }
    });

    it('tracks sinEstimateRR through the cofunction identity', () => {
        // cos(x) = sin(x + pi/2). Both estimates carry their own minimax
        // error, so the two bounds add.
        for (const [cd, sd] of [[4, 5], [8, 9], [10, 11]] as const) {
            check(fc.double({ min: -30, max: 30, noNaN: true }), x => {
                const a = cosEstimateRR(x, cd);
                const b = sinEstimateRR(x + GTE_C_HALF_PI, sd);
                return Math.abs(a - b)
                    <= MAX_ERROR[cd] + SIN_MAX_ERROR[sd] + 1e-13;
            });
        }
    });
});
