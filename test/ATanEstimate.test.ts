import { describe, it, expect } from 'vitest';
import {
    atanEstimate, atanEstimateRR, getATanEstimateMaxError
} from '../src/ATanEstimate.js';
import { GTE_C_HALF_PI } from '../src/Constants.js';
import { tanEstimate } from '../src/TanEstimate.js';
import { check, fc } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream ATanEstimate.h.
// ---------------------------------------------------------------------------

describe('ATanEstimate verification', () => {
    const inDomain = fc.double({ min: -1, max: 1, noNaN: true });
    const anyMagnitude = fc.tuple(
        fc.double({ min: 1, max: 2, noNaN: true, noDefaultInfinity: true }),
        fc.integer({ min: -300, max: 300 }),
        fc.boolean()
    ).map(([m, e, neg]) => (neg ? -m : m) * Math.pow(2, e));

    it('is within the documented bound on fast-check samples of [-1,1]', () => {
        for (const degree of DEGREES) {
            check(inDomain, x =>
                Math.abs(atanEstimate(x, degree) - Math.atan(x))
                    <= MAX_ERROR[degree]);
        }
    });

    it('has a max error over [-1,1] that decreases with the degree', () => {
        // Measured independently of the published table.
        const samples = 40000;
        let previous = Number.POSITIVE_INFINITY;
        for (const degree of DEGREES) {
            let observed = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = -1 + (2 * i) / samples;
                observed = Math.max(observed,
                    Math.abs(atanEstimate(x, degree) - Math.atan(x)));
            }
            expect(observed).toBeLessThanOrEqual(MAX_ERROR[degree]);
            expect(observed).toBeGreaterThan(0.9 * MAX_ERROR[degree]);
            expect(observed).toBeLessThan(previous);
            previous = observed;
        }
    }, 30000);

    it('is exactly odd, for every input and degree', () => {
        // Only odd powers appear and the polynomial sees x*x, so the sign
        // symmetry is bit exact.
        for (const degree of DEGREES) {
            check(fc.double({ min: -1e3, max: 1e3, noNaN: true }), x =>
                Object.is(atanEstimate(-x, degree) + 0,
                    -atanEstimate(x, degree) + 0));
        }
    });

    it('honours the interpolation constraint p(1) = pi/4', () => {
        // Upstream constrains the minimax fit to pass through atan(1). The
        // published coefficients meet it exactly at degrees 3 and 5 and to
        // 4-8 ulps at degrees 7 and 9, but only to 76 and 560 ulps at
        // degrees 11 and 13 (6.2e-14 at degree 13) - the residual of the
        // constrained solve, six orders of magnitude below the degree's own
        // 3.6e-7 max error, so it is a curiosity rather than a defect.
        for (const degree of DEGREES) {
            expect(Math.abs(atanEstimate(1, degree) - Math.PI / 4))
                .toBeLessThanOrEqual(1e-13);
            expect(Math.abs(atanEstimate(-1, degree) + Math.PI / 4))
                .toBeLessThanOrEqual(1e-13);
        }
    });

    it('applies the reciprocal identity outside [-1,1] exactly', () => {
        // Upstream: atan(x) = +-pi/2 - atan(1/x) for |x| > 1. The branch on
        // the sign is what the estimate would get wrong if transcribed with
        // a single sign.
        for (const degree of DEGREES) {
            check(anyMagnitude, x => {
                const expected = Math.abs(x) <= 1
                    ? atanEstimate(x, degree)
                    : (x > 1 ? GTE_C_HALF_PI - atanEstimate(1 / x, degree)
                        : -GTE_C_HALF_PI - atanEstimate(1 / x, degree));
                return atanEstimateRR(x, degree) + 0 === expected + 0;
            });
        }
    });

    it('stays within the max error over the whole real line', () => {
        // The identity transfers the bound from [-1,1] to every x, up to the
        // rounding of the pi/2 constant and of 1/x.
        for (const degree of DEGREES) {
            check(anyMagnitude, x =>
                Math.abs(atanEstimateRR(x, degree) - Math.atan(x))
                    <= MAX_ERROR[degree] + 1e-15);
        }
    });

    it('is exactly odd in the range-reduced form as well', () => {
        for (const degree of DEGREES) {
            check(anyMagnitude, x =>
                Object.is(atanEstimateRR(-x, degree) + 0,
                    -atanEstimateRR(x, degree) + 0));
        }
    });

    it('is strictly increasing over the whole real line', () => {
        // atan is increasing everywhere; the reciprocal branch must not
        // introduce a step at |x| = 1.
        const samples = 20000;
        for (const degree of DEGREES) {
            let previous = Number.NEGATIVE_INFINITY;
            for (let i = 0; i <= samples; ++i) {
                // A tangent-spaced sweep covers both branches densely.
                const u = -1.5 + (3 * i) / samples;
                const x = Math.tan(u);
                const v = atanEstimateRR(x, degree);
                expect(v).toBeGreaterThan(previous);
                previous = v;
            }
        }
    }, 30000);

    it('inverts tanEstimate to the combined accuracy on [-pi/4,pi/4]', () => {
        // atan(tan(u)) = u, cross-checked against another file of this group.
        for (const degree of [11, 13]) {
            check(fc.double({ min: -0.78, max: 0.78, noNaN: true }), u => {
                const t = tanEstimate(u, 13);
                const back = atanEstimateRR(t, degree);
                // d(atan)/dx = 1/(1+t^2) <= 1, so the tan error transfers
                // with a factor of at most one.
                return Math.abs(back - u)
                    <= MAX_ERROR[degree] + 1.06e-8 + 1e-13;
            });
        }
    });
});
