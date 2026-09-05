import { describe, it, expect } from 'vitest';
import {
    tanEstimate, tanEstimateRR, getTanEstimateMaxError
} from '../src/TanEstimate.js';
import { GTE_C_HALF_PI, GTE_C_PI, GTE_C_QUARTER_PI } from '../src/Constants.js';
import { sinEstimate } from '../src/SinEstimate.js';
import { cosEstimate } from '../src/CosEstimate.js';
import { check, fc } from './helpers/arbitraries.js';
import { exactRemainder } from './helpers/exact.js';

// Upstream-documented maximum errors of SinEstimate and CosEstimate,
// used by the tan = sin/cos cross-check below.
const SIN_MAX_ERROR: Record<number, number> = {
    9: 5.2010783457846e-9,
    11: 1.9323431743601e-11
};
const COS_MAX_ERROR: Record<number, number> = {
    10: 2.7008567604626e-10
};

const DEGREES = [3, 5, 7, 9, 11, 13] as const;

// Upstream-documented maximum errors of the estimates on [-pi/4,pi/4].
const MAX_ERROR: Record<number, number> = {
    3: 1.1661892256205e-2,
    5: 5.8431854390146e-4,
    7: 3.5418688397793e-5,
    9: 2.2988173248307e-6,
    11: 1.5426258070939e-7,
    13: 1.0550265105991e-8
};

describe('tanEstimate', () => {
    it('stays within the documented max error on a dense grid of [-pi/4,pi/4]', () => {
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree];
            let maxObserved = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = -GTE_C_QUARTER_PI + (2 * GTE_C_QUARTER_PI * i) / samples;
                const error = Math.abs(tanEstimate(x, degree) - Math.tan(x));
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
            expect(tanEstimate(0, degree)).toBe(0);
        }
    });

    it('is an odd function (exact sign symmetry)', () => {
        // Only odd powers appear: p(x) = x*q(x^2), so p(-x) = -p(x) exactly.
        for (const degree of DEGREES) {
            for (const x of [0.1, 0.3, 0.5, GTE_C_QUARTER_PI]) {
                expect(tanEstimate(-x, degree)).toBe(-tanEstimate(x, degree));
            }
        }
    });

    it('satisfies the constraint p(pi/4) = 1 up to roundoff', () => {
        for (const degree of DEGREES) {
            expect(Math.abs(tanEstimate(GTE_C_QUARTER_PI, degree) - 1))
                .toBeLessThanOrEqual(1e-14);
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => tanEstimate(0.5, 1)).toThrow('Invalid degree.');
        expect(() => tanEstimate(0.5, 4)).toThrow('Invalid degree.');
        expect(() => tanEstimate(0.5, 15)).toThrow('Invalid degree.');
        expect(() => tanEstimate(0.5, 3.5)).toThrow('Invalid degree.');
    });
});

describe('tanEstimateRR', () => {
    it('agrees with tanEstimate on [-pi/4,pi/4]', () => {
        // remainder(x, pi) = x exactly there, so the same branch is taken.
        for (const degree of DEGREES) {
            for (const x of [-GTE_C_QUARTER_PI, -0.5, 0, 0.5, GTE_C_QUARTER_PI]) {
                expect(tanEstimateRR(x, degree)).toBe(tanEstimate(x, degree));
            }
        }
    });

    it('stays within a derivative-scaled error bound over a wide range', () => {
        // For reduced |y| <= pi/4 the polynomial bound applies directly.
        // For |y| in (pi/4,pi/2) the identity tan(y) = (1+p)/(1-p) with
        // p = tan(y - pi/4) + e amplifies the polynomial error by
        // 2/(1-p)^2 = (1+tan(y)^2)/(1+p^2) <= 1 + tan(y)^2 = sec(y)^2, so
        // |est - tan(x)| <= bound * sec(x)^2 up to roundoff. Points too
        // close to the poles are skipped.
        const samples = 4096;
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree];
            for (let i = 0; i <= samples; ++i) {
                const x = -10 + (20 * i) / samples;
                if (Math.abs(Math.cos(x)) < 0.01) {
                    continue;
                }
                const tanx = Math.tan(x);
                const tolerance = bound * (1 + tanx * tanx) * (1 + 1e-6) + 1e-12;
                expect(Math.abs(tanEstimateRR(x, degree) - tanx))
                    .toBeLessThanOrEqual(tolerance);
            }
        }
    });

    it('uses the shift identities accurately at sample angles beyond pi/4', () => {
        // tan(pi/3) = sqrt(3), tan(-pi/3) = -sqrt(3).
        for (const degree of DEGREES) {
            const bound = MAX_ERROR[degree];
            expect(Math.abs(tanEstimateRR(Math.PI / 3, degree) - Math.sqrt(3)))
                .toBeLessThanOrEqual(4 * bound * (1 + 1e-6) + 1e-13);
            expect(Math.abs(tanEstimateRR(-Math.PI / 3, degree) + Math.sqrt(3)))
                .toBeLessThanOrEqual(4 * bound * (1 + 1e-6) + 1e-13);
        }
    });

    it('respects periodicity: matches values shifted by pi within tolerance', () => {
        for (const degree of DEGREES) {
            for (const x of [0.3, 0.7, -0.5]) {
                const base = tanEstimateRR(x, degree);
                const shifted = tanEstimateRR(x + Math.PI, degree);
                expect(Math.abs(shifted - base)).toBeLessThanOrEqual(1e-12);
            }
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => tanEstimateRR(1, 2)).toThrow('Invalid degree.');
        expect(() => tanEstimateRR(1, 15)).toThrow('Invalid degree.');
    });
});

describe('getTanEstimateMaxError', () => {
    it('returns the documented bounds', () => {
        for (const degree of DEGREES) {
            expect(getTanEstimateMaxError(degree)).toBe(MAX_ERROR[degree]);
        }
    });

    it('is strictly decreasing in the degree', () => {
        for (let i = 1; i < DEGREES.length; ++i) {
            expect(getTanEstimateMaxError(DEGREES[i]))
                .toBeLessThan(getTanEstimateMaxError(DEGREES[i - 1]));
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => getTanEstimateMaxError(1)).toThrow('Invalid degree.');
        expect(() => getTanEstimateMaxError(15)).toThrow('Invalid degree.');
    });
});

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream TanEstimate.h.
// ---------------------------------------------------------------------------

describe('TanEstimate verification', () => {
    const anyMagnitude = fc.tuple(
        fc.double({ min: 1, max: 2, noNaN: true, noDefaultInfinity: true }),
        fc.integer({ min: -40, max: 62 }),
        fc.boolean()
    ).map(([m, e, neg]) => (neg ? -m : m) * Math.pow(2, e));

    // The upstream reduction, with the argument reduction done exactly.
    const referenceRR = (x: number, degree: number): number => {
        const r = exactRemainder(x, GTE_C_PI);
        let y: number;
        if (r > GTE_C_HALF_PI) { y = r - GTE_C_PI; }
        else if (r < -GTE_C_HALF_PI) { y = r + GTE_C_PI; }
        else { y = r; }
        if (Math.abs(y) <= GTE_C_QUARTER_PI) { return tanEstimate(y, degree); }
        if (y > GTE_C_QUARTER_PI) {
            const poly = tanEstimate(y - GTE_C_QUARTER_PI, degree);
            return (1 + poly) / (1 - poly);
        }
        const poly = tanEstimate(y + GTE_C_QUARTER_PI, degree);
        return (-1 + poly) / (1 + poly);
    };

    it('reduces the argument exactly, as std::remainder does', () => {
        // std::remainder is computed as if in infinite precision; the rounded
        // quotient x - Math.round(x/pi)*pi is not the same function and is
        // already wrong by 1e-6 near |x| = 1e10.
        for (const degree of DEGREES) {
            check(anyMagnitude, x =>
                tanEstimateRR(x, degree) + 0 === referenceRR(x, degree) + 0);
        }
        for (const x of [1e10, 1e12, 1e14, 1e16, -9769346117865922,
            -26789498295071.3, 123456789012345, 4.5e15]) {
            for (const degree of DEGREES) {
                expect(tanEstimateRR(x, degree) + 0)
                    .toBe(referenceRR(x, degree) + 0);
            }
        }
    });

    it('never returns NaN, whatever the magnitude of x', () => {
        // remainder(x, pi) lands in [-pi/2,pi/2]; the identity branches then
        // divide by 1 -+ tan(z) with |z| <= pi/4, which is bounded away from
        // zero. A reduction that overshoots can drive the divisor to zero.
        for (const degree of DEGREES) {
            check(anyMagnitude, x => Number.isFinite(tanEstimateRR(x, degree)));
        }
    });

    it('is within the documented bound on fast-check samples of [-pi/4,pi/4]', () => {
        const inDomain = fc.double(
            { min: -GTE_C_QUARTER_PI, max: GTE_C_QUARTER_PI, noNaN: true });
        for (const degree of DEGREES) {
            check(inDomain, x =>
                Math.abs(tanEstimate(x, degree) - Math.tan(x))
                    <= MAX_ERROR[degree]);
        }
    });

    it('has a max error over the domain that decreases with the degree', () => {
        const samples = 20000;
        let previous = Number.POSITIVE_INFINITY;
        for (const degree of DEGREES) {
            let observed = 0;
            for (let i = 0; i <= samples; ++i) {
                const x = -GTE_C_QUARTER_PI
                    + (2 * GTE_C_QUARTER_PI * i) / samples;
                observed = Math.max(observed,
                    Math.abs(tanEstimate(x, degree) - Math.tan(x)));
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
                Object.is(tanEstimate(-x, degree) + 0,
                    -tanEstimate(x, degree) + 0));
        }
    });

    it('satisfies tan = sin/cos on the reduced domain', () => {
        // Cross-check against two independent estimates from this group. On
        // [-pi/4,pi/4] the cosine is at least sqrt(2)/2, so the quotient's
        // error is bounded by (errSin + |tan| * errCos) / cos.
        const inDomain = fc.double(
            { min: -GTE_C_QUARTER_PI, max: GTE_C_QUARTER_PI, noNaN: true });
        for (const [td, sd, cd] of [[9, 9, 10], [13, 11, 10]] as const) {
            check(inDomain, x => {
                const s = sinEstimate(x, sd);
                const c = cosEstimate(x, cd);
                const bound = (SIN_MAX_ERROR[sd] + Math.abs(s / c)
                    * COS_MAX_ERROR[cd]) / Math.abs(c) + MAX_ERROR[td];
                return Math.abs(tanEstimate(x, td) - s / c) <= bound;
            });
        }
    });

    it('matches the pi-periodicity of the tangent', () => {
        for (const degree of DEGREES) {
            check(fc.double({ min: -1.5, max: 1.5, noNaN: true }), x => {
                if (Math.abs(Math.tan(x)) > 20) { return true; }
                const a = tanEstimateRR(x, degree);
                const b = tanEstimateRR(x + 8 * GTE_C_PI, degree);
                return Math.abs(a - b) <= 1e-11 * (1 + Math.abs(a));
            });
        }
    });
});
