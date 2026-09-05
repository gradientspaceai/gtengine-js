import { describe, it, expect } from 'vitest';
import {
    chebyshevRatio, chebyshevRatioUsingCosAngle,
    chebyshevRatios, chebyshevRatiosUsingCosAngle
} from '../src/ChebyshevRatio.js';
import { GTE_C_PI } from '../src/Constants.js';
import { check, fc } from './helpers/arbitraries.js';

describe('chebyshevRatio', () => {
    it('computes sin(t*A)/sin(A) for angles in (0,pi)', () => {
        for (const angle of [0.1, 0.5, 1, Math.PI / 2, 2, 3, 3.14]) {
            for (const t of [0, 0.25, 0.5, 0.75, 1]) {
                expect(chebyshevRatio(t, angle))
                    .toBe(Math.sin(t * angle) / Math.sin(angle));
            }
        }
    });

    it('returns t for angle 0 (removable singularity, l\'Hospital)', () => {
        for (const t of [0, 0.25, 0.5, 0.75, 1]) {
            expect(chebyshevRatio(t, 0)).toBe(t);
        }
    });

    it('approaches t continuously as the angle approaches 0', () => {
        for (const t of [0.25, 0.5, 0.75]) {
            expect(Math.abs(chebyshevRatio(t, 1e-8) - t)).toBeLessThan(1e-15);
        }
    });

    it('is exact at the endpoints of t', () => {
        for (const angle of [0.5, 1.5, 3]) {
            expect(chebyshevRatio(0, angle)).toBe(0);
            expect(chebyshevRatio(1, angle)).toBe(1);
        }
    });

    it('throws for angles outside [0,pi)', () => {
        expect(() => chebyshevRatio(0.5, -0.1)).toThrow('Invalid angle.');
        expect(() => chebyshevRatio(0.5, GTE_C_PI)).toThrow('Invalid angle.');
        expect(() => chebyshevRatio(0.5, 4)).toThrow('Invalid angle.');
    });
});

describe('chebyshevRatioUsingCosAngle', () => {
    it('agrees with chebyshevRatio for the angle acos(cosAngle)', () => {
        for (const cosAngle of [-0.99, -0.5, 0, 0.5, 0.99]) {
            const angle = Math.acos(cosAngle);
            for (const t of [0, 0.25, 0.5, 0.75, 1]) {
                expect(chebyshevRatioUsingCosAngle(t, cosAngle))
                    .toBe(Math.sin(t * angle) / Math.sin(angle));
            }
        }
    });

    it('returns t for cosAngle = 1 (angle 0)', () => {
        for (const t of [0, 0.3, 1]) {
            expect(chebyshevRatioUsingCosAngle(t, 1)).toBe(t);
            // Values beyond 1 also take the angle-0 branch, as upstream.
            expect(chebyshevRatioUsingCosAngle(t, 1.5)).toBe(t);
        }
    });

    it('throws for cosAngle <= -1 (angle pi)', () => {
        expect(() => chebyshevRatioUsingCosAngle(0.5, -1)).toThrow('Invalid angle.');
        expect(() => chebyshevRatioUsingCosAngle(0.5, -2)).toThrow('Invalid angle.');
    });
});

describe('chebyshevRatios', () => {
    it('returns the pair {f(1-t,A), f(t,A)}', () => {
        for (const angle of [0.25, 1, 2.5]) {
            for (const t of [0, 0.25, 0.5, 0.9, 1]) {
                const [f0, f1] = chebyshevRatios(t, angle);
                expect(f0).toBe(chebyshevRatio(1 - t, angle));
                expect(f1).toBe(chebyshevRatio(t, angle));
            }
        }
    });

    it('returns {1-t, t} for angle 0', () => {
        expect(chebyshevRatios(0.25, 0)).toEqual([0.75, 0.25]);
        expect(chebyshevRatios(1, 0)).toEqual([0, 1]);
    });

    it('gives slerp-style barycentric weights whose sum is 1/cos(A/2)', () => {
        // f(1/2,A) + f(1/2,A) = 2*sin(A/2)/sin(A) = 1/cos(A/2), which
        // approaches 1 as the angle approaches 0.
        for (const angle of [0.01, 0.1, 1]) {
            const [f0, f1] = chebyshevRatios(0.5, angle);
            expect(f0).toBe(f1);
            expect(Math.abs(f0 + f1 - 1 / Math.cos(angle / 2))).toBeLessThan(1e-14);
        }
    });

    it('throws for angles outside [0,pi)', () => {
        expect(() => chebyshevRatios(0.5, -1)).toThrow('Invalid angle.');
        expect(() => chebyshevRatios(0.5, GTE_C_PI)).toThrow('Invalid angle.');
    });
});

describe('chebyshevRatiosUsingCosAngle', () => {
    it('agrees with chebyshevRatios for the angle acos(cosAngle)', () => {
        for (const cosAngle of [-0.9, 0, 0.9]) {
            const angle = Math.acos(cosAngle);
            for (const t of [0, 0.25, 0.5, 1]) {
                expect(chebyshevRatiosUsingCosAngle(t, cosAngle))
                    .toEqual(chebyshevRatios(t, angle));
            }
        }
    });

    it('returns {1-t, t} for cosAngle >= 1 (angle 0)', () => {
        expect(chebyshevRatiosUsingCosAngle(0.25, 1)).toEqual([0.75, 0.25]);
        expect(chebyshevRatiosUsingCosAngle(0.25, 2)).toEqual([0.75, 0.25]);
    });

    it('throws for cosAngle <= -1 (angle pi)', () => {
        expect(() => chebyshevRatiosUsingCosAngle(0.5, -1)).toThrow('Invalid angle.');
    });
});

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream ChebyshevRatio.h.
// ---------------------------------------------------------------------------

describe('ChebyshevRatio verification', () => {
    // Angles strictly inside (0,pi); the endpoints are error cases.
    const angleIn = fc.double({ min: 1e-6, max: Math.PI - 1e-6, noNaN: true });
    const anyT = fc.double({ min: -2, max: 3, noNaN: true });

    it('returns the pair {f(1-t,A), f(t,A)} with bit-identical divisions', () => {
        // Upstream hoists sin(angle) out of the two divisions but states that
        // the results must still equal the single-ratio function exactly; a
        // port that multiplied by 1/sin(angle) instead would break this.
        check(fc.tuple(anyT, angleIn), ([t, angle]) => {
            const pair = chebyshevRatios(t, angle);
            return pair[0] === chebyshevRatio(1 - t, angle)
                && pair[1] === chebyshevRatio(t, angle);
        });
        check(fc.tuple(anyT, fc.double({ min: -0.999, max: 0.999, noNaN: true })),
            ([t, cosAngle]) => {
                const pair = chebyshevRatiosUsingCosAngle(t, cosAngle);
                return pair[0] === chebyshevRatioUsingCosAngle(1 - t, cosAngle)
                    && pair[1] === chebyshevRatioUsingCosAngle(t, cosAngle);
            });
    });

    it('satisfies the Chebyshev three-term recurrence', () => {
        // f(t+1,A) = 2*cos(A)*f(t,A) - f(t-1,A) is the identity the ratio is
        // named after, and it is independent of the implementation.
        check(fc.tuple(fc.double({ min: -1, max: 1, noNaN: true }),
            fc.double({ min: 1e-2, max: Math.PI - 1e-2, noNaN: true })),
            ([t, angle]) => {
                const next = chebyshevRatio(t + 1, angle);
                const recurrence = 2 * Math.cos(angle) * chebyshevRatio(t, angle)
                    - chebyshevRatio(t - 1, angle);
                // sin(angle) in the denominator is bounded below by sin(0.01),
                // so the relative conditioning is about 100.
                return Math.abs(next - recurrence)
                    <= 1e-12 * (1 + Math.abs(next));
            });
    });

    it('is exact at t = 0 and t = 1 for every angle', () => {
        check(angleIn, angle => {
            return chebyshevRatio(0, angle) === 0
                && chebyshevRatio(1, angle) === 1
                && chebyshevRatios(0, angle)[1] === 0
                && chebyshevRatios(1, angle)[0] === 0
                && chebyshevRatios(0, angle)[0] === 1
                && chebyshevRatios(1, angle)[1] === 1;
        });
    });

    it('agrees between the angle and cosine entry points', () => {
        // acos(cos(A)) reproduces A to within the conditioning of acos near
        // the endpoints, so the two forms agree to a relative tolerance.
        check(fc.tuple(anyT, fc.double({ min: 0.05, max: Math.PI - 0.05, noNaN: true })),
            ([t, angle]) => {
                const a = chebyshevRatio(t, angle);
                const b = chebyshevRatioUsingCosAngle(t, Math.cos(angle));
                return Math.abs(a - b) <= 1e-9 * (1 + Math.abs(a));
            });
    });

    it('interpolates the endpoints of the arc for t in [0,1]', () => {
        // Both weights are non-negative there, and the sum-to-product form
        //   sin((1-t)A) + sin(tA) = 2*sin(A/2)*cos((t - 1/2)A)
        // gives f(1-t)+f(t) = cos((t - 1/2)A)/cos(A/2), which is at least 1
        // and peaks at 1/cos(A/2) in the middle of the arc.
        check(fc.tuple(fc.double({ min: 0, max: 1, noNaN: true }), angleIn),
            ([t, angle]) => {
                const [f0, f1] = chebyshevRatios(t, angle);
                const sum = f0 + f1;
                const expected = Math.cos((t - 0.5) * angle)
                    / Math.cos(angle / 2);
                return f0 >= 0 && f1 >= 0
                    && sum >= 1 - 1e-12
                    && Math.abs(sum - expected) <= 1e-9 * (1 + sum);
            });
    });

    it('degenerates to linear weights as the angle vanishes', () => {
        check(fc.double({ min: 0, max: 1, noNaN: true }), t => {
            expect(chebyshevRatio(t, 0)).toBe(t);
            expect(chebyshevRatios(t, 0)).toEqual([1 - t, t]);
            expect(chebyshevRatioUsingCosAngle(t, 1)).toBe(t);
            expect(chebyshevRatiosUsingCosAngle(t, 1)).toEqual([1 - t, t]);
            // Continuity: a tiny angle is within rounding of the limit.
            return Math.abs(chebyshevRatio(t, 1e-8) - t) <= 1e-12;
        });
    });

    it('rejects every angle outside [0,pi), NaN included', () => {
        // Upstream falls through both branches and calls LogError; NaN fails
        // the > 0 and == 0 tests, so it lands there too.
        for (const bad of [-1e-8, -1, Math.PI, 4, Number.NaN,
            Number.POSITIVE_INFINITY]) {
            expect(() => chebyshevRatio(0.5, bad)).toThrow('Invalid angle.');
            expect(() => chebyshevRatios(0.5, bad)).toThrow('Invalid angle.');
        }
        for (const bad of [-1, -1.5, Number.NEGATIVE_INFINITY]) {
            expect(() => chebyshevRatioUsingCosAngle(0.5, bad))
                .toThrow('Invalid angle.');
            expect(() => chebyshevRatiosUsingCosAngle(0.5, bad))
                .toThrow('Invalid angle.');
        }
    });

    it('treats cosAngle >= 1 and NaN as the zero angle, as upstream does', () => {
        // The cosine entry points test `cosAngle < 1` first, so anything not
        // less than one - including out-of-range values and NaN - takes the
        // l'Hospital branch instead of raising. Preserved from upstream.
        for (const c of [1, 1.5, Number.NaN, Number.POSITIVE_INFINITY]) {
            expect(chebyshevRatioUsingCosAngle(0.25, c)).toBe(0.25);
            expect(chebyshevRatiosUsingCosAngle(0.25, c)).toEqual([0.75, 0.25]);
        }
    });
});
