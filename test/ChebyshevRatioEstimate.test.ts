import { describe, it, expect } from 'vitest';
import {
    chebyshevRatioEstimate, getChebyshevRatioEstimateMaxError,
    chebyshevRatioEstimateR, getChebyshevRatioEstimateRMaxError
} from '../src/ChebyshevRatioEstimate';

// Upstream-documented maximum errors for degrees 1..16 (angle in [0,pi/2]).
const MAX_ERROR = [
    1.8249897492955e-2,
    5.2760601519022e-3,
    1.8055057987877e-3,
    6.7244299646175e-4,
    2.6386437427495e-4,
    1.0731422197408e-4,
    4.4805894183764e-5,
    1.9088088593749e-5,
    8.2629028074211e-6,
    3.6237273527418e-6,
    1.6064797200289e-6,
    7.1872518425665e-7,
    3.2407757655229e-7,
    1.4712279927665e-7,
    6.7187475472075e-8,
    3.0844086507110e-8
];

// Upstream-documented maximum errors for degrees 1..12 (angle in [0,pi/4]).
const MAX_ERROR_R = [
    8.6832275204274e-4,
    6.6040175097815e-5,
    6.1949661303018e-6,
    6.4578503422564e-7,
    7.1792162659179e-8,
    8.3364721792379e-9,
    9.9903230132981e-10,
    1.2262002524466e-10,
    1.5335510639148e-11,
    1.9472201628901e-12,
    2.5046631435544e-13,
    3.2751579226443e-14
];

// Exact Chebyshev ratio f(t,A) = sin(t*A)/sin(A) with x = cos(A); f(t,0) = t.
function exactRatio(t: number, x: number): number {
    if (x >= 1) {
        return t;
    }
    const angle = Math.acos(x);
    return Math.sin(t * angle) / Math.sin(angle);
}

describe('chebyshevRatioEstimate', () => {
    it('stays within the documented max error for A in [0,pi/2]', () => {
        const tSamples = 64;
        const aSamples = 64;
        for (let degree = 1; degree <= 16; ++degree) {
            const bound = MAX_ERROR[degree - 1];
            let maxObserved = 0;
            for (let i = 0; i <= tSamples; ++i) {
                const t = i / tSamples;
                for (let j = 0; j <= aSamples; ++j) {
                    const angle = (Math.PI / 2) * (j / aSamples);
                    const x = Math.cos(angle);
                    const [f0, f1] = chebyshevRatioEstimate(t, x, degree);
                    const e0 = Math.abs(f0 - exactRatio(1 - t, x));
                    const e1 = Math.abs(f1 - exactRatio(t, x));
                    maxObserved = Math.max(maxObserved, e0, e1);
                }
            }
            expect(maxObserved).toBeLessThanOrEqual(bound * (1 + 1e-6) + 1e-14);
        }
    });

    it('approaches the bound on a dense grid for low degrees (tight minimax)', () => {
        const tSamples = 128;
        const aSamples = 128;
        for (let degree = 1; degree <= 4; ++degree) {
            const bound = MAX_ERROR[degree - 1];
            let maxObserved = 0;
            for (let i = 0; i <= tSamples; ++i) {
                const t = i / tSamples;
                for (let j = 0; j <= aSamples; ++j) {
                    const x = Math.cos((Math.PI / 2) * (j / aSamples));
                    const [, f1] = chebyshevRatioEstimate(t, x, degree);
                    maxObserved = Math.max(maxObserved, Math.abs(f1 - exactRatio(t, x)));
                }
            }
            expect(maxObserved).toBeGreaterThan(0.25 * bound);
        }
    });

    it('is exact for x = 1 (angle 0): returns {1-t, t}', () => {
        for (const degree of [1, 4, 8, 16]) {
            for (const t of [0, 0.25, 0.5, 0.75, 1]) {
                expect(chebyshevRatioEstimate(t, 1, degree)).toEqual([1 - t, t]);
            }
        }
    });

    it('is exact for t = 0 and t = 1 in the corresponding component', () => {
        // For t = 0 the f1 recursion starts at 0 and every term retains the
        // zero factor; symmetrically for t = 1 and f0.
        for (const degree of [1, 8, 16]) {
            for (const x of [0, 0.5, 1]) {
                expect(chebyshevRatioEstimate(0, x, degree)[1]).toBe(0);
                expect(chebyshevRatioEstimate(1, x, degree)[0]).toBe(0);
            }
        }
    });

    it('has the symmetry f(t) <-> f(1-t) between the two components', () => {
        for (const degree of [3, 10]) {
            for (const x of [0.2, 0.7]) {
                for (const t of [0.125, 0.375]) {
                    const f = chebyshevRatioEstimate(t, x, degree);
                    const g = chebyshevRatioEstimate(1 - t, x, degree);
                    expect(f[0]).toBe(g[1]);
                    expect(f[1]).toBe(g[0]);
                }
            }
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => chebyshevRatioEstimate(0.5, 0.5, 0)).toThrow('Invalid degree.');
        expect(() => chebyshevRatioEstimate(0.5, 0.5, 17)).toThrow('Invalid degree.');
        expect(() => chebyshevRatioEstimate(0.5, 0.5, 1.5)).toThrow('Invalid degree.');
    });
});

describe('chebyshevRatioEstimateR', () => {
    it('stays within the documented max error for A in [0,pi/4]', () => {
        const tSamples = 64;
        const aSamples = 64;
        for (let degree = 1; degree <= 12; ++degree) {
            const bound = MAX_ERROR_R[degree - 1];
            let maxObserved = 0;
            for (let i = 0; i <= tSamples; ++i) {
                const t = i / tSamples;
                for (let j = 0; j <= aSamples; ++j) {
                    const angle = (Math.PI / 4) * (j / aSamples);
                    const x = Math.cos(angle);
                    const [f0, f1] = chebyshevRatioEstimateR(t, x, degree);
                    const e0 = Math.abs(f0 - exactRatio(1 - t, x));
                    const e1 = Math.abs(f1 - exactRatio(t, x));
                    maxObserved = Math.max(maxObserved, e0, e1);
                }
            }
            // The additive cushion absorbs roundoff of the evaluation and of
            // the exact reference, which matters for the smallest bounds.
            expect(maxObserved).toBeLessThanOrEqual(bound * (1 + 1e-6) + 5e-15);
        }
    });

    it('is more accurate than the full-range estimate at the same degree', () => {
        for (let degree = 1; degree <= 12; ++degree) {
            expect(getChebyshevRatioEstimateRMaxError(degree))
                .toBeLessThan(getChebyshevRatioEstimateMaxError(degree));
        }
    });

    it('is exact for x = 1 (angle 0): returns {1-t, t}', () => {
        for (const degree of [1, 6, 12]) {
            for (const t of [0, 0.5, 1]) {
                expect(chebyshevRatioEstimateR(t, 1, degree)).toEqual([1 - t, t]);
            }
        }
    });

    it('throws for invalid degrees', () => {
        expect(() => chebyshevRatioEstimateR(0.5, 0.9, 0)).toThrow('Invalid degree.');
        expect(() => chebyshevRatioEstimateR(0.5, 0.9, 13)).toThrow('Invalid degree.');
    });
});

describe('max error queries', () => {
    it('return the documented bounds and decrease with degree', () => {
        for (let degree = 1; degree <= 16; ++degree) {
            expect(getChebyshevRatioEstimateMaxError(degree)).toBe(MAX_ERROR[degree - 1]);
        }
        for (let degree = 1; degree <= 12; ++degree) {
            expect(getChebyshevRatioEstimateRMaxError(degree)).toBe(MAX_ERROR_R[degree - 1]);
        }
        for (let degree = 2; degree <= 16; ++degree) {
            expect(getChebyshevRatioEstimateMaxError(degree))
                .toBeLessThan(getChebyshevRatioEstimateMaxError(degree - 1));
        }
        for (let degree = 2; degree <= 12; ++degree) {
            expect(getChebyshevRatioEstimateRMaxError(degree))
                .toBeLessThan(getChebyshevRatioEstimateRMaxError(degree - 1));
        }
    });

    it('throw for invalid degrees', () => {
        expect(() => getChebyshevRatioEstimateMaxError(0)).toThrow('Invalid degree.');
        expect(() => getChebyshevRatioEstimateMaxError(17)).toThrow('Invalid degree.');
        expect(() => getChebyshevRatioEstimateRMaxError(0)).toThrow('Invalid degree.');
        expect(() => getChebyshevRatioEstimateRMaxError(13)).toThrow('Invalid degree.');
    });
});
