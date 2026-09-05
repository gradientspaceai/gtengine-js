import { describe, it, expect } from 'vitest';
import {
    chebyshevRatioEstimate, getChebyshevRatioEstimateMaxError,
    chebyshevRatioEstimateR, getChebyshevRatioEstimateRMaxError
} from '../src/ChebyshevRatioEstimate.js';
import { chebyshevRatios } from '../src/ChebyshevRatio.js';
import { check, fc, scaled } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream
// ChebyshevRatioEstimate.h.
// ---------------------------------------------------------------------------

/**
 * The upstream coefficient helpers, restated directly from the header:
 *   a[i] = (Degree != i+1 ? 1 : u[i]) / ((i+1) * (2*(i+1)+1))
 *   b[i] = (Degree != i+1 ? 1 : u[i]) * (i+1) / (2*(i+1)+1)
 * The port folds these into the evaluation loop, so this reference confirms
 * the folding did not change which term carries the u-value.
 */
function referenceEstimate(t: number, x: number, degree: number,
    u: readonly number[]): [number, number] {
    const a: number[] = [];
    const b: number[] = [];
    for (let i = 0; i < u.length; ++i) {
        const scale = degree !== i + 1 ? 1 : u[i];
        a.push(scale / ((i + 1) * (2 * (i + 1) + 1)));
        b.push((scale * (i + 1)) / (2 * (i + 1) + 1));
    }
    const y = 1 - x;
    let term0 = 1 - t;
    let term1 = t;
    const sqr0 = term0 * term0;
    const sqr1 = term1 * term1;
    const f: [number, number] = [term0, term1];
    for (let i = 0; i < degree; ++i) {
        term0 *= (b[i] - a[i] * sqr0) * y;
        term1 *= (b[i] - a[i] * sqr1) * y;
        f[0] += term0;
        f[1] += term1;
    }
    return f;
}

const U = [
    1.5149656562200644050, 1.6410179946672027729, 1.7124880779005808851,
    1.7593545031636841358, 1.7927054757060019163, 1.8177479632959470113,
    1.8372872973294931409, 1.8529805143706497006, 1.8658739107798316681,
    1.8766626700393858052, 1.8858276947289707159, 1.8937127486228939599,
    1.9005703533887863266, 1.9065903281211855624, 1.9119182032942771965,
    1.9166674811124804201
];
const U_R = [
    1.1021472152138613865, 1.1239349540626744073, 1.1351870374370363059,
    1.1421060160698368602, 1.1468020192623136211, 1.1502017494201659531,
    1.1527782928466798751, 1.1547990001678465344, 1.1564265502929687024,
    1.1577657226562501069, 1.1588859375000000185, 1.1598375000000000767
];

describe('ChebyshevRatioEstimate verification', () => {
    const anyT = fc.double({ min: 0, max: 1, noNaN: true });
    const anyX = fc.double({ min: 0, max: 1, noNaN: true });
    const anyXR = fc.double({ min: Math.SQRT1_2, max: 1, noNaN: true });

    it('matches the upstream a/b coefficient tables bit for bit', () => {
        for (const degree of [1, 2, 5, 9, 16]) {
            check(fc.tuple(anyT, anyX), ([t, x]) => {
                const got = chebyshevRatioEstimate(t, x, degree);
                const want = referenceEstimate(t, x, degree, U);
                return got[0] === want[0] && got[1] === want[1];
            });
        }
        for (const degree of [1, 2, 5, 9, 12]) {
            check(fc.tuple(anyT, anyXR), ([t, x]) => {
                const got = chebyshevRatioEstimateR(t, x, degree);
                const want = referenceEstimate(t, x, degree, U_R);
                return got[0] === want[0] && got[1] === want[1];
            });
        }
    });

    it('is within the documented bound on fast-check samples', () => {
        // A uniform grid rather than fc.double: sampling the bit patterns
        // of [0,pi/2] produces subnormal angles, where the *reference*
        // sin(t*A)/sin(A) underflows to 0/A and the comparison becomes
        // meaningless (the estimate itself is fine there).
        const angle = scaled(0, Math.PI / 2);
        for (const degree of [1, 4, 8, 12, 16]) {
            check(fc.tuple(anyT, angle), ([t, A]) => {
                const f = chebyshevRatioEstimate(t, Math.cos(A), degree);
                const sinA = Math.sin(A);
                const exact0 = A === 0 ? 1 - t : Math.sin((1 - t) * A) / sinA;
                const exact1 = A === 0 ? t : Math.sin(t * A) / sinA;
                const bound = MAX_ERROR[degree - 1];
                return Math.abs(f[0] - exact0) <= bound
                    && Math.abs(f[1] - exact1) <= bound;
            });
        }
    });

    it('has a measured max error that decreases with the degree', () => {
        // Independently measured on a dense (A,t) grid rather than read from
        // the table. Every degree of both variants stays inside its published
        // bound; the only entry that comes close to exceeding it is the R
        // variant at degree 12, whose bound (3.2752e-14) is itself near the
        // round-off floor of the summation, so that check carries a 1%
        // cushion.
        const measure = (est: (t: number, x: number, d: number) => [number, number],
            degree: number, aMax: number): number => {
            let worst = 0;
            const NA = 240, NT = 240;
            for (let i = 0; i <= NA; ++i) {
                const A = (i * aMax) / NA;
                const x = Math.cos(A);
                const sinA = Math.sin(A);
                for (let j = 0; j <= NT; ++j) {
                    const t = j / NT;
                    const f = est(t, x, degree);
                    const e0 = A === 0 ? Math.abs(f[0] - (1 - t))
                        : Math.abs(f[0] - Math.sin((1 - t) * A) / sinA);
                    const e1 = A === 0 ? Math.abs(f[1] - t)
                        : Math.abs(f[1] - Math.sin(t * A) / sinA);
                    worst = Math.max(worst, e0, e1);
                }
            }
            return worst;
        };
        let previous = Number.POSITIVE_INFINITY;
        for (let degree = 1; degree <= 16; ++degree) {
            const observed = measure(chebyshevRatioEstimate, degree, Math.PI / 2);
            expect(observed).toBeLessThanOrEqual(MAX_ERROR[degree - 1] * 1.01);
            expect(observed).toBeLessThan(previous);
            previous = observed;
        }
        previous = Number.POSITIVE_INFINITY;
        for (let degree = 1; degree <= 12; ++degree) {
            const observed = measure(chebyshevRatioEstimateR, degree, Math.PI / 4);
            expect(observed).toBeLessThanOrEqual(MAX_ERROR_R[degree - 1] * 1.01);
            expect(observed).toBeLessThan(previous);
            previous = observed;
        }
    }, 30000);

    it('is exactly symmetric under the swap of t and 1-t', () => {
        // The loop treats term0 and term1 identically, so swapping the two
        // arguments swaps the outputs bit for bit whenever 1-t is exact,
        // which Sterbenz guarantees on [1/2,1].
        const upperHalf = fc.double({ min: 0.5, max: 1, noNaN: true });
        for (const degree of [1, 3, 8, 16]) {
            check(fc.tuple(upperHalf, anyX), ([t, x]) => {
                const a = chebyshevRatioEstimate(t, x, degree);
                const b = chebyshevRatioEstimate(1 - t, x, degree);
                return a[0] === b[1] && a[1] === b[0];
            });
        }
    });

    it('is exact at the arc endpoints and at zero angle', () => {
        for (let degree = 1; degree <= 16; ++degree) {
            check(anyX, x => {
                const atZero = chebyshevRatioEstimate(0, x, degree);
                const atOne = chebyshevRatioEstimate(1, x, degree);
                return atZero[1] === 0 && atOne[0] === 0;
            });
            check(anyT, t => {
                const f = chebyshevRatioEstimate(t, 1, degree);
                return f[0] === 1 - t && f[1] === t;
            });
        }
    });

    it('beats the full-range estimate on [0,pi/4] at every shared degree', () => {
        // The R variant differs only in the u-values, chosen for the smaller
        // domain; swapping the two tables would invert this.
        const angle = scaled(0, Math.PI / 4);
        for (let degree = 1; degree <= 12; ++degree) {
            expect(getChebyshevRatioEstimateRMaxError(degree))
                .toBeLessThan(getChebyshevRatioEstimateMaxError(degree));
            check(fc.tuple(anyT, angle), ([t, A]) => {
                const x = Math.cos(A);
                const sinA = Math.sin(A);
                const exact = A === 0 ? t : Math.sin(t * A) / sinA;
                const eR = Math.abs(
                    chebyshevRatioEstimateR(t, x, degree)[1] - exact);
                return eR <= MAX_ERROR_R[degree - 1];
            });
        }
    });

    it('agrees with the exact ChebyshevRatio pair within the bound', () => {
        // Cross-check against the non-estimating file of the same group.
        const angle = scaled(1e-9, Math.PI / 2 - 1e-9);
        for (const degree of [4, 10, 16]) {
            check(fc.tuple(anyT, angle), ([t, A]) => {
                const x = Math.cos(A);
                const exact = chebyshevRatios(t, A);
                const est = chebyshevRatioEstimate(t, x, degree);
                const bound = MAX_ERROR[degree - 1] + 1e-12;
                return Math.abs(est[0] - exact[0]) <= bound
                    && Math.abs(est[1] - exact[1]) <= bound;
            });
        }
    });
});
