import { describe, it, expect } from 'vitest';
import { Integration } from '../src/Integration.js';
import { check, expectClose, fc, scaled } from './helpers/arbitraries.js';

describe('Integration.trapezoidRule', () => {
    it('is exact for a linear integrand', () => {
        // The trapezoid rule reproduces linear functions exactly.
        const value = Integration.trapezoidRule(2, 0, 1, t => 3 * t + 1);
        expect(value).toBeCloseTo(2.5, 12);
        const refined = Integration.trapezoidRule(17, 0, 1, t => 3 * t + 1);
        expect(refined).toBeCloseTo(2.5, 12);
    });

    it('approximates a quadratic with the expected O(h^2) error', () => {
        const exact = 1 / 3;
        const coarse = Integration.trapezoidRule(5, 0, 1, t => t * t);
        const fine = Integration.trapezoidRule(9, 0, 1, t => t * t);
        expect(coarse).toBeCloseTo(exact, 1);
        expect(fine).toBeCloseTo(exact, 2);
        // Halving h reduces the error by a factor of about four.
        const ratio = Math.abs(coarse - exact) / Math.abs(fine - exact);
        expect(ratio).toBeGreaterThan(3.5);
        expect(ratio).toBeLessThan(4.5);
    });

    it('converges for a transcendental integrand', () => {
        // integral of sin over [0, pi] is 2.
        const value = Integration.trapezoidRule(1001, 0, Math.PI, Math.sin);
        expect(value).toBeCloseTo(2, 5);
    });

    it('changes sign when the endpoints are swapped', () => {
        const forward = Integration.trapezoidRule(101, 0, 1, t => t * t);
        const backward = Integration.trapezoidRule(101, 1, 0, t => t * t);
        expect(backward).toBeCloseTo(-forward, 12);
    });
});

describe('Integration.romberg', () => {
    it('reduces to the two-sample trapezoid rule for order 1', () => {
        const f = (t: number) => Math.exp(t);
        expect(Integration.romberg(1, 0, 1, f))
            .toBeCloseTo(Integration.trapezoidRule(2, 0, 1, f), 14);
    });

    it('is exact for polynomials of low degree', () => {
        // Romberg with order k is exact for polynomials of degree < 2k.
        const cubic = (t: number) => 4 * t * t * t - 3 * t * t + 2 * t - 1;
        // integral over [0, 2] is 16 - 8 + 4 - 2 = 10.
        expect(Integration.romberg(3, 0, 2, cubic)).toBeCloseTo(10, 10);
    });

    it('converges rapidly for smooth integrands', () => {
        expect(Integration.romberg(10, 0, Math.PI, Math.sin))
            .toBeCloseTo(2, 12);
        expect(Integration.romberg(12, 0, 1, Math.exp))
            .toBeCloseTo(Math.E - 1, 12);
        expect(Integration.romberg(12, 0, 1, t => 4 / (1 + t * t)))
            .toBeCloseTo(Math.PI, 12);
    });

    it('beats the trapezoid rule at a comparable sample count', () => {
        const exact = Math.E - 1;
        // Romberg of order 8 uses about 2^7 subintervals.
        const romberg = Integration.romberg(8, 0, 1, Math.exp);
        const trapezoid = Integration.trapezoidRule(129, 0, 1, Math.exp);
        expect(Math.abs(romberg - exact))
            .toBeLessThan(Math.abs(trapezoid - exact));
    });
});

describe('Integration.computeQuadratureInfo', () => {
    it('produces the degree-2 Gauss-Legendre rule', () => {
        const { roots, coefficients } = Integration.computeQuadratureInfo(2);
        const sorted = roots.slice().sort((a, b) => a - b);
        expect(sorted[0]).toBeCloseTo(-1 / Math.sqrt(3), 10);
        expect(sorted[1]).toBeCloseTo(1 / Math.sqrt(3), 10);
        expect(coefficients[0]).toBeCloseTo(1, 10);
        expect(coefficients[1]).toBeCloseTo(1, 10);
    });

    it('produces the degree-3 Gauss-Legendre rule', () => {
        const { roots, coefficients } = Integration.computeQuadratureInfo(3);
        const pairs = roots.map((r, i) => [r, coefficients[i]] as [number, number])
            .sort((a, b) => a[0] - b[0]);
        const expectedRoot = Math.sqrt(3 / 5);
        expect(pairs[0][0]).toBeCloseTo(-expectedRoot, 10);
        expect(pairs[1][0]).toBeCloseTo(0, 10);
        expect(pairs[2][0]).toBeCloseTo(expectedRoot, 10);
        expect(pairs[0][1]).toBeCloseTo(5 / 9, 10);
        expect(pairs[1][1]).toBeCloseTo(8 / 9, 10);
        expect(pairs[2][1]).toBeCloseTo(5 / 9, 10);
    });

    it('produces roots that are zeros of the Legendre polynomial', () => {
        // P[k](x) = ((2k-1)*x*P[k-1](x) - (k-1)*P[k-2](x))/k.
        const legendre = (degree: number, x: number): number => {
            let pkm2 = 1;
            let pkm1 = x;
            if (degree === 0) {
                return pkm2;
            }
            for (let k = 2; k <= degree; ++k) {
                const pk = ((2 * k - 1) * x * pkm1 - (k - 1) * pkm2) / k;
                pkm2 = pkm1;
                pkm1 = pk;
            }
            return pkm1;
        };

        for (let degree = 2; degree <= 7; ++degree) {
            const { roots, coefficients } =
                Integration.computeQuadratureInfo(degree);
            expect(roots.length).toBe(degree);
            for (const root of roots) {
                expect(Math.abs(root)).toBeLessThan(1);
                expect(legendre(degree, root)).toBeCloseTo(0, 9);
            }
            // The weights are positive and sum to the length of [-1, 1].
            let sum = 0;
            for (const c of coefficients) {
                expect(c).toBeGreaterThan(0);
                sum += c;
            }
            expect(sum).toBeCloseTo(2, 9);
        }
    });
});

describe('Integration.gaussianQuadrature', () => {
    it('integrates polynomials of degree at most 2n-1 exactly', () => {
        for (let degree = 2; degree <= 6; ++degree) {
            const { roots, coefficients } =
                Integration.computeQuadratureInfo(degree);
            for (let p = 0; p <= 2 * degree - 1; ++p) {
                // integral over [0, 1] of t^p is 1/(p+1).
                const value = Integration.gaussianQuadrature(
                    roots, coefficients, 0, 1, t => Math.pow(t, p));
                expect(value).toBeCloseTo(1 / (p + 1), 8);
            }
        }
    });

    it('integrates over a general interval', () => {
        const { roots, coefficients } = Integration.computeQuadratureInfo(6);
        expect(Integration.gaussianQuadrature(
            roots, coefficients, 0, Math.PI, Math.sin)).toBeCloseTo(2, 8);
        expect(Integration.gaussianQuadrature(
            roots, coefficients, -1, 3, t => t * t * t))
            .toBeCloseTo(20, 8);
    });

    it('changes sign when the endpoints are swapped', () => {
        const { roots, coefficients } = Integration.computeQuadratureInfo(4);
        const forward = Integration.gaussianQuadrature(
            roots, coefficients, 0, 2, Math.exp);
        const backward = Integration.gaussianQuadrature(
            roots, coefficients, 2, 0, Math.exp);
        expect(backward).toBeCloseTo(-forward, 12);
    });

    it('agrees with Romberg integration on a smooth integrand', () => {
        const f = (t: number) => Math.cos(t) * Math.exp(-t);
        const { roots, coefficients } = Integration.computeQuadratureInfo(8);
        const gauss = Integration.gaussianQuadrature(
            roots, coefficients, 0, 2, f);
        const romberg = Integration.romberg(12, 0, 2, f);
        expect(gauss).toBeCloseTo(romberg, 9);
    });
});

// ---------------------------------------------------------------------------
// Verification (group V38).
// ---------------------------------------------------------------------------
describe('Integration verification', () => {
    // Horner evaluation of c[0] + c[1] t + ... and the exact antiderivative.
    function polyEval(c: readonly number[], t: number): number {
        let s = 0;
        for (let i = c.length - 1; i >= 0; --i) { s = s * t + c[i]; }
        return s;
    }

    function polyIntegral(c: readonly number[], a: number, b: number): number {
        let s = 0;
        for (let i = 0; i < c.length; ++i) {
            s += c[i] * (Math.pow(b, i + 1) - Math.pow(a, i + 1)) / (i + 1);
        }
        return s;
    }

    // The Legendre polynomial P_n evaluated by its three-term recurrence,
    // independent of the coefficient array Integration builds.
    function legendre(n: number, x: number): number {
        let pkm1 = 1;
        let pk = x;
        if (n === 0) { return pkm1; }
        for (let k = 2; k <= n; ++k) {
            const next = ((2 * k - 1) * x * pk - (k - 1) * pkm1) / k;
            pkm1 = pk;
            pk = next;
        }
        return pk;
    }

    const coeffs = (maxDegree: number) => fc.array(scaled(-3, 3),
        { minLength: 1, maxLength: maxDegree + 1 });

    const bounds = fc.tuple(scaled(-3, 3), scaled(-3, 3));

    it('the trapezoid rule is exact for linear integrands', () => {
        check(fc.tuple(coeffs(1), bounds, fc.integer({ min: 2, max: 20 })),
            ([c, [a, b], numSamples]) => {
                const got = Integration.trapezoidRule(numSamples, a, b,
                    t => polyEval(c, t));
                const exact = polyIntegral(c, a, b);
                expectClose(got, exact, 1e-12, 1e-12);
            });
    });

    it('the trapezoid rule reduces to the two-point formula for numSamples 2',
        () => {
            check(fc.tuple(bounds, coeffs(4)), ([[a, b], c]) => {
                const f = (t: number) => polyEval(c, t);
                expectClose(Integration.trapezoidRule(2, a, b, f),
                    0.5 * (b - a) * (f(a) + f(b)), 1e-12, 1e-12);
            });
        });

    it('Romberg of order k is exact for polynomials of degree at most 2k-1',
        () => {
            check(fc.tuple(fc.integer({ min: 1, max: 6 }), bounds, fc.nat())
                .chain(([order, ab]) => fc.tuple(fc.constant(order),
                    fc.constant(ab), coeffs(2 * order - 1))),
                ([order, [a, b], c]) => {
                    const exact = polyIntegral(c, a, b);
                    const got = Integration.romberg(order, a, b,
                        t => polyEval(c, t));
                    // Richardson extrapolation of the trapezoid rule cancels
                    // the error terms through h^(2k-2), so a polynomial of
                    // degree 2k-1 is integrated exactly up to round-off.
                    expectClose(got, exact, 1e-11, 1e-11);
                });
        });

    it('Gauss-Legendre of degree n is exact for polynomials of degree at most '
        + '2n-1', () => {
            const infos = new Map<number, { roots: number[];
                coefficients: number[] }>();
            for (let n = 2; n <= 8; ++n) {
                infos.set(n, Integration.computeQuadratureInfo(n));
            }
            check(fc.tuple(fc.integer({ min: 2, max: 8 }), bounds, fc.nat())
                .chain(([n, ab]) => fc.tuple(fc.constant(n), fc.constant(ab),
                    coeffs(2 * n - 1))),
                ([n, [a, b], c]) => {
                    const { roots, coefficients } = infos.get(n)!;
                    const got = Integration.gaussianQuadrature(roots,
                        coefficients, a, b, t => polyEval(c, t));
                    expectClose(got, polyIntegral(c, a, b), 1e-10, 1e-10);
                }, 150);
        });

    it('the quadrature roots are zeros of the Legendre polynomial and the '
        + 'coefficients are a positive symmetric partition of 2', () => {
            for (let n = 2; n <= 8; ++n) {
                const { roots, coefficients } =
                    Integration.computeQuadratureInfo(n);
                expect(roots.length).toBe(n);
                expect(coefficients.length).toBe(n);
                let sum = 0;
                for (let i = 0; i < n; ++i) {
                    expect(Math.abs(roots[i])).toBeLessThan(1);
                    expect(Math.abs(legendre(n, roots[i])))
                        .toBeLessThanOrEqual(1e-12);
                    // Gauss-Legendre weights are strictly positive.
                    expect(coefficients[i]).toBeGreaterThan(0);
                    sum += coefficients[i];
                    // The nodes and weights are symmetric about the origin,
                    // and Find returns them in increasing order.
                    expectClose(roots[i], -roots[n - 1 - i], 1e-13, 1e-13);
                    expectClose(coefficients[i], coefficients[n - 1 - i],
                        1e-12, 1e-12);
                    if (i + 1 < n) {
                        expect(roots[i]).toBeLessThan(roots[i + 1]);
                    }
                }
                // The rule integrates the constant 1 over [-1,1] exactly.
                expectClose(sum, 2, 1e-12, 1e-12);
            }
        });

    it('all three rules are linear and reverse sign with the interval', () => {
        const info = Integration.computeQuadratureInfo(5);
        check(fc.tuple(coeffs(4), coeffs(4), bounds, scaled(-3, 3)),
            ([c0, c1, [a, b], s]) => {
                const f = (t: number) => polyEval(c0, t);
                const g = (t: number) => polyEval(c1, t);
                const h = (t: number) => f(t) + s * g(t);
                const rules: Array<(fn: (t: number) => number,
                    lo: number, hi: number) => number> = [
                        (fn, lo, hi) => Integration.trapezoidRule(9, lo, hi, fn),
                        (fn, lo, hi) => Integration.romberg(5, lo, hi, fn),
                        (fn, lo, hi) => Integration.gaussianQuadrature(
                            info.roots, info.coefficients, lo, hi, fn)
                    ];
                for (const rule of rules) {
                    const scale = 1e-10 * (1 + Math.abs(s))
                        * (1 + Math.abs(b - a));
                    expectClose(rule(h, a, b),
                        rule(f, a, b) + s * rule(g, a, b), scale, 1e-10);
                    expectClose(rule(f, b, a), -rule(f, a, b), scale, 1e-10);
                }
            }, 100);
    });

    it('Romberg is additive over a split interval for polynomials', () => {
        check(fc.tuple(coeffs(5), scaled(-3, 3), scaled(-3, 3),
            scaled(-3, 3)), ([c, a, m, b]) => {
                const f = (t: number) => polyEval(c, t);
                const whole = Integration.romberg(6, a, b, f);
                const left = Integration.romberg(6, a, m, f);
                const right = Integration.romberg(6, m, b, f);
                expectClose(whole, left + right, 1e-10, 1e-10);
            });
    });

    it('the three rules agree on smooth transcendental integrands', () => {
        const info = Integration.computeQuadratureInfo(8);
        // The last column is the error an 8-point Gauss-Legendre rule can
        // reach on that integrand: it is spectrally accurate for entire
        // functions but only algebraically so for 1/(1+t^2), whose poles at
        // +-i sit close to the interval.
        const cases: Array<[(t: number) => number, number, number, number,
            number]> = [
                [Math.sin, 0, Math.PI, 2, 1e-8],
                [Math.exp, 0, 1, Math.E - 1, 1e-8],
                [t => 1 / (1 + t * t), -1, 1, Math.PI / 2, 1e-5],
                [t => t * Math.cos(t), 0, 1,
                    Math.cos(1) + Math.sin(1) - 1, 1e-8]
            ];
        for (const [f, a, b, exact, gaussTol] of cases) {
            expectClose(Integration.romberg(12, a, b, f), exact, 1e-11, 1e-11);
            expectClose(Integration.gaussianQuadrature(info.roots,
                info.coefficients, a, b, f), exact, gaussTol, gaussTol);
            // The trapezoid rule converges as O(h^2).
            const t1 = Integration.trapezoidRule(101, a, b, f);
            const t2 = Integration.trapezoidRule(201, a, b, f);
            expect(Math.abs(t2 - exact)).toBeLessThan(
                Math.abs(t1 - exact) * 0.3 + 1e-14);
        }
    });

    it('a degenerate interval integrates to zero', () => {
        const info = Integration.computeQuadratureInfo(4);
        check(fc.tuple(scaled(-3, 3), coeffs(4)), ([a, c]) => {
            const f = (t: number) => polyEval(c, t);
            expect(Integration.trapezoidRule(7, a, a, f) + 0).toBe(0);
            expect(Integration.romberg(4, a, a, f) + 0).toBe(0);
            expect(Integration.gaussianQuadrature(info.roots,
                info.coefficients, a, a, f) + 0).toBe(0);
        });
    });
});
