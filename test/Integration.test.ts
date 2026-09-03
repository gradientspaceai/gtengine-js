import { describe, it, expect } from 'vitest';
import { Integration } from '../src/Integration.js';

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
