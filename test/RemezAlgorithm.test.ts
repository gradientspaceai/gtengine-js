import { describe, expect, it } from 'vitest';
import { REMEZ_FAILURE, RemezAlgorithm } from '../src/RemezAlgorithm';

// Sample the error F(x) - P(x) densely and return its extreme values.
function sampleError(F: (x: number) => number, coefficients: readonly number[],
    xMin: number, xMax: number, numSamples: number = 20001):
    { maxError: number; minSigned: number; maxSigned: number } {
    let maxError = 0;
    let minSigned = Number.MAX_VALUE;
    let maxSigned = -Number.MAX_VALUE;
    for (let i = 0; i < numSamples; ++i) {
        const x = xMin + (xMax - xMin) * i / (numSamples - 1);
        let p = 0;
        for (let j = coefficients.length - 1; j >= 0; --j) {
            p = coefficients[j] + x * p;
        }
        const e = F(x) - p;
        maxError = Math.max(maxError, Math.abs(e));
        minSigned = Math.min(minSigned, e);
        maxSigned = Math.max(maxSigned, e);
    }
    return { maxError, minSigned, maxSigned };
}

describe('RemezAlgorithm', () => {
    describe('input validation', () => {
        const F = (x: number): number => Math.exp(x);
        const FDer = (x: number): number => Math.exp(x);

        it('throws when xMin >= xMax', () => {
            const remez = new RemezAlgorithm();
            expect(() => remez.execute(F, FDer, 1, 1, 2, 8, 64, 64)).toThrow('Invalid input.');
            expect(() => remez.execute(F, FDer, 1, 0, 2, 8, 64, 64)).toThrow('Invalid input.');
        });

        it('throws when the degree is zero', () => {
            const remez = new RemezAlgorithm();
            expect(() => remez.execute(F, FDer, 0, 1, 0, 8, 64, 64)).toThrow('Invalid input.');
        });

        it('throws when an iteration count is zero', () => {
            const remez = new RemezAlgorithm();
            expect(() => remez.execute(F, FDer, 0, 1, 2, 0, 64, 64)).toThrow('Invalid input.');
            expect(() => remez.execute(F, FDer, 0, 1, 2, 8, 0, 64)).toThrow('Invalid input.');
            expect(() => remez.execute(F, FDer, 0, 1, 2, 8, 64, 0)).toThrow('Invalid input.');
        });
    });

    describe('degree-1 minimax approximation of exp(x) on [0,1]', () => {
        // The best linear approximation a + b*x of e^x on [0,1] is known in
        // closed form. The slope is the divided difference
        //   b = (e^1 - e^0) / (1 - 0) = e - 1,
        // the interior equioscillation point is where the derivative matches
        // the slope,
        //   x* = ln(e - 1),
        // and the intercept is
        //   a = (f(0) + f(x*))/2 - b*x*/2 = e/2 - (e-1)*ln(e-1)/2.
        // The equioscillation error is E = f(0) - a = 1 - a, attained with
        // signs +, -, + at x = 0, x*, 1.
        const b = Math.E - 1;
        const xStar = Math.log(Math.E - 1);
        const a = 0.5 * Math.E - 0.5 * b * xStar;
        const E = 1 - a;

        const F = (x: number): number => Math.exp(x);
        const FDer = (x: number): number => Math.exp(x);

        it('reproduces the closed-form coefficients', () => {
            // Sanity check on the closed-form values themselves.
            expect(b).toBeCloseTo(1.718281828459045, 12);
            expect(a).toBeCloseTo(0.8940665837422168, 12);
            expect(E).toBeCloseTo(0.1059334162577832, 12);

            const remez = new RemezAlgorithm();
            const iterations = remez.execute(F, FDer, 0, 1, 1, 16, 1024, 128);
            expect(iterations).not.toBe(REMEZ_FAILURE);
            expect(iterations).toBeGreaterThan(0);

            const c = remez.getCoefficients();
            expect(c).toHaveLength(2);
            expect(c[0]).toBeCloseTo(a, 12);
            expect(c[1]).toBeCloseTo(b, 12);
            expect(Math.abs(remez.getEstimatedMaxError())).toBeCloseTo(E, 12);
        });

        it('has an equioscillating error with the documented extrema', () => {
            const remez = new RemezAlgorithm();
            remez.execute(F, FDer, 0, 1, 1, 16, 1024, 128);

            const xNodes = remez.getXNodes();
            expect(xNodes).toHaveLength(3);
            expect(xNodes[0]).toBe(0);
            expect(xNodes[1]).toBeCloseTo(xStar, 10);
            expect(xNodes[2]).toBe(1);

            const errors = remez.getErrors();
            expect(errors).toHaveLength(3);
            // The signs alternate + - +.
            expect(errors[0]).toBeGreaterThan(0);
            expect(errors[1]).toBeLessThan(0);
            expect(errors[2]).toBeGreaterThan(0);
            // All three have the same magnitude E (the equioscillation
            // property that characterizes the minimax polynomial).
            for (const e of errors) {
                expect(Math.abs(e)).toBeCloseTo(E, 10);
            }
        });

        it('produces the true minimax error over a dense sampling', () => {
            const remez = new RemezAlgorithm();
            remez.execute(F, FDer, 0, 1, 1, 16, 1024, 128);
            const sampled = sampleError(F, remez.getCoefficients(), 0, 1);
            expect(sampled.maxError).toBeCloseTo(E, 8);
            expect(sampled.maxSigned).toBeCloseTo(E, 8);
            expect(sampled.minSigned).toBeCloseTo(-E, 8);
        });

        it('beats the Taylor and interpolating polynomials of the same degree', () => {
            const remez = new RemezAlgorithm();
            remez.execute(F, FDer, 0, 1, 1, 16, 1024, 128);
            const minimax = sampleError(F, remez.getCoefficients(), 0, 1).maxError;

            // Degree-1 Taylor polynomial about x = 0.
            const taylor = sampleError(F, [1, 1], 0, 1).maxError;
            // Linear interpolant through the endpoints.
            const interp = sampleError(F, [1, Math.E - 1], 0, 1).maxError;

            expect(minimax).toBeLessThan(taylor);
            expect(minimax).toBeLessThan(interp);
            // The interpolant's error is exactly twice the minimax error for
            // a convex function on an interval.
            expect(interp).toBeCloseTo(2 * minimax, 6);
        });
    });

    describe('higher-degree approximations', () => {
        it('equioscillates for a degree-2 approximation of sin on [0, pi/2]', () => {
            const F = (x: number): number => Math.sin(x);
            const FDer = (x: number): number => Math.cos(x);
            const remez = new RemezAlgorithm();
            const iterations = remez.execute(F, FDer, 0, 0.5 * Math.PI, 2, 32, 1024, 128);
            expect(iterations).not.toBe(REMEZ_FAILURE);

            const errors = remez.getErrors();
            expect(errors).toHaveLength(4);
            const magnitude = Math.abs(errors[0]);
            expect(magnitude).toBeGreaterThan(0);
            for (let i = 0; i < errors.length; ++i) {
                expect(Math.abs(errors[i])).toBeCloseTo(magnitude, 8);
                // Alternating signs starting positive at x = xMin.
                expect(Math.sign(errors[i])).toBe(i % 2 === 0 ? 1 : -1);
            }
            expect(Math.abs(remez.getEstimatedMaxError())).toBeCloseTo(magnitude, 8);

            // No point of the interval has a larger error.
            const sampled = sampleError(F, remez.getCoefficients(), 0, 0.5 * Math.PI);
            expect(sampled.maxError).toBeCloseTo(magnitude, 8);
        });

        it('improves monotonically as the degree increases', () => {
            const F = (x: number): number => Math.exp(x);
            const FDer = (x: number): number => Math.exp(x);
            let previous = Number.MAX_VALUE;
            for (let degree = 1; degree <= 6; ++degree) {
                const remez = new RemezAlgorithm();
                const iterations = remez.execute(F, FDer, 0, 1, degree, 32, 1024, 128);
                expect(iterations).not.toBe(REMEZ_FAILURE);
                const error = Math.abs(remez.getEstimatedMaxError());
                expect(error).toBeGreaterThan(0);
                expect(error).toBeLessThan(previous);
                expect(remez.getCoefficients()).toHaveLength(degree + 1);
                expect(remez.getXNodes()).toHaveLength(degree + 2);
                previous = error;
            }
            // Degree 6 approximates exp on [0,1] to about 4e-8, which is the
            // magnitude predicted by the Chebyshev coefficient
            // 2*exp(1/2)*I_7(1/2) of exp((t+1)/2) on [-1,1].
            expect(previous).toBeLessThan(1e-07);
            expect(previous).toBeGreaterThan(1e-09);
        });

        it('recovers the low-order coefficients of exp near the origin', () => {
            // On a narrow interval the minimax polynomial is close to the
            // Taylor polynomial, so the coefficients approach 1, 1, 1/2.
            const F = (x: number): number => Math.exp(x);
            const FDer = (x: number): number => Math.exp(x);
            const remez = new RemezAlgorithm();
            remez.execute(F, FDer, -0.01, 0.01, 2, 32, 1024, 128);
            const c = remez.getCoefficients();
            expect(c[0]).toBeCloseTo(1, 6);
            expect(c[1]).toBeCloseTo(1, 4);
            expect(c[2]).toBeCloseTo(0.5, 2);
            // The leading neglected term is x^3/6 on |x| <= 0.01, reduced by
            // the Chebyshev factor 1/4, i.e. about 4.2e-8.
            expect(Math.abs(remez.getEstimatedMaxError())).toBeLessThan(1e-07);
        });
    });

    describe('exact-fit and termination behavior', () => {
        it('reports failure when F is itself a polynomial of the given degree', () => {
            // The approximation is exact, so the errors are all zero and the
            // oscillation test fails immediately; upstream signals this by
            // returning the size_t sentinel.
            const F = (x: number): number => 2 * x * x - 3 * x + 1;
            const FDer = (x: number): number => 4 * x - 3;
            const remez = new RemezAlgorithm();
            const iterations = remez.execute(F, FDer, -1, 1, 2, 16, 512, 128);
            expect(iterations).toBe(REMEZ_FAILURE);

            // The exact coefficients were still computed before the test.
            const c = remez.getCoefficients();
            expect(c[0]).toBeCloseTo(1, 10);
            expect(c[1]).toBeCloseTo(-3, 10);
            expect(c[2]).toBeCloseTo(2, 10);
            expect(Math.abs(remez.getEstimatedMaxError())).toBeLessThan(1e-12);
        });

        it('returns the iteration cap when the budget is exhausted', () => {
            const F = (x: number): number => Math.exp(x);
            const FDer = (x: number): number => Math.exp(x);
            const remez = new RemezAlgorithm();
            const iterations = remez.execute(F, FDer, 0, 1, 3, 2, 1024, 128);
            expect(iterations).toBe(2);
        });

        it('converges in few iterations and is stable under extra ones', () => {
            const F = (x: number): number => Math.exp(x);
            const FDer = (x: number): number => Math.exp(x);

            const few = new RemezAlgorithm();
            few.execute(F, FDer, 0, 1, 3, 8, 1024, 128);
            const many = new RemezAlgorithm();
            many.execute(F, FDer, 0, 1, 3, 64, 1024, 128);

            expect(many.getCoefficients().length).toBe(few.getCoefficients().length);
            for (let i = 0; i < few.getCoefficients().length; ++i) {
                expect(many.getCoefficients()[i]).toBeCloseTo(few.getCoefficients()[i], 12);
            }
            expect(many.getEstimatedMaxError()).toBeCloseTo(few.getEstimatedMaxError(), 14);
        });

        it('can be reused for a second problem', () => {
            const remez = new RemezAlgorithm();
            remez.execute((x) => Math.exp(x), (x) => Math.exp(x), 0, 1, 2, 32, 1024, 128);
            const first = remez.getCoefficients().slice();

            remez.execute((x) => Math.sin(x), (x) => Math.cos(x), 0, 1, 4, 32, 1024, 128);
            expect(remez.getCoefficients()).toHaveLength(5);
            const sampled = sampleError(Math.sin, remez.getCoefficients(), 0, 1);
            expect(sampled.maxError).toBeLessThan(5e-05);

            remez.execute((x) => Math.exp(x), (x) => Math.exp(x), 0, 1, 2, 32, 1024, 128);
            for (let i = 0; i < first.length; ++i) {
                expect(remez.getCoefficients()[i]).toBeCloseTo(first[i], 12);
            }
        });
    });

    describe('accessors before execute', () => {
        it('has empty outputs on a fresh object', () => {
            const remez = new RemezAlgorithm();
            expect(remez.getCoefficients()).toEqual([]);
            expect(remez.getXNodes()).toEqual([]);
            expect(remez.getErrors()).toEqual([]);
            expect(remez.getEstimatedMaxError()).toBe(0);
        });
    });

    describe('interval handling', () => {
        it('produces the same approximation under a shift of the interval', () => {
            // Approximating G(x) = F(x - shift) on [xMin+shift, xMax+shift]
            // yields the same error magnitude as F on [xMin, xMax].
            const shift = 3;
            const direct = new RemezAlgorithm();
            direct.execute((x) => Math.exp(x), (x) => Math.exp(x), 0, 1, 3, 32, 1024, 128);

            const shifted = new RemezAlgorithm();
            shifted.execute((x) => Math.exp(x - shift), (x) => Math.exp(x - shift),
                shift, 1 + shift, 3, 32, 1024, 128);

            expect(Math.abs(shifted.getEstimatedMaxError()))
                .toBeCloseTo(Math.abs(direct.getEstimatedMaxError()), 12);
            const directNodes = direct.getXNodes();
            const shiftedNodes = shifted.getXNodes();
            for (let i = 0; i < directNodes.length; ++i) {
                expect(shiftedNodes[i] - shift).toBeCloseTo(directNodes[i], 8);
            }
        });

        it('handles an interval that straddles the origin', () => {
            const F = (x: number): number => Math.cos(x);
            const FDer = (x: number): number => -Math.sin(x);
            const remez = new RemezAlgorithm();
            const iterations = remez.execute(F, FDer, -1, 1, 5, 32, 1024, 128);
            expect(iterations).not.toBe(REMEZ_FAILURE);

            // cos is even, so the odd-order coefficients are essentially zero.
            const c = remez.getCoefficients();
            expect(c).toHaveLength(6);
            expect(Math.abs(c[1])).toBeLessThan(1e-10);
            expect(Math.abs(c[3])).toBeLessThan(1e-10);
            expect(Math.abs(c[5])).toBeLessThan(1e-10);
            expect(c[0]).toBeCloseTo(1, 4);

            const sampled = sampleError(F, c, -1, 1);
            expect(sampled.maxError).toBeLessThan(1e-04);
            expect(1 - c[0]).toBeCloseTo(sampled.maxError, 6);
        });

        it('signals failure for an even function, even degree, symmetric interval', () => {
            // For an even F on [-a,a] the error of the even minimax polynomial
            // of even degree n satisfies E(-a) = E(a), but the algorithm's
            // n+2 nodes require alternating signs at the two endpoints when
            // n is even. No oscillatory node set exists, so the sentinel is
            // returned. (The degree-n and degree-(n+1) minimax polynomials
            // coincide in this configuration.)
            const F = (x: number): number => Math.cos(x);
            const FDer = (x: number): number => -Math.sin(x);
            const remez = new RemezAlgorithm();
            expect(remez.execute(F, FDer, -1, 1, 4, 32, 1024, 128)).toBe(REMEZ_FAILURE);
        });
    });
});
