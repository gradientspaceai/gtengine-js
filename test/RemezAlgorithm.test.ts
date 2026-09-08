import { describe, expect, it } from 'vitest';
import { REMEZ_FAILURE, RemezAlgorithm } from '../src/RemezAlgorithm.js';
import {
    check, expectClose, fc, scaled
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V39): independent review against upstream RemezAlgorithm.h.
//
// Upstream limitation found here and preserved by the port: ComputeXExtremes
// pins the first and last x-nodes to xMin and xMax, so the alternation set is
// forced to contain both endpoints. When the true alternation set does not -
// which happens as soon as the error function has an interior extremum just
// inside an endpoint - the exchange converges to a fixed point that is not
// the minimax polynomial, and the errors at those (wrong) nodes still
// alternate, so the algorithm reports success. The last test in this block
// pins the numbers for one such case.
// ---------------------------------------------------------------------------

describe('RemezAlgorithm verification', () => {
    type Sample = {
        F: (x: number) => number;
        FDer: (x: number) => number;
    };

    // Functions with a derivative of constant sign and no inflection point on
    // [-1,1]: their minimax error function attains its extreme values at the
    // interval endpoints, which is the configuration upstream's node update
    // assumes.
    const convexFamily = fc.tuple(fc.boolean(), scaled(0.6, 2))
        .map(([isExp, a]): Sample => (isExp
            ? {
                F: (x: number) => Math.exp(a * x),
                FDer: (x: number) => a * Math.exp(a * x)
            }
            : {
                F: (x: number) => 1 / (x + a + 1.5),
                FDer: (x: number) => -1 / ((x + a + 1.5) * (x + a + 1.5))
            }));

    // The convex family plus a sinusoid, whose inflection point inside the
    // interval is what defeats the pinned end nodes.
    const family: fc.Arbitrary<Sample> = fc.oneof(convexFamily,
        scaled(0.6, 2).map((a): Sample => ({
            F: (x: number) => Math.sin(a * x + 0.3),
            FDer: (x: number) => a * Math.cos(a * x + 0.3)
        })));

    // The largest |F(x) - P(x)| on [xMin,xMax], sampled densely.
    function maxError(F: (x: number) => number,
        coefficients: readonly number[], xMin: number, xMax: number): number {
        const samples = 4001;
        let worst = 0;
        for (let i = 0; i <= samples; ++i) {
            const x = xMin + ((xMax - xMin) * i) / samples;
            let p = 0;
            for (let k = coefficients.length - 1; k >= 0; --k) {
                p = coefficients[k] + x * p;
            }
            worst = Math.max(worst, Math.abs(F(x) - p));
        }
        return worst;
    }

    // The degree-n polynomial interpolating F at the n+1 Chebyshev points of
    // [xMin,xMax], in Lagrange form. It is a near-best approximation, so its
    // error is a meaningful upper bound for the minimax error, and it is not
    // an implementation of the algorithm under test.
    function chebyshevInterpolantMaxError(F: (x: number) => number, n: number,
        xMin: number, xMax: number): number {
        const center = 0.5 * (xMax + xMin), radius = 0.5 * (xMax - xMin);
        const nodes: number[] = [], values: number[] = [];
        for (let k = 0; k <= n; ++k) {
            const x = center + radius
                * Math.cos((Math.PI * (2 * k + 1)) / (2 * n + 2));
            nodes.push(x);
            values.push(F(x));
        }
        const samples = 4001;
        let worst = 0;
        for (let i = 0; i <= samples; ++i) {
            const x = xMin + ((xMax - xMin) * i) / samples;
            let p = 0;
            for (let k = 0; k <= n; ++k) {
                let basis = 1;
                for (let j = 0; j <= n; ++j) {
                    if (j !== k) {
                        basis *= (x - nodes[j]) / (nodes[k] - nodes[j]);
                    }
                }
                p += values[k] * basis;
            }
            worst = Math.max(worst, Math.abs(F(x) - p));
        }
        return worst;
    }

    it('equioscillates at the x-nodes with the estimated maximum error',
        () => {
            check(fc.tuple(family, fc.integer({ min: 1, max: 5 })),
                ([sample, degree]) => {
                    const remez = new RemezAlgorithm();
                    const iterations = remez.execute(sample.F, sample.FDer,
                        -1, 1, degree, 32, 1024, 128);
                    // The sentinel means no oscillatory node set was found;
                    // that case has its own tests above.
                    fc.pre(iterations !== REMEZ_FAILURE);

                    const errors = remez.getErrors();
                    const estimate = Math.abs(remez.getEstimatedMaxError());
                    expect(errors.length).toBe(degree + 2);
                    expect(estimate).toBeGreaterThan(0);
                    for (let i = 0; i < errors.length; ++i) {
                        // Equal magnitude at every node ...
                        expectClose(Math.abs(errors[i]), estimate, 0, 1e-4);
                        // ... with alternating signs.
                        if (i > 0) {
                            expect(errors[i] * errors[i - 1])
                                .toBeLessThanOrEqual(0);
                        }
                    }

                    // The level at the nodes is a lower bound for the true
                    // maximum error of the returned polynomial, since the
                    // nodes lie in the interval.
                    const sampled = maxError(sample.F, remez.getCoefficients(),
                        -1, 1);
                    expect(sampled).toBeGreaterThanOrEqual(estimate
                        * (1 - 1e-9));
                }, 100);
        });

    it('is the minimax polynomial when the extremes are at the endpoints',
        () => {
            check(fc.tuple(convexFamily, fc.integer({ min: 1, max: 5 })),
                ([sample, degree]) => {
                    const remez = new RemezAlgorithm();
                    const iterations = remez.execute(sample.F, sample.FDer,
                        -1, 1, degree, 32, 1024, 128);
                    fc.pre(iterations !== REMEZ_FAILURE);

                    const estimate = Math.abs(remez.getEstimatedMaxError());
                    const sampled = maxError(sample.F,
                        remez.getCoefficients(), -1, 1);
                    // The equioscillation level is the true maximum error: on
                    // this family the observed gap over 200 values of the
                    // shape parameter and five degrees never exceeds 1e-10
                    // relative.
                    expectClose(sampled, estimate, 0, 1e-8);

                    // And no polynomial of the same degree does better, in
                    // particular the Chebyshev interpolant (observed ratio at
                    // most 0.92).
                    const reference = chebyshevInterpolantMaxError(sample.F,
                        degree, -1, 1);
                    expect(sampled).toBeLessThanOrEqual(reference
                        * (1 + 1e-9));
                }, 100);
        });

    it('is invariant under an affine change of the interval', () => {
        check(fc.tuple(family, fc.integer({ min: 1, max: 4 }),
            scaled(-3, 3), scaled(0.5, 4)), ([sample, degree, c, r]) => {
            // G(y) = F((y - c)/r) on [c - r, c + r] is F on [-1,1]
            // reparameterized, so the minimax error is the same.
            const G = (y: number) => sample.F((y - c) / r);
            const GDer = (y: number) => sample.FDer((y - c) / r) / r;
            const base = new RemezAlgorithm();
            const it0 = base.execute(sample.F, sample.FDer, -1, 1, degree, 32,
                1024, 128);
            const moved = new RemezAlgorithm();
            const it1 = moved.execute(G, GDer, c - r, c + r, degree, 32, 1024,
                128);
            fc.pre(it0 !== REMEZ_FAILURE && it1 !== REMEZ_FAILURE);
            // The two runs perform different arithmetic, so only the value of
            // the minimax error is comparable, not the digits.
            expectClose(Math.abs(moved.getEstimatedMaxError()),
                Math.abs(base.getEstimatedMaxError()), 0, 1e-6);
        }, 100);
    });

    it('does no worse as the degree increases', () => {
        check(convexFamily, sample => {
            let previous = Number.POSITIVE_INFINITY;
            for (let degree = 1; degree <= 5; ++degree) {
                const remez = new RemezAlgorithm();
                const iterations = remez.execute(sample.F, sample.FDer, -1, 1,
                    degree, 32, 1024, 128);
                if (iterations === REMEZ_FAILURE) { continue; }
                const estimate = Math.abs(remez.getEstimatedMaxError());
                expect(estimate).toBeLessThanOrEqual(previous);
                previous = estimate;
            }
        }, 100);
    });

    it('reproduces the textbook degree-1 minimax approximation of exp on '
        + '[0,1]', () => {
        // The best linear approximation of a strictly convex F on [a,b] has
        // slope m = (F(b)-F(a))/(b-a), touches F where F'(c) = m, and its
        // error is E = (F(a) - F(c) + m*(c-a))/2. For F = exp on [0,1] this
        // is the classical P(x) = 0.89406 + 1.71828 x with E = 0.10593.
        const m = Math.E - 1;
        const c = Math.log(m);
        const E = (1 - m + m * c) / 2;
        const p0 = (1 + m) / 2 - (m * c) / 2;
        const remez = new RemezAlgorithm();
        const iterations = remez.execute(x => Math.exp(x), x => Math.exp(x),
            0, 1, 1, 32, 1024, 128);
        expect(iterations).not.toBe(REMEZ_FAILURE);
        const coefficients = remez.getCoefficients();
        expectClose(coefficients[0], p0, 1e-12, 1e-12);
        expectClose(coefficients[1], m, 1e-12, 1e-12);
        expectClose(Math.abs(remez.getEstimatedMaxError()), E, 1e-12, 1e-12);
    });

    it('is at least as accurate as the upstream SinEstimate table', () => {
        // SinEstimate.h approximates sin(x) on [-pi/2, pi/2] by an odd
        // polynomial that interpolates sin at the origin; the Remez
        // polynomial of the same total degree is the unconstrained minimax
        // approximation, so it must be at least as accurate. The bounds are
        // the C_SIN_EST_MAX_ERROR entries of src/SinEstimate.ts.
        //
        // The interval is [0, pi/2]: the SinEstimate polynomial is odd, so
        // its maximum error over [-pi/2, pi/2] is attained on [0, pi/2] as
        // well. Running Remez on the symmetric interval instead hits the
        // no-alternation case, because for an odd F the error function is odd,
        // its extrema come in +/- pairs, and an odd degree needs an odd
        // number of alternation points.
        const upstream: [number, number][] = [[3, 1.3481903639146e-2],
        [5, 1.4001209384651e-4], [7, 1.0205878939740e-6],
        [9, 5.2010783457846e-9], [11, 1.9323431743601e-11]];
        for (const [degree, bound] of upstream) {
            const remez = new RemezAlgorithm();
            const iterations = remez.execute(x => Math.sin(x),
                x => Math.cos(x), 0, Math.PI / 2, degree, 32, 1024, 128);
            expect(iterations).not.toBe(REMEZ_FAILURE);
            expect(Math.abs(remez.getEstimatedMaxError()))
                .toBeLessThanOrEqual(bound);
        }
    });

    it('records the upstream node-pinning limitation', () => {
        // F(x) = sin(2x + 0.3) on [-1,1] has an inflection point at
        // x = -0.15, and the alternation set of its degree-1 minimax
        // polynomial does not contain x = -1. ComputeXExtremes pins the first
        // node there anyway, so the exchange converges to a different fixed
        // point: the errors at the three nodes alternate with magnitude
        // 0.33346, the algorithm reports success, and yet the polynomial it
        // returns has a true maximum error of 0.49382 - while the actual
        // minimax error, found by a direct search over the two coefficients,
        // is 0.36192. The reported estimate therefore understates the error
        // of the returned polynomial by a third.
        //
        // The port preserves this: a fix means letting the outermost nodes
        // move into the interior, which is a different node-exchange step
        // rather than a local correction, and no file in the library consumes
        // RemezAlgorithm. See the PR's upstream-bug section.
        const F = (x: number) => Math.sin(2 * x + 0.3);
        const FDer = (x: number) => 2 * Math.cos(2 * x + 0.3);
        const remez = new RemezAlgorithm();
        const iterations = remez.execute(F, FDer, -1, 1, 1, 32, 1024, 128);
        expect(iterations).not.toBe(REMEZ_FAILURE);

        const errors = remez.getErrors();
        for (let i = 0; i < errors.length; ++i) {
            expectClose(Math.abs(errors[i]), 0.3334601918959971, 1e-12, 1e-12);
            if (i > 0) { expect(errors[i] * errors[i - 1]).toBeLessThan(0); }
        }
        expect(remez.getXNodes()[0]).toBe(-1);
        expect(remez.getXNodes()[2]).toBe(1);

        const sampled = maxError(F, remez.getCoefficients(), -1, 1);
        expectClose(sampled, 0.4938154705205917, 1e-6, 1e-6);
        expect(sampled / Math.abs(remez.getEstimatedMaxError()))
            .toBeGreaterThan(1.4);
    });
});
