import { describe, it, expect } from 'vitest';
import { FastGaussianBlur1 } from '../src/FastGaussianBlur1.js';
import { check, expectClose, fc, finite } from './helpers/arbitraries.js';

function sum(values: ArrayLike<number>): number {
    let s = 0;
    for (let i = 0; i < values.length; ++i) {
        s += values[i];
    }
    return s;
}

describe('FastGaussianBlur1', () => {
    it('leaves a constant image unchanged', () => {
        const xBound = 17;
        const input = new Array<number>(xBound).fill(3.25);
        const output = new Array<number>(xBound).fill(0);
        const blur = new FastGaussianBlur1();
        blur.execute(xBound, input, output, 1.5, Math.log(1.2));
        for (let x = 0; x < xBound; ++x) {
            expect(output[x]).toBeCloseTo(3.25, 12);
        }
    });

    it('preserves total mass and smooths an impulse monotonically', () => {
        const xBound = 65;
        const base = 1.05;  // stability requires 1 < b < exp(0.5) for d = 1
        const logBase = Math.log(base);
        let scale = 1.0;
        let input = new Array<number>(xBound).fill(0);
        input[32] = 1;
        const mass0 = sum(input);
        let prevMax = 1;
        const blur = new FastGaussianBlur1();

        for (let iteration = 0; iteration < 10; ++iteration) {
            const output = new Array<number>(xBound).fill(0);
            blur.execute(xBound, input, output, scale, logBase);

            // Total mass is preserved (impulse support stays away from the
            // boundary, and linear interpolation weights sum to one).
            expect(sum(output)).toBeCloseTo(mass0, 8);

            // The peak decays monotonically (heat equation maximum
            // principle).
            const max = Math.max(...output);
            expect(max).toBeLessThan(prevMax);
            expect(max).toBeGreaterThan(0);
            prevMax = max;

            // Symmetric input stays symmetric.
            for (let x = 0; x < xBound; ++x) {
                expect(output[x]).toBeCloseTo(output[xBound - 1 - x], 10);
            }

            input = output;
            scale *= base;
        }

        // After several iterations the impulse has visibly spread.
        expect(input[32]).toBeLessThan(0.5);
        expect(input[30]).toBeGreaterThan(0);
        expect(input[34]).toBeGreaterThan(0);
    });

    it('works with typed array input and output', () => {
        const xBound = 33;
        const input = new Float64Array(xBound);
        input[16] = 8;
        const output = new Float64Array(xBound);
        const blur = new FastGaussianBlur1();
        blur.execute(xBound, input, output, 1.0, Math.log(1.1));
        expect(sum(output)).toBeCloseTo(8, 10);
        expect(output[16]).toBeLessThan(8);
    });

    it('uses boundary values for samples that reach outside the image', () => {
        const xBound = 5;
        // Large scale forces the boundary-value branches at every pixel.
        const input = [1, 0, 0, 0, 2];
        const output = new Array<number>(xBound).fill(0);
        const logBase = Math.log(1.2);
        const blur = new FastGaussianBlur1();
        blur.execute(xBound, input, output, 10.0, logBase);
        for (let x = 0; x < xBound; ++x) {
            // xsum = -2 input[x] + input[xBound-1] + input[0] = -2 input[x] + 3.
            const expected = input[x] + logBase * (-2 * input[x] + 3);
            expect(output[x]).toBeCloseTo(expected, 12);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V24): properties cross-checking the port against the
// upstream FastGaussianBlur1.h finite-difference stencil.
// ---------------------------------------------------------------------------

describe('FastGaussianBlur1 verification', () => {
    // Small integer images: every stencil weight below is a dyadic rational
    // for the scales used, so the identities hold exactly in binary64.
    const image = (n: number) => fc.array(fc.integer({ min: -16, max: 16 }),
        { minLength: n, maxLength: n });
    const bound = fc.integer({ min: 2, max: 12 });
    // Scales whose fractional part is a negative power of two keep the
    // linear-interpolation weights exact.
    const dyadicScale = fc.constantFrom(0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3);
    const dyadicLogBase = fc.constantFrom(0.0625, 0.125, 0.25, 0.5);

    function blur(input: readonly number[], scale: number, logBase: number): number[] {
        const output = new Array<number>(input.length).fill(0);
        new FastGaussianBlur1().execute(input.length, input, output, scale, logBase);
        return output;
    }

    it('is linear in the input image', () => {
        check(fc.tuple(bound, dyadicScale, dyadicLogBase,
            fc.integer({ min: -3, max: 3 }), fc.integer({ min: -3, max: 3 })),
            ([n, scale, logBase, a, b]) => {
                fc.assert(fc.property(image(n), image(n), (u, v) => {
                    const combined = u.map((ui, i) => a * ui + b * v[i]);
                    const lhs = blur(combined, scale, logBase);
                    const bu = blur(u, scale, logBase);
                    const bv = blur(v, scale, logBase);
                    for (let i = 0; i < n; ++i) {
                        expectClose(lhs[i], a * bu[i] + b * bv[i], 1e-12, 1e-12);
                    }
                }), { numRuns: 5 });
            }, 40);
    });

    it('leaves every constant image unchanged (zero discrete Laplacian)', () => {
        check(fc.tuple(bound, dyadicScale, dyadicLogBase, finite(-100, 100)),
            ([n, scale, logBase, c]) => {
                const output = blur(new Array<number>(n).fill(c), scale, logBase);
                for (const value of output) {
                    expectClose(value, c, 1e-12, 1e-12);
                }
            });
    });

    it('commutes with reflection of the image', () => {
        // The stencil samples x + scale and x - scale symmetrically and uses
        // the near boundary value on both ends, so reversing the input
        // reverses the output exactly.
        check(fc.tuple(bound, dyadicScale, dyadicLogBase), ([n, scale, logBase]) => {
            fc.assert(fc.property(image(n), input => {
                const forward = blur(input, scale, logBase);
                const reversed = blur([...input].reverse(), scale, logBase);
                for (let i = 0; i < n; ++i) {
                    expectClose(reversed[i], forward[n - 1 - i], 1e-12, 1e-12);
                }
            }), { numRuns: 8 });
        }, 40);
    });

    it('matches a direct evaluation of the upstream stencil', () => {
        // Independent re-derivation: sample the image at x +- scale with
        // linear interpolation, clamped to the end samples, and form
        // center + logBase * (f(x+s) - 2 f(x) + f(x-s)).
        check(fc.tuple(bound, dyadicScale, dyadicLogBase), ([n, scale, logBase]) => {
            fc.assert(fc.property(image(n), input => {
                const output = blur(input, scale, logBase);
                for (let x = 0; x < n; ++x) {
                    const plus = x + scale;
                    const minus = x - scale;
                    const xp1 = Math.floor(plus);
                    const xm1 = Math.ceil(minus);
                    const fPlus = xp1 >= n - 1
                        ? input[n - 1]
                        : input[xp1] + (plus - xp1) * (input[xp1 + 1] - input[xp1]);
                    const fMinus = xm1 <= 0
                        ? input[0]
                        : input[xm1] + (minus - xm1) * (input[xm1] - input[xm1 - 1]);
                    const expected = input[x]
                        + logBase * (fPlus - 2 * input[x] + fMinus);
                    expectClose(output[x], expected, 1e-12, 1e-12);
                }
            }), { numRuns: 8 });
        }, 40);
    });

    it('reduces to the end samples when the scale spans the image', () => {
        // For scale >= xBound - 1 both branches take the boundary values, so
        // the update is input[0] + input[n-1] - 2 * input[x].
        check(fc.tuple(bound, dyadicLogBase), ([n, logBase]) => {
            fc.assert(fc.property(image(n), input => {
                const output = blur(input, n, logBase);
                for (let x = 0; x < n; ++x) {
                    expectClose(output[x],
                        input[x] + logBase * (input[0] + input[n - 1] - 2 * input[x]),
                        1e-12, 1e-12);
                }
            }), { numRuns: 8 });
        }, 40);
    });

    it('writes integer results through an Int32Array output', () => {
        // Upstream casts to T on store; Int32Array assignment truncates
        // toward zero exactly as static_cast<int32_t> does.
        check(fc.tuple(bound, dyadicScale), ([n, scale]) => {
            fc.assert(fc.property(image(n), input => {
                const asDoubles = blur(input, scale, 0.25);
                const asInts = new Int32Array(n);
                new FastGaussianBlur1().execute(n, input, asInts, scale, 0.25);
                for (let x = 0; x < n; ++x) {
                    // Normalize -0 to 0: Int32Array stores +0 for a truncated -0.3.
                    expect(asInts[x]).toBe(Math.trunc(asDoubles[x]) + 0);
                }
            }), { numRuns: 8 });
        }, 40);
    });
});
