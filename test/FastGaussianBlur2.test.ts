import { describe, it, expect } from 'vitest';
import { FastGaussianBlur1 } from '../src/FastGaussianBlur1.js';
import { FastGaussianBlur2 } from '../src/FastGaussianBlur2.js';
import { check, expectClose, fc, finite } from './helpers/arbitraries.js';

function sum(values: ArrayLike<number>): number {
    let s = 0;
    for (let i = 0; i < values.length; ++i) {
        s += values[i];
    }
    return s;
}

// Simple deterministic pseudorandom generator for the cross-check.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('FastGaussianBlur2', () => {
    it('leaves a constant image unchanged', () => {
        const xBound = 9, yBound = 7;
        const input = new Array<number>(xBound * yBound).fill(-2.5);
        const output = new Array<number>(xBound * yBound).fill(0);
        const blur = new FastGaussianBlur2();
        blur.execute(xBound, yBound, input, output, 1.5, Math.log(1.2));
        for (let i = 0; i < output.length; ++i) {
            expect(output[i]).toBeCloseTo(-2.5, 12);
        }
    });

    it('preserves total mass and smooths an impulse monotonically', () => {
        const xBound = 33, yBound = 33;
        const base = 1.05;  // stability requires 1 < b < exp(1) for d = 2
        const logBase = Math.log(base);
        let scale = 1.0;
        let input = new Array<number>(xBound * yBound).fill(0);
        const cx = 16, cy = 16;
        input[cx + xBound * cy] = 1;
        const mass0 = sum(input);
        let prevMax = 1;
        const blur = new FastGaussianBlur2();

        for (let iteration = 0; iteration < 8; ++iteration) {
            const output = new Array<number>(xBound * yBound).fill(0);
            blur.execute(xBound, yBound, input, output, scale, logBase);

            expect(sum(output)).toBeCloseTo(mass0, 8);

            const max = Math.max(...output);
            expect(max).toBeLessThan(prevMax);
            expect(max).toBeGreaterThan(0);
            prevMax = max;

            // Reflection symmetry about both axes through the center.
            for (let y = 0; y < yBound; ++y) {
                for (let x = 0; x < xBound; ++x) {
                    const v = output[x + xBound * y];
                    expect(v).toBeCloseTo(output[(xBound - 1 - x) + xBound * y], 10);
                    expect(v).toBeCloseTo(output[x + xBound * (yBound - 1 - y)], 10);
                }
            }

            input = output;
            scale *= base;
        }

        expect(input[cx + xBound * cy]).toBeLessThan(0.5);
        expect(input[cx + 2 + xBound * cy]).toBeGreaterThan(0);
    });

    it('matches per-axis FastGaussianBlur1 passes on random input', () => {
        // The 2D update is output = center + logBase * (xsum + ysum), so it
        // equals rowBlur(input) + columnBlur(input) - input, where each
        // one-dimensional blur is computed by FastGaussianBlur1.
        const xBound = 12, yBound = 10;
        const scale = 1.75;
        const logBase = Math.log(1.3);
        const rand = makeRandom(0x2d2d2d);
        const input = Array.from({ length: xBound * yBound }, () => rand());

        const output = new Array<number>(xBound * yBound).fill(0);
        new FastGaussianBlur2().execute(xBound, yBound, input, output, scale, logBase);

        const blur1 = new FastGaussianBlur1();
        const rowBlur = new Array<number>(xBound * yBound).fill(0);
        for (let y = 0; y < yBound; ++y) {
            const row = input.slice(xBound * y, xBound * (y + 1));
            const rowOut = new Array<number>(xBound).fill(0);
            blur1.execute(xBound, row, rowOut, scale, logBase);
            for (let x = 0; x < xBound; ++x) {
                rowBlur[x + xBound * y] = rowOut[x];
            }
        }
        const colBlur = new Array<number>(xBound * yBound).fill(0);
        for (let x = 0; x < xBound; ++x) {
            const col = new Array<number>(yBound);
            for (let y = 0; y < yBound; ++y) {
                col[y] = input[x + xBound * y];
            }
            const colOut = new Array<number>(yBound).fill(0);
            blur1.execute(yBound, col, colOut, scale, logBase);
            for (let y = 0; y < yBound; ++y) {
                colBlur[x + xBound * y] = colOut[y];
            }
        }

        for (let i = 0; i < output.length; ++i) {
            expect(output[i]).toBeCloseTo(rowBlur[i] + colBlur[i] - input[i], 12);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V24): properties cross-checking the port against the
// upstream FastGaussianBlur2.h finite-difference stencil.
// ---------------------------------------------------------------------------

describe('FastGaussianBlur2 verification', () => {
    const bounds = fc.tuple(fc.integer({ min: 2, max: 7 }), fc.integer({ min: 2, max: 7 }));
    const dyadicScale = fc.constantFrom(0.25, 0.5, 0.75, 1, 1.5, 2);
    const dyadicLogBase = fc.constantFrom(0.0625, 0.125, 0.25);
    const pixels = (n: number) => fc.array(fc.integer({ min: -16, max: 16 }),
        { minLength: n, maxLength: n });

    function blur2(xBound: number, yBound: number, input: readonly number[],
        scale: number, logBase: number): number[] {
        const output = new Array<number>(xBound * yBound).fill(0);
        new FastGaussianBlur2().execute(xBound, yBound, input, output, scale, logBase);
        return output;
    }

    it('equals the sum of the two separable 1-dimensional passes minus the image', () => {
        // out2 = c + logBase * (xsum + ysum) while a row pass gives
        // c + logBase * xsum and a column pass gives c + logBase * ysum, so
        // out2 = rowPass + columnPass - input, exactly.
        check(fc.tuple(bounds, dyadicScale, dyadicLogBase),
            ([[xBound, yBound], scale, logBase]) => {
                fc.assert(fc.property(pixels(xBound * yBound), input => {
                    const output = blur2(xBound, yBound, input, scale, logBase);
                    const blur1 = new FastGaussianBlur1();

                    const rowPass = new Array<number>(xBound * yBound).fill(0);
                    for (let y = 0; y < yBound; ++y) {
                        const row = input.slice(xBound * y, xBound * (y + 1));
                        const out = new Array<number>(xBound).fill(0);
                        blur1.execute(xBound, row, out, scale, logBase);
                        for (let x = 0; x < xBound; ++x) {
                            rowPass[x + xBound * y] = out[x];
                        }
                    }

                    const columnPass = new Array<number>(xBound * yBound).fill(0);
                    for (let x = 0; x < xBound; ++x) {
                        const column = new Array<number>(yBound);
                        for (let y = 0; y < yBound; ++y) {
                            column[y] = input[x + xBound * y];
                        }
                        const out = new Array<number>(yBound).fill(0);
                        blur1.execute(yBound, column, out, scale, logBase);
                        for (let y = 0; y < yBound; ++y) {
                            columnPass[x + xBound * y] = out[y];
                        }
                    }

                    for (let i = 0; i < output.length; ++i) {
                        expectClose(output[i], rowPass[i] + columnPass[i] - input[i],
                            1e-12, 1e-12);
                    }
                }), { numRuns: 6 });
            }, 40);
    });

    it('is linear in the input image', () => {
        check(fc.tuple(bounds, dyadicScale, dyadicLogBase,
            fc.integer({ min: -3, max: 3 }), fc.integer({ min: -3, max: 3 })),
            ([[xBound, yBound], scale, logBase, a, b]) => {
                const n = xBound * yBound;
                fc.assert(fc.property(pixels(n), pixels(n), (u, v) => {
                    const combined = u.map((ui, i) => a * ui + b * v[i]);
                    const lhs = blur2(xBound, yBound, combined, scale, logBase);
                    const bu = blur2(xBound, yBound, u, scale, logBase);
                    const bv = blur2(xBound, yBound, v, scale, logBase);
                    for (let i = 0; i < n; ++i) {
                        expectClose(lhs[i], a * bu[i] + b * bv[i], 1e-12, 1e-12);
                    }
                }), { numRuns: 5 });
            }, 40);
    });

    it('leaves every constant image unchanged', () => {
        check(fc.tuple(bounds, dyadicScale, dyadicLogBase, finite(-100, 100)),
            ([[xBound, yBound], scale, logBase, c]) => {
                const input = new Array<number>(xBound * yBound).fill(c);
                for (const value of blur2(xBound, yBound, input, scale, logBase)) {
                    expectClose(value, c, 1e-12, 1e-12);
                }
            });
    });

    it('commutes with reflection in x and in y', () => {
        check(fc.tuple(bounds, dyadicScale, dyadicLogBase),
            ([[xBound, yBound], scale, logBase]) => {
                fc.assert(fc.property(pixels(xBound * yBound), input => {
                    const output = blur2(xBound, yBound, input, scale, logBase);

                    const flipX = new Array<number>(xBound * yBound);
                    const flipY = new Array<number>(xBound * yBound);
                    for (let y = 0; y < yBound; ++y) {
                        for (let x = 0; x < xBound; ++x) {
                            flipX[x + xBound * y] = input[(xBound - 1 - x) + xBound * y];
                            flipY[x + xBound * y] = input[x + xBound * (yBound - 1 - y)];
                        }
                    }
                    const outX = blur2(xBound, yBound, flipX, scale, logBase);
                    const outY = blur2(xBound, yBound, flipY, scale, logBase);
                    for (let y = 0; y < yBound; ++y) {
                        for (let x = 0; x < xBound; ++x) {
                            expectClose(outX[x + xBound * y],
                                output[(xBound - 1 - x) + xBound * y], 1e-12, 1e-12);
                            expectClose(outY[x + xBound * y],
                                output[x + xBound * (yBound - 1 - y)], 1e-12, 1e-12);
                        }
                    }
                }), { numRuns: 6 });
            }, 40);
    });

    it('transposing the image transposes the result', () => {
        // The x and y sweeps use the same stencil, so the filter commutes
        // with the (x, y) -> (y, x) relabeling; this catches a swapped bound
        // or a swapped index in the linearization i = x + xBound * y.
        check(fc.tuple(bounds, dyadicScale, dyadicLogBase),
            ([[xBound, yBound], scale, logBase]) => {
                fc.assert(fc.property(pixels(xBound * yBound), input => {
                    const output = blur2(xBound, yBound, input, scale, logBase);
                    const transposed = new Array<number>(xBound * yBound);
                    for (let y = 0; y < yBound; ++y) {
                        for (let x = 0; x < xBound; ++x) {
                            transposed[y + yBound * x] = input[x + xBound * y];
                        }
                    }
                    const outT = blur2(yBound, xBound, transposed, scale, logBase);
                    for (let y = 0; y < yBound; ++y) {
                        for (let x = 0; x < xBound; ++x) {
                            expectClose(outT[y + yBound * x], output[x + xBound * y],
                                1e-12, 1e-12);
                        }
                    }
                }), { numRuns: 6 });
            }, 40);
    });
});
