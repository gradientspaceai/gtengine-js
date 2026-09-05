import { describe, it, expect } from 'vitest';
import { FastGaussianBlur1 } from '../src/FastGaussianBlur1.js';
import { FastGaussianBlur3 } from '../src/FastGaussianBlur3.js';
import { check, expectClose, fc, finite } from './helpers/arbitraries.js';

function sum(values: ArrayLike<number>): number {
    let s = 0;
    for (let i = 0; i < values.length; ++i) {
        s += values[i];
    }
    return s;
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('FastGaussianBlur3', () => {
    it('leaves a constant image unchanged', () => {
        const xBound = 6, yBound = 5, zBound = 4;
        const n = xBound * yBound * zBound;
        const input = new Array<number>(n).fill(1.75);
        const output = new Array<number>(n).fill(0);
        const blur = new FastGaussianBlur3();
        blur.execute(xBound, yBound, zBound, input, output, 1.5, Math.log(1.2));
        for (let i = 0; i < n; ++i) {
            expect(output[i]).toBeCloseTo(1.75, 12);
        }
    });

    it('preserves total mass and smooths an impulse monotonically', () => {
        const xBound = 17, yBound = 17, zBound = 17;
        const n = xBound * yBound * zBound;
        const base = 1.05;  // stability requires 1 < b < exp(1.5) for d = 3
        const logBase = Math.log(base);
        let scale = 1.0;
        let input = new Array<number>(n).fill(0);
        const c = 8;
        const center = c + xBound * (c + yBound * c);
        input[center] = 1;
        const mass0 = sum(input);
        let prevMax = 1;
        const blur = new FastGaussianBlur3();

        for (let iteration = 0; iteration < 5; ++iteration) {
            const output = new Array<number>(n).fill(0);
            blur.execute(xBound, yBound, zBound, input, output, scale, logBase);

            expect(sum(output)).toBeCloseTo(mass0, 8);

            const max = Math.max(...output);
            expect(max).toBeLessThan(prevMax);
            expect(max).toBeGreaterThan(0);
            prevMax = max;

            input = output;
            scale *= base;
        }

        // Spread away from the impulse in all three axes.
        expect(input[center]).toBeLessThan(0.5);
        expect(input[(c + 1) + xBound * (c + yBound * c)]).toBeGreaterThan(0);
        expect(input[c + xBound * ((c + 1) + yBound * c)]).toBeGreaterThan(0);
        expect(input[c + xBound * (c + yBound * (c + 1))]).toBeGreaterThan(0);

        // Axis-permutation symmetry of the impulse response.
        for (let z = 0; z < zBound; ++z) {
            for (let y = 0; y < yBound; ++y) {
                for (let x = 0; x < xBound; ++x) {
                    const v = input[x + xBound * (y + yBound * z)];
                    expect(v).toBeCloseTo(input[y + xBound * (x + yBound * z)], 10);
                    expect(v).toBeCloseTo(input[x + xBound * (z + yBound * y)], 10);
                }
            }
        }
    });

    it('matches per-axis FastGaussianBlur1 passes on random input', () => {
        // The 3D update is output = center + logBase * (xsum + ysum + zsum),
        // so it equals xBlur + yBlur + zBlur - 2 * input with each axis blur
        // computed by FastGaussianBlur1.
        const xBound = 7, yBound = 6, zBound = 5;
        const n = xBound * yBound * zBound;
        const scale = 1.4;
        const logBase = Math.log(1.25);
        const rand = makeRandom(0xb07);
        const input = Array.from({ length: n }, () => rand());

        const output = new Array<number>(n).fill(0);
        new FastGaussianBlur3().execute(xBound, yBound, zBound, input, output, scale, logBase);

        const blur1 = new FastGaussianBlur1();
        const axisBlur = (bound: number, stride: number, lineStart: (j: number) => number,
            numLines: number): number[] => {
            const result = new Array<number>(n).fill(0);
            for (let j = 0; j < numLines; ++j) {
                const start = lineStart(j);
                const line = new Array<number>(bound);
                for (let i = 0; i < bound; ++i) {
                    line[i] = input[start + stride * i];
                }
                const lineOut = new Array<number>(bound).fill(0);
                blur1.execute(bound, line, lineOut, scale, logBase);
                for (let i = 0; i < bound; ++i) {
                    result[start + stride * i] = lineOut[i];
                }
            }
            return result;
        };

        const xBlur = axisBlur(xBound, 1,
            (j) => xBound * j, yBound * zBound);
        const yBlur = axisBlur(yBound, xBound,
            (j) => (j % xBound) + xBound * yBound * Math.floor(j / xBound), xBound * zBound);
        const zBlur = axisBlur(zBound, xBound * yBound,
            (j) => j, xBound * yBound);

        for (let i = 0; i < n; ++i) {
            expect(output[i]).toBeCloseTo(xBlur[i] + yBlur[i] + zBlur[i] - 2 * input[i], 12);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V24): properties cross-checking the port against the
// upstream FastGaussianBlur3.h finite-difference stencil.
// ---------------------------------------------------------------------------

describe('FastGaussianBlur3 verification', () => {
    const bounds = fc.tuple(fc.integer({ min: 2, max: 5 }),
        fc.integer({ min: 2, max: 5 }), fc.integer({ min: 2, max: 5 }));
    const dyadicScale = fc.constantFrom(0.25, 0.5, 0.75, 1, 1.5);
    const dyadicLogBase = fc.constantFrom(0.0625, 0.125, 0.25);
    const voxels = (n: number) => fc.array(fc.integer({ min: -16, max: 16 }),
        { minLength: n, maxLength: n });

    function blur3(xBound: number, yBound: number, zBound: number,
        input: readonly number[], scale: number, logBase: number): number[] {
        const output = new Array<number>(xBound * yBound * zBound).fill(0);
        new FastGaussianBlur3().execute(xBound, yBound, zBound, input, output,
            scale, logBase);
        return output;
    }

    it('equals the three separable 1-dimensional passes minus twice the image', () => {
        // out3 = c + logBase * (xsum + ysum + zsum) and each axis pass gives
        // c + logBase * (that axis sum), so out3 = px + py + pz - 2 * input.
        check(fc.tuple(bounds, dyadicScale, dyadicLogBase),
            ([[xBound, yBound, zBound], scale, logBase]) => {
                const n = xBound * yBound * zBound;
                const index = (x: number, y: number, z: number) =>
                    x + xBound * (y + yBound * z);
                fc.assert(fc.property(voxels(n), input => {
                    const output = blur3(xBound, yBound, zBound, input, scale, logBase);
                    const blur1 = new FastGaussianBlur1();
                    const pass = (bound: number, gather: (i: number, a: number, b: number) => number,
                        scatter: (i: number, a: number, b: number, value: number) => void,
                        outerA: number, outerB: number) => {
                        for (let a = 0; a < outerA; ++a) {
                            for (let b = 0; b < outerB; ++b) {
                                const line = new Array<number>(bound);
                                for (let i = 0; i < bound; ++i) { line[i] = gather(i, a, b); }
                                const out = new Array<number>(bound).fill(0);
                                blur1.execute(bound, line, out, scale, logBase);
                                for (let i = 0; i < bound; ++i) { scatter(i, a, b, out[i]); }
                            }
                        }
                    };

                    const px = new Array<number>(n).fill(0);
                    pass(xBound, (i, y, z) => input[index(i, y, z)],
                        (i, y, z, v) => { px[index(i, y, z)] = v; }, yBound, zBound);
                    const py = new Array<number>(n).fill(0);
                    pass(yBound, (i, x, z) => input[index(x, i, z)],
                        (i, x, z, v) => { py[index(x, i, z)] = v; }, xBound, zBound);
                    const pz = new Array<number>(n).fill(0);
                    pass(zBound, (i, x, y) => input[index(x, y, i)],
                        (i, x, y, v) => { pz[index(x, y, i)] = v; }, xBound, yBound);

                    for (let i = 0; i < n; ++i) {
                        expectClose(output[i], px[i] + py[i] + pz[i] - 2 * input[i],
                            1e-12, 1e-12);
                    }
                }), { numRuns: 5 });
            }, 40);
    });

    it('is linear in the input image', () => {
        check(fc.tuple(bounds, dyadicScale, dyadicLogBase,
            fc.integer({ min: -3, max: 3 }), fc.integer({ min: -3, max: 3 })),
            ([[xBound, yBound, zBound], scale, logBase, a, b]) => {
                const n = xBound * yBound * zBound;
                fc.assert(fc.property(voxels(n), voxels(n), (u, v) => {
                    const combined = u.map((ui, i) => a * ui + b * v[i]);
                    const lhs = blur3(xBound, yBound, zBound, combined, scale, logBase);
                    const bu = blur3(xBound, yBound, zBound, u, scale, logBase);
                    const bv = blur3(xBound, yBound, zBound, v, scale, logBase);
                    for (let i = 0; i < n; ++i) {
                        expectClose(lhs[i], a * bu[i] + b * bv[i], 1e-12, 1e-12);
                    }
                }), { numRuns: 4 });
            }, 40);
    });

    it('leaves every constant image unchanged', () => {
        check(fc.tuple(bounds, dyadicScale, dyadicLogBase, finite(-100, 100)),
            ([[xBound, yBound, zBound], scale, logBase, c]) => {
                const input = new Array<number>(xBound * yBound * zBound).fill(c);
                for (const value of blur3(xBound, yBound, zBound, input, scale, logBase)) {
                    expectClose(value, c, 1e-12, 1e-12);
                }
            });
    });

    it('is invariant under a permutation of the axes', () => {
        // Relabeling (x, y, z) -> (y, z, x) must permute the result the same
        // way; this catches a swapped bound in i = x + xBound * (y + yBound * z).
        check(fc.tuple(bounds, dyadicScale, dyadicLogBase),
            ([[xBound, yBound, zBound], scale, logBase]) => {
                const n = xBound * yBound * zBound;
                fc.assert(fc.property(voxels(n), input => {
                    const output = blur3(xBound, yBound, zBound, input, scale, logBase);
                    const permuted = new Array<number>(n);
                    for (let z = 0; z < zBound; ++z) {
                        for (let y = 0; y < yBound; ++y) {
                            for (let x = 0; x < xBound; ++x) {
                                permuted[y + yBound * (z + zBound * x)] =
                                    input[x + xBound * (y + yBound * z)];
                            }
                        }
                    }
                    const outP = blur3(yBound, zBound, xBound, permuted, scale, logBase);
                    for (let z = 0; z < zBound; ++z) {
                        for (let y = 0; y < yBound; ++y) {
                            for (let x = 0; x < xBound; ++x) {
                                expectClose(outP[y + yBound * (z + zBound * x)],
                                    output[x + xBound * (y + yBound * z)], 1e-12, 1e-12);
                            }
                        }
                    }
                }), { numRuns: 5 });
            }, 40);
    });

    it('commutes with reflection along each axis', () => {
        check(fc.tuple(bounds, dyadicScale, dyadicLogBase),
            ([[xBound, yBound, zBound], scale, logBase]) => {
                const n = xBound * yBound * zBound;
                const index = (x: number, y: number, z: number) =>
                    x + xBound * (y + yBound * z);
                fc.assert(fc.property(voxels(n), input => {
                    const output = blur3(xBound, yBound, zBound, input, scale, logBase);
                    for (const axis of [0, 1, 2]) {
                        const flipped = new Array<number>(n);
                        for (let z = 0; z < zBound; ++z) {
                            for (let y = 0; y < yBound; ++y) {
                                for (let x = 0; x < xBound; ++x) {
                                    const sx = axis === 0 ? xBound - 1 - x : x;
                                    const sy = axis === 1 ? yBound - 1 - y : y;
                                    const sz = axis === 2 ? zBound - 1 - z : z;
                                    flipped[index(x, y, z)] = input[index(sx, sy, sz)];
                                }
                            }
                        }
                        const out = blur3(xBound, yBound, zBound, flipped, scale, logBase);
                        for (let z = 0; z < zBound; ++z) {
                            for (let y = 0; y < yBound; ++y) {
                                for (let x = 0; x < xBound; ++x) {
                                    const sx = axis === 0 ? xBound - 1 - x : x;
                                    const sy = axis === 1 ? yBound - 1 - y : y;
                                    const sz = axis === 2 ? zBound - 1 - z : z;
                                    expectClose(out[index(x, y, z)],
                                        output[index(sx, sy, sz)], 1e-12, 1e-12);
                                }
                            }
                        }
                    }
                }), { numRuns: 4 });
            }, 40);
    });
});
