import { describe, it, expect } from 'vitest';
import { FastGaussianBlur1 } from '../src/FastGaussianBlur1';
import { FastGaussianBlur2 } from '../src/FastGaussianBlur2';

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
