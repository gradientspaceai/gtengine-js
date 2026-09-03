import { describe, it, expect } from 'vitest';
import { FastGaussianBlur1 } from '../src/FastGaussianBlur1.js';
import { FastGaussianBlur3 } from '../src/FastGaussianBlur3.js';

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
