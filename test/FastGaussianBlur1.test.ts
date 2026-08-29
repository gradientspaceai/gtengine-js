import { describe, it, expect } from 'vitest';
import { FastGaussianBlur1 } from '../src/FastGaussianBlur1';

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
