import { describe, expect, it } from 'vitest';
import { ApprPolynomialSpecial2 } from '../src/ApprPolynomialSpecial2.js';
import { ApprQuery } from '../src/ApprQuery.js';

// Solve the 2x2 system A*X = B by Cramer's rule; an independent check.
function solve2(A: number[][], B: number[]): number[] {
    const det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    return [
        (B[0] * A[1][1] - A[0][1] * B[1]) / det,
        (A[0][0] * B[1] - B[0] * A[1][0]) / det
    ];
}

describe('ApprPolynomialSpecial2', () => {
    it('validates the degrees in the constructor', () => {
        expect(() => new ApprPolynomialSpecial2([]))
            .toThrowError('The input array must have elements.');
        expect(() => new ApprPolynomialSpecial2([0, 2, 2]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial2([2, 1]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial2([-1, 3]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial2([0, 1, 5])).not.toThrow();
    });

    it('initializes the parameters to zero', () => {
        const fitter = new ApprPolynomialSpecial2([0, 2, 5]);
        expect(fitter.getParameters()).toEqual([0, 0, 0]);
        expect(fitter.getMinimumRequired()).toBe(3);
        expect(fitter.getXDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
    });

    it('reproduces data from w = 3 + 5*x^2 with the term set {1, x^2}', () => {
        // The internal transform maps x in [-1,1] to itself and w affinely,
        // so the model {1, x^2} represents the data exactly.
        const observations = [-1, -0.5, 0, 0.5, 1].map(
            (x) => [x, 3 + 5 * x * x]);

        const fitter = new ApprPolynomialSpecial2([0, 2]);
        expect(fitter.fit(observations)).toBe(true);
        expect(fitter.getXDomain()).toEqual([-1, 1]);

        // In the transformed space, w' = -1 + 2*(x')^2.
        const p = fitter.getParameters();
        expect(p[0]).toBeCloseTo(-1, 10);
        expect(p[1]).toBeCloseTo(2, 10);

        // evaluate() undoes the transform, so it reproduces the original
        // data and interpolates/extrapolates the underlying parabola.
        for (const [x, w] of observations) {
            expect(fitter.evaluate(x)).toBeCloseTo(w, 10);
            expect(fitter.error([x, w])).toBeLessThan(1e-10);
        }
        for (const x of [0.25, -0.8, 2, -3]) {
            expect(fitter.evaluate(x)).toBeCloseTo(3 + 5 * x * x, 8);
        }
    });

    it('restricts the fit to the requested powers', () => {
        // With the term set {1, x^3} an odd cubic through the origin is
        // recovered but a quadratic component cannot be.
        const observations = [-1, -0.5, 0, 0.5, 1].map(
            (x) => [x, 4 * x * x * x]);
        const fitter = new ApprPolynomialSpecial2([0, 3]);
        expect(fitter.fit(observations)).toBe(true);
        for (const [x, w] of observations) {
            expect(fitter.evaluate(x)).toBeCloseTo(w, 10);
        }
        expect(fitter.evaluate(0.75)).toBeCloseTo(4 * 0.75 ** 3, 8);
    });

    it('matches the transformed normal equations for the term set {1, x}', () => {
        const observations = [[0, 1], [1, 1.9], [2, 3.2], [3, 3.9], [4, 5.4]];

        // Reproduce the upstream transform to [-1,1]^2.
        const xs = observations.map((o) => o[0]);
        const ws = observations.map((o) => o[1]);
        const xmin = Math.min(...xs), xmax = Math.max(...xs);
        const wmin = Math.min(...ws), wmax = Math.max(...ws);
        const sx = 1 / (xmax - xmin), sw = 1 / (wmax - wmin);
        const tx = xs.map((x) => -1 + 2 * sx * (x - xmin));
        const tw = ws.map((w) => -1 + 2 * sw * (w - wmin));

        // A and B are the normalized sums over the terms (x^0, x^1).
        const n = observations.length;
        const A = [[0, 0], [0, 0]];
        const B = [0, 0];
        for (let i = 0; i < n; ++i) {
            for (let r = 0; r < 2; ++r) {
                for (let c = 0; c < 2; ++c) {
                    A[r][c] += Math.pow(tx[i], r + c);
                }
                B[r] += Math.pow(tx[i], r) * tw[i];
            }
        }
        for (let r = 0; r < 2; ++r) {
            B[r] /= n;
            for (let c = 0; c < 2; ++c) {
                A[r][c] /= n;
            }
        }
        const expected = solve2(A, B);

        const fitter = new ApprPolynomialSpecial2([0, 1]);
        expect(fitter.fit(observations)).toBe(true);
        const p = fitter.getParameters();
        expect(p[0]).toBeCloseTo(expected[0], 10);
        expect(p[1]).toBeCloseTo(expected[1], 10);

        // Because the transform is affine in both x and w, the {1, x} fit
        // agrees with the ordinary least-squares line in the original space.
        let s1 = 0, sxs = 0, sx2 = 0, sws = 0, sxw = 0;
        for (const [x, w] of observations) {
            s1 += 1; sxs += x; sx2 += x * x; sws += w; sxw += x * w;
        }
        const det = s1 * sx2 - sxs * sxs;
        const a0 = (sws * sx2 - sxs * sxw) / det;
        const a1 = (s1 * sxw - sws * sxs) / det;
        for (const [x] of observations) {
            expect(fitter.evaluate(x)).toBeCloseTo(a0 + a1 * x, 8);
        }
    });

    it('reports failure for invalid indices when validation is enabled', () => {
        const fitter = new ApprPolynomialSpecial2([0, 1]);
        ApprQuery.validateIndices = true;
        try {
            expect(fitter.fit([[0, 0], [1, 1]], [0, 9])).toBe(false);
            expect(fitter.getParameters()).toEqual([0, 0]);
        }
        finally {
            ApprQuery.validateIndices = false;
        }
    });

    it('copies the parameters between models', () => {
        const source = new ApprPolynomialSpecial2([0, 2]);
        source.fit([-1, -0.5, 0, 0.5, 1].map((x) => [x, 3 + 5 * x * x]));
        const target = new ApprPolynomialSpecial2([0, 2]);
        target.copyParameters(source);
        expect(target.getParameters()).toEqual(source.getParameters());
        expect(target.getXDomain()).toEqual(source.getXDomain());
        expect(target.evaluate(0.3)).toBe(source.evaluate(0.3));

        // The copy is deep, including the transform state.
        const copied = target.getParameters().slice();
        const value = target.evaluate(0.3);
        source.fit([-2, -1, 0, 1, 2].map((x) => [x, 1 - x * x]));
        expect(target.getParameters()).toEqual(copied);
        expect(target.evaluate(0.3)).toBe(value);
    });
});
