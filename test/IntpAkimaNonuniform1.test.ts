import { describe, it, expect } from 'vitest';
import { IntpAkimaNonuniform1 } from '../src/IntpAkimaNonuniform1.js';
import { IntpAkimaUniform1 } from '../src/IntpAkimaUniform1.js';

describe('IntpAkimaNonuniform1', () => {
    it('throws for invalid inputs', () => {
        expect(() => new IntpAkimaNonuniform1(2, [0, 1], [0, 1]))
            .toThrow('Invalid input to IntpAkima1 constructor.');
        // X is shorter than the sample count.
        expect(() => new IntpAkimaNonuniform1(4, [0, 1, 2], [0, 1, 2, 3]))
            .toThrow('Invalid input.');
        // X is not strictly increasing.
        expect(() => new IntpAkimaNonuniform1(3, [0, 1, 1], [0, 1, 2]))
            .toThrow('Invalid input.');
        expect(() => new IntpAkimaNonuniform1(3, [0, 2, 1], [0, 1, 2]))
            .toThrow('Invalid input.');
    });

    it('provides member access', () => {
        const X = [-1, 0.25, 3, 3.5, 10];
        const F = [0, 1, 4, 9, 16];
        const interp = new IntpAkimaNonuniform1(5, X, F);
        expect(interp.getQuantity()).toBe(5);
        expect(interp.getX()).toBe(X);
        expect(interp.getF()).toBe(F);
        expect(interp.getXMin()).toBe(-1);
        expect(interp.getXMax()).toBe(10);
    });

    it('interpolates the samples exactly', () => {
        const X = [-2, -1.5, 0, 0.25, 3, 7, 7.5];
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const interp = new IntpAkimaNonuniform1(X.length, X, F);
        for (let i = 0; i < X.length; ++i) {
            expect(interp.evaluate(X[i])).toBeCloseTo(F[i], 12);
            expect(interp.evaluate(0, X[i])).toBeCloseTo(F[i], 12);
        }
    });

    it('reproduces linear data exactly', () => {
        const X = [-2, -1.5, 0, 0.25, 3, 7, 7.5];
        const a = 1.5, b = -0.75;
        const F = X.map(x => a + b * x);
        const interp = new IntpAkimaNonuniform1(X.length, X, F);
        for (let k = 0; k <= 40; ++k) {
            const x = X[0] + (k / 40) * (X[X.length - 1] - X[0]);
            expect(interp.evaluate(x)).toBeCloseTo(a + b * x, 11);
            expect(interp.evaluate(1, x)).toBeCloseTo(b, 11);
            expect(interp.evaluate(2, x)).toBeCloseTo(0, 11);
            expect(interp.evaluate(3, x)).toBeCloseTo(0, 11);
        }
    });

    it('agrees with the uniform interpolator on uniformly spaced data', () => {
        const xMin = -1, dx = 0.5;
        const F = [1, -3, 2.5, 7, 0, -1.25, 4, 4, 6];
        const X: number[] = [];
        for (let i = 0; i < F.length; ++i) {
            X.push(xMin + i * dx);
        }
        const nonuniform = new IntpAkimaNonuniform1(F.length, X, F);
        const uniform = new IntpAkimaUniform1(F.length, xMin, dx, F);
        expect(nonuniform.getXMin()).toBeCloseTo(uniform.getXMin(), 14);
        expect(nonuniform.getXMax()).toBeCloseTo(uniform.getXMax(), 14);
        for (let k = 0; k <= 60; ++k) {
            const x = xMin + (k / 60) * (F.length - 1) * dx;
            for (let order = 0; order <= 3; ++order) {
                expect(nonuniform.evaluate(order, x))
                    .toBeCloseTo(uniform.evaluate(order, x), 10);
            }
        }
    });

    it('averages the two slopes when both slope pairs are equal', () => {
        // The x-values are nonuniform but the slopes are (1,1,1,2,2,2), so
        // at sample index 3 the estimated derivative is 0.5*(1 + 2) = 1.5.
        const X = [0, 1, 3, 4, 6, 6.5, 9];
        const slopes = [1, 1, 1, 2, 2, 2];
        const F = [0];
        for (let i = 0; i < slopes.length; ++i) {
            F.push(F[i] + slopes[i] * (X[i + 1] - X[i]));
        }
        const interp = new IntpAkimaNonuniform1(X.length, X, F);
        expect(interp.evaluate(1, X[3])).toBeCloseTo(1.5, 12);
        expect(interp.evaluate(1, X[1])).toBeCloseTo(1, 12);
        expect(interp.evaluate(1, X[5])).toBeCloseTo(2, 12);
    });

    it('is C1 continuous at the interior knots', () => {
        const X = [-2, -1.5, 0, 0.25, 3, 7, 7.5];
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const interp = new IntpAkimaNonuniform1(X.length, X, F);
        const h = 1e-9;
        for (let i = 1; i < X.length - 1; ++i) {
            expect(interp.evaluate(1, X[i] - h))
                .toBeCloseTo(interp.evaluate(1, X[i] + h), 6);
            expect(interp.evaluate(X[i] - h))
                .toBeCloseTo(interp.evaluate(X[i] + h), 6);
        }
    });

    it('matches finite differences of the function', () => {
        const X = [-2, -1.5, 0, 0.25, 3, 7, 7.5];
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const interp = new IntpAkimaNonuniform1(X.length, X, F);
        const h = 1e-5;
        for (const x of [-1.9, -0.4, 0.1, 1.7, 5.5, 7.2]) {
            const d1 = (interp.evaluate(x + h) - interp.evaluate(x - h)) / (2 * h);
            expect(interp.evaluate(1, x)).toBeCloseTo(d1, 5);
        }
    });

    it('clamps inputs outside the domain', () => {
        const X = [-2, -1.5, 0, 0.25, 3];
        const F = [1, -3, 2.5, 7, 0];
        const interp = new IntpAkimaNonuniform1(X.length, X, F);
        expect(interp.evaluate(-100)).toBeCloseTo(F[0], 12);
        expect(interp.evaluate(100)).toBeCloseTo(F[F.length - 1], 12);
        expect(interp.evaluate(1, -100)).toBe(interp.evaluate(1, X[0]));
        expect(interp.evaluate(1, 100)).toBe(interp.evaluate(1, X[X.length - 1]));
    });

    it('returns zero for derivative orders larger than three', () => {
        const X = [0, 1, 3, 4, 6];
        const F = [1, -3, 2.5, 7, 0];
        const interp = new IntpAkimaNonuniform1(5, X, F);
        expect(interp.evaluate(4, 2)).toBe(0);
    });

    it('handles the minimum sample count', () => {
        const X = [0, 1.5, 4];
        const F = [0, 2, 1];
        const interp = new IntpAkimaNonuniform1(3, X, F);
        for (let i = 0; i < 3; ++i) {
            expect(interp.evaluate(X[i])).toBeCloseTo(F[i], 12);
        }
    });

    it('is invariant under an affine change of the x-variable', () => {
        // Akima interpolation of {(x[i], f[i])} and {(a*x[i]+b, f[i])} give
        // the same function values at corresponding points.
        const X = [-2, -1.5, 0, 0.25, 3, 7, 7.5];
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const a = 2.5, b = -1;
        const Y = X.map(x => a * x + b);
        const interp0 = new IntpAkimaNonuniform1(X.length, X, F);
        const interp1 = new IntpAkimaNonuniform1(Y.length, Y, F);
        for (let k = 0; k <= 50; ++k) {
            const x = X[0] + (k / 50) * (X[X.length - 1] - X[0]);
            expect(interp1.evaluate(a * x + b)).toBeCloseTo(interp0.evaluate(x), 10);
            expect(a * interp1.evaluate(1, a * x + b))
                .toBeCloseTo(interp0.evaluate(1, x), 10);
        }
    });
});
