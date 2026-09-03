import { describe, expect, it } from 'vitest';
import { ApprPolynomialSpecial3 } from '../src/ApprPolynomialSpecial3.js';
import { ApprQuery } from '../src/ApprQuery.js';

const grid = [-1, -0.5, 0, 0.5, 1];

describe('ApprPolynomialSpecial3', () => {
    it('validates the degrees in the constructor', () => {
        expect(() => new ApprPolynomialSpecial3([0, 1], [0]))
            .toThrowError('The input arrays must have the same size.');
        expect(() => new ApprPolynomialSpecial3([], []))
            .toThrowError('The input array must have elements.');
        expect(() => new ApprPolynomialSpecial3([0, 0], [0, 1]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial3([0, 1], [1, 0]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial3([0, 1, 2], [0, 2, 3]))
            .not.toThrow();
    });

    it('rejects the natural term set {(0,0),(1,0),(0,1)}, as upstream does', () => {
        // Upstream requires each degree list to be strictly increasing,
        // which is stronger than "distinct <p,q> pairs"; the port preserves
        // the restriction.
        expect(() => new ApprPolynomialSpecial3([0, 1, 0], [0, 0, 1]))
            .toThrowError('Degrees must be increasing.');
    });

    it('initializes the parameters to zero', () => {
        const fitter = new ApprPolynomialSpecial3([0, 1], [0, 2]);
        expect(fitter.getParameters()).toEqual([0, 0]);
        expect(fitter.getMinimumRequired()).toBe(2);
        expect(fitter.getXDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
        expect(fitter.getYDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
    });

    it('reproduces data from w = 5 + 2*x*y with the term set {1, x*y}', () => {
        const observations: number[][] = [];
        for (const x of grid) {
            for (const y of grid) {
                observations.push([x, y, 5 + 2 * x * y]);
            }
        }

        const fitter = new ApprPolynomialSpecial3([0, 1], [0, 1]);
        expect(fitter.fit(observations)).toBe(true);
        expect(fitter.getXDomain()).toEqual([-1, 1]);
        expect(fitter.getYDomain()).toEqual([-1, 1]);

        // The transform is the identity on x and y here, and w' = x'*y', so
        // the transformed coefficients are (0, 1).
        const p = fitter.getParameters();
        expect(p[0]).toBeCloseTo(0, 10);
        expect(p[1]).toBeCloseTo(1, 10);

        for (const obs of observations) {
            expect(fitter.evaluate(obs[0], obs[1])).toBeCloseTo(obs[2], 10);
            expect(fitter.error(obs)).toBeLessThan(1e-10);
        }
        // Interpolation and extrapolation follow the same surface.
        for (const [x, y] of [[0.3, -0.4], [0.75, 0.75], [2, 3]]) {
            expect(fitter.evaluate(x, y)).toBeCloseTo(5 + 2 * x * y, 8);
        }
    });

    it('reproduces data from a three-term set {1, x*y, x^2*y^2}', () => {
        const observations: number[][] = [];
        for (const x of grid) {
            for (const y of grid) {
                const t = x * y;
                observations.push([x, y, 1 + 3 * t + 4 * t * t]);
            }
        }

        const fitter = new ApprPolynomialSpecial3([0, 1, 2], [0, 1, 2]);
        expect(fitter.fit(observations)).toBe(true);
        for (const obs of observations) {
            expect(fitter.evaluate(obs[0], obs[1])).toBeCloseTo(obs[2], 8);
        }
        for (const [x, y] of [[0.3, -0.4], [0.9, 0.6]]) {
            const t = x * y;
            expect(fitter.evaluate(x, y)).toBeCloseTo(1 + 3 * t + 4 * t * t, 6);
        }
    });

    it('produces a genuine least-squares fit for data off the model', () => {
        // w = 5 + 2*x*y plus a symmetric perturbation at each sample, so the
        // fit returns the unperturbed surface.
        const observations: number[][] = [];
        for (const x of grid) {
            for (const y of grid) {
                observations.push([x, y, 5 + 2 * x * y + 0.4]);
                observations.push([x, y, 5 + 2 * x * y - 0.4]);
            }
        }

        const fitter = new ApprPolynomialSpecial3([0, 1], [0, 1]);
        expect(fitter.fit(observations)).toBe(true);

        // The residuals are orthogonal to both basis terms in the
        // transformed space; equivalently, they cancel in pairs here.
        for (const x of grid) {
            for (const y of grid) {
                expect(fitter.evaluate(x, y)).toBeCloseTo(5 + 2 * x * y, 8);
            }
        }
    });

    it('reports failure for invalid indices when validation is enabled', () => {
        const fitter = new ApprPolynomialSpecial3([0, 1], [0, 1]);
        ApprQuery.validateIndices = true;
        try {
            expect(fitter.fit([[0, 0, 0], [1, 1, 1]], [0, 4])).toBe(false);
            expect(fitter.getParameters()).toEqual([0, 0]);
        }
        finally {
            ApprQuery.validateIndices = false;
        }
    });

    it('copies the parameters between models', () => {
        const observations: number[][] = [];
        for (const x of grid) {
            for (const y of grid) {
                observations.push([x, y, 5 + 2 * x * y]);
            }
        }

        const source = new ApprPolynomialSpecial3([0, 1], [0, 1]);
        source.fit(observations);
        const target = new ApprPolynomialSpecial3([0, 1], [0, 1]);
        target.copyParameters(source);
        expect(target.getParameters()).toEqual(source.getParameters());
        expect(target.getXDomain()).toEqual(source.getXDomain());
        expect(target.getYDomain()).toEqual(source.getYDomain());
        expect(target.evaluate(0.3, 0.4)).toBe(source.evaluate(0.3, 0.4));

        const copied = target.getParameters().slice();
        const value = target.evaluate(0.3, 0.4);
        source.fit(observations.map(([x, y]) => [x, y, 1 - x * y]));
        expect(target.getParameters()).toEqual(copied);
        expect(target.evaluate(0.3, 0.4)).toBe(value);
    });
});
