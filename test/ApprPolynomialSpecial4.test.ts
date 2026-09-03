import { describe, expect, it } from 'vitest';
import { ApprPolynomialSpecial4 } from '../src/ApprPolynomialSpecial4.js';
import { ApprQuery } from '../src/ApprQuery.js';

const grid = [-1, -0.5, 0, 0.5, 1];

function cube(f: (x: number, y: number, z: number) => number): number[][] {
    const observations: number[][] = [];
    for (const x of grid) {
        for (const y of grid) {
            for (const z of grid) {
                observations.push([x, y, z, f(x, y, z)]);
            }
        }
    }
    return observations;
}

describe('ApprPolynomialSpecial4', () => {
    it('validates the degrees in the constructor', () => {
        expect(() => new ApprPolynomialSpecial4([0, 1], [0, 1], [0]))
            .toThrowError('The input arrays must have the same size.');
        expect(() => new ApprPolynomialSpecial4([], [], []))
            .toThrowError('The input array must have elements.');
        expect(() => new ApprPolynomialSpecial4([0, 1], [0, 1], [1, 1]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial4([0, 2], [0, 1], [3, 1]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial4([0, 1, 2], [0, 1, 3], [0, 2, 4]))
            .not.toThrow();
    });

    it('initializes the parameters to zero', () => {
        const fitter = new ApprPolynomialSpecial4([0, 1], [0, 2], [0, 3]);
        expect(fitter.getParameters()).toEqual([0, 0]);
        expect(fitter.getMinimumRequired()).toBe(2);
        expect(fitter.getXDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
        expect(fitter.getYDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
        expect(fitter.getZDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
    });

    it('reproduces data from w = 7 + 3*x*y*z with the term set {1, x*y*z}', () => {
        const f = (x: number, y: number, z: number): number =>
            7 + 3 * x * y * z;
        const observations = cube(f);

        const fitter = new ApprPolynomialSpecial4([0, 1], [0, 1], [0, 1]);
        expect(fitter.fit(observations)).toBe(true);
        expect(fitter.getXDomain()).toEqual([-1, 1]);
        expect(fitter.getYDomain()).toEqual([-1, 1]);
        expect(fitter.getZDomain()).toEqual([-1, 1]);

        // With the identity transform on x, y and z, w' = x'*y'*z', so the
        // transformed coefficients are (0, 1).
        const p = fitter.getParameters();
        expect(p[0]).toBeCloseTo(0, 10);
        expect(p[1]).toBeCloseTo(1, 10);

        for (const obs of observations) {
            expect(fitter.evaluate(obs[0], obs[1], obs[2]))
                .toBeCloseTo(obs[3], 10);
            expect(fitter.error(obs)).toBeLessThan(1e-10);
        }
        for (const [x, y, z] of [[0.3, -0.4, 0.7], [2, 1, -3]]) {
            expect(fitter.evaluate(x, y, z)).toBeCloseTo(f(x, y, z), 8);
        }
    });

    it('reproduces data from the three-term set {1, x*y*z, x^2*y^2*z^2}', () => {
        const f = (x: number, y: number, z: number): number => {
            const t = x * y * z;
            return 2 - t + 5 * t * t;
        };
        const observations = cube(f);

        const fitter = new ApprPolynomialSpecial4([0, 1, 2], [0, 1, 2], [0, 1, 2]);
        expect(fitter.fit(observations)).toBe(true);
        for (const obs of observations) {
            expect(fitter.evaluate(obs[0], obs[1], obs[2]))
                .toBeCloseTo(obs[3], 8);
        }
        for (const [x, y, z] of [[0.3, -0.4, 0.9], [0.8, 0.8, 0.8]]) {
            expect(fitter.evaluate(x, y, z)).toBeCloseTo(f(x, y, z), 6);
        }
    });

    it('produces a genuine least-squares fit for data off the model', () => {
        const observations: number[][] = [];
        for (const [x, y, z] of cube(() => 0).map((o) => [o[0], o[1], o[2]])) {
            const w = 7 + 3 * x * y * z;
            observations.push([x, y, z, w + 0.6]);
            observations.push([x, y, z, w - 0.6]);
        }

        const fitter = new ApprPolynomialSpecial4([0, 1], [0, 1], [0, 1]);
        expect(fitter.fit(observations)).toBe(true);
        for (const x of grid) {
            for (const y of grid) {
                for (const z of grid) {
                    expect(fitter.evaluate(x, y, z))
                        .toBeCloseTo(7 + 3 * x * y * z, 8);
                }
            }
        }
    });

    it('reports failure for invalid indices when validation is enabled', () => {
        const fitter = new ApprPolynomialSpecial4([0, 1], [0, 1], [0, 1]);
        ApprQuery.validateIndices = true;
        try {
            expect(fitter.fit([[0, 0, 0, 0], [1, 1, 1, 1]], [0, 6])).toBe(false);
            expect(fitter.getParameters()).toEqual([0, 0]);
        }
        finally {
            ApprQuery.validateIndices = false;
        }
    });

    it('copies the parameters between models', () => {
        const observations = cube((x, y, z) => 7 + 3 * x * y * z);

        const source = new ApprPolynomialSpecial4([0, 1], [0, 1], [0, 1]);
        source.fit(observations);
        const target = new ApprPolynomialSpecial4([0, 1], [0, 1], [0, 1]);
        target.copyParameters(source);
        expect(target.getParameters()).toEqual(source.getParameters());
        expect(target.getXDomain()).toEqual(source.getXDomain());
        expect(target.getYDomain()).toEqual(source.getYDomain());
        expect(target.getZDomain()).toEqual(source.getZDomain());
        expect(target.evaluate(0.3, 0.4, 0.5))
            .toBe(source.evaluate(0.3, 0.4, 0.5));

        const copied = target.getParameters().slice();
        const value = target.evaluate(0.3, 0.4, 0.5);
        source.fit(cube((x, y, z) => 1 - x * y * z));
        expect(target.getParameters()).toEqual(copied);
        expect(target.evaluate(0.3, 0.4, 0.5)).toBe(value);
    });
});
