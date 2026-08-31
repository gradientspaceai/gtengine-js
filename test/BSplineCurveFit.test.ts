import { describe, it, expect } from 'vitest';
import { BSplineCurveFit } from '../src/BSplineCurveFit';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Flat sample data of 'numSamples' points of the given dimension, generated
// by 'f' evaluated at the uniform parameters t = i/(numSamples-1).
function sample(dimension: number, numSamples: number,
    f: (t: number) => number[]): number[] {
    const data: number[] = [];
    for (let i = 0; i < numSamples; ++i) {
        const point = f(i / (numSamples - 1));
        for (let j = 0; j < dimension; ++j) {
            data.push(point[j]);
        }
    }
    return data;
}

describe('BSplineCurveFit construction', () => {
    it('validates its preconditions', () => {
        const data = sample(2, 20, (t) => [t, 2 * t]);
        expect(() => new BSplineCurveFit(0, 20, data, 3, 8)).toThrow();
        // degree must satisfy 1 <= degree < numControls.
        expect(() => new BSplineCurveFit(2, 20, data, 0, 8)).toThrow();
        expect(() => new BSplineCurveFit(2, 20, data, 8, 8)).toThrow();
        // numControls <= numSamples - degree - 1.
        expect(() => new BSplineCurveFit(2, 20, data, 3, 17)).toThrow();
        expect(() => new BSplineCurveFit(2, 20, [1, 2], 3, 8)).toThrow();
        // The boundary case is allowed.
        expect(() => new BSplineCurveFit(2, 20, data, 3, 16)).not.toThrow();
    });

    it('reports the input and output information', () => {
        const data = sample(3, 25, (t) => [t, t * t, 1]);
        const fit = new BSplineCurveFit(3, 25, data, 3, 10);
        expect(fit.getDimension()).toBe(3);
        expect(fit.getNumSamples()).toBe(25);
        expect(fit.getDegree()).toBe(3);
        expect(fit.getNumControls()).toBe(10);
        expect(fit.getSampleData()).toBe(data);
        expect(fit.getControlData().length).toBe(30);
        // The fitted basis is open and uniform with numControls control
        // points.
        const basis = fit.getBasis();
        expect(basis.getDegree()).toBe(3);
        expect(basis.getNumControls()).toBe(10);
        expect(basis.isOpen()).toBe(true);
        expect(basis.isUniform()).toBe(true);
        expect(basis.getMinDomain()).toBe(0);
        expect(basis.getMaxDomain()).toBe(1);
    });
});

describe('BSplineCurveFit fitting', () => {
    it('matches the first and last samples exactly', () => {
        const data = sample(2, 30, (t) => [Math.cos(3 * t), Math.sin(3 * t)]);
        const fit = new BSplineCurveFit(2, 30, data, 3, 10);
        const controls = fit.getControlData();
        expect(controls[0]).toBe(data[0]);
        expect(controls[1]).toBe(data[1]);
        expect(controls[2 * 9]).toBe(data[2 * 29]);
        expect(controls[2 * 9 + 1]).toBe(data[2 * 29 + 1]);
        // An open spline interpolates the end control points.
        expect(fit.getPosition(0)).toEqual([data[0], data[1]]);
        expect(fit.getPosition(1)).toEqual([data[58], data[59]]);
    });

    it('recovers constant data exactly', () => {
        const data = sample(3, 40, () => [2, -3, 0.5]);
        const fit = new BSplineCurveFit(3, 40, data, 3, 12);
        for (let i = 0; i < 12; ++i) {
            expect(fit.getControlData()[3 * i]).toBeCloseTo(2, 10);
            expect(fit.getControlData()[3 * i + 1]).toBeCloseTo(-3, 10);
            expect(fit.getControlData()[3 * i + 2]).toBeCloseTo(0.5, 10);
        }
        for (const t of [0, 0.13, 0.5, 0.87, 1]) {
            const position = fit.getPosition(t);
            expect(position[0]).toBeCloseTo(2, 10);
            expect(position[1]).toBeCloseTo(-3, 10);
            expect(position[2]).toBeCloseTo(0.5, 10);
        }
    });

    it('recovers a straight line, which lies in the spline space', () => {
        const p0 = [1, -2, 3], p1 = [4, 6, -5];
        const numSamples = 50;
        const data = sample(3, numSamples,
            (t) => [0, 1, 2].map((k) => p0[k] + t * (p1[k] - p0[k])));
        const fit = new BSplineCurveFit(3, numSamples, data, 3, 12);
        for (const t of [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1]) {
            const position = fit.getPosition(t);
            for (let k = 0; k < 3; ++k) {
                expect(position[k]).toBeCloseTo(p0[k] + t * (p1[k] - p0[k]), 8);
            }
        }
    });

    it('recovers a cubic polynomial, which lies in the degree-3 spline space', () => {
        const f = (t: number): number[] =>
            [t * t * t - 2 * t + 1, 3 * t * t - t, -t * t * t + 5];
        const numSamples = 60;
        const data = sample(3, numSamples, f);
        const fit = new BSplineCurveFit(3, numSamples, data, 3, 15);
        for (const t of [0.05, 0.2, 0.4, 0.61, 0.8, 0.95]) {
            const position = fit.getPosition(t);
            const expected = f(t);
            for (let k = 0; k < 3; ++k) {
                expect(position[k]).toBeCloseTo(expected[k], 7);
            }
        }
    });

    it('approximates noisy data with a small residual', () => {
        const rand = makeRandom(99);
        const truth = (t: number): number[] => [Math.cos(6 * t), Math.sin(6 * t)];
        const numSamples = 200;
        const data: number[] = [];
        for (let i = 0; i < numSamples; ++i) {
            const t = i / (numSamples - 1);
            const point = truth(t);
            data.push(point[0] + 0.01 * (2 * rand() - 1));
            data.push(point[1] + 0.01 * (2 * rand() - 1));
        }
        const fit = new BSplineCurveFit(2, numSamples, data, 3, 20);

        // The fit smooths the noise: the residual against the noise-free
        // curve is far smaller than the noise amplitude at interior samples.
        let maxError = 0;
        for (let i = 5; i < numSamples - 5; ++i) {
            const t = i / (numSamples - 1);
            const position = fit.getPosition(t);
            const point = truth(t);
            maxError = Math.max(maxError,
                Math.abs(position[0] - point[0]), Math.abs(position[1] - point[1]));
        }
        expect(maxError).toBeLessThan(0.01);
    });

    it('evaluates derivatives consistently with finite differences', () => {
        const data = sample(2, 40, (t) => [t, Math.sin(4 * t)]);
        const fit = new BSplineCurveFit(2, 40, data, 3, 12);
        const h = 1e-6;
        for (const t of [0.2, 0.5, 0.8]) {
            const derivative = fit.evaluate(t, 1);
            const plus = fit.getPosition(t + h);
            const minus = fit.getPosition(t - h);
            for (let k = 0; k < 2; ++k) {
                expect(derivative[k]).toBeCloseTo((plus[k] - minus[k]) / (2 * h), 4);
            }
        }
    });

    it('clamps parameters outside [0,1]', () => {
        const data = sample(2, 30, (t) => [t, t * t]);
        const fit = new BSplineCurveFit(2, 30, data, 3, 10);
        expect(fit.getPosition(-5)).toEqual(fit.getPosition(0));
        expect(fit.getPosition(5)).toEqual(fit.getPosition(1));
    });
});
