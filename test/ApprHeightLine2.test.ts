import { describe, it, expect } from 'vitest';
import { ApprHeightLine2 } from '../src/ApprHeightLine2.js';
import { Vector } from '../src/Vector.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

describe('ApprHeightLine2', () => {
    it('initializes the parameters to zero', () => {
        const fitter = new ApprHeightLine2();
        const p = fitter.getParameters();
        expect(p.average.values).toEqual([0, 0]);
        expect(p.coefficients.values).toEqual([0, 0]);
        expect(fitter.getMinimumRequired()).toBe(2);
    });

    it('recovers a line from points that lie exactly on it', () => {
        // y = 2x + 1 sampled at x = 0..4.
        const points: Vector[] = [];
        for (let x = 0; x <= 4; ++x) {
            points.push(v2(x, 2 * x + 1));
        }
        const fitter = new ApprHeightLine2();
        expect(fitter.fit(points)).toBe(true);
        const p = fitter.getParameters();
        expect(p.average.values[0]).toBeCloseTo(2, 12);
        expect(p.average.values[1]).toBeCloseTo(5, 12);
        expect(p.coefficients.values[0]).toBeCloseTo(2, 12);
        expect(p.coefficients.values[1]).toBe(-1);

        for (const point of points) {
            expect(fitter.error(point)).toBeCloseTo(0, 20);
        }
    });

    it('recovers a horizontal line', () => {
        const points = [v2(-5, 7), v2(0, 7), v2(3, 7), v2(11, 7)];
        const fitter = new ApprHeightLine2();
        expect(fitter.fit(points)).toBe(true);
        expect(fitter.getParameters().coefficients.values[0]).toBeCloseTo(0, 12);
        expect(fitter.getParameters().average.values[1]).toBeCloseTo(7, 12);
    });

    it('matches the closed-form least-squares slope', () => {
        const random = makeRandom(4242);
        for (let trial = 0; trial < 20; ++trial) {
            const points: Vector[] = [];
            const n = 12;
            for (let i = 0; i < n; ++i) {
                const x = 10 * random() - 5;
                points.push(v2(x, 3 * x - 4 + (2 * random() - 1)));
            }
            const fitter = new ApprHeightLine2();
            expect(fitter.fit(points)).toBe(true);
            const p = fitter.getParameters();

            // Independent computation of the slope, sum(dx*dy)/sum(dx*dx).
            let sx = 0, sy = 0;
            for (const q of points) {
                sx += q.values[0];
                sy += q.values[1];
            }
            const xAvr = sx / n, yAvr = sy / n;
            let sxy = 0, sxx = 0;
            for (const q of points) {
                sxy += (q.values[0] - xAvr) * (q.values[1] - yAvr);
                sxx += (q.values[0] - xAvr) * (q.values[0] - xAvr);
            }
            expect(p.average.values[0]).toBeCloseTo(xAvr, 10);
            expect(p.average.values[1]).toBeCloseTo(yAvr, 10);
            expect(p.coefficients.values[0]).toBeCloseTo(sxy / sxx, 8);

            // The residuals of a least-squares fit sum to zero.
            let residual = 0;
            for (const q of points) {
                residual += (q.values[1] - yAvr)
                    - p.coefficients.values[0] * (q.values[0] - xAvr);
            }
            expect(residual).toBeCloseTo(0, 8);
        }
    });

    it('computes the squared vertical error', () => {
        const points = [v2(0, 0), v2(1, 2), v2(2, 4)];
        const fitter = new ApprHeightLine2();
        fitter.fit(points);
        // The fitted line is y = 2x; at (1,5) the vertical deviation is 3.
        expect(fitter.error(v2(1, 5))).toBeCloseTo(9, 10);
        expect(fitter.error(v2(3, 6))).toBeCloseTo(0, 10);
    });

    it('fails for points with a single x value', () => {
        const points = [v2(2, 0), v2(2, 1), v2(2, 5)];
        const fitter = new ApprHeightLine2();
        expect(fitter.fit(points)).toBe(false);
        expect(fitter.getParameters().average.values).toEqual([0, 0]);
        expect(fitter.getParameters().coefficients.values).toEqual([0, 0]);
    });

    it('fails for coincident points', () => {
        const points = [v2(1, 1), v2(1, 1), v2(1, 1)];
        const fitter = new ApprHeightLine2();
        expect(fitter.fit(points)).toBe(false);
        expect(fitter.getParameters().coefficients.values).toEqual([0, 0]);
    });

    it('fits an indexed subset and a contiguous range', () => {
        const points = [v2(0, 0), v2(1, 2), v2(2, 4), v2(3, 100)];
        const fitter = new ApprHeightLine2();
        expect(fitter.fit(points, [0, 1, 2])).toBe(true);
        expect(fitter.getParameters().coefficients.values[0]).toBeCloseTo(2, 12);
        expect(fitter.fit(points, 0, 2)).toBe(true);
        expect(fitter.getParameters().coefficients.values[0]).toBeCloseTo(2, 12);
    });

    it('deep-copies the parameters', () => {
        const source = new ApprHeightLine2();
        source.fit([v2(0, 0), v2(1, 2)]);
        const target = new ApprHeightLine2();
        target.copyParameters(source);
        expect(target.getParameters().coefficients.values[0]).toBeCloseTo(2, 12);
        source.getParameters().average.values[0] = 77;
        expect(target.getParameters().average.values[0]).toBeCloseTo(0.5, 12);
    });
});
