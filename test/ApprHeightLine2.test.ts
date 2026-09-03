import { describe, it, expect } from 'vitest';
import { ApprHeightLine2 } from '../src/ApprHeightLine2.js';
import { Vector } from '../src/Vector.js';
import { check, expectClose, expectVectorClose, fc, finite, vector } from './helpers/arbitraries.js';

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

describe('ApprHeightLine2 verification', () => {
    // Distinct abscissas so the fit is well posed.
    const xsArb = fc.array(finite(-10, 10), { minLength: 2, maxLength: 12 })
        .map(xs => [...new Set(xs)])
        .filter(xs => xs.length >= 2
            && Math.max(...xs) - Math.min(...xs) > 0.1);

    it('recovers the exact line of height data', () => {
        check(fc.tuple(xsArb, finite(-5, 5), finite(-5, 5)),
            ([xs, a, b]) => {
                const points = xs.map(
                    x => Vector.fromArray([x, a * x + b]));
                const fitter = new ApprHeightLine2();
                expect(fitter.fit(points)).toBe(true);
                const p = fitter.getParameters();

                // The parameters are ((xAvr,yAvr),(a,-1)).
                expect(p.coefficients.get(1)).toBe(-1);
                expectClose(p.coefficients.get(0), a, 1e-8, 1e-8);

                const xAvr = xs.reduce((u, v) => u + v, 0) / xs.length;
                expectClose(p.average.get(0), xAvr, 1e-9, 1e-9);
                expectClose(p.average.get(1), a * xAvr + b, 1e-8, 1e-8);

                // The error of every sample vanishes.
                for (const q of points) {
                    expectClose(fitter.error(q), 0, 1e-14, 0);
                }
            });
    });

    it('satisfies the least-squares normal equation', () => {
        // The minimizer of sum_i [a*(x_i-xAvr) - (y_i-yAvr)]^2 satisfies
        // sum_i [a*dx_i - dy_i] * dx_i = 0.
        check(fc.tuple(xsArb, fc.array(finite(-10, 10),
            { minLength: 2, maxLength: 12 })), ([xs, ys]) => {
                const n = Math.min(xs.length, ys.length);
                if (n < 2) { return; }
                const points = Array.from({ length: n },
                    (_, i) => Vector.fromArray([xs[i], ys[i]]));
                const fitter = new ApprHeightLine2();
                if (!fitter.fit(points)) { return; }
                const p = fitter.getParameters();
                const a = p.coefficients.get(0);

                let residual = 0, scale = 0;
                for (const q of points) {
                    const dx = q.get(0) - p.average.get(0);
                    const dy = q.get(1) - p.average.get(1);
                    residual += (a * dx - dy) * dx;
                    scale += Math.abs(a * dx * dx) + Math.abs(dy * dx);
                }
                expect(Math.abs(residual))
                    .toBeLessThanOrEqual(1e-9 + 1e-9 * scale);
            });
    });

    it('is equivariant under a translation of the samples', () => {
        check(fc.tuple(xsArb, finite(-5, 5), finite(-5, 5), vector(2, -20, 20)),
            ([xs, a, b, t]) => {
                const points = xs.map(x => Vector.fromArray([x, a * x + b]));
                const shifted = points.map(q => Vector.fromArray(
                    [q.get(0) + t.get(0), q.get(1) + t.get(1)]));

                const f0 = new ApprHeightLine2();
                const f1 = new ApprHeightLine2();
                expect(f0.fit(points)).toBe(true);
                expect(f1.fit(shifted)).toBe(true);

                // The slope is unchanged and the average translates.
                expectClose(f1.getParameters().coefficients.get(0),
                    f0.getParameters().coefficients.get(0), 1e-8, 1e-8);
                expectClose(f1.getParameters().average.get(0),
                    f0.getParameters().average.get(0) + t.get(0), 1e-8, 1e-8);
                expectClose(f1.getParameters().average.get(1),
                    f0.getParameters().average.get(1) + t.get(1), 1e-8, 1e-8);
            });
    });

    it('fails and zeroes the parameters for vertical or empty data', () => {
        // covar00 == 0 is the documented failure condition: every sample has
        // the same abscissa. The abscissa is an integer and the sample count
        // a power of two so that the average is exactly x -- upstream
        // divides by multiplying with 1/n, which is inexact for other counts
        // and leaves a covariance of the order of one ulp.
        check(fc.tuple(fc.integer({ min: -8, max: 8 }),
            fc.constantFrom(1, 2, 4, 8),
            fc.array(finite(-10, 10), { minLength: 8, maxLength: 8 })),
            ([x, n, ys]) => {
                const points = ys.slice(0, n).map(
                    y => Vector.fromArray([x, y]));
                const fitter = new ApprHeightLine2();
                expect(fitter.fit(points)).toBe(false);
                const p = fitter.getParameters();
                expect(p.average.values).toEqual([0, 0]);
                expect(p.coefficients.values).toEqual([0, 0]);
            });

        const empty = new ApprHeightLine2();
        expect(empty.fit([])).toBe(false);
        expect(empty.getParameters().average.values).toEqual([0, 0]);
    });

    it('copyParameters deep-copies the parameters', () => {
        check(fc.tuple(xsArb, finite(-5, 5), finite(-5, 5)), ([xs, a, b]) => {
            const source = new ApprHeightLine2();
            expect(source.fit(xs.map(
                x => Vector.fromArray([x, a * x + b])))).toBe(true);
            const target = new ApprHeightLine2();
            target.copyParameters(source);

            const s = source.getParameters();
            const t = target.getParameters();
            expect(t.average).not.toBe(s.average);
            expect(t.coefficients).not.toBe(s.coefficients);
            expectVectorClose(t.average, s.average, 0, 0);
            expectVectorClose(t.coefficients, s.coefficients, 0, 0);

            s.average.set(0, 4321);
            expect(t.average.get(0)).not.toBe(4321);
        });
    });
});
