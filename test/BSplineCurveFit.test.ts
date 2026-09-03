import { describe, it, expect } from 'vitest';
import { check, fc, finite, expectClose, seededRandom } from './helpers/arbitraries.js';
import { BSplineCurveFit } from '../src/BSplineCurveFit.js';
import { BSplineCurve } from '../src/BSplineCurve.js';
import { BasisFunctionInput } from '../src/BasisFunction.js';
import { Vector } from '../src/Vector.js';

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

// ---------------------------------------------------------------------------
// Verification pass (V05). The fitted control points are cross-checked
// against the separately ported BSplineCurve class (which shares no code with
// BSplineCurveFit::Evaluate) and against the least-squares normal equations
// solved with an independent dense Gaussian elimination.
// ---------------------------------------------------------------------------
describe('BSplineCurveFit verification', () => {
    interface Case {
        dimension: number; numSamples: number; degree: number; numControls: number;
        samples: number[];
    }

    const caseArb: fc.Arbitrary<Case> = fc.tuple(
        fc.integer({ min: 1, max: 3 }),     // dimension
        fc.integer({ min: 1, max: 3 }),     // degree
        fc.integer({ min: 0, max: 4 }),     // numControls - degree - 1
        fc.integer({ min: 0, max: 8 }),     // extra samples
        fc.array(finite(-10, 10), { minLength: 3 * 40, maxLength: 3 * 40 })
    ).map(([dimension, degree, extraControls, extraSamples, data]) => {
        const numControls = degree + 1 + extraControls;
        const numSamples = numControls + degree + 1 + extraSamples;
        const samples: number[] = [];
        for (let i = 0; i < numSamples * dimension; ++i) {
            samples.push(data[i % data.length]);
        }
        return { dimension, numSamples, degree, numControls, samples };
    });

    const fitOf = (c: Case): BSplineCurveFit =>
        new BSplineCurveFit(c.dimension, c.numSamples, c.samples, c.degree, c.numControls);

    /** The fitted curve rebuilt as a BSplineCurve over the same open uniform basis. */
    function asCurve(fit: BSplineCurveFit, dimension: number): BSplineCurve {
        const controls: Vector[] = [];
        const data = fit.getControlData();
        for (let i = 0; i < fit.getNumControls(); ++i) {
            const p = new Vector(dimension);
            for (let k = 0; k < dimension; ++k) { p.values[k] = data[i * dimension + k]; }
            controls.push(p);
        }
        const input = new BasisFunctionInput(fit.getNumControls(), fit.getDegree());
        return new BSplineCurve(dimension, input, controls);
    }

    it('evaluates the same curve as BSplineCurve over the fitted controls', () => {
        check(fc.tuple(caseArb, finite(0, 1)), ([c, t]) => {
            const fit = fitOf(c);
            const curve = asCurve(fit, c.dimension);
            const got = fit.getPosition(t);
            const want = curve.getPosition(t);
            for (let k = 0; k < c.dimension; ++k) {
                expectClose(got[k], want.values[k], 1e-9, 1e-10);
            }
            return true;
        }, 60);
    });

    it('forces the end control points to the end samples', () => {
        check(caseArb, c => {
            const fit = fitOf(c);
            const data = fit.getControlData();
            expect(data.length).toBe(c.dimension * c.numControls);
            const lastControl = c.dimension * (c.numControls - 1);
            const lastSample = c.dimension * (c.numSamples - 1);
            for (let k = 0; k < c.dimension; ++k) {
                expect(data[k]).toBe(c.samples[k]);
                expect(data[lastControl + k]).toBe(c.samples[lastSample + k]);
            }
            // An open uniform B-spline interpolates its end control points,
            // so the fitted curve passes through the end samples.
            const p0 = fit.getPosition(0);
            const p1 = fit.getPosition(1);
            for (let k = 0; k < c.dimension; ++k) {
                expectClose(p0[k], c.samples[k], 1e-11, 1e-11);
                expectClose(p1[k], c.samples[lastSample + k], 1e-11, 1e-11);
            }
            return true;
        }, 60);
    });

    // Any polynomial of degree <= the spline degree lies in the spline space,
    // so the least-squares fit reproduces it. The end-point override is a
    // no-op in that case because the exact fit already interpolates the ends.
    it('reproduces polynomials that lie in the spline space', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 3 }),
            fc.array(finite(-5, 5), { minLength: 4, maxLength: 4 }),
            finite(0, 1)), ([degree, coeff, t]) => {
                const numControls = degree + 3;
                const numSamples = numControls + degree + 6;
                const poly = (x: number): number => {
                    let value = 0;
                    for (let k = degree; k >= 0; --k) { value = value * x + coeff[k]; }
                    return value;
                };
                const samples: number[] = [];
                for (let i = 0; i < numSamples; ++i) {
                    samples.push(poly(i / (numSamples - 1)));
                }
                const fit = new BSplineCurveFit(1, numSamples, samples, degree, numControls);
                // The normal equations are mildly ill conditioned, hence the
                // 1e-8 relative tolerance rather than machine precision.
                expectClose(fit.getPosition(t)[0], poly(t), 1e-8, 1e-8);
                return true;
            }, 80);
    });

    // The fit solves (A^T A) Q = A^T P. Rebuild that dense system from the
    // basis functions and solve it with plain Gaussian elimination, then
    // compare with the fitted controls (excluding the two overridden ends).
    it('solves the least-squares normal equations', () => {
        const degree = 2;
        const numControls = 6;
        const numSamples = 20;
        const rand = seededRandom(0x8f1c2);
        const samples: number[] = [];
        for (let i = 0; i < numSamples; ++i) { samples.push(2 * rand() - 1); }
        const fit = new BSplineCurveFit(1, numSamples, samples, degree, numControls);

        const basis = fit.getBasis();
        // A is numSamples x numControls with A[i][j] = N_j(t_i).
        const a: number[][] = [];
        for (let i = 0; i < numSamples; ++i) {
            const t = i / (numSamples - 1);
            const { minIndex, maxIndex } = basis.evaluate(t, 0);
            const row = new Array<number>(numControls).fill(0);
            for (let j = minIndex; j <= maxIndex; ++j) {
                row[j] = basis.getValue(0, j);
            }
            a.push(row);
        }
        // Normal equations.
        const m: number[][] = [];
        for (let r = 0; r < numControls; ++r) {
            const row = new Array<number>(numControls + 1).fill(0);
            for (let cIdx = 0; cIdx < numControls; ++cIdx) {
                let s = 0;
                for (let i = 0; i < numSamples; ++i) { s += a[i][r] * a[i][cIdx]; }
                row[cIdx] = s;
            }
            let rhs = 0;
            for (let i = 0; i < numSamples; ++i) { rhs += a[i][r] * samples[i]; }
            row[numControls] = rhs;
            m.push(row);
        }
        // Gaussian elimination with partial pivoting.
        for (let col = 0; col < numControls; ++col) {
            let pivot = col;
            for (let r = col + 1; r < numControls; ++r) {
                if (Math.abs(m[r][col]) > Math.abs(m[pivot][col])) { pivot = r; }
            }
            const tmp = m[col]; m[col] = m[pivot]; m[pivot] = tmp;
            const d = m[col][col];
            for (let cIdx = col; cIdx <= numControls; ++cIdx) { m[col][cIdx] /= d; }
            for (let r = 0; r < numControls; ++r) {
                if (r === col) { continue; }
                const f = m[r][col];
                if (f === 0) { continue; }
                for (let cIdx = col; cIdx <= numControls; ++cIdx) {
                    m[r][cIdx] -= f * m[col][cIdx];
                }
            }
        }

        const controls = fit.getControlData();
        // Skip the first and last, which the constructor overrides with the
        // end samples.
        for (let i = 1; i < numControls - 1; ++i) {
            expectClose(controls[i], m[i][numControls], 1e-9, 1e-9);
        }
        // And the ends really were overridden.
        expect(controls[0]).toBe(samples[0]);
        expect(controls[numControls - 1]).toBe(samples[numSamples - 1]);
    });

    it('rejects inputs that violate the documented preconditions', () => {
        const samples = new Array<number>(64).fill(0);
        expect(() => new BSplineCurveFit(0, 20, samples, 2, 6)).toThrow(/Invalid dimension/);
        expect(() => new BSplineCurveFit(1, 20, samples, 0, 6)).toThrow(/Invalid degree/);
        expect(() => new BSplineCurveFit(1, 20, samples, 6, 6)).toThrow(/Invalid degree/);
        expect(() => new BSplineCurveFit(1, 8, samples, 2, 6)).toThrow(/Invalid number of controls/);
        expect(() => new BSplineCurveFit(1, 20, [0, 1], 2, 6)).toThrow(/Invalid sample data/);
    });

    it('clamps parameters outside the domain of the open spline', () => {
        check(fc.tuple(caseArb, finite(1, 5)), ([c, over]) => {
            const fit = fitOf(c);
            const below = fit.getPosition(-over);
            const above = fit.getPosition(1 + over);
            for (let k = 0; k < c.dimension; ++k) {
                expectClose(below[k], fit.getPosition(0)[k], 1e-12, 1e-12);
                expectClose(above[k], fit.getPosition(1)[k], 1e-12, 1e-12);
            }
            return true;
        }, 40);
    });
});
