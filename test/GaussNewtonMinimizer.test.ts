import { describe, it, expect } from 'vitest';
import { GaussNewtonMinimizer } from '../src/GaussNewtonMinimizer.js';
import { Matrix } from '../src/Matrix.js';
import { Vector } from '../src/Vector.js';

// A deterministic pseudorandom generator so the randomized cross-checks are
// reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state + 0x6D2B79F5) >>> 0;
        let t = state;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

// Solve the 2x2 system [a b; c d] x = (e, f) directly.
function solve2(a: number, b: number, c: number, d: number, e: number,
    f: number): [number, number] {
    const det = a * d - b * c;
    return [(e * d - b * f) / det, (a * f - e * c) / det];
}

// The linear least-squares problem F(p) = A*p - b, where A is numF-by-numP.
// The Gauss-Newton normal equations are exact for this F, so one step from
// any starting point lands on the least-squares minimizer.
function makeLinearProblem(A: Matrix, b: Vector) {
    const numF = A.numRows, numP = A.numCols;
    const fFunction = (p: Vector, f: Vector): void => {
        for (let r = 0; r < numF; ++r) {
            let sum = -b.get(r);
            for (let c = 0; c < numP; ++c) {
                sum += A.get(r, c) * p.get(c);
            }
            f.set(r, sum);
        }
    };
    const jFunction = (_p: Vector, j: Matrix): void => {
        for (let r = 0; r < numF; ++r) {
            for (let c = 0; c < numP; ++c) {
                j.set(r, c, A.get(r, c));
            }
        }
    };
    return { fFunction, jFunction, numF, numP };
}

// The model y(t) = p0 * exp(p1 * t), fit to samples of an exact instance of
// the model, so the global minimum has zero error at the true parameters.
const exponentialTrue = [2, -0.5];
const exponentialT = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5];
const exponentialY = exponentialT.map(
    t => exponentialTrue[0] * Math.exp(exponentialTrue[1] * t));

function exponentialF(p: Vector, f: Vector): void {
    for (let i = 0; i < exponentialT.length; ++i) {
        f.set(i, p.get(0) * Math.exp(p.get(1) * exponentialT[i])
            - exponentialY[i]);
    }
}

function exponentialJ(p: Vector, j: Matrix): void {
    for (let i = 0; i < exponentialT.length; ++i) {
        const e = Math.exp(p.get(1) * exponentialT[i]);
        j.set(i, 0, e);
        j.set(i, 1, p.get(0) * exponentialT[i] * e);
    }
}

// The equivalent "J plus" callback: J^T*J and -J^T*F computed directly.
function exponentialJPlus(p: Vector, jtj: Matrix, negJTF: Vector): void {
    const j = new Matrix(exponentialT.length, 2);
    exponentialJ(p, j);
    const f = new Vector(exponentialT.length);
    exponentialF(p, f);
    for (let r = 0; r < 2; ++r) {
        for (let c = 0; c < 2; ++c) {
            let sum = 0;
            for (let i = 0; i < exponentialT.length; ++i) {
                sum += j.get(i, r) * j.get(i, c);
            }
            jtj.set(r, c, sum);
        }
        let sum = 0;
        for (let i = 0; i < exponentialT.length; ++i) {
            sum += j.get(i, r) * f.get(i);
        }
        negJTF.set(r, -sum);
    }
}

describe('GaussNewtonMinimizer construction', () => {
    it('reports the dimensions and rejects invalid ones', () => {
        const minimizer = GaussNewtonMinimizer.fromJFunction(2, 5,
            () => { }, () => { });
        expect(minimizer.getNumPDimensions()).toBe(2);
        expect(minimizer.getNumFDimensions()).toBe(5);

        expect(() => GaussNewtonMinimizer.fromJFunction(0, 5,
            () => { }, () => { })).toThrow('Invalid dimensions.');
        expect(() => GaussNewtonMinimizer.fromJFunction(2, 0,
            () => { }, () => { })).toThrow('Invalid dimensions.');
        expect(() => GaussNewtonMinimizer.fromJPlusFunction(-1, 5,
            () => { }, () => { })).toThrow('Invalid dimensions.');
    });
});

describe('GaussNewtonMinimizer on linear least squares', () => {
    it('lands on the normal-equations solution in a single step', () => {
        // Fit the line y = p0 + p1*t to the samples (0,1), (1,2), (2,2),
        // (3,5). The least-squares solution comes from the normal equations
        // [n sum(t); sum(t) sum(t^2)] p = [sum(y); sum(t*y)].
        const t = [0, 1, 2, 3];
        const y = [1, 2, 2, 5];
        const A = Matrix.fromArray(4, 2, [
            1, t[0],
            1, t[1],
            1, t[2],
            1, t[3]
        ]);
        const b = Vector.fromArray(y);
        const sumT = t.reduce((a, v) => a + v, 0);
        const sumTT = t.reduce((a, v) => a + v * v, 0);
        const sumY = y.reduce((a, v) => a + v, 0);
        const sumTY = t.reduce((a, v, i) => a + v * y[i], 0);
        const expected = solve2(t.length, sumT, sumT, sumTT, sumY, sumTY);

        const { fFunction, jFunction } = makeLinearProblem(A, b);
        const minimizer = GaussNewtonMinimizer.fromJFunction(2, 4, fFunction,
            jFunction);

        // The starting point is irrelevant for a linear model.
        for (const p0 of [Vector.fromArray([0, 0]),
        Vector.fromArray([10, -7])]) {
            const result = minimizer.minimize(p0, 8, 1e-12, 1e-12);
            expect(result.minLocation.get(0)).toBeCloseTo(expected[0], 10);
            expect(result.minLocation.get(1)).toBeCloseTo(expected[1], 10);

            // The minimum error is the sum of squared residuals at the
            // least-squares solution.
            let expectedError = 0;
            for (let i = 0; i < t.length; ++i) {
                const r = expected[0] + expected[1] * t[i] - y[i];
                expectedError += r * r;
            }
            expect(result.minError).toBeCloseTo(expectedError, 10);
        }
    });

    it('solves a consistent linear system exactly', () => {
        // A square, invertible A with b = A*pTrue, so the minimum error is 0.
        const A = Matrix.fromArray(2, 2, [3, 1, 1, 2]);
        const pTrue = Vector.fromArray([-1, 4]);
        const b = Vector.fromArray([
            3 * pTrue.get(0) + 1 * pTrue.get(1),
            1 * pTrue.get(0) + 2 * pTrue.get(1)
        ]);
        const { fFunction, jFunction } = makeLinearProblem(A, b);
        const minimizer = GaussNewtonMinimizer.fromJFunction(2, 2, fFunction,
            jFunction);
        const result = minimizer.minimize(Vector.fromArray([5, 5]), 10, 1e-14,
            1e-14);
        expect(result.minLocation.get(0)).toBeCloseTo(-1, 10);
        expect(result.minLocation.get(1)).toBeCloseTo(4, 10);
        expect(result.minError).toBeLessThan(1e-20);
        // Note: 'converged' stays false here. For an exactly linear model
        // the second step reproduces the same location, so the error does
        // not strictly decrease and the tolerance tests - which upstream
        // performs only inside the strict-decrease branch - are never
        // reached.
        expect(result.converged).toBe(false);
    });
});

describe('GaussNewtonMinimizer on a nonlinear model', () => {
    it('recovers the parameters of an exponential model', () => {
        const minimizer = GaussNewtonMinimizer.fromJFunction(2,
            exponentialT.length, exponentialF, exponentialJ);
        const result = minimizer.minimize(Vector.fromArray([1.5, -0.3]), 64,
            1e-14, 1e-16);
        expect(result.minLocation.get(0)).toBeCloseTo(exponentialTrue[0], 6);
        expect(result.minLocation.get(1)).toBeCloseTo(exponentialTrue[1], 6);
        expect(result.minError).toBeLessThan(1e-14);
        expect(result.converged).toBe(true);
        expect(result.numIterations).toBeLessThanOrEqual(64);
    });

    it('gives the same result through the J and J-plus callbacks', () => {
        const rng = makeRng(1234567);
        for (let trial = 0; trial < 20; ++trial) {
            const p0 = Vector.fromArray([
                exponentialTrue[0] + 0.5 * (2 * rng() - 1),
                exponentialTrue[1] + 0.2 * (2 * rng() - 1)
            ]);
            const viaJ = GaussNewtonMinimizer.fromJFunction(2,
                exponentialT.length, exponentialF, exponentialJ)
                .minimize(p0, 32, 1e-14, 1e-16);
            const viaJPlus = GaussNewtonMinimizer.fromJPlusFunction(2,
                exponentialT.length, exponentialF, exponentialJPlus)
                .minimize(p0, 32, 1e-14, 1e-16);

            expect(viaJPlus.converged).toBe(viaJ.converged);
            expect(viaJPlus.numIterations).toBe(viaJ.numIterations);
            expect(viaJPlus.minLocation.get(0))
                .toBeCloseTo(viaJ.minLocation.get(0), 9);
            expect(viaJPlus.minLocation.get(1))
                .toBeCloseTo(viaJ.minLocation.get(1), 9);
            expect(viaJPlus.minError).toBeCloseTo(viaJ.minError, 12);
        }
    });
});

describe('GaussNewtonMinimizer degenerate behavior', () => {
    it('stops when the Cholesky factorization fails', () => {
        // F(p) = p0 + p1 - 1 has the rank-1 Jacobian [1 1], so J^T*J is
        // singular and cannot be factored.
        const fFunction = (p: Vector, f: Vector): void => {
            f.set(0, p.get(0) + p.get(1) - 1);
        };
        const jFunction = (_p: Vector, j: Matrix): void => {
            j.set(0, 0, 1);
            j.set(0, 1, 1);
        };
        const minimizer = GaussNewtonMinimizer.fromJFunction(2, 1, fFunction,
            jFunction);
        const p0 = Vector.fromArray([0, 0]);
        const result = minimizer.minimize(p0, 10, 1e-12, 1e-12);
        expect(result.converged).toBe(false);
        expect(result.numIterations).toBe(1);
        expect(result.minLocation.values).toEqual([0, 0]);
        expect(result.minError).toBeCloseTo(1, 12);
    });

    it('does no iterations when maxIterations is zero', () => {
        const minimizer = GaussNewtonMinimizer.fromJFunction(2,
            exponentialT.length, exponentialF, exponentialJ);
        const p0 = Vector.fromArray([1.5, -0.3]);
        const result = minimizer.minimize(p0, 0, 1e-12, 1e-12);
        expect(result.converged).toBe(false);
        expect(result.minLocation.values).toEqual(p0.values);
        // The initial error is |F(p0)|^2.
        const f = new Vector(exponentialT.length);
        exponentialF(p0, f);
        let error = 0;
        for (let i = 0; i < f.size; ++i) {
            error += f.get(i) * f.get(i);
        }
        expect(result.minError).toBeCloseTo(error, 14);
        // Upstream reports the loop counter, which is one more than the
        // number of iterations performed when the loop is exhausted.
        expect(result.numIterations).toBe(1);
    });

    it('clamps negative tolerances to zero', () => {
        // With negative tolerances no convergence test can succeed, so the
        // minimizer runs all iterations; the result must still be the best
        // location found.
        const minimizer = GaussNewtonMinimizer.fromJFunction(2,
            exponentialT.length, exponentialF, exponentialJ);
        const result = minimizer.minimize(Vector.fromArray([1.5, -0.3]), 16,
            -1, -1);
        expect(result.minLocation.get(0)).toBeCloseTo(exponentialTrue[0], 6);
        expect(result.minLocation.get(1)).toBeCloseTo(exponentialTrue[1], 6);
        // The loop was exhausted, so the counter is maxIterations + 1.
        expect(result.numIterations).toBe(17);
    });

    it('does not mutate the caller supplied initial point', () => {
        const minimizer = GaussNewtonMinimizer.fromJFunction(2,
            exponentialT.length, exponentialF, exponentialJ);
        const p0 = Vector.fromArray([1.5, -0.3]);
        const result = minimizer.minimize(p0, 16, 1e-14, 1e-16);
        expect(p0.values).toEqual([1.5, -0.3]);
        expect(result.minLocation).not.toBe(p0);
    });
});
