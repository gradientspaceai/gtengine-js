import { describe, it, expect } from 'vitest';
import {
    LevenbergMarquardtMinimizer
} from '../src/LevenbergMarquardtMinimizer.js';
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

// The model y(t) = p0 * exp(p1 * t), fit to samples of an exact instance of
// the model, so the global minimum has zero error at the true parameters.
const expTrue = [2, -0.5];
const expT = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5];
const expY = expT.map(t => expTrue[0] * Math.exp(expTrue[1] * t));

function expF(p: Vector, f: Vector): void {
    for (let i = 0; i < expT.length; ++i) {
        f.set(i, p.get(0) * Math.exp(p.get(1) * expT[i]) - expY[i]);
    }
}

function expJ(p: Vector, j: Matrix): void {
    for (let i = 0; i < expT.length; ++i) {
        const e = Math.exp(p.get(1) * expT[i]);
        j.set(i, 0, e);
        j.set(i, 1, p.get(0) * expT[i] * e);
    }
}

function expJPlus(p: Vector, jtj: Matrix, negJTF: Vector): void {
    const j = new Matrix(expT.length, 2);
    expJ(p, j);
    const f = new Vector(expT.length);
    expF(p, f);
    for (let r = 0; r < 2; ++r) {
        for (let c = 0; c < 2; ++c) {
            let sum = 0;
            for (let i = 0; i < expT.length; ++i) {
                sum += j.get(i, r) * j.get(i, c);
            }
            jtj.set(r, c, sum);
        }
        let sum = 0;
        for (let i = 0; i < expT.length; ++i) {
            sum += j.get(i, r) * f.get(i);
        }
        negJTF.set(r, -sum);
    }
}

// F(x) = atan(x) is the classic example for which the undamped Gauss-Newton
// (Newton) step overshoots and increases the error when |x| is large enough.
// It forces the Levenberg-Marquardt inner loop to increase lambda several
// times before an update reduces the error, which is exactly the situation
// in which the upstream stale-residual bug shows up.
function atanF(p: Vector, f: Vector): void {
    f.set(0, Math.atan(p.get(0)));
}

function atanJ(p: Vector, j: Matrix): void {
    const x = p.get(0);
    j.set(0, 0, 1 / (1 + x * x));
}

function atanJPlus(p: Vector, jtj: Matrix, negJTF: Vector): void {
    const x = p.get(0);
    const d = 1 / (1 + x * x);
    jtj.set(0, 0, d * d);
    negJTF.set(0, -d * Math.atan(x));
}

describe('LevenbergMarquardtMinimizer construction', () => {
    it('reports the dimensions and rejects invalid ones', () => {
        const minimizer = LevenbergMarquardtMinimizer.fromJFunction(3, 7,
            () => { }, () => { });
        expect(minimizer.getNumPDimensions()).toBe(3);
        expect(minimizer.getNumFDimensions()).toBe(7);

        expect(() => LevenbergMarquardtMinimizer.fromJFunction(0, 5,
            () => { }, () => { })).toThrow('Invalid dimensions.');
        expect(() => LevenbergMarquardtMinimizer.fromJPlusFunction(2, 0,
            () => { }, () => { })).toThrow('Invalid dimensions.');
    });
});

describe('LevenbergMarquardtMinimizer on least-squares problems', () => {
    it('recovers the parameters of an exponential model', () => {
        const minimizer = LevenbergMarquardtMinimizer.fromJFunction(2,
            expT.length, expF, expJ);
        const result = minimizer.minimize(Vector.fromArray([1.0, -0.1]), 128,
            1e-12, 1e-14, 0.001, 10, 8);
        expect(result.minLocation.get(0)).toBeCloseTo(expTrue[0], 5);
        expect(result.minLocation.get(1)).toBeCloseTo(expTrue[1], 5);
        expect(result.minError).toBeLessThan(1e-10);
    });

    it('finds the least-squares line through noisy samples', () => {
        // Fit y = p0 + p1*t. The Jacobian is constant, so the minimum is the
        // normal-equations solution.
        const t = [0, 1, 2, 3, 4];
        const y = [1.1, 1.9, 3.2, 3.8, 5.1];
        const n = t.length;
        const sumT = t.reduce((a, v) => a + v, 0);
        const sumTT = t.reduce((a, v) => a + v * v, 0);
        const sumY = y.reduce((a, v) => a + v, 0);
        const sumTY = t.reduce((a, v, i) => a + v * y[i], 0);
        const det = n * sumTT - sumT * sumT;
        const expected = [
            (sumY * sumTT - sumT * sumTY) / det,
            (n * sumTY - sumY * sumT) / det
        ];

        const fFunction = (p: Vector, f: Vector): void => {
            for (let i = 0; i < n; ++i) {
                f.set(i, p.get(0) + p.get(1) * t[i] - y[i]);
            }
        };
        const jFunction = (_p: Vector, j: Matrix): void => {
            for (let i = 0; i < n; ++i) {
                j.set(i, 0, 1);
                j.set(i, 1, t[i]);
            }
        };

        const minimizer = LevenbergMarquardtMinimizer.fromJFunction(2, n,
            fFunction, jFunction);
        const result = minimizer.minimize(Vector.fromArray([-5, 5]), 200,
            1e-14, 1e-16, 0.001, 10, 8);
        expect(result.minLocation.get(0)).toBeCloseTo(expected[0], 6);
        expect(result.minLocation.get(1)).toBeCloseTo(expected[1], 6);
    });

    it('succeeds where Gauss-Newton fails, because lambda regularizes a '
        + 'singular normal-equations matrix', () => {
            // F(p) = p0 + p1 - 1 has the rank-1 Jacobian [1 1], so J^T*J is
            // singular. Adding lambda*average(diagonal) makes it positive
            // definite, so the Levenberg-Marquardt step is computable.
            const fFunction = (p: Vector, f: Vector): void => {
                f.set(0, p.get(0) + p.get(1) - 1);
            };
            const jFunction = (_p: Vector, j: Matrix): void => {
                j.set(0, 0, 1);
                j.set(0, 1, 1);
            };
            const minimizer = LevenbergMarquardtMinimizer.fromJFunction(2, 1,
                fFunction, jFunction);
            const result = minimizer.minimize(Vector.fromArray([0, 0]), 200,
                1e-14, 1e-16, 0.001, 10, 8);
            // Any point on the line p0 + p1 = 1 is a global minimizer.
            expect(result.minError).toBeLessThan(1e-12);
            expect(result.minLocation.get(0) + result.minLocation.get(1))
                .toBeCloseTo(1, 6);
        });
});

describe('LevenbergMarquardtMinimizer lambda adjustment', () => {
    it('increases lambda until an update reduces the error', () => {
        // Count the callbacks so the inner adjustment loop is observably
        // exercised: an adjustment repeats the iteration at the same
        // pCurrent, which requires re-evaluating F there (see the
        // upstream-bug note in the port).
        let numF = 0, numJ = 0;
        const countingF = (p: Vector, f: Vector): void => {
            ++numF;
            atanF(p, f);
        };
        const countingJ = (p: Vector, j: Matrix): void => {
            ++numJ;
            atanJ(p, j);
        };

        const minimizer = LevenbergMarquardtMinimizer.fromJFunction(1, 1,
            countingF, countingJ);
        const result = minimizer.minimize(Vector.fromArray([2]), 64, 1e-14,
            1e-16, 0.001, 10, 12);

        // Each call of the iteration evaluates J once and F once (at the
        // candidate); the extra F evaluations are the re-evaluations at a
        // repeated pCurrent, so their presence proves lambda was adjusted.
        expect(numF).toBeGreaterThan(numJ + 1);
        // The global minimum of |atan(x)|^2 is at x = 0.
        expect(Math.abs(result.minLocation.get(0))).toBeLessThan(1e-3);
        expect(result.minError).toBeLessThan(1e-6);
    });

    it('gives the same result through the J and J-plus callbacks', () => {
        // The J-plus callback computes J^T*J and -J^T*F from pCurrent alone,
        // so it is immune to the upstream stale-residual bug. Agreement of
        // the two paths on a problem that adjusts lambda is the regression
        // test for the fix.
        const viaJ = LevenbergMarquardtMinimizer.fromJFunction(1, 1, atanF,
            atanJ).minimize(Vector.fromArray([2]), 64, 1e-14, 1e-16, 0.001,
                10, 12);
        const viaJPlus = LevenbergMarquardtMinimizer.fromJPlusFunction(1, 1,
            atanF, atanJPlus).minimize(Vector.fromArray([2]), 64, 1e-14,
                1e-16, 0.001, 10, 12);
        expect(viaJPlus.minLocation.get(0))
            .toBeCloseTo(viaJ.minLocation.get(0), 9);
        expect(viaJPlus.minError).toBeCloseTo(viaJ.minError, 12);

        // The same agreement must hold for the intermediate iterates, which
        // is where the stale-residual bug is visible: with the upstream
        // code, the first iterate of the J path is 1.369... instead of
        // -0.7678....
        for (let iterations = 1; iterations <= 6; ++iterations) {
            const a = LevenbergMarquardtMinimizer.fromJFunction(1, 1, atanF,
                atanJ).minimize(Vector.fromArray([2]), iterations, 1e-14,
                    1e-16, 0.001, 10, 12);
            const b = LevenbergMarquardtMinimizer.fromJPlusFunction(1, 1,
                atanF, atanJPlus).minimize(Vector.fromArray([2]), iterations,
                    1e-14, 1e-16, 0.001, 10, 12);
            expect(b.minLocation.get(0))
                .toBeCloseTo(a.minLocation.get(0), 12);
            expect(b.minError).toBeCloseTo(a.minError, 14);
        }

        // The same agreement on the exponential fit, from several starts.
        const rng = makeRng(24680);
        for (let trial = 0; trial < 12; ++trial) {
            const p0 = Vector.fromArray([
                expTrue[0] + 1.0 * (2 * rng() - 1),
                expTrue[1] + 0.4 * (2 * rng() - 1)
            ]);
            const a = LevenbergMarquardtMinimizer.fromJFunction(2, expT.length,
                expF, expJ).minimize(p0, 64, 1e-12, 1e-14, 0.001, 10, 8);
            const b = LevenbergMarquardtMinimizer.fromJPlusFunction(2,
                expT.length, expF, expJPlus)
                .minimize(p0, 64, 1e-12, 1e-14, 0.001, 10, 8);
            expect(b.minLocation.get(0)).toBeCloseTo(a.minLocation.get(0), 7);
            expect(b.minLocation.get(1)).toBeCloseTo(a.minLocation.get(1), 7);
            expect(b.minError).toBeCloseTo(a.minError, 12);
        }
    });

    it('falls back to Gauss-Newton when the lambda inputs are invalid',
        () => {
            // Nonpositive lambdaFactor or lambdaAdjust forces
            // maxAdjustments = 1 and lambda = 0, which is a plain
            // Gauss-Newton iteration. The exponential problem is well
            // behaved from a nearby start, so it still converges.
            const p0 = Vector.fromArray([1.8, -0.45]);
            const damped = LevenbergMarquardtMinimizer.fromJFunction(2,
                expT.length, expF, expJ)
                .minimize(p0, 64, 1e-12, 1e-14, 0, 10, 8);
            const alsoDamped = LevenbergMarquardtMinimizer.fromJFunction(2,
                expT.length, expF, expJ)
                .minimize(p0, 64, 1e-12, 1e-14, 0.001, -1, 8);
            for (const result of [damped, alsoDamped]) {
                expect(result.minLocation.get(0)).toBeCloseTo(expTrue[0], 6);
                expect(result.minLocation.get(1)).toBeCloseTo(expTrue[1], 6);
                expect(result.minError).toBeLessThan(1e-12);
            }
        });
});

describe('LevenbergMarquardtMinimizer degenerate behavior', () => {
    it('stops when the Cholesky factorization fails', () => {
        // A zero Jacobian makes J^T*J zero, and no lambda adjustment can
        // help because the diagonal adjustment is proportional to the
        // (zero) diagonal sum.
        const fFunction = (p: Vector, f: Vector): void => {
            f.set(0, p.get(0) - 3);
        };
        const jFunction = (_p: Vector, j: Matrix): void => {
            j.set(0, 0, 0);
        };
        const minimizer = LevenbergMarquardtMinimizer.fromJFunction(1, 1,
            fFunction, jFunction);
        const result = minimizer.minimize(Vector.fromArray([0]), 10, 1e-12,
            1e-12, 0.001, 10, 4);
        expect(result.converged).toBe(false);
        expect(result.numIterations).toBe(1);
        expect(result.minLocation.values).toEqual([0]);
        expect(result.minError).toBeCloseTo(9, 12);
    });

    it('does no iterations when maxIterations is zero', () => {
        const minimizer = LevenbergMarquardtMinimizer.fromJFunction(2,
            expT.length, expF, expJ);
        const p0 = Vector.fromArray([1.5, -0.3]);
        const result = minimizer.minimize(p0, 0, 1e-12, 1e-12, 0.001, 10, 8);
        expect(result.converged).toBe(false);
        expect(result.minLocation.values).toEqual(p0.values);
        const f = new Vector(expT.length);
        expF(p0, f);
        let error = 0;
        for (let i = 0; i < f.size; ++i) {
            error += f.get(i) * f.get(i);
        }
        expect(result.minError).toBeCloseTo(error, 14);
        // Upstream reports the loop counter, which is one more than the
        // number of iterations performed when the loop is exhausted.
        expect(result.numIterations).toBe(1);
    });

    it('does not mutate the caller supplied initial point', () => {
        const minimizer = LevenbergMarquardtMinimizer.fromJFunction(2,
            expT.length, expF, expJ);
        const p0 = Vector.fromArray([1.5, -0.3]);
        const result = minimizer.minimize(p0, 32, 1e-12, 1e-14, 0.001, 10, 8);
        expect(p0.values).toEqual([1.5, -0.3]);
        expect(result.minLocation).not.toBe(p0);
    });
});
