import { describe, it, expect } from 'vitest';
import { OdeImplicitEuler } from '../src/OdeImplicitEuler.js';
import { Matrix, determinant, lInfinityNorm, mulMatrix } from '../src/Matrix.js';
import { OdeSolver } from '../src/OdeSolver.js';
import { Vector } from '../src/Vector.js';
import { LinearSystem } from '../src/LinearSystem.js';
import {
    check, expectClose, fc, scaled, wellScaledMatrix, wellScaledVector
} from './helpers/arbitraries.js';

describe('OdeImplicitEuler', () => {
    it('is an OdeSolver whose step size is readable and writable', () => {
        const solver = new OdeImplicitEuler(0.25, (_t, x) => x, () => Matrix.fromArray(1, 1, [1]));
        expect(solver).toBeInstanceOf(OdeSolver);
        expect(solver.getTDelta()).toBe(0.25);
        solver.setTDelta(0.5);
        expect(solver.getTDelta()).toBe(0.5);
    });

    it('takes the exact backward-Euler step for dx/dt = -k*x', () => {
        // The update is x + h * (I - h*DF)^{-1} * F, which for F = -k*x and
        // DF = -k is x + h*(-k*x)/(1 + h*k) = x/(1 + h*k). That is the exact
        // solution of the implicit equation y = x + h*(-k*y).
        const k = 3;
        const h = 0.5;
        const solver = new OdeImplicitEuler(h,
            (_t, x) => Vector.fromArray([-k * x.get(0)]),
            () => Matrix.fromArray(1, 1, [-k]));
        const { tOut, xOut } = solver.update(0, Vector.fromArray([1]));
        expect(tOut).toBe(h);
        expect(xOut.get(0)).toBeCloseTo(1 / (1 + h * k), 14);
    });

    it('solves the implicit equation for a 2x2 linear system', () => {
        // x' = A*x with A = {{0,1},{-2,-3}}. The update solves
        // (I - h*A)*y = A*x for y, then x + h*y. Equivalently the returned
        // value is (I - h*A)^{-1} * x, the standard backward-Euler iterate.
        const h = 0.25;
        const A = [0, 1, -2, -3];
        const F = (_t: number, x: Vector) => Vector.fromArray([
            A[0] * x.get(0) + A[1] * x.get(1),
            A[2] * x.get(0) + A[3] * x.get(1)
        ]);
        const solver = new OdeImplicitEuler(h, F,
            () => Matrix.fromArray(2, 2, A));
        const xIn = Vector.fromArray([1, -1]);
        const { xOut } = solver.update(0, xIn);

        // Direct solve of (I - h*A)*y = x.
        const m = [1 - h * A[0], -h * A[1], -h * A[2], 1 - h * A[3]];
        const det = m[0] * m[3] - m[1] * m[2];
        const y0 = (m[3] * xIn.get(0) - m[1] * xIn.get(1)) / det;
        const y1 = (-m[2] * xIn.get(0) + m[0] * xIn.get(1)) / det;
        expect(xOut.get(0)).toBeCloseTo(y0, 12);
        expect(xOut.get(1)).toBeCloseTo(y1, 12);
    });

    it('reproduces forward Euler when the derivative matrix is zero', () => {
        // With DF = 0 the matrix I - h*DF is the identity, so the update
        // reduces to xOut = xIn + h*F.
        const h = 0.3;
        const solver = new OdeImplicitEuler(h,
            (t, _x) => Vector.fromArray([t + 1, 2]),
            () => Matrix.fromArray(2, 2, [0, 0, 0, 0]));
        const { xOut } = solver.update(1, Vector.fromArray([4, 5]));
        expect(xOut.get(0)).toBeCloseTo(4 + h * 2, 14);
        expect(xOut.get(1)).toBeCloseTo(5 + h * 2, 14);
    });

    it('does not modify the input vector', () => {
        const solver = new OdeImplicitEuler(0.5,
            (_t, x) => x, () => Matrix.fromArray(2, 2, [1, 0, 0, 1]));
        const xIn = Vector.fromArray([1, 2]);
        const { xOut } = solver.update(0, xIn);
        expect(xIn.values).toEqual([1, 2]);
        expect(xOut).not.toBe(xIn);
    });

    it('is unconditionally stable on a stiff decay where Euler diverges', () => {
        // x' = -100*x. Backward Euler with h = 0.1 damps to zero, whereas
        // forward Euler with the same step size has amplification factor
        // 1 - 100*h = -9 and diverges.
        const k = 100;
        const h = 0.1;
        const solver = new OdeImplicitEuler(h,
            (_t, x) => Vector.fromArray([-k * x.get(0)]),
            () => Matrix.fromArray(1, 1, [-k]));
        let t = 0;
        let x = Vector.fromArray([1]);
        for (let i = 0; i < 20; ++i) {
            const result = solver.update(t, x);
            t = result.tOut;
            x = result.xOut;
        }
        expect(Math.abs(x.get(0))).toBeLessThan(1e-10);
        // The iterate is (1 + h*k)^{-n}.
        expect(x.get(0)).toBeCloseTo(Math.pow(1 + h * k, -20), 20);
    });

    it('converges to the exact solution with first-order accuracy', () => {
        // x' = -x with x(0) = 1 has the solution exp(-t). Backward Euler is
        // first order, so halving h halves the error at t = 1.
        const F = (_t: number, x: Vector) => Vector.fromArray([-x.get(0)]);
        const DF = () => Matrix.fromArray(1, 1, [-1]);

        const errorForSteps = (numSteps: number): number => {
            const solver = new OdeImplicitEuler(1 / numSteps, F, DF);
            let t = 0;
            let x = Vector.fromArray([1]);
            for (let i = 0; i < numSteps; ++i) {
                const result = solver.update(t, x);
                t = result.tOut;
                x = result.xOut;
            }
            return Math.abs(x.get(0) - Math.exp(-1));
        };

        const e0 = errorForSteps(200);
        const e1 = errorForSteps(400);
        expect(e0 / e1).toBeGreaterThan(1.7);
        expect(e0 / e1).toBeLessThan(2.3);
    });

    it('produces the zero vector when I - h*DF is singular', () => {
        // Upstream Inverse(...) returns the zero matrix for a noninvertible
        // input, so the port does the same and the update degenerates to
        // xOut = xIn.
        const h = 1;
        const solver = new OdeImplicitEuler(h,
            (_t, x) => Vector.fromArray([x.get(0), x.get(1)]),
            () => Matrix.fromArray(2, 2, [1, 0, 0, 1]));
        const { xOut } = solver.update(0, Vector.fromArray([3, 4]));
        expect(xOut.get(0)).toBe(3);
        expect(xOut.get(1)).toBe(4);
    });

    it('handles a 3x3 nonsymmetric derivative matrix', () => {
        const h = 0.2;
        const A = [
            2, -1, 0,
            1, 3, -2,
            0, 1, 1
        ];
        const F = (_t: number, x: Vector) => {
            const r = new Vector(3);
            for (let i = 0; i < 3; ++i) {
                r.values[i] = A[0 + 3 * i] * x.get(0) + A[1 + 3 * i] * x.get(1) +
                    A[2 + 3 * i] * x.get(2);
            }
            return r;
        };
        const solver = new OdeImplicitEuler(h, F,
            () => Matrix.fromArray(3, 3, A));
        const xIn = Vector.fromArray([1, 2, -1]);
        const { xOut } = solver.update(0, xIn);

        // Verify (I - h*A) * (xOut - xIn)/h = A*xIn, the equation solved by
        // the update.
        const rhs = F(0, xIn);
        for (let r = 0; r < 3; ++r) {
            let sum = 0;
            for (let c = 0; c < 3; ++c) {
                const identity = (r === c ? 1 : 0);
                sum += (identity - h * A[c + 3 * r]) * (xOut.get(c) - xIn.get(c)) / h;
            }
            expect(sum).toBeCloseTo(rhs.get(r), 10);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V39): independent review against upstream OdeImplicitEuler.h.
// ---------------------------------------------------------------------------

describe('OdeImplicitEuler verification', () => {
    const rate = fc.oneof(scaled(-2, -0.5), scaled(0.5, 2));

    // A dimension, a derivative matrix and a step size. Upstream places no
    // condition on DF, but the update inverts I - h*DF, so draws whose
    // I - h*DF is nearly singular are skipped: the residual bounds below
    // would be meaningless there. The rejection test is scale relative, as
    // in invertibleMatrix.
    const derivativeSetup = fc.integer({ min: 1, max: 5 }).chain(n =>
        fc.tuple(fc.constant(n), wellScaledMatrix(n, n, -2, 2),
            wellScaledVector(n, -5, 5), scaled(0.05, 0.5)));

    function conditioned(n: number, DF: Matrix, h: number): Matrix {
        const dg = new Matrix(n, n);
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                dg.set(r, c, (r === c ? 1 : 0) - h * DF.get(r, c));
            }
        }
        const norm = lInfinityNorm(dg);
        fc.pre(Math.abs(determinant(dg)) >= 1e-2 * Math.pow(norm, n));
        return dg;
    }

    it('solves (I - h*DF) * (xOut - xIn) = h * F(xIn)', () => {
        // This is the linear system that defines the update, whatever F and
        // DF are; upstream never assumes DF is the actual derivative of F.
        check(fc.tuple(derivativeSetup, wellScaledVector(5, -5, 5)),
            ([[n, DF, xIn, h], Fvalues]) => {
                const dg = conditioned(n, DF, h);
                const F = new Vector(n);
                for (let i = 0; i < n; ++i) { F.set(i, Fvalues.get(i)); }
                const solver = new OdeImplicitEuler(h, () => F, () => DF);
                const { tOut, xOut } = solver.update(0, xIn);
                expect(tOut).toBe(h);
                for (let r = 0; r < n; ++r) {
                    let sum = 0;
                    for (let c = 0; c < n; ++c) {
                        sum += dg.get(r, c) * (xOut.get(c) - xIn.get(c));
                    }
                    // The residual of a solve is bounded by the round-off of
                    // the products it is assembled from, amplified by the
                    // condition number of I - h*DF.
                    expectClose(sum, h * F.get(r), 1e-9, 1e-9);
                }
            });
    });

    it('is the backward-Euler iterate (I - h*A)^{-1} for a linear field',
        () => {
            check(derivativeSetup, ([n, A, xIn, h]) => {
                // For F(t,x) = A*x the update simplifies to
                // xOut = (I - h*A)^{-1} * xIn, which is what the general
                // linear solver computes from the same matrix.
                const dg = conditioned(n, A, h);
                const solver = new OdeImplicitEuler(h,
                    (_t, x) => mulMatrix(A, x), () => A);
                const { xOut } = solver.update(0, xIn);
                const reference = LinearSystem.solve(n, dg.values,
                    xIn.values);
                expect(reference.invertible).toBe(true);
                for (let i = 0; i < n; ++i) {
                    expectClose(xOut.get(i), reference.X[i], 1e-9, 1e-9);
                }
            });
        });

    it('is first-order accurate on dx/dt = a*x', () => {
        check(rate, a => {
            const F = (_t: number, x: Vector) =>
                Vector.fromArray([a * x.get(0)]);
            const DF = () => Matrix.fromArray(1, 1, [a]);
            const errorForSteps = (numSteps: number): number => {
                const solver = new OdeImplicitEuler(1 / numSteps, F, DF);
                let t = 0;
                let x = Vector.fromArray([1]);
                for (let i = 0; i < numSteps; ++i) {
                    const r = solver.update(t, x);
                    t = r.tOut;
                    x = r.xOut;
                }
                return Math.abs(x.get(0) - Math.exp(a));
            };
            // The observed ratios over this range of a lie in [1.99, 2.04].
            const ratio = errorForSteps(64) / errorForSteps(128);
            expect(ratio).toBeGreaterThan(1.8);
            expect(ratio).toBeLessThan(2.2);
        });
    });

    it('is unconditionally stable on stiff decay', () => {
        check(fc.tuple(scaled(10, 1e6), scaled(0.01, 2)), ([k, h]) => {
            // dx/dt = -k*x with k*h arbitrarily large: the implicit iterate
            // is x/(1 + k*h), which never overshoots, whereas the explicit
            // Euler factor 1 - k*h diverges for k*h > 2.
            const solver = new OdeImplicitEuler(h,
                (_t, x) => Vector.fromArray([-k * x.get(0)]),
                () => Matrix.fromArray(1, 1, [-k]));
            let t = 0;
            let x = Vector.fromArray([1]);
            for (let i = 0; i < 16; ++i) {
                const r = solver.update(t, x);
                t = r.tOut;
                x = r.xOut;
                expect(x.get(0)).toBeGreaterThan(0);
                expect(x.get(0)).toBeLessThan(1);
            }
            // The update evaluates x - x*h*k/(1 + h*k) rather than
            // x/(1 + h*k), which cancels about log10(h*k) digits per step.
            expectClose(x.get(0), Math.pow(1 + k * h, -16), 0, 1e-6);
        });
    });

    it('leaves the state unchanged when I - h*DF is singular', () => {
        // Upstream Inverse() returns the zero matrix for a noninvertible
        // input, so the increment h * (dgInverse * F) is zero.
        check(fc.tuple(wellScaledVector(2, -5, 5),
            fc.constantFrom(2, 1, 0.5, 0.25, 0.125)),
        ([xIn, h]) => {
                // DF = I/h makes I - h*DF the zero matrix. The step sizes
                // are exact powers of two, so h*(1/h) is exactly 1 and the
                // matrix is exactly singular rather than merely tiny.
                const solver = new OdeImplicitEuler(h,
                    (_t, x) => x,
                    () => Matrix.fromArray(2, 2, [1 / h, 0, 0, 1 / h]));
                const { xOut } = solver.update(0, xIn);
                expect(xOut.get(0)).toBe(xIn.get(0));
                expect(xOut.get(1)).toBe(xIn.get(1));
            });
    });

    it('never aliases its input or its output', () => {
        check(fc.tuple(wellScaledVector(3, -5, 5), scaled(0.05, 2)),
            ([xIn, h]) => {
                const before = xIn.clone();
                const solver = new OdeImplicitEuler(h, (_t, x) => x,
                    () => new Matrix(3, 3));
                const { xOut } = solver.update(0, xIn);
                expect(xOut).not.toBe(xIn);
                for (let i = 0; i < 3; ++i) {
                    expect(xIn.get(i)).toBe(before.get(i));
                    // With DF = 0 the step is the explicit Euler step.
                    expect(xOut.get(i)).toBe(xIn.get(i) + h * xIn.get(i));
                }
            });
    });

    it('takes a Matrix from the derivative callback (V39 API alignment)',
        () => {
            // The callback used to return a row-major number[] and the file
            // carried a private Gauss-Jordan inverse, which is neither the
            // upstream algorithm (Inverse(GMatrix): GaussianElimination with
            // full pivoting, which is inverse() in Matrix.ts) nor the
            // convention the rest of the library uses.
            const solver = new OdeImplicitEuler(0.5,
                (_t, x) => x, () => Matrix.fromArray(1, 1, [2]));
            const { xOut } = solver.update(0, Vector.fromArray([1]));
            // I - 0.5*2 = 0 is singular, so the increment vanishes.
            expect(xOut.get(0)).toBe(1);
        });
});
