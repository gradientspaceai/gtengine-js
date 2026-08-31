import { describe, it, expect } from 'vitest';
import { OdeImplicitEuler } from '../src/OdeImplicitEuler';
import { OdeSolver } from '../src/OdeSolver';
import { Vector } from '../src/Vector';

describe('OdeImplicitEuler', () => {
    it('is an OdeSolver whose step size is readable and writable', () => {
        const solver = new OdeImplicitEuler(0.25, (_t, x) => x, () => [1]);
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
            () => [-k]);
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
        const solver = new OdeImplicitEuler(h, F, () => A);
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
            () => [0, 0, 0, 0]);
        const { xOut } = solver.update(1, Vector.fromArray([4, 5]));
        expect(xOut.get(0)).toBeCloseTo(4 + h * 2, 14);
        expect(xOut.get(1)).toBeCloseTo(5 + h * 2, 14);
    });

    it('does not modify the input vector', () => {
        const solver = new OdeImplicitEuler(0.5,
            (_t, x) => x, () => [1, 0, 0, 1]);
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
            () => [-k]);
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
        const DF = () => [-1];

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
            () => [1, 0, 0, 1]);
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
        const solver = new OdeImplicitEuler(h, F, () => A);
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
