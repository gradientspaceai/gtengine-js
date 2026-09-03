import { describe, it, expect } from 'vitest';
import { OdeRungeKutta4 } from '../src/OdeRungeKutta4.js';
import { OdeSolver } from '../src/OdeSolver.js';
import { Vector } from '../src/Vector.js';

describe('OdeRungeKutta4', () => {
    it('is an OdeSolver whose step size is readable and writable', () => {
        const solver = new OdeRungeKutta4(0.25, (_t, x) => x);
        expect(solver).toBeInstanceOf(OdeSolver);
        expect(solver.getTDelta()).toBe(0.25);
        solver.setTDelta(0.5);
        expect(solver.getTDelta()).toBe(0.5);
    });

    it('takes the exact classical RK4 step for dx/dt = x', () => {
        // The four stages produce the truncated exponential series
        // x*(1 + h + h^2/2 + h^3/6 + h^4/24).
        const h = 0.25;
        const solver = new OdeRungeKutta4(h, (_t, x) => x);
        const { tOut, xOut } = solver.update(0, Vector.fromArray([1]));
        const expected = 1 + h + h * h / 2 + h * h * h / 6 + h * h * h * h / 24;
        expect(tOut).toBe(h);
        expect(xOut.get(0)).toBeCloseTo(expected, 15);
        // The truncated series underestimates exp(h).
        expect(xOut.get(0)).toBeLessThan(Math.exp(h));
    });

    it('evaluates F four times at the Runge-Kutta nodes', () => {
        const times: number[] = [];
        const solver = new OdeRungeKutta4(0.4, (t, _x) => {
            times.push(t);
            return Vector.fromArray([0]);
        });
        const { tOut } = solver.update(1, Vector.fromArray([7]));
        expect(times.length).toBe(4);
        expect(times[0]).toBe(1);
        expect(times[1]).toBeCloseTo(1.2, 15);
        expect(times[2]).toBeCloseTo(1.2, 15);
        expect(times[3]).toBeCloseTo(1.4, 15);
        expect(tOut).toBeCloseTo(1.4, 15);
    });

    it('integrates polynomials of degree <= 3 exactly (Simpson rule)', () => {
        // x' = 4*t^3 with x(0) = 0 has the solution t^4. The RK4 quadrature
        // for a right-hand side independent of x is Simpson's rule, which is
        // exact for cubics.
        const solver = new OdeRungeKutta4(2, (t, _x) =>
            Vector.fromArray([4 * t * t * t]));
        const { xOut } = solver.update(0, Vector.fromArray([0]));
        expect(xOut.get(0)).toBeCloseTo(16, 10);
    });

    it('does not modify the input vector', () => {
        const solver = new OdeRungeKutta4(0.5, (_t, x) => x);
        const xIn = Vector.fromArray([1, 2, 3]);
        const { xOut } = solver.update(0, xIn);
        expect(xIn.values).toEqual([1, 2, 3]);
        expect(xOut).not.toBe(xIn);
    });

    it('solves the harmonic oscillator with fourth-order accuracy', () => {
        const F = (_t: number, x: Vector) =>
            Vector.fromArray([x.get(1), -x.get(0)]);

        const errorForSteps = (numSteps: number): number => {
            const tFinal = 1;
            const h = tFinal / numSteps;
            const solver = new OdeRungeKutta4(h, F);
            let t = 0;
            let x = Vector.fromArray([1, 0]);
            for (let i = 0; i < numSteps; ++i) {
                const result = solver.update(t, x);
                t = result.tOut;
                x = result.xOut;
            }
            return Math.max(
                Math.abs(x.get(0) - Math.cos(tFinal)),
                Math.abs(x.get(1) + Math.sin(tFinal)));
        };

        const e0 = errorForSteps(8);
        const e1 = errorForSteps(16);
        // The ratio of errors is approximately 2^4 for a fourth-order method.
        expect(e0 / e1).toBeGreaterThan(13);
        expect(e0 / e1).toBeLessThan(19);
    });

    it('is much more accurate than the exact solution requires for e', () => {
        // Integrate x' = x from 0 to 1 with 10 steps; the result approximates
        // e with error on the order of h^4.
        const numSteps = 10;
        const solver = new OdeRungeKutta4(1 / numSteps, (_t, x) => x);
        let t = 0;
        let x = Vector.fromArray([1]);
        for (let i = 0; i < numSteps; ++i) {
            const result = solver.update(t, x);
            t = result.tOut;
            x = result.xOut;
        }
        expect(t).toBeCloseTo(1, 12);
        expect(x.get(0)).toBeCloseTo(Math.E, 5);
    });
});
