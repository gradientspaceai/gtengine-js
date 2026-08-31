import { describe, it, expect } from 'vitest';
import { OdeEuler } from '../src/OdeEuler';
import { OdeSolver } from '../src/OdeSolver';
import { Vector } from '../src/Vector';

describe('OdeEuler', () => {
    it('is an OdeSolver whose step size is readable and writable', () => {
        const solver = new OdeEuler(0.25, (_t, x) => x);
        expect(solver).toBeInstanceOf(OdeSolver);
        expect(solver.getTDelta()).toBe(0.25);
        solver.setTDelta(0.5);
        expect(solver.getTDelta()).toBe(0.5);
    });

    it('takes the exact Euler step for dx/dt = x', () => {
        // One step of x' = x from x(0) = 1 is x + h*x = 1 + h.
        const h = 0.125;
        const solver = new OdeEuler(h, (_t, x) => x);
        const { tOut, xOut } = solver.update(0, Vector.fromArray([1]));
        expect(tOut).toBe(h);
        expect(xOut.get(0)).toBe(1 + h);
    });

    it('evaluates F at (tIn, xIn) and advances t by tDelta', () => {
        const calls: Array<{ t: number; x: number }> = [];
        const solver = new OdeEuler(0.1, (t, x) => {
            calls.push({ t, x: x.get(0) });
            return Vector.fromArray([t * x.get(0)]);
        });
        const { tOut, xOut } = solver.update(2, Vector.fromArray([3]));
        expect(calls).toEqual([{ t: 2, x: 3 }]);
        expect(tOut).toBeCloseTo(2.1, 15);
        // x + h * (t * x) = 3 + 0.1 * 6.
        expect(xOut.get(0)).toBeCloseTo(3.6, 12);
    });

    it('does not modify the input vector', () => {
        const solver = new OdeEuler(0.5, (_t, x) => x);
        const xIn = Vector.fromArray([1, 2, 3]);
        const { xOut } = solver.update(0, xIn);
        expect(xIn.values).toEqual([1, 2, 3]);
        expect(xOut).not.toBe(xIn);
        expect(xOut.values).toEqual([1.5, 3, 4.5]);
    });

    it('solves the harmonic oscillator with first-order accuracy', () => {
        // x'' = -x written as the system (x0,x1)' = (x1,-x0) has the exact
        // solution (cos t, -sin t) for the initial value (1,0). The global
        // error of Euler's method is O(h), so halving h halves the error.
        const F = (_t: number, x: Vector) =>
            Vector.fromArray([x.get(1), -x.get(0)]);

        const errorForSteps = (numSteps: number): number => {
            const tFinal = 1;
            const h = tFinal / numSteps;
            const solver = new OdeEuler(h, F);
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

        const e0 = errorForSteps(200);
        const e1 = errorForSteps(400);
        expect(e1).toBeLessThan(e0);
        // The ratio of errors is approximately 2^1 for a first-order method.
        expect(e0 / e1).toBeGreaterThan(1.7);
        expect(e0 / e1).toBeLessThan(2.3);
    });

    it('reproduces the closed-form Euler iterate for dx/dt = x', () => {
        // n steps of size h produce (1 + h)^n exactly (up to rounding).
        const h = 0.01;
        const numSteps = 100;
        const solver = new OdeEuler(h, (_t, x) => x);
        let t = 0;
        let x = Vector.fromArray([1]);
        for (let i = 0; i < numSteps; ++i) {
            const result = solver.update(t, x);
            t = result.tOut;
            x = result.xOut;
        }
        expect(x.get(0)).toBeCloseTo(Math.pow(1 + h, numSteps), 12);
        // The Euler iterate underestimates e = 2.71828...
        expect(x.get(0)).toBeLessThan(Math.E);
    });
});
