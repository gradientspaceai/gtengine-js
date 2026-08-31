import { describe, it, expect } from 'vitest';
import { OdeMidpoint } from '../src/OdeMidpoint';
import { OdeSolver } from '../src/OdeSolver';
import { Vector } from '../src/Vector';

describe('OdeMidpoint', () => {
    it('is an OdeSolver whose step size is readable and writable', () => {
        const solver = new OdeMidpoint(0.25, (_t, x) => x);
        expect(solver).toBeInstanceOf(OdeSolver);
        expect(solver.getTDelta()).toBe(0.25);
        solver.setTDelta(0.5);
        expect(solver.getTDelta()).toBe(0.5);
    });

    it('takes the exact midpoint step for dx/dt = x', () => {
        // xTemp = x*(1 + h/2), so xOut = x*(1 + h + h^2/2).
        const h = 0.125;
        const solver = new OdeMidpoint(h, (_t, x) => x);
        const { tOut, xOut } = solver.update(0, Vector.fromArray([1]));
        expect(tOut).toBe(h);
        expect(xOut.get(0)).toBeCloseTo(1 + h + 0.5 * h * h, 15);
    });

    it('evaluates F at the interval midpoint', () => {
        const calls: Array<{ t: number; x: number }> = [];
        const solver = new OdeMidpoint(0.4, (t, x) => {
            calls.push({ t, x: x.get(0) });
            return Vector.fromArray([1]);
        });
        const { tOut, xOut } = solver.update(1, Vector.fromArray([5]));
        expect(calls.length).toBe(2);
        expect(calls[0]).toEqual({ t: 1, x: 5 });
        // The second evaluation is at t + h/2 and x + (h/2)*F.
        expect(calls[1].t).toBeCloseTo(1.2, 15);
        expect(calls[1].x).toBeCloseTo(5.2, 15);
        expect(tOut).toBeCloseTo(1.4, 15);
        expect(xOut.get(0)).toBeCloseTo(5.4, 15);
    });

    it('does not modify the input vector', () => {
        const solver = new OdeMidpoint(0.5, (_t, x) => x);
        const xIn = Vector.fromArray([1, 2]);
        const { xOut } = solver.update(0, xIn);
        expect(xIn.values).toEqual([1, 2]);
        expect(xOut).not.toBe(xIn);
    });

    it('integrates a time-dependent right-hand side exactly for quadratics', () => {
        // x' = 2*t with x(0) = 0 has the solution t^2. The midpoint rule is
        // exact for a linear integrand, so a single step is exact.
        const solver = new OdeMidpoint(2, (t, _x) => Vector.fromArray([2 * t]));
        const { tOut, xOut } = solver.update(0, Vector.fromArray([0]));
        expect(tOut).toBe(2);
        expect(xOut.get(0)).toBeCloseTo(4, 12);
    });

    it('solves the harmonic oscillator with second-order accuracy', () => {
        const F = (_t: number, x: Vector) =>
            Vector.fromArray([x.get(1), -x.get(0)]);

        const errorForSteps = (numSteps: number): number => {
            const tFinal = 1;
            const h = tFinal / numSteps;
            const solver = new OdeMidpoint(h, F);
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

        const e0 = errorForSteps(50);
        const e1 = errorForSteps(100);
        // The ratio of errors is approximately 2^2 for a second-order method.
        expect(e0 / e1).toBeGreaterThan(3.5);
        expect(e0 / e1).toBeLessThan(4.5);
    });
});
