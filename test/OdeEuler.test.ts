import { describe, it, expect } from 'vitest';
import { OdeEuler } from '../src/OdeEuler.js';
import { OdeSolver } from '../src/OdeSolver.js';
import { Vector } from '../src/Vector.js';
import {
    check, expectClose, fc, scaled, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V39): independent review against upstream OdeEuler.h.
// ---------------------------------------------------------------------------

describe('OdeEuler verification', () => {
    // The decay/growth rate of dx/dt = a*x, kept away from zero so that the
    // discretization error - and therefore the measured order - is not lost
    // in round-off.
    const rate = fc.oneof(scaled(-2, -0.5), scaled(0.5, 2));

    it('is first-order accurate on dx/dt = a*x', () => {
        check(rate, a => {
            const F = (_t: number, x: Vector) =>
                Vector.fromArray([a * x.get(0)]);
            const errorForSteps = (numSteps: number): number => {
                const solver = new OdeEuler(1 / numSteps, F);
                let t = 0;
                let x = Vector.fromArray([1]);
                for (let i = 0; i < numSteps; ++i) {
                    const r = solver.update(t, x);
                    t = r.tOut;
                    x = r.xOut;
                }
                return Math.abs(x.get(0) - Math.exp(a));
            };
            // Halving the step halves the global error of a first-order
            // scheme. The observed ratios over this range of a lie in
            // [1.96, 2.01]; the band allows for the O(h) correction to the
            // asymptotic value 2.
            const ratio = errorForSteps(64) / errorForSteps(128);
            expect(ratio).toBeGreaterThan(1.8);
            expect(ratio).toBeLessThan(2.2);
        });
    });

    it('reproduces the closed-form iterate (1 + a*h)^n for dx/dt = a*x',
        () => {
            check(fc.tuple(rate, fc.integer({ min: 1, max: 24 })),
                ([a, numSteps]) => {
                    const h = 1 / numSteps;
                    const solver = new OdeEuler(h,
                        (_t, x) => Vector.fromArray([a * x.get(0)]));
                    let t = 0;
                    let x = Vector.fromArray([1]);
                    for (let i = 0; i < numSteps; ++i) {
                        const r = solver.update(t, x);
                        t = r.tOut;
                        x = r.xOut;
                    }
                    // Each step multiplies by (1 + a*h) up to one rounding,
                    // so n steps agree with the power to within n ulps.
                    expectClose(x.get(0), Math.pow(1 + a * h, numSteps), 0,
                        numSteps * 4 * Number.EPSILON);
                });
        });

    it('integrates a right-hand side independent of x exactly', () => {
        check(fc.tuple(wellScaledVector(4, -5, 5), scaled(-2, 2),
            scaled(-3, 3)), ([c, tDelta, tIn]) => {
            // F(t,x) = c makes the Euler step exact: x(t + h) = x + h*c.
            const solver = new OdeEuler(tDelta, () => c);
            const xIn = Vector.fromArray([1, 2, 3, 4]);
            const { tOut, xOut } = solver.update(tIn, xIn);
            expect(tOut).toBe(tIn + tDelta);
            for (let i = 0; i < 4; ++i) {
                expect(xOut.get(i)).toBe(xIn.get(i) + tDelta * c.get(i));
            }
        });
    });

    it('never aliases its input or its output', () => {
        check(fc.tuple(wellScaledVector(3, -5, 5), scaled(-2, 2)),
            ([xIn, tDelta]) => {
                const before = xIn.clone();
                const solver = new OdeEuler(tDelta, (_t, x) => x);
                const { xOut } = solver.update(0, xIn);
                expect(xOut).not.toBe(xIn);
                for (let i = 0; i < 3; ++i) {
                    expect(xIn.get(i)).toBe(before.get(i));
                }
            });
    });
});
