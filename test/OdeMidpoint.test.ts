import { describe, it, expect } from 'vitest';
import { OdeMidpoint } from '../src/OdeMidpoint.js';
import { OdeSolver } from '../src/OdeSolver.js';
import { Vector } from '../src/Vector.js';
import {
    check, expectClose, fc, scaled, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V39): independent review against upstream OdeMidpoint.h.
// ---------------------------------------------------------------------------

describe('OdeMidpoint verification', () => {
    const rate = fc.oneof(scaled(-2, -0.5), scaled(0.5, 2));

    it('is second-order accurate on dx/dt = a*x', () => {
        check(rate, a => {
            const F = (_t: number, x: Vector) =>
                Vector.fromArray([a * x.get(0)]);
            const errorForSteps = (numSteps: number): number => {
                const solver = new OdeMidpoint(1 / numSteps, F);
                let t = 0;
                let x = Vector.fromArray([1]);
                for (let i = 0; i < numSteps; ++i) {
                    const r = solver.update(t, x);
                    t = r.tOut;
                    x = r.xOut;
                }
                return Math.abs(x.get(0) - Math.exp(a));
            };
            // Halving the step quarters the global error of a second-order
            // scheme; the observed ratios lie in [3.91, 4.10].
            const ratio = errorForSteps(32) / errorForSteps(64);
            expect(ratio).toBeGreaterThan(3.6);
            expect(ratio).toBeLessThan(4.4);
        });
    });

    it('reproduces the closed-form iterate (1 + a*h + (a*h)^2/2)^n', () => {
        check(fc.tuple(rate, fc.integer({ min: 1, max: 24 })),
            ([a, numSteps]) => {
                // One midpoint step of dx/dt = a*x multiplies x by
                // 1 + a*h*(1 + a*h/2), the degree-2 Taylor polynomial of
                // exp(a*h).
                const h = 1 / numSteps;
                const solver = new OdeMidpoint(h,
                    (_t, x) => Vector.fromArray([a * x.get(0)]));
                let t = 0;
                let x = Vector.fromArray([1]);
                for (let i = 0; i < numSteps; ++i) {
                    const r = solver.update(t, x);
                    t = r.tOut;
                    x = r.xOut;
                }
                const factor = 1 + a * h * (1 + 0.5 * a * h);
                expectClose(x.get(0), Math.pow(factor, numSteps), 0,
                    numSteps * 8 * Number.EPSILON);
            });
    });

    it('applies the midpoint quadrature rule when F depends only on t', () => {
        check(fc.tuple(scaled(-3, 3), scaled(-3, 3), scaled(-2, 2),
            scaled(-3, 3)), ([c0, c1, tDelta, tIn]) => {
            // For dx/dt = c0 + c1*t the update is x + h*(c0 + c1*(t + h/2)),
            // which is the exact integral of a linear integrand.
            const solver = new OdeMidpoint(tDelta,
                (t, _x) => Vector.fromArray([c0 + c1 * t]));
            const xIn = Vector.fromArray([7]);
            const { tOut, xOut } = solver.update(tIn, xIn);
            expect(tOut).toBe(tIn + tDelta);
            const exact = 7 + c0 * tDelta
                + 0.5 * c1 * ((tIn + tDelta) * (tIn + tDelta) - tIn * tIn);
            expectClose(xOut.get(0), exact, 1e-12, 1e-12);
        });
    });

    it('evaluates F at the midpoint of the step, not at either end', () => {
        check(fc.tuple(scaled(-2, 2), scaled(-3, 3)), ([tDelta, tIn]) => {
            const times: number[] = [];
            const solver = new OdeMidpoint(tDelta, (t, x) => {
                times.push(t);
                return x;
            });
            solver.update(tIn, Vector.fromArray([1]));
            expect(times.length).toBe(2);
            expect(times[0]).toBe(tIn);
            expect(times[1]).toBe(tIn + 0.5 * tDelta);
        });
    });

    it('never aliases its input or its output', () => {
        check(fc.tuple(wellScaledVector(3, -5, 5), scaled(-2, 2)),
            ([xIn, tDelta]) => {
                const before = xIn.clone();
                const solver = new OdeMidpoint(tDelta, (_t, x) => x);
                const { xOut } = solver.update(0, xIn);
                expect(xOut).not.toBe(xIn);
                for (let i = 0; i < 3; ++i) {
                    expect(xIn.get(i)).toBe(before.get(i));
                }
            });
    });
});
