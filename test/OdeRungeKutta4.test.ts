import { describe, it, expect } from 'vitest';
import { OdeRungeKutta4 } from '../src/OdeRungeKutta4.js';
import { OdeSolver } from '../src/OdeSolver.js';
import { Vector } from '../src/Vector.js';
import {
    check, expectClose, fc, scaled, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V39): independent review against upstream OdeRungeKutta4.h.
// ---------------------------------------------------------------------------

describe('OdeRungeKutta4 verification', () => {
    const rate = fc.oneof(scaled(-2, -0.5), scaled(0.5, 2));

    it('is fourth-order accurate on dx/dt = a*x', () => {
        check(rate, a => {
            const F = (_t: number, x: Vector) =>
                Vector.fromArray([a * x.get(0)]);
            const errorForSteps = (numSteps: number): number => {
                const solver = new OdeRungeKutta4(1 / numSteps, F);
                let t = 0;
                let x = Vector.fromArray([1]);
                for (let i = 0; i < numSteps; ++i) {
                    const r = solver.update(t, x);
                    t = r.tOut;
                    x = r.xOut;
                }
                return Math.abs(x.get(0) - Math.exp(a));
            };
            // Halving the step divides the global error of a fourth-order
            // scheme by 16; the observed ratios lie in [15.2, 16.9]. The
            // step counts are large enough for the asymptotic regime and
            // small enough that the errors (about 1e-6) stay far above the
            // round-off floor.
            const ratio = errorForSteps(16) / errorForSteps(32);
            expect(ratio).toBeGreaterThan(14);
            expect(ratio).toBeLessThan(18);
        });
    });

    it('reproduces the degree-4 Taylor amplification for dx/dt = a*x', () => {
        check(fc.tuple(rate, fc.integer({ min: 1, max: 24 })),
            ([a, numSteps]) => {
                // The classical RK4 step applied to a linear equation
                // multiplies x by the degree-4 Taylor polynomial of
                // exp(a*h).
                const h = 1 / numSteps;
                const solver = new OdeRungeKutta4(h,
                    (_t, x) => Vector.fromArray([a * x.get(0)]));
                let t = 0;
                let x = Vector.fromArray([1]);
                for (let i = 0; i < numSteps; ++i) {
                    const r = solver.update(t, x);
                    t = r.tOut;
                    x = r.xOut;
                }
                const z = a * h;
                const factor = 1 + z * (1 + z * (0.5 + z * (1 / 6 + z / 24)));
                expectClose(x.get(0), Math.pow(factor, numSteps), 0,
                    numSteps * 16 * Number.EPSILON);
            });
    });

    it('applies Simpson quadrature when F depends only on t', () => {
        check(fc.tuple(fc.array(scaled(-3, 3), { minLength: 4, maxLength: 4 }),
            scaled(-2, 2), scaled(-3, 3)), ([c, tDelta, tIn]) => {
            // Both middle stages evaluate F at t + h/2, so the update is
            // x + (h/6)*(f(t) + 4*f(t + h/2) + f(t + h)): Simpson's rule,
            // which is exact for integrands of degree <= 3.
            const f = (t: number) =>
                c[0] + t * (c[1] + t * (c[2] + t * c[3]));
            const antiderivative = (t: number) =>
                t * (c[0] + t * (c[1] / 2 + t * (c[2] / 3 + t * c[3] / 4)));
            const solver = new OdeRungeKutta4(tDelta,
                (t, _x) => Vector.fromArray([f(t)]));
            const { tOut, xOut } = solver.update(tIn, Vector.fromArray([7]));
            expect(tOut).toBe(tIn + tDelta);
            const exact = 7 + antiderivative(tIn + tDelta)
                - antiderivative(tIn);
            // Cancellation in the antiderivative difference costs a few
            // digits when tIn is far from zero and tDelta is small.
            expectClose(xOut.get(0), exact, 1e-11, 1e-11);
        });
    });

    it('evaluates F at the four Runge-Kutta nodes in order', () => {
        check(fc.tuple(scaled(-2, 2), scaled(-3, 3)), ([tDelta, tIn]) => {
            const times: number[] = [];
            const solver = new OdeRungeKutta4(tDelta, (t, x) => {
                times.push(t);
                return x;
            });
            solver.update(tIn, Vector.fromArray([1]));
            expect(times.length).toBe(4);
            expect(times[0]).toBe(tIn);
            expect(times[1]).toBe(tIn + 0.5 * tDelta);
            expect(times[2]).toBe(tIn + 0.5 * tDelta);
            expect(times[3]).toBe(tIn + tDelta);
        });
    });

    it('nearly conserves the energy of a harmonic oscillator', () => {
        check(scaled(0.5, 2), omega => {
            // x'' = -omega^2 * x written as a first-order system. The energy
            // omega^2*x^2 + v^2 is conserved exactly by the flow, and RK4
            // drifts by O(h^4) per unit time.
            const h = 0.01;
            const solver = new OdeRungeKutta4(h, (_t, s) =>
                Vector.fromArray([s.get(1), -omega * omega * s.get(0)]));
            const energy = (s: Vector) =>
                omega * omega * s.get(0) * s.get(0) + s.get(1) * s.get(1);
            let t = 0;
            let s = Vector.fromArray([1, 0]);
            const e0 = energy(s);
            for (let i = 0; i < 500; ++i) {
                const r = solver.update(t, s);
                t = r.tOut;
                s = r.xOut;
            }
            // The amplification factor per step is |P4(i*omega*h)|, which is
            // 1 - (omega*h)^6/72 + ...; over 500 steps of h = 0.01 with
            // omega <= 2 the relative drift stays below 1e-6.
            expect(Math.abs(energy(s) / e0 - 1)).toBeLessThan(1e-6);
        });
    });

    it('never aliases its input or its output', () => {
        check(fc.tuple(wellScaledVector(3, -5, 5), scaled(-2, 2)),
            ([xIn, tDelta]) => {
                const before = xIn.clone();
                const solver = new OdeRungeKutta4(tDelta, (_t, x) => x);
                const { xOut } = solver.update(0, xIn);
                expect(xOut).not.toBe(xIn);
                for (let i = 0; i < 3; ++i) {
                    expect(xIn.get(i)).toBe(before.get(i));
                }
            });
    });
});
