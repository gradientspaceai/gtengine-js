import { describe, it, expect } from 'vitest';
import { OdeSolver, type OdeFunction } from '../src/OdeSolver.js';
import { OdeEuler } from '../src/OdeEuler.js';
import { OdeMidpoint } from '../src/OdeMidpoint.js';
import { OdeRungeKutta4 } from '../src/OdeRungeKutta4.js';
import { Vector, mul } from '../src/Vector.js';
import { check, expectClose, fc, finite as vfinite } from './helpers/arbitraries.js';

// A minimal concrete subclass (an explicit Euler step) used to exercise the
// abstract base class. The real OdeEuler port arrives in a later batch.
class TestEulerSolver extends OdeSolver<number[]> {
    constructor(tDelta: number, F: OdeFunction<number[]>) {
        super(tDelta, F);
    }

    update(tIn: number, xIn: number[]): { tOut: number; xOut: number[] } {
        const fValue = this.mFunction(tIn, xIn);
        const xOut = xIn.map((x, i) => x + this.mTDelta * fValue[i]);
        return { tOut: tIn + this.mTDelta, xOut };
    }
}

describe('OdeSolver', () => {
    it('stores and updates tDelta via the accessors', () => {
        const solver = new TestEulerSolver(0.25, (_t, x) => x);
        expect(solver.getTDelta()).toBe(0.25);
        solver.setTDelta(0.125);
        expect(solver.getTDelta()).toBe(0.125);
    });

    it('passes t and x to the stored function F', () => {
        let seenT = NaN;
        let seenX: number[] = [];
        const solver = new TestEulerSolver(0.5, (t, x) => {
            seenT = t;
            seenX = x;
            return [0, 0];
        });
        const { tOut, xOut } = solver.update(2, [3, 4]);
        expect(seenT).toBe(2);
        expect(seenX).toEqual([3, 4]);
        expect(tOut).toBe(2.5);
        expect(xOut).toEqual([3, 4]);
    });

    it('a subclass update approximates dx/dt = x (exponential growth)', () => {
        // Solve dx/dt = x with x(0) = 1 to t = 1 using many small Euler
        // steps; the result must approach e.
        const n = 100000;
        const solver = new TestEulerSolver(1 / n, (_t, x) => [x[0]]);
        let t = 0;
        let x = [1];
        for (let i = 0; i < n; ++i) {
            ({ tOut: t, xOut: x } = solver.update(t, x));
        }
        expect(t).toBeCloseTo(1, 10);
        expect(x[0]).toBeCloseTo(Math.E, 4);
    });

    it('a subclass update approximates a time-dependent equation dx/dt = 2t', () => {
        // x(t) = t^2 with x(0) = 0.
        const n = 20000;
        const solver = new TestEulerSolver(1 / n, (t, _x) => [2 * t]);
        let t = 0;
        let x = [0];
        for (let i = 0; i < n; ++i) {
            ({ tOut: t, xOut: x } = solver.update(t, x));
        }
        // Explicit Euler on dx/dt = 2t has global error exactly h = 1/n.
        expect(Math.abs(x[0] - 1)).toBeLessThanOrEqual(1.5 / n);
    });
});

// ---------------------------------------------------------------------------
// Verification (group V38): OdeSolver.h is an abstract base class holding
// tDelta and F, so the properties below pin the base-class contract and
// exercise it through the concrete solvers that derive from it.
// ---------------------------------------------------------------------------
describe('OdeSolver verification', () => {
    it('setTDelta/getTDelta round-trip and the stored value is used by update',
        () => {
            check(fc.tuple(vfinite(-5, 5), vfinite(-5, 5), vfinite(-5, 5)),
                ([h0, h1, x0]) => {
                    // dx/dt = 1, so an explicit Euler step advances x by
                    // exactly tDelta and t by exactly tDelta.
                    const solver = new TestEulerSolver(h0, () => [1]);
                    expect(solver.getTDelta()).toBe(h0);
                    solver.setTDelta(h1);
                    expect(solver.getTDelta()).toBe(h1);
                    const { tOut, xOut } = solver.update(0, [x0]);
                    // '+ 0' normalizes the -0/+0 tie that toBe (Object.is)
                    // would otherwise reject when h1 is 0.
                    expect(tOut + 0).toBe(0 + h1 + 0);
                    expect(xOut[0] + 0).toBe(x0 + h1 + 0);
                });
        });

    it('a derived solver receives (t, x) with the documented arity', () => {
        // The OdeFunction type is (t, x) => TVector. TypeScript accepts a
        // lower-arity callback where a higher-arity one is declared, so a
        // one-argument lambda would silently bind to t rather than to x.
        // Pin that F is called with two arguments in the (t, x) order.
        const seen: Array<[number, number[]]> = [];
        const F: OdeFunction<number[]> = function (t, x) {
            expect(arguments.length).toBe(2);
            seen.push([t, x.slice()]);
            return [0];
        };
        const solver = new TestEulerSolver(0.5, F);
        solver.update(7, [11]);
        expect(seen).toEqual([[7, [11]]]);
    });

    it('the concrete solvers of the library satisfy the base-class contract',
        () => {
            const F: OdeFunction<Vector> = (_t, x) => x.clone();
            for (const solver of [
                new OdeEuler(0.25, F),
                new OdeMidpoint(0.25, F),
                new OdeRungeKutta4(0.25, F)
            ]) {
                expect(solver).toBeInstanceOf(OdeSolver);
                expect(solver.getTDelta()).toBe(0.25);
                const xIn = Vector.fromArray([1, 2, 3]);
                const { tOut, xOut } = solver.update(1, xIn);
                expect(tOut).toBeCloseTo(1.25, 12);
                // Upstream states that a derived Update must tolerate xIn
                // and xOut being the same object; the ported solvers build a
                // new Vector, so xIn must be untouched.
                expect(xIn.values).toEqual([1, 2, 3]);
                expect(xOut).not.toBe(xIn);
            }
        });

    it('one step of the derived solvers reproduces their Taylor expansions',
        () => {
            // dx/dt = lambda*x has the exact solution x(t) = x0*exp(lambda*t).
            // Euler reproduces the degree-1 Taylor polynomial of exp(lambda*h),
            // midpoint the degree-2 one and RK4 the degree-4 one.
            check(fc.tuple(vfinite(-2, 2), vfinite(-2, 2), vfinite(-3, 3)),
                ([lambda, h, x0]) => {
                    const F: OdeFunction<Vector> = (_t, x) => mul(lambda, x);
                    const z = lambda * h;
                    const cases: Array<[OdeSolver<Vector>, number]> = [
                        [new OdeEuler(h, F), 1 + z],
                        [new OdeMidpoint(h, F), 1 + z + z * z / 2],
                        [new OdeRungeKutta4(h, F),
                            1 + z + z * z / 2 + z ** 3 / 6 + z ** 4 / 24]
                    ];
                    for (const [solver, factor] of cases) {
                        const { xOut } = solver.update(0, Vector.fromArray([x0]));
                        expectClose(xOut.get(0), x0 * factor, 1e-12, 1e-12);
                    }
                });
        });

    it('the derived solvers have the documented order of accuracy', () => {
        // Integrate dx/dt = x from 0 to 1 with n steps and compare with e.
        // Halving h must reduce the global error by 2^p with p the order.
        const F: OdeFunction<Vector> = (_t, x) => x.clone();
        const errorAt = (make: (h: number) => OdeSolver<Vector>, n: number) => {
            const solver = make(1 / n);
            let t = 0;
            let x = Vector.fromArray([1]);
            for (let i = 0; i < n; ++i) {
                ({ tOut: t, xOut: x } = solver.update(t, x));
            }
            return Math.abs(x.get(0) - Math.E);
        };
        const orders: Array<[(h: number) => OdeSolver<Vector>, number]> = [
            [h => new OdeEuler(h, F), 1],
            [h => new OdeMidpoint(h, F), 2],
            [h => new OdeRungeKutta4(h, F), 4]
        ];
        for (const [make, p] of orders) {
            const e0 = errorAt(make, 64);
            const e1 = errorAt(make, 128);
            // The ratio tends to 2^p from below/above; allow 20%.
            expect(e0 / e1).toBeGreaterThan(2 ** p * 0.8);
            expect(e0 / e1).toBeLessThan(2 ** p * 1.25);
        }
    });
});
