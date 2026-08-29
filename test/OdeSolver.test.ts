import { describe, it, expect } from 'vitest';
import { OdeSolver, type OdeFunction } from '../src/OdeSolver';

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
