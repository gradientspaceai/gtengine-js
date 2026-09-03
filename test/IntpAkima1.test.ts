import { describe, it, expect } from 'vitest';
import { IntpAkima1, IntpAkima1Polynomial } from '../src/IntpAkima1.js';

// Minimal concrete subclass on the uniform lattice x = 0, 1, ..., q-1 whose
// segment polynomials are Hermite cubics built from prescribed derivatives.
// It exercises the abstract base machinery (polynomial storage, clamped
// evaluation, lookup dispatch, computeDerivative).
class TestInterp extends IntpAkima1 {
    constructor(F: number[], FD: number[]) {
        super(F.length, F);
        for (let i = 0; i < F.length - 1; ++i) {
            // Cubic Hermite on [0,1] local coordinates: matches F[i], FD[i]
            // at dx=0 and F[i+1], FD[i+1] at dx=1.
            const p = this.mPoly[i];
            const df = F[i + 1] - F[i];
            p.setCoeff(0, F[i]);
            p.setCoeff(1, FD[i]);
            p.setCoeff(2, 3 * df - 2 * FD[i] - FD[i + 1]);
            p.setCoeff(3, -2 * df + FD[i] + FD[i + 1]);
        }
    }

    getXMin(): number {
        return 0;
    }

    getXMax(): number {
        return this.mQuantity - 1;
    }

    protected lookup(x: number): { index: number, dx: number } {
        const index = Math.min(Math.floor(x), this.mQuantity - 2);
        return { index, dx: x - index };
    }

    // Expose the protected derivative estimator for testing.
    derivative(slope: readonly number[], offset: number): number {
        return this.computeDerivative(slope, offset);
    }
}

describe('IntpAkima1Polynomial', () => {
    it('evaluates the cubic and its derivatives', () => {
        // P(x) = 1 + 2x + 3x^2 + 4x^3.
        const p = new IntpAkima1Polynomial();
        p.setCoeff(0, 1);
        p.setCoeff(1, 2);
        p.setCoeff(2, 3);
        p.setCoeff(3, 4);
        expect(p.getCoeff(2)).toBe(3);

        expect(p.evaluate(0.5)).toBeCloseTo(1 + 1 + 0.75 + 0.5, 14);
        expect(p.evaluate(0, 0.5)).toBe(p.evaluate(0.5));
        // P'(x) = 2 + 6x + 12x^2.
        expect(p.evaluate(1, 0.5)).toBeCloseTo(2 + 3 + 3, 14);
        // P''(x) = 6 + 24x.
        expect(p.evaluate(2, 0.5)).toBeCloseTo(18, 14);
        // P'''(x) = 24.
        expect(p.evaluate(3, 0.5)).toBe(24);
        // Orders beyond the degree return zero.
        expect(p.evaluate(4, 0.5)).toBe(0);
    });

    it('defaults to the zero polynomial', () => {
        const p = new IntpAkima1Polynomial();
        for (let order = 0; order <= 3; ++order) {
            expect(p.evaluate(order, 0.7)).toBe(0);
        }
    });
});

describe('IntpAkima1', () => {
    it('throws when fewer than three points are supplied', () => {
        expect(() => new TestInterp([1, 2], [0, 0])).toThrow(
            'Invalid input to IntpAkima1 constructor.');
    });

    it('provides member access', () => {
        const F = [1, 4, 9, 16];
        const interp = new TestInterp(F, [0, 0, 0, 0]);
        expect(interp.getQuantity()).toBe(4);
        expect(interp.getF()).toBe(F);
        expect(interp.getXMin()).toBe(0);
        expect(interp.getXMax()).toBe(3);
    });

    it('interpolates the samples and derivatives', () => {
        // f(x) = x^2 on x = 0..3 with exact derivatives 2x; the piecewise
        // Hermite cubic then reproduces x^2 exactly.
        const F = [0, 1, 4, 9];
        const FD = [0, 2, 4, 6];
        const interp = new TestInterp(F, FD);

        for (let i = 0; i < 4; ++i) {
            expect(interp.evaluate(i)).toBeCloseTo(F[i], 13);
        }
        for (const x of [0.25, 0.5, 1.5, 2.75, 3]) {
            expect(interp.evaluate(x)).toBeCloseTo(x * x, 13);
            expect(interp.evaluate(0, x)).toBeCloseTo(x * x, 13);
            expect(interp.evaluate(1, x)).toBeCloseTo(2 * x, 13);
            expect(interp.evaluate(2, x)).toBeCloseTo(2, 13);
            expect(interp.evaluate(3, x)).toBeCloseTo(0, 13);
        }
        expect(interp.evaluate(4, 1.5)).toBe(0);
    });

    it('clamps evaluations to [xmin, xmax]', () => {
        const F = [0, 1, 4, 9];
        const FD = [0, 2, 4, 6];
        const interp = new TestInterp(F, FD);
        expect(interp.evaluate(-5)).toBeCloseTo(interp.evaluate(0), 14);
        expect(interp.evaluate(100)).toBeCloseTo(interp.evaluate(3), 14);
        expect(interp.evaluate(1, -5)).toBeCloseTo(interp.evaluate(1, 0), 14);
        expect(interp.evaluate(1, 100)).toBeCloseTo(interp.evaluate(1, 3), 14);
    });

    describe('computeDerivative', () => {
        const interp = new TestInterp([0, 1, 4], [0, 2, 4]);

        it('uses the Akima weighted average when all slopes differ', () => {
            // ad0 = |8-4| = 4, ad1 = |1-2| = 1,
            // result = (4*2 + 1*4) / 5 = 12/5.
            expect(interp.derivative([1, 2, 4, 8], 0)).toBeCloseTo(2.4, 14);
        });

        it('returns slope[1] when the middle slopes are equal', () => {
            expect(interp.derivative([0, 2, 2, 9], 0)).toBe(2);
        });

        it('returns slope[2] when only the right pair is equal', () => {
            // slope0 != slope1, slope2 == slope3 -> slope[2].
            expect(interp.derivative([0, 1, 3, 3], 0)).toBe(3);
        });

        it('returns slope[1] when the left pair is equal and right pair differs', () => {
            expect(interp.derivative([2, 2, 3, 5], 0)).toBe(2);
        });

        it('averages the middle slopes when both outer pairs are equal', () => {
            expect(interp.derivative([1, 1, 3, 3], 0)).toBe(2);
        });

        it('respects the window offset', () => {
            expect(interp.derivative([99, 42, 1, 2, 4, 8, 77], 2)).toBeCloseTo(2.4, 14);
        });
    });
});
