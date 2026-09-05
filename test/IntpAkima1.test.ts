import { describe, it, expect } from 'vitest';
import { IntpAkima1, IntpAkima1Polynomial } from '../src/IntpAkima1.js';
import { check, expectClose, fc, finite, scaled, wellScaled } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V28): property-based cross-checks against IntpAkima1.h.
// ---------------------------------------------------------------------------
describe('IntpAkima1 verification', () => {
    // computeDerivative is the Akima slope estimator; the concrete uniform and
    // nonuniform subclasses feed it a sliding window of four slopes.
    const cd = (s: number[]): number => {
        const t = new TestInterp([0, 1, 2], [0, 0, 0]);
        return t.derivative(s, 0);
    };
    // wellScaled snaps |x| < 1e-3 to exactly zero, so no slope is subnormal;
    // a subnormal would underflow when multiplied by the power-of-two scale
    // below and break the exact-scaling argument.
    const slopes = () => fc.array(wellScaled(-6, 6), { minLength: 4, maxLength: 4 });
    // +/- 2^m, so multiplying a double by it is exact.
    const signedPowerOfTwo = () =>
        fc.tuple(fc.integer({ min: -4, max: 4 }), fc.boolean())
            .map(([m, neg]) => (neg ? -1 : 1) * Math.pow(2, m));

    it('always returns a value between the two middle slopes', () => {
        // This is the property that makes the Akima spline overshoot-free:
        // the estimate is a convex combination of slope[1] and slope[2] in the
        // generic branch and one of them (or their midpoint) otherwise.
        check(slopes(), s => {
            const d = cd(s);
            const lo = Math.min(s[1], s[2]);
            const hi = Math.max(s[1], s[2]);
            expect(d).toBeGreaterThanOrEqual(lo - 1e-12 * (1 + Math.abs(lo)));
            expect(d).toBeLessThanOrEqual(hi + 1e-12 * (1 + Math.abs(hi)));
        });
    });

    it('reproduces the common slope of linear data exactly', () => {
        check(finite(-6, 6), m => {
            expect(cd([m, m, m, m])).toBe(m);
        });
    });

    it('is invariant under reversal of the slope window', () => {
        // Reversing the sample order negates the abscissa and the ordinate
        // slopes together, so the estimate must be unchanged. Every branch of
        // upstream's decision tree maps onto its mirror image.
        check(slopes(), s => {
            expect(cd(s)).toBe(cd([s[3], s[2], s[1], s[0]]));
        });
    });

    it('is homogeneous of degree one in the slopes', () => {
        // The scale is a signed power of two so that k*slope, the differences
        // |s3-s2| and |s0-s1| and their quotient all scale exactly. With a
        // general k the weights are formed from differences of nearly equal
        // slopes, and that cancellation makes the ratio -- though not the
        // result's enclosure between s1 and s2 -- arbitrarily ill-conditioned.
        check(fc.tuple(slopes(), signedPowerOfTwo()), ([s, k]) => {
            // '+ 0' normalizes -0 to 0; toBe uses Object.is, which separates
            // them, and a zero slope makes the two sides differ only in sign
            // of zero.
            expect(cd(s.map(si => si * k)) + 0).toBe(k * cd(s) + 0);
        });
    });

    it('is exact when only one side has curvature', () => {
        // slope[0] == slope[1] != slope[2] == slope[3] is the "corner" case;
        // upstream averages the two middle slopes there.
        check(fc.tuple(finite(-5, 5), finite(-5, 5)), ([a, b]) => {
            expectClose(cd([a, a, b, b]), 0.5 * (a + b), 1e-15, 1e-15);
            // One-sided flatness picks the slope on the flat side's far end.
            if (a !== b) {
                expect(cd([a, a, b, 2 * b - a])).toBe(a);
                expect(cd([2 * a - b, a, b, b])).toBe(b);
            }
        });
    });

    it('clamps the query to [xMin, xMax] before evaluating', () => {
        check(fc.tuple(fc.array(finite(-5, 5), { minLength: 5, maxLength: 5 }),
            fc.array(finite(-3, 3), { minLength: 5, maxLength: 5 }),
            finite(-40, 40)), ([F, FD, x]) => {
            const t = new TestInterp(F, FD);
            const clamped = Math.min(Math.max(x, t.getXMin()), t.getXMax());
            expect(t.evaluate(x)).toBe(t.evaluate(clamped));
            for (let order = 0; order <= 3; ++order) {
                expect(t.evaluate(order, x)).toBe(t.evaluate(order, clamped));
            }
        });
    });

    it('interpolates the data and slopes it was built from', () => {
        check(fc.tuple(fc.array(finite(-5, 5), { minLength: 5, maxLength: 5 }),
            fc.array(finite(-3, 3), { minLength: 5, maxLength: 5 })), ([F, FD]) => {
            const t = new TestInterp(F, FD);
            for (let i = 0; i < F.length; ++i) {
                expectClose(t.evaluate(i), F[i], 1e-12, 1e-12);
                expectClose(t.evaluate(1, i), FD[i], 1e-11, 1e-11);
            }
        });
    });

    it('reports the first derivative consistently with a central difference', () => {
        const h = 1e-5;
        check(fc.tuple(fc.array(wellScaled(-5, 5), { minLength: 5, maxLength: 5 }),
            fc.array(wellScaled(-3, 3), { minLength: 5, maxLength: 5 }),
            scaled(0.05, 3.95)), ([F, FD, x]) => {
            const t = new TestInterp(F, FD);
            const fd = (t.evaluate(0, x + h) - t.evaluate(0, x - h)) / (2 * h);
            // Central differences are exact for cubics up to the h^2/6 * P'''
            // term (P''' <= ~160 here) plus roundoff O(eps * |P| / h) ~ 1e-10.
            // scaled() draws x from a 4096-point lattice whose step is 9.5e-4,
            // so x is never within h = 1e-5 of a knot and both offsets lie in
            // the same cell.
            expectClose(t.evaluate(1, x), fd, 1e-4, 1e-6);
        });
    });

    it('polynomial evaluate(order, x) differentiates the cubic exactly', () => {
        check(fc.tuple(fc.array(finite(-5, 5), { minLength: 4, maxLength: 4 }),
            scaled(-2, 2)), ([c, x]) => {
            const p = new IntpAkima1Polynomial();
            for (let i = 0; i < 4; ++i) { p.setCoeff(i, c[i]); }
            expect(p.evaluate(0, x)).toBe(p.evaluate(x));
            expectClose(p.evaluate(1, x), c[1] + 2 * c[2] * x + 3 * c[3] * x * x,
                1e-12, 1e-13);
            expectClose(p.evaluate(2, x), 2 * c[2] + 6 * c[3] * x, 1e-12, 1e-13);
            expect(p.evaluate(3, x)).toBe(6 * c[3]);
            expect(p.evaluate(4, x)).toBe(0);
        });
    });

    it('rejects fewer than three samples', () => {
        expect(() => new TestInterp([1, 2], [0, 0])).toThrow();
        expect(() => new TestInterp([1], [0])).toThrow();
    });
});
