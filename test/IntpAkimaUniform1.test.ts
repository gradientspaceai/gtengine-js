import { describe, it, expect } from 'vitest';
import { IntpAkimaUniform1 } from '../src/IntpAkimaUniform1.js';
import { check, expectClose, fc, finite, seededRandom }
    from './helpers/arbitraries.js';

describe('IntpAkimaUniform1', () => {
    it('throws for invalid inputs', () => {
        expect(() => new IntpAkimaUniform1(2, 0, 1, [0, 1]))
            .toThrow('Invalid input to IntpAkima1 constructor.');
        expect(() => new IntpAkimaUniform1(3, 0, 0, [0, 1, 2]))
            .toThrow('Spacing must be positive.');
        expect(() => new IntpAkimaUniform1(3, 0, -1, [0, 1, 2]))
            .toThrow('Spacing must be positive.');
    });

    it('provides member access', () => {
        const F = [0, 1, 4, 9, 16];
        const interp = new IntpAkimaUniform1(5, -1, 0.5, F);
        expect(interp.getQuantity()).toBe(5);
        expect(interp.getF()).toBe(F);
        expect(interp.getXMin()).toBe(-1);
        expect(interp.getXSpacing()).toBe(0.5);
        expect(interp.getXMax()).toBeCloseTo(1, 15);
    });

    it('interpolates the samples exactly', () => {
        const xMin = -2;
        const dx = 0.75;
        const F = [1, -3, 2.5, 7, 0, -1.25, 4, 4];
        const interp = new IntpAkimaUniform1(F.length, xMin, dx, F);
        for (let i = 0; i < F.length; ++i) {
            expect(interp.evaluate(xMin + i * dx)).toBeCloseTo(F[i], 12);
            expect(interp.evaluate(0, xMin + i * dx)).toBeCloseTo(F[i], 12);
        }
    });

    it('reproduces linear data exactly', () => {
        // All slopes are equal, so the estimated derivatives are the common
        // slope and the quadratic and cubic coefficients vanish.
        const xMin = 1;
        const dx = 0.25;
        const a = -0.5, b = 3;
        const F: number[] = [];
        for (let i = 0; i < 6; ++i) {
            F.push(a + b * (xMin + i * dx));
        }
        const interp = new IntpAkimaUniform1(F.length, xMin, dx, F);
        for (let k = 0; k <= 40; ++k) {
            const x = xMin + (k / 40) * (F.length - 1) * dx;
            expect(interp.evaluate(x)).toBeCloseTo(a + b * x, 12);
            expect(interp.evaluate(1, x)).toBeCloseTo(b, 12);
            expect(interp.evaluate(2, x)).toBeCloseTo(0, 12);
            expect(interp.evaluate(3, x)).toBeCloseTo(0, 12);
        }
    });

    it('reproduces constant data exactly', () => {
        const F = [2, 2, 2, 2, 2];
        const interp = new IntpAkimaUniform1(5, 0, 1, F);
        for (let k = 0; k <= 20; ++k) {
            const x = (k / 20) * 4;
            expect(interp.evaluate(x)).toBeCloseTo(2, 14);
            expect(interp.evaluate(1, x)).toBeCloseTo(0, 14);
        }
    });

    it('averages the two slopes when both slope pairs are equal', () => {
        // The slopes are (1,1,1,2,2,2). At sample index 3 the four-slope
        // window is (1,1,2,2), so slope[0] == slope[1] and
        // slope[2] == slope[3] and the estimated derivative is
        // 0.5*(1 + 2) = 1.5.
        const F = [0, 1, 2, 3, 5, 7, 9];
        const interp = new IntpAkimaUniform1(F.length, 0, 1, F);
        expect(interp.evaluate(1, 3)).toBeCloseTo(1.5, 14);

        // Away from the corner the data is linear, so the derivative is the
        // local slope.
        expect(interp.evaluate(1, 1)).toBeCloseTo(1, 14);
        expect(interp.evaluate(1, 5)).toBeCloseTo(2, 14);
    });

    it('uses the weighted-average derivative formula', () => {
        // The slopes are (0,1,3,6). Prepending the extrapolated boundary
        // slopes, the window centered at sample index 2 is
        // (0, 1, 3, 6), so ad0 = |6 - 3| = 3, ad1 = |0 - 1| = 1 and the
        // derivative is (3*1 + 1*3)/(3 + 1) = 1.5.
        const F = [0, 0, 1, 4, 10];
        const interp = new IntpAkimaUniform1(F.length, 0, 1, F);
        expect(interp.evaluate(1, 2)).toBeCloseTo(1.5, 14);
    });

    it('is C1 continuous at the interior knots', () => {
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const interp = new IntpAkimaUniform1(F.length, 0, 1, F);
        const h = 1e-9;
        for (let i = 1; i < F.length - 1; ++i) {
            const left = interp.evaluate(1, i - h);
            const right = interp.evaluate(1, i + h);
            expect(left).toBeCloseTo(right, 6);
            expect(interp.evaluate(i - h)).toBeCloseTo(interp.evaluate(i + h), 6);
        }
    });

    it('matches finite differences of the function', () => {
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const interp = new IntpAkimaUniform1(F.length, 0, 1, F);
        const h = 1e-5;
        for (const x of [0.3, 1.7, 2.5, 3.1, 4.8, 5.4]) {
            const d1 = (interp.evaluate(x + h) - interp.evaluate(x - h)) / (2 * h);
            expect(interp.evaluate(1, x)).toBeCloseTo(d1, 6);
            const d2 = (interp.evaluate(1, x + h) - interp.evaluate(1, x - h)) / (2 * h);
            expect(interp.evaluate(2, x)).toBeCloseTo(d2, 5);
        }
    });

    it('clamps inputs outside the domain', () => {
        const F = [1, -3, 2.5, 7, 0];
        const xMin = -1, dx = 0.5;
        const interp = new IntpAkimaUniform1(F.length, xMin, dx, F);
        expect(interp.evaluate(-100)).toBeCloseTo(F[0], 14);
        expect(interp.evaluate(100)).toBeCloseTo(F[F.length - 1], 12);
        expect(interp.evaluate(1, -100)).toBe(interp.evaluate(1, xMin));
        expect(interp.evaluate(1, 100)).toBe(interp.evaluate(1, interp.getXMax()));
    });

    it('returns zero for derivative orders larger than three', () => {
        const F = [1, -3, 2.5, 7, 0];
        const interp = new IntpAkimaUniform1(5, 0, 1, F);
        expect(interp.evaluate(4, 1.5)).toBe(0);
        expect(interp.evaluate(7, 1.5)).toBe(0);
    });

    it('is cubic on each interval, matching the third derivative', () => {
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const interp = new IntpAkimaUniform1(F.length, 0, 1, F);
        // Within a single interval the polynomial has a constant third
        // derivative.
        const c3 = interp.evaluate(3, 2.25);
        expect(interp.evaluate(3, 2.75)).toBeCloseTo(c3, 14);
        // The second derivative varies linearly with x on the interval.
        const s0 = interp.evaluate(2, 2.25);
        const s1 = interp.evaluate(2, 2.75);
        expect(s1 - s0).toBeCloseTo(c3 * 0.5, 12);
    });

    it('handles the minimum sample count', () => {
        const F = [0, 2, 1];
        const interp = new IntpAkimaUniform1(3, 0, 1, F);
        expect(interp.evaluate(0)).toBeCloseTo(0, 14);
        expect(interp.evaluate(1)).toBeCloseTo(2, 14);
        expect(interp.evaluate(2)).toBeCloseTo(1, 14);
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V29): property-based checks against upstream
// IntpAkimaUniform1.h and the shared IntpAkima1.h base.
// ---------------------------------------------------------------------------

// A uniform Akima interpolator built from generated parameters. The spacing is
// bounded away from zero so that the divided differences and the finite
// differences below keep their significant digits.
const uniformAkima = fc.tuple(
    fc.integer({ min: -4, max: 4 }),                       // xMin
    fc.integer({ min: 1, max: 4 }).map(k => k / 2),        // xSpacing in [0.5,2]
    fc.array(fc.integer({ min: -20, max: 20 }),
        { minLength: 3, maxLength: 9 })                    // samples
).map(([xMin, xSpacing, F]) =>
    ({ xMin, xSpacing, F, intp: new IntpAkimaUniform1(F.length, xMin, xSpacing, F) }));

describe('IntpAkimaUniform1 verification', () => {
    // The polynomial of cell i has constant term F[i] and Lookup returns
    // dx = 0 there, so every sample but the last is reproduced bit-exactly.
    // The last sample is the right end of the last cell and is only
    // reproduced up to rounding.
    it('interpolates every sample', () => {
        check(uniformAkima, ({ xMin, xSpacing, F, intp }) => {
            for (let i = 0; i + 1 < F.length; ++i) {
                expect(intp.evaluate(xMin + xSpacing * i)).toBe(F[i]);
            }
            expectClose(intp.evaluate(intp.getXMax()), F[F.length - 1],
                1e-9, 1e-12);
            return true;
        });
    });

    // All divided differences of a linear sample set are equal, so
    // ComputeDerivative takes the slope[1] == slope[2] branch and the cubic
    // terms vanish.
    it('reproduces linear data, with zero second and third derivatives', () => {
        check(fc.tuple(
            fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: 1, max: 4 }).map(k => k / 2),
            fc.integer({ min: -9, max: 9 }),
            fc.integer({ min: -9, max: 9 }),
            fc.integer({ min: 3, max: 9 })
        ), ([xMin, xSpacing, a, b, n]) => {
            const F: number[] = [];
            for (let i = 0; i < n; ++i) {
                F.push(a + b * (xMin + xSpacing * i));
            }
            const intp = new IntpAkimaUniform1(n, xMin, xSpacing, F);
            const rand = seededRandom(n * 131 + 17);
            for (let k = 0; k < 8; ++k) {
                const x = intp.getXMin()
                    + rand() * (intp.getXMax() - intp.getXMin());
                expectClose(intp.evaluate(x), a + b * x, 1e-9, 1e-12);
                expectClose(intp.evaluate(1, x), b, 1e-9, 1e-12);
                expectClose(intp.evaluate(2, x), 0, 1e-9, 1e-12);
                expectClose(intp.evaluate(3, x), 0, 1e-9, 1e-12);
            }
            return true;
        });
    });

    // The interpolant is a cubic on each cell, so the fourth derivative is
    // zero and the second central difference is exact for it: the only error
    // is the round-off of the difference quotient, which is bounded by
    // eps*|F|/h^2 with h = spacing/4.
    it('matches the exact second central difference inside a cell', () => {
        check(uniformAkima, ({ xMin, xSpacing, F, intp }) => {
            const h = xSpacing / 4;
            for (let i = 0; i + 1 < F.length; ++i) {
                const x = xMin + xSpacing * (i + 0.5);
                const d2 = (intp.evaluate(x + h) - 2 * intp.evaluate(x)
                    + intp.evaluate(x - h)) / (h * h);
                expectClose(intp.evaluate(2, x), d2, 1e-8, 1e-8);
            }
            return true;
        });
    });

    // The fourth-order five-point stencil is exact for polynomials of degree
    // at most four, so it reproduces the first derivative of the cell cubic.
    it('matches the fourth-order first-derivative stencil inside a cell', () => {
        check(uniformAkima, ({ xMin, xSpacing, F, intp }) => {
            const h = xSpacing / 8;
            for (let i = 0; i + 1 < F.length; ++i) {
                const x = xMin + xSpacing * (i + 0.5);
                const d1 = (-intp.evaluate(x + 2 * h) + 8 * intp.evaluate(x + h)
                    - 8 * intp.evaluate(x - h) + intp.evaluate(x - 2 * h))
                    / (12 * h);
                expectClose(intp.evaluate(1, x), d1, 1e-8, 1e-8);
            }
            return true;
        });
    });

    // The third derivative of a cubic is constant, so it does not vary inside
    // a cell.
    it('has a constant third derivative inside each cell', () => {
        check(uniformAkima, ({ xMin, xSpacing, F, intp }) => {
            for (let i = 0; i + 1 < F.length; ++i) {
                const base = intp.evaluate(3, xMin + xSpacing * (i + 0.25));
                expect(intp.evaluate(3, xMin + xSpacing * (i + 0.75)))
                    .toBe(base);
            }
            return true;
        });
    });

    // The base class clamps x into [xMin, xMax] before the lookup, so an
    // out-of-domain query returns exactly the boundary value.
    it('clamps queries outside the domain', () => {
        check(fc.tuple(uniformAkima, finite(0, 50), fc.integer({ min: 0, max: 3 })),
            ([{ intp }, delta, order]) => {
                expect(intp.evaluate(order, intp.getXMin() - delta))
                    .toBe(intp.evaluate(order, intp.getXMin()));
                expect(intp.evaluate(order, intp.getXMax() + delta))
                    .toBe(intp.evaluate(order, intp.getXMax()));
                return true;
            });
    });

    // Polynomial::operator()(order, x) falls through its switch and returns
    // zero for orders above three (and, as a JavaScript-only case, for
    // negative orders).
    it('returns zero for derivative orders outside [0,3]', () => {
        check(fc.tuple(uniformAkima, fc.integer({ min: 4, max: 12 })),
            ([{ intp }, order]) => {
                const x = 0.5 * (intp.getXMin() + intp.getXMax());
                expect(intp.evaluate(order, x)).toBe(0);
                return true;
            });
    });

    // The samples are aliased, not copied (upstream stores Real const* mF),
    // but the polynomial coefficients are computed once in the constructor,
    // so a later mutation of F does not change the evaluations.
    it('exposes the aliased sample array through getF', () => {
        check(uniformAkima, ({ F, intp }) => {
            expect(intp.getF()).toBe(F);
            expect(intp.getQuantity()).toBe(F.length);
            return true;
        });
    });
});
