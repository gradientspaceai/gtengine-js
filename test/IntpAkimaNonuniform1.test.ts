import { describe, it, expect } from 'vitest';
import { IntpAkimaNonuniform1 } from '../src/IntpAkimaNonuniform1.js';
import { IntpAkimaUniform1 } from '../src/IntpAkimaUniform1.js';
import { check, expectClose, fc, finite, seededRandom }
    from './helpers/arbitraries.js';

describe('IntpAkimaNonuniform1', () => {
    it('throws for invalid inputs', () => {
        expect(() => new IntpAkimaNonuniform1(2, [0, 1], [0, 1]))
            .toThrow('Invalid input to IntpAkima1 constructor.');
        // X is shorter than the sample count.
        expect(() => new IntpAkimaNonuniform1(4, [0, 1, 2], [0, 1, 2, 3]))
            .toThrow('Invalid input.');
        // X is not strictly increasing.
        expect(() => new IntpAkimaNonuniform1(3, [0, 1, 1], [0, 1, 2]))
            .toThrow('Invalid input.');
        expect(() => new IntpAkimaNonuniform1(3, [0, 2, 1], [0, 1, 2]))
            .toThrow('Invalid input.');
    });

    it('provides member access', () => {
        const X = [-1, 0.25, 3, 3.5, 10];
        const F = [0, 1, 4, 9, 16];
        const interp = new IntpAkimaNonuniform1(5, X, F);
        expect(interp.getQuantity()).toBe(5);
        expect(interp.getX()).toBe(X);
        expect(interp.getF()).toBe(F);
        expect(interp.getXMin()).toBe(-1);
        expect(interp.getXMax()).toBe(10);
    });

    it('interpolates the samples exactly', () => {
        const X = [-2, -1.5, 0, 0.25, 3, 7, 7.5];
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const interp = new IntpAkimaNonuniform1(X.length, X, F);
        for (let i = 0; i < X.length; ++i) {
            expect(interp.evaluate(X[i])).toBeCloseTo(F[i], 12);
            expect(interp.evaluate(0, X[i])).toBeCloseTo(F[i], 12);
        }
    });

    it('reproduces linear data exactly', () => {
        const X = [-2, -1.5, 0, 0.25, 3, 7, 7.5];
        const a = 1.5, b = -0.75;
        const F = X.map(x => a + b * x);
        const interp = new IntpAkimaNonuniform1(X.length, X, F);
        for (let k = 0; k <= 40; ++k) {
            const x = X[0] + (k / 40) * (X[X.length - 1] - X[0]);
            expect(interp.evaluate(x)).toBeCloseTo(a + b * x, 11);
            expect(interp.evaluate(1, x)).toBeCloseTo(b, 11);
            expect(interp.evaluate(2, x)).toBeCloseTo(0, 11);
            expect(interp.evaluate(3, x)).toBeCloseTo(0, 11);
        }
    });

    it('agrees with the uniform interpolator on uniformly spaced data', () => {
        const xMin = -1, dx = 0.5;
        const F = [1, -3, 2.5, 7, 0, -1.25, 4, 4, 6];
        const X: number[] = [];
        for (let i = 0; i < F.length; ++i) {
            X.push(xMin + i * dx);
        }
        const nonuniform = new IntpAkimaNonuniform1(F.length, X, F);
        const uniform = new IntpAkimaUniform1(F.length, xMin, dx, F);
        expect(nonuniform.getXMin()).toBeCloseTo(uniform.getXMin(), 14);
        expect(nonuniform.getXMax()).toBeCloseTo(uniform.getXMax(), 14);
        for (let k = 0; k <= 60; ++k) {
            const x = xMin + (k / 60) * (F.length - 1) * dx;
            for (let order = 0; order <= 3; ++order) {
                expect(nonuniform.evaluate(order, x))
                    .toBeCloseTo(uniform.evaluate(order, x), 10);
            }
        }
    });

    it('averages the two slopes when both slope pairs are equal', () => {
        // The x-values are nonuniform but the slopes are (1,1,1,2,2,2), so
        // at sample index 3 the estimated derivative is 0.5*(1 + 2) = 1.5.
        const X = [0, 1, 3, 4, 6, 6.5, 9];
        const slopes = [1, 1, 1, 2, 2, 2];
        const F = [0];
        for (let i = 0; i < slopes.length; ++i) {
            F.push(F[i] + slopes[i] * (X[i + 1] - X[i]));
        }
        const interp = new IntpAkimaNonuniform1(X.length, X, F);
        expect(interp.evaluate(1, X[3])).toBeCloseTo(1.5, 12);
        expect(interp.evaluate(1, X[1])).toBeCloseTo(1, 12);
        expect(interp.evaluate(1, X[5])).toBeCloseTo(2, 12);
    });

    it('is C1 continuous at the interior knots', () => {
        const X = [-2, -1.5, 0, 0.25, 3, 7, 7.5];
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const interp = new IntpAkimaNonuniform1(X.length, X, F);
        const h = 1e-9;
        for (let i = 1; i < X.length - 1; ++i) {
            expect(interp.evaluate(1, X[i] - h))
                .toBeCloseTo(interp.evaluate(1, X[i] + h), 6);
            expect(interp.evaluate(X[i] - h))
                .toBeCloseTo(interp.evaluate(X[i] + h), 6);
        }
    });

    it('matches finite differences of the function', () => {
        const X = [-2, -1.5, 0, 0.25, 3, 7, 7.5];
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const interp = new IntpAkimaNonuniform1(X.length, X, F);
        const h = 1e-5;
        for (const x of [-1.9, -0.4, 0.1, 1.7, 5.5, 7.2]) {
            const d1 = (interp.evaluate(x + h) - interp.evaluate(x - h)) / (2 * h);
            expect(interp.evaluate(1, x)).toBeCloseTo(d1, 5);
        }
    });

    it('clamps inputs outside the domain', () => {
        const X = [-2, -1.5, 0, 0.25, 3];
        const F = [1, -3, 2.5, 7, 0];
        const interp = new IntpAkimaNonuniform1(X.length, X, F);
        expect(interp.evaluate(-100)).toBeCloseTo(F[0], 12);
        expect(interp.evaluate(100)).toBeCloseTo(F[F.length - 1], 12);
        expect(interp.evaluate(1, -100)).toBe(interp.evaluate(1, X[0]));
        expect(interp.evaluate(1, 100)).toBe(interp.evaluate(1, X[X.length - 1]));
    });

    it('returns zero for derivative orders larger than three', () => {
        const X = [0, 1, 3, 4, 6];
        const F = [1, -3, 2.5, 7, 0];
        const interp = new IntpAkimaNonuniform1(5, X, F);
        expect(interp.evaluate(4, 2)).toBe(0);
    });

    it('handles the minimum sample count', () => {
        const X = [0, 1.5, 4];
        const F = [0, 2, 1];
        const interp = new IntpAkimaNonuniform1(3, X, F);
        for (let i = 0; i < 3; ++i) {
            expect(interp.evaluate(X[i])).toBeCloseTo(F[i], 12);
        }
    });

    it('is invariant under an affine change of the x-variable', () => {
        // Akima interpolation of {(x[i], f[i])} and {(a*x[i]+b, f[i])} give
        // the same function values at corresponding points.
        const X = [-2, -1.5, 0, 0.25, 3, 7, 7.5];
        const F = [1, -3, 2.5, 7, 0, -1.25, 4];
        const a = 2.5, b = -1;
        const Y = X.map(x => a * x + b);
        const interp0 = new IntpAkimaNonuniform1(X.length, X, F);
        const interp1 = new IntpAkimaNonuniform1(Y.length, Y, F);
        for (let k = 0; k <= 50; ++k) {
            const x = X[0] + (k / 50) * (X[X.length - 1] - X[0]);
            expect(interp1.evaluate(a * x + b)).toBeCloseTo(interp0.evaluate(x), 10);
            expect(a * interp1.evaluate(1, a * x + b))
                .toBeCloseTo(interp0.evaluate(1, x), 10);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V29): property-based checks against upstream
// IntpAkimaNonuniform1.h and the shared IntpAkima1.h base.
// ---------------------------------------------------------------------------

// Strictly increasing abscissae with gaps of at least 1/2, so that the
// divided differences keep their significant digits, paired with samples.
const nonuniformAkima = fc.tuple(
    fc.integer({ min: -6, max: 6 }),
    fc.array(fc.integer({ min: 1, max: 6 }), { minLength: 2, maxLength: 8 }),
    fc.array(fc.integer({ min: -20, max: 20 }), { minLength: 9, maxLength: 9 })
).map(([x0, gaps, samples]) => {
    const X = [x0];
    for (const g of gaps) {
        X.push(X[X.length - 1] + g / 2);
    }
    const F = samples.slice(0, X.length);
    return { X, F, intp: new IntpAkimaNonuniform1(X.length, X, F) };
});

describe('IntpAkimaNonuniform1 verification', () => {
    // The polynomial of cell i has constant term F[i] and Lookup returns
    // dx = 0 there, so every sample but the last is reproduced bit-exactly.
    it('interpolates every sample', () => {
        check(nonuniformAkima, ({ X, F, intp }) => {
            for (let i = 0; i + 1 < X.length; ++i) {
                expect(intp.evaluate(X[i])).toBe(F[i]);
            }
            expectClose(intp.evaluate(intp.getXMax()), F[F.length - 1],
                1e-9, 1e-12);
            expect(intp.getXMin()).toBe(X[0]);
            expect(intp.getXMax()).toBe(X[X.length - 1]);
            expect(intp.getX()).toBe(X);
            return true;
        });
    });

    // On a uniform grid the two interpolators must agree. They differ only in
    // the rounding of the divided differences: the nonuniform one divides by
    // X[i+1]-X[i] while the uniform one multiplies by 1/spacing.
    it('agrees with IntpAkimaUniform1 on a uniform grid', () => {
        check(fc.tuple(
            fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: 1, max: 4 }).map(k => k / 2),
            fc.array(fc.integer({ min: -20, max: 20 }),
                { minLength: 3, maxLength: 9 })
        ), ([xMin, xSpacing, F]) => {
            const n = F.length;
            const X: number[] = [];
            for (let i = 0; i < n; ++i) {
                X.push(xMin + xSpacing * i);
            }
            const nonuniform = new IntpAkimaNonuniform1(n, X, F);
            const uniform = new IntpAkimaUniform1(n, xMin, xSpacing, F);
            const rand = seededRandom(n * 977 + 5);
            for (let k = 0; k < 10; ++k) {
                const x = X[0] + rand() * (X[n - 1] - X[0]);
                for (let order = 0; order <= 2; ++order) {
                    expectClose(nonuniform.evaluate(order, x),
                        uniform.evaluate(order, x), 1e-9, 1e-10);
                }
            }
            return true;
        });
    });

    // All divided differences of a linear sample set are equal, so
    // ComputeDerivative takes the slope[1] == slope[2] branch and the cubic
    // terms vanish.
    it('reproduces linear data, with zero second and third derivatives', () => {
        check(fc.tuple(nonuniformAkima, fc.integer({ min: -9, max: 9 }),
            fc.integer({ min: -9, max: 9 })), ([{ X }, a, b]) => {
                const F = X.map(x => a + b * x);
                const intp = new IntpAkimaNonuniform1(X.length, X, F);
                const rand = seededRandom(X.length * 313 + 29);
                for (let k = 0; k < 8; ++k) {
                    const x = X[0] + rand() * (X[X.length - 1] - X[0]);
                    expectClose(intp.evaluate(x), a + b * x, 1e-9, 1e-12);
                    expectClose(intp.evaluate(1, x), b, 1e-9, 1e-12);
                    expectClose(intp.evaluate(2, x), 0, 1e-9, 1e-12);
                }
                return true;
            });
    });

    // The interpolant is a cubic on each cell, so the second central
    // difference is exact for it up to the round-off of the quotient.
    it('matches the exact second central difference inside a cell', () => {
        check(nonuniformAkima, ({ X, intp }) => {
            for (let i = 0; i + 1 < X.length; ++i) {
                const width = X[i + 1] - X[i];
                const h = width / 4;
                const x = X[i] + width / 2;
                const d2 = (intp.evaluate(x + h) - 2 * intp.evaluate(x)
                    + intp.evaluate(x - h)) / (h * h);
                expectClose(intp.evaluate(2, x), d2, 1e-7, 1e-7);
            }
            return true;
        });
    });

    // Akima interpolation commutes with an increasing affine change of the
    // independent variable: the divided differences scale by 1/alpha and the
    // derivative estimates are weighted averages of them, so the interpolant
    // is the same function reparameterized.
    it('is equivariant under an increasing affine change of x', () => {
        check(fc.tuple(nonuniformAkima,
            fc.integer({ min: 1, max: 8 }).map(k => k / 2),
            fc.integer({ min: -5, max: 5 })), ([{ X, F }, alpha, beta]) => {
                const Y = X.map(x => alpha * x + beta);
                const mapped = new IntpAkimaNonuniform1(Y.length, Y, F);
                const base = new IntpAkimaNonuniform1(X.length, X, F);
                const rand = seededRandom(X.length * 71 + 13);
                for (let k = 0; k < 8; ++k) {
                    const x = X[0] + rand() * (X[X.length - 1] - X[0]);
                    expectClose(mapped.evaluate(alpha * x + beta),
                        base.evaluate(x), 1e-9, 1e-11);
                    expectClose(alpha * mapped.evaluate(1, alpha * x + beta),
                        base.evaluate(1, x), 1e-9, 1e-11);
                }
                return true;
            });
    });

    // The base class clamps x into [xMin, xMax] before the lookup, so an
    // out-of-domain query returns exactly the boundary value.
    it('clamps queries outside the domain', () => {
        check(fc.tuple(nonuniformAkima, finite(0, 50),
            fc.integer({ min: 0, max: 3 })), ([{ intp }, delta, order]) => {
                expect(intp.evaluate(order, intp.getXMin() - delta))
                    .toBe(intp.evaluate(order, intp.getXMin()));
                expect(intp.evaluate(order, intp.getXMax() + delta))
                    .toBe(intp.evaluate(order, intp.getXMax()));
                return true;
            });
    });

    // Upstream asserts X[i+1] > X[i] for every i.
    it('rejects abscissae that are not strictly increasing', () => {
        check(fc.tuple(nonuniformAkima, fc.integer({ min: 0, max: 7 })),
            ([{ X, F }, which]) => {
                if (X.length < 4) {
                    return true;
                }
                const i = which % (X.length - 1);
                const bad = X.slice();
                bad[i + 1] = bad[i];
                expect(() => new IntpAkimaNonuniform1(bad.length, bad, F))
                    .toThrow('Invalid input.');
                return true;
            });
    });
});
