import { describe, it, expect } from 'vitest';
import { BSRational } from '../src/BSRational.js';
import {
    PolynomialRootRational, RootsLinear, rootsNumberOps, rootsRationalOps, sortRoots
} from '../src/RootsLinear.js';
import {
    check, expectClose, fc, nonzero, wellScaled
} from './helpers/arbitraries.js';

describe('PolynomialRootRational', () => {
    it('defaults to the invalid root (0, 0)', () => {
        const r = new PolynomialRootRational();
        expect(r.x.toNumber()).toBe(0);
        expect(r.m).toBe(0);
    });

    it('compares using only the root estimate', () => {
        const a = new PolynomialRootRational(BSRational.fromNumber(1), 1);
        const b = new PolynomialRootRational(BSRational.fromNumber(1), 3);
        const c = new PolynomialRootRational(BSRational.fromNumber(2), 1);
        expect(a.equals(b)).toBe(true);
        expect(a.lessThan(b)).toBe(false);
        expect(a.lessThan(c)).toBe(true);
        expect(c.lessThan(a)).toBe(false);
    });
});

describe('sortRoots', () => {
    it('orders by the root estimate for both scalar types', () => {
        const numeric = [
            rootsNumberOps.makeRoot(3, 1),
            rootsNumberOps.makeRoot(-1, 2),
            rootsNumberOps.makeRoot(0, 1)
        ];
        sortRoots(rootsNumberOps, numeric);
        expect(numeric.map(r => r.x)).toEqual([-1, 0, 3]);
        expect(numeric.map(r => r.m)).toEqual([2, 1, 1]);

        const rational = [
            rootsRationalOps.makeRoot(BSRational.fromNumber(3), 1),
            rootsRationalOps.makeRoot(BSRational.fromNumber(-1), 2),
            rootsRationalOps.makeRoot(BSRational.fromNumber(0), 1)
        ];
        sortRoots(rootsRationalOps, rational);
        expect(rational.map(r => r.x.toNumber())).toEqual([-1, 0, 3]);
    });
});

describe('RootsLinear.solve', () => {
    it('reports no root when the degree is smaller than 1', () => {
        // g(x) = 0 has every x as a solution, g(x) = 5 has none. Both report
        // zero roots.
        expect(RootsLinear.solve(0, 0)).toEqual([]);
        expect(RootsLinear.solve(5, 0)).toEqual([]);
    });

    it('reports the zero-valued root exactly', () => {
        const roots = RootsLinear.solve(0, 7);
        expect(roots.length).toBe(1);
        expect(roots[0].x).toBe(0);
        expect(roots[0].m).toBe(1);
        // The zero is exact, not a rounding artifact of -0/7.
        expect(Object.is(roots[0].x, 0)).toBe(true);
    });

    it('solves g0 + g1*x = 0 as -g0/g1', () => {
        const cases: [number, number, number][] = [
            [-6, 3, 2],
            [6, 3, -2],
            [1, 4, -0.25],
            [-7, -2, -3.5],
            [1e-8, 1e8, -1e-16]
        ];
        for (const [g0, g1, expected] of cases) {
            const roots = RootsLinear.solve(g0, g1);
            expect(roots.length).toBe(1);
            expect(roots[0].m).toBe(1);
            expect(roots[0].x).toBeCloseTo(expected, 12);
            // The root satisfies the polynomial to machine precision.
            expect(Math.abs(g0 + g1 * roots[0].x)).toBeLessThanOrEqual(
                4 * Number.EPSILON * Math.abs(g0));
        }
    });
});

describe('RootsLinear.solveMonic', () => {
    it('solves m0 + x = 0 as -m0', () => {
        for (const m0 of [0, 1, -1, 3.5, -2.75, 1e100, -1e-100]) {
            const roots = RootsLinear.solveMonic(m0);
            expect(roots.length).toBe(1);
            expect(roots[0].m).toBe(1);
            expect(roots[0].x).toBe(-m0);
        }
    });
});

describe('RootsLinear rational instantiation', () => {
    const rat = (x: number, y?: number) => BSRational.fromNumber(x, y);

    it('mirrors the floating-point results exactly for representable inputs', () => {
        expect(RootsLinear.solveRational(rat(0), rat(0))).toEqual([]);
        expect(RootsLinear.solveRational(rat(5), rat(0))).toEqual([]);

        const zero = RootsLinear.solveRational(rat(0), rat(7));
        expect(zero.length).toBe(1);
        expect(zero[0].x.getSign()).toBe(0);
        expect(zero[0].m).toBe(1);

        const roots = RootsLinear.solveRational(rat(-6), rat(3));
        expect(roots.length).toBe(1);
        expect(roots[0].x.equals(rat(2))).toBe(true);
    });

    it('is exact where the floating-point path rounds', () => {
        // -1/3 is not representable in binary floating point, but the
        // rational path stores it exactly, so g0 + g1*x is exactly zero.
        const roots = RootsLinear.solveRational(rat(1), rat(3));
        expect(roots.length).toBe(1);
        expect(roots[0].x.equals(rat(-1, 3))).toBe(true);
        const residual = rat(1).add(rat(3).mul(roots[0].x));
        expect(residual.getSign()).toBe(0);

        // The floating-point path cannot be exact here.
        const fpRoot = RootsLinear.solve(1, 3)[0].x;
        expect(fpRoot).toBeCloseTo(-1 / 3, 15);
    });

    it('solves the monic polynomial with rational coefficients', () => {
        const roots = RootsLinear.solveMonicRational(rat(2, 5));
        expect(roots.length).toBe(1);
        expect(roots[0].m).toBe(1);
        expect(roots[0].x.equals(rat(-2, 5))).toBe(true);
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V44): property-based comparison against upstream
// RootsLinear.h.
// ---------------------------------------------------------------------------

describe('RootsLinear verification', () => {
    it('reports no root exactly when the leading coefficient is zero', () => {
        check(fc.tuple(wellScaled(-10, 10), wellScaled(-10, 10)), ([g0, g1]) => {
            const roots = RootsLinear.solve(g0, g1);
            if (g1 === 0) {
                // Upstream reports no roots whether g0 is zero (every x is a
                // solution) or not (no solution).
                expect(roots).toEqual([]);
            } else {
                expect(roots.length).toBe(1);
                expect(roots[0].m).toBe(1);
            }
        });
    });

    it('solves g0 + g1*x = 0 with a residual at the rounding level', () => {
        check(fc.tuple(wellScaled(-10, 10), nonzero(-10, 10, 1e-2)), ([g0, g1]) => {
            const x = RootsLinear.solve(g0, g1)[0].x;
            // -g0/g1 is a single correctly rounded division, so the residual
            // is bounded by one ulp of |g0| (plus the ulp of the product).
            expect(Math.abs(g0 + g1 * x))
                .toBeLessThanOrEqual(4 * Number.EPSILON * Math.abs(g0));
        });
    });

    it('returns the positive zero root when g0 is zero', () => {
        check(nonzero(-10, 10, 1e-2), g1 => {
            const roots = RootsLinear.solve(0, g1);
            expect(roots.length).toBe(1);
            expect(roots[0].m).toBe(1);
            // Upstream writes the literal T(0), not -g0/g1, so the sign of
            // the zero never depends on the sign of g1.
            expect(Object.is(roots[0].x, 0)).toBe(true);
        });
    });

    it('solves the monic polynomial exactly', () => {
        check(wellScaled(-10, 10), m0 => {
            const roots = RootsLinear.solveMonic(m0);
            expect(roots.length).toBe(1);
            expect(roots[0].m).toBe(1);
            expect(roots[0].x === -m0).toBe(true);
        });
    });

    it('solves the rational instantiation exactly', () => {
        check(fc.tuple(fc.integer({ min: -20, max: 20 }),
            fc.integer({ min: -20, max: 20 }).filter(v => v !== 0)),
            ([g0, g1]) => {
                const roots = RootsLinear.solveRational(BSRational.fromNumber(g0),
                    BSRational.fromNumber(g1));
                expect(roots.length).toBe(1);
                // BSRational division is exact, so g0 + g1*x is exactly zero.
                const residual = BSRational.fromNumber(g0).add(
                    BSRational.fromNumber(g1).mul(roots[0].x));
                expect(residual.getSign()).toBe(0);
            });
    });

    it('sorts roots by the estimate for both scalar types', () => {
        check(fc.array(fc.integer({ min: -9, max: 9 }),
            { minLength: 1, maxLength: 6 }), values => {
                const numeric = values.map(v => rootsNumberOps.makeRoot(v, 1));
                sortRoots(rootsNumberOps, numeric);
                expect(numeric.map(r => r.x))
                    .toEqual([...values].sort((a, b) => a - b));

                const rational = values.map(v =>
                    rootsRationalOps.makeRoot(BSRational.fromNumber(v), 1));
                sortRoots(rootsRationalOps, rational);
                expect(rational.map(r => r.x.toNumber()))
                    .toEqual([...values].sort((a, b) => a - b));
            }, 100);
    });

    it('keeps the two RootsScalarOps records mutually consistent', () => {
        // Both instantiations of the upstream template must agree on every
        // primitive the shared Solve bodies use.
        check(fc.tuple(wellScaled(-8, 8), nonzero(-8, 8, 1e-2)), ([a, b]) => {
            const ra = BSRational.fromNumber(a);
            const rb = BSRational.fromNumber(b);
            expect(rootsRationalOps.isZero(ra)).toBe(rootsNumberOps.isZero(a));
            expect(rootsRationalOps.lessThan(ra, rb))
                .toBe(rootsNumberOps.lessThan(a, b));
            // '+ 0' normalizes the signed zero: negating the double 0 gives
            // -0 while BSRational has no signed zero, so the rational
            // instantiation of upstream's 'roots[0] = { -m0, 1 }' reports +0
            // where the double instantiation reports -0. Both are faithful to
            // their upstream instantiation.
            expect(rootsRationalOps.negate(ra).toNumber() + 0)
                .toBe(rootsNumberOps.negate(a) + 0);
            expectClose(rootsRationalOps.div(ra, rb).toNumber(),
                rootsNumberOps.div(a, b), 1e-15, 1e-15);
            expect(rootsRationalOps.zero().getSign()).toBe(0);
            // toRational/fromRational round-trips a double exactly.
            expect(rootsNumberOps.fromRational(rootsNumberOps.toRational(a)) === a)
                .toBe(true);
            // rootsRationalOps clones rather than aliases, so a later mutation
            // of the result cannot corrupt the input.
            const clone = rootsRationalOps.toRational(ra);
            expect(clone).not.toBe(ra);
            expect(clone.equals(ra)).toBe(true);
        });
    });
});
