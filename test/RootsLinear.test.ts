import { describe, it, expect } from 'vitest';
import { BSRational } from '../src/BSRational.js';
import {
    PolynomialRootRational, RootsLinear, rootsNumberOps, rootsRationalOps, sortRoots
} from '../src/RootsLinear.js';

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
