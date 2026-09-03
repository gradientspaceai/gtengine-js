import { describe, it, expect, afterEach } from 'vitest';
import { RootsPolynomial, type RootMultiplicity } from '../src/RootsPolynomial.js';

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

function evalPoly(c: readonly number[], t: number): number {
    let result = 0;
    for (let i = c.length - 1; i >= 0; --i) {
        result = t * result + c[i];
    }
    return result;
}

// Coefficients (ascending powers) of the monic polynomial with the given roots.
function coeffFromRoots(roots: readonly number[]): number[] {
    let c = [1];
    for (const r of roots) {
        const next = new Array<number>(c.length + 1).fill(0);
        for (let i = 0; i < c.length; ++i) {
            next[i] -= r * c[i];
            next[i + 1] += c[i];
        }
        c = next;
    }
    return c;
}

function expectRootsNear(actual: RootMultiplicity[],
    expected: Array<{ root: number; multiplicity: number }>, tol: number): void {
    expect(actual.length).toBe(expected.length);
    for (let i = 0; i < expected.length; ++i) {
        expect(Math.abs(actual[i].root - expected[i].root)).toBeLessThanOrEqual(tol);
        expect(actual[i].multiplicity).toBe(expected[i].multiplicity);
    }
}

describe('RootsPolynomial', () => {
    afterEach(() => {
        RootsPolynomial.rootsLowDegreeBlock = null;
    });

    describe('solveQuadratic', () => {
        it('finds two simple roots of (x-2)(x-3), sorted ascending', () => {
            const rm = RootsPolynomial.solveQuadratic(6, -5, 1);
            expectRootsNear(rm, [
                { root: 2, multiplicity: 1 },
                { root: 3, multiplicity: 1 }
            ], 1e-14);
        });

        it('finds the double root of (x-2)^2', () => {
            const rm = RootsPolynomial.solveQuadratic(4, -4, 1);
            expectRootsNear(rm, [{ root: 2, multiplicity: 2 }], 0);
        });

        it('returns no real roots for x^2 + 1', () => {
            expect(RootsPolynomial.solveQuadratic(1, 0, 1)).toEqual([]);
        });

        it('handles a non-monic quadratic 2x^2 - 8', () => {
            const rm = RootsPolynomial.solveQuadratic(-8, 0, 2);
            expectRootsNear(rm, [
                { root: -2, multiplicity: 1 },
                { root: 2, multiplicity: 1 }
            ], 0);
        });
    });

    describe('solveCubic', () => {
        it('finds three simple roots of (x-1)(x-2)(x-3), sorted ascending', () => {
            const rm = RootsPolynomial.solveCubic(-6, 11, -6, 1);
            expectRootsNear(rm, [
                { root: 1, multiplicity: 1 },
                { root: 2, multiplicity: 1 },
                { root: 3, multiplicity: 1 }
            ], 1e-12);
        });

        it('finds the double root of (x-1)^2 (x+2)', () => {
            // x^3 - 3x + 2, already depressed; the classifier delta is
            // exactly zero in double precision.
            const rm = RootsPolynomial.solveCubic(2, -3, 0, 1);
            expectRootsNear(rm, [
                { root: -2, multiplicity: 1 },
                { root: 1, multiplicity: 2 }
            ], 0);
        });

        it('finds the triple root of (x-1)^3', () => {
            const rm = RootsPolynomial.solveCubic(-1, 3, -3, 1);
            expectRootsNear(rm, [{ root: 1, multiplicity: 3 }], 0);
        });

        it('finds the single real root of x^3 + x + 1', () => {
            const rm = RootsPolynomial.solveCubic(1, 1, 0, 1);
            expectRootsNear(rm, [{ root: -0.6823278038280193, multiplicity: 1 }], 1e-12);
        });

        it('finds the single real root of x^3 - 8 (c1 = 0 branch)', () => {
            const rm = RootsPolynomial.solveCubic(-8, 0, 0, 1);
            expectRootsNear(rm, [{ root: 2, multiplicity: 1 }], 1e-14);
        });

        it('handles a root of zero (c0 = 0 branch)', () => {
            // x(x-1)(x+1) = x^3 - x.
            const rm = RootsPolynomial.solveCubic(0, -1, 0, 1);
            expectRootsNear(rm, [
                { root: -1, multiplicity: 1 },
                { root: 0, multiplicity: 1 },
                { root: 1, multiplicity: 1 }
            ], 0);
        });

        it('bumps the multiplicity of a zero root (x^3 branch)', () => {
            const rm = RootsPolynomial.solveCubic(0, 0, 0, 1);
            expectRootsNear(rm, [{ root: 0, multiplicity: 3 }], 0);
        });
    });

    describe('solveQuartic', () => {
        it('finds four simple roots of the biquadratic (x^2-1)(x^2-4)', () => {
            // x^4 - 5x^2 + 4.
            const rm = RootsPolynomial.solveQuartic(4, 0, -5, 0, 1);
            expectRootsNear(rm, [
                { root: -2, multiplicity: 1 },
                { root: -1, multiplicity: 1 },
                { root: 1, multiplicity: 1 },
                { root: 2, multiplicity: 1 }
            ], 1e-14);
        });

        it('finds two double roots of (x^2-1)^2', () => {
            // x^4 - 2x^2 + 1.
            const rm = RootsPolynomial.solveQuartic(1, 0, -2, 0, 1);
            expectRootsNear(rm, [
                { root: -1, multiplicity: 2 },
                { root: 1, multiplicity: 2 }
            ], 0);
        });

        it('finds the triple root of (x-1)^3 (x+3)', () => {
            // x^4 - 6x^2 + 8x - 3, already depressed with delta = 0, a0 = 0.
            const rm = RootsPolynomial.solveQuartic(-3, 8, -6, 0, 1);
            expectRootsNear(rm, [
                { root: -3, multiplicity: 1 },
                { root: 1, multiplicity: 3 }
            ], 0);
        });

        it('finds four simple roots of the general quartic (x-1)(x-2)(x-3)(x-5)', () => {
            // x^4 - 11x^3 + 41x^2 - 61x + 30.
            const rm = RootsPolynomial.solveQuartic(30, -61, 41, -11, 1);
            expectRootsNear(rm, [
                { root: 1, multiplicity: 1 },
                { root: 2, multiplicity: 1 },
                { root: 3, multiplicity: 1 },
                { root: 5, multiplicity: 1 }
            ], 1e-8);
        });

        it('finds the two real roots of x^4 - 1', () => {
            const rm = RootsPolynomial.solveQuartic(-1, 0, 0, 0, 1);
            expectRootsNear(rm, [
                { root: -1, multiplicity: 1 },
                { root: 1, multiplicity: 1 }
            ], 0);
        });

        it('returns no real roots for x^4 + 1', () => {
            expect(RootsPolynomial.solveQuartic(1, 0, 0, 0, 1)).toEqual([]);
        });

        it('finds the quadruple root of x^4', () => {
            const rm = RootsPolynomial.solveQuartic(0, 0, 0, 0, 1);
            expectRootsNear(rm, [{ root: 0, multiplicity: 4 }], 0);
        });
    });

    describe('getRootInfo*', () => {
        it('reports quadratic multiplicities', () => {
            expect(RootsPolynomial.getRootInfoQuadratic(6, -5, 1)).toEqual([1, 1]);
            expect(RootsPolynomial.getRootInfoQuadratic(4, -4, 1)).toEqual([2]);
            expect(RootsPolynomial.getRootInfoQuadratic(1, 0, 1)).toEqual([]);
        });

        it('reports cubic multiplicities', () => {
            expect(RootsPolynomial.getRootInfoCubic(-6, 11, -6, 1)).toEqual([1, 1, 1]);
            expect(RootsPolynomial.getRootInfoCubic(2, -3, 0, 1)).toEqual([1, 2]);
            expect(RootsPolynomial.getRootInfoCubic(-1, 3, -3, 1)).toEqual([3]);
            expect(RootsPolynomial.getRootInfoCubic(1, 1, 0, 1)).toEqual([1]);
        });

        it('reports quartic multiplicities', () => {
            expect(RootsPolynomial.getRootInfoQuartic(4, 0, -5, 0, 1)).toEqual([1, 1, 1, 1]);
            expect(RootsPolynomial.getRootInfoQuartic(1, 0, -2, 0, 1)).toEqual([2, 2]);
            expect(RootsPolynomial.getRootInfoQuartic(-3, 8, -6, 0, 1)).toEqual([3, 1]);
            expect(RootsPolynomial.getRootInfoQuartic(30, -61, 41, -11, 1)).toEqual([1, 1, 1, 1]);
            expect(RootsPolynomial.getRootInfoQuartic(-1, 0, 0, 0, 1)).toEqual([1, 1]);
            expect(RootsPolynomial.getRootInfoQuartic(1, 0, 0, 0, 1)).toEqual([]);
            expect(RootsPolynomial.getRootInfoQuartic(0, 0, 0, 0, 1)).toEqual([4]);
        });
    });

    describe('find (general degree)', () => {
        it('finds the roots of a quartic in ascending order', () => {
            const c = [30, -61, 41, -11, 1];
            const roots = RootsPolynomial.find(4, c, 2048);
            expect(roots.length).toBe(4);
            const expected = [1, 2, 3, 5];
            for (let i = 0; i < 4; ++i) {
                expect(Math.abs(roots[i] - expected[i])).toBeLessThanOrEqual(1e-8);
            }
        });

        it('handles degenerate degrees', () => {
            // Identically zero polynomial.
            expect(RootsPolynomial.find(2, [0, 0, 0], 64)).toEqual([0]);
            // Nonzero constant.
            expect(RootsPolynomial.find(0, [5], 64)).toEqual([]);
            // Linear.
            expect(RootsPolynomial.find(1, [-2, 1], 64)).toEqual([2]);
            // Degree reduction by trailing zero coefficients: 1*t + 0*t^2.
            // (The linear solver computes -c[0]/c[1] = -0, which equals 0.)
            const reduced = RootsPolynomial.find(2, [0, 1, 0], 64);
            expect(reduced.length).toBe(1);
            expect(reduced[0] === 0).toBe(true);
        });

        it('reports roots of odd multiplicity > 1 (residuals are zero)', () => {
            // (x-1)^3: every reported root must satisfy p(root) ~ 0. As in
            // upstream, a repeated root may be reported more than once
            // because subinterval endpoints have zero function values.
            const c = [-1, 3, -3, 1];
            const roots = RootsPolynomial.find(3, c, 2048);
            expect(roots.length).toBeGreaterThanOrEqual(1);
            for (const r of roots) {
                expect(Math.abs(r - 1)).toBeLessThanOrEqual(1e-5);
                expect(Math.abs(evalPoly(c, r))).toBeLessThanOrEqual(1e-14);
            }
        });

        it('finds a bracketed root by bisection (interval overload)', () => {
            // p(t) = t^2 - 2 on [1,2].
            const r = RootsPolynomial.find(2, [-2, 0, 1], 1, 2, 2048);
            expect(r.found).toBe(true);
            expect(Math.abs(r.root - Math.SQRT2)).toBeLessThanOrEqual(Number.EPSILON);
        });

        it('interval overload detects roots at the interval endpoints', () => {
            const c = [-4, 0, 1];  // t^2 - 4, exact root at t = 2
            let r = RootsPolynomial.find(2, c, 2, 3, 2048);
            expect(r).toEqual({ found: true, root: 2 });
            r = RootsPolynomial.find(2, c, 0, 2, 2048);
            expect(r).toEqual({ found: true, root: 2 });
        });

        it('interval overload rejects a non-bracketing interval', () => {
            const r = RootsPolynomial.find(2, [-2, 0, 1], 2, 3, 2048);
            expect(r.found).toBe(false);
        });

        it('interval overload rejects invalid inputs', () => {
            // Reversed interval that brackets a root sign-wise.
            expect(RootsPolynomial.find(2, [-2, 0, 1], 2, 1, 100).found).toBe(false);
            // Zero iteration budget.
            expect(RootsPolynomial.find(2, [-2, 0, 1], 1, 2, 0).found).toBe(false);
        });
    });

    describe('rootsLowDegreeBlock test hook', () => {
        it('reports the classification blocks taken by the solvers', () => {
            const blocks: number[] = [];
            RootsPolynomial.rootsLowDegreeBlock = (b) => blocks.push(b);

            RootsPolynomial.solveQuadratic(-1, 0, 1);   // two simple roots
            expect(blocks).toEqual([0]);

            blocks.length = 0;
            RootsPolynomial.solveQuadratic(0, 0, 1);    // double root
            expect(blocks).toEqual([1]);

            blocks.length = 0;
            RootsPolynomial.solveQuadratic(1, 0, 1);    // complex pair
            expect(blocks).toEqual([2]);

            blocks.length = 0;
            // (x-1)(x-2)(x-3) depresses to c0 = 0 exactly, so it takes the
            // depressed-quadratic path (block 0) then inserts the zero root
            // for the cubic (block 4).
            RootsPolynomial.solveCubic(-6, 11, -6, 1);
            expect(blocks).toEqual([0, 4]);

            blocks.length = 0;
            // (x-1)(x-2)(x-4): asymmetric roots, delta > 0 general path.
            RootsPolynomial.solveCubic(-8, 14, -7, 1);
            expect(blocks).toEqual([7]);

            blocks.length = 0;
            RootsPolynomial.solveQuartic(1, 0, -2, 0, 1);  // two double roots
            expect(blocks).toEqual([27]);
        });
    });

    describe('randomized cross-checks', () => {
        it('recovers randomized cubic roots', () => {
            const rng = makeRng(0x51a7e3);
            for (let trial = 0; trial < 200; ++trial) {
                let r: number[];
                do {
                    r = [rng() * 10 - 5, rng() * 10 - 5, rng() * 10 - 5];
                    r.sort((a, b) => a - b);
                } while (r[1] - r[0] < 0.25 || r[2] - r[1] < 0.25);

                const c = coeffFromRoots(r);
                const rm = RootsPolynomial.solveCubic(c[0], c[1], c[2], c[3]);
                expect(rm.length).toBe(3);
                for (let i = 0; i < 3; ++i) {
                    expect(rm[i].multiplicity).toBe(1);
                    expect(Math.abs(rm[i].root - r[i])).toBeLessThanOrEqual(1e-7);
                }
            }
        });

        it('recovers randomized quartic roots', () => {
            const rng = makeRng(0xbee5);
            for (let trial = 0; trial < 200; ++trial) {
                let r: number[];
                do {
                    r = [rng() * 8 - 4, rng() * 8 - 4, rng() * 8 - 4, rng() * 8 - 4];
                    r.sort((a, b) => a - b);
                } while (r[1] - r[0] < 0.5 || r[2] - r[1] < 0.5 || r[3] - r[2] < 0.5);

                const c = coeffFromRoots(r);
                const rm = RootsPolynomial.solveQuartic(c[0], c[1], c[2], c[3], c[4]);
                expect(rm.length).toBe(4);
                for (let i = 0; i < 4; ++i) {
                    expect(rm[i].multiplicity).toBe(1);
                    expect(Math.abs(rm[i].root - r[i])).toBeLessThanOrEqual(1e-6);
                }
            }
        });

        it('recovers randomized degree-5 roots with find', () => {
            const rng = makeRng(0xdead01);
            for (let trial = 0; trial < 50; ++trial) {
                let r: number[];
                do {
                    r = [0, 0, 0, 0, 0].map(() => rng() * 8 - 4);
                    r.sort((a, b) => a - b);
                } while (r[1] - r[0] < 0.5 || r[2] - r[1] < 0.5 ||
                    r[3] - r[2] < 0.5 || r[4] - r[3] < 0.5);

                const c = coeffFromRoots(r);
                const roots = RootsPolynomial.find(5, c, 2048);
                expect(roots.length).toBe(5);
                for (let i = 0; i < 5; ++i) {
                    expect(Math.abs(roots[i] - r[i])).toBeLessThanOrEqual(1e-8);
                }
            }
        });
    });
});
