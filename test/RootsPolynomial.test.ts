import { describe, it, expect, afterEach } from 'vitest';
import { RootsCubic } from '../src/RootsCubic.js';
import { RootsPolynomial, type RootMultiplicity } from '../src/RootsPolynomial.js';
import { RootsQuadratic } from '../src/RootsQuadratic.js';
import { RootsQuartic } from '../src/RootsQuartic.js';
import { check, fc, nonzero } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (V44): property-based comparison against upstream
// RootsPolynomial.h.
//
// Upstream templates the Solve*/GetRootInfo* family on a Rational type so the
// branch tests (delta, a0, a1, c0, c1, ...) are exact. The port instantiates
// the floating-point path, so the classification is subject to rounding. The
// properties below therefore pin the two things that must hold regardless of
// classification -- every reported root is a root of the polynomial, and the
// reported structure is self-consistent with GetRootInfo -- and characterize
// the rounding sensitivity instead of asserting exactness.
// ---------------------------------------------------------------------------

const rpLattice = fc.integer({ min: -5, max: 5 });
const rpLead = fc.integer({ min: -4, max: 4 }).filter(v => v !== 0);

function multiplicitySum(rm: readonly RootMultiplicity[]): number {
    return rm.reduce((s, r) => s + r.multiplicity, 0);
}

describe('RootsPolynomial verification', () => {
    it('solves quadratics with roots sorted, distinct and residual-free', () => {
        check(fc.tuple(rpLattice, rpLattice, rpLead), ([a, b, lead]) => {
            const c = coeffFromRoots([a, b]).map(v => v * lead);
            const rm = RootsPolynomial.solveQuadratic(c[0], c[1], c[2]);
            expect(multiplicitySum(rm)).toBeLessThanOrEqual(2);
            for (let i = 0; i < rm.length; ++i) {
                if (i > 0) { expect(rm[i - 1].root < rm[i].root).toBe(true); }
                const err = Math.min(Math.abs(rm[i].root - a),
                    Math.abs(rm[i].root - b));
                // A double root resolves only to about sqrt(eps).
                expect(err).toBeLessThanOrEqual(1e-7 * (1 + Math.abs(rm[i].root)));
            }
            // Integer roots make the classifier c0 = q0 - (q1/2)^2 exact, so
            // both roots are always found.
            expect(rm.length).toBe(a === b ? 1 : 2);
        }, 100);
    });

    it('solves cubics and quartics with every reported root a genuine root', () => {
        check(fc.tuple(fc.array(rpLattice, { minLength: 3, maxLength: 4 }), rpLead),
            ([rs, lead]) => {
                const c = coeffFromRoots(rs).map(v => v * lead);
                const rm = rs.length === 3
                    ? RootsPolynomial.solveCubic(c[0], c[1], c[2], c[3])
                    : RootsPolynomial.solveQuartic(c[0], c[1], c[2], c[3], c[4]);
                expect(multiplicitySum(rm)).toBeLessThanOrEqual(rs.length);
                for (let i = 0; i < rm.length; ++i) {
                    if (i > 0) { expect(rm[i - 1].root < rm[i].root).toBe(true); }
                    const err = Math.min(...rs.map(r => Math.abs(rm[i].root - r)));
                    // Repeated roots lose up to eps^(1/3) of accuracy, and the
                    // floating-point classifier can split a multiple root into
                    // nearby simple ones.
                    expect(err).toBeLessThanOrEqual(1e-4 * (1 + Math.abs(rm[i].root)));
                }
            }, 100);
    });

    it('agrees with getRootInfo on the number and multiplicities of the roots', () => {
        check(fc.tuple(fc.array(rpLattice, { minLength: 2, maxLength: 4 }), rpLead),
            ([rs, lead]) => {
                const c = coeffFromRoots(rs).map(v => v * lead);
                let rm: RootMultiplicity[];
                let info: number[];
                if (rs.length === 2) {
                    rm = RootsPolynomial.solveQuadratic(c[0], c[1], c[2]);
                    info = RootsPolynomial.getRootInfoQuadratic(c[0], c[1], c[2]);
                } else if (rs.length === 3) {
                    rm = RootsPolynomial.solveCubic(c[0], c[1], c[2], c[3]);
                    info = RootsPolynomial.getRootInfoCubic(c[0], c[1], c[2], c[3]);
                } else {
                    rm = RootsPolynomial.solveQuartic(c[0], c[1], c[2], c[3], c[4]);
                    info = RootsPolynomial.getRootInfoQuartic(c[0], c[1], c[2], c[3],
                        c[4]);
                }
                // The two families walk identical branch trees, so they agree
                // on the multiplicity multiset. The only way they can differ
                // is upstream's std::map::insert collapsing two roots that
                // round to the same double, which can only lose entries.
                expect(rm.length).toBeLessThanOrEqual(info.length);
                if (rm.length === info.length) {
                    expect(rm.map(r => r.multiplicity).sort((x, y) => x - y))
                        .toEqual([...info].sort((x, y) => x - y));
                }
                expect(info.reduce((s, m) => s + m, 0))
                    .toBeLessThanOrEqual(rs.length);
            }, 100);
    });

    it('agrees with the exactly classified low-degree solvers on distinct roots', () => {
        // Distinct integer roots keep the classifier quantities bounded away
        // from zero; a repeated root makes them exactly zero in rational
        // arithmetic and merely small in double, which is where the two
        // families part company (see the characterization test below).
        check(fc.tuple(fc.uniqueArray(rpLattice, { minLength: 2, maxLength: 4 }),
            rpLead),
            ([rs, lead]) => {
                const c = coeffFromRoots(rs).map(v => v * lead);
                const exact = rs.length === 2
                    ? RootsQuadratic.solve(true, c[0], c[1], c[2])
                    : (rs.length === 3
                        ? RootsCubic.solve(true, c[0], c[1], c[2], c[3])
                        : RootsQuartic.solve(true, c[0], c[1], c[2], c[3], c[4]));
                const rm = rs.length === 2
                    ? RootsPolynomial.solveQuadratic(c[0], c[1], c[2])
                    : (rs.length === 3
                        ? RootsPolynomial.solveCubic(c[0], c[1], c[2], c[3])
                        : RootsPolynomial.solveQuartic(c[0], c[1], c[2], c[3], c[4]));
                // Integer coefficients keep the depressed classifiers exact
                // enough that the deprecated solver matches the rational
                // classification of RootsQuadratic/Cubic/Quartic.
                expect(rm.length).toBe(exact.length);
                for (let i = 0; i < rm.length; ++i) {
                    expect(rm[i].multiplicity).toBe(exact[i].m);
                    expect(Math.abs(rm[i].root - exact[i].x))
                        .toBeLessThanOrEqual(1e-6 * (1 + Math.abs(exact[i].x)));
                }
            }, 100);
    });

    it('reports no real roots for a product of two irreducible quadratics', () => {
        check(fc.tuple(fc.integer({ min: -3, max: 3 }), fc.integer({ min: 1, max: 8 }),
            fc.integer({ min: -3, max: 3 }), fc.integer({ min: 1, max: 8 }), rpLead)
            .filter(([p1, p0, q1, q0]) =>
                p1 * p1 - 4 * p0 < 0 && q1 * q1 - 4 * q0 < 0),
            ([p1, p0, q1, q0, lead]) => {
                // (x^2 + p1 x + p0)(x^2 + q1 x + q0), scaled.
                const c = [p0 * q0, p1 * q0 + p0 * q1, p0 + q0 + p1 * q1, p1 + q1, 1]
                    .map(v => v * lead);
                expect(RootsPolynomial.solveQuartic(c[0], c[1], c[2], c[3], c[4]))
                    .toEqual([]);
                expect(RootsPolynomial.getRootInfoQuartic(c[0], c[1], c[2], c[3], c[4]))
                    .toEqual([]);
            }, 100);
    });

    it('finds the roots of a general polynomial within its Cauchy bound', () => {
        check(fc.tuple(fc.uniqueArray(rpLattice, { minLength: 1, maxLength: 5 }),
            rpLead),
            ([rs, lead]) => {
                const c = coeffFromRoots(rs).map(v => v * lead);
                const degree = c.length - 1;
                const roots = RootsPolynomial.find(degree, c, 1024) as number[];
                // The Cauchy bound must contain every root.
                const invLeading = 1 / c[degree];
                let maxValue = 0;
                for (let i = 0; i < degree; ++i) {
                    maxValue = Math.max(maxValue, Math.abs(c[i] * invLeading));
                }
                const bound = 1 + maxValue;
                for (const x of roots) {
                    expect(Math.abs(x)).toBeLessThanOrEqual(bound);
                    const err = Math.min(...rs.map(r => Math.abs(x - r)));
                    expect(err).toBeLessThanOrEqual(1e-6 * (1 + Math.abs(x)));
                }
                // Distinct roots all change the sign of p, so all are found.
                for (const r of rs) {
                    const err = Math.min(...roots.map(x => Math.abs(x - r)),
                        Number.POSITIVE_INFINITY);
                    expect(err).toBeLessThanOrEqual(1e-6 * (1 + Math.abs(r)));
                }
            }, 100);
    });

    it('handles the documented degenerate inputs of find', () => {
        check(fc.tuple(nonzero(-8, 8, 1e-2), fc.integer({ min: 0, max: 4 })),
            ([a, degree]) => {
                // A nonzero constant has no roots.
                expect(RootsPolynomial.find(0, [a], 64)).toEqual([]);
                // The identically zero polynomial reports the single root 0.
                expect(RootsPolynomial.find(degree,
                    new Array<number>(degree + 1).fill(0), 64)).toEqual([0]);
                // A negative degree, or too few coefficients, reports nothing.
                expect(RootsPolynomial.find(-1, [a], 64)).toEqual([]);
                expect(RootsPolynomial.find(degree + 2, [a], 64)).toEqual([]);
                // Trailing zero coefficients lower the effective degree.
                expect(RootsPolynomial.find(degree + 1, [-a, a,
                    ...new Array<number>(degree).fill(0)], 64)).toEqual([1]);
            }, 60);
    });

    it('respects the bracket contract of the bounded find', () => {
        check(fc.tuple(fc.uniqueArray(rpLattice, { minLength: 1, maxLength: 3 }),
            rpLead, fc.double({ min: 0.05, max: 0.9, noNaN: true })),
            ([rs, lead, frac]) => {
                const c = coeffFromRoots(rs).map(v => v * lead);
                const degree = c.length - 1;
                const r = rs[0];
                const lo = r - (1 - frac);
                const hi = r + frac;
                const out = RootsPolynomial.find(degree, c, lo, hi, 1024) as
                    { found: boolean; root: number };
                if (out.found) {
                    expect(Math.abs(evalPoly(c, out.root)))
                        .toBeLessThan(1e-8 * (1 + Math.abs(out.root)) ** degree);
                }
                // Reversed endpoints are rejected, but only after the
                // endpoint-value tests, which upstream performs first.
                const reversed = RootsPolynomial.find(degree, c, hi, lo, 1024) as
                    { found: boolean; root: number };
                if (evalPoly(c, hi) !== 0 && evalPoly(c, lo) !== 0) {
                    expect(reversed.found).toBe(false);
                    expect(reversed.root === 0).toBe(true);
                }
                // A zero iteration budget is rejected the same way.
                const noIterations = RootsPolynomial.find(degree, c, lo, hi, 0) as
                    { found: boolean; root: number };
                if (evalPoly(c, lo) !== 0 && evalPoly(c, hi) !== 0) {
                    expect(noIterations.found).toBe(false);
                }
            }, 100);
    });

    it('invokes the low-degree block callback for the branch it takes', () => {
        const seen: number[] = [];
        RootsPolynomial.rootsLowDegreeBlock = (block: number) => { seen.push(block); };
        try {
            // Two simple roots of the depressed quadratic: block 0.
            RootsPolynomial.solveQuadratic(-1, 0, 1);
            expect(seen).toContain(0);
            seen.length = 0;
            // One double root: block 1.
            RootsPolynomial.solveQuadratic(1, -2, 1);
            expect(seen).toContain(1);
            seen.length = 0;
            // A complex-conjugate pair: block 2.
            RootsPolynomial.solveQuadratic(1, 0, 1);
            expect(seen).toContain(2);
            seen.length = 0;
            // Three simple roots of the depressed cubic: block 7. The cubic
            // x^3 - 3x + 1 is already depressed, so c0 and c1 survive the
            // transformation unchanged and delta = 81 > 0.
            RootsPolynomial.solveCubic(1, -3, 0, 1);
            expect(seen).toContain(7);
            seen.length = 0;
            // x^3 - 3x^2 + 2x = x(x-1)(x-2) depresses to c0 = 0, so the
            // depressed quadratic runs (block 0) and the cubic adds the zero
            // root itself (block 4).
            RootsPolynomial.solveCubic(0, 2, -3, 1);
            expect(seen).toEqual([0, 4]);
        } finally {
            RootsPolynomial.rootsLowDegreeBlock = null;
        }
    });

    it('classifies repeated roots by rounding, not exactly (port deviation)', () => {
        // Upstream templates the classification on a Rational type, so the
        // branch tests are exact. This port instantiates Rational = double.
        //
        // (1) -x (x + 1)^2 has the exact classification { (-1, 2), (0, 1) }.
        // Depressing it gives c1 = -1/3 and c0 = -2/27 whose rounded double
        // values make delta = -(4 c1^3 + 27 c0^2) a tiny nonzero instead of
        // exactly zero, so the solver takes the "one simple root" branch and
        // loses the double root at -1 completely.
        const cubic = RootsPolynomial.solveCubic(0, -1, -2, -1);
        expect(cubic.length).toBe(1);
        expect(cubic[0].multiplicity).toBe(1);
        expect(Math.abs(cubic[0].root)).toBeLessThan(1e-15);
        const exactCubic = RootsCubic.solve(true, 0, -1, -2, -1);
        expect(exactCubic.map(r => [r.x, r.m])).toEqual([[-1, 2], [0, 1]]);

        // (2) The reverse error. For the coefficients (1, -2/3, 1/9) as
        // doubles, the exact discriminant is negative -- there are no real
        // roots -- but the depressed constant c0 = q0 - (q1/2)^2 rounds to
        // exactly zero, so the solver reports a double root at 3.
        const quadratic = RootsPolynomial.solveQuadratic(1, -2 / 3, 1 / 9);
        expect(quadratic.length).toBe(1);
        expect(quadratic[0].multiplicity).toBe(2);
        expect(quadratic[0].root).toBeCloseTo(3, 9);
        expect(RootsQuadratic.solve(true, 1, -2 / 3, 1 / 9)).toEqual([]);

        // This is the conditioning that dependent files inherit from the
        // deprecated solver; RootsLinear/Quadratic/Cubic/Quartic classify
        // with exact rational arithmetic and do not have it.
    });
});
