import { describe, it, expect } from 'vitest';
import { BSRational } from '../src/BSRational.js';
import type { PolynomialRoot } from '../src/PolynomialRoot.js';
import { RootsGeneralPolynomial } from '../src/RootsGeneralPolynomial.js';
import { RootsPolynomial } from '../src/RootsPolynomial.js';
import { RootsCubic } from '../src/RootsCubic.js';
import { RootsQuartic } from '../src/RootsQuartic.js';
import { check, expectClose, fc, nonzero, wellScaled } from './helpers/arbitraries.js';

const rat = (x: number, y?: number) => BSRational.fromNumber(x, y);
const BOTH = [true, false];

function fromRoots(roots: readonly number[]): number[] {
    let p = [1];
    for (const r of roots) {
        const q = new Array<number>(p.length + 1).fill(0);
        for (let i = 0; i < p.length; ++i) {
            q[i + 1] += p[i];
            q[i] -= r * p[i];
        }
        p = q;
    }
    return p;
}

// Multiply the polynomial p by the quadratic x^2 + b*x + c.
function timesQuadratic(p: readonly number[], b: number, c: number): number[] {
    const q = new Array<number>(p.length + 2).fill(0);
    for (let i = 0; i < p.length; ++i) {
        q[i + 2] += p[i];
        q[i + 1] += b * p[i];
        q[i] += c * p[i];
    }
    return q;
}

function evaluate(p: readonly number[], x: number): number {
    let result = p[p.length - 1];
    for (let i = p.length - 2; i >= 0; --i) {
        result = x * result + p[i];
    }
    return result;
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function summarize(roots: readonly PolynomialRoot[]): [number, number][] {
    return roots.map(r => [r.x, r.m]);
}

function solveBoth(p: readonly number[]): PolynomialRoot[][] {
    return BOTH.map(useBisection =>
        RootsQuartic.solve(useBisection, p[0], p[1], p[2], p[3], p[4]));
}

describe('RootsQuartic.solve degree reduction and zero-valued roots', () => {
    it('falls through to the cubic solver when g4 is zero', () => {
        const p = fromRoots([1, 2, 3]);
        for (const useBisection of BOTH) {
            expect(summarize(RootsQuartic.solve(useBisection, p[0], p[1], p[2], p[3], 0)))
                .toEqual([[1, 1], [2, 1], [3, 1]]);
        }
    });

    it('handles zero-valued roots of every multiplicity', () => {
        for (const useBisection of BOTH) {
            // x^4.
            expect(summarize(RootsQuartic.solve(useBisection, 0, 0, 0, 0, 3)))
                .toEqual([[0, 4]]);
            // x^3 (x - 2).
            expect(summarize(RootsQuartic.solve(useBisection, 0, 0, 0, -2, 1)))
                .toEqual([[0, 3], [2, 1]]);
            // x^2 (x - 1)(x + 1) = x^4 - x^2.
            expect(summarize(RootsQuartic.solve(useBisection, 0, 0, -1, 0, 1)))
                .toEqual([[-1, 1], [0, 2], [1, 1]]);
            // x^2 (x^2 + 1): only the double zero is real.
            expect(summarize(RootsQuartic.solve(useBisection, 0, 0, 1, 0, 1)))
                .toEqual([[0, 2]]);
            // x (x - 1)(x - 2)(x - 3).
            const p = fromRoots([0, 1, 2, 3]);
            expect(summarize(RootsQuartic.solve(useBisection, p[0], p[1], p[2], p[3], p[4])))
                .toEqual([[0, 1], [1, 1], [2, 1], [3, 1]]);
        }
    });
});

describe('RootsQuartic.solve classification', () => {
    it('finds four distinct real roots in increasing order', () => {
        const p = fromRoots([1, 2, 3, 4]);
        for (const roots of solveBoth(p)) {
            expect(summarize(roots)).toEqual([[1, 1], [2, 1], [3, 1], [4, 1]]);
        }
    });

    it('finds two double roots', () => {
        const p = fromRoots([2, 2, -3, -3]);
        for (const roots of solveBoth(p)) {
            expect(summarize(roots)).toEqual([[-3, 2], [2, 2]]);
        }
    });

    it('finds a triple root plus a simple root', () => {
        for (const rs of [[1, 1, 1, 5], [-2, -2, -2, 1], [3, 3, 3, -4]]) {
            const p = fromRoots(rs);
            const expected = [...new Set(rs)].sort((a, b) => a - b);
            for (const roots of solveBoth(p)) {
                expect(roots.length).toBe(2);
                expect(roots.map(r => r.x)).toEqual(expected);
                const triple = roots.find(r => r.m === 3);
                expect(triple).toBeDefined();
                expect(triple!.x).toBe(rs[0]);
            }
        }
    });

    it('finds a quadruple root', () => {
        for (const r of [1, -5, 0.25]) {
            const p = fromRoots([r, r, r, r]);
            for (const roots of solveBoth(p)) {
                expect(summarize(roots)).toEqual([[r, 4]]);
            }
        }
    });

    it('finds a double root plus two simple roots', () => {
        const p = fromRoots([1, 1, 3, 7]);
        for (const roots of solveBoth(p)) {
            expect(summarize(roots)).toEqual([[1, 2], [3, 1], [7, 1]]);
        }
    });

    it('finds a double root when the other pair is complex', () => {
        // (x - 1)^2 (x^2 + 1).
        const p = timesQuadratic(fromRoots([1, 1]), 0, 1);
        for (const roots of solveBoth(p)) {
            expect(summarize(roots)).toEqual([[1, 2]]);
        }
    });

    it('finds two simple real roots when the other pair is complex', () => {
        // (x - 1)(x - 2)(x^2 + 1).
        const p = timesQuadratic(fromRoots([1, 2]), 0, 1);
        for (const roots of solveBoth(p)) {
            expect(roots.length).toBe(2);
            expect(roots[0].x).toBeCloseTo(1, 10);
            expect(roots[1].x).toBeCloseTo(2, 10);
            expect(roots.map(r => r.m)).toEqual([1, 1]);
        }
    });

    it('reports no roots for two complex-conjugate pairs', () => {
        // (x^2 + 1)(x^2 + 4) = x^4 + 5x^2 + 4.
        for (const roots of solveBoth([4, 0, 5, 0, 1])) {
            expect(roots).toEqual([]);
        }
        // (x^2 + 1)^2.
        for (const roots of solveBoth([1, 0, 2, 0, 1])) {
            expect(roots).toEqual([]);
        }
        // (x^2 - 2x + 2)(x^2 + 4x + 5), no real roots and d1 != 0.
        const p = timesQuadratic(timesQuadratic([1], -2, 2), 4, 5);
        for (const roots of solveBoth(p)) {
            expect(roots).toEqual([]);
        }
    });
});

describe('RootsQuartic biquadratic paths', () => {
    it('solves x^4 + d2*x^2 + d0 across all sign combinations', () => {
        for (const useBisection of BOTH) {
            // Four real roots: (x^2 - 1)(x^2 - 4).
            expect(summarize(RootsQuartic.solveDepressed(useBisection, 4, 0, -5)))
                .toEqual([[-2, 1], [-1, 1], [1, 1], [2, 1]]);
            // Two real roots, one complex pair: (x^2 - 2)(x^2 + 3).
            const two = RootsQuartic.solveDepressed(useBisection, -6, 0, 1);
            expect(two.length).toBe(2);
            expect(two[0].x).toBeCloseTo(-Math.SQRT2, 12);
            expect(two[1].x).toBeCloseTo(Math.SQRT2, 12);
            // Two double real roots: (x^2 - 4)^2.
            expect(summarize(RootsQuartic.solveDepressed(useBisection, 16, 0, -8)))
                .toEqual([[-2, 2], [2, 2]]);
            // One complex pair of multiplicity 2: (x^2 + 4)^2.
            expect(RootsQuartic.solveDepressed(useBisection, 16, 0, 8)).toEqual([]);
            // Two complex pairs, negative T: x^4 + x^2 + 1.
            expect(RootsQuartic.solveDepressed(useBisection, 1, 0, 1)).toEqual([]);
        }
    });
});

describe('RootsQuartic.solveMonic and solveDepressed', () => {
    it('agrees with the general solver on monic and depressed forms', () => {
        for (const rs of [[1, 2, 3, 4], [-1, -1, 2, 5], [0.5, 0.5, 0.5, 0.5]]) {
            const p = fromRoots(rs);
            for (const useBisection of BOTH) {
                const general = RootsQuartic.solve(useBisection, p[0], p[1], p[2], p[3], p[4]);
                const monic = RootsQuartic.solveMonic(useBisection, p[0], p[1], p[2], p[3]);
                expect(summarize(monic)).toEqual(summarize(general));
            }
        }

        // A genuinely depressed quartic: sum of roots is zero.
        const rs = [-3, -1, 1, 3];
        const p = fromRoots(rs);
        expect(p[3]).toBe(0);
        for (const useBisection of BOTH) {
            const depressed = RootsQuartic.solveDepressed(useBisection, p[0], p[1], p[2]);
            expect(depressed.map(r => r.m)).toEqual([1, 1, 1, 1]);
            for (let i = 0; i < 4; ++i) {
                expect(depressed[i].x).toBeCloseTo(rs[i], 12);
            }
        }
    });
});

describe('RootsQuartic rational instantiation', () => {
    it('classifies an exactly quadruple rational root', () => {
        // (x - 1/3)^4.
        const roots = RootsQuartic.solveMonicRational(true,
            rat(1, 81), rat(-4, 27), rat(6, 9), rat(-4, 3));
        expect(roots.length).toBe(1);
        expect(roots[0].m).toBe(4);
        expect(roots[0].x.equals(rat(1, 3))).toBe(true);
    });

    it('classifies an exact double root with two simple roots', () => {
        // (x - 1/2)^2 (x - 2)(x + 3) = the monic quartic below.
        const roots = RootsQuartic.solveMonicRational(false,
            rat(-3, 2), rat(25, 4), rat(-27, 4), rat(0));
        expect(roots.map(r => r.m)).toEqual([1, 2, 1]);
        expect(roots[0].x.equals(rat(-3))).toBe(true);
        expect(roots[1].x.equals(rat(1, 2))).toBe(true);
        expect(roots[2].x.equals(rat(2))).toBe(true);
    });

    it('agrees with the floating-point path on representable coefficients', () => {
        for (const useBisection of BOTH) {
            for (const rs of [[1, 2, 3, 4], [-2, -2, 1, 6], [0, 0, 1, -1], [1, 1, 1, 1]]) {
                const p = fromRoots(rs);
                const fp = RootsQuartic.solve(useBisection, p[0], p[1], p[2], p[3], p[4]);
                const r = RootsQuartic.solveRational(useBisection,
                    rat(p[0]), rat(p[1]), rat(p[2]), rat(p[3]), rat(p[4]));
                expect(r.map(v => v.m)).toEqual(fp.map(v => v.m));
                for (let i = 0; i < fp.length; ++i) {
                    expect(r[i].x.toNumber()).toBeCloseTo(fp[i].x, 12);
                }
            }
        }
    });

    it('exposes ComputeDepressedRoots directly', () => {
        for (const useBisection of BOTH) {
            // x^4 - 5x^2 + 4 = (x^2 - 1)(x^2 - 4).
            const roots = RootsQuartic.computeDepressedRoots(useBisection,
                rat(4), rat(0), rat(-5));
            expect(roots.map(r => [r.x.toNumber(), r.m]))
                .toEqual([[-2, 1], [-1, 1], [1, 1], [2, 1]]);
        }
    });
});

describe('RootsQuartic randomized cross-checks', () => {
    it('matches RootsGeneralPolynomial on random real-rooted quartics', () => {
        const random = makeRandom(0x9e3779b9);
        for (let trial = 0; trial < 40; ++trial) {
            const rs = [0, 0, 0, 0].map(() => 16 * random() - 8);
            rs.sort((a, b) => a - b);
            const p = fromRoots(rs);
            const general = RootsGeneralPolynomial.solve(p);
            expect(general.length).toBe(4);
            for (const useBisection of BOTH) {
                const mine = RootsQuartic.solve(useBisection, p[0], p[1], p[2], p[3], p[4]);
                expect(mine.length).toBe(4);
                for (let i = 0; i < 4; ++i) {
                    expect(mine[i].m).toBe(1);
                    expect(mine[i].x).toBeCloseTo(rs[i], 5);
                    expect(mine[i].x).toBeCloseTo(general[i], 5);
                }
            }
        }
    });

    it('matches RootsPolynomial.solveQuartic and keeps residuals small on random coefficients', () => {
        const random = makeRandom(0x13579bdf);
        let checked = 0;
        for (let trial = 0; trial < 80; ++trial) {
            const g4 = (random() < 0.5 ? 1 : -1) * (0.5 + random());
            const g3 = 6 * random() - 3;
            const g2 = 6 * random() - 3;
            const g1 = 6 * random() - 3;
            const g0 = 6 * random() - 3;
            const p = [g0, g1, g2, g3, g4];
            const legacy = RootsPolynomial.solveQuartic(g0, g1, g2, g3, g4);
            const bisect = RootsQuartic.solve(true, g0, g1, g2, g3, g4);
            const closed = RootsQuartic.solve(false, g0, g1, g2, g3, g4);
            expect(bisect.length).toBe(closed.length);
            expect(bisect.length).toBe(legacy.length);
            const derivative = [g1, 2 * g2, 3 * g3, 4 * g4];
            for (let i = 0; i < bisect.length; ++i) {
                expect(bisect[i].m).toBe(legacy[i].multiplicity);
                expect(bisect[i].x).toBeCloseTo(legacy[i].root, 5);
                expect(closed[i].x).toBeCloseTo(bisect[i].x, 8);
                const slope = Math.abs(evaluate(derivative, bisect[i].x));
                if (slope > 1e-2) {
                    expect(Math.abs(evaluate(p, bisect[i].x)) / slope).toBeLessThan(1e-9);
                    ++checked;
                }
            }
        }
        expect(checked).toBeGreaterThan(20);
    });

    it('finds exactly the real roots of random products of two quadratics', () => {
        const random = makeRandom(0x2468ace0);
        for (let trial = 0; trial < 40; ++trial) {
            // One real pair and one complex pair.
            const r0 = 8 * random() - 4;
            const r1 = 8 * random() - 4;
            const b = 4 * random() - 2;
            const c = b * b / 4 + 0.5 + random();
            const p = timesQuadratic(fromRoots([r0, r1]), b, c);
            const expected = [Math.min(r0, r1), Math.max(r0, r1)];
            for (const useBisection of BOTH) {
                const roots = RootsQuartic.solve(useBisection, p[0], p[1], p[2], p[3], p[4]);
                expect(roots.length).toBe(2);
                expect(roots[0].x).toBeCloseTo(expected[0], 5);
                expect(roots[1].x).toBeCloseTo(expected[1], 5);
            }
        }
    });

    it('resolves a near-double root pair without collapsing the multiplicities', () => {
        for (const e of [1e-6, 1e-8]) {
            const p = fromRoots([1, 1 + e, 3, 5]);
            for (const useBisection of BOTH) {
                const roots = RootsQuartic.solve(useBisection, p[0], p[1], p[2], p[3], p[4]);
                expect(roots.length).toBe(4);
                expect(roots.map(r => r.m)).toEqual([1, 1, 1, 1]);
                expect(roots[0].x).toBeCloseTo(1, 5);
                expect(roots[1].x).toBeCloseTo(1 + e, 5);
                expect(roots[2].x).toBeCloseTo(3, 8);
                expect(roots[3].x).toBeCloseTo(5, 8);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V45): property-based comparison against upstream
// RootsQuartic.h.
// ---------------------------------------------------------------------------

// Distinct roots of `rs` with their multiplicities, in increasing order: the
// exact classification upstream computes with rational arithmetic.
function classify(rs: readonly number[]): [number, number][] {
    const sorted = [...rs].sort((a, b) => a - b);
    const out: [number, number][] = [];
    for (const r of sorted) {
        if (out.length > 0 && out[out.length - 1][0] === r) {
            ++out[out.length - 1][1];
        } else {
            out.push([r, 1]);
        }
    }
    return out;
}

const smallInt = fc.integer({ min: -6, max: 6 });
const nonzeroSmallInt = fc.integer({ min: -3, max: 3 }).filter(v => v !== 0);
// An irreducible real quadratic x^2 + b*x + c (b^2 - 4c < 0) with small
// integer coefficients, so every product stays exact in binary64.
const irreducibleQuadratic = fc.tuple(fc.integer({ min: -4, max: 4 }),
    fc.integer({ min: 1, max: 10 })).filter(([b, c]) => b * b - 4 * c < 0);
// A power of two scales a double exactly, leaving ComputeClassifiers' exact
// rational output unchanged.
const powerOfTwo = fc.integer({ min: -20, max: 20 }).map(k => Math.pow(2, k));

describe('RootsQuartic verification', () => {
    it('reproduces the exact classification of integer-rooted quartics', () => {
        check(fc.tuple(smallInt, smallInt, smallInt, smallInt, nonzeroSmallInt),
            ([a, b, c, d, lead]) => {
                const rs = [a, b, c, d];
                // |coefficient| <= 3 * 6^4 is an exact integer in binary64.
                const p = fromRoots(rs).map(v => v * lead);
                const expected = classify(rs);
                for (const useBisection of BOTH) {
                    const roots = RootsQuartic.solve(useBisection, p[0], p[1], p[2],
                        p[3], p[4]);
                    expect(roots.length).toBe(expected.length);
                    for (let i = 0; i < roots.length; ++i) {
                        // The multiplicity comes from exact rational sign
                        // tests, so it must be exactly right.
                        expect(roots[i].m).toBe(expected[i][1]);
                        expectClose(roots[i].x, expected[i][0], 1e-8, 1e-8);
                    }
                }
            });
    });

    it('reports no real roots for a product of two irreducible quadratics', () => {
        // Upstream's positive-discriminant early-out tests only d2 > 0, which
        // lets quartics with two complex-conjugate pairs and d2 <= 0 fall
        // through to the four-real-root extraction (upstream bug, fixed in the
        // port). Every product of two irreducible quadratics is such a
        // quartic, so this property fails on the upstream branch condition.
        check(fc.tuple(irreducibleQuadratic, irreducibleQuadratic, nonzeroSmallInt,
            fc.boolean()),
            ([[b0, c0], [b1, c1], lead, useBisection]) => {
                const p = timesQuadratic(timesQuadratic([1], b0, c0), b1, c1)
                    .map(v => v * lead);
                expect(RootsQuartic.solve(useBisection, p[0], p[1], p[2], p[3], p[4]))
                    .toEqual([]);
            });
    });

    it('finds exactly the real pair when the other pair is complex', () => {
        check(fc.tuple(smallInt, smallInt, irreducibleQuadratic, nonzeroSmallInt,
            fc.boolean()),
            ([r0, r1, [b, c], lead, useBisection]) => {
                const p = timesQuadratic(fromRoots([r0, r1]), b, c).map(v => v * lead);
                const roots = RootsQuartic.solve(useBisection, p[0], p[1], p[2], p[3],
                    p[4]);
                const expected = classify([r0, r1]);
                expect(roots.length).toBe(expected.length);
                for (let i = 0; i < roots.length; ++i) {
                    expect(roots[i].m).toBe(expected[i][1]);
                    expectClose(roots[i].x, expected[i][0], 1e-8, 1e-8);
                }
            });
    });

    it('is invariant under exact (power-of-two) scaling of the coefficients', () => {
        check(fc.tuple(smallInt, smallInt, smallInt, smallInt, nonzeroSmallInt,
            powerOfTwo, fc.boolean()),
            ([a, b, c, d, lead, scale, useBisection]) => {
                const p = fromRoots([a, b, c, d]).map(v => v * lead);
                const base = RootsQuartic.solve(useBisection, p[0], p[1], p[2], p[3],
                    p[4]);
                const scaled = RootsQuartic.solve(useBisection, scale * p[0],
                    scale * p[1], scale * p[2], scale * p[3], scale * p[4]);
                expect(scaled.length).toBe(base.length);
                for (let i = 0; i < base.length; ++i) {
                    expect(scaled[i].m).toBe(base[i].m);
                    // ComputeClassifiers divides by the leading coefficient
                    // exactly, so the depressed quartic is identical.
                    expect(scaled[i].x === base[i].x).toBe(true);
                }
            });
    });

    it('agrees between the bisection and closed-form solvers', () => {
        check(fc.tuple(wellScaled(-4, 4), wellScaled(-4, 4), wellScaled(-4, 4),
            wellScaled(-4, 4), nonzero(-4, 4, 0.25)),
            ([g0, g1, g2, g3, g4]) => {
                const bisect = RootsQuartic.solve(true, g0, g1, g2, g3, g4);
                const closed = RootsQuartic.solve(false, g0, g1, g2, g3, g4);
                // The classification is exact and shared, so the count and the
                // multiplicities cannot differ between the two solvers.
                expect(closed.length).toBe(bisect.length);
                for (let i = 0; i < bisect.length; ++i) {
                    expect(closed[i].m).toBe(bisect[i].m);
                    // A root of multiplicity m is conditioned like eps^(1/m).
                    const tol = bisect[i].m === 1 ? 1e-6 : 1e-3;
                    expectClose(closed[i].x, bisect[i].x, tol, tol);
                }
            });
    });

    it('returns roots in increasing order, multiplicities summing to at most 4', () => {
        check(fc.tuple(wellScaled(-4, 4), wellScaled(-4, 4), wellScaled(-4, 4),
            wellScaled(-4, 4), wellScaled(-4, 4), fc.boolean()),
            ([g0, g1, g2, g3, g4, useBisection]) => {
                const roots = RootsQuartic.solve(useBisection, g0, g1, g2, g3, g4);
                let sum = 0;
                for (let i = 0; i < roots.length; ++i) {
                    expect(Number.isFinite(roots[i].x)).toBe(true);
                    expect(roots[i].m).toBeGreaterThanOrEqual(1);
                    sum += roots[i].m;
                    if (i > 0) { expect(roots[i - 1].x < roots[i].x).toBe(true); }
                }
                expect(sum).toBeLessThanOrEqual(4);
            });
    });

    it('keeps the backward error of every simple root near machine precision', () => {
        check(fc.tuple(nonzero(-3, 3, 0.1), nonzero(-3, 3, 0.1), nonzero(-3, 3, 0.1),
            nonzero(-3, 3, 0.1), nonzero(-3, 3, 0.5), fc.boolean()),
            ([g0, g1, g2, g3, g4, useBisection]) => {
                const p = [g0, g1, g2, g3, g4];
                const dp = [g1, 2 * g2, 3 * g3, 4 * g4];
                for (const r of RootsQuartic.solve(useBisection, g0, g1, g2, g3, g4)) {
                    if (r.m !== 1) { continue; }
                    const slope = Math.abs(evaluate(dp, r.x));
                    // |p(x)| / |p'(x)| is the first-order distance to the true
                    // root; only meaningful away from a nearly multiple root.
                    if (slope > 1e-2) {
                        expect(Math.abs(evaluate(p, r.x)) / slope).toBeLessThan(1e-8);
                    }
                }
            });
    });

    it('agrees across the general, monic and depressed entry points', () => {
        check(fc.tuple(wellScaled(-4, 4), wellScaled(-4, 4), wellScaled(-4, 4),
            wellScaled(-4, 4), fc.boolean()),
            ([m0, m1, m2, m3, useBisection]) => {
                const general = RootsQuartic.solve(useBisection, m0, m1, m2, m3, 1);
                const monic = RootsQuartic.solveMonic(useBisection, m0, m1, m2, m3);
                expect(summarize(monic)).toEqual(summarize(general));
                const depressed = RootsQuartic.solveDepressed(useBisection, m0, m1, m2);
                expect(summarize(depressed)).toEqual(
                    summarize(RootsQuartic.solveMonic(useBisection, m0, m1, m2, 0)));
            });
    });

    it('drops to the cubic solver when the leading coefficient is zero', () => {
        check(fc.tuple(wellScaled(-4, 4), wellScaled(-4, 4), wellScaled(-4, 4),
            wellScaled(-4, 4), fc.boolean()),
            ([g0, g1, g2, g3, useBisection]) => {
                expect(summarize(RootsQuartic.solve(useBisection, g0, g1, g2, g3, 0)))
                    .toEqual(summarize(RootsCubic.solve(useBisection, g0, g1, g2, g3)));
            });
    });

    it('matches the floating-point path on the rational instantiation', () => {
        check(fc.tuple(smallInt, smallInt, smallInt, smallInt, nonzeroSmallInt,
            fc.boolean()),
            ([a, b, c, d, lead, useBisection]) => {
                const p = fromRoots([a, b, c, d]).map(v => v * lead);
                const fp = RootsQuartic.solve(useBisection, p[0], p[1], p[2], p[3], p[4]);
                const r = RootsQuartic.solveRational(useBisection, rat(p[0]), rat(p[1]),
                    rat(p[2]), rat(p[3]), rat(p[4]));
                expect(r.map(v => v.m)).toEqual(fp.map(v => v.m));
                for (let i = 0; i < fp.length; ++i) {
                    expectClose(r[i].x.toNumber(), fp[i].x, 1e-12, 1e-12);
                }
            });
    });

    it('agrees with RootsGeneralPolynomial on square-free real-rooted quartics', () => {
        // RootsGeneralPolynomial brackets roots by sign changes and therefore
        // cannot see a root of even multiplicity, so the cross-check needs
        // four distinct roots.
        check(fc.tuple(smallInt, smallInt, smallInt, smallInt, fc.boolean())
            .filter(([a, b, c, d]) => new Set([a, b, c, d]).size === 4),
            ([a, b, c, d, useBisection]) => {
                const rs = [a, b, c, d].sort((x, y) => x - y);
                const p = fromRoots(rs);
                const general = RootsGeneralPolynomial.solve(p);
                const mine = RootsQuartic.solve(useBisection, p[0], p[1], p[2], p[3],
                    p[4]);
                expect(mine.length).toBe(4);
                expect(general.length).toBe(4);
                for (let i = 0; i < 4; ++i) {
                    expectClose(mine[i].x, general[i], 1e-6, 1e-6);
                }
            }, 100);
    });
});
