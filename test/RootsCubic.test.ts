import { describe, it, expect } from 'vitest';
import { BSRational } from '../src/BSRational.js';
import type { PolynomialRoot } from '../src/PolynomialRoot.js';
import { RootsCubic } from '../src/RootsCubic.js';
import { RootsGeneralPolynomial } from '../src/RootsGeneralPolynomial.js';
import { RootsPolynomial } from '../src/RootsPolynomial.js';
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

// Multiply the polynomial p by the irreducible quadratic x^2 + b*x + c.
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

describe('RootsCubic.solve degree reduction and zero-valued roots', () => {
    it('falls through to the quadratic solver when g3 is zero', () => {
        for (const useBisection of BOTH) {
            expect(summarize(RootsCubic.solve(useBisection, -6, 1, 1, 0)))
                .toEqual([[-3, 1], [2, 1]]);
        }
    });

    it('handles zero-valued roots of every multiplicity', () => {
        for (const useBisection of BOTH) {
            // x^3.
            expect(summarize(RootsCubic.solve(useBisection, 0, 0, 0, 2)))
                .toEqual([[0, 3]]);
            // x^2 (x - 4).
            expect(summarize(RootsCubic.solve(useBisection, 0, 0, -4, 1)))
                .toEqual([[0, 2], [4, 1]]);
            // x (x - 1)(x - 2) = x^3 - 3x^2 + 2x.
            expect(summarize(RootsCubic.solve(useBisection, 0, 2, -3, 1)))
                .toEqual([[0, 1], [1, 1], [2, 1]]);
            // x (x^2 + 1): one real root at zero.
            expect(summarize(RootsCubic.solve(useBisection, 0, 1, 0, 1)))
                .toEqual([[0, 1]]);
        }
    });
});

describe('RootsCubic.solve classification', () => {
    it('finds three distinct real roots in increasing order', () => {
        for (const useBisection of BOTH) {
            const p = fromRoots([1, 2, 3]);
            expect(summarize(RootsCubic.solve(useBisection, p[0], p[1], p[2], p[3])))
                .toEqual([[1, 1], [2, 1], [3, 1]]);
            const q = fromRoots([-5, -1, 4]);
            const qRoots = RootsCubic.solve(useBisection, q[0], q[1], q[2], q[3]);
            expect(qRoots.map(r => r.m)).toEqual([1, 1, 1]);
            expect(qRoots[0].x).toBeCloseTo(-5, 12);
            expect(qRoots[1].x).toBeCloseTo(-1, 12);
            expect(qRoots[2].x).toBeCloseTo(4, 12);
        }
    });

    it('finds a double root plus a simple root', () => {
        for (const useBisection of BOTH) {
            // (x - 1)^2 (x + 2).
            const p = fromRoots([1, 1, -2]);
            expect(summarize(RootsCubic.solve(useBisection, p[0], p[1], p[2], p[3])))
                .toEqual([[-2, 1], [1, 2]]);
            // (x + 3)^2 (x - 1): the double root sorts first.
            const q = fromRoots([-3, -3, 1]);
            expect(summarize(RootsCubic.solve(useBisection, q[0], q[1], q[2], q[3])))
                .toEqual([[-3, 2], [1, 1]]);
        }
    });

    it('finds a triple root', () => {
        for (const useBisection of BOTH) {
            const p = fromRoots([2, 2, 2]);
            expect(summarize(RootsCubic.solve(useBisection, p[0], p[1], p[2], p[3])))
                .toEqual([[2, 3]]);
            const q = fromRoots([-0.5, -0.5, -0.5]);
            expect(summarize(RootsCubic.solve(useBisection, q[0], q[1], q[2], q[3])))
                .toEqual([[-0.5, 3]]);
        }
    });

    it('finds the single real root when a complex pair is present', () => {
        for (const useBisection of BOTH) {
            // (x - 3)(x^2 + 1).
            const p = timesQuadratic(fromRoots([3]), 0, 1);
            const roots = RootsCubic.solve(useBisection, p[0], p[1], p[2], p[3]);
            expect(roots.length).toBe(1);
            expect(roots[0].m).toBe(1);
            expect(roots[0].x).toBeCloseTo(3, 12);
        }
    });
});

describe('RootsCubic.solveDepressed', () => {
    it('solves the pure cases d1 = 0 and d0 = 0', () => {
        for (const useBisection of BOTH) {
            // x^3 + d0 with d0 > 0 and d0 < 0.
            expect(RootsCubic.solveDepressed(useBisection, 8, 0)[0].x)
                .toBeCloseTo(-2, 12);
            expect(RootsCubic.solveDepressed(useBisection, -27, 0)[0].x)
                .toBeCloseTo(3, 12);
            // x^3 + d1*x = x (x^2 + d1).
            expect(summarize(RootsCubic.solveDepressed(useBisection, 0, 4)))
                .toEqual([[0, 1]]);
            const roots = RootsCubic.solveDepressed(useBisection, 0, -4);
            expect(roots.map(r => r.m)).toEqual([1, 1, 1]);
            expect(roots[0].x).toBeCloseTo(-2, 12);
            expect(roots[1].x).toBe(0);
            expect(roots[2].x).toBeCloseTo(2, 12);
            // x^3.
            expect(summarize(RootsCubic.solveDepressed(useBisection, 0, 0)))
                .toEqual([[0, 3]]);
        }
    });

    it('solves the Cardano example x^3 - 2x - 5', () => {
        for (const useBisection of BOTH) {
            const roots = RootsCubic.solveDepressed(useBisection, -5, -2);
            expect(roots.length).toBe(1);
            expect(roots[0].x).toBeCloseTo(2.0945514815423265, 12);
        }
    });

    it('bounds the single real root correctly outside [-1,1] (upstream bug fix)', () => {
        // Upstream ComputeDepressedRootsBisection uses b = max(1,|d0|,|d1|)
        // as a root bound. For x^3 - 0.9x - 0.9 that is b = 1, but the only
        // real root is about 1.2686, so upstream's bisection would return the
        // interval endpoint 1 as the "root". The port uses Cauchy's bound
        // 1 + max(|d0|,|d1|).
        const expected = 1.2686314944439885;
        for (const useBisection of BOTH) {
            const roots = RootsCubic.solveDepressed(useBisection, -0.9, -0.9);
            expect(roots.length).toBe(1);
            expect(roots[0].x).toBeCloseTo(expected, 12);
            expect(Math.abs(evaluate([-0.9, -0.9, 0, 1], roots[0].x)))
                .toBeLessThan(1e-14);
        }

        // Same failure mode mirrored about the origin.
        for (const useBisection of BOTH) {
            const roots = RootsCubic.solveDepressed(useBisection, 0.9, -0.9);
            expect(roots.length).toBe(1);
            expect(roots[0].x).toBeCloseTo(-expected, 12);
        }
    });
});

describe('RootsCubic rational instantiation', () => {
    it('classifies an exactly triple rational root that floating point misses', () => {
        // (x - 1/3)^3 = x^3 - x^2 + x/3 - 1/27.
        const roots = RootsCubic.solveMonicRational(true,
            rat(-1, 27), rat(1, 3), rat(-1));
        expect(roots.length).toBe(1);
        expect(roots[0].m).toBe(3);
        expect(roots[0].x.equals(rat(1, 3))).toBe(true);
    });

    it('classifies an exact double root with a rational simple root', () => {
        // (x - 1/2)^2 (x + 5/2).
        const roots = RootsCubic.solveMonicRational(false,
            rat(5, 8), rat(-9, 4), rat(3, 2));
        expect(roots.length).toBe(2);
        expect(roots.map(r => r.m)).toEqual([1, 2]);
        expect(roots[0].x.equals(rat(-5, 2))).toBe(true);
        expect(roots[1].x.equals(rat(1, 2))).toBe(true);
    });

    it('agrees with the floating-point path on representable coefficients', () => {
        for (const useBisection of BOTH) {
            for (const rs of [[1, 2, 3], [-4, -4, 2], [0.5, 0.5, 0.5], [0, 1, -1]]) {
                const p = fromRoots(rs);
                const fp = RootsCubic.solve(useBisection, p[0], p[1], p[2], p[3]);
                const r = RootsCubic.solveRational(useBisection,
                    rat(p[0]), rat(p[1]), rat(p[2]), rat(p[3]));
                expect(r.map(v => v.m)).toEqual(fp.map(v => v.m));
                for (let i = 0; i < fp.length; ++i) {
                    expect(r[i].x.toNumber()).toBeCloseTo(fp[i].x, 12);
                }
            }
        }
    });

    it('exposes ComputeDepressedRoots directly', () => {
        for (const useBisection of BOTH) {
            // x^3 - 3x + 2 = (x - 1)^2 (x + 2), delta = 0.
            const roots = RootsCubic.computeDepressedRoots(useBisection, rat(2), rat(-3));
            expect(roots.map(r => [r.x.toNumber(), r.m])).toEqual([[-2, 1], [1, 2]]);
        }
    });
});

describe('RootsCubic randomized cross-checks', () => {
    it('matches RootsGeneralPolynomial on random real-rooted cubics', () => {
        const random = makeRandom(0x0c0ffee1);
        for (let trial = 0; trial < 60; ++trial) {
            const rs = [20 * random() - 10, 20 * random() - 10, 20 * random() - 10];
            rs.sort((a, b) => a - b);
            const p = fromRoots(rs);
            const general = RootsGeneralPolynomial.solve(p);
            expect(general.length).toBe(3);
            for (const useBisection of BOTH) {
                const mine = RootsCubic.solve(useBisection, p[0], p[1], p[2], p[3]);
                expect(mine.length).toBe(3);
                for (let i = 0; i < 3; ++i) {
                    expect(mine[i].m).toBe(1);
                    expect(mine[i].x).toBeCloseTo(rs[i], 6);
                    expect(mine[i].x).toBeCloseTo(general[i], 6);
                }
            }
        }
    });

    it('matches RootsPolynomial.solveCubic and keeps residuals small on random coefficients', () => {
        const random = makeRandom(0x1234abcd);
        for (let trial = 0; trial < 120; ++trial) {
            const g3 = (random() < 0.5 ? 1 : -1) * (0.5 + random());
            const g2 = 8 * random() - 4;
            const g1 = 8 * random() - 4;
            const g0 = 8 * random() - 4;
            const p = [g0, g1, g2, g3];
            const legacy = RootsPolynomial.solveCubic(g0, g1, g2, g3);
            const bisect = RootsCubic.solve(true, g0, g1, g2, g3);
            const closed = RootsCubic.solve(false, g0, g1, g2, g3);
            expect(bisect.length).toBe(closed.length);
            // The rational classification is exact, so the number of distinct
            // real roots must agree with the legacy floating-point solver on
            // these well-separated random coefficients.
            expect(bisect.length).toBe(legacy.length);
            const derivative = [g1, 2 * g2, 3 * g3];
            for (let i = 0; i < bisect.length; ++i) {
                expect(bisect[i].m).toBe(legacy[i].multiplicity);
                expect(bisect[i].x).toBeCloseTo(legacy[i].root, 6);
                expect(closed[i].x).toBeCloseTo(bisect[i].x, 8);
                // Residual scaled by |p'(x)|: this is the backward error of a
                // simple root and must be near machine precision.
                const slope = Math.abs(evaluate(derivative, bisect[i].x));
                if (slope > 1e-3) {
                    expect(Math.abs(evaluate(p, bisect[i].x)) / slope).toBeLessThan(1e-10);
                }
            }
        }
    });

    it('finds one real root for random cubics with a complex-conjugate pair', () => {
        const random = makeRandom(0x77aabb11);
        for (let trial = 0; trial < 60; ++trial) {
            const r = 10 * random() - 5;
            // x^2 + b*x + c with b^2 - 4c < 0.
            const b = 4 * random() - 2;
            const c = b * b / 4 + 0.25 + random();
            const p = timesQuadratic(fromRoots([r]), b, c);
            for (const useBisection of BOTH) {
                const roots = RootsCubic.solve(useBisection, p[0], p[1], p[2], p[3]);
                expect(roots.length).toBe(1);
                expect(roots[0].m).toBe(1);
                expect(roots[0].x).toBeCloseTo(r, 6);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V45): property-based comparison against upstream
// RootsCubic.h.
// ---------------------------------------------------------------------------

// Distinct roots of `rs` with their multiplicities, in increasing order. This
// is the exact classification upstream computes with rational arithmetic, so
// the solver must reproduce it for every integer-rooted cubic.
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

const smallInt = fc.integer({ min: -8, max: 8 });
const nonzeroSmallInt = fc.integer({ min: -4, max: 4 }).filter(v => v !== 0);
// Exact scale factors: a power of two multiplies a double without rounding,
// so scaling the coefficients leaves the exact rational classifiers computed
// by ComputeClassifiers unchanged and the results must be bit-identical.
const powerOfTwo = fc.integer({ min: -20, max: 20 }).map(k => Math.pow(2, k));

describe('RootsCubic verification', () => {
    it('reproduces the exact classification of integer-rooted cubics', () => {
        check(fc.tuple(smallInt, smallInt, smallInt, nonzeroSmallInt),
            ([a, b, c, lead]) => {
                const rs = [a, b, c];
                // Integer coefficients of magnitude <= 4 * 8^3 are exact.
                const p = fromRoots(rs).map(v => v * lead);
                const expected = classify(rs);
                for (const useBisection of BOTH) {
                    const roots = RootsCubic.solve(useBisection, p[0], p[1], p[2], p[3]);
                    expect(roots.length).toBe(expected.length);
                    for (let i = 0; i < roots.length; ++i) {
                        // The multiplicity comes from an exact rational sign
                        // test, so it must be exactly right.
                        expect(roots[i].m).toBe(expected[i][1]);
                        expectClose(roots[i].x, expected[i][0], 1e-9, 1e-9);
                    }
                }
            });
    });

    it('is invariant under exact (power-of-two) scaling of the coefficients', () => {
        check(fc.tuple(smallInt, smallInt, smallInt, nonzeroSmallInt, powerOfTwo,
            fc.boolean()),
            ([a, b, c, lead, scale, useBisection]) => {
                const p = fromRoots([a, b, c]).map(v => v * lead);
                const base = RootsCubic.solve(useBisection, p[0], p[1], p[2], p[3]);
                const scaled = RootsCubic.solve(useBisection, scale * p[0], scale * p[1],
                    scale * p[2], scale * p[3]);
                expect(scaled.length).toBe(base.length);
                for (let i = 0; i < base.length; ++i) {
                    expect(scaled[i].m).toBe(base[i].m);
                    // ComputeClassifiers divides by the leading coefficient in
                    // exact rational arithmetic, so the depressed cubic and
                    // hence every root is identical, not merely close.
                    expect(scaled[i].x === base[i].x).toBe(true);
                }
            });
    });

    it('agrees between the bisection and closed-form solvers', () => {
        check(fc.tuple(wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
            nonzero(-5, 5, 0.25)),
            ([g0, g1, g2, g3]) => {
                const bisect = RootsCubic.solve(true, g0, g1, g2, g3);
                const closed = RootsCubic.solve(false, g0, g1, g2, g3);
                // Both paths share the exact rational classification, so the
                // count and the multiplicities cannot differ; only the root
                // estimates come from floating-point arithmetic.
                expect(closed.length).toBe(bisect.length);
                for (let i = 0; i < bisect.length; ++i) {
                    expect(closed[i].m).toBe(bisect[i].m);
                    // A root of multiplicity m is conditioned like eps^(1/m),
                    // so only simple roots get a tight tolerance.
                    const tol = bisect[i].m === 1 ? 1e-7 : 1e-4;
                    expectClose(closed[i].x, bisect[i].x, tol, tol);
                }
            });
    });

    it('returns roots in increasing order, multiplicities summing to at most 3', () => {
        check(fc.tuple(wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
            wellScaled(-5, 5), fc.boolean()),
            ([g0, g1, g2, g3, useBisection]) => {
                const roots = RootsCubic.solve(useBisection, g0, g1, g2, g3);
                let sum = 0;
                for (let i = 0; i < roots.length; ++i) {
                    expect(Number.isFinite(roots[i].x)).toBe(true);
                    expect(roots[i].m).toBeGreaterThanOrEqual(1);
                    sum += roots[i].m;
                    if (i > 0) { expect(roots[i - 1].x < roots[i].x).toBe(true); }
                }
                expect(sum).toBeLessThanOrEqual(3);
                // g3 = 0 degenerates to the quadratic (or linear) solver; a
                // genuine cubic always has at least one real root.
                if (g3 !== 0) { expect(roots.length).toBeGreaterThanOrEqual(1); }
            });
    });

    it('keeps the backward error of every simple root near machine precision', () => {
        check(fc.tuple(nonzero(-4, 4, 0.1), nonzero(-4, 4, 0.1), nonzero(-4, 4, 0.1),
            nonzero(-4, 4, 0.5), fc.boolean()),
            ([g0, g1, g2, g3, useBisection]) => {
                const p = [g0, g1, g2, g3];
                const dp = [g1, 2 * g2, 3 * g3];
                for (const r of RootsCubic.solve(useBisection, g0, g1, g2, g3)) {
                    if (r.m !== 1) { continue; }
                    const slope = Math.abs(evaluate(dp, r.x));
                    // |p(x)| / |p'(x)| is the first-order distance to the true
                    // root; it is meaningful only when the root is not nearly
                    // multiple (slope bounded away from zero).
                    if (slope > 1e-2) {
                        expect(Math.abs(evaluate(p, r.x)) / slope).toBeLessThan(1e-9);
                    }
                }
            });
    });

    it('agrees with RootsGeneralPolynomial on real-rooted cubics', () => {
        // RootsGeneralPolynomial brackets roots by sign changes, so it cannot
        // see a root of even multiplicity; restrict the cross-check to
        // square-free cubics.
        check(fc.tuple(smallInt, smallInt, smallInt, fc.boolean())
            .filter(([a, b, c]) => a !== b && b !== c && a !== c),
            ([a, b, c, useBisection]) => {
                const rs = [a, b, c].sort((x, y) => x - y);
                const p = fromRoots(rs);
                const general = RootsGeneralPolynomial.solve(p);
                const mine = RootsCubic.solve(useBisection, p[0], p[1], p[2], p[3]);
                // Both report distinct roots only, so the lists must have the
                // same length.
                expect(mine.length).toBe(general.length);
                for (let i = 0; i < mine.length; ++i) {
                    expectClose(mine[i].x, general[i], 1e-6, 1e-6);
                }
            }, 100);
    });

    it('finds exactly one real root when the other two are a complex pair', () => {
        check(fc.tuple(smallInt, fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: 1, max: 12 }), nonzeroSmallInt, fc.boolean())
            .filter(([, b, c]) => b * b - 4 * c < 0),
            ([r, b, c, lead, useBisection]) => {
                const p = timesQuadratic(fromRoots([r]), b, c).map(v => v * lead);
                const roots = RootsCubic.solve(useBisection, p[0], p[1], p[2], p[3]);
                expect(roots.length).toBe(1);
                expect(roots[0].m).toBe(1);
                expectClose(roots[0].x, r, 1e-9, 1e-9);
            });
    });

    it('agrees across the general, monic and depressed entry points', () => {
        check(fc.tuple(wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
            fc.boolean()),
            ([m0, m1, m2, useBisection]) => {
                const general = RootsCubic.solve(useBisection, m0, m1, m2, 1);
                const monic = RootsCubic.solveMonic(useBisection, m0, m1, m2);
                expect(summarize(monic)).toEqual(summarize(general));
                const depressed = RootsCubic.solveDepressed(useBisection, m0, m1);
                expect(summarize(depressed))
                    .toEqual(summarize(RootsCubic.solveMonic(useBisection, m0, m1, 0)));
            });
    });

    it('matches the floating-point path on the rational instantiation', () => {
        check(fc.tuple(smallInt, smallInt, smallInt, nonzeroSmallInt, fc.boolean()),
            ([a, b, c, lead, useBisection]) => {
                const p = fromRoots([a, b, c]).map(v => v * lead);
                const fp = RootsCubic.solve(useBisection, p[0], p[1], p[2], p[3]);
                const r = RootsCubic.solveRational(useBisection,
                    rat(p[0]), rat(p[1]), rat(p[2]), rat(p[3]));
                expect(r.map(v => v.m)).toEqual(fp.map(v => v.m));
                for (let i = 0; i < fp.length; ++i) {
                    // rootsRationalOps keeps the Rational root, rootsNumberOps
                    // rounds it once to a double, so the two agree to an ulp.
                    expectClose(r[i].x.toNumber(), fp[i].x, 1e-12, 1e-12);
                }
            });
    });

    it('never lets the bisection bound exclude the single real root', () => {
        // Upstream's b = max(1,|d0|,|d1|) is not a root bound for
        // x^3 + d1*x + d0; the port uses Cauchy's 1 + max(|d0|,|d1|). Sweep
        // the whole one-real-root region and check the residual rather than a
        // reference root, so the property is independent of any other solver.
        check(fc.tuple(wellScaled(-1, 1), wellScaled(-1, 0), fc.boolean())
            .filter(([d0, d1]) => -27 * d0 * d0 - 4 * d1 * d1 * d1 < 0),
            ([d0, d1, useBisection]) => {
                const roots = RootsCubic.solveDepressed(useBisection, d0, d1);
                expect(roots.length).toBe(1);
                const x = roots[0].x;
                const slope = Math.abs(3 * x * x + d1);
                if (slope > 1e-2) {
                    expect(Math.abs(evaluate([d0, d1, 0, 1], x)) / slope)
                        .toBeLessThan(1e-9);
                }
            });
    });
});
