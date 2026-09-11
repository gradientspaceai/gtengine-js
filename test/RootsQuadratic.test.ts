import { describe, it, expect } from 'vitest';
import { BSRational } from '../src/BSRational.js';
import type { PolynomialRoot } from '../src/PolynomialRoot.js';
import { RootsPolynomial } from '../src/RootsPolynomial.js';
import { RootsLinear } from '../src/RootsLinear.js';
import { RootsQuadratic, rationalSqrtViaQuadratic } from '../src/RootsQuadratic.js';
import { check, fc, nonzero, wellScaled } from './helpers/arbitraries.js';

const rat = (x: number, y?: number) => BSRational.fromNumber(x, y);

// The coefficients, in increasing order of power, of the monic polynomial
// with the given roots.
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

function evaluate(p: readonly number[], x: number): number {
    let result = p[p.length - 1];
    for (let i = p.length - 2; i >= 0; --i) {
        result = x * result + p[i];
    }
    return result;
}

// A deterministic linear congruential generator so the randomized tests are
// reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

const BOTH = [true, false];

function summarize(roots: readonly PolynomialRoot[]): [number, number][] {
    return roots.map(r => [r.x, r.m]);
}

describe('RootsQuadratic.solve (general)', () => {
    it('falls through to the linear solver when g2 is zero', () => {
        for (const useBisection of BOTH) {
            expect(RootsQuadratic.solve(useBisection, -6, 3, 0)).toEqual(
                [{ x: 2, m: 1 }]);
            expect(RootsQuadratic.solve(useBisection, 1, 0, 0)).toEqual([]);
        }
    });

    it('finds two distinct real roots with the correct order', () => {
        for (const useBisection of BOTH) {
            // (x - 2)(x + 3) = x^2 + x - 6, scaled by 5.
            const roots = RootsQuadratic.solve(useBisection, -30, 5, 5);
            expect(summarize(roots)).toEqual([[-3, 1], [2, 1]]);
        }
    });

    it('finds a double root with multiplicity 2', () => {
        for (const useBisection of BOTH) {
            // (x - 2)^2 = x^2 - 4x + 4.
            expect(summarize(RootsQuadratic.solve(useBisection, 4, -4, 1)))
                .toEqual([[2, 2]]);
            // (x + 1.5)^2, scaled by -2.
            expect(summarize(RootsQuadratic.solve(useBisection, -4.5, -6, -2)))
                .toEqual([[-1.5, 2]]);
        }
    });

    it('reports no roots for a complex-conjugate pair', () => {
        for (const useBisection of BOTH) {
            expect(RootsQuadratic.solve(useBisection, 1, 0, 1)).toEqual([]);
            expect(RootsQuadratic.solve(useBisection, 5, 2, 1)).toEqual([]);
        }
    });

    it('handles zero-valued roots exactly', () => {
        for (const useBisection of BOTH) {
            // x^2 = 0.
            expect(summarize(RootsQuadratic.solve(useBisection, 0, 0, 3)))
                .toEqual([[0, 2]]);
            // x^2 - 5x = x (x - 5).
            expect(summarize(RootsQuadratic.solve(useBisection, 0, -5, 1)))
                .toEqual([[0, 1], [5, 1]]);
            // x^2 + 5x, so the nonzero root sorts before zero.
            expect(summarize(RootsQuadratic.solve(useBisection, 0, 5, 1)))
                .toEqual([[-5, 1], [0, 1]]);
        }
    });

    it('is scale invariant', () => {
        for (const useBisection of BOTH) {
            const reference = RootsQuadratic.solve(useBisection, -6, 1, 1);
            for (const s of [1e-6, 0.5, 3, 1e6]) {
                const roots = RootsQuadratic.solve(useBisection, -6 * s, s, s);
                expect(roots.length).toBe(2);
                expect(roots[0].x).toBeCloseTo(reference[0].x, 12);
                expect(roots[1].x).toBeCloseTo(reference[1].x, 12);
            }
        }
    });
});

describe('RootsQuadratic.solveMonic and solveDepressed', () => {
    it('agrees with the general solver on the monic polynomial', () => {
        for (const useBisection of BOTH) {
            for (const roots of [[1, 4], [-3, -3], [0, 2], [0, 0], [-1.5, 7.25]]) {
                const p = fromRoots(roots);
                const general = RootsQuadratic.solve(useBisection, p[0], p[1], p[2]);
                const monic = RootsQuadratic.solveMonic(useBisection, p[0], p[1]);
                expect(summarize(monic)).toEqual(summarize(general));
            }
        }
    });

    it('solves the depressed quadratic d0 + x^2 = 0', () => {
        for (const useBisection of BOTH) {
            expect(summarize(RootsQuadratic.solveDepressed(useBisection, -4)))
                .toEqual([[-2, 1], [2, 1]]);
            expect(summarize(RootsQuadratic.solveDepressed(useBisection, 0)))
                .toEqual([[0, 2]]);
            expect(RootsQuadratic.solveDepressed(useBisection, 4)).toEqual([]);

            // sqrt(2) to machine precision, and the roots are symmetric.
            const roots = RootsQuadratic.solveDepressed(useBisection, -2);
            expect(roots[1].x).toBeCloseTo(Math.SQRT2, 15);
            expect(roots[0].x).toBe(-roots[1].x);
        }
    });

    it('produces near-machine-precision square roots for small and large d0', () => {
        for (const useBisection of BOTH) {
            for (const d0 of [-1e-12, -0.3, -1, -7, -1e12]) {
                const roots = RootsQuadratic.solveDepressed(useBisection, d0);
                expect(roots.length).toBe(2);
                const expected = Math.sqrt(-d0);
                expect(Math.abs(roots[1].x - expected)).toBeLessThanOrEqual(
                    4 * Number.EPSILON * expected);
            }
        }
    });
});

describe('RootsQuadratic rational instantiation', () => {
    it('classifies the discriminant exactly where floating point cannot', () => {
        // x^2 - 2*(1/3)*x + 1/9 = (x - 1/3)^2. In floating point the
        // coefficients are inexact, but the rational path sees the exact
        // double root.
        const roots = RootsQuadratic.solveMonicRational(true, rat(1, 9), rat(-2, 3));
        expect(roots.length).toBe(1);
        expect(roots[0].m).toBe(2);
        expect(roots[0].x.equals(rat(1, 3))).toBe(true);
    });

    it('agrees with the floating-point path on representable coefficients', () => {
        for (const useBisection of BOTH) {
            for (const roots of [[2, -3], [1.5, 1.5], [0, 4], [-8, 8]]) {
                const p = fromRoots(roots);
                const fp = RootsQuadratic.solve(useBisection, p[0], p[1], p[2]);
                const r = RootsQuadratic.solveRational(useBisection,
                    rat(p[0]), rat(p[1]), rat(p[2]));
                expect(r.map(v => v.m)).toEqual(fp.map(v => v.m));
                for (let i = 0; i < fp.length; ++i) {
                    expect(r[i].x.toNumber()).toBeCloseTo(fp[i].x, 12);
                }
            }
        }
    });

    it('computes depressed roots identically through the class and free function', () => {
        for (const useBisection of BOTH) {
            const viaClass = RootsQuadratic.computeDepressedRoots(useBisection, rat(-9));
            expect(viaClass.map(v => v.x.toNumber())).toEqual([-3, 3]);
            expect(viaClass.map(v => v.m)).toEqual([1, 1]);
        }
    });
});

describe('rationalSqrtViaQuadratic', () => {
    it('returns the non-negative square root, or zero for non-positive input', () => {
        for (const useBisection of BOTH) {
            expect(rationalSqrtViaQuadratic(useBisection, rat(0)).getSign()).toBe(0);
            expect(rationalSqrtViaQuadratic(useBisection, rat(-3)).getSign()).toBe(0);
            expect(rationalSqrtViaQuadratic(useBisection, rat(16)).equals(rat(4))).toBe(true);
            expect(rationalSqrtViaQuadratic(useBisection, rat(2)).toNumber())
                .toBeCloseTo(Math.SQRT2, 15);
        }
    });
});

describe('RootsQuadratic randomized cross-checks', () => {
    it('matches RootsPolynomial.solveQuadratic on random coefficients', () => {
        const random = makeRandom(0x5eed1234);
        for (let trial = 0; trial < 200; ++trial) {
            const g2 = 4 * random() - 2 + (random() < 0.5 ? 0.5 : -0.5);
            const g1 = 8 * random() - 4;
            const g0 = 8 * random() - 4;
            const mine = RootsQuadratic.solve(true, g0, g1, g2);
            const closed = RootsQuadratic.solve(false, g0, g1, g2);
            const legacy = RootsPolynomial.solveQuadratic(g0, g1, g2);
            expect(mine.length).toBe(legacy.length);
            expect(closed.length).toBe(mine.length);
            for (let i = 0; i < mine.length; ++i) {
                expect(mine[i].m).toBe(legacy[i].multiplicity);
                expect(mine[i].x).toBeCloseTo(legacy[i].root, 8);
                expect(closed[i].x).toBeCloseTo(mine[i].x, 12);
            }
        }
    });

    it('produces residuals near machine precision on random real-rooted quadratics', () => {
        const random = makeRandom(0xbeef4321);
        for (let trial = 0; trial < 200; ++trial) {
            const r0 = 20 * random() - 10;
            const r1 = 20 * random() - 10;
            const p = fromRoots([r0, r1]);
            for (const useBisection of BOTH) {
                const roots = RootsQuadratic.solve(useBisection, p[0], p[1], p[2]);
                expect(roots.length).toBeGreaterThanOrEqual(1);
                const scale = 1 + Math.abs(p[0]) + Math.abs(p[1]);
                for (const root of roots) {
                    // The residual is scaled by the derivative magnitude,
                    // which is the natural conditioning of a simple root.
                    expect(Math.abs(evaluate(p, root.x))).toBeLessThan(1e-9 * scale);
                }
                const sorted = [Math.min(r0, r1), Math.max(r0, r1)];
                expect(roots[0].x).toBeCloseTo(sorted[0], 6);
                expect(roots[roots.length - 1].x).toBeCloseTo(sorted[1], 6);
            }
        }
    });

    it('resolves near-multiple roots better than the naive discriminant', () => {
        // (x - 1)(x - (1 + e)) for tiny e. The rational classification sees
        // the two distinct roots even though the double-precision
        // discriminant b^2 - 4ac underflows to zero.
        for (const e of [1e-9, 1e-10]) {
            const p = fromRoots([1, 1 + e]);
            for (const useBisection of BOTH) {
                const roots = RootsQuadratic.solve(useBisection, p[0], p[1], p[2]);
                expect(roots.length).toBe(2);
                expect(roots[0].m).toBe(1);
                expect(roots[1].m).toBe(1);
                expect(roots[0].x).toBeLessThan(roots[1].x);
                expect(roots[0].x).toBeCloseTo(1, 8);
                expect(roots[1].x).toBeCloseTo(1 + e, 8);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V44): property-based comparison against upstream
// RootsQuadratic.h.
// ---------------------------------------------------------------------------

// Distinct roots of `rs` with their multiplicities, in increasing order.
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

const smallInt = fc.integer({ min: -12, max: 12 });
const nonzeroSmallInt = fc.integer({ min: -5, max: 5 }).filter(v => v !== 0);
// Power-of-two scale factors multiply a double exactly.
const powerOfTwo = fc.integer({ min: -20, max: 20 }).map(k => Math.pow(2, k));

describe('RootsQuadratic verification', () => {
    it('reproduces the exact classification of integer-rooted quadratics', () => {
        check(fc.tuple(smallInt, smallInt, nonzeroSmallInt), ([a, b, lead]) => {
            const p = fromRoots([a, b]).map(v => v * lead);
            const expected = classify([a, b]);
            for (const useBisection of BOTH) {
                const roots = RootsQuadratic.solve(useBisection, p[0], p[1], p[2]);
                expect(roots.length).toBe(expected.length);
                for (let i = 0; i < roots.length; ++i) {
                    // The multiplicity comes from the exact rational sign of
                    // the depressed constant, so it must be exactly right.
                    expect(roots[i].m).toBe(expected[i][1]);
                    expect(Math.abs(roots[i].x - expected[i][0]))
                        .toBeLessThanOrEqual(1e-12 * (1 + Math.abs(expected[i][0])));
                }
            }
        });
    });

    it('reports no real roots exactly when the discriminant is negative', () => {
        check(fc.tuple(smallInt, smallInt, nonzeroSmallInt, fc.boolean()),
            ([g0, g1, g2, useBisection]) => {
                const roots = RootsQuadratic.solve(useBisection, g0, g1, g2);
                // The integer discriminant is exact in binary64.
                const disc = g1 * g1 - 4 * g0 * g2;
                if (disc < 0) {
                    expect(roots).toEqual([]);
                } else if (disc === 0) {
                    expect(roots.length).toBe(1);
                    expect(roots[0].m).toBe(2);
                } else {
                    expect(roots.length).toBe(2);
                    expect(roots.map(r => r.m)).toEqual([1, 1]);
                }
            });
    });

    it('is invariant under exact (power-of-two) scaling of the coefficients', () => {
        check(fc.tuple(smallInt, smallInt, nonzeroSmallInt, powerOfTwo, fc.boolean()),
            ([g0, g1, g2, scale, useBisection]) => {
                const base = RootsQuadratic.solve(useBisection, g0, g1, g2);
                const scaled = RootsQuadratic.solve(useBisection, scale * g0,
                    scale * g1, scale * g2);
                expect(scaled.length).toBe(base.length);
                for (let i = 0; i < base.length; ++i) {
                    expect(scaled[i].m).toBe(base[i].m);
                    // ComputeClassifiers divides by g2 in exact rational
                    // arithmetic, so the depressed quadratic is identical.
                    expect(scaled[i].x === base[i].x).toBe(true);
                }
            });
    });

    it('agrees between the bisection and closed-form solvers', () => {
        check(fc.tuple(wellScaled(-8, 8), wellScaled(-8, 8), nonzero(-8, 8, 0.25)),
            ([g0, g1, g2]) => {
                const bisect = RootsQuadratic.solve(true, g0, g1, g2);
                const closed = RootsQuadratic.solve(false, g0, g1, g2);
                expect(closed.length).toBe(bisect.length);
                for (let i = 0; i < bisect.length; ++i) {
                    expect(closed[i].m).toBe(bisect[i].m);
                    const tol = 1e-9 * (1 + Math.abs(bisect[i].x));
                    expect(Math.abs(closed[i].x - bisect[i].x))
                        .toBeLessThanOrEqual(tol);
                }
            });
    });

    it('returns roots in increasing order with multiplicities summing to at most 2', () => {
        check(fc.tuple(wellScaled(-8, 8), wellScaled(-8, 8), wellScaled(-8, 8),
            fc.boolean()),
            ([g0, g1, g2, useBisection]) => {
                const roots = RootsQuadratic.solve(useBisection, g0, g1, g2);
                let sum = 0;
                for (let i = 0; i < roots.length; ++i) {
                    expect(Number.isFinite(roots[i].x)).toBe(true);
                    sum += roots[i].m;
                    if (i > 0) { expect(roots[i - 1].x < roots[i].x).toBe(true); }
                }
                expect(sum).toBeLessThanOrEqual(2);
            });
    });

    it('keeps the backward error of every simple root near machine precision', () => {
        check(fc.tuple(nonzero(-8, 8, 0.1), nonzero(-8, 8, 0.1), nonzero(-8, 8, 0.5),
            fc.boolean()),
            ([g0, g1, g2, useBisection]) => {
                const p = [g0, g1, g2];
                for (const r of RootsQuadratic.solve(useBisection, g0, g1, g2)) {
                    if (r.m !== 1) { continue; }
                    const slope = Math.abs(g1 + 2 * g2 * r.x);
                    if (slope > 1e-2) {
                        expect(Math.abs(evaluate(p, r.x)) / slope).toBeLessThan(1e-10);
                    }
                }
            });
    });

    it('drops to the linear solver when the leading coefficient is zero', () => {
        check(fc.tuple(wellScaled(-8, 8), wellScaled(-8, 8), fc.boolean()),
            ([g0, g1, useBisection]) => {
                expect(summarize(RootsQuadratic.solve(useBisection, g0, g1, 0)))
                    .toEqual(summarize(RootsLinear.solve(g0, g1)));
            });
    });

    it('agrees across the general, monic and depressed entry points', () => {
        check(fc.tuple(wellScaled(-8, 8), wellScaled(-8, 8), fc.boolean()),
            ([m0, m1, useBisection]) => {
                expect(summarize(RootsQuadratic.solveMonic(useBisection, m0, m1)))
                    .toEqual(summarize(RootsQuadratic.solve(useBisection, m0, m1, 1)));
                expect(summarize(RootsQuadratic.solveDepressed(useBisection, m0)))
                    .toEqual(summarize(
                        RootsQuadratic.solveMonic(useBisection, m0, 0)));
            });
    });

    it('matches the floating-point path on the rational instantiation', () => {
        check(fc.tuple(smallInt, smallInt, nonzeroSmallInt, fc.boolean()),
            ([g0, g1, g2, useBisection]) => {
                const fp = RootsQuadratic.solve(useBisection, g0, g1, g2);
                const r = RootsQuadratic.solveRational(useBisection,
                    rat(g0), rat(g1), rat(g2));
                expect(r.map(v => v.m)).toEqual(fp.map(v => v.m));
                for (let i = 0; i < fp.length; ++i) {
                    expect(Math.abs(r[i].x.toNumber() - fp[i].x))
                        .toBeLessThanOrEqual(1e-12 * (1 + Math.abs(fp[i].x)));
                }
            });
    });

    it('extracts square roots through rationalSqrtViaQuadratic', () => {
        // The port's stand-in for upstream's repeated
        //   ComputeDepressedRoots(useBisection, -rV, rQRoots); rQRoots[1].x
        // idiom, which reads a stale array entry when rV <= 0.
        check(fc.tuple(fc.integer({ min: -20, max: 60 }), fc.boolean()),
            ([v, useBisection]) => {
                const s = rationalSqrtViaQuadratic(useBisection, rat(v));
                if (v <= 0) {
                    // No real square root (or zero): the port returns zero
                    // instead of whatever the previous call left behind.
                    expect(s.getSign()).toBe(0);
                } else {
                    expect(s.getSign()).toBeGreaterThan(0);
                    const approx = s.toNumber();
                    expect(Math.abs(approx * approx - v))
                        .toBeLessThanOrEqual(1e-12 * v);
                }
            }, 100);
    });

    it('bounds every depressed root by max(1,|d0|), the upstream bisection bound', () => {
        // The bisection path searches [0,b] with b = max(1,|d0|). The root of
        // x^2 + d0 is sqrt(-d0) <= b for every d0 <= 0, so unlike the cubic
        // the bound is valid and the search never returns an endpoint.
        check(fc.tuple(wellScaled(-50, 0), fc.boolean()), ([d0, useBisection]) => {
            const roots = RootsQuadratic.solveDepressed(useBisection, d0);
            if (d0 === 0) {
                expect(summarize(roots)).toEqual([[0, 2]]);
                return;
            }
            expect(roots.length).toBe(2);
            const b = Math.max(1, Math.abs(d0));
            expect(roots[1].x).toBeLessThanOrEqual(b);
            expect(roots[0].x === -roots[1].x).toBe(true);
            expect(Math.abs(roots[1].x * roots[1].x + d0))
                .toBeLessThanOrEqual(1e-12 * Math.abs(d0));
        });
    });
});
