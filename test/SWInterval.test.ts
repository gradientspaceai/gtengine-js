import { describe, it, expect } from 'vitest';
import { check, fc } from './helpers/arbitraries.js';
import { SWInterval, nextDown, nextUp } from '../src/SWInterval.js';

// Deterministic pseudorandom generator (LCG) for reproducible tests.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function expectContains(w: SWInterval, value: number): void {
    expect(w.get(0)).toBeLessThanOrEqual(value);
    expect(w.get(1)).toBeGreaterThanOrEqual(value);
}

describe('nextDown/nextUp (std::nextafter toward -max/+max)', () => {
    it('steps one ulp around 1', () => {
        expect(nextUp(1)).toBe(1 + Math.pow(2, -52));
        expect(nextDown(1)).toBe(1 - Math.pow(2, -53));
        expect(nextDown(nextUp(1))).toBe(1);
        expect(nextUp(nextDown(1))).toBe(1);
    });

    it('steps one ulp for negative values', () => {
        expect(nextUp(-1)).toBe(-(1 - Math.pow(2, -53)));
        expect(nextDown(-1)).toBe(-(1 + Math.pow(2, -52)));
    });

    it('handles zeros and subnormals', () => {
        expect(nextUp(0)).toBe(Number.MIN_VALUE);
        expect(nextDown(0)).toBe(-Number.MIN_VALUE);
        expect(nextUp(-0)).toBe(Number.MIN_VALUE);
        expect(nextDown(-0)).toBe(-Number.MIN_VALUE);
        expect(nextDown(Number.MIN_VALUE)).toBe(0);
        expect(nextUp(-Number.MIN_VALUE)).toBe(-0);
        expect(nextUp(Number.MIN_VALUE)).toBe(2 * Number.MIN_VALUE);
    });

    it('handles infinities and MAX_VALUE like std::nextafter(x, +/-max)', () => {
        // Stepping from an infinity toward +/-max returns +/-MAX_VALUE.
        expect(nextDown(Number.POSITIVE_INFINITY)).toBe(Number.MAX_VALUE);
        expect(nextUp(Number.NEGATIVE_INFINITY)).toBe(-Number.MAX_VALUE);
        expect(nextUp(Number.POSITIVE_INFINITY)).toBe(Number.MAX_VALUE);
        expect(nextDown(Number.NEGATIVE_INFINITY)).toBe(-Number.MAX_VALUE);
        // from == to returns to.
        expect(nextUp(Number.MAX_VALUE)).toBe(Number.MAX_VALUE);
        expect(nextDown(-Number.MAX_VALUE)).toBe(-Number.MAX_VALUE);
        // Stepping away from +/-MAX_VALUE in the other direction is normal.
        expect(nextDown(Number.MAX_VALUE)).toBeLessThan(Number.MAX_VALUE);
        expect(nextUp(-Number.MAX_VALUE)).toBeGreaterThan(-Number.MAX_VALUE);
    });

    it('propagates NaN', () => {
        expect(Number.isNaN(nextUp(Number.NaN))).toBe(true);
        expect(Number.isNaN(nextDown(Number.NaN))).toBe(true);
    });

    it('is monotone: nextDown(x) < x < nextUp(x) for finite nonextreme x', () => {
        const rand = makeRandom(0x1234);
        for (let trial = 0; trial < 200; ++trial) {
            const x = (rand() - 0.5) * Math.pow(10, Math.floor(rand() * 40) - 20);
            expect(nextDown(x)).toBeLessThan(x);
            expect(nextUp(x)).toBeGreaterThan(x);
        }
    });
});

describe('SWInterval construction', () => {
    it('default-constructs to [0, 0]', () => {
        const v = new SWInterval();
        expect(v.getEndpoints()).toEqual([0, 0]);
    });

    it('creates the degenerate interval [e, e]', () => {
        const v = new SWInterval(3.5);
        expect(v.getEndpoints()).toEqual([3.5, 3.5]);
    });

    it('creates [e0, e1] and reads endpoints via get', () => {
        const v = new SWInterval(-1, 2);
        expect(v.get(0)).toBe(-1);
        expect(v.get(1)).toBe(2);
    });

    it('fromEndpoints and clone copy the endpoints', () => {
        const v = SWInterval.fromEndpoints([1, 4]);
        expect(v.getEndpoints()).toEqual([1, 4]);
        const c = v.clone();
        expect(c.getEndpoints()).toEqual([1, 4]);
        expect(c).not.toBe(v);
    });

    it('throwOnInvalid traps e0 > e1 when enabled', () => {
        expect(() => new SWInterval(2, 1)).not.toThrow();
        SWInterval.throwOnInvalid = true;
        try {
            expect(() => new SWInterval(2, 1)).toThrow('Invalid SWInterval.');
            expect(() => new SWInterval(1, 2)).not.toThrow();
        } finally {
            SWInterval.throwOnInvalid = false;
        }
    });

    it('exposes the upstream constants', () => {
        expect(SWInterval.zero).toBe(0);
        expect(SWInterval.one).toBe(1);
        expect(SWInterval.max).toBe(Number.MAX_VALUE);
        expect(SWInterval.inf).toBe(Number.POSITIVE_INFINITY);
    });
});

describe('SWInterval leaf-node operations', () => {
    it('add(u, v) encloses the exact sum', () => {
        const w = SWInterval.add(0.1, 0.2);
        expect(w.get(0)).toBe(nextDown(0.1 + 0.2));
        expect(w.get(1)).toBe(nextUp(0.1 + 0.2));
        expectContains(w, 0.3);
    });

    it('sub/mul/div leaves widen the rounded result by one ulp each way', () => {
        const s = SWInterval.sub(1.0, 0.3);
        expect(s.getEndpoints()).toEqual([nextDown(1.0 - 0.3), nextUp(1.0 - 0.3)]);
        expectContains(s, 0.7);

        const m = SWInterval.mul(0.1, 0.1);
        expect(m.getEndpoints()).toEqual([nextDown(0.1 * 0.1), nextUp(0.1 * 0.1)]);
        expectContains(m, 0.01);

        const d = SWInterval.div(1.0, 3.0);
        expect(d.getEndpoints()).toEqual([nextDown(1 / 3), nextUp(1 / 3)]);
        expectContains(d, 1 / 3);
    });

    it('div(u, 0) returns the entire set of real numbers', () => {
        const w = SWInterval.div(1.0, 0.0);
        expect(w.getEndpoints()).toEqual(
            [Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY]);
    });

    it('reciprocal helpers produce the upstream half-infinite intervals', () => {
        const rd = SWInterval.reciprocalDown(4);
        expect(rd.get(0)).toBe(nextDown(0.25));
        expect(rd.get(1)).toBe(Number.POSITIVE_INFINITY);
        const ru = SWInterval.reciprocalUp(-4);
        expect(ru.get(0)).toBe(Number.NEGATIVE_INFINITY);
        expect(ru.get(1)).toBe(nextUp(-0.25));
        expect(SWInterval.reals().getEndpoints()).toEqual(
            [Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY]);
    });
});

describe('SWInterval unary and additive operations', () => {
    it('negate maps [e0, e1] to [-e1, -e0] with no widening', () => {
        const v = new SWInterval(-1, 2);
        expect(v.negate().getEndpoints()).toEqual([-2, 1]);
    });

    it('add and sub produce the upstream endpoints', () => {
        const u = new SWInterval(2, 3);
        const v = new SWInterval(4, 5);
        expect(u.add(v).getEndpoints()).toEqual([nextDown(6), nextUp(8)]);
        expect(u.sub(v).getEndpoints()).toEqual([nextDown(-3), nextUp(-1)]);
        expect(u.add(10).getEndpoints()).toEqual([nextDown(12), nextUp(13)]);
        expect(u.sub(1).getEndpoints()).toEqual([nextDown(1), nextUp(2)]);
        expect(SWInterval.scalarSub(10, u).getEndpoints()).toEqual(
            [nextDown(7), nextUp(8)]);
    });
});

describe('SWInterval multiplication sign handling', () => {
    const cases: Array<[[number, number], [number, number], number, number]> = [
        // [u, v, exact lo, exact hi] for each sign configuration.
        [[2, 3], [4, 5], 8, 15],        // u >= 0, v >= 0
        [[2, 3], [-5, -4], -15, -8],    // u >= 0, v <= 0
        [[2, 3], [-4, 5], -12, 15],     // u >= 0, v mixed
        [[-3, -2], [4, 5], -15, -8],    // u <= 0, v >= 0
        [[-3, -2], [-5, -4], 8, 15],    // u <= 0, v <= 0
        [[-3, -2], [-4, 5], -15, 12],   // u <= 0, v mixed
        [[-2, 3], [4, 5], -10, 15],     // u mixed, v >= 0
        [[-2, 3], [-5, -4], -15, 10],   // u mixed, v <= 0
        [[-2, 3], [-5, 4], -15, 12]     // u mixed, v mixed (Mul2)
    ];

    it('produces the exact product bounds widened by one ulp', () => {
        for (const [ue, ve, lo, hi] of cases) {
            const u = new SWInterval(ue[0], ue[1]);
            const v = new SWInterval(ve[0], ve[1]);
            const w = u.mul(v);
            expect(w.get(0)).toBe(nextDown(lo));
            expect(w.get(1)).toBe(nextUp(hi));
        }
    });

    it('handles scalar multiplication of both signs', () => {
        const u = new SWInterval(-2, 3);
        expect(u.mul(2).getEndpoints()).toEqual([nextDown(-4), nextUp(6)]);
        expect(u.mul(-2).getEndpoints()).toEqual([nextDown(-6), nextUp(4)]);
        expect(u.mul(0).getEndpoints()).toEqual([nextDown(0), nextUp(0)]);
    });
});

describe('SWInterval division', () => {
    it('divides by an interval excluding zero', () => {
        const u = new SWInterval(1, 2);
        const v = new SWInterval(4, 5);
        const w = u.div(v);
        expectContains(w, 1 / 5);
        expectContains(w, 2 / 4);
        expect(w.get(0)).toBeGreaterThan(0.19);
        expect(w.get(1)).toBeLessThan(0.51);
    });

    it('divides by scalars of both signs and returns reals for zero', () => {
        const u = new SWInterval(1, 2);
        expect(u.div(2).getEndpoints()).toEqual([nextDown(0.5), nextUp(1)]);
        expect(u.div(-2).getEndpoints()).toEqual([nextDown(-1), nextUp(-0.5)]);
        expect(u.div(0).getEndpoints()).toEqual(
            [Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY]);
    });

    it('handles divisor intervals with a zero endpoint', () => {
        const u = new SWInterval(1, 2);
        // [1,2] / [0,4] = [1,2] * [nextDown(1/4), +inf). The upper bound
        // 2 * inf = inf is pulled back to MAX_VALUE by the nextafter
        // widening toward +max, exactly as in upstream Mul.
        const w0 = u.div(new SWInterval(0, 4));
        expect(w0.get(0)).toBeLessThanOrEqual(0.25);
        expect(w0.get(0)).toBeGreaterThan(0.24);
        expect(w0.get(1)).toBe(Number.MAX_VALUE);
        // [1,2] / [-4,0] = (~-MAX, ~-1/4], mirrored.
        const w1 = u.div(new SWInterval(-4, 0));
        expect(w1.get(0)).toBe(-Number.MAX_VALUE);
        expect(w1.get(1)).toBeGreaterThanOrEqual(-0.25);
        expect(w1.get(1)).toBeLessThan(-0.24);
    });

    it('returns reals when the divisor interval strictly contains zero', () => {
        const u = new SWInterval(1, 2);
        expect(u.div(new SWInterval(-1, 1)).getEndpoints()).toEqual(
            [Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY]);
        expect(SWInterval.scalarDiv(1, new SWInterval(-1, 1)).getEndpoints())
            .toEqual([Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY]);
    });

    it('scalarDiv matches the reciprocal formulation', () => {
        const w = SWInterval.scalarDiv(1, new SWInterval(2, 4));
        expectContains(w, 0.25);
        expectContains(w, 0.5);
        const wn = SWInterval.scalarDiv(-3, new SWInterval(2, 4));
        expectContains(wn, -1.5);
        expectContains(wn, -0.75);
    });
});

describe('SWInterval enclosure property (randomized)', () => {
    const rand = makeRandom(0xb13);

    function randomInterval(): SWInterval {
        const a = (rand() - 0.5) * 200;
        const b = (rand() - 0.5) * 200;
        return new SWInterval(Math.min(a, b), Math.max(a, b));
    }

    function sample(v: SWInterval): number {
        const t = rand();
        const s = v.get(0) + t * (v.get(1) - v.get(0));
        return Math.min(Math.max(s, v.get(0)), v.get(1));
    }

    it('add/sub/mul results always contain the pointwise double result', () => {
        for (let trial = 0; trial < 500; ++trial) {
            const u = randomInterval();
            const v = randomInterval();
            const su = sample(u);
            const sv = sample(v);
            expectContains(u.add(v), su + sv);
            expectContains(u.sub(v), su - sv);
            expectContains(u.mul(v), su * sv);
            expectContains(u.negate(), -su);
        }
    });

    it('div results contain the pointwise double result for nonzero divisors', () => {
        for (let trial = 0; trial < 500; ++trial) {
            const u = randomInterval();
            // Divisor with both endpoints of the same sign, bounded away
            // from zero.
            const sign = rand() < 0.5 ? -1 : 1;
            const a = sign * (0.5 + rand() * 100);
            const b = sign * (0.5 + rand() * 100);
            const v = new SWInterval(Math.min(a, b), Math.max(a, b));
            const su = sample(u);
            const sv = sample(v);
            expectContains(u.div(v), su / sv);
            expectContains(SWInterval.scalarDiv(su, v), su / sv);
        }
    });

    it('scalar operations always contain the pointwise double result', () => {
        for (let trial = 0; trial < 500; ++trial) {
            const u = randomInterval();
            const su = sample(u);
            const s = (rand() - 0.5) * 200;
            expectContains(u.add(s), su + s);
            expectContains(u.sub(s), su - s);
            expectContains(u.mul(s), su * s);
            expectContains(SWInterval.scalarSub(s, u), s - su);
            if (s !== 0) {
                expectContains(u.div(s), su / s);
            }
        }
    });

    it('composite expressions preserve enclosure', () => {
        for (let trial = 0; trial < 200; ++trial) {
            const a = (rand() - 0.5) * 20;
            const b = (rand() - 0.5) * 20;
            const c = (rand() - 0.5) * 20;
            const d = 0.5 + rand() * 10;
            // Leaf intervals for the raw floating-point variables, combined
            // at interior nodes: a*b + c/d - a.
            const w = SWInterval.mul(a, b)
                .add(SWInterval.div(c, d))
                .sub(new SWInterval(a));
            expectContains(w, a * b + c / d - a);
        }
    });

    it('degenerate zero-width intervals stay tight', () => {
        const v = new SWInterval(5);
        const diff = v.sub(v);
        expectContains(diff, 0);
        // The result is widened by exactly one ulp on each side of 0.
        expect(diff.get(0)).toBe(-Number.MIN_VALUE);
        expect(diff.get(1)).toBe(Number.MIN_VALUE);

        const z = new SWInterval(0);
        expect(z.mul(z).getEndpoints()).toEqual(
            [-Number.MIN_VALUE, Number.MIN_VALUE]);
        expectContains(z.add(z), 0);
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V05). The existing suite checks enclosure of the
// *rounded* double result; these properties check enclosure of the *exact*
// real result, computed with bigint rationals, which is what interval
// arithmetic actually promises.
// ---------------------------------------------------------------------------
describe('SWInterval verification', () => {
    // ---- exact rational arithmetic on doubles ------------------------------

    interface Rat { n: bigint; d: bigint; }   // d > 0

    /** The exact value of a finite double as a rational number. */
    function ratOf(v: number): Rat {
        const buffer = new DataView(new ArrayBuffer(8));
        buffer.setFloat64(0, v);
        const bits = buffer.getBigUint64(0);
        const s = (bits >> 63n) === 1n ? -1n : 1n;
        const e = Number((bits >> 52n) & 0x7FFn);
        const t = bits & 0xFFFFFFFFFFFFFn;
        const m = e === 0 ? t : (t | 0x10000000000000n);
        const p = (e === 0 ? -1074 : e - 1075);
        return p >= 0
            ? { n: s * m * (1n << BigInt(p)), d: 1n }
            : { n: s * m, d: 1n << BigInt(-p) };
    }

    const ratAdd = (a: Rat, b: Rat): Rat => ({ n: a.n * b.d + b.n * a.d, d: a.d * b.d });
    const ratSub = (a: Rat, b: Rat): Rat => ({ n: a.n * b.d - b.n * a.d, d: a.d * b.d });
    const ratMul = (a: Rat, b: Rat): Rat => ({ n: a.n * b.n, d: a.d * b.d });
    const ratDiv = (a: Rat, b: Rat): Rat => (b.n < 0n
        ? { n: -a.n * b.d, d: a.d * -b.n }
        : { n: a.n * b.d, d: a.d * b.n });
    /** sign(a - b) for two rationals with positive denominators. */
    const ratCmp = (a: Rat, b: Rat): number => {
        const l = a.n * b.d;
        const r = b.n * a.d;
        return l < r ? -1 : (l > r ? 1 : 0);
    };

    /** The exact real 'value' lies in the closed interval w. */
    function expectContainsExact(w: SWInterval, value: Rat): void {
        const lo = w.get(0);
        const hi = w.get(1);
        if (lo !== Number.NEGATIVE_INFINITY) {
            expect(ratCmp(ratOf(lo), value)).toBeLessThanOrEqual(0);
        }
        if (hi !== Number.POSITIVE_INFINITY) {
            expect(ratCmp(value, ratOf(hi))).toBeLessThanOrEqual(0);
        }
    }

    /** Endpoint equality that treats -0 and +0 as equal, as C++ == does. */
    function expectSameEndpoints(a: SWInterval, b: SWInterval): void {
        expect(a.get(0) === b.get(0)).toBe(true);
        expect(a.get(1) === b.get(1)).toBe(true);
    }

    // ---- arbitraries -------------------------------------------------------

    // Moderate magnitudes so that no bound overflows to an infinity, which
    // the nextafter widening would pull back to +/-MAX_VALUE and break
    // enclosure (a caveat inherited from upstream).
    const endpoint = fc.double({ min: -1e6, max: 1e6, noNaN: true, noDefaultInfinity: true });
    const awayFromZero = endpoint.filter(v => Math.abs(v) >= 1e-3);

    const intervalArb: fc.Arbitrary<SWInterval> = fc.tuple(endpoint, endpoint)
        .map(([a, b]) => new SWInterval(Math.min(a, b), Math.max(a, b)));

    // A divisor interval bounded away from zero, so the reciprocal cannot
    // overflow.
    const divisorArb: fc.Arbitrary<SWInterval> = fc.tuple(awayFromZero, awayFromZero)
        .filter(([a, b]) => (a > 0 && b > 0) || (a < 0 && b < 0))
        .map(([a, b]) => new SWInterval(Math.min(a, b), Math.max(a, b)));

    /** Points of the interval whose images bound the result of a monotone op. */
    function samplePoints(u: SWInterval): number[] {
        const a = u.get(0);
        const b = u.get(1);
        return [a, b, a + (b - a) * 0.5, a + (b - a) * 0.25, a + (b - a) * 0.75];
    }

    // ---- nextDown / nextUp -------------------------------------------------

    // Independent characterisation: nextUp(x) is the least double greater
    // than x, so nothing lies strictly between; the midpoint of x and its
    // successor must round back to one of the two.
    it('nextUp and nextDown are the adjacent doubles', () => {
        check(fc.double({ noNaN: true, noDefaultInfinity: true }), x => {
            const up = nextUp(x);
            const down = nextDown(x);
            if (x !== Number.MAX_VALUE) {
                expect(up).toBeGreaterThan(x);
                const mid = x / 2 + up / 2;
                expect(mid === x || mid === up).toBe(true);
                expect(nextDown(up) === x).toBe(true);
            }
            if (x !== -Number.MAX_VALUE) {
                expect(down).toBeLessThan(x);
                const mid = x / 2 + down / 2;
                expect(mid === x || mid === down).toBe(true);
                expect(nextUp(down) === x).toBe(true);
            }
            // Odd symmetry, which the interval operations rely on.
            expect(nextUp(-x) === -nextDown(x)).toBe(true);
            return true;
        });
    });

    // ---- enclosure of the exact real result --------------------------------

    it('add and sub enclose the exact real result for every pair of points', () => {
        check(fc.tuple(intervalArb, intervalArb), ([u, v]) => {
            const sum = u.add(v);
            const diff = u.sub(v);
            expect(sum.get(0)).toBeLessThanOrEqual(sum.get(1));
            expect(diff.get(0)).toBeLessThanOrEqual(diff.get(1));
            for (const a of samplePoints(u)) {
                for (const b of samplePoints(v)) {
                    expectContainsExact(sum, ratAdd(ratOf(a), ratOf(b)));
                    expectContainsExact(diff, ratSub(ratOf(a), ratOf(b)));
                }
            }
            return true;
        }, 100);
    });

    it('mul encloses the exact real product in every sign configuration', () => {
        check(fc.tuple(intervalArb, intervalArb), ([u, v]) => {
            const prod = u.mul(v);
            expect(prod.get(0)).toBeLessThanOrEqual(prod.get(1));
            for (const a of samplePoints(u)) {
                for (const b of samplePoints(v)) {
                    expectContainsExact(prod, ratMul(ratOf(a), ratOf(b)));
                }
            }
            return true;
        }, 100);
    });

    it('div encloses the exact real quotient when the divisor excludes zero', () => {
        check(fc.tuple(intervalArb, divisorArb), ([u, v]) => {
            const quot = u.div(v);
            expect(quot.get(0)).toBeLessThanOrEqual(quot.get(1));
            for (const a of samplePoints(u)) {
                for (const b of samplePoints(v)) {
                    expectContainsExact(quot, ratDiv(ratOf(a), ratOf(b)));
                }
            }
            return true;
        }, 100);
    });

    it('the scalar overloads enclose the exact real result', () => {
        check(fc.tuple(intervalArb, awayFromZero), ([u, s]) => {
            const rs = ratOf(s);
            for (const a of samplePoints(u)) {
                const ra = ratOf(a);
                expectContainsExact(u.add(s), ratAdd(ra, rs));
                expectContainsExact(u.sub(s), ratSub(ra, rs));
                expectContainsExact(u.mul(s), ratMul(ra, rs));
                expectContainsExact(SWInterval.scalarSub(s, u), ratSub(rs, ra));
                expectContainsExact(u.div(s), ratDiv(ra, rs));
            }
            return true;
        }, 100);
    });

    it('scalarDiv encloses the exact real quotient', () => {
        check(fc.tuple(awayFromZero, divisorArb), ([s, v]) => {
            const w = SWInterval.scalarDiv(s, v);
            for (const b of samplePoints(v)) {
                expectContainsExact(w, ratDiv(ratOf(s), ratOf(b)));
            }
            return true;
        }, 100);
    });

    it('the leaf-node operations enclose the exact real result', () => {
        check(fc.tuple(endpoint, endpoint), ([a, b]) => {
            const ra = ratOf(a), rb = ratOf(b);
            expectContainsExact(SWInterval.add(a, b), ratAdd(ra, rb));
            expectContainsExact(SWInterval.sub(a, b), ratSub(ra, rb));
            expectContainsExact(SWInterval.mul(a, b), ratMul(ra, rb));
            if (b === 0) {
                expect(SWInterval.div(a, b).getEndpoints())
                    .toEqual([-Infinity, Infinity]);
            } else if (Number.isFinite(a / b)) {
                expectContainsExact(SWInterval.div(a, b), ratDiv(ra, rb));
            }
            // else the quotient overflows; nextafter pulls the infinite bound
            // back to +/-MAX_VALUE and enclosure is lost. That is an upstream
            // caveat of the nextafter widening, exercised below.
            return true;
        });
    });

    // Upstream caveat, preserved: a bound that overflows to an infinity is
    // pulled back to +/-MAX_VALUE by nextafter(value, +/-max), so the result
    // no longer encloses the exact quotient.
    it('an overflowing leaf quotient loses enclosure at MAX_VALUE', () => {
        const w = SWInterval.div(8.881784197001252e-16, -Number.MIN_VALUE);
        expect(w.getEndpoints()).toEqual([-Number.MAX_VALUE, -Number.MAX_VALUE]);
    });

    // ---- algebraic identities of the widened operations --------------------

    it('interval addition and multiplication are commutative endpoint for endpoint', () => {
        check(fc.tuple(intervalArb, intervalArb), ([u, v]) => {
            expectSameEndpoints(u.mul(v), v.mul(u));
            expectSameEndpoints(u.add(v), v.add(u));
            return true;
        });
    });

    // nextDown(-x) = -nextUp(x), so negating one operand of a subtraction is
    // the same as negating the mirrored result. This pins the direction of
    // every widening in add/sub/negate at once.
    it('negation commutes with the additive operations', () => {
        check(fc.tuple(intervalArb, intervalArb), ([u, v]) => {
            expectSameEndpoints(u.negate().negate(), u);
            expectSameEndpoints(u.negate().add(v), v.sub(u));
            expectSameEndpoints(u.sub(v).negate(), v.sub(u));
            expectSameEndpoints(u.negate().mul(v), u.mul(v).negate());
            return true;
        });
    });

    it('a scalar operand behaves like the degenerate interval [s, s]', () => {
        check(fc.tuple(intervalArb, endpoint), ([u, s]) => {
            const deg = new SWInterval(s);
            expectSameEndpoints(u.add(s), u.add(deg));
            expectSameEndpoints(u.sub(s), u.sub(deg));
            expectSameEndpoints(u.mul(s), u.mul(deg));
            expectSameEndpoints(SWInterval.scalarSub(s, u), deg.sub(u));
            return true;
        });
    });

    it('a degenerate product is exactly one ulp wide on each side', () => {
        check(fc.tuple(endpoint, endpoint), ([a, b]) => {
            const w = new SWInterval(a).mul(new SWInterval(b));
            expect(w.get(0) === nextDown(a * b)).toBe(true);
            expect(w.get(1) === nextUp(a * b)).toBe(true);
            return true;
        });
    });

    // ---- divisor intervals containing zero ---------------------------------

    it('a divisor interval straddling zero yields the whole real line', () => {
        check(fc.tuple(intervalArb, fc.double({ min: -1e3, max: -1e-6, noNaN: true }),
            fc.double({ min: 1e-6, max: 1e3, noNaN: true })), ([u, lo, hi]) => {
                const v = new SWInterval(lo, hi);
                expect(u.div(v).getEndpoints()).toEqual([-Infinity, Infinity]);
                expect(SWInterval.scalarDiv(1, v).getEndpoints())
                    .toEqual([-Infinity, Infinity]);
                return true;
            });
    });

    // Upstream quirk, preserved: the half-infinite reciprocal of [0, b] is
    // multiplied by the dividend, and nextafter pulls the resulting infinite
    // bound back to +/-MAX_VALUE instead of leaving it infinite.
    it('a divisor interval with a zero endpoint yields a bound at MAX_VALUE', () => {
        check(fc.tuple(fc.double({ min: 1, max: 100, noNaN: true }),
            fc.double({ min: 1, max: 100, noNaN: true })), ([a, b]) => {
                const u = new SWInterval(a, a + b);
                const w = u.div(new SWInterval(0, b));
                expect(w.get(1)).toBe(Number.MAX_VALUE);
                expectContainsExact(w, ratDiv(ratOf(a), ratOf(b)));
                const w2 = u.div(new SWInterval(-b, 0));
                expect(w2.get(0)).toBe(-Number.MAX_VALUE);
                expectContainsExact(w2, ratDiv(ratOf(a), ratOf(-b)));
                return true;
            });
    });

    // ---- value semantics ---------------------------------------------------

    it('results and getEndpoints never alias the operands', () => {
        const u = new SWInterval(1, 2);
        const v = new SWInterval(3, 4);
        const e = u.getEndpoints();
        e[0] = 99;
        expect(u.get(0)).toBe(1);
        for (const w of [u.clone(), u.negate(), u.add(v), u.sub(v), u.mul(v),
            u.div(v), u.add(1), SWInterval.scalarSub(1, u), SWInterval.scalarDiv(1, u)]) {
            expect(w).not.toBe(u);
            expect(w).not.toBe(v);
        }
        expect(u.getEndpoints()).not.toBe(u.getEndpoints());
    });
});
