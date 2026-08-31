import { describe, it, expect } from 'vitest';
import { FPInterval } from '../src/FPInterval';

// Deterministic pseudorandom generator (LCG) for reproducible tests.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// ---------------------------------------------------------------------------
// Exact (arbitrary precision) reference arithmetic used to verify enclosure.
// Every finite double is exactly m * 2^e for integers m and e, so sums and
// products of doubles are compared without any rounding by using BigInt.
// ---------------------------------------------------------------------------

interface Dyadic { m: bigint; e: number; }

const dyadicBuffer = new ArrayBuffer(8);
const dyadicF = new Float64Array(dyadicBuffer);
const dyadicU = new BigUint64Array(dyadicBuffer);

function toDyadic(x: number): Dyadic {
    dyadicF[0] = x;
    const bits = dyadicU[0];
    const negative = ((bits >> 63n) & 1n) === 1n;
    const biased = Number((bits >> 52n) & 0x7ffn);
    const frac = bits & 0xfffffffffffffn;
    let m: bigint;
    let e: number;
    if (biased === 0) {
        m = frac;
        e = -1074;
    } else {
        m = frac | (1n << 52n);
        e = biased - 1075;
    }
    return { m: negative ? -m : m, e };
}

// Scale a dyadic to the exponent 'e' (which must be <= d.e).
function scaleTo(d: Dyadic, e: number): bigint {
    return d.m << BigInt(d.e - e);
}

function addDyadic(a: Dyadic, b: Dyadic): Dyadic {
    const e = Math.min(a.e, b.e);
    return { m: scaleTo(a, e) + scaleTo(b, e), e };
}

function mulDyadic(a: Dyadic, b: Dyadic): Dyadic {
    return { m: a.m * b.m, e: a.e + b.e };
}

// Returns -1, 0 or +1 as a < b, a == b or a > b.
function cmpDyadic(a: Dyadic, b: Dyadic): number {
    const e = Math.min(a.e, b.e);
    const ma = scaleTo(a, e);
    const mb = scaleTo(b, e);
    return ma < mb ? -1 : (ma > mb ? 1 : 0);
}

// The interval [w0, w1] must contain the exact value 'exact'. Infinite
// endpoints are accepted as bounds on either side.
function expectEnclosesDyadic(w: FPInterval, exact: Dyadic): void {
    const lo = w.get(0);
    const hi = w.get(1);
    if (lo !== Number.NEGATIVE_INFINITY) {
        expect(Number.isFinite(lo)).toBe(true);
        expect(cmpDyadic(toDyadic(lo), exact)).toBeLessThanOrEqual(0);
    }
    if (hi !== Number.POSITIVE_INFINITY) {
        expect(Number.isFinite(hi)).toBe(true);
        expect(cmpDyadic(toDyadic(hi), exact)).toBeGreaterThanOrEqual(0);
    }
}

function expectContains(w: FPInterval, value: number): void {
    expect(w.get(0)).toBeLessThanOrEqual(value);
    expect(w.get(1)).toBeGreaterThanOrEqual(value);
}

// A random double in [-1024, 1024] with a full 53-bit significand.
function randomValue(rand: () => number): number {
    return (2 * rand() - 1) * 1024 * (1 + rand());
}

describe('FPInterval construction and access', () => {
    it('default-constructs the degenerate interval [0, 0]', () => {
        const u = new FPInterval();
        expect(u.get(0)).toBe(0);
        expect(u.get(1)).toBe(0);
    });

    it('constructs the degenerate interval [e, e] from one endpoint', () => {
        const u = new FPInterval(1.25);
        expect(u.getEndpoints()).toEqual([1.25, 1.25]);
    });

    it('constructs [e0, e1] and clones by value', () => {
        const u = new FPInterval(-2, 3);
        const v = u.clone();
        expect(v.getEndpoints()).toEqual([-2, 3]);
        expect(v).not.toBe(u);
    });

    it('constructs from an endpoint pair', () => {
        const u = FPInterval.fromEndpoints([-1, 4]);
        expect(u.getEndpoints()).toEqual([-1, 4]);
    });

    it('optionally traps e0 > e1 (GTE_THROW_ON_INVALID_INTERVAL)', () => {
        expect(() => new FPInterval(3, 2)).not.toThrow();
        FPInterval.throwOnInvalid = true;
        try {
            expect(() => new FPInterval(3, 2)).toThrow(/Invalid FPInterval/);
            expect(() => new FPInterval(2, 3)).not.toThrow();
            // The one-argument form is always valid.
            expect(() => new FPInterval(3)).not.toThrow();
        } finally {
            FPInterval.throwOnInvalid = false;
        }
    });
});

describe('FPInterval emulated directed rounding', () => {
    it('rounds one representable value outward', () => {
        expect(FPInterval.roundUp(1)).toBe(1 + Math.pow(2, -52));
        expect(FPInterval.roundDown(1)).toBe(1 - Math.pow(2, -53));
        expect(FPInterval.roundUp(-1)).toBe(-(1 - Math.pow(2, -53)));
        expect(FPInterval.roundDown(-1)).toBe(-(1 + Math.pow(2, -52)));
    });

    it('crosses zero through the subnormals', () => {
        expect(FPInterval.roundUp(0)).toBe(Number.MIN_VALUE);
        expect(FPInterval.roundDown(0)).toBe(-Number.MIN_VALUE);
    });

    // This is the property that distinguishes FPInterval from the sibling
    // SWInterval port: directed rounding must leave an infinite bound alone
    // (upstream issue: std::nextafter on an infinite bound breaks enclosure).
    it('leaves an infinite bound infinite in the outward direction', () => {
        expect(FPInterval.roundUp(Number.POSITIVE_INFINITY))
            .toBe(Number.POSITIVE_INFINITY);
        expect(FPInterval.roundDown(Number.NEGATIVE_INFINITY))
            .toBe(Number.NEGATIVE_INFINITY);
    });

    it('pulls an infinity back toward the finite range in the inward direction', () => {
        // An overflowed upper computation still bounds the exact value from
        // below by MAX_VALUE.
        expect(FPInterval.roundDown(Number.POSITIVE_INFINITY))
            .toBe(Number.MAX_VALUE);
        expect(FPInterval.roundUp(Number.NEGATIVE_INFINITY))
            .toBe(-Number.MAX_VALUE);
    });

    it('propagates NaN', () => {
        expect(Number.isNaN(FPInterval.roundUp(NaN))).toBe(true);
        expect(Number.isNaN(FPInterval.roundDown(NaN))).toBe(true);
    });
});

describe('FPInterval leaf-node operations', () => {
    it('produces degenerate intervals for exact sums and differences', () => {
        expect(FPInterval.add(0.5, 0.25).getEndpoints()).toEqual([0.75, 0.75]);
        expect(FPInterval.sub(1, 0.25).getEndpoints()).toEqual([0.75, 0.75]);
        expect(FPInterval.add(1, -1).getEndpoints()).toEqual([0, 0]);
        expect(FPInterval.add(3, 4).getEndpoints()).toEqual([7, 7]);
    });

    it('produces degenerate intervals for exact products and quotients', () => {
        expect(FPInterval.mul(1.5, 4).getEndpoints()).toEqual([6, 6]);
        expect(FPInterval.mul(3, 0).getEndpoints()).toEqual([0, 0]);
        expect(FPInterval.div(1, 4).getEndpoints()).toEqual([0.25, 0.25]);
        expect(FPInterval.div(7, 2).getEndpoints()).toEqual([3.5, 3.5]);
    });

    it('widens by one ulp for inexact results', () => {
        const w = FPInterval.add(1, Math.pow(2, -60));
        expect(w.get(0)).toBe(FPInterval.roundDown(1));
        expect(w.get(1)).toBe(FPInterval.roundUp(1));
        expect(w.get(0)).toBeLessThan(w.get(1));

        const q = FPInterval.div(1, 3);
        expect(q.get(0)).toBeLessThan(1 / 3);
        expect(q.get(1)).toBeGreaterThan(1 / 3);
    });

    it('returns the reals for division by zero', () => {
        expect(FPInterval.div(1, 0).getEndpoints())
            .toEqual([Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY]);
        expect(FPInterval.reals().getEndpoints())
            .toEqual([Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY]);
    });

    it('encloses the exact sum, difference and product (randomized)', () => {
        const rand = makeRandom(20260830);
        for (let trial = 0; trial < 4000; ++trial) {
            const a = randomValue(rand);
            const b = randomValue(rand);
            const da = toDyadic(a);
            const db = toDyadic(b);

            expectEnclosesDyadic(FPInterval.add(a, b), addDyadic(da, db));
            expectEnclosesDyadic(FPInterval.sub(a, b),
                addDyadic(da, { m: -db.m, e: db.e }));
            expectEnclosesDyadic(FPInterval.mul(a, b), mulDyadic(da, db));
        }
    });

    it('is degenerate exactly when the operation is exact (randomized)', () => {
        const rand = makeRandom(777);
        for (let trial = 0; trial < 4000; ++trial) {
            const a = randomValue(rand);
            const b = randomValue(rand);

            const sum = FPInterval.add(a, b);
            const exactSum =
                cmpDyadic(toDyadic(a + b), addDyadic(toDyadic(a), toDyadic(b))) === 0;
            expect(sum.get(0) === sum.get(1)).toBe(exactSum);
            if (exactSum) {
                expect(sum.get(0)).toBe(a + b);
            }

            const product = FPInterval.mul(a, b);
            const exactProduct =
                cmpDyadic(toDyadic(a * b), mulDyadic(toDyadic(a), toDyadic(b))) === 0;
            // The Dekker guards can report "not proven exact" for a product
            // that is in fact exact, which only widens the interval; they
            // never report exactness for an inexact product.
            if (product.get(0) === product.get(1)) {
                expect(exactProduct).toBe(true);
                expect(product.get(0)).toBe(a * b);
            }
        }
    });

    it('encloses the exact quotient (randomized)', () => {
        const rand = makeRandom(31337);
        for (let trial = 0; trial < 4000; ++trial) {
            const a = randomValue(rand);
            const b = randomValue(rand);
            if (b === 0) {
                continue;
            }
            const w = FPInterval.div(a, b);
            // exact(a/b) is between lo and hi iff lo*b and hi*b bracket a
            // (with the sign of b taken into account). Use exact products.
            const da = toDyadic(a);
            const dlo = mulDyadic(toDyadic(w.get(0)), toDyadic(b));
            const dhi = mulDyadic(toDyadic(w.get(1)), toDyadic(b));
            if (b > 0) {
                expect(cmpDyadic(dlo, da)).toBeLessThanOrEqual(0);
                expect(cmpDyadic(dhi, da)).toBeGreaterThanOrEqual(0);
            } else {
                expect(cmpDyadic(dlo, da)).toBeGreaterThanOrEqual(0);
                expect(cmpDyadic(dhi, da)).toBeLessThanOrEqual(0);
            }
        }
    });
});

describe('FPInterval arithmetic on intervals', () => {
    it('negates by reversing the endpoints', () => {
        expect(new FPInterval(-2, 3).negate().getEndpoints()).toEqual([-3, 2]);
    });

    it('adds and subtracts intervals with exact endpoints', () => {
        const u = new FPInterval(1, 2);
        const v = new FPInterval(0.25, 0.5);
        expect(u.add(v).getEndpoints()).toEqual([1.25, 2.5]);
        expect(u.sub(v).getEndpoints()).toEqual([0.5, 1.75]);
        expect(u.add(1).getEndpoints()).toEqual([2, 3]);
        expect(u.sub(1).getEndpoints()).toEqual([0, 1]);
        expect(FPInterval.scalarSub(4, u).getEndpoints()).toEqual([2, 3]);
    });

    it('multiplies intervals over all nine sign cases', () => {
        const cases: Array<[number, number]> = [[1, 2], [-2, -1], [-1, 3]];
        for (const [u0, u1] of cases) {
            for (const [v0, v1] of cases) {
                const w = new FPInterval(u0, u1).mul(new FPInterval(v0, v1));
                const products = [u0 * v0, u0 * v1, u1 * v0, u1 * v1];
                expect(w.get(0)).toBeLessThanOrEqual(Math.min(...products));
                expect(w.get(1)).toBeGreaterThanOrEqual(Math.max(...products));
                // The bounds are attained (all endpoints here are exact).
                expect(w.get(0)).toBe(Math.min(...products));
                expect(w.get(1)).toBe(Math.max(...products));
            }
        }
    });

    it('multiplies by a scalar on either side of zero', () => {
        const u = new FPInterval(-1, 3);
        expect(u.mul(2).getEndpoints()).toEqual([-2, 6]);
        expect(u.mul(-2).getEndpoints()).toEqual([-6, 2]);
        // The product with the scalar zero is the degenerate interval at
        // zero (the lower endpoint is the signed zero -0, as in C++).
        const zeroProduct = u.mul(0);
        expect(zeroProduct.get(0) === 0).toBe(true);
        expect(zeroProduct.get(1) === 0).toBe(true);
    });

    it('encloses every product of sampled operands (randomized)', () => {
        const rand = makeRandom(99);
        for (let trial = 0; trial < 500; ++trial) {
            const a = randomValue(rand);
            const b = randomValue(rand);
            const c = randomValue(rand);
            const d = randomValue(rand);
            const u = new FPInterval(Math.min(a, b), Math.max(a, b));
            const v = new FPInterval(Math.min(c, d), Math.max(c, d));
            const w = u.mul(v);
            for (let s = 0; s <= 4; ++s) {
                for (let t = 0; t <= 4; ++t) {
                    // The interpolation rounds, so clamp back into the
                    // interval before using the point as a witness.
                    const x = Math.min(u.get(1), Math.max(u.get(0),
                        u.get(0) + (s / 4) * (u.get(1) - u.get(0))));
                    const y = Math.min(v.get(1), Math.max(v.get(0),
                        v.get(0) + (t / 4) * (v.get(1) - v.get(0))));
                    expectEnclosesDyadic(w, mulDyadic(toDyadic(x), toDyadic(y)));
                }
            }
            // productLowerBound/productUpperBound agree with mul whenever
            // at least one of the two intervals does not straddle zero.
            // (Upstream's both-straddle branch is wrong; see the dedicated
            // test below.)
            const uStraddles = u.get(0) < 0 && u.get(1) > 0;
            const vStraddles = v.get(0) < 0 && v.get(1) > 0;
            if (!uStraddles || !vStraddles) {
                const lo = FPInterval.productLowerBound(
                    u.getEndpoints(), v.getEndpoints());
                const hi = FPInterval.productUpperBound(
                    u.getEndpoints(), v.getEndpoints());
                expect(lo).toBe(w.get(0));
                expect(hi).toBe(w.get(1));
            }
        }
    });

    it('preserves the upstream productLowerBound/productUpperBound defect', () => {
        // Upstream's both-straddle branch returns u[0]*v[0] as the lower
        // bound and u[1]*v[1] as the upper bound. Neither is correct; the
        // port preserves the behavior and this test documents it. See the
        // comment on productLowerBound.
        const u = [-4, 1];
        const v = [-4, 1];
        expect(FPInterval.productLowerBound(u, v)).toBe(16);
        expect(FPInterval.productUpperBound(u, v)).toBe(1);

        // The correct product interval, from mul, is [-4, 16].
        const w = new FPInterval(u[0], u[1]).mul(new FPInterval(v[0], v[1]));
        expect(w.getEndpoints()).toEqual([-4, 16]);
    });

    it('divides by an interval that excludes zero', () => {
        const u = new FPInterval(1, 2);
        const v = new FPInterval(4, 8);
        const w = u.div(v);
        expectContains(w, 1 / 8);
        expectContains(w, 2 / 4);
        expect(w.get(0)).toBeLessThanOrEqual(0.125);
        expect(w.get(1)).toBeGreaterThanOrEqual(0.5);

        const negative = u.div(new FPInterval(-8, -4));
        expectContains(negative, -0.5);
        expectContains(negative, -0.125);
    });

    it('divides by a scalar', () => {
        const u = new FPInterval(1, 2);
        expect(u.div(4).getEndpoints()).toEqual([0.25, 0.5]);
        expect(u.div(-4).getEndpoints()).toEqual([-0.5, -0.25]);
        expect(u.div(0).getEndpoints())
            .toEqual([Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY]);
    });

    it('returns the reals when the divisor straddles zero', () => {
        const u = new FPInterval(1, 2);
        expect(u.div(new FPInterval(-1, 1)).getEndpoints())
            .toEqual([Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY]);
        expect(FPInterval.scalarDiv(1, new FPInterval(-1, 1)).getEndpoints())
            .toEqual([Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY]);
    });

    it('handles a divisor with a zero endpoint via the half-line reciprocals', () => {
        // [0, 2] has reciprocal [1/2, +inf).
        const w = new FPInterval(1, 1).div(new FPInterval(0, 2));
        expect(w.get(1)).toBe(Number.POSITIVE_INFINITY);
        expect(w.get(0)).toBeLessThanOrEqual(0.5);

        // [-2, 0] has reciprocal (-inf, -1/2].
        const v = new FPInterval(1, 1).div(new FPInterval(-2, 0));
        expect(v.get(0)).toBe(Number.NEGATIVE_INFINITY);
        expect(v.get(1)).toBeGreaterThanOrEqual(-0.5);
    });

    it('agrees with the free-function scalar division', () => {
        const v = new FPInterval(4, 8);
        expect(FPInterval.scalarDiv(2, v).getEndpoints())
            .toEqual(new FPInterval(2, 2).div(v).getEndpoints());
    });

    it('keeps infinite endpoints infinite through arithmetic', () => {
        // The regression guard for the nextafter-on-infinity defect: an
        // interval with a +inf upper bound must keep it after adding or
        // multiplying by a positive quantity.
        const halfLine = FPInterval.reciprocalDown(2);
        expect(halfLine.get(0)).toBeLessThanOrEqual(0.5);
        expect(halfLine.get(1)).toBe(Number.POSITIVE_INFINITY);

        expect(halfLine.add(1).get(1)).toBe(Number.POSITIVE_INFINITY);
        expect(halfLine.mul(3).get(1)).toBe(Number.POSITIVE_INFINITY);
        expect(halfLine.sub(1).get(1)).toBe(Number.POSITIVE_INFINITY);

        const lowerHalfLine = FPInterval.reciprocalUp(-2);
        expect(lowerHalfLine.get(0)).toBe(Number.NEGATIVE_INFINITY);
        expect(lowerHalfLine.get(1)).toBeGreaterThanOrEqual(-0.5);
        expect(lowerHalfLine.add(1).get(0)).toBe(Number.NEGATIVE_INFINITY);
        expect(lowerHalfLine.mul(3).get(0)).toBe(Number.NEGATIVE_INFINITY);
    });

    it('overflows outward rather than reporting an unbounded interval bound', () => {
        // MAX_VALUE + MAX_VALUE overflows in round-to-nearest. The lower
        // bound must stay finite (it still bounds the exact sum from below)
        // and the upper bound must become +inf.
        const w = FPInterval.add(Number.MAX_VALUE, Number.MAX_VALUE);
        expect(w.get(0)).toBe(Number.MAX_VALUE);
        expect(w.get(1)).toBe(Number.POSITIVE_INFINITY);

        const v = FPInterval.add(-Number.MAX_VALUE, -Number.MAX_VALUE);
        expect(v.get(0)).toBe(Number.NEGATIVE_INFINITY);
        expect(v.get(1)).toBe(-Number.MAX_VALUE);
    });

    it('encloses an evaluated expression tree', () => {
        // f(x, y) = (x + y) * (x - y) / (x * y) evaluated at nearby values.
        const rand = makeRandom(4242);
        for (let trial = 0; trial < 200; ++trial) {
            const x = 1 + rand();
            const y = 2 + rand();
            const ix = new FPInterval(x, x);
            const iy = new FPInterval(y, y);
            const w = ix.add(iy).mul(ix.sub(iy)).div(ix.mul(iy));
            const expected = ((x + y) * (x - y)) / (x * y);
            expectContains(w, expected);
            // The interval is tight: a handful of ulps wide at most.
            expect(w.get(1) - w.get(0)).toBeLessThan(1e-12);
        }
    });
});
