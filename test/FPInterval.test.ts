import { describe, it, expect } from 'vitest';
import { FPInterval } from '../src/FPInterval.js';
import { BSRational } from '../src/BSRational.js';
import { check, fc, finite } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (group V38). The defining property of an interval type is
// enclosure: whatever the emulated directed rounding does, the exact real
// result of every operation must lie between the two endpoints. These
// properties check that against BSRational, which represents every dyadic
// rational (and every quotient of two of them) exactly.
// ---------------------------------------------------------------------------
describe('FPInterval verification', () => {
    const rat = (x: number): BSRational => BSRational.fromNumber(x);

    // The exact value must satisfy lo <= exact <= hi. Infinite endpoints are
    // vacuously satisfied (BSRational cannot represent them).
    function expectEncloses(w: FPInterval, exact: BSRational,
        label: string): void {
        const lo = w.get(0);
        const hi = w.get(1);
        expect(lo, `${label}: lo=${lo} must not be NaN`).not.toBeNaN();
        expect(hi, `${label}: hi=${hi} must not be NaN`).not.toBeNaN();
        expect(lo <= hi, `${label}: [${lo}, ${hi}] is not an interval`)
            .toBe(true);
        if (Number.isFinite(lo)) {
            expect(rat(lo).lessThanOrEqual(exact),
                `${label}: lo=${lo} exceeds the exact value`).toBe(true);
        } else {
            expect(lo).toBe(Number.NEGATIVE_INFINITY);
        }
        if (Number.isFinite(hi)) {
            expect(rat(hi).greaterThanOrEqual(exact),
                `${label}: hi=${hi} is below the exact value`).toBe(true);
        } else {
            expect(hi).toBe(Number.POSITIVE_INFINITY);
        }
    }

    // Operands spanning the whole exponent range, subnormals included: the
    // enclosure is an exact statement, so there is no tolerance to protect.
    const value = fc.oneof(
        finite(-10, 10),
        finite(-1, 1),
        fc.constantFrom(0, -0, 1, -1, 0.5, -0.5, 3, Number.MIN_VALUE,
            -Number.MIN_VALUE, Number.MAX_VALUE, -Number.MAX_VALUE,
            2.2250738585072014e-308, 1e300, -1e300, 1e-300),
        fc.integer({ min: -1000, max: 1000 }).map(i => i / 8)
    );

    const interval = fc.tuple(value, value).map(([a, b]) =>
        new FPInterval(Math.min(a, b), Math.max(a, b)));

    it('roundDown and roundUp step one representable value outward', () => {
        check(value, x => {
            const down = FPInterval.roundDown(x);
            const up = FPInterval.roundUp(x);
            expect(down).toBeLessThanOrEqual(x);
            expect(up).toBeGreaterThanOrEqual(x);
            if (Number.isFinite(x)) {
                // Exactly one step, so no representable value is skipped.
                // '+ 0' normalizes the -0/+0 tie that toBe (Object.is) would
                // otherwise reject at the zero crossing.
                expect(FPInterval.roundUp(down) + 0).toBe(x + 0);
                expect(FPInterval.roundDown(up) + 0).toBe(x + 0);
            }
        });
    });

    it('the infinities round outward the way enclosure requires', () => {
        // An upper computation that overflowed to +inf must stay +inf, and a
        // lower one must fall back to a finite bound below the true value.
        expect(FPInterval.roundUp(Number.POSITIVE_INFINITY))
            .toBe(Number.POSITIVE_INFINITY);
        expect(FPInterval.roundDown(Number.POSITIVE_INFINITY))
            .toBe(Number.MAX_VALUE);
        expect(FPInterval.roundDown(Number.NEGATIVE_INFINITY))
            .toBe(Number.NEGATIVE_INFINITY);
        expect(FPInterval.roundUp(Number.NEGATIVE_INFINITY))
            .toBe(-Number.MAX_VALUE);
    });

    it('the leaf operations enclose the exact sum, difference and product',
        () => {
            check(fc.tuple(value, value), ([u, v]) => {
                expectEncloses(FPInterval.add(u, v), rat(u).add(rat(v)),
                    `add(${u},${v})`);
                expectEncloses(FPInterval.sub(u, v), rat(u).sub(rat(v)),
                    `sub(${u},${v})`);
                expectEncloses(FPInterval.mul(u, v), rat(u).mul(rat(v)),
                    `mul(${u},${v})`);
                if (v !== 0) {
                    expectEncloses(FPInterval.div(u, v), rat(u).div(rat(v)),
                        `div(${u},${v})`);
                } else {
                    expect(FPInterval.div(u, v).getEndpoints()).toEqual(
                        [Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY]);
                }
            });
        });

    it('interval addition and subtraction enclose the exact endpoint results',
        () => {
            check(fc.tuple(interval, interval), ([u, v]) => {
                const sum = u.add(v);
                expectEncloses(sum, rat(u.get(0)).add(rat(v.get(0))),
                    'add lower');
                expectEncloses(sum, rat(u.get(1)).add(rat(v.get(1))),
                    'add upper');
                const diff = u.sub(v);
                expectEncloses(diff, rat(u.get(0)).sub(rat(v.get(1))),
                    'sub lower');
                expectEncloses(diff, rat(u.get(1)).sub(rat(v.get(0))),
                    'sub upper');
            });
        });

    it('interval multiplication encloses every exact corner product', () => {
        // A product over a box attains its extremes at the corners, so the
        // interval must contain all four exact corner products.
        check(fc.tuple(interval, interval), ([u, v]) => {
            const w = u.mul(v);
            for (const a of u.getEndpoints()) {
                for (const b of v.getEndpoints()) {
                    expectEncloses(w, rat(a).mul(rat(b)),
                        `mul corner ${a}*${b}`);
                }
            }
        });
    });

    it('scalar multiplication agrees with multiplication by a point interval',
        () => {
            check(fc.tuple(interval, value), ([u, s]) => {
                expect(u.mul(s).getEndpoints())
                    .toEqual(u.mul(new FPInterval(s)).getEndpoints());
                for (const a of u.getEndpoints()) {
                    expectEncloses(u.mul(s), rat(a).mul(rat(s)),
                        `scalar mul ${a}*${s}`);
                }
            });
        });

    it('interval division encloses every exact corner quotient when the divisor '
        + 'excludes zero', () => {
            check(fc.tuple(interval, interval), ([u, v]) => {
                if (!(v.get(0) > 0 || v.get(1) < 0)) {
                    return;
                }
                // A subnormal divisor endpoint makes its reciprocal overflow
                // to an infinity, and an infinity times a zero endpoint of the
                // dividend is NaN. That is upstream behaviour (directed
                // rounding overflows to the infinity as well); it is pinned by
                // its own test below rather than swept into this property.
                if (!Number.isFinite(1 / v.get(0))
                    || !Number.isFinite(1 / v.get(1))) {
                    return;
                }
                const w = u.div(v);
                for (const a of u.getEndpoints()) {
                    for (const b of v.getEndpoints()) {
                        expectEncloses(w, rat(a).div(rat(b)),
                            `div corner ${a}/${b}`);
                    }
                }
            });
        });

    it('a divisor that straddles zero produces the reals, and a zero endpoint '
        + 'produces a half line', () => {
            check(fc.tuple(interval, interval), ([u, v]) => {
                const v0 = v.get(0);
                const v1 = v.get(1);
                if (v0 < 0 && v1 > 0) {
                    expect(u.div(v).getEndpoints()).toEqual(
                        [Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY]);
                    expect(FPInterval.scalarDiv(1, v).getEndpoints()).toEqual(
                        [Number.NEGATIVE_INFINITY, Number.POSITIVE_INFINITY]);
                } else if ((v0 === 0 || v1 === 0) && v0 !== v1
                    && u.get(0) !== 0 && u.get(1) !== 0) {
                    // Upstream replaces 1/[0,v1] by [1/v1, +inf) and
                    // 1/[v0,0] by (-inf, 1/v0], which is a superset of the
                    // true reciprocal set, then multiplies. The half line is
                    // unbounded on the side the zero endpoint produces.
                    const w = u.div(v);
                    expect(w.get(0)).toBeLessThanOrEqual(w.get(1));
                    expect(Number.isFinite(w.get(0))
                        && Number.isFinite(w.get(1))).toBe(false);
                }
            });
        });

    it('upstream (preserved): an infinite reciprocal times a zero endpoint is '
        + 'NaN, so the enclosure is lost', () => {
            // 1/Number.MIN_VALUE overflows to +inf. Upstream's directed
            // rounding overflows to the infinity in exactly the same way, and
            // operator* then evaluates 0 * (+inf).
            const w = new FPInterval(0, 0).div(
                new FPInterval(Number.MIN_VALUE, 1));
            expect(Number.isNaN(w.get(0)) || Number.isNaN(w.get(1))).toBe(true);
            // The half-line reciprocal of a divisor with a zero endpoint does
            // the same when the dividend is the degenerate zero interval.
            const h = new FPInterval(0, 0).div(new FPInterval(0, 2));
            expect(Number.isNaN(h.get(0)) || Number.isNaN(h.get(1))).toBe(true);
            // With no zero endpoint in the dividend the result is the
            // (unbounded but valid) half line.
            expect(new FPInterval(1, 2).div(new FPInterval(0, 2))
                .getEndpoints()).toEqual([0.5, Number.POSITIVE_INFINITY]);
            // The same loss happens without any division: an intermediate
            // that overflows to an infinity, multiplied by a zero endpoint.
            const overflowed = new FPInterval(Number.MAX_VALUE)
                .add(new FPInterval(Number.MAX_VALUE));
            expect(overflowed.get(1)).toBe(Number.POSITIVE_INFINITY);
            const lost = overflowed.mul(new FPInterval(0));
            expect(Number.isNaN(lost.get(0)) || Number.isNaN(lost.get(1)))
                .toBe(true);
        });

    it('scalarSub and scalarDiv agree with the point-interval operations', () => {
        check(fc.tuple(value, interval), ([s, v]) => {
            expect(FPInterval.scalarSub(s, v).getEndpoints())
                .toEqual(new FPInterval(s).sub(v).getEndpoints());
            expect(FPInterval.scalarDiv(s, v).getEndpoints())
                .toEqual(new FPInterval(s).div(v).getEndpoints());
        });
    });

    it('negation reverses the endpoints and commutes with the operations', () => {
        check(fc.tuple(interval, interval), ([u, v]) => {
            // '+ 0' on both sides normalizes the -0/+0 ties that toEqual
            // (Object.is) would otherwise reject.
            const norm = (w: FPInterval) => w.getEndpoints().map(x => x + 0);
            expect(norm(u.negate())).toEqual([-u.get(1) + 0, -u.get(0) + 0]);
            // -(-u) = u.
            expect(norm(u.negate().negate()))
                .toEqual([u.get(0) + 0, u.get(1) + 0]);
            // (-u)*v = -(u*v) as sets: the endpoints must match after
            // reversal, since every branch of operator* uses the same
            // rounding on the mirrored corners.
            const a = u.negate().mul(v).getEndpoints();
            const b = u.mul(v).negate().getEndpoints();
            expect(a[0] + 0).toBe(b[0] + 0);
            expect(a[1] + 0).toBe(b[1] + 0);
        });
    });

    it('a degenerate interval stays degenerate exactly when the result is '
        + 'representable', () => {
            check(fc.tuple(value, value), ([u, v]) => {
                const sum = FPInterval.add(u, v);
                if (sum.get(0) === sum.get(1)) {
                    // Proven exact, so the exact rational sum is the double.
                    expect(rat(sum.get(0)).equals(rat(u).add(rat(v)))).toBe(true);
                }
                const product = FPInterval.mul(u, v);
                if (product.get(0) === product.get(1)) {
                    expect(rat(product.get(0)).equals(rat(u).mul(rat(v))))
                        .toBe(true);
                }
            });
        });

    it('an expression tree evaluated on intervals encloses its exact value',
        () => {
            // w = (a + b) * c - a / d, with the same expression evaluated
            // exactly in BSRational.
            // Moderate operands only: an intermediate that overflows to an
            // infinity and is then multiplied by a zero endpoint gives NaN,
            // which is upstream behaviour and is pinned above. 1/d must not
            // overflow either, for the same reason.
            const moderate = fc.oneof(finite(-10, 10),
                fc.integer({ min: -1000, max: 1000 }).map(i => i / 8),
                fc.constantFrom(0, -0, 1, -1, 0.5, 3));
            check(fc.tuple(moderate, moderate, moderate,
                moderate.filter(x => x !== 0 && Number.isFinite(1 / x))),
                ([a, b, c, d]) => {
                    const w = new FPInterval(a).add(new FPInterval(b))
                        .mul(new FPInterval(c))
                        .sub(new FPInterval(a).div(new FPInterval(d)));
                    const exact = rat(a).add(rat(b)).mul(rat(c))
                        .sub(rat(a).div(rat(d)));
                    expectEncloses(w, exact, 'expression tree');
                }, 100);
        });

    it('upstream defect (preserved): productLowerBound/productUpperBound are '
        + 'wrong when both intervals straddle zero', () => {
            // Every other branch agrees with operator*, so pin both the
            // agreement and the one branch that does not agree.
            check(fc.tuple(interval, interval), ([u, v]) => {
                const ue = u.getEndpoints();
                const ve = v.getEndpoints();
                const lo = FPInterval.productLowerBound(ue, ve);
                const hi = FPInterval.productUpperBound(ue, ve);
                const w = u.mul(v);
                if (ue[0] < 0 && ue[1] > 0 && ve[0] < 0 && ve[1] > 0) {
                    // The defective branch: upstream returns u0*v0 (a product
                    // of two negatives, hence positive) as a lower bound and
                    // u1*v1 as an upper bound.
                    expect(lo).toBe(FPInterval.mul(ue[0], ve[0]).get(0));
                    expect(hi).toBe(FPInterval.mul(ue[1], ve[1]).get(1));
                } else {
                    expect(lo + 0).toBe(w.get(0) + 0);
                    expect(hi + 0).toBe(w.get(1) + 0);
                }
            });
            // The reported witnesses.
            expect(FPInterval.productLowerBound([-1, 1], [-1, 1])).toBe(1);
            expect(new FPInterval(-1, 1).mul(new FPInterval(-1, 1))
                .getEndpoints()).toEqual([-1, 1]);
            expect(FPInterval.productUpperBound([-4, 1], [-4, 1])).toBe(1);
            expect(new FPInterval(-4, 1).mul(new FPInterval(-4, 1))
                .getEndpoints()).toEqual([-4, 16]);
        });
});
