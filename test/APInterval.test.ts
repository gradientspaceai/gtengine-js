import { describe, it, expect } from 'vitest';
import { APInterval, isInfinite } from '../src/APInterval.js';
import { BSRational } from '../src/BSRational.js';
import { check, fc } from './helpers/arbitraries.js';

// Convenience constructors for exact rationals.
function r(numerator: number, denominator: number = 1): BSRational {
    return BSRational.fromNumber(numerator, denominator);
}

const zero = r(0);

function num(x: BSRational): number {
    return x.toNumber();
}

// [w0, w1] as a pair of doubles, for compact comparisons of exact results
// whose values are representable in binary64.
function pair(w: APInterval): [number, number] {
    const e = w.getEndpoints();
    return [e[0].toNumber(), e[1].toNumber()];
}

function contains(w: APInterval, x: BSRational): boolean {
    return w.get(0).lessThanOrEqual(x) && x.lessThanOrEqual(w.get(1));
}

// The exact product (or quotient) interval, computed by brute force from the
// four endpoint combinations. This is the mathematical definition that the
// upstream case analysis implements.
function bruteForce(u: APInterval, v: APInterval,
    op: (a: BSRational, b: BSRational) => BSRational): [number, number] {
    const [u0, u1] = u.getEndpoints();
    const [v0, v1] = v.getEndpoints();
    const values = [op(u0, v0), op(u0, v1), op(u1, v0), op(u1, v1)];
    let lo = values[0], hi = values[0];
    for (const value of values) {
        if (value.lessThan(lo)) { lo = value; }
        if (hi.lessThan(value)) { hi = value; }
    }
    return [lo.toNumber(), hi.toNumber()];
}

describe('APInterval construction and access', () => {
    it('creates [0,0], [e,e] and [e0,e1]', () => {
        expect(pair(new APInterval())).toEqual([0, 0]);
        expect(pair(new APInterval(r(3, 4)))).toEqual([0.75, 0.75]);
        expect(pair(new APInterval(r(-1, 2), r(3, 8)))).toEqual([-0.5, 0.375]);
        expect(pair(APInterval.fromEndpoints([r(1, 4), r(1, 2)]))).toEqual([0.25, 0.5]);
    });

    it('is immutable: endpoints are copied in and out', () => {
        const e0 = r(1, 2);
        const e1 = r(3, 2);
        const u = new APInterval(e0, e1);

        // Mutating the inputs after construction must not change u.
        e0.negate();
        e1.setSign(-1);
        expect(pair(u)).toEqual([0.5, 1.5]);

        // Mutating the outputs must not change u either.
        const out = u.get(0);
        out.negate();
        expect(num(u.get(0))).toBe(0.5);

        const clone = u.clone();
        expect(pair(clone)).toEqual([0.5, 1.5]);
    });

    it('optionally traps an invalid interval', () => {
        expect(() => new APInterval(r(2), r(1))).not.toThrow();
        APInterval.throwOnInvalid = true;
        try {
            expect(() => new APInterval(r(2), r(1))).toThrow(/Invalid interval/);
            expect(() => new APInterval(r(1), r(2))).not.toThrow();
        } finally {
            APInterval.throwOnInvalid = false;
        }
    });
});

describe('APInterval leaf-node operations', () => {
    it('produces degenerate intervals for exact scalar arithmetic', () => {
        // Unlike SWInterval, the arbitrary-precision endpoints are exact, so
        // the leaf-node results are the degenerate intervals [w,w].
        expect(pair(APInterval.add(r(1, 3), r(1, 6)))).toEqual([0.5, 0.5]);
        expect(pair(APInterval.sub(r(1, 3), r(1, 6)))).toEqual([1 / 6, 1 / 6]);
        expect(pair(APInterval.mul(r(2, 3), r(3, 4)))).toEqual([0.5, 0.5]);
        expect(pair(APInterval.div(r(3, 8), r(3, 4)))).toEqual([0.5, 0.5]);
    });

    it('returns the reals for division by zero', () => {
        const w = APInterval.div(r(1), zero);
        expect(isInfinite(w.get(0))).toBe(true);
        expect(isInfinite(w.get(1))).toBe(true);
        expect(w.get(0).getSign()).toBe(-2);
        expect(w.get(1).getSign()).toBe(+2);
        expect(isInfinite(r(1, 3))).toBe(false);
        expect(isInfinite(zero)).toBe(false);
    });

    it('computes the four-argument internal forms', () => {
        // Add(u0,u1,v0,v1) = [u0+v0, u1+v1].
        expect(pair(APInterval.add(r(1), r(2), r(10), r(20)))).toEqual([11, 22]);
        // Sub(u0,u1,v0,v1) = [u0-v1, u1-v0].
        expect(pair(APInterval.sub(r(1), r(2), r(10), r(20)))).toEqual([-19, -8]);
        // Mul(u0,u1,v0,v1) = [u0*v0, u1*v1].
        expect(pair(APInterval.mul(r(1), r(2), r(10), r(20)))).toEqual([10, 40]);
        // Div(u0,u1,v0,v1) = [u0/v1, u1/v0].
        expect(pair(APInterval.div(r(1), r(2), r(10), r(20)))).toEqual([0.05, 0.2]);
        // Mul2 for the both-straddling-zero case.
        expect(pair(APInterval.mul2(r(-2), r(3), r(-5), r(7))))
            .toEqual([-15, 21]);
        // Reciprocal(v0,v1) = [1/v1, 1/v0].
        expect(pair(APInterval.reciprocal(r(2), r(4)))).toEqual([0.25, 0.5]);
    });

    it('builds half-infinite intervals from the reciprocal helpers', () => {
        const down = APInterval.reciprocalDown(r(4));
        expect(num(down.get(0))).toBe(0.25);
        expect(isInfinite(down.get(1))).toBe(true);

        const up = APInterval.reciprocalUp(r(-4));
        expect(isInfinite(up.get(0))).toBe(true);
        expect(num(up.get(1))).toBe(-0.25);
    });
});

describe('APInterval arithmetic', () => {
    it('negates by reversing and negating the endpoints', () => {
        expect(pair(new APInterval(r(-1, 2), r(3, 4)).negate())).toEqual([-0.75, 0.5]);
    });

    it('adds and subtracts exactly, for intervals and scalars', () => {
        const u = new APInterval(r(1, 3), r(1, 2));
        const v = new APInterval(r(1, 6), r(5, 6));
        expect(pair(u.add(v))).toEqual([0.5, 4 / 3]);
        expect(pair(u.sub(v))).toEqual([-0.5, 1 / 3]);
        expect(pair(u.add(r(1, 6)))).toEqual([0.5, 2 / 3]);
        expect(pair(u.sub(r(1, 6)))).toEqual([1 / 6, 1 / 3]);
        expect(pair(APInterval.scalarSub(r(1), u))).toEqual([0.5, 2 / 3]);

        // Scalar addition is commutative: s + v has the same endpoints as
        // v.add(s).
        expect(pair(v.add(r(-2)))).toEqual([num(r(-11, 6)), num(r(-7, 6))]);
    });

    it('multiplies correctly in all nine sign cases', () => {
        const pos = new APInterval(r(2), r(5));
        const neg = new APInterval(r(-5), r(-2));
        const mix = new APInterval(r(-3), r(7));
        const mul = (a: BSRational, b: BSRational) => a.mul(b);
        for (const u of [pos, neg, mix]) {
            for (const v of [pos, neg, mix]) {
                expect(pair(u.mul(v))).toEqual(bruteForce(u, v, mul));
            }
        }
        // Zero-touching endpoints exercise the >= and <= boundaries.
        const touch = new APInterval(zero, r(4));
        expect(pair(touch.mul(mix))).toEqual(bruteForce(touch, mix, mul));
        expect(pair(mix.mul(touch))).toEqual(bruteForce(mix, touch, mul));
    });

    it('multiplies by a scalar on either side', () => {
        const u = new APInterval(r(-3), r(7));
        expect(pair(u.mul(r(2)))).toEqual([-6, 14]);
        expect(pair(u.mul(r(-2)))).toEqual([-14, 6]);
        expect(pair(u.mul(zero))).toEqual([0, 0]);
    });

    it('divides by intervals and scalars that exclude zero', () => {
        const div = (a: BSRational, b: BSRational) => a.div(b);
        const us = [new APInterval(r(2), r(5)), new APInterval(r(-5), r(-2)),
            new APInterval(r(-3), r(7))];
        const vs = [new APInterval(r(2), r(4)), new APInterval(r(-4), r(-2))];
        for (const u of us) {
            for (const v of vs) {
                expect(pair(u.div(v))).toEqual(bruteForce(u, v, div));
                expect(pair(APInterval.scalarDiv(r(3), v)))
                    .toEqual(bruteForce(new APInterval(r(3)), v, div));
            }
        }
        expect(pair(us[0].div(r(2)))).toEqual([1, 2.5]);
        expect(pair(us[0].div(r(-2)))).toEqual([-2.5, -1]);
    });

    it('returns indeterminate results when the divisor straddles zero', () => {
        const u = new APInterval(r(1), r(2));
        const straddle = new APInterval(r(-1), r(3));
        for (const w of [u.div(zero), u.div(straddle),
            APInterval.scalarDiv(r(1), straddle)]) {
            expect(isInfinite(w.get(0))).toBe(true);
            expect(isInfinite(w.get(1))).toBe(true);
        }

        // A divisor with a zero endpoint gives a half-infinite result whose
        // finite endpoint is exact.
        const upper = u.div(new APInterval(zero, r(4)));
        expect(num(upper.get(0))).toBe(0.25);
        expect(isInfinite(upper.get(1))).toBe(true);

        const lower = u.div(new APInterval(r(-4), zero));
        expect(isInfinite(lower.get(0))).toBe(true);
        expect(num(lower.get(1))).toBe(-0.25);
    });
});

describe('APInterval containment invariants', () => {
    it('encloses the exact result for randomized exact inputs', () => {
        // A deterministic linear congruential generator.
        let seed = 987654321;
        const next = (bound: number): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed % bound;
        };
        const randomRational = (): BSRational => {
            return r(next(21) - 10, next(7) + 1);
        };
        const randomInterval = (): APInterval => {
            const a = randomRational();
            const b = randomRational();
            return a.lessThanOrEqual(b) ? new APInterval(a, b) : new APInterval(b, a);
        };
        // The endpoints and the midpoint are exact members of the interval.
        const samples = (w: APInterval): BSRational[] => {
            const [e0, e1] = w.getEndpoints();
            return [e0, e1, BSRational.ldexp(e0.add(e1), -1)];
        };

        for (let trial = 0; trial < 40; ++trial) {
            const u = randomInterval();
            const v = randomInterval();
            const sum = u.add(v);
            const difference = u.sub(v);
            const product = u.mul(v);
            const divisorExcludesZero =
                v.get(0).getSign() > 0 || v.get(1).getSign() < 0;
            const quotient = divisorExcludesZero ? u.div(v) : null;
            for (const x of samples(u)) {
                for (const y of samples(v)) {
                    expect(contains(sum, x.add(y))).toBe(true);
                    expect(contains(difference, x.sub(y))).toBe(true);
                    expect(contains(product, x.mul(y))).toBe(true);
                    if (quotient !== null) {
                        expect(contains(quotient, x.div(y))).toBe(true);
                    }
                }
            }
            // The endpoints of a computed interval are ordered.
            for (const w of [sum, difference, product]) {
                expect(w.get(0).lessThanOrEqual(w.get(1))).toBe(true);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md): the upstream case analysis of
// operator* and operator/ is checked against the exact hull of the endpoint
// combinations, and every result is checked to contain the exact result of
// the operation applied to interior points.
// ---------------------------------------------------------------------------

describe('APInterval verification', () => {
    const rational: fc.Arbitrary<BSRational> =
        fc.tuple(fc.integer({ min: -30, max: 30 }),
            fc.integer({ min: 1, max: 16 }))
            .map(([n, d]) => BSRational.fromNumber(n, d));

    const interval: fc.Arbitrary<APInterval> =
        fc.tuple(rational, rational).map(([a, b]) =>
            a.lessThanOrEqual(b) ? new APInterval(a, b) : new APInterval(b, a));

    // Intervals that do not contain zero, so that division is determinate.
    const nonzeroInterval: fc.Arbitrary<APInterval> = interval.filter(w => {
        const e = w.getEndpoints();
        return e[0].getSign() > 0 || e[1].getSign() < 0;
    });

    // The exact hull of the four endpoint combinations, which is the
    // mathematical range of a bilinear operation on the two intervals.
    function hull(u: APInterval, v: APInterval,
        op: (a: BSRational, b: BSRational) => BSRational): [BSRational, BSRational] {
        const [u0, u1] = u.getEndpoints();
        const [v0, v1] = v.getEndpoints();
        const values = [op(u0, v0), op(u0, v1), op(u1, v0), op(u1, v1)];
        let lo = values[0], hi = values[0];
        for (const value of values) {
            if (value.lessThan(lo)) { lo = value; }
            if (hi.lessThan(value)) { hi = value; }
        }
        return [lo, hi];
    }

    function expectSame(w: APInterval, expected: [BSRational, BSRational]): void {
        const e = w.getEndpoints();
        expect(e[0].equals(expected[0])).toBe(true);
        expect(e[1].equals(expected[1])).toBe(true);
    }

    // A point of the interval: u0 + t*(u1 - u0) with t in [0,1] rational.
    function interior(w: APInterval, t: BSRational): BSRational {
        const e = w.getEndpoints();
        return e[0].add(t.mul(e[1].sub(e[0])));
    }

    const fraction: fc.Arbitrary<BSRational> =
        fc.integer({ min: 0, max: 8 }).map(i => BSRational.fromNumber(i, 8));

    it('add and sub are exactly the hull of the endpoint combinations', () => {
        check(fc.tuple(interval, interval), ([u, v]) => {
            expectSame(u.add(v), hull(u, v, (a, b) => a.add(b)));
            expectSame(u.sub(v), hull(u, v, (a, b) => a.sub(b)));
            // Scalar overloads agree with the degenerate interval.
            const s = v.get(0);
            expectSame(u.add(s), hull(u, new APInterval(s), (a, b) => a.add(b)));
            expectSame(u.sub(s), hull(u, new APInterval(s), (a, b) => a.sub(b)));
            expectSame(APInterval.scalarSub(s, u),
                hull(new APInterval(s), u, (a, b) => a.sub(b)));
        });
    });

    it('the nine sign cases of operator* reproduce the exact hull', () => {
        check(fc.tuple(interval, interval), ([u, v]) => {
            const w = u.mul(v);
            expectSame(w, hull(u, v, (a, b) => a.mul(b)));
            // Multiplication of intervals is commutative.
            expectSame(v.mul(u), w.getEndpoints());
            // The scalar overload agrees with the degenerate interval.
            const s = v.get(0);
            expectSame(u.mul(s), hull(u, new APInterval(s), (a, b) => a.mul(b)));
        });
    });

    it('operator/ reproduces the exact hull when the divisor excludes zero',
        () => {
            check(fc.tuple(interval, nonzeroInterval), ([u, v]) => {
                expectSame(u.div(v), hull(u, v, (a, b) => a.div(b)));
                expectSame(APInterval.scalarDiv(u.get(0), v),
                    hull(new APInterval(u.get(0)), v, (a, b) => a.div(b)));
                const s = v.get(0);
                expectSame(u.div(s),
                    hull(u, new APInterval(s), (a, b) => a.div(b)));
            });
        });

    it('contains the exact result for every pair of interior points', () => {
        check(fc.tuple(interval, nonzeroInterval, fraction, fraction),
            ([u, v, s, t]) => {
                const x = interior(u, s);
                const y = interior(v, t);
                for (const [w, value] of [
                    [u.add(v), x.add(y)],
                    [u.sub(v), x.sub(y)],
                    [u.mul(v), x.mul(y)],
                    [u.div(v), x.div(y)]
                ] as [APInterval, BSRational][]) {
                    const e = w.getEndpoints();
                    expect(e[0].lessThanOrEqual(value)).toBe(true);
                    expect(value.lessThanOrEqual(e[1])).toBe(true);
                }
            });
    });

    it('negation reverses and negates the endpoints', () => {
        check(interval, u => {
            const e = u.getEndpoints();
            const n = u.negate();
            expect(n.get(0).equals(e[1].negated())).toBe(true);
            expect(n.get(1).equals(e[0].negated())).toBe(true);
            // Negation is an involution and turns * into its negation.
            expect(n.negate().get(0).equals(e[0])).toBe(true);
            expect(n.negate().get(1).equals(e[1])).toBe(true);
        });
    });

    it('the leaf-node operations are the degenerate exact results', () => {
        check(fc.tuple(rational, rational), ([a, b]) => {
            for (const [w, value] of [
                [APInterval.add(a, b), a.add(b)],
                [APInterval.sub(a, b), a.sub(b)],
                [APInterval.mul(a, b), a.mul(b)]
            ] as [APInterval, BSRational][]) {
                expect(w.get(0).equals(value)).toBe(true);
                expect(w.get(1).equals(value)).toBe(true);
            }
            const q = APInterval.div(a, b);
            if (b.getSign() !== 0) {
                expect(q.get(0).equals(a.div(b))).toBe(true);
                expect(q.get(1).equals(a.div(b))).toBe(true);
            } else {
                expect(isInfinite(q.get(0))).toBe(true);
                expect(isInfinite(q.get(1))).toBe(true);
            }
        });
    });

    it('yields the indeterminate interval exactly where upstream does', () => {
        check(fc.tuple(interval, interval), ([u, v]) => {
            const e = v.getEndpoints();
            const straddles = e[0].getSign() < 0 && e[1].getSign() > 0;
            const w = u.div(v);
            if (straddles) {
                expect(isInfinite(w.get(0))).toBe(true);
                expect(isInfinite(w.get(1))).toBe(true);
            } else if (e[0].getSign() === 0 && e[1].getSign() === 0) {
                // [0,0]: reciprocalDown of 0 divides by zero.
                expect(() => u.div(v)).toThrow();
            } else {
                // A divisor with one zero endpoint gives a half-line. For a
                // nonnegative numerator interval with a positive upper bound
                // the finite endpoint is the exact quotient by the nonzero
                // divisor endpoint and the other endpoint is a sentinel.
                const uu = u.getEndpoints();
                if (uu[0].getSign() >= 0 && uu[1].getSign() > 0) {
                    if (e[0].getSign() === 0 && e[1].getSign() > 0) {
                        // v = [0,b]: reciprocal is [1/b, +infinity).
                        expect(w.get(0).equals(uu[0].div(e[1]))).toBe(true);
                        expect(isInfinite(w.get(1))).toBe(true);
                    } else if (e[1].getSign() === 0 && e[0].getSign() < 0) {
                        // v = [a,0]: reciprocal is (-infinity, 1/a].
                        expect(isInfinite(w.get(0))).toBe(true);
                        expect(w.get(1).equals(uu[0].div(e[0]))).toBe(true);
                    }
                }
            }
            // Division by the zero scalar is the whole real line.
            const z = u.div(new BSRational());
            expect(isInfinite(z.get(0))).toBe(true);
            expect(isInfinite(z.get(1))).toBe(true);
        });
    });

    it('is immutable: no accessor exposes the internal endpoints', () => {
        check(fc.tuple(rational, rational), ([a, b]) => {
            const e0 = a.lessThanOrEqual(b) ? a : b;
            const e1 = a.lessThanOrEqual(b) ? b : a;
            const u = new APInterval(e0, e1);
            const before = u.getEndpoints();

            // Mutating the constructor inputs, a get() result and a
            // getEndpoints() result must all leave u unchanged.
            e0.negate();
            e1.setSign(-1);
            u.get(0).negate();
            const taken = u.getEndpoints();
            taken[0].setSign(1);
            taken[1].negate();

            const after = u.getEndpoints();
            expect(after[0].equals(before[0])).toBe(true);
            expect(after[1].equals(before[1])).toBe(true);
        });
    });
});
