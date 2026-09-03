import { describe, it, expect } from 'vitest';
import { check, fc } from './helpers/arbitraries.js';
import { QFNumber, type QFCoefficient } from '../src/QFNumber.js';

// Numeric evaluation of a quadratic field number (recursively for nested
// coefficients), used as an independent cross-check.
function evalC(c: QFCoefficient): number {
    return typeof c === 'number' ? c : evalQF(c);
}

function evalQF(q: QFNumber): number {
    return evalC(q.x[0]) + evalC(q.x[1]) * Math.sqrt(q.d);
}

// Deterministic pseudorandom generator (LCG) for reproducible tests.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function randomInt(rand: () => number, lo: number, hi: number): number {
    return lo + Math.floor(rand() * (hi - lo + 1));
}

// Numeric coefficient equality (=== treats -0 and +0 as equal, matching the
// C++ operator== on double).
function expectSameCoeffs(p: QFNumber, q: QFNumber): void {
    expect(p.x[0] === q.x[0]).toBe(true);
    expect(p.x[1] === q.x[1]).toBe(true);
}

describe('QFNumber construction', () => {
    it('default-constructs to 0 + 0 * sqrt(0)', () => {
        const q = new QFNumber();
        expect(q.x[0]).toBe(0);
        expect(q.x[1]).toBe(0);
        expect(q.d).toBe(0);
    });

    it('fromD creates 0 + 0 * sqrt(d)', () => {
        const q = QFNumber.fromD(5);
        expect(q.x[0]).toBe(0);
        expect(q.x[1]).toBe(0);
        expect(q.d).toBe(5);
    });

    it('stores coefficients and d', () => {
        const q = new QFNumber(3, -2, 7);
        expect(q.x[0]).toBe(3);
        expect(q.x[1]).toBe(-2);
        expect(q.d).toBe(7);
    });

    it('clone produces an independent deep copy', () => {
        const inner = new QFNumber(1, 2, 2);
        const q = new QFNumber(inner, new QFNumber(3, 4, 2), 3);
        const c = q.clone();
        expect(evalQF(c)).toBe(evalQF(q));
        (c.x[0] as QFNumber).x[0] = 99;
        expect((q.x[0] as QFNumber).x[0]).toBe(1);
    });

    it('implements the ArbitraryPrecisionNumber marker with division', () => {
        const q = new QFNumber(1, 1, 2);
        expect(q.isArbitraryPrecision).toBe(true);
        expect(q.hasDivisionOperator).toBe(true);
    });
});

describe('QFNumber known values', () => {
    it('(1 + sqrt(2)) * (1 - sqrt(2)) = -1', () => {
        const a = new QFNumber(1, 1, 2);
        const b = new QFNumber(1, -1, 2);
        const p = a.mul(b);
        expect(p.x[0]).toBe(-1);
        expect(p.x[1]).toBe(0);
        expect(p.d).toBe(2);
    });

    it('sqrt(2) * sqrt(2) = 2', () => {
        const r2 = new QFNumber(0, 1, 2);
        const p = r2.mul(r2);
        expect(p.x[0]).toBe(2);
        expect(p.x[1]).toBe(0);
    });

    it('golden ratio satisfies phi^2 = phi + 1', () => {
        const phi = new QFNumber(0.5, 0.5, 5);
        const sq = phi.mul(phi);
        const inc = phi.add(1);
        expect(sq.x[0]).toBe(inc.x[0]);
        expect(sq.x[1]).toBe(inc.x[1]);
        expect(evalQF(phi)).toBeCloseTo((1 + Math.sqrt(5)) / 2, 12);
    });

    it('q / q = 1 exactly', () => {
        const q = new QFNumber(3, -5, 7);
        const r = q.div(q);
        expect(r.x[0]).toBe(1);
        // The x[1] numerator is exactly zero; the sign of the zero depends
        // on the sign of the denominator (IEEE 0 / -166 = -0), as in C++.
        expect(r.x[1] === 0).toBe(true);
    });

    it('scalar operations match upstream formulas', () => {
        const q = new QFNumber(1, 2, 3);
        const a = q.add(2);
        expect([a.x[0], a.x[1], a.d]).toEqual([3, 2, 3]);
        const s = q.sub(3);
        expect([s.x[0], s.x[1], s.d]).toEqual([-2, 2, 3]);
        const m = q.mul(2);
        expect([m.x[0], m.x[1], m.d]).toEqual([2, 4, 3]);
        const v = q.div(2);
        expect([v.x[0], v.x[1], v.d]).toEqual([0.5, 1, 3]);
        const ss = QFNumber.scalarSub(5, q);
        expect([ss.x[0], ss.x[1], ss.d]).toEqual([4, -2, 3]);
        // 1 / (1 + 2*sqrt(3)): denom = 1 - 4*3 = -11.
        const sd = QFNumber.scalarDiv(1, q);
        expect(sd.x[0]).toBeCloseTo(-1 / 11, 15);
        expect(sd.x[1]).toBeCloseTo(2 / 11, 15);
        expect(evalQF(q) * evalQF(sd)).toBeCloseTo(1, 12);
    });

    it('negate flips both coefficients', () => {
        const q = new QFNumber(1, -2, 3);
        const n = q.negate();
        expect([n.x[0], n.x[1], n.d]).toEqual([-1, 2, 3]);
        const z = q.add(n);
        expect(z.x[0]).toBe(0);
        expect(z.x[1]).toBe(0);
    });
});

describe('QFNumber field axioms (randomized, exact integer coefficients)', () => {
    const rand = makeRandom(0x0b13);
    const trials = 100;

    it('satisfies commutativity, associativity and distributivity', () => {
        for (let trial = 0; trial < trials; ++trial) {
            const d = [2, 3, 5][randomInt(rand, 0, 2)];
            const a = new QFNumber(randomInt(rand, -10, 10), randomInt(rand, -10, 10), d);
            const b = new QFNumber(randomInt(rand, -10, 10), randomInt(rand, -10, 10), d);
            const c = new QFNumber(randomInt(rand, -10, 10), randomInt(rand, -10, 10), d);

            // Commutativity (exact for integer coefficients).
            expectSameCoeffs(a.add(b), b.add(a));
            expectSameCoeffs(a.mul(b), b.mul(a));

            // Associativity.
            expectSameCoeffs(a.add(b).add(c), a.add(b.add(c)));
            expectSameCoeffs(a.mul(b).mul(c), a.mul(b.mul(c)));

            // Distributivity.
            expectSameCoeffs(a.mul(b.add(c)), a.mul(b).add(a.mul(c)));
        }
    });

    it('has additive and multiplicative identities and inverses', () => {
        for (let trial = 0; trial < trials; ++trial) {
            const d = [2, 3, 5][randomInt(rand, 0, 2)];
            let a = new QFNumber(randomInt(rand, -10, 10), randomInt(rand, -10, 10), d);
            if (a.x[0] === 0 && a.x[1] === 0) {
                a = new QFNumber(1, 1, d);
            }
            const zero = QFNumber.fromD(d);
            const one = new QFNumber(1, 0, d);

            expect(a.add(zero).x).toEqual(a.x);
            expect(a.mul(one).x).toEqual(a.x);

            const negA = a.negate();
            const sumZero = a.add(negA);
            expect(sumZero.x[0]).toBe(0);
            expect(sumZero.x[1]).toBe(0);

            // Multiplicative inverse (nonzero since sqrt(d) is irrational
            // and the coefficients are integers).
            const inv = QFNumber.scalarDiv(1, a);
            const prod = a.mul(inv);
            expect(prod.x[0] as number).toBeCloseTo(1, 10);
            expect(prod.x[1] as number).toBeCloseTo(0, 10);
        }
    });

    it('div is the inverse of mul', () => {
        for (let trial = 0; trial < trials; ++trial) {
            const d = [2, 3, 5][randomInt(rand, 0, 2)];
            const a = new QFNumber(randomInt(rand, -10, 10), randomInt(rand, -10, 10), d);
            let b = new QFNumber(randomInt(rand, -10, 10), randomInt(rand, -10, 10), d);
            if (b.x[0] === 0 && b.x[1] === 0) {
                b = new QFNumber(2, 1, d);
            }
            const q = a.mul(b).div(b);
            expect(q.x[0] as number).toBeCloseTo(a.x[0] as number, 10);
            expect(q.x[1] as number).toBeCloseTo(a.x[1] as number, 10);
        }
    });
});

describe('QFNumber comparisons', () => {
    it('agrees with numeric evaluation for random operands', () => {
        const rand = makeRandom(0x517e);
        for (let trial = 0; trial < 200; ++trial) {
            const d = [2, 3, 5][randomInt(rand, 0, 2)];
            const a = new QFNumber(randomInt(rand, -10, 10), randomInt(rand, -10, 10), d);
            const b = new QFNumber(randomInt(rand, -10, 10), randomInt(rand, -10, 10), d);
            const va = evalQF(a), vb = evalQF(b);
            expect(a.equals(b)).toBe(va === vb);
            expect(a.notEquals(b)).toBe(va !== vb);
            expect(a.lessThan(b)).toBe(va < vb);
            expect(a.greaterThan(b)).toBe(va > vb);
            expect(a.lessThanEqual(b)).toBe(va <= vb);
            expect(a.greaterThanEqual(b)).toBe(va >= vb);
        }
    });

    it('detects equality across representations when d is a perfect square', () => {
        // 1 + 1 * sqrt(4) = 3 + 0 * sqrt(4).
        const a = new QFNumber(1, 1, 4);
        const b = new QFNumber(3, 0, 4);
        expect(a.equals(b)).toBe(true);
        expect(a.lessThan(b)).toBe(false);
        expect(b.lessThan(a)).toBe(false);
    });

    it('compares by x[0] alone when d = 0', () => {
        const a = new QFNumber(1, 7, 0);
        const b = new QFNumber(1, -9, 0);
        expect(a.equals(b)).toBe(true);
        expect(new QFNumber(0, 5, 0).lessThan(a)).toBe(true);
    });

    it('orders known irrational values correctly', () => {
        const two = new QFNumber(2, 0, 2);
        const onePlusRoot2 = new QFNumber(1, 1, 2); // ~2.414
        expect(two.lessThan(onePlusRoot2)).toBe(true);
        expect(onePlusRoot2.greaterThan(two)).toBe(true);
        expect(onePlusRoot2.equals(onePlusRoot2)).toBe(true);
    });
});

describe('QFNumber nested fields (N = 2)', () => {
    it('(sqrt(2) + sqrt(3))^2 = 5 + 2*sqrt(6)', () => {
        // z = sqrt(2) + sqrt(3) as (0 + 1*sqrt(2)) + (1 + 0*sqrt(2))*sqrt(3).
        const z = new QFNumber(new QFNumber(0, 1, 2), new QFNumber(1, 0, 2), 3);
        expect(evalQF(z)).toBeCloseTo(Math.sqrt(2) + Math.sqrt(3), 12);
        const zz = z.mul(z);
        expect(evalQF(zz)).toBeCloseTo(5 + 2 * Math.sqrt(6), 12);
        // z^2 = (2 + 3) + (2*sqrt(2))*sqrt(3): check the exact coefficients.
        const x0 = zz.x[0] as QFNumber;
        const x1 = zz.x[1] as QFNumber;
        expect([x0.x[0], x0.x[1]]).toEqual([5, 0]);
        expect([x1.x[0], x1.x[1]]).toEqual([0, 2]);
    });

    it('nested arithmetic agrees with numeric evaluation', () => {
        const rand = makeRandom(0xbeef);
        for (let trial = 0; trial < 50; ++trial) {
            const mk = () => new QFNumber(
                new QFNumber(randomInt(rand, -5, 5), randomInt(rand, -5, 5), 2),
                new QFNumber(randomInt(rand, -5, 5), randomInt(rand, -5, 5), 2),
                3);
            const a = mk();
            const b = mk();
            const va = evalQF(a), vb = evalQF(b);
            expect(evalQF(a.add(b))).toBeCloseTo(va + vb, 10);
            expect(evalQF(a.sub(b))).toBeCloseTo(va - vb, 10);
            expect(evalQF(a.mul(b))).toBeCloseTo(va * vb, 8);
            if (Math.abs(vb) > 0.01) {
                expect(evalQF(a.div(b))).toBeCloseTo(va / vb, 8);
            }
            expect(evalQF(a.negate())).toBeCloseTo(-va, 12);
            expect(evalQF(a.add(3))).toBeCloseTo(va + 3, 12);
            expect(evalQF(QFNumber.scalarSub(3, a))).toBeCloseTo(3 - va, 12);
        }
    });
});

describe('QFNumber mismatched-d trapping', () => {
    it('does not throw by default', () => {
        const a = new QFNumber(1, 1, 2);
        const b = new QFNumber(1, 1, 3);
        expect(() => a.add(b)).not.toThrow();
    });

    it('throws when assertOnMismatchedD is enabled', () => {
        const a = new QFNumber(1, 1, 2);
        const b = new QFNumber(1, 1, 3);
        QFNumber.assertOnMismatchedD = true;
        try {
            expect(() => a.add(b)).toThrow('Mismatched d-value.');
            expect(() => a.sub(b)).toThrow('Mismatched d-value.');
            expect(() => a.mul(b)).toThrow('Mismatched d-value.');
            expect(() => a.div(b)).toThrow('Mismatched d-value.');
            expect(() => a.equals(b)).toThrow('Mismatched d-value.');
            expect(() => a.lessThan(b)).toThrow('Mismatched d-value.');
            expect(() => a.add(a.clone())).not.toThrow();
        } finally {
            QFNumber.assertOnMismatchedD = false;
        }
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V05). Quadratic field elements with integer coefficients
// are compared against exact bigint symbolic arithmetic on a + b*sqrt(d).
// ---------------------------------------------------------------------------
describe('QFNumber verification', () => {
    // Small integer coefficients keep every product exact in binary64.
    const coeff = fc.integer({ min: -12, max: 12 });
    // d = 0 exercises the degenerate branch, 4 and 9 are perfect squares (so
    // sqrt(d) is rational and distinct coefficient pairs can be equal), the
    // rest are squarefree.
    const dValue = fc.constantFrom(0, 2, 3, 4, 5, 6, 7, 9);

    interface Sym { a: bigint; b: bigint; d: bigint; }

    // C++ operator== on double treats -0 and +0 as equal; the coefficient
    // arithmetic legitimately produces -0 (for example -x with x = 0), so
    // coefficients are compared with === rather than Object.is.
    const expectCoeffs = (q: QFNumber, a: number, b: number): void => {
        expect(q.x[0] === a).toBe(true);
        expect(q.x[1] === b).toBe(true);
    };

    const expectSameX = (p0: QFNumber, p1: QFNumber): void => {
        expect(p0.x[0] === p1.x[0]).toBe(true);
        expect(p0.x[1] === p1.x[1]).toBe(true);
    };

    const symOf = (q: QFNumber): Sym => ({
        a: BigInt(q.x[0] as number), b: BigInt(q.x[1] as number), d: BigInt(q.d)
    });

    /** Exact sign of a + b*sqrt(d) for integer a, b and d >= 0. */
    function exactSign(a: bigint, b: bigint, d: bigint): number {
        if (d === 0n || b === 0n) {
            return a < 0n ? -1 : (a > 0n ? 1 : 0);
        }
        if (a === 0n) {
            return b < 0n ? -1 : 1;
        }
        if (a > 0n && b > 0n) { return 1; }
        if (a < 0n && b < 0n) { return -1; }
        // Opposite signs: compare the squares of the two terms.
        const lhs = a * a;
        const rhs = b * b * d;
        if (lhs === rhs) { return 0; }
        return (lhs > rhs) === (a > 0n) ? 1 : -1;
    }

    const qfArb = fc.tuple(coeff, coeff, dValue)
        .map(([a, b, d]) => new QFNumber(a, b, d));

    /** Two numbers of the same quadratic field. */
    const pairArb = fc.tuple(coeff, coeff, coeff, coeff, dValue)
        .map(([a0, b0, a1, b1, d]) =>
            [new QFNumber(a0, b0, d), new QFNumber(a1, b1, d)] as const);

    // ---- exact arithmetic --------------------------------------------------

    it('add, sub, negate and mul match exact bigint symbolic arithmetic', () => {
        check(pairArb, ([q0, q1]) => {
            const s0 = symOf(q0);
            const s1 = symOf(q1);

            const sum = q0.add(q1);
            expectCoeffs(sum, Number(s0.a + s1.a), Number(s0.b + s1.b));
            expect(sum.d).toBe(q0.d);

            const diff = q0.sub(q1);
            expectCoeffs(diff, Number(s0.a - s1.a), Number(s0.b - s1.b));

            const prod = q0.mul(q1);
            expectCoeffs(prod, Number(s0.a * s1.a + s0.b * s1.b * s0.d),
                Number(s0.a * s1.b + s0.b * s1.a));
            expect(prod.d).toBe(q0.d);

            const neg = q0.negate();
            expectCoeffs(neg, Number(-s0.a), Number(-s0.b));
            expect(neg.d).toBe(q0.d);
            return true;
        });
    });

    it('the scalar overloads match the upstream formulas', () => {
        check(fc.tuple(qfArb, coeff), ([q, s]) => {
            const sym = symOf(q);
            const bs = BigInt(s);

            expectCoeffs(q.add(s), Number(sym.a + bs), Number(sym.b));
            expectCoeffs(q.sub(s), Number(sym.a - bs), Number(sym.b));
            expectCoeffs(q.mul(s), Number(sym.a * bs), Number(sym.b * bs));
            // s - q keeps the sqrt coefficient negated.
            expectCoeffs(QFNumber.scalarSub(s, q), Number(bs - sym.a), Number(-sym.b));
            return true;
        });
    });

    it('multiplication is commutative, associative and distributes over addition', () => {
        check(fc.tuple(coeff, coeff, coeff, coeff, coeff, coeff, dValue),
            ([a0, b0, a1, b1, a2, b2, d]) => {
                const p = new QFNumber(a0, b0, d);
                const q = new QFNumber(a1, b1, d);
                const r = new QFNumber(a2, b2, d);
                expectSameX(p.mul(q), q.mul(p));
                expectSameX(p.add(q), q.add(p));
                expectSameX(p.mul(q).mul(r), p.mul(q.mul(r)));
                expectSameX(p.add(q).add(r), p.add(q.add(r)));
                expectSameX(p.mul(q.add(r)), p.mul(q).add(p.mul(r)));
                return true;
            });
    });

    // ---- division ----------------------------------------------------------

    it('div inverts mul and scalarDiv agrees with the QFNumber division', () => {
        check(pairArb, ([q0, q1]) => {
            const s1 = symOf(q1);
            // The divisor must have a nonzero norm a^2 - b^2 d.
            const norm = s1.a * s1.a - s1.b * s1.b * s1.d;
            if (norm === 0n) { return true; }

            const quotient = q0.mul(q1).div(q1);
            expect(quotient.x[0] as number).toBeCloseTo(q0.x[0] as number, 8);
            expect(quotient.x[1] as number).toBeCloseTo(q0.x[1] as number, 8);

            // s / q equals (s + 0*sqrt(d)) / q.
            const viaDiv = new QFNumber(3, 0, q1.d).div(q1);
            const viaScalar = QFNumber.scalarDiv(3, q1);
            expect(viaScalar.x[0] as number).toBeCloseTo(viaDiv.x[0] as number, 10);
            expect(viaScalar.x[1] as number).toBeCloseTo(viaDiv.x[1] as number, 10);
            return true;
        });
    });

    // Upstream uses q0.d (the dividend's d), not q1.d, in the denominator and
    // in numer0. With mismatched d-values that choice is observable; it is
    // pinned here because the trap is off by default.
    it('div uses the dividend d-value, as upstream does', () => {
        const q0 = new QFNumber(1, 1, 2);
        const q1 = new QFNumber(3, 1, 5);
        const result = q0.div(q1);
        // denom = 3*3 - 1*1*2 = 7, numer0 = 1*3 - 1*1*2 = 1, numer1 = 1*3 - 1*1 = 2.
        expect(result.x[0]).toBe(1 / 7);
        expect(result.x[1]).toBe(2 / 7);
        expect(result.d).toBe(2);
    });

    // ---- comparisons -------------------------------------------------------

    it('the comparisons realise the exact order of a + b*sqrt(d)', () => {
        check(pairArb, ([q0, q1]) => {
            const s0 = symOf(q0);
            const s1 = symOf(q1);
            const c = exactSign(s0.a - s1.a, s0.b - s1.b, s0.d);

            expect(q0.lessThan(q1)).toBe(c < 0);
            expect(q0.greaterThan(q1)).toBe(c > 0);
            expect(q0.equals(q1)).toBe(c === 0);
            expect(q0.notEquals(q1)).toBe(c !== 0);
            expect(q0.lessThanEqual(q1)).toBe(c <= 0);
            expect(q0.greaterThanEqual(q1)).toBe(c >= 0);
            return true;
        });
    });

    it('the order is a strict weak ordering consistent with numeric evaluation', () => {
        check(fc.tuple(coeff, coeff, coeff, coeff, coeff, coeff, dValue),
            ([a0, b0, a1, b1, a2, b2, d]) => {
                const p = new QFNumber(a0, b0, d);
                const q = new QFNumber(a1, b1, d);
                const r = new QFNumber(a2, b2, d);
                // Irreflexive and antisymmetric.
                expect(p.lessThan(p)).toBe(false);
                expect(p.lessThan(q) && q.lessThan(p)).toBe(false);
                // Transitive.
                if (p.lessThan(q) && q.lessThan(r)) {
                    expect(p.lessThan(r)).toBe(true);
                }
                // Consistent with the floating-point evaluation, except in the
                // near-tie cases that the exact comparison resolves and the
                // double evaluation cannot.
                const vp = evalQF(p), vq = evalQF(q);
                if (Math.abs(vp - vq) > 1e-9) {
                    expect(p.lessThan(q)).toBe(vp < vq);
                }
                return true;
            });
    });

    it('perfect-square d identifies numerically equal representations', () => {
        // 1 + 2*sqrt(9) = 7 = 7 + 0*sqrt(9).
        expect(new QFNumber(1, 2, 9).equals(new QFNumber(7, 0, 9))).toBe(true);
        // 5 - 1*sqrt(4) = 3 = 1 + 1*sqrt(4).
        expect(new QFNumber(5, -1, 4).equals(new QFNumber(1, 1, 4))).toBe(true);
        // A squarefree d never allows that.
        expect(new QFNumber(1, 2, 5).equals(new QFNumber(7, 0, 5))).toBe(false);
    });

    // ---- value semantics ---------------------------------------------------

    // Regression: upstream stores coefficients by value, so no result may
    // share a nested coefficient object with the operands that produced it.
    it('results never alias the nested coefficients of their operands', () => {
        const inner = (a: number, b: number): QFNumber => new QFNumber(a, b, 2);
        const q0 = new QFNumber(inner(1, 2), inner(3, 4), 5);
        const q1 = new QFNumber(inner(5, 6), inner(7, 8), 5);

        const results: QFNumber[] = [
            q0.clone(), q0.negate(), q0.add(1), q0.sub(1), q0.mul(2), q0.div(2),
            q0.add(q1), q0.sub(q1), q0.mul(q1), QFNumber.scalarSub(1, q0)
        ];
        for (const r of results) {
            for (const c of r.x) {
                expect(c).not.toBe(q0.x[0]);
                expect(c).not.toBe(q0.x[1]);
                expect(c).not.toBe(q1.x[0]);
                expect(c).not.toBe(q1.x[1]);
            }
        }

        // Mutating a result must not disturb the operands.
        const sum = q0.add(1);
        (sum.x[1] as QFNumber).x[0] = 99;
        expect((q0.x[1] as QFNumber).x[0]).toBe(3);

        // The constructor copies too.
        const c0 = inner(10, 11);
        const held = new QFNumber(c0, inner(12, 13), 5);
        expect(held.x[0]).not.toBe(c0);
        c0.x[0] = -1;
        expect((held.x[0] as QFNumber).x[0]).toBe(10);
    });

    // ---- nested fields (N = 2) ---------------------------------------------

    it('nested arithmetic agrees with numeric evaluation', () => {
        const nested = fc.tuple(coeff, coeff, coeff, coeff)
            .map(([a, b, c, e]) =>
                new QFNumber(new QFNumber(a, b, 2), new QFNumber(c, e, 2), 3));
        check(fc.tuple(nested, nested), ([p, q]) => {
            expect(evalQF(p.add(q))).toBeCloseTo(evalQF(p) + evalQF(q), 9);
            expect(evalQF(p.sub(q))).toBeCloseTo(evalQF(p) - evalQF(q), 9);
            expect(evalQF(p.mul(q))).toBeCloseTo(evalQF(p) * evalQF(q), 8);
            expect(evalQF(p.negate())).toBeCloseTo(-evalQF(p), 9);
            expect(evalQF(p.mul(3))).toBeCloseTo(3 * evalQF(p), 9);
            return true;
        });
    });

    it('(sqrt(2) + sqrt(3))^2 = 5 + 2 sqrt(6) in the nested field', () => {
        // z = 0 + 1*sqrt(2) + (1 + 0*sqrt(2))*sqrt(3).
        const z = new QFNumber(new QFNumber(0, 1, 2), new QFNumber(1, 0, 2), 3);
        const zz = z.mul(z);
        // x[0] = 2 + 3 = 5, x[1] = 0 + 2*sqrt(2)*... evaluate numerically.
        expect(evalQF(zz)).toBeCloseTo(5 + 2 * Math.sqrt(6), 12);
        expect((zz.x[0] as QFNumber).x[0]).toBe(5);
        expect((zz.x[0] as QFNumber).x[1]).toBe(0);
        expect((zz.x[1] as QFNumber).x[0]).toBe(0);
        expect((zz.x[1] as QFNumber).x[1]).toBe(2);
    });

    it('mixing coefficient depths is rejected by the comparisons', () => {
        const shallow = new QFNumber(1, 2, 3);
        const deep = new QFNumber(new QFNumber(1, 0, 2), new QFNumber(0, 0, 2), 3);
        expect(() => shallow.equals(deep)).toThrow(/Mismatched coefficient depth/);
        expect(() => shallow.lessThan(deep)).toThrow(/Mismatched coefficient depth/);
    });
});
