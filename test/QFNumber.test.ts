import { describe, it, expect } from 'vitest';
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
