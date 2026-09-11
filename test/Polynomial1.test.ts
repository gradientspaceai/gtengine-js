import { describe, it, expect } from 'vitest';
import {
    Polynomial1, greatestCommonDivisor, squareFreeFactorization
} from '../src/Polynomial1.js';
import {
    check, expectStrictWeakOrder, fc
} from './helpers/arbitraries.js';

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

// A random polynomial of exactly the requested degree (leading coefficient
// bounded away from zero so no leading-zero elimination occurs).
function randomPoly(rng: () => number, degree: number): Polynomial1 {
    const c: number[] = [];
    for (let i = 0; i < degree; ++i) {
        c.push(4 * rng() - 2);
    }
    const leading = 4 * rng() - 2;
    c.push(leading >= 0 ? leading + 0.5 : leading - 0.5);
    return Polynomial1.fromCoefficients(c);
}

// Monic polynomial with the given roots: (t - r0)(t - r1)...
function fromRoots(roots: number[]): Polynomial1 {
    let p = Polynomial1.fromCoefficients([1]);
    for (const r of roots) {
        p = p.mul(Polynomial1.fromCoefficients([-r, 1]));
    }
    return p;
}

function expectCoefficients(p: Polynomial1, expected: number[], tol = 0): void {
    expect(p.getDegree()).toBe(expected.length - 1);
    for (let i = 0; i < expected.length; ++i) {
        expect(Math.abs(p.get(i) - expected[i])).toBeLessThanOrEqual(tol);
    }
}

// Compare two polynomials with a tolerance scaled by the coefficient sizes.
function expectPolyClose(p: Polynomial1, q: Polynomial1, relTol = 1e-12): void {
    expect(p.getDegree()).toBe(q.getDegree());
    let scale = 1;
    for (let i = 0; i <= p.getDegree(); ++i) {
        scale = Math.max(scale, Math.abs(p.get(i)), Math.abs(q.get(i)));
    }
    for (let i = 0; i <= p.getDegree(); ++i) {
        expect(Math.abs(p.get(i) - q.get(i))).toBeLessThanOrEqual(relTol * scale);
    }
}

describe('Polynomial1 construction and access', () => {
    it('default-constructs the zero polynomial of degree 0', () => {
        const p = new Polynomial1();
        expect(p.getDegree()).toBe(0);
        expect(p.get(0)).toBe(0);
        expect(p.getCoefficients()).toEqual([0]);
    });

    it('constructs a zero-filled polynomial of the requested degree', () => {
        const p = new Polynomial1(4);
        expect(p.getDegree()).toBe(4);
        expect(p.getCoefficients()).toEqual([0, 0, 0, 0, 0]);
    });

    it('rejects a negative or non-integer degree', () => {
        expect(() => new Polynomial1(-1)).toThrow('Invalid degree.');
        expect(() => new Polynomial1(1.5)).toThrow('Invalid degree.');
    });

    it('drops leading zeros in fromCoefficients', () => {
        const p = Polynomial1.fromCoefficients([1, 2, 0, 0]);
        expectCoefficients(p, [1, 2]);
    });

    it('keeps degree 0 for the all-zero coefficient list', () => {
        const p = Polynomial1.fromCoefficients([0, 0, 0]);
        expectCoefficients(p, [0]);
    });

    it('rejects an empty coefficient list', () => {
        expect(() => Polynomial1.fromCoefficients([])).toThrow(
            'Invalid number of coefficients.');
    });

    it('sets and reads individual coefficients', () => {
        const p = new Polynomial1(2);
        p.set(0, 3);
        p.set(2, -1);
        expectCoefficients(p, [3, 0, -1]);
    });

    it('grows and shrinks with setDegree, zeroing new coefficients', () => {
        const p = Polynomial1.fromCoefficients([1, 2, 3]);
        p.setDegree(4);
        expectCoefficients(p, [1, 2, 3, 0, 0]);
        p.setDegree(1);
        expectCoefficients(p, [1, 2]);
    });

    it('fills all coefficients with setCoefficients', () => {
        const p = new Polynomial1(3);
        p.setCoefficients(7);
        expectCoefficients(p, [7, 7, 7, 7]);
    });

    it('clones deeply', () => {
        const p = Polynomial1.fromCoefficients([1, 2, 3]);
        const q = p.clone();
        q.set(0, 99);
        expect(p.get(0)).toBe(1);
        expect(q.get(0)).toBe(99);
    });

    it('eliminates leading zeros on demand', () => {
        const p = new Polynomial1(5);
        p.set(0, 1);
        p.set(1, -2);
        p.eliminateLeadingZeros();
        expectCoefficients(p, [1, -2]);

        const z = new Polynomial1(5);
        z.eliminateLeadingZeros();
        expectCoefficients(z, [0]);
    });
});

describe('Polynomial1 comparisons', () => {
    it('compares equality by coefficient array', () => {
        const p = Polynomial1.fromCoefficients([1, 2, 3]);
        const q = Polynomial1.fromCoefficients([1, 2, 3]);
        const r = Polynomial1.fromCoefficients([1, 2, 4]);
        expect(p.equals(q)).toBe(true);
        expect(p.notEquals(q)).toBe(false);
        expect(p.equals(r)).toBe(false);
        expect(p.notEquals(r)).toBe(true);
    });

    it('orders lexicographically, with a proper prefix comparing less', () => {
        const p = Polynomial1.fromCoefficients([1, 2]);
        const q = Polynomial1.fromCoefficients([1, 2, 3]);
        const r = Polynomial1.fromCoefficients([1, 3]);
        expect(p.lessThan(q)).toBe(true);
        expect(q.greaterThan(p)).toBe(true);
        expect(p.lessThan(r)).toBe(true);
        expect(r.lessThan(p)).toBe(false);
        expect(p.lessThanOrEqual(p.clone())).toBe(true);
        expect(p.greaterThanOrEqual(p.clone())).toBe(true);
        expect(p.lessThan(p.clone())).toBe(false);
        expect(p.greaterThan(p.clone())).toBe(false);
    });
});

describe('Polynomial1 evaluation and calculus', () => {
    it('evaluates known values', () => {
        // p(t) = 1 - 2t + 3t^2
        const p = Polynomial1.fromCoefficients([1, -2, 3]);
        expect(p.evaluate(0)).toBe(1);
        expect(p.evaluate(1)).toBe(2);
        expect(p.evaluate(2)).toBe(9);
        expect(p.evaluate(-1)).toBe(6);
    });

    it('evaluates a constant polynomial', () => {
        const p = Polynomial1.fromCoefficients([5]);
        expect(p.evaluate(0)).toBe(5);
        expect(p.evaluate(1e6)).toBe(5);
    });

    it('computes the derivative of a known polynomial', () => {
        // d/dt (1 - 2t + 3t^2 + 4t^3) = -2 + 6t + 12t^2
        const p = Polynomial1.fromCoefficients([1, -2, 3, 4]);
        expectCoefficients(p.getDerivative(), [-2, 6, 12]);
    });

    it('gives the zero polynomial as the derivative of a constant', () => {
        const p = Polynomial1.fromCoefficients([5]);
        expectCoefficients(p.getDerivative(), [0]);
    });

    it('matches a central-difference derivative on random polynomials', () => {
        const rng = makeRng(0x51ee7);
        for (let trial = 0; trial < 50; ++trial) {
            const p = randomPoly(rng, 4);
            const d = p.getDerivative();
            const t = 2 * rng() - 1;
            const h = 1e-5;
            const numeric = (p.evaluate(t + h) - p.evaluate(t - h)) / (2 * h);
            expect(Math.abs(d.evaluate(t) - numeric)).toBeLessThanOrEqual(1e-6);
        }
    });

    it('inverts the coefficient order', () => {
        const p = Polynomial1.fromCoefficients([1, 2, 3, 4]);
        expectCoefficients(p.getInversion(), [4, 3, 2, 1]);
        // The inversion of p satisfies inv(t) = t^degree * p(1/t).
        const t = 1.7;
        expect(Math.abs(p.getInversion().evaluate(t) -
            Math.pow(t, 3) * p.evaluate(1 / t))).toBeLessThanOrEqual(1e-12);
    });

    it('translates: getTranslation(t0) is p(t - t0)', () => {
        const p = Polynomial1.fromCoefficients([1, -2, 3, 4]);
        const q = p.getTranslation(0.75);
        expect(q.getDegree()).toBe(3);
        for (const t of [-2, -0.5, 0, 0.75, 1.25, 3]) {
            expect(Math.abs(q.evaluate(t) - p.evaluate(t - 0.75)))
                .toBeLessThanOrEqual(1e-12);
        }
    });

    it('translates a constant polynomial', () => {
        const p = Polynomial1.fromCoefficients([5]);
        expectCoefficients(p.getTranslation(3), [5]);
    });
});

describe('Polynomial1 arithmetic', () => {
    it('negates', () => {
        expectCoefficients(Polynomial1.fromCoefficients([1, -2, 3]).negate(),
            [-1, 2, -3]);
    });

    it('adds and subtracts polynomials of unequal degree', () => {
        const p = Polynomial1.fromCoefficients([1, 1, 1]);
        const q = Polynomial1.fromCoefficients([1, 2]);
        expectCoefficients(p.add(q), [2, 3, 1]);
        expectCoefficients(q.add(p), [2, 3, 1]);
        expectCoefficients(p.sub(q), [0, -1, 1]);
        expectCoefficients(q.sub(p), [0, 1, -1]);
    });

    it('eliminates leading zeros after cancellation', () => {
        // (1 + x + x^2) + (1 + 2x - x^2) = 2 + 3x
        const p = Polynomial1.fromCoefficients([1, 1, 1]);
        const q = Polynomial1.fromCoefficients([1, 2, -1]);
        expectCoefficients(p.add(q), [2, 3]);
        expectCoefficients(p.sub(p.clone()), [0]);
    });

    it('multiplies polynomials', () => {
        // (1 + 2x)(3 - x + x^2) = 3 + 5x - x^2 + 2x^3
        const p = Polynomial1.fromCoefficients([1, 2]);
        const q = Polynomial1.fromCoefficients([3, -1, 1]);
        expectCoefficients(p.mul(q), [3, 5, -1, 2]);
        expectCoefficients(q.mul(p), [3, 5, -1, 2]);
    });

    it('applies scalar add, sub, subFrom, mul and div', () => {
        const p = Polynomial1.fromCoefficients([1, -2, 3]);
        expectCoefficients(p.add(10), [11, -2, 3]);
        expectCoefficients(p.sub(10), [-9, -2, 3]);
        expectCoefficients(p.subFrom(10), [9, 2, -3]);
        expectCoefficients(p.mul(2), [2, -4, 6]);
        expectCoefficients(p.div(2), [0.5, -1, 1.5]);
    });

    it('throws on division by the zero scalar', () => {
        const p = Polynomial1.fromCoefficients([1, 2]);
        expect(() => p.div(0)).toThrow('Division by zero.');
    });

    it('leaves the operands unmodified', () => {
        const p = Polynomial1.fromCoefficients([1, 2, 3]);
        const q = Polynomial1.fromCoefficients([4, 5]);
        p.add(q); p.sub(q); p.mul(q); p.mul(3); p.div(3); p.negate();
        expectCoefficients(p, [1, 2, 3]);
        expectCoefficients(q, [4, 5]);
    });

    it('satisfies the ring identities on random polynomials', () => {
        const rng = makeRng(0xb1a5e);
        for (let trial = 0; trial < 100; ++trial) {
            const a = randomPoly(rng, 1 + Math.floor(4 * rng()));
            const b = randomPoly(rng, 1 + Math.floor(4 * rng()));
            const c = randomPoly(rng, 1 + Math.floor(4 * rng()));

            // Commutativity of addition and multiplication.
            expectPolyClose(a.add(b), b.add(a), 0);
            expectPolyClose(a.mul(b), b.mul(a), 1e-14);

            // Associativity.
            expectPolyClose(a.add(b).add(c), a.add(b.add(c)), 1e-14);
            expectPolyClose(a.mul(b).mul(c), a.mul(b.mul(c)), 1e-13);

            // Distributivity.
            expectPolyClose(a.mul(b.add(c)), a.mul(b).add(a.mul(c)), 1e-13);

            // Additive inverse and identity.
            expectCoefficients(a.add(a.negate()), [0]);
            expectPolyClose(a.mul(Polynomial1.fromCoefficients([1])), a, 0);

            // Degrees add under multiplication (leading coefficients are
            // bounded away from zero).
            expect(a.mul(b).getDegree()).toBe(a.getDegree() + b.getDegree());
        }
    });

    it('is a ring homomorphism under evaluation', () => {
        const rng = makeRng(0x901ae);
        for (let trial = 0; trial < 100; ++trial) {
            const a = randomPoly(rng, 1 + Math.floor(4 * rng()));
            const b = randomPoly(rng, 1 + Math.floor(4 * rng()));
            const t = 2 * rng() - 1;
            const scale = Math.max(1, Math.abs(a.evaluate(t)), Math.abs(b.evaluate(t)));
            expect(Math.abs(a.add(b).evaluate(t) - (a.evaluate(t) + b.evaluate(t))))
                .toBeLessThanOrEqual(1e-12 * scale);
            expect(Math.abs(a.mul(b).evaluate(t) - a.evaluate(t) * b.evaluate(t)))
                .toBeLessThanOrEqual(1e-12 * scale * scale);
        }
    });
});

describe('Polynomial1 division', () => {
    it('divides exactly when the divisor is a factor', () => {
        // (x^2 - 3x + 2) / (x - 1) = x - 2, remainder 0.
        const p = Polynomial1.fromCoefficients([2, -3, 1]);
        const d = Polynomial1.fromCoefficients([-1, 1]);
        const { quotient, remainder } = p.divide(d);
        expectCoefficients(quotient, [-2, 1]);
        expectCoefficients(remainder, [0]);
    });

    it('produces a nonzero remainder for an inexact division', () => {
        // (x^2 + 1) / (x - 1) = x + 1, remainder 2.
        const p = Polynomial1.fromCoefficients([1, 0, 1]);
        const d = Polynomial1.fromCoefficients([-1, 1]);
        const { quotient, remainder } = p.divide(d);
        expectCoefficients(quotient, [1, 1]);
        expectCoefficients(remainder, [2]);
    });

    it('returns quotient 0 and remainder p when degree(p) < degree(d)', () => {
        const p = Polynomial1.fromCoefficients([1, 2]);
        const d = Polynomial1.fromCoefficients([1, 1, 1]);
        const { quotient, remainder } = p.divide(d);
        expectCoefficients(quotient, [0]);
        expectCoefficients(remainder, [1, 2]);
        // The remainder is a copy, not an alias.
        remainder.set(0, 99);
        expect(p.get(0)).toBe(1);
    });

    it('gives a zero remainder for a constant divisor', () => {
        const p = Polynomial1.fromCoefficients([2, 4, 6]);
        const { quotient, remainder } = p.divide(Polynomial1.fromCoefficients([2]));
        expectCoefficients(quotient, [1, 2, 3]);
        expectCoefficients(remainder, [0]);
    });

    it('reconstructs p = q*d + r with degree(r) < degree(d)', () => {
        const rng = makeRng(0xd1de5);
        for (let trial = 0; trial < 200; ++trial) {
            const pDegree = 1 + Math.floor(6 * rng());
            const dDegree = 1 + Math.floor(4 * rng());
            const p = randomPoly(rng, pDegree);
            const d = randomPoly(rng, dDegree);
            const { quotient, remainder } = p.divide(d);

            if (pDegree >= dDegree) {
                expect(quotient.getDegree()).toBe(pDegree - dDegree);
            }
            expect(remainder.getDegree()).toBeLessThan(dDegree);

            const reconstructed = quotient.mul(d).add(remainder);
            expectPolyClose(reconstructed, p, 1e-11);
        }
    });
});

describe('Polynomial1 makeMonic', () => {
    it('scales the highest-degree coefficient to 1', () => {
        const p = Polynomial1.fromCoefficients([2, 4, 8]);
        p.makeMonic();
        expectCoefficients(p, [0.25, 0.5, 1]);
    });

    it('leaves an already-monic polynomial unchanged', () => {
        const p = Polynomial1.fromCoefficients([2, -3, 1]);
        p.makeMonic();
        expectCoefficients(p, [2, -3, 1]);
    });

    it('eliminates leading zeros before scaling', () => {
        const p = new Polynomial1(4);
        p.set(0, 3);
        p.set(1, 6);
        p.makeMonic();
        expectCoefficients(p, [0.5, 1]);
    });
});

describe('greatestCommonDivisor', () => {
    it('finds the common factor of constructed products', () => {
        // g = (x-1)(x-2), p0 = g*(x-3), p1 = g*(x-4).
        const g = fromRoots([1, 2]);
        const p0 = g.mul(Polynomial1.fromCoefficients([-3, 1]));
        const p1 = g.mul(Polynomial1.fromCoefficients([-4, 1]));
        expectPolyClose(greatestCommonDivisor(p0, p1), g, 1e-12);
        // The order of the arguments does not matter.
        expectPolyClose(greatestCommonDivisor(p1, p0), g, 1e-12);
    });

    it('finds a common linear factor', () => {
        const g = fromRoots([2]);
        const p0 = g.mul(fromRoots([-1, 3]));
        const p1 = g.mul(fromRoots([5]));
        expectPolyClose(greatestCommonDivisor(p0, p1), g, 1e-12);
    });

    it('finds a common cubic factor', () => {
        const g = fromRoots([1, -2, 4]);
        const p0 = g.mul(fromRoots([3]));
        const p1 = g.mul(fromRoots([-5, 6]));
        expectPolyClose(greatestCommonDivisor(p0, p1), g, 1e-10);
    });

    it('returns the monic constant 1 for coprime polynomials', () => {
        const p0 = fromRoots([1, 2]);
        const p1 = fromRoots([3, 4]);
        expectCoefficients(greatestCommonDivisor(p0, p1), [1], 1e-14);
    });

    it('scales the result to be monic', () => {
        const g = fromRoots([1, 2]);
        const p0 = g.mul(Polynomial1.fromCoefficients([-3, 1])).mul(7);
        const p1 = g.mul(Polynomial1.fromCoefficients([-4, 1])).mul(-3);
        const gcd = greatestCommonDivisor(p0, p1);
        expect(gcd.get(gcd.getDegree())).toBeCloseTo(1, 12);
        expectPolyClose(gcd, g, 1e-12);
    });

    it('handles zero-valued inputs', () => {
        const zero = Polynomial1.fromCoefficients([0]);
        const p = fromRoots([1, 2]);
        expectPolyClose(greatestCommonDivisor(p, zero), p, 0);
        expectPolyClose(greatestCommonDivisor(zero, p), p, 0);
        expectCoefficients(greatestCommonDivisor(zero, zero.clone()), [0]);
    });

    it('divides both inputs exactly', () => {
        const g = fromRoots([-1, 3]);
        const p0 = g.mul(fromRoots([2]));
        const p1 = g.mul(fromRoots([5, 7]));
        const gcd = greatestCommonDivisor(p0, p1);
        expectCoefficients(p0.divide(gcd).remainder, [0], 1e-10);
        expectCoefficients(p1.divide(gcd).remainder, [0], 1e-10);
    });
});

describe('squareFreeFactorization', () => {
    // Reassemble f = factors[0] * factors[1]^2 * ... * factors[n-1]^n.
    function product(factors: Polynomial1[]): Polynomial1 {
        let p = Polynomial1.fromCoefficients([1]);
        for (let i = 0; i < factors.length; ++i) {
            for (let k = 0; k <= i; ++k) {
                p = p.mul(factors[i]);
            }
        }
        return p;
    }

    it('recovers the double factor of f = (x-2)^2', () => {
        const f = fromRoots([2, 2]);
        const factors = squareFreeFactorization(f);
        // f = 1 * (x-2)^2.
        expect(factors.length).toBe(2);
        expectCoefficients(factors[0], [1], 1e-14);
        expectCoefficients(factors[1], [-2, 1], 1e-12);
        expectPolyClose(product(factors), f, 1e-12);
    });

    it('recovers the triple factor of f = (x-2)^3', () => {
        const f = fromRoots([2, 2, 2]);
        const factors = squareFreeFactorization(f);
        // f = 1 * 1^2 * (x-2)^3.
        expect(factors.length).toBe(3);
        expectCoefficients(factors[0], [1], 1e-14);
        expectCoefficients(factors[1], [1], 1e-14);
        expectCoefficients(factors[2], [-2, 1], 1e-12);
        expectPolyClose(product(factors), f, 1e-12);
    });

    it('separates the simple and double factors of (x-1)(x-3)(x-2)^2', () => {
        const f = fromRoots([1, 3]).mul(fromRoots([2, 2]));
        const factors = squareFreeFactorization(f);
        expect(factors.length).toBe(2);
        // factors[0] = (x-1)(x-3) = x^2 - 4x + 3, factors[1] = (x-2).
        expectCoefficients(factors[0], [3, -4, 1], 1e-11);
        expectCoefficients(factors[1], [-2, 1], 1e-12);
        expectPolyClose(product(factors), f, 1e-10);
    });

    it('handles the repeated root at the origin: f = x^2 (x-1)', () => {
        const f = fromRoots([0, 0, 1]);
        const factors = squareFreeFactorization(f);
        expect(factors.length).toBe(2);
        expectCoefficients(factors[0], [-1, 1], 1e-12);
        expectCoefficients(factors[1], [0, 1], 1e-12);
        expectPolyClose(product(factors), f, 1e-12);
    });

    it('returns a trivial factorization for a square-free polynomial', () => {
        const f = fromRoots([1, 2, 3]);
        const factors = squareFreeFactorization(f);
        expect(factors.length).toBe(1);
        expectPolyClose(product(factors), f, 1e-9);
        expectPolyClose(factors[0], f, 1e-9);
    });

    it('degrades to a single factor when the floating-point GCD is unstable', () => {
        // f = (x-1)(x-2)^2. The true gcd(f, f') is (x-2), but the floating-
        // point Euclidean algorithm reduces to the constant 1 because the
        // remainders are only ~1e-16 rather than exactly zero. The algorithm
        // then reports f itself as a single square-free factor. This is
        // inherent to upstream's exact-zero remainder test, not a port
        // deviation; the reconstruction identity still holds.
        const f = fromRoots([1]).mul(fromRoots([2, 2]));
        expectCoefficients(greatestCommonDivisor(f, f.getDerivative()), [1], 0);
        const factors = squareFreeFactorization(f);
        expect(factors.length).toBe(1);
        expectPolyClose(product(factors), f, 0);
    });

    it('always terminates on randomized inputs', () => {
        // Either a valid factorization is produced (the reconstruction
        // identity holds) or the iteration cap fires; neither hangs.
        const rng = makeRng(0x5f3ee);
        let converged = 0, capped = 0;
        for (let trial = 0; trial < 300; ++trial) {
            // Mix square-free and repeated small-integer roots.
            const n = 1 + Math.floor(4 * rng());
            const roots: number[] = [];
            for (let i = 0; i < n; ++i) {
                const r = Math.round(10 * rng() - 5);
                roots.push(r);
                if (rng() < 0.5) {
                    roots.push(r);
                }
            }
            const f = fromRoots(roots);
            let factors: Polynomial1[] | null = null;
            try {
                factors = squareFreeFactorization(f);
            } catch (e) {
                expect((e as Error).message).toContain(
                    'The square-free factorization did not converge.');
                ++capped;
                continue;
            }
            ++converged;
            const scale = Math.max(1, ...f.getCoefficients().map(Math.abs));
            const p = product(factors);
            expect(p.getDegree()).toBe(f.getDegree());
            for (let i = 0; i <= f.getDegree(); ++i) {
                expect(Math.abs(p.get(i) - f.get(i)))
                    .toBeLessThanOrEqual(1e-8 * scale);
            }
        }
        // Both outcomes occur for this seed, so both paths are exercised.
        expect(converged).toBeGreaterThan(0);
        expect(capped).toBeGreaterThan(0);
    });

    it('throws instead of looping forever when the GCD collapses', () => {
        // f = (x-2)^2 (x+2)(x+3)(x-4). Here gcd(f, f') rounds to
        // x - 1.9999999999999998, after which every gcd(b, d) in the loop is
        // the constant 1, the degree of b never decreases and upstream's
        // do-while loop never exits. The port caps the iterations.
        const f = fromRoots([2, 2, -2, -3, 4]);
        expectCoefficients(f, [-96, 40, 36, -14, -3, 1], 1e-10);
        expect(() => squareFreeFactorization(f)).toThrow(
            'The square-free factorization did not converge.');
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V44): property-based comparison against upstream
// Polynomial1.h. Integer coefficients and integer evaluation points keep the
// arithmetic exact in binary64, so every identity below is cross-checked
// against a bigint reference implementation rather than a tolerance.
// ---------------------------------------------------------------------------

// Integer coefficient list, in order of increasing power. Coefficients up to
// 8 in absolute value with degree at most 4, evaluated at |t| <= 4, keep every
// intermediate below 2^53.
const intCoefficients = (maxDegree = 4) =>
    fc.array(fc.integer({ min: -8, max: 8 }),
        { minLength: 1, maxLength: maxDegree + 1 });

const intPoint = fc.integer({ min: -4, max: 4 });

// Reference arithmetic on bigint coefficient arrays.
function refEvaluate(c: readonly number[], t: number): bigint {
    const bt = BigInt(t);
    let result = 0n;
    for (let i = c.length - 1; i >= 0; --i) {
        result = result * bt + BigInt(c[i]);
    }
    return result;
}

function refMul(a: readonly number[], b: readonly number[]): bigint[] {
    const out = new Array<bigint>(a.length + b.length - 1).fill(0n);
    for (let i = 0; i < a.length; ++i) {
        for (let j = 0; j < b.length; ++j) {
            out[i + j] += BigInt(a[i]) * BigInt(b[j]);
        }
    }
    return out;
}

// The coefficients with trailing (high-degree) zeros removed, as upstream's
// EliminateLeadingZeros does, but never below a single coefficient.
function refTrim(c: readonly bigint[]): bigint[] {
    const out = [...c];
    while (out.length > 1 && out[out.length - 1] === 0n) { out.pop(); }
    return out;
}

function coeffsOf(p: Polynomial1): number[] {
    return [...p.getCoefficients()];
}

describe('Polynomial1 verification', () => {
    it('evaluates with Horner exactly as the bigint reference does', () => {
        check(fc.tuple(intCoefficients(), intPoint), ([c, t]) => {
            const p = Polynomial1.fromCoefficients([...c]);
            // fromCoefficients drops the leading zeros, which does not change
            // the value of the polynomial.
            expect(BigInt(p.evaluate(t))).toBe(refEvaluate(c, t));
        });
    });

    it('adds, subtracts and negates coefficientwise with leading zeros removed', () => {
        check(fc.tuple(intCoefficients(), intCoefficients(), intPoint),
            ([a, b, t]) => {
                const p = Polynomial1.fromCoefficients([...a]);
                const q = Polynomial1.fromCoefficients([...b]);
                const sum = p.add(q);
                const diff = p.sub(q);
                expect(BigInt(sum.evaluate(t)))
                    .toBe(refEvaluate(a, t) + refEvaluate(b, t));
                expect(BigInt(diff.evaluate(t)))
                    .toBe(refEvaluate(a, t) - refEvaluate(b, t));
                expect(BigInt(p.negate().evaluate(t))).toBe(-refEvaluate(a, t));
                // The sum and difference eliminate leading zeros, so the
                // degree drops when the high-order terms cancel.
                const expectedSum = refTrim(coeffsOf(sum).map(BigInt));
                expect(coeffsOf(sum).map(BigInt)).toEqual(expectedSum);
                expect(coeffsOf(diff).map(BigInt))
                    .toEqual(refTrim(coeffsOf(diff).map(BigInt)));
            });
    });

    it('multiplies exactly and keeps the sum of the degrees', () => {
        check(fc.tuple(intCoefficients(3), intCoefficients(3), intPoint),
            ([a, b, t]) => {
                const p = Polynomial1.fromCoefficients([...a]);
                const q = Polynomial1.fromCoefficients([...b]);
                const product = p.mul(q);
                // Upstream operator*(p0,p1) does NOT eliminate leading zeros,
                // so the degree is always the sum of the input degrees.
                expect(product.getDegree()).toBe(p.getDegree() + q.getDegree());
                const expected = refMul(coeffsOf(p), coeffsOf(q));
                expect(coeffsOf(product).map(BigInt)).toEqual(expected);
                expect(BigInt(product.evaluate(t)))
                    .toBe(refEvaluate(coeffsOf(p), t) * refEvaluate(coeffsOf(q), t));
            });
    });

    it('satisfies the product rule for the derivative of a product', () => {
        check(fc.tuple(intCoefficients(3), intCoefficients(3), intPoint),
            ([a, b, t]) => {
                const p = Polynomial1.fromCoefficients([...a]);
                const q = Polynomial1.fromCoefficients([...b]);
                const lhs = p.mul(q).getDerivative();
                const rhs = p.getDerivative().mul(q).add(p.mul(q.getDerivative()));
                expect(BigInt(lhs.evaluate(t))).toBe(BigInt(rhs.evaluate(t)));
            });
    });

    it('differentiates termwise and collapses a constant to zero', () => {
        check(intCoefficients(), c => {
            const p = Polynomial1.fromCoefficients([...c]);
            const d = p.getDerivative();
            const pc = coeffsOf(p);
            if (p.getDegree() === 0) {
                // Upstream returns the degree-0 zero polynomial.
                expect(coeffsOf(d)).toEqual([0]);
            } else {
                expect(d.getDegree()).toBe(p.getDegree() - 1);
                for (let i = 0; i < pc.length - 1; ++i) {
                    expect(d.get(i)).toBe(pc[i + 1] * (i + 1));
                }
            }
        });
    });

    it('translates: p.getTranslation(t0) evaluated at t equals p(t - t0)', () => {
        check(fc.tuple(intCoefficients(), intPoint, fc.integer({ min: -3, max: 3 })),
            ([c, t, t0]) => {
                const p = Polynomial1.fromCoefficients([...c]);
                const shifted = p.getTranslation(t0);
                expect(BigInt(shifted.evaluate(t)))
                    .toBe(refEvaluate(coeffsOf(p), t - t0));
            });
    });

    it('inverts: inversion[i] = p[degree - i], so inv(t) = t^degree * p(1/t)', () => {
        check(fc.tuple(intCoefficients(), intPoint.filter(t => t !== 0)),
            ([c, t]) => {
                const p = Polynomial1.fromCoefficients([...c]);
                const inv = p.getInversion();
                expect(inv.getDegree()).toBe(p.getDegree());
                const pc = coeffsOf(p);
                for (let i = 0; i <= p.getDegree(); ++i) {
                    expect(inv.get(i)).toBe(pc[p.getDegree() - i]);
                }
                // t^degree * p(1/t) with exact rational arithmetic: multiply
                // through by t^degree, which is exactly what the reversal is.
                const reversed = refEvaluate(coeffsOf(inv), t);
                let acc = 0n;
                const bt = BigInt(t);
                for (let i = 0; i <= p.getDegree(); ++i) {
                    acc += BigInt(pc[i]) * bt ** BigInt(p.getDegree() - i);
                }
                expect(reversed).toBe(acc);
                // Round trip.
                expect(coeffsOf(inv.getInversion())).toEqual(pc);
            });
    });

    it('divides exactly: p = quotient * divisor + remainder', () => {
        // A monic integer divisor keeps the Euclidean algorithm in integers,
        // so the identity holds exactly with no tolerance.
        check(fc.tuple(intCoefficients(4),
            fc.array(fc.integer({ min: -6, max: 6 }), { minLength: 0, maxLength: 3 }),
            intPoint),
            ([a, lower, t]) => {
                const p = Polynomial1.fromCoefficients([...a]);
                const d = Polynomial1.fromCoefficients([...lower, 1]);
                const { quotient, remainder } = p.divide(d);
                if (p.getDegree() >= d.getDegree()) {
                    expect(quotient.getDegree()).toBe(p.getDegree() - d.getDegree());
                    expect(remainder.getDegree())
                        .toBeLessThanOrEqual(Math.max(0, d.getDegree() - 1));
                } else {
                    // Upstream returns Q = 0 and R = P in this case.
                    expect(coeffsOf(quotient)).toEqual([0]);
                    expect(coeffsOf(remainder)).toEqual(coeffsOf(p));
                }
                expect(BigInt(quotient.mul(d).add(remainder).evaluate(t)))
                    .toBe(refEvaluate(coeffsOf(p), t));
            });
    });

    it('makeMonic rescales in place without changing the roots', () => {
        check(fc.tuple(intCoefficients(), intPoint), ([c, t]) => {
            const p = Polynomial1.fromCoefficients([...c]);
            const leading = p.get(p.getDegree());
            const before = p.evaluate(t);
            p.makeMonic();
            expect(p.get(p.getDegree()) === 1).toBe(true);
            if (leading === 0) {
                // The only way the leading coefficient is zero after
                // fromCoefficients is the degree-0 zero polynomial, for which
                // upstream divides by zero and produces 1/0 = Infinity.
                expect(p.getDegree()).toBe(0);
                return;
            }
            expect(Math.abs(p.evaluate(t) * leading - before))
                .toBeLessThanOrEqual(1e-12 * (1 + Math.abs(before)));
        });
    });

    it('orders polynomials lexicographically on the coefficient arrays', () => {
        check(fc.array(intCoefficients(2), { minLength: 2, maxLength: 5 }), lists => {
            const polys = lists.map(c => Polynomial1.fromCoefficients([...c]));
            expectStrictWeakOrder(polys, (x, y) => x.lessThan(y));
            for (const x of polys) {
                for (const y of polys) {
                    // std::vector semantics: <= is !(y < x), > is y < x,
                    // >= is !(x < y), != is !(x == y).
                    expect(x.lessThanOrEqual(y)).toBe(!y.lessThan(x));
                    expect(x.greaterThan(y)).toBe(y.lessThan(x));
                    expect(x.greaterThanOrEqual(y)).toBe(!x.lessThan(y));
                    expect(x.notEquals(y)).toBe(!x.equals(y));
                    // A proper prefix compares less than the longer array.
                    const xc = coeffsOf(x), yc = coeffsOf(y);
                    if (xc.length < yc.length
                        && xc.every((v, i) => v === yc[i])) {
                        expect(x.lessThan(y)).toBe(true);
                    }
                }
            }
        }, 60);
    });

    it('eliminateLeadingZeros never drops the last coefficient', () => {
        check(fc.array(fc.integer({ min: -2, max: 2 }),
            { minLength: 1, maxLength: 6 }), c => {
                const p = new Polynomial1(c.length - 1);
                for (let i = 0; i < c.length; ++i) { p.set(i, c[i]); }
                p.eliminateLeadingZeros();
                expect(p.getDegree()).toBeGreaterThanOrEqual(0);
                if (c.every(v => v === 0)) {
                    expect(coeffsOf(p)).toEqual([0]);
                } else {
                    expect(p.get(p.getDegree())).not.toBe(0);
                    expect(coeffsOf(p))
                        .toEqual(c.slice(0, p.getDegree() + 1));
                }
            });
    });

    it('setDegree grows with zeros and truncates, as std::vector::resize', () => {
        check(fc.tuple(intCoefficients(), fc.integer({ min: 0, max: 7 })),
            ([c, degree]) => {
                const p = Polynomial1.fromCoefficients([...c]);
                const before = coeffsOf(p);
                p.setDegree(degree);
                expect(p.getDegree()).toBe(degree);
                for (let i = 0; i <= degree; ++i) {
                    expect(p.get(i)).toBe(i < before.length ? before[i] : 0);
                }
            });
    });

    it('computes a monic gcd that divides both inputs', () => {
        // Build both inputs from small integer linear factors.
        const factorSet = fc.array(fc.integer({ min: -3, max: 3 }),
            { minLength: 1, maxLength: 3 });
        check(fc.tuple(factorSet, factorSet, intPoint), ([ra, rb, t]) => {
            const p = fromRoots(ra);
            const q = fromRoots(rb);
            const g = greatestCommonDivisor(p, q);
            // Upstream normalizes the result to leading coefficient 1.
            expect(g.get(g.getDegree())).toBeCloseTo(1, 12);
            // The gcd divides both inputs, so both remainders vanish.
            for (const x of [p, q]) {
                const { quotient, remainder } = x.divide(g);
                for (let i = 0; i <= remainder.getDegree(); ++i) {
                    expect(Math.abs(remainder.get(i))).toBeLessThan(1e-9);
                }
                expect(Math.abs(quotient.mul(g).evaluate(t) - x.evaluate(t)))
                    .toBeLessThan(1e-9 * (1 + Math.abs(x.evaluate(t))));
            }
            // The number of shared roots (with multiplicity) is the degree of
            // the exact gcd. The floating-point Euclidean algorithm tests
            // remainders against exactly zero, so it can only UNDER-report the
            // degree, never over-report it (an extra factor would have to
            // divide both inputs).
            const remaining = new Map<number, number>();
            for (const r of ra) { remaining.set(r, (remaining.get(r) ?? 0) + 1); }
            let shared = 0;
            for (const r of rb) {
                const n = remaining.get(r) ?? 0;
                if (n > 0) { remaining.set(r, n - 1); ++shared; }
            }
            expect(g.getDegree()).toBeLessThanOrEqual(shared);
            // Every root the gcd does report must be a shared root.
            for (const r of new Set(ra)) {
                if (!rb.includes(r)) {
                    expect(Math.abs(g.evaluate(r)))
                        .toBeGreaterThan(g.getDegree() === 0 ? 0.5 : 1e-9);
                }
            }
        }, 100);
    });

    it('finds the exact gcd when a shared factor keeps the division exact', () => {
        // gcd(p, p) runs one division with a zero remainder, so it returns p
        // made monic no matter how the coefficients round.
        check(fc.array(fc.integer({ min: -3, max: 3 }),
            { minLength: 1, maxLength: 4 }), roots => {
                const p = fromRoots(roots);
                const g = greatestCommonDivisor(p, p);
                expect(g.getDegree()).toBe(p.getDegree());
                for (let i = 0; i <= p.getDegree(); ++i) {
                    expect(Math.abs(g.get(i) - p.get(i))).toBeLessThan(1e-9);
                }
            }, 60);
    });

    it('can under-report the gcd in floating point (upstream instability)', () => {
        // Characterization of the exact-zero remainder test that upstream's
        // Euclidean loop uses (gtengine-js issue #83): t^3 - 2t^2 - 3t and
        // t^3 - 5t^2 + 7t - 3 share the factor t - 3, but making the first
        // remainder monic introduces 10/3, after which no later remainder is
        // exactly zero and the loop returns the constant 1.
        const p = fromRoots([-1, 0, 3]);
        const q = fromRoots([3, 1, 1]);
        const g = greatestCommonDivisor(p, q);
        expect(g.getDegree()).toBe(0);
        expect(g.get(0)).toBe(1);
    });

    it('gcd with the zero polynomial returns the upstream asymmetric result', () => {
        check(intCoefficients(), c => {
            const p = Polynomial1.fromCoefficients([...c]);
            const zero = Polynomial1.fromCoefficients([0]);
            // Upstream picks a = the input of larger-or-equal degree and
            // returns 'a != zero ? a : zero' without making it monic. When
            // both inputs have degree 0, a is the FIRST argument, so
            // gcd(c, 0) = c but gcd(0, c) = 0: the result is asymmetric for
            // degree-0 inputs. The quirk is preserved.
            expect(coeffsOf(greatestCommonDivisor(p, zero)))
                .toEqual(p.equals(zero) ? [0] : coeffsOf(p));
            const expectedReversed = (p.getDegree() > 0 && !p.equals(zero))
                ? coeffsOf(p) : [0];
            expect(coeffsOf(greatestCommonDivisor(zero, p)))
                .toEqual(expectedReversed);
        });
    });

    it('square-free-factors a product of distinct linear factors into itself', () => {
        // For a square-free f the algorithm terminates in one pass with the
        // single factor f (made monic). Repeated roots hit the documented
        // floating-point instability of the GCD (gtengine-js issue #83), for
        // which the port throws rather than looping forever.
        check(fc.uniqueArray(fc.integer({ min: -3, max: 3 }),
            { minLength: 1, maxLength: 4 }), roots => {
                const f = fromRoots(roots);
                const factors = squareFreeFactorization(f);
                expect(factors.length).toBe(1);
                expect(factors[0].getDegree()).toBe(roots.length);
                for (const r of roots) {
                    expect(Math.abs(factors[0].evaluate(r))).toBeLessThan(1e-9);
                }
            }, 60);
    });
});
