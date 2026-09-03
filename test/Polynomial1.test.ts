import { describe, it, expect } from 'vitest';
import {
    Polynomial1, greatestCommonDivisor, squareFreeFactorization
} from '../src/Polynomial1.js';

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
