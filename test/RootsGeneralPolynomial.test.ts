import { describe, it, expect } from 'vitest';
import { BSRational } from '../src/BSRational';
import { RootsGeneralPolynomial } from '../src/RootsGeneralPolynomial';

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

// The multiplication of the polynomial p by the quadratic x^2 + b*x + c,
// which has no real roots when b^2 < 4*c.
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
    let result = 0;
    for (let i = p.length - 1; i >= 0; --i) {
        result = x * result + p[i];
    }
    return result;
}

// A deterministic pseudorandom generator, so failures reproduce.
function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

describe('RootsGeneralPolynomial.solve degenerate inputs', () => {
    it('reports no roots for constants, including the zero polynomial', () => {
        expect(RootsGeneralPolynomial.solve([])).toEqual([]);
        expect(RootsGeneralPolynomial.solve([0])).toEqual([]);
        expect(RootsGeneralPolynomial.solve([5])).toEqual([]);
        // High-order zero coefficients are removed, leaving a constant.
        expect(RootsGeneralPolynomial.solve([5, 0, 0, 0])).toEqual([]);
        expect(RootsGeneralPolynomial.solve([0, 0, 0])).toEqual([]);
    });

    it('solves linear polynomials exactly', () => {
        expect(RootsGeneralPolynomial.solve([-6, 3])).toEqual([2]);
        expect(RootsGeneralPolynomial.solve([1, 4])).toEqual([-0.25]);
        // Trailing zero coefficients reduce the degree to 1.
        expect(RootsGeneralPolynomial.solve([-6, 3, 0, 0])).toEqual([2]);
    });

    it('removes high-order zero coefficients of a higher-degree polynomial', () => {
        // The port fixes an upstream allocation bug here; see the file header.
        expect(RootsGeneralPolynomial.solve([-1, 0, 1, 0, 0, 0]))
            .toEqual(RootsGeneralPolynomial.solve([-1, 0, 1]));
        expect(RootsGeneralPolynomial.solve([-6, 11, -6, 1, 0]))
            .toEqual(RootsGeneralPolynomial.solve([-6, 11, -6, 1]));
    });

    it('ignores the useThreading argument', () => {
        const p = fromRoots([-3, -1, 2, 5]);
        expect(RootsGeneralPolynomial.solve(p, true))
            .toEqual(RootsGeneralPolynomial.solve(p, false));
    });
});

describe('RootsGeneralPolynomial.solve known factorizations', () => {
    it('finds the roots of quadratics', () => {
        expect(RootsGeneralPolynomial.solve([-1, 0, 1])).toEqual([-1, 1]);
        // x^2 - 5x + 6 = (x-2)(x-3)
        const roots = RootsGeneralPolynomial.solve([6, -5, 1]);
        expect(roots.length).toBe(2);
        expect(roots[0]).toBeCloseTo(2, 14);
        expect(roots[1]).toBeCloseTo(3, 14);
        // x^2 + 1 has no real roots.
        expect(RootsGeneralPolynomial.solve([1, 0, 1])).toEqual([]);
        // x^2 = 0 has the double root 0, reported once.
        expect(RootsGeneralPolynomial.solve([0, 0, 1])).toEqual([0]);
    });

    it('finds the roots of cubics and quartics', () => {
        // (x+3)(x-1)(x-4)
        let roots = RootsGeneralPolynomial.solve(fromRoots([-3, 1, 4]));
        expect(roots.length).toBe(3);
        expect(roots[0]).toBeCloseTo(-3, 12);
        expect(roots[1]).toBeCloseTo(1, 12);
        expect(roots[2]).toBeCloseTo(4, 12);

        // (x-2)(x^2+1): one real root.
        roots = RootsGeneralPolynomial.solve(timesQuadratic(fromRoots([2]), 0, 1));
        expect(roots.length).toBe(1);
        expect(roots[0]).toBeCloseTo(2, 12);

        // (x+1)(x-0.5)(x-2)(x-7)
        roots = RootsGeneralPolynomial.solve(fromRoots([-1, 0.5, 2, 7]));
        expect(roots.length).toBe(4);
        expect(roots[0]).toBeCloseTo(-1, 12);
        expect(roots[1]).toBeCloseTo(0.5, 12);
        expect(roots[2]).toBeCloseTo(2, 12);
        expect(roots[3]).toBeCloseTo(7, 12);

        // (x^2+1)(x^2+x+1): no real roots.
        expect(RootsGeneralPolynomial.solve(
            timesQuadratic(timesQuadratic([1], 0, 1), 1, 1))).toEqual([]);
    });

    it('reports a repeated root once, whatever the multiplicity', () => {
        // (x-1)^2
        expect(RootsGeneralPolynomial.solve([1, -2, 1])).toEqual([1]);
        // (x-1)^3
        expect(RootsGeneralPolynomial.solve([-1, 3, -3, 1])).toEqual([1]);
        // (x-1)^2 (x-3): the double root is an exact root of the derivative.
        let roots = RootsGeneralPolynomial.solve([-3, 7, -5, 1]);
        expect(roots.length).toBe(2);
        expect(roots[0]).toBeCloseTo(1, 12);
        expect(roots[1]).toBeCloseTo(3, 12);
        // (x-1)^3 (x-2): a triple root loses about two thirds of the
        // significant digits, as expected for an ill-conditioned root.
        roots = RootsGeneralPolynomial.solve([2, -7, 9, -5, 1]);
        expect(roots.length).toBe(2);
        expect(roots[0]).toBeCloseTo(1, 4);
        expect(roots[1]).toBeCloseTo(2, 12);
    });

    it('handles tiny and huge roots', () => {
        // (x-1e-3)(x-2e-3)(x-3e-3)
        const small = RootsGeneralPolynomial.solve([-6e-9, 1.1e-5, -0.006, 1]);
        expect(small.length).toBe(3);
        expect(small[0]).toBeCloseTo(0.001, 12);
        expect(small[1]).toBeCloseTo(0.002, 12);
        expect(small[2]).toBeCloseTo(0.003, 12);

        // x^2 - 1e12
        const large = RootsGeneralPolynomial.solve([-1e12, 0, 1]);
        expect(large).toEqual([-1e6, 1e6]);
    });

    it('returns the roots in increasing order', () => {
        const roots = RootsGeneralPolynomial.solve(fromRoots([5, -2, 3, -7, 0]));
        for (let i = 1; i < roots.length; ++i) {
            expect(roots[i]).toBeGreaterThan(roots[i - 1]);
        }
        expect(roots.length).toBe(5);
    });
});

describe('RootsGeneralPolynomial.solve ill-conditioned inputs', () => {
    it('solves the Wilkinson-style polynomial of degree 6', () => {
        // (x-1)(x-2)...(x-6)
        const p = [720, -1764, 1624, -735, 175, -21, 1];
        const roots = RootsGeneralPolynomial.solve(p);
        expect(roots.length).toBe(6);
        for (let i = 0; i < 6; ++i) {
            expect(roots[i]).toBeCloseTo(i + 1, 8);
        }
    });

    it('solves the Wilkinson-style polynomial of degree 8', () => {
        // (x-1)(x-2)...(x-8)
        const p = [40320, -109584, 118124, -67284, 22449, -4536, 546, -36, 1];
        const roots = RootsGeneralPolynomial.solve(p);
        expect(roots.length).toBe(8);
        for (let i = 0; i < 8; ++i) {
            expect(Math.abs(roots[i] - (i + 1))).toBeLessThan(1e-8);
        }
    });

    it('solves a polynomial with widely separated scales', () => {
        // (x - 1e-4)(x - 1)(x - 1e4)
        const roots = RootsGeneralPolynomial.solve(fromRoots([1e-4, 1, 1e4]));
        expect(roots.length).toBe(3);
        expect(roots[0]).toBeCloseTo(1e-4, 12);
        expect(roots[1]).toBeCloseTo(1, 8);
        expect(roots[2] / 1e4).toBeCloseTo(1, 10);
    });
});

describe('RootsGeneralPolynomial.solve randomized cross-checks', () => {
    it('recovers well-separated random roots and leaves small residuals', () => {
        const rand = makeRandom(20260902);
        for (let trial = 0; trial < 60; ++trial) {
            const degree = 2 + (trial % 4);

            // Well-separated roots in [-5,5], sorted increasingly.
            const roots: number[] = [];
            let x = -5;
            for (let i = 0; i < degree; ++i) {
                x += 0.6 + 2 * rand();
                roots.push(x);
            }

            const p = fromRoots(roots);
            const computed = RootsGeneralPolynomial.solve(p);
            expect(computed.length).toBe(degree);
            for (let i = 0; i < degree; ++i) {
                expect(Math.abs(computed[i] - roots[i]))
                    .toBeLessThan(1e-9 * (1 + Math.abs(roots[i])));
            }

            // The residuals are small relative to the coefficient scale.
            const scale = p.reduce((a, c) => Math.max(a, Math.abs(c)), 0);
            for (const r of computed) {
                expect(Math.abs(evaluate(p, r))).toBeLessThan(1e-8 * scale);
            }
        }
    });

    it('ignores complex-conjugate factors while keeping the real roots', () => {
        const rand = makeRandom(987654321);
        for (let trial = 0; trial < 30; ++trial) {
            const roots = [-3 - rand(), 0.5 + rand(), 4 + rand()];
            // x^2 + b*x + c with b^2 < 4c has no real roots.
            const c = 1 + 3 * rand();
            const b = (2 * rand() - 1) * Math.sqrt(4 * c) * 0.9;
            const p = timesQuadratic(fromRoots(roots), b, c);
            const computed = RootsGeneralPolynomial.solve(p);
            expect(computed.length).toBe(3);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(computed[i] - roots[i]))
                    .toBeLessThan(1e-9 * (1 + Math.abs(roots[i])));
            }
        }
    });

    it('is invariant under scaling of the polynomial', () => {
        const rand = makeRandom(13579);
        for (let trial = 0; trial < 20; ++trial) {
            const roots = [-2 - rand(), rand(), 3 + rand()];
            const p = fromRoots(roots);
            const scale = trial % 2 === 0 ? -7.5 : 1e3;
            const scaled = p.map(c => scale * c);
            const a = RootsGeneralPolynomial.solve(p);
            const b = RootsGeneralPolynomial.solve(scaled);
            expect(b.length).toBe(a.length);
            for (let i = 0; i < a.length; ++i) {
                expect(Math.abs(a[i] - b[i])).toBeLessThan(1e-9 * (1 + Math.abs(a[i])));
            }
        }
    });
});

describe('RootsGeneralPolynomial.solveRational', () => {
    const rat = (x: number): BSRational => BSRational.fromNumber(x);

    it('reports no roots for constants and solves linear polynomials exactly', () => {
        expect(RootsGeneralPolynomial.solveRational([])).toEqual([]);
        expect(RootsGeneralPolynomial.solveRational([rat(3)])).toEqual([]);
        expect(RootsGeneralPolynomial.solveRational(
            [rat(3), rat(0), rat(0)])).toEqual([]);

        // 3x - 6 = 0 has the exact rational root 2.
        let roots = RootsGeneralPolynomial.solveRational([rat(-6), rat(3)]);
        expect(roots.length).toBe(1);
        expect(roots[0].toNumber()).toBe(2);

        // 3x + 1 = 0 has the exact rational root -1/3, which is not a
        // floating-point number.
        roots = RootsGeneralPolynomial.solveRational([rat(1), rat(3)]);
        expect(roots.length).toBe(1);
        expect(roots[0].mul(rat(3)).add(rat(1)).getSign()).toBe(0);
        expect(roots[0].toNumber()).toBeCloseTo(-1 / 3, 15);
    });

    it('matches the floating-point solver on quadratics and cubics', () => {
        const cases: number[][] = [
            [-1, 0, 1],
            [6, -5, 1],
            [1, 0, 1],
            fromRoots([-3, 1, 4]),
            fromRoots([-1, 0.5, 2, 7])
        ];
        for (const p of cases) {
            const expected = RootsGeneralPolynomial.solve(p);
            const roots = RootsGeneralPolynomial.solveRational(p.map(rat));
            expect(roots.map(r => r.toNumber())).toEqual(expected);
        }
    });

    it('finds exact rational roots of a monic polynomial with rational data', () => {
        // (x - 1/2)(x - 3/4) = x^2 - (5/4) x + 3/8, using exact rationals.
        const half = BSRational.fromNumber(1).div(rat(2));
        const threeQuarters = BSRational.fromNumber(3).div(rat(4));
        const p = [half.mul(threeQuarters),
            half.add(threeQuarters).negated(), rat(1)];
        const roots = RootsGeneralPolynomial.solveRational(p);
        expect(roots.length).toBe(2);
        // The roots are found by floating-point bisection and converted back
        // to rationals, so they are the nearest floating-point numbers (or a
        // neighbor of one) rather than the exact rationals.
        expect(roots[0].toNumber()).toBeCloseTo(0.5, 15);
        expect(roots[1].toNumber()).toBeCloseTo(0.75, 15);
    });

    it('removes high-order zero coefficients', () => {
        const roots = RootsGeneralPolynomial.solveRational(
            [rat(-1), rat(0), rat(1), rat(0), rat(0)]);
        expect(roots.map(r => r.toNumber())).toEqual([-1, 1]);
    });

    it('does not modify its input array of rationals', () => {
        const p = [rat(-6), rat(11), rat(-6), rat(2)];
        const before = p.map(r => r.toNumber());
        RootsGeneralPolynomial.solveRational(p);
        expect(p.map(r => r.toNumber())).toEqual(before);
    });
});
