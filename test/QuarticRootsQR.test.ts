import { describe, it, expect } from 'vitest';
import { CubicRootsQR } from '../src/CubicRootsQR.js';
import { QuarticRootsQR, type QuarticRootsQRMatrix } from '../src/QuarticRootsQR.js';
import { RootsQuartic } from '../src/RootsQuartic.js';
import { check, fc, nonzero, wellScaled } from './helpers/arbitraries.js';

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

// p(x) = c0 + c1*x + c2*x^2 + c3*x^3 + x^4.
function evalQuartic(c0: number, c1: number, c2: number, c3: number, x: number): number {
    return c0 + x * (c1 + x * (c2 + x * (c3 + x)));
}

// Monic quartic coefficients from roots r0, r1, r2, r3.
function coeffFromRoots(r: number[]): [number, number, number, number] {
    const [a, b, c, d] = r;
    const c3 = -(a + b + c + d);
    const c2 = a * b + a * c + a * d + b * c + b * d + c * d;
    const c1 = -(a * b * c + a * b * d + a * c * d + b * c * d);
    const c0 = a * b * c * d;
    return [c0, c1, c2, c3];
}

// Monic quartic coefficients from the two quadratic factors
// (x^2 + p1*x + p0)(x^2 + q1*x + q0).
function coeffFromQuadratics(p1: number, p0: number, q1: number, q0: number):
    [number, number, number, number] {
    return [p0 * q0, p1 * q0 + p0 * q1, p0 + q0 + p1 * q1, p1 + q1];
}

const MAX_ITERATIONS = 1024;

describe('QuarticRootsQR', () => {
    it('finds four simple roots of (x-1)(x-2)(x-3)(x-4)', () => {
        const [c0, c1, c2, c3] = coeffFromRoots([1, 2, 3, 4]);
        const qr = new QuarticRootsQR();
        const { iterations, numRoots, roots } = qr.solve(MAX_ITERATIONS, c0, c1, c2, c3);
        expect(iterations).toBeLessThan(MAX_ITERATIONS);
        expect(numRoots).toBe(4);
        const sorted = roots.slice(0, numRoots).sort((a, b) => a - b);
        for (let i = 0; i < 4; ++i) {
            expect(Math.abs(sorted[i] - (i + 1))).toBeLessThanOrEqual(1e-8);
        }
    });

    it('finds symmetric roots of (x^2-1)(x^2-4)', () => {
        // x^4 - 5x^2 + 4
        const qr = new QuarticRootsQR();
        const { numRoots, roots } = qr.solve(MAX_ITERATIONS, 4, 0, -5, 0);
        expect(numRoots).toBe(4);
        const sorted = roots.slice(0, numRoots).sort((a, b) => a - b);
        const expected = [-2, -1, 1, 2];
        for (let i = 0; i < 4; ++i) {
            expect(Math.abs(sorted[i] - expected[i])).toBeLessThanOrEqual(1e-8);
        }
    });

    it('finds the two real roots of x^4 - 1 (c1 = c2 = 0 cycle avoidance)', () => {
        const qr = new QuarticRootsQR();
        const { numRoots, roots } = qr.solve(MAX_ITERATIONS, -1, 0, 0, 0);
        expect(numRoots).toBe(2);
        const sorted = roots.slice(0, numRoots).sort((a, b) => a - b);
        expect(Math.abs(sorted[0] - (-1))).toBeLessThanOrEqual(1e-10);
        expect(Math.abs(sorted[1] - 1)).toBeLessThanOrEqual(1e-10);
    });

    it('reports no real roots for x^4 + 1', () => {
        const qr = new QuarticRootsQR();
        const { numRoots } = qr.solve(MAX_ITERATIONS, 1, 0, 0, 0);
        expect(numRoots).toBe(0);
    });

    it('reports no real roots for (x^2+1)(x^2+2x+5)', () => {
        const [c0, c1, c2, c3] = coeffFromQuadratics(0, 1, 2, 5);
        const qr = new QuarticRootsQR();
        const { numRoots } = qr.solve(MAX_ITERATIONS, c0, c1, c2, c3);
        expect(numRoots).toBe(0);
    });

    it('finds only the two real roots of (x-1)(x+2)(x^2+x+1)', () => {
        const [c0, c1, c2, c3] = coeffFromQuadratics(1, -2, 1, 1);
        const qr = new QuarticRootsQR();
        const { numRoots, roots } = qr.solve(MAX_ITERATIONS, c0, c1, c2, c3);
        expect(numRoots).toBe(2);
        const sorted = roots.slice(0, numRoots).sort((a, b) => a - b);
        expect(Math.abs(sorted[0] - (-2))).toBeLessThanOrEqual(1e-8);
        expect(Math.abs(sorted[1] - 1)).toBeLessThanOrEqual(1e-8);
    });

    it('handles the quadruple root of x^4 = 0', () => {
        const qr = new QuarticRootsQR();
        const { numRoots, roots } = qr.solve(MAX_ITERATIONS, 0, 0, 0, 0);
        expect(numRoots).toBeGreaterThanOrEqual(1);
        for (let i = 0; i < numRoots; ++i) {
            expect(Math.abs(roots[i])).toBeLessThanOrEqual(1e-3);
        }
    });

    it('handles a double root: (x-2)^2 (x+1)(x-5)', () => {
        const [c0, c1, c2, c3] = coeffFromRoots([2, 2, -1, 5]);
        const qr = new QuarticRootsQR();
        const { numRoots, roots } = qr.solve(MAX_ITERATIONS, c0, c1, c2, c3);
        expect(numRoots).toBeGreaterThanOrEqual(2);
        const sorted = roots.slice(0, numRoots).sort((a, b) => a - b);
        // The simple roots are well-conditioned; a double root is accurate
        // only to about sqrt(eps).
        expect(Math.abs(sorted[0] - (-1))).toBeLessThanOrEqual(1e-8);
        expect(Math.abs(sorted[numRoots - 1] - 5)).toBeLessThanOrEqual(1e-8);
        for (let i = 1; i < numRoots - 1; ++i) {
            expect(Math.abs(sorted[i] - 2)).toBeLessThanOrEqual(1e-6);
        }
    });

    it('has small residuals at every reported root', () => {
        const [c0, c1, c2, c3] = coeffFromRoots([-3.25, -0.5, 1.75, 6]);
        const qr = new QuarticRootsQR();
        const { numRoots, roots } = qr.solve(MAX_ITERATIONS, c0, c1, c2, c3);
        expect(numRoots).toBe(4);
        for (let i = 0; i < numRoots; ++i) {
            expect(Math.abs(evalQuartic(c0, c1, c2, c3, roots[i])))
                .toBeLessThanOrEqual(1e-8);
        }
    });

    it('computes eigenvalues from a caller-supplied companion matrix', () => {
        // Companion matrix of (x-1)(x-2)(x-3)(x-4)
        //   = x^4 - 10x^3 + 35x^2 - 50x + 24.
        const [c0, c1, c2, c3] = coeffFromRoots([1, 2, 3, 4]);
        const A: QuarticRootsQRMatrix = [
            [0, 0, 0, -c0],
            [1, 0, 0, -c1],
            [0, 1, 0, -c2],
            [0, 0, 1, -c3]
        ];
        const qr = new QuarticRootsQR();
        const { numRoots, roots } = qr.solveMatrix(MAX_ITERATIONS, A);
        expect(numRoots).toBe(4);
        const sorted = roots.slice(0, numRoots).sort((a, b) => a - b);
        for (let i = 0; i < 4; ++i) {
            expect(Math.abs(sorted[i] - (i + 1))).toBeLessThanOrEqual(1e-7);
        }
    });

    it('reports maxIterations when given no iteration budget', () => {
        const qr = new QuarticRootsQR();
        const [c0, c1, c2, c3] = coeffFromRoots([1, 2, 3, 4]);
        const { iterations, numRoots } = qr.solve(0, c0, c1, c2, c3);
        expect(iterations).toBe(0);
        expect(numRoots).toBe(0);
    });

    it('recovers randomized well-separated real roots', () => {
        const rng = makeRng(0x4a17c);
        const qr = new QuarticRootsQR();
        for (let trial = 0; trial < 100; ++trial) {
            // Choose four roots in [-5,5] separated by at least 0.25.
            let r: number[];
            do {
                r = [rng() * 10 - 5, rng() * 10 - 5, rng() * 10 - 5, rng() * 10 - 5];
                r.sort((a, b) => a - b);
            } while (r[1] - r[0] < 0.25 || r[2] - r[1] < 0.25 || r[3] - r[2] < 0.25);

            const [c0, c1, c2, c3] = coeffFromRoots(r);
            const { numRoots, roots } = qr.solve(MAX_ITERATIONS, c0, c1, c2, c3);
            expect(numRoots).toBe(4);
            const sorted = roots.slice(0, 4).sort((a, b) => a - b);
            for (let i = 0; i < 4; ++i) {
                expect(Math.abs(sorted[i] - r[i])).toBeLessThanOrEqual(1e-6);
            }
        }
    });

    it('recovers two real roots when the other two are complex', () => {
        const rng = makeRng(0x2c0e1);
        const qr = new QuarticRootsQR();
        for (let trial = 0; trial < 100; ++trial) {
            // Real roots a < b with a gap, plus an irreducible quadratic
            // x^2 + q1*x + q0 with q1^2 - 4*q0 < 0.
            const a = rng() * 6 - 5;
            const b = a + 0.5 + rng() * 4;
            const q1 = rng() * 4 - 2;
            const q0 = q1 * q1 / 4 + 0.25 + rng();

            const [c0, c1, c2, c3] = coeffFromQuadratics(-(a + b), a * b, q1, q0);
            const { numRoots, roots } = qr.solve(MAX_ITERATIONS, c0, c1, c2, c3);
            expect(numRoots).toBe(2);
            const sorted = roots.slice(0, 2).sort((x, y) => x - y);
            expect(Math.abs(sorted[0] - a)).toBeLessThanOrEqual(1e-6);
            expect(Math.abs(sorted[1] - b)).toBeLessThanOrEqual(1e-6);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V44): property-based comparison against upstream
// QuarticRootsQR.h.
// ---------------------------------------------------------------------------

const qrLattice = fc.integer({ min: -5, max: 5 });

function derivative(c0: number, c1: number, c2: number, c3: number,
    x: number): number {
    return c1 + x * (2 * c2 + x * (3 * c3 + x * 4));
}

describe('QuarticRootsQR verification', () => {
    it('returns the requested iteration count and no roots when it cannot uncouple', () => {
        check(fc.tuple(wellScaled(-4, 4), wellScaled(-4, 4), wellScaled(-4, 4),
            wellScaled(-4, 4)),
            ([c0, c1, c2, c3]) => {
                const out = new QuarticRootsQR().solve(0, c0, c1, c2, c3);
                expect(out.iterations).toBe(0);
                expect(out.numRoots).toBe(0);
                expect(out.roots).toEqual([0, 0, 0, 0]);
            });
    });

    it('finds all four roots of a quartic with distinct integer roots', () => {
        check(fc.tuple(qrLattice, qrLattice, qrLattice, qrLattice)
            .filter(rs => new Set(rs).size === 4),
            rs => {
                const [c0, c1, c2, c3] = coeffFromRoots([...rs]);
                const out = new QuarticRootsQR().solve(MAX_ITERATIONS, c0, c1, c2, c3);
                expect(out.numRoots).toBe(4);
                const found = out.roots.slice(0, 4).sort((x, y) => x - y);
                const expected = [...rs].sort((x, y) => x - y);
                for (let i = 0; i < 4; ++i) {
                    expect(Math.abs(found[i] - expected[i]))
                        .toBeLessThanOrEqual(1e-7 * (1 + Math.abs(expected[i])));
                }
            }, 100);
    });

    it('reports only genuine roots for any integer-rooted quartic', () => {
        check(fc.tuple(qrLattice, qrLattice, qrLattice, qrLattice), rs => {
            const [c0, c1, c2, c3] = coeffFromRoots([...rs]);
            const out = new QuarticRootsQR().solve(MAX_ITERATIONS, c0, c1, c2, c3);
            expect(out.numRoots).toBeGreaterThanOrEqual(0);
            expect(out.numRoots).toBeLessThanOrEqual(4);
            for (let i = 0; i < out.numRoots; ++i) {
                const f = out.roots[i];
                const err = Math.min(...rs.map(e => Math.abs(f - e)));
                // A root of multiplicity m carries about eps^(1/m) accuracy.
                expect(err).toBeLessThanOrEqual(1e-3 * (1 + Math.abs(f)));
            }
            for (let i = out.numRoots; i < 4; ++i) {
                expect(out.roots[i] === 0).toBe(true);
            }
        }, 100);
    });

    it('reports no real roots for a product of two irreducible quadratics', () => {
        check(fc.tuple(fc.integer({ min: -3, max: 3 }), fc.integer({ min: 1, max: 8 }),
            fc.integer({ min: -3, max: 3 }), fc.integer({ min: 1, max: 8 }))
            .filter(([p1, p0, q1, q0]) =>
                p1 * p1 - 4 * p0 < 0 && q1 * q1 - 4 * q0 < 0),
            ([p1, p0, q1, q0]) => {
                const [c0, c1, c2, c3] = coeffFromQuadratics(p1, p0, q1, q0);
                const out = new QuarticRootsQR().solve(MAX_ITERATIONS, c0, c1, c2, c3);
                expect(out.numRoots).toBe(0);
            }, 100);
    });

    it('reports residuals at the level of the root conditioning', () => {
        check(fc.tuple(nonzero(-3, 3, 0.1), nonzero(-3, 3, 0.1), nonzero(-3, 3, 0.1),
            nonzero(-3, 3, 0.1)),
            ([c0, c1, c2, c3]) => {
                const out = new QuarticRootsQR().solve(MAX_ITERATIONS, c0, c1, c2, c3);
                for (let i = 0; i < out.numRoots; ++i) {
                    const x = out.roots[i];
                    expect(Number.isFinite(x)).toBe(true);
                    const slope = Math.abs(derivative(c0, c1, c2, c3, x));
                    if (slope > 1e-2) {
                        expect(Math.abs(evalQuartic(c0, c1, c2, c3, x)) / slope)
                            .toBeLessThan(1e-7);
                    }
                }
            });
    });

    it('agrees with the exact RootsQuartic classifier on well-separated roots', () => {
        check(fc.tuple(nonzero(-3, 3, 0.2), nonzero(-3, 3, 0.2), nonzero(-3, 3, 0.2),
            nonzero(-3, 3, 0.2)),
            ([c0, c1, c2, c3]) => {
                const qr = new QuarticRootsQR().solve(MAX_ITERATIONS, c0, c1, c2, c3);
                const exact = RootsQuartic.solveMonic(true, c0, c1, c2, c3);
                const found = qr.roots.slice(0, qr.numRoots);
                for (const f of found) {
                    const err = Math.min(...exact.map(r => Math.abs(f - r.x)),
                        Number.POSITIVE_INFINITY);
                    expect(err).toBeLessThanOrEqual(1e-3 * (1 + Math.abs(f)));
                }
                // The converse holds only for simple, well-separated roots:
                // see the characterization test about even multiplicities.
                const wellSeparated = exact.length > 0 && exact.every(r => r.m === 1)
                    && exact.every((r, i) => i === 0 || r.x - exact[i - 1].x > 1e-2);
                if (wellSeparated) {
                    for (const r of exact) {
                        const err = Math.min(...found.map(f => Math.abs(f - r.x)),
                            Number.POSITIVE_INFINITY);
                        expect(err).toBeLessThanOrEqual(1e-5 * (1 + Math.abs(r.x)));
                    }
                }
            }, 100);
    });

    it('solves the companion matrix in place', () => {
        check(fc.tuple(wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3),
            wellScaled(-3, 3)),
            ([c0, c1, c2, c3]) => {
                const A: QuarticRootsQRMatrix = [
                    [0, 0, 0, -c0], [1, 0, 0, -c1], [0, 1, 0, -c2], [0, 0, 1, -c3]];
                const before = A.map(row => [...row]);
                const out = new QuarticRootsQR().solveMatrix(MAX_ITERATIONS, A);
                // solveMatrix mutates its argument, as upstream documents.
                expect(A).not.toEqual(before);
                for (let i = 0; i < out.numRoots; ++i) {
                    const x = out.roots[i];
                    const slope = Math.abs(derivative(c0, c1, c2, c3, x));
                    if (slope > 1e-2) {
                        expect(Math.abs(evalQuartic(c0, c1, c2, c3, x)) / slope)
                            .toBeLessThan(1e-7);
                    }
                }
            });
    });
});

describe('QuarticRootsQR / CubicRootsQR even-multiplicity characterization', () => {
    it('drops real roots of even multiplicity (upstream behavior, preserved)', () => {
        // The deflated 2x2 block of a repeated real root has a zero
        // discriminant in exact arithmetic; rounding makes it slightly
        // negative and upstream's GetQuadraticRoots then emits nothing for
        // that block. The roots are lost silently: the caller cannot tell
        // this case from a genuine complex-conjugate pair. Recorded as an
        // upstream suspect and preserved rather than "fixed" with a
        // tolerance, since upstream defines the behavior.
        const cubic = new CubicRootsQR().solve(MAX_ITERATIONS, 36, 24, -11);
        // (x - 6)^2 (x + 1) = x^3 - 11x^2 + 24x + 36.
        expect(cubic.numRoots).toBe(1);
        expect(cubic.roots[0]).toBeCloseTo(-1, 10);

        // (x - 1)^2 (x - 2)^2 = x^4 - 6x^3 + 13x^2 - 12x + 4: only the pair
        // near 2 survives, split into two nearby simple roots.
        const quartic = new QuarticRootsQR().solve(MAX_ITERATIONS, 4, -12, 13, -6);
        expect(quartic.numRoots).toBe(2);
        for (let i = 0; i < 2; ++i) {
            expect(Math.abs(quartic.roots[i] - 2)).toBeLessThan(1e-6);
        }
    });
});
