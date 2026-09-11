import { describe, it, expect } from 'vitest';
import { CubicRootsQR, type CubicRootsQRMatrix } from '../src/CubicRootsQR.js';
import { RootsCubic } from '../src/RootsCubic.js';
import { check, fc, nonzero, wellScaled } from './helpers/arbitraries.js';

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

// p(x) = c0 + c1*x + c2*x^2 + x^3.
function evalCubic(c0: number, c1: number, c2: number, x: number): number {
    return c0 + x * (c1 + x * (c2 + x));
}

// Monic cubic coefficients from roots r0, r1, r2.
function coeffFromRoots(r0: number, r1: number, r2: number): [number, number, number] {
    return [-r0 * r1 * r2, r0 * r1 + r0 * r2 + r1 * r2, -(r0 + r1 + r2)];
}

const MAX_ITERATIONS = 1024;

describe('CubicRootsQR', () => {
    it('finds three simple roots of (x-1)(x-2)(x-3)', () => {
        const qr = new CubicRootsQR();
        const { iterations, numRoots, roots } = qr.solve(MAX_ITERATIONS, -6, 11, -6);
        expect(iterations).toBeLessThan(MAX_ITERATIONS);
        expect(numRoots).toBe(3);
        const sorted = roots.slice(0, numRoots).sort((a, b) => a - b);
        expect(Math.abs(sorted[0] - 1)).toBeLessThanOrEqual(1e-8);
        expect(Math.abs(sorted[1] - 2)).toBeLessThanOrEqual(1e-8);
        expect(Math.abs(sorted[2] - 3)).toBeLessThanOrEqual(1e-8);
    });

    it('finds the single real root of x^3 + x + 1', () => {
        const qr = new CubicRootsQR();
        const { numRoots, roots } = qr.solve(MAX_ITERATIONS, 1, 1, 0);
        expect(numRoots).toBe(1);
        expect(Math.abs(roots[0] - (-0.6823278038280193))).toBeLessThanOrEqual(1e-10);
    });

    it('finds the single real root of x^3 - 1 (c1 = c2 = 0 cycle avoidance)', () => {
        const qr = new CubicRootsQR();
        const { numRoots, roots } = qr.solve(MAX_ITERATIONS, -1, 0, 0);
        expect(numRoots).toBe(1);
        expect(Math.abs(roots[0] - 1)).toBeLessThanOrEqual(1e-10);
    });

    it('handles the triple root of x^3 = 0', () => {
        const qr = new CubicRootsQR();
        const { numRoots, roots } = qr.solve(MAX_ITERATIONS, 0, 0, 0);
        expect(numRoots).toBeGreaterThanOrEqual(1);
        for (let i = 0; i < numRoots; ++i) {
            expect(Math.abs(roots[i])).toBeLessThanOrEqual(1e-5);
        }
    });

    it('handles the triple root of (x-1)^3', () => {
        // Multiple roots are ill-conditioned (perturbation ~ eps^(1/3)),
        // so only require the residual and the root cluster to be small.
        const qr = new CubicRootsQR();
        const { numRoots, roots } = qr.solve(MAX_ITERATIONS, -1, 3, -3);
        expect(numRoots).toBeGreaterThanOrEqual(1);
        for (let i = 0; i < numRoots; ++i) {
            expect(Math.abs(roots[i] - 1)).toBeLessThanOrEqual(1e-4);
            expect(Math.abs(evalCubic(-1, 3, -3, roots[i]))).toBeLessThanOrEqual(1e-12);
        }
    });

    it('handles a double root: (x-2)^2 (x+1)', () => {
        const [c0, c1, c2] = coeffFromRoots(2, 2, -1);
        const qr = new CubicRootsQR();
        const { numRoots, roots } = qr.solve(MAX_ITERATIONS, c0, c1, c2);
        expect(numRoots).toBeGreaterThanOrEqual(1);
        const sorted = roots.slice(0, numRoots).sort((a, b) => a - b);
        // The simple root -1 is well-conditioned; a double root is accurate
        // only to about sqrt(eps).
        expect(Math.abs(sorted[0] - (-1))).toBeLessThanOrEqual(1e-8);
        for (let i = 1; i < numRoots; ++i) {
            expect(Math.abs(sorted[i] - 2)).toBeLessThanOrEqual(1e-6);
        }
    });

    it('computes eigenvalues from a caller-supplied companion matrix', () => {
        // Companion matrix of (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6:
        // c0 = -6, c1 = 11, c2 = -6.
        const A: CubicRootsQRMatrix = [
            [0, 0, 6],
            [1, 0, -11],
            [0, 1, 6]
        ];
        const qr = new CubicRootsQR();
        const { numRoots, roots } = qr.solveMatrix(MAX_ITERATIONS, A);
        expect(numRoots).toBe(3);
        const sorted = roots.slice(0, numRoots).sort((a, b) => a - b);
        expect(Math.abs(sorted[0] - 1)).toBeLessThanOrEqual(1e-8);
        expect(Math.abs(sorted[1] - 2)).toBeLessThanOrEqual(1e-8);
        expect(Math.abs(sorted[2] - 3)).toBeLessThanOrEqual(1e-8);
    });

    it('reports maxIterations when given no iteration budget', () => {
        const qr = new CubicRootsQR();
        const { iterations, numRoots } = qr.solve(0, -6, 11, -6);
        expect(iterations).toBe(0);
        expect(numRoots).toBe(0);
    });

    it('recovers randomized well-separated real roots', () => {
        const rng = makeRng(0x2f6e2b1);
        const qr = new CubicRootsQR();
        for (let trial = 0; trial < 100; ++trial) {
            // Choose three roots in [-5,5] separated by at least 0.25.
            let r: number[];
            do {
                r = [rng() * 10 - 5, rng() * 10 - 5, rng() * 10 - 5];
                r.sort((a, b) => a - b);
            } while (r[1] - r[0] < 0.25 || r[2] - r[1] < 0.25);

            const [c0, c1, c2] = coeffFromRoots(r[0], r[1], r[2]);
            const { numRoots, roots } = qr.solve(MAX_ITERATIONS, c0, c1, c2);
            expect(numRoots).toBe(3);
            const sorted = roots.slice(0, 3).sort((a, b) => a - b);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(sorted[i] - r[i])).toBeLessThanOrEqual(1e-7);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V44): property-based comparison against upstream
// CubicRootsQR.h.
// ---------------------------------------------------------------------------

const qrLattice = fc.integer({ min: -6, max: 6 });

describe('CubicRootsQR verification', () => {
    it('returns the requested iteration count and no roots when it cannot uncouple', () => {
        check(fc.tuple(wellScaled(-4, 4), wellScaled(-4, 4), wellScaled(-4, 4)),
            ([c0, c1, c2]) => {
                // With no iterations allowed the loop never runs, so upstream
                // returns maxIterations with the zero-filled roots array.
                const out = new CubicRootsQR().solve(0, c0, c1, c2);
                expect(out.iterations).toBe(0);
                expect(out.numRoots).toBe(0);
                expect(out.roots).toEqual([0, 0, 0]);
            });
    });

    it('finds all three roots of a cubic with distinct integer roots', () => {
        // Distinct integer roots keep the deflated 2x2 block well away from a
        // zero discriminant. For a repeated root the discriminant is zero in
        // exact arithmetic and can round negative, in which case upstream's
        // GetQuadraticRoots reports nothing at all for that block; the next
        // property covers what survives in that case.
        check(fc.tuple(qrLattice, qrLattice, qrLattice)
            .filter(([a, b, c]) => a !== b && b !== c && a !== c),
            ([a, b, c]) => {
                const [c0, c1, c2] = coeffFromRoots(a, b, c);
                const out = new CubicRootsQR().solve(MAX_ITERATIONS, c0, c1, c2);
                expect(out.numRoots).toBe(3);
                const found = out.roots.slice(0, 3).sort((x, y) => x - y);
                const expected = [a, b, c].sort((x, y) => x - y);
                for (let i = 0; i < 3; ++i) {
                    expect(Math.abs(found[i] - expected[i]))
                        .toBeLessThanOrEqual(1e-8 * (1 + Math.abs(expected[i])));
                }
            }, 100);
    });

    it('reports only genuine roots for any integer-rooted cubic', () => {
        check(fc.tuple(qrLattice, qrLattice, qrLattice), ([a, b, c]) => {
            const [c0, c1, c2] = coeffFromRoots(a, b, c);
            const out = new CubicRootsQR().solve(MAX_ITERATIONS, c0, c1, c2);
            // A cubic always has at least one real root and the QR iteration
            // always uncouples, so at least one is reported.
            expect(out.numRoots).toBeGreaterThanOrEqual(1);
            const expected = [a, b, c];
            for (let i = 0; i < out.numRoots; ++i) {
                const f = out.roots[i];
                const err = Math.min(...expected.map(e => Math.abs(f - e)));
                // A double root is only accurate to about sqrt(eps).
                expect(err).toBeLessThanOrEqual(1e-4 * (1 + Math.abs(f)));
            }
        }, 100);
    });

    it('reports residuals at the level of the root conditioning', () => {
        check(fc.tuple(nonzero(-4, 4, 0.1), nonzero(-4, 4, 0.1), nonzero(-4, 4, 0.1)),
            ([c0, c1, c2]) => {
                const out = new CubicRootsQR().solve(MAX_ITERATIONS, c0, c1, c2);
                for (let i = 0; i < out.numRoots; ++i) {
                    const x = out.roots[i];
                    const slope = Math.abs(c1 + 2 * c2 * x + 3 * x * x);
                    if (slope > 1e-2) {
                        expect(Math.abs(evalCubic(c0, c1, c2, x)) / slope)
                            .toBeLessThan(1e-8);
                    }
                }
            });
    });

    it('agrees with the exact RootsCubic classifier on which roots are real', () => {
        check(fc.tuple(nonzero(-4, 4, 0.2), nonzero(-4, 4, 0.2), nonzero(-4, 4, 0.2)),
            ([c0, c1, c2]) => {
                const qr = new CubicRootsQR().solve(MAX_ITERATIONS, c0, c1, c2);
                const exact = RootsCubic.solveMonic(true, c0, c1, c2);
                const found = qr.roots.slice(0, qr.numRoots).sort((x, y) => x - y);
                // Every QR root must be one of the exactly classified roots.
                for (const f of found) {
                    const err = Math.min(...exact.map(r => Math.abs(f - r.x)));
                    expect(err).toBeLessThanOrEqual(1e-4 * (1 + Math.abs(f)));
                }
                // The converse holds only when the roots are well separated:
                // a nearly double root leaves the deflated 2x2 block with a
                // discriminant that can round negative, and upstream then
                // reports no root for that block at all.
                const wellSeparated = exact.every(r => r.m === 1)
                    && exact.every((r, i) =>
                        i === 0 || r.x - exact[i - 1].x > 1e-2);
                if (wellSeparated) {
                    for (const r of exact) {
                        const err = Math.min(...found.map(f => Math.abs(f - r.x)));
                        expect(err).toBeLessThanOrEqual(1e-6 * (1 + Math.abs(r.x)));
                    }
                }
            }, 100);
    });

    it('solves the companion matrix in place and matches the coefficient entry point', () => {
        check(fc.tuple(wellScaled(-4, 4), wellScaled(-4, 4), wellScaled(-4, 4)),
            ([c0, c1, c2]) => {
                // Reproduce what solve() builds: the companion matrix, without
                // the fixed first Householder step that breaks the QR cycle.
                const A: CubicRootsQRMatrix = [[0, 0, -c0], [1, 0, -c1], [0, 1, -c2]];
                const before = A.map(row => [...row]);
                const viaMatrix = new CubicRootsQR().solveMatrix(MAX_ITERATIONS, A);
                // solveMatrix mutates its argument, as upstream documents.
                expect(A).not.toEqual(before);

                const viaCoefficients =
                    new CubicRootsQR().solve(MAX_ITERATIONS, c0, c1, c2);
                // The two entry points differ only by that extra iteration, so
                // both must report real roots of the same polynomial.
                for (const out of [viaCoefficients, viaMatrix]) {
                    for (let i = 0; i < out.numRoots; ++i) {
                        const x = out.roots[i];
                        const slope = Math.abs(c1 + 2 * c2 * x + 3 * x * x);
                        if (slope > 1e-2) {
                            expect(Math.abs(evalCubic(c0, c1, c2, x)) / slope)
                                .toBeLessThan(1e-8);
                        }
                    }
                }
            });
    });

    it('handles the QR-cycle coefficients c1 = c2 = 0', () => {
        // Upstream applies one fixed Householder reflection before the loop
        // precisely so that x^3 + c0 does not cycle.
        check(nonzero(-100, 100, 1e-2), c0 => {
            const out = new CubicRootsQR().solve(MAX_ITERATIONS, c0, 0, 0);
            expect(out.numRoots).toBeGreaterThanOrEqual(1);
            const expected = -Math.cbrt(c0);
            const err = Math.min(...out.roots.slice(0, out.numRoots)
                .map(f => Math.abs(f - expected)));
            expect(err).toBeLessThanOrEqual(1e-8 * (1 + Math.abs(expected)));
        });
    });

    it('never reports a NaN root and zero-fills the unused entries', () => {
        check(fc.tuple(wellScaled(-6, 6), wellScaled(-6, 6), wellScaled(-6, 6)),
            ([c0, c1, c2]) => {
                const out = new CubicRootsQR().solve(MAX_ITERATIONS, c0, c1, c2);
                expect(out.roots.length).toBe(3);
                expect(out.numRoots).toBeGreaterThanOrEqual(0);
                expect(out.numRoots).toBeLessThanOrEqual(3);
                for (let i = 0; i < out.numRoots; ++i) {
                    expect(Number.isFinite(out.roots[i])).toBe(true);
                }
                for (let i = out.numRoots; i < 3; ++i) {
                    expect(out.roots[i] === 0).toBe(true);
                }
            });
    });
});
