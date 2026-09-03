import { describe, it, expect } from 'vitest';
import { CubicRootsQR, type CubicRootsQRMatrix } from '../src/CubicRootsQR.js';

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
