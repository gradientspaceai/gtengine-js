import { describe, it, expect } from 'vitest';
import { QuarticRootsQR, type QuarticRootsQRMatrix } from '../src/QuarticRootsQR.js';

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
