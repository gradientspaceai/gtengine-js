import { describe, it, expect } from 'vitest';
import { UnsymmetricEigenvalues } from '../src/UnsymmetricEigenvalues.js';

// Build the (row-major) companion matrix of the monic polynomial
// x^n + c[n-1]*x^{n-1} + ... + c[1]*x + c[0]. Its eigenvalues are exactly
// the roots of the polynomial.
function companionMatrix(c: number[]): number[] {
    const n = c.length;
    const m = new Array<number>(n * n).fill(0);
    for (let r = 1; r < n; ++r) {
        m[(r - 1) + r * n] = 1;  // subdiagonal A(r, r-1) = 1
    }
    for (let r = 0; r < n; ++r) {
        m[(n - 1) + r * n] = -c[r];  // last column A(r, n-1) = -c[r]
    }
    return m;
}

// Expand monic polynomial coefficients from its roots.
function coefficientsFromRoots(roots: number[]): number[] {
    let coeff = [1];
    for (const root of roots) {
        const next = new Array<number>(coeff.length + 1).fill(0);
        for (let i = 0; i < coeff.length; ++i) {
            next[i] += coeff[i] * (-root);
            next[i + 1] += coeff[i];
        }
        coeff = next;
    }
    // coeff[i] is the coefficient of x^i; drop the leading 1.
    return coeff.slice(0, coeff.length - 1);
}

function solveCompanion(roots: number[], maxIterations = 1024): { numEigenvalues: number; eigenvalues: number[]; iterations: number } {
    const n = roots.length;
    const solver = new UnsymmetricEigenvalues(n, maxIterations);
    const iterations = solver.solve(companionMatrix(coefficientsFromRoots(roots)), +1);
    const { numEigenvalues, eigenvalues } = solver.getEigenvalues();
    return { numEigenvalues, eigenvalues, iterations };
}

describe('UnsymmetricEigenvalues', () => {
    it('computes the roots of (x-1)(x-2)(x-3) from its companion matrix', () => {
        const { numEigenvalues, eigenvalues, iterations } = solveCompanion([1, 2, 3]);
        expect(iterations).toBeLessThan(1024);
        expect(numEigenvalues).toBe(3);
        expect(eigenvalues[0]).toBeCloseTo(1, 8);
        expect(eigenvalues[1]).toBeCloseTo(2, 8);
        expect(eigenvalues[2]).toBeCloseTo(3, 8);
    });

    it('computes the roots of a quartic with distinct real roots', () => {
        const roots = [-3, -0.5, 1.25, 4];
        const { numEigenvalues, eigenvalues } = solveCompanion(roots);
        expect(numEigenvalues).toBe(4);
        for (let i = 0; i < 4; ++i) {
            expect(eigenvalues[i]).toBeCloseTo(roots[i], 7);
        }
    });

    it('computes the roots of a degree-6 polynomial with distinct real roots', () => {
        const roots = [-5, -2, -1, 0.5, 3, 7];
        const { numEigenvalues, eigenvalues } = solveCompanion(roots);
        expect(numEigenvalues).toBe(6);
        for (let i = 0; i < roots.length; ++i) {
            expect(eigenvalues[i]).toBeCloseTo(roots[i], 6);
        }
    });

    it('reports only the real eigenvalues when complex pairs exist', () => {
        // p(x) = (x-1)(x-2)(x^2+1) has real roots {1, 2} and a complex
        // conjugate pair {i, -i}.
        // x^2+1 contributes coefficients via multiplication.
        const realRoots = [1, 2];
        // Multiply (x^2 + 1) into the expansion of (x-1)(x-2).
        // (x-1)(x-2) = x^2 - 3x + 2, so
        // p(x) = x^4 - 3x^3 + 3x^2 - 3x + 2 with c = [2, -3, 3, -3].
        const solver = new UnsymmetricEigenvalues(4, 1024);
        solver.solve(companionMatrix([2, -3, 3, -3]), +1);
        const { numEigenvalues, eigenvalues } = solver.getEigenvalues();
        expect(numEigenvalues).toBe(2);
        for (let i = 0; i < numEigenvalues; ++i) {
            expect(eigenvalues[i]).toBeCloseTo(realRoots[i], 8);
        }
    });

    it('supports decreasing sort order', () => {
        const solver = new UnsymmetricEigenvalues(3, 1024);
        solver.solve(companionMatrix(coefficientsFromRoots([1, 2, 3])), -1);
        const { numEigenvalues, eigenvalues } = solver.getEigenvalues();
        expect(numEigenvalues).toBe(3);
        expect(eigenvalues[0]).toBeCloseTo(3, 8);
        expect(eigenvalues[1]).toBeCloseTo(2, 8);
        expect(eigenvalues[2]).toBeCloseTo(1, 8);
    });

    it('supports no sorting (sortType 0) and returns all real eigenvalues', () => {
        const roots = [2, -4, 0.75];
        const solver = new UnsymmetricEigenvalues(3, 1024);
        solver.solve(companionMatrix(coefficientsFromRoots(roots)), 0);
        const { numEigenvalues, eigenvalues } = solver.getEigenvalues();
        expect(numEigenvalues).toBe(3);
        const sorted = eigenvalues.slice().sort((x, y) => x - y);
        const expected = roots.slice().sort((x, y) => x - y);
        for (let i = 0; i < 3; ++i) {
            expect(sorted[i]).toBeCloseTo(expected[i], 8);
        }
    });

    it('solves a general (non-companion) unsymmetric matrix', () => {
        // A 3x3 matrix with known eigenvalues {1, 2, 3}:
        // A = [[2, 0, 0], [1, 3, -1], [1, 1, 1]].
        // det(A - x I) = (2-x)((3-x)(1-x)+1) = (2-x)(x^2-4x+4) = (2-x)(x-2)^2.
        // Use instead a matrix with distinct eigenvalues: upper triangular
        // plus a similarity transform would do, but simplest is a matrix
        // whose characteristic polynomial is known. Take
        // A = [[6, -1, 0], [2, 3, 0], [0, 0, -2]] with eigenvalues
        // {4, 5, -2} (the 2x2 block has trace 9, det 20).
        const a = [
            6, -1, 0,
            2, 3, 0,
            0, 0, -2
        ];
        const solver = new UnsymmetricEigenvalues(3, 1024);
        solver.solve(a, +1);
        const { numEigenvalues, eigenvalues } = solver.getEigenvalues();
        expect(numEigenvalues).toBe(3);
        expect(eigenvalues[0]).toBeCloseTo(-2, 8);
        expect(eigenvalues[1]).toBeCloseTo(4, 8);
        expect(eigenvalues[2]).toBeCloseTo(5, 8);
    });

    it('handles repeated roots (reports the well-separated root)', () => {
        // Repeated roots are ill-conditioned for the QR iteration: the 2x2
        // diagonal block for the double root 2 converges with a slightly
        // negative discriminant, so upstream classifies it as a complex
        // pair and omits it. Every eigenvalue that is reported must match
        // a true root, and the simple root 5 must be found.
        const { numEigenvalues, eigenvalues } = solveCompanion([2, 2, 5]);
        expect(numEigenvalues).toBeGreaterThanOrEqual(1);
        for (let i = 0; i < numEigenvalues; ++i) {
            const nearestDist = Math.min(Math.abs(eigenvalues[i] - 2), Math.abs(eigenvalues[i] - 5));
            expect(nearestDist).toBeLessThanOrEqual(1e-4);
        }
        expect(Math.abs(eigenvalues[numEigenvalues - 1] - 5)).toBeLessThanOrEqual(1e-8);
    });

    it('returns 0 iterations and no eigenvalues for invalid construction', () => {
        const solver = new UnsymmetricEigenvalues(2, 1024);  // size < 3
        const iterations = solver.solve([1, 0, 0, 1], +1);
        expect(iterations).toBe(0);
        const { numEigenvalues, eigenvalues } = solver.getEigenvalues();
        expect(numEigenvalues).toBe(0);
        expect(eigenvalues).toEqual([]);

        const solver2 = new UnsymmetricEigenvalues(3, 0);  // maxIterations < 1
        expect(solver2.solve(companionMatrix([1, 1, 1]), +1)).toBe(0);
    });
});
