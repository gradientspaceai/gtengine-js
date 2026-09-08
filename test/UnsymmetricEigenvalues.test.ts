import { describe, it, expect } from 'vitest';
import { UnsymmetricEigenvalues } from '../src/UnsymmetricEigenvalues.js';
import { SymmetricEigensolver3x3 } from '../src/SymmetricEigensolver3x3.js';
import { check, expectClose, fc, nonzero, scaled, wellScaled } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (group V38).
// ---------------------------------------------------------------------------
describe('UnsymmetricEigenvalues verification', () => {
    // Row-major companion matrix of the monic polynomial with the given roots.
    // Its characteristic polynomial is exactly prod(x - roots[i]).
    function companionOf(roots: readonly number[]): number[] {
        const n = roots.length;
        // Expand prod(x - r) = c[0] + c[1] x + ... + x^n.
        let c = [1];
        for (const r of roots) {
            const next = new Array<number>(c.length + 1).fill(0);
            for (let i = 0; i < c.length; ++i) {
                next[i] -= r * c[i];
                next[i + 1] += c[i];
            }
            c = next;
        }
        const M = new Array<number>(n * n).fill(0);
        for (let i = 1; i < n; ++i) {
            M[i * n + (i - 1)] = 1;
        }
        for (let i = 0; i < n; ++i) {
            M[i * n + (n - 1)] = -c[i];
        }
        return M;
    }

    // |det(A - lambda*I)| by plain Gaussian elimination with partial pivoting.
    function charPolyValue(n: number, M: readonly number[], lambda: number): number {
        const a: number[] = M.slice();
        for (let i = 0; i < n; ++i) {
            a[i * n + i] -= lambda;
        }
        let det = 1;
        for (let col = 0; col < n; ++col) {
            let piv = col;
            for (let r = col + 1; r < n; ++r) {
                if (Math.abs(a[r * n + col]) > Math.abs(a[piv * n + col])) {
                    piv = r;
                }
            }
            if (a[piv * n + col] === 0) {
                return 0;
            }
            if (piv !== col) {
                for (let c = 0; c < n; ++c) {
                    const t = a[piv * n + c];
                    a[piv * n + c] = a[col * n + c];
                    a[col * n + c] = t;
                }
                det = -det;
            }
            det *= a[col * n + col];
            for (let r = col + 1; r < n; ++r) {
                const f = a[r * n + col] / a[col * n + col];
                for (let c = col; c < n; ++c) {
                    a[r * n + c] -= f * a[col * n + c];
                }
            }
        }
        return det;
    }

    // Distinct integer-lattice roots, so the companion matrix is well
    // conditioned and every root is well separated from the others.
    const distinctRoots = (n: number): fc.Arbitrary<number[]> =>
        fc.uniqueArray(fc.integer({ min: -6, max: 6 }),
            { minLength: n, maxLength: n });

    it('recovers the roots of a polynomial from its companion matrix', () => {
        check(fc.tuple(fc.integer({ min: 3, max: 6 }), fc.nat()).chain(
            ([n, seed]) => fc.tuple(fc.constant(n), fc.constant(seed),
                distinctRoots(n))).map(([n, , roots]) => ({ n, roots })),
            ({ n, roots }) => {
                const M = companionOf(roots);
                const solver = new UnsymmetricEigenvalues(n, 8192);
                solver.solve(M, +1);
                const { numEigenvalues, eigenvalues } = solver.getEigenvalues();
                const sortedRoots = roots.slice().sort((x, y) => x - y);
                // Every eigenvalue reported must be one of the true roots.
                // The companion matrix of a degree-6 polynomial with roots up
                // to 6 has entries in the thousands, so scale the tolerance.
                const scale = Math.max(...M.map(Math.abs), 1);
                for (let i = 0; i < numEigenvalues; ++i) {
                    const nearest = Math.min(...sortedRoots.map(
                        r => Math.abs(r - eigenvalues[i])));
                    expect(nearest).toBeLessThanOrEqual(1e-6 * scale);
                }
                // Increasing order was requested.
                for (let i = 0; i + 1 < numEigenvalues; ++i) {
                    expect(eigenvalues[i]).toBeLessThanOrEqual(eigenvalues[i + 1]);
                }
                // Up to n = 4 the Francis iteration always decouples every
                // block for these companion matrices (measured over 400
                // random root sets per size). For n = 5 and n = 6 it can
                // exhaust its iterations on one block and report fewer
                // eigenvalues -- the case upstream flags in its own comment
                // ("This happened with root finding using the companion
                // matrix of a polynomial") -- so only the subset property
                // above is asserted there.
                if (n <= 4) {
                    expect(numEigenvalues).toBe(n);
                    for (let i = 0; i < n; ++i) {
                        expectClose(eigenvalues[i], sortedRoots[i],
                            1e-6 * scale, 1e-9);
                    }
                }
            }, 100);
    });

    it('every reported eigenvalue is a root of the characteristic polynomial',
        () => {
            check(fc.tuple(fc.integer({ min: 3, max: 5 }), fc.nat({ max: 1e9 }))
                .chain(([n]) => fc.tuple(fc.constant(n),
                    fc.array(scaled(-4, 4), { minLength: n * n, maxLength: n * n })))
                .map(([n, entries]) => ({ n, M: entries })),
                ({ n, M }) => {
                    const solver = new UnsymmetricEigenvalues(n, 8192);
                    solver.solve(M, 0);
                    const { numEigenvalues, eigenvalues } =
                        solver.getEigenvalues();
                    expect(numEigenvalues).toBeLessThanOrEqual(n);
                    const scale = Math.max(...M.map(Math.abs), 1);
                    for (let i = 0; i < numEigenvalues; ++i) {
                        const lambda = eigenvalues[i];
                        // det(A - lambda*I) = 0 to within the conditioning of
                        // the determinant: its derivative is bounded by
                        // n*scale^(n-1), so use that as the scale factor.
                        const value = Math.abs(charPolyValue(n, M, lambda));
                        const grad = n * Math.pow(
                            Math.max(scale, Math.abs(lambda)) + 1, n - 1);
                        expect(value).toBeLessThanOrEqual(1e-7 * grad);
                    }
                }, 100);
        });

    it('regression (#42): a decoupled trailing 1x1 block is reported', () => {
        // For an upper-triangular matrix every subdiagonal entry is zero, so
        // the matrix is already quasi-triangular and every diagonal entry is a
        // 1x1 block. Upstream's packing loop runs only to mSizeM1 and drops
        // A(N-1,N-1); the port reports it.
        check(fc.tuple(fc.integer({ min: 3, max: 7 }), fc.nat()).chain(
            ([n]) => fc.tuple(fc.constant(n),
                fc.uniqueArray(fc.integer({ min: -9, max: 9 }),
                    { minLength: n, maxLength: n }),
                fc.array(scaled(-5, 5), { minLength: 49, maxLength: 49 })))
            .map(([n, diag, upper]) => ({ n, diag, upper })),
            ({ n, diag, upper }) => {
                const M = new Array<number>(n * n).fill(0);
                for (let r = 0; r < n; ++r) {
                    M[r * n + r] = diag[r];
                    for (let c = r + 1; c < n; ++c) {
                        M[r * n + c] = upper[r * n + c];
                    }
                }
                const solver = new UnsymmetricEigenvalues(n, 1024);
                const iterations = solver.solve(M, +1);
                // Nothing to do: the matrix is already quasi-triangular.
                expect(iterations).toBe(0);
                const { numEigenvalues, eigenvalues } = solver.getEigenvalues();
                expect(numEigenvalues).toBe(n);
                const sorted = diag.slice().sort((x, y) => x - y);
                for (let i = 0; i < n; ++i) {
                    expect(eigenvalues[i]).toBe(sorted[i]);
                }
                // In particular the last diagonal entry, the one upstream
                // drops, is present.
                expect(eigenvalues).toContain(diag[n - 1]);
            }, 100);
    });

    it('agrees with the symmetric 3x3 eigensolver on symmetric input', () => {
        // The diagonal entries are kept away from zero: see the
        // non-convergence property below for what happens when they are not.
        check(fc.tuple(nonzero(-8, 8, 0.25), wellScaled(-8, 8),
            wellScaled(-8, 8), nonzero(-8, 8, 0.25), wellScaled(-8, 8),
            nonzero(-8, 8, 0.25)),
            ([a00, a01, a02, a11, a12, a22]) => {
                const M = [a00, a01, a02, a01, a11, a12, a02, a12, a22];
                const solver = new UnsymmetricEigenvalues(3, 8192);
                solver.solve(M, +1);
                const { numEigenvalues, eigenvalues } = solver.getEigenvalues();
                const ref = new SymmetricEigensolver3x3().solve(a00, a01, a02,
                    a11, a12, a22, false, +1);
                const scale = Math.max(...M.map(Math.abs), 1);
                // A real symmetric matrix has three real eigenvalues; the QR
                // iteration finds them unless it exhausts its iterations.
                expect(numEigenvalues).toBe(3);
                for (let i = 0; i < 3; ++i) {
                    expectClose(eigenvalues[i], ref.evals[i], 1e-8 * scale,
                        1e-8);
                }
            }, 100);
    });

    it('the sort types permute one set of eigenvalues', () => {
        check(fc.array(scaled(-4, 4), { minLength: 16, maxLength: 16 }), M => {
            const results = [-1, 0, 1].map(sortType => {
                const solver = new UnsymmetricEigenvalues(4, 8192);
                solver.solve(M, sortType);
                return solver.getEigenvalues();
            });
            const [dec, none, inc] = results;
            expect(inc.numEigenvalues).toBe(none.numEigenvalues);
            expect(dec.numEigenvalues).toBe(none.numEigenvalues);
            const k = none.numEigenvalues;
            for (let i = 0; i + 1 < k; ++i) {
                expect(inc.eigenvalues[i])
                    .toBeLessThanOrEqual(inc.eigenvalues[i + 1]);
                expect(dec.eigenvalues[i])
                    .toBeGreaterThanOrEqual(dec.eigenvalues[i + 1]);
            }
            const sortedNone = none.eigenvalues.slice().sort((x, y) => x - y);
            for (let i = 0; i < k; ++i) {
                expect(inc.eigenvalues[i] + 0).toBe(sortedNone[i] + 0);
                expect(dec.eigenvalues[i] + 0).toBe(sortedNone[k - 1 - i] + 0);
            }
        }, 100);
    });

    it('upstream: the Francis iteration cycles for a zero-diagonal symmetric '
        + 'matrix and reports no eigenvalues', () => {
            // A = [[0,8,0],[8,0,8],[0,8,0]] is symmetric with characteristic
            // polynomial -L*(L^2 - 128), so its eigenvalues are the well
            // separated reals 0 and +-8*sqrt(2). Upstream's FrancisQRStep has
            // no exceptional shift, and the double shift taken from a
            // trailing 2x2 block of trace zero makes the iteration cycle: the
            // block never decouples, Solve consumes every iteration and the
            // packing loop reports nothing. Preserved (the port matches
            // upstream); pinned here so a future change is deliberate.
            const M = [0, 8, 0, 8, 0, 8, 0, 8, 0];
            for (const maxIterations of [16, 256, 4096]) {
                const solver = new UnsymmetricEigenvalues(3, maxIterations);
                expect(solver.solve(M, +1)).toBe(maxIterations);
                expect(solver.getEigenvalues())
                    .toEqual({ numEigenvalues: 0, eigenvalues: [] });
            }
            // The true eigenvalues, for the record.
            const s2 = Math.sqrt(128);
            for (const lambda of [-s2, 0, s2]) {
                expect(Math.abs(charPolyValue(3, M, lambda)))
                    .toBeLessThanOrEqual(1e-11);
            }
        });

    it('getEigenvalues returns a copy that does not alias solver state', () => {
        const solver = new UnsymmetricEigenvalues(3, 512);
        solver.solve([2, 1, 0, 0, 3, 1, 0, 0, 5], +1);
        const first = solver.getEigenvalues();
        expect(first.eigenvalues).toEqual([2, 3, 5]);
        first.eigenvalues[0] = 12345;
        const second = solver.getEigenvalues();
        expect(second.eigenvalues).toEqual([2, 3, 5]);
        expect(second.eigenvalues).not.toBe(first.eigenvalues);
    });

    it('rejects invalid constructor arguments the way upstream does', () => {
        check(fc.tuple(fc.integer({ min: -3, max: 2 }),
            fc.integer({ min: 0, max: 4 })), ([size, maxIterations]) => {
                // Upstream keeps mSize = 0 unless size >= 3 and
                // maxIterations > 0; Solve then returns 0 and reports no
                // eigenvalues.
                const solver = new UnsymmetricEigenvalues(size, maxIterations);
                expect(solver.solve([1, 0, 0, 0, 1, 0, 0, 0, 1], +1)).toBe(0);
                expect(solver.getEigenvalues())
                    .toEqual({ numEigenvalues: 0, eigenvalues: [] });
            });
    });
});
