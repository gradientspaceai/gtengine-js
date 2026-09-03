import { describe, it, expect } from 'vitest';
import { SymmetricEigensolver } from '../src/SymmetricEigensolver.js';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Element (row, col) of an NxN row-major matrix.
function at(n: number, m: readonly number[], row: number, col: number): number {
    return m[col + n * row];
}

// The i-th eigenvector is column i of the eigenvector matrix.
function column(n: number, m: readonly number[], col: number): number[] {
    const v = new Array<number>(n).fill(0);
    for (let row = 0; row < n; ++row) {
        v[row] = m[col + n * row];
    }
    return v;
}

// Q * D * Q^T where D is the diagonal matrix of the given values.
function reconstruct(n: number, q: readonly number[],
    d: readonly number[]): number[] {
    const a = new Array<number>(n * n).fill(0);
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            let sum = 0;
            for (let k = 0; k < n; ++k) {
                sum += at(n, q, r, k) * d[k] * at(n, q, c, k);
            }
            a[c + n * r] = sum;
        }
    }
    return a;
}

// The largest absolute difference between two equal-length arrays.
function maxDiff(a: readonly number[], b: readonly number[]): number {
    let m = 0;
    for (let i = 0; i < a.length; ++i) {
        m = Math.max(m, Math.abs(a[i] - b[i]));
    }
    return m;
}

// Q^T * Q compared against the identity.
function orthogonalityError(n: number, q: readonly number[]): number {
    let m = 0;
    for (let i = 0; i < n; ++i) {
        for (let j = 0; j < n; ++j) {
            let sum = 0;
            for (let r = 0; r < n; ++r) {
                sum += at(n, q, r, i) * at(n, q, r, j);
            }
            m = Math.max(m, Math.abs(sum - (i === j ? 1 : 0)));
        }
    }
    return m;
}

// The determinant by Gaussian elimination with partial pivoting.
function determinant(n: number, m: readonly number[]): number {
    const a = m.slice();
    let det = 1;
    for (let col = 0; col < n; ++col) {
        let pivot = col;
        for (let row = col + 1; row < n; ++row) {
            if (Math.abs(a[col + n * row]) > Math.abs(a[col + n * pivot])) {
                pivot = row;
            }
        }
        if (a[col + n * pivot] === 0) {
            return 0;
        }
        if (pivot !== col) {
            for (let c = 0; c < n; ++c) {
                const t = a[c + n * col];
                a[c + n * col] = a[c + n * pivot];
                a[c + n * pivot] = t;
            }
            det = -det;
        }
        det *= a[col + n * col];
        for (let row = col + 1; row < n; ++row) {
            const factor = a[col + n * row] / a[col + n * col];
            for (let c = col; c < n; ++c) {
                a[c + n * row] -= factor * a[c + n * col];
            }
        }
    }
    return det;
}

// A random symmetric NxN matrix with entries in [-1,1].
function randomSymmetric(n: number, rand: () => number): number[] {
    const a = new Array<number>(n * n).fill(0);
    for (let r = 0; r < n; ++r) {
        for (let c = r; c < n; ++c) {
            const value = 2 * rand() - 1;
            a[c + n * r] = value;
            a[r + n * c] = value;
        }
    }
    return a;
}

describe('SymmetricEigensolver', () => {
    it('is inert when the size is too small', () => {
        for (const size of [1, 0, -3]) {
            const solver = new SymmetricEigensolver(size, 32);
            expect(solver.solve([1], 0)).toBe(0);
            expect(solver.getEigenvalues()).toEqual([]);
            expect(solver.getEigenvectors()).toEqual([]);
            expect(solver.getEigenvector(0)).toEqual([]);
            expect(solver.getEigenvalue(0)).toBe(Number.MAX_VALUE);
            expect(solver.getEigenvectorMatrixType()).toBe(-1);
        }
    });

    it('is inert when the iteration budget is not positive', () => {
        const solver = new SymmetricEigensolver(3, 0);
        expect(solver.solve([1, 0, 0, 0, 1, 0, 0, 0, 1], 0)).toBe(0);
    });

    it('converges immediately for a diagonal matrix', () => {
        // The tridiagonalization of a diagonal matrix leaves the
        // superdiagonal zero, so the very first convergence test succeeds.
        const a = [3, 0, 0, 0, -1, 0, 0, 0, 5];
        const solver = new SymmetricEigensolver(3, 32);
        expect(solver.solve(a, 0)).toBe(0);
        expect(solver.getEigenvalues()).toEqual([3, -1, 5]);

        expect(solver.solve(a, +1)).toBe(0);
        expect(solver.getEigenvalues()).toEqual([-1, 3, 5]);

        expect(solver.solve(a, -1)).toBe(0);
        expect(solver.getEigenvalues()).toEqual([5, 3, -1]);
    });

    it('handles the zero matrix', () => {
        const n = 4;
        const solver = new SymmetricEigensolver(n, 32);
        expect(solver.solve(new Array<number>(n * n).fill(0), +1)).toBe(0);
        expect(solver.getEigenvalues()).toEqual([0, 0, 0, 0]);
        const q = solver.getEigenvectors();
        expect(orthogonalityError(n, q)).toBeLessThan(1e-14);
    });

    it('computes the known decomposition of a 2x2 matrix', () => {
        // A = {{2,1},{1,2}} has eigenvalues 1 and 3 with eigenvectors
        // (1,-1)/sqrt(2) and (1,1)/sqrt(2).
        const n = 2;
        const a = [2, 1, 1, 2];
        const solver = new SymmetricEigensolver(n, 32);
        const iterations = solver.solve(a, +1);
        expect(iterations).toBeGreaterThan(0);
        expect(iterations).toBeLessThan(32);

        const values = solver.getEigenvalues();
        expect(values[0]).toBeCloseTo(1, 14);
        expect(values[1]).toBeCloseTo(3, 14);

        const q = solver.getEigenvectors();
        const v0 = column(n, q, 0);
        const v1 = column(n, q, 1);
        const invSqrt2 = 1 / Math.sqrt(2);
        expect(Math.abs(v0[0])).toBeCloseTo(invSqrt2, 13);
        expect(v0[0] + v0[1]).toBeCloseTo(0, 13);
        expect(Math.abs(v1[0])).toBeCloseTo(invSqrt2, 13);
        expect(v1[0] - v1[1]).toBeCloseTo(0, 13);
    });

    it('computes the known decomposition of the 3x3 second-difference matrix', () => {
        // A = tridiag(-1,2,-1) of size 3 has eigenvalues 2 - sqrt(2), 2 and
        // 2 + sqrt(2).
        const n = 3;
        const a = [
            2, -1, 0,
            -1, 2, -1,
            0, -1, 2
        ];
        const solver = new SymmetricEigensolver(n, 64);
        solver.solve(a, +1);
        const values = solver.getEigenvalues();
        expect(values[0]).toBeCloseTo(2 - Math.SQRT2, 13);
        expect(values[1]).toBeCloseTo(2, 13);
        expect(values[2]).toBeCloseTo(2 + Math.SQRT2, 13);

        // The eigenvector for the middle eigenvalue is (1,0,-1)/sqrt(2).
        const q = solver.getEigenvectors();
        const v = column(n, q, 1);
        expect(Math.abs(v[0])).toBeCloseTo(1 / Math.SQRT2, 12);
        expect(v[1]).toBeCloseTo(0, 12);
        expect(v[0] + v[2]).toBeCloseTo(0, 12);
    });

    it('satisfies A*v = lambda*v for every eigenpair', () => {
        const n = 5;
        const a = [
            4, 1, -2, 0, 3,
            1, 2, 0, 1, -1,
            -2, 0, 3, 2, 0,
            0, 1, 2, -5, 1,
            3, -1, 0, 1, 6
        ];
        const solver = new SymmetricEigensolver(n, 64);
        solver.solve(a, -1);
        const values = solver.getEigenvalues();
        const q = solver.getEigenvectors();

        for (let i = 0; i < n; ++i) {
            const v = column(n, q, i);
            for (let r = 0; r < n; ++r) {
                let sum = 0;
                for (let c = 0; c < n; ++c) {
                    sum += at(n, a, r, c) * v[c];
                }
                expect(sum).toBeCloseTo(values[i] * v[r], 11);
            }
            // getEigenvector(i) computes the same column incrementally.
            expect(maxDiff(solver.getEigenvector(i), v)).toBeLessThan(1e-12);
            expect(solver.getEigenvalue(i)).toBe(values[i]);
        }

        // The eigenvalues are sorted in decreasing order.
        for (let i = 1; i < n; ++i) {
            expect(values[i - 1]).toBeGreaterThanOrEqual(values[i]);
        }
    });

    it('preserves the trace and the determinant', () => {
        const n = 4;
        const a = [
            1, 2, 3, 4,
            2, 5, 6, 7,
            3, 6, 8, 9,
            4, 7, 9, 10
        ];
        const solver = new SymmetricEigensolver(n, 64);
        solver.solve(a, 0);
        const values = solver.getEigenvalues();

        let trace = 0;
        for (let i = 0; i < n; ++i) {
            trace += at(n, a, i, i);
        }
        const sum = values.reduce((s, v) => s + v, 0);
        expect(sum).toBeCloseTo(trace, 11);

        const product = values.reduce((s, v) => s * v, 1);
        expect(product).toBeCloseTo(determinant(n, a), 9);
    });

    it('handles repeated eigenvalues', () => {
        // A = 2*I + u*u^T with u = (1,1,1,1) has eigenvalue 6 (multiplicity
        // 1) and eigenvalue 2 (multiplicity 3).
        const n = 4;
        const a = new Array<number>(n * n).fill(0);
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                a[c + n * r] = (r === c ? 2 : 0) + 1;
            }
        }

        const solver = new SymmetricEigensolver(n, 64);
        solver.solve(a, +1);
        const values = solver.getEigenvalues();
        expect(values[0]).toBeCloseTo(2, 12);
        expect(values[1]).toBeCloseTo(2, 12);
        expect(values[2]).toBeCloseTo(2, 12);
        expect(values[3]).toBeCloseTo(6, 12);

        const q = solver.getEigenvectors();
        expect(orthogonalityError(n, q)).toBeLessThan(1e-13);
        expect(maxDiff(reconstruct(n, q, values), a)).toBeLessThan(1e-12);
    });

    it('handles a matrix that decouples into two blocks', () => {
        // The block-diagonal matrix diag({{2,1},{1,2}}, {{5,2},{2,5}}) has
        // eigenvalues 1, 3, 3 and 7. The tridiagonal reduction has a zero
        // superdiagonal entry between the blocks, exercising the code that
        // isolates the lower-right unreduced block.
        const n = 4;
        const a = [
            2, 1, 0, 0,
            1, 2, 0, 0,
            0, 0, 5, 2,
            0, 0, 2, 5
        ];
        const solver = new SymmetricEigensolver(n, 64);
        solver.solve(a, +1);
        const values = solver.getEigenvalues();
        expect(values[0]).toBeCloseTo(1, 12);
        expect(values[1]).toBeCloseTo(3, 12);
        expect(values[2]).toBeCloseTo(3, 12);
        expect(values[3]).toBeCloseTo(7, 12);

        const q = solver.getEigenvectors();
        expect(maxDiff(reconstruct(n, q, values), a)).toBeLessThan(1e-12);
    });

    it('builds a correct Q when a Householder step degenerates (upstream fix)', () => {
        // Upstream Tridiagonalize stores 2/Dot(v,v) = 2 for a Householder
        // step whose subcolumn is already zero even though the identity was
        // applied, so the rebuilt reflection in GetEigenvectors flips the
        // sign of row i+1 of Q. See the comment in src/SymmetricEigensolver.
        // Every matrix below has at least one such degenerate step.
        const cases: Array<{ n: number; a: number[] }> = [
            {
                // Block diagonal 2+2; the step for i = 1 degenerates.
                n: 4,
                a: [
                    2, 1, 0, 0,
                    1, 2, 0, 0,
                    0, 0, 5, 2,
                    0, 0, 2, 5
                ]
            },
            {
                // Already tridiagonal, so every Householder step degenerates.
                n: 5,
                a: [
                    1, 2, 0, 0, 0,
                    2, 3, 1, 0, 0,
                    0, 1, -4, 5, 0,
                    0, 0, 5, 2, 1,
                    0, 0, 0, 1, 7
                ]
            },
            {
                // Block diagonal 1+4 with a leading isolated eigenvalue.
                n: 5,
                a: [
                    8, 0, 0, 0, 0,
                    0, 1, 2, 3, 4,
                    0, 2, 5, 6, 7,
                    0, 3, 6, 9, 1,
                    0, 4, 7, 1, 2
                ]
            }
        ];

        for (const { n, a } of cases) {
            for (const sortType of [0, +1, -1]) {
                const solver = new SymmetricEigensolver(n, 256);
                solver.solve(a, sortType);
                const values = solver.getEigenvalues();
                const q = solver.getEigenvectors();
                expect(orthogonalityError(n, q)).toBeLessThan(1e-12);
                expect(maxDiff(reconstruct(n, q, values), a)).toBeLessThan(1e-11);
                for (let i = 0; i < n; ++i) {
                    expect(maxDiff(solver.getEigenvector(i), column(n, q, i)))
                        .toBeLessThan(1e-12);
                }
            }
        }
    });

    it('does not drop the trailing eigenvalue when the last block decouples', () => {
        // A tridiagonal matrix whose last superdiagonal entry is already
        // zero, so the trailing 1x1 block is decoupled from the start. All N
        // eigenvalues must still be reported.
        const n = 5;
        const a = [
            4, 1, 0, 0, 0,
            1, 4, 1, 0, 0,
            0, 1, 4, 1, 0,
            0, 0, 1, 4, 0,
            0, 0, 0, 0, 9
        ];
        const solver = new SymmetricEigensolver(n, 64);
        solver.solve(a, +1);
        const values = solver.getEigenvalues();
        expect(values.length).toBe(n);

        // The leading 4x4 block is 4*I + tridiag(1,0,1) whose eigenvalues are
        // 4 + 2*cos(k*pi/5) for k = 1..4; the trailing block contributes 9.
        const expected = [9];
        for (let k = 1; k <= 4; ++k) {
            expected.push(4 + 2 * Math.cos(k * Math.PI / 5));
        }
        expected.sort((p, q2) => p - q2);
        for (let i = 0; i < n; ++i) {
            expect(values[i]).toBeCloseTo(expected[i], 12);
        }

        const q = solver.getEigenvectors();
        expect(orthogonalityError(n, q)).toBeLessThan(1e-13);
        expect(maxDiff(reconstruct(n, q, values), a)).toBeLessThan(1e-12);
    });

    it('sorts eigenvalues and eigenvectors consistently', () => {
        const n = 5;
        const rand = makeRandom(20260830);
        const a = randomSymmetric(n, rand);

        const unsorted = new SymmetricEigensolver(n, 64);
        unsorted.solve(a, 0);
        const valuesU = unsorted.getEigenvalues();
        const vectorsU = unsorted.getEigenvectors();

        const increasing = new SymmetricEigensolver(n, 64);
        increasing.solve(a, +1);
        const valuesI = increasing.getEigenvalues();
        const vectorsI = increasing.getEigenvectors();

        const decreasing = new SymmetricEigensolver(n, 64);
        decreasing.solve(a, -1);
        const valuesD = decreasing.getEigenvalues();

        expect(valuesI.slice().sort((p, q) => p - q)).toEqual(valuesI);
        expect(valuesD.slice().sort((p, q) => q - p)).toEqual(valuesD);
        expect(valuesD.slice().reverse()).toEqual(valuesI);
        expect(valuesU.slice().sort((p, q) => p - q)).toEqual(valuesI);

        // The sorted eigenvector matrix is the unsorted one with permuted
        // columns; both reconstruct A.
        expect(maxDiff(reconstruct(n, vectorsU, valuesU), a)).toBeLessThan(1e-13);
        expect(maxDiff(reconstruct(n, vectorsI, valuesI), a)).toBeLessThan(1e-13);
        for (let i = 0; i < n; ++i) {
            const j = valuesU.indexOf(valuesI[i]);
            expect(maxDiff(column(n, vectorsI, i), column(n, vectorsU, j)))
                .toBeLessThan(1e-15);
        }
    });

    it('reports whether the eigenvector matrix is a rotation or a reflection', () => {
        const rand = makeRandom(777);
        for (let n = 2; n <= 6; ++n) {
            for (const sortType of [0, +1, -1]) {
                const a = randomSymmetric(n, rand);
                const solver = new SymmetricEigensolver(n, 128);
                solver.solve(a, sortType);
                // The type is only final after getEigenvectors().
                const q = solver.getEigenvectors();
                const type = solver.getEigenvectorMatrixType();
                expect(type === 0 || type === 1).toBe(true);
                const det = determinant(n, q);
                expect(det).toBeCloseTo(type === 1 ? 1 : -1, 10);
            }
        }
    });

    it('reconstructs random symmetric matrices of several sizes', () => {
        const rand = makeRandom(13579);
        for (let n = 2; n <= 10; ++n) {
            for (let trial = 0; trial < 5; ++trial) {
                const a = randomSymmetric(n, rand);
                const solver = new SymmetricEigensolver(n, 4096);
                const iterations = solver.solve(a, +1);
                expect(iterations).not.toBe(SymmetricEigensolver.noConvergence);

                const values = solver.getEigenvalues();
                const q = solver.getEigenvectors();
                expect(orthogonalityError(n, q)).toBeLessThan(1e-12);
                expect(maxDiff(reconstruct(n, q, values), a)).toBeLessThan(1e-12);

                // Each single-eigenvector query matches the full matrix.
                for (let i = 0; i < n; ++i) {
                    expect(maxDiff(solver.getEigenvector(i), column(n, q, i)))
                        .toBeLessThan(1e-12);
                }
            }
        }
    });

    it('reports failure to converge within the iteration budget', () => {
        const n = 6;
        const a = randomSymmetric(n, makeRandom(2468));
        const solver = new SymmetricEigensolver(n, 1);
        expect(solver.solve(a, +1)).toBe(SymmetricEigensolver.noConvergence);
    });

    it('can be reused for several matrices', () => {
        const n = 3;
        const solver = new SymmetricEigensolver(n, 64);

        const a0 = [2, 0, 0, 0, 4, 0, 0, 0, 6];
        solver.solve(a0, +1);
        expect(solver.getEigenvalues()[0]).toBeCloseTo(2, 14);

        const a1 = [
            2, -1, 0,
            -1, 2, -1,
            0, -1, 2
        ];
        solver.solve(a1, +1);
        const values = solver.getEigenvalues();
        expect(values[0]).toBeCloseTo(2 - Math.SQRT2, 13);
        expect(values[2]).toBeCloseTo(2 + Math.SQRT2, 13);

        solver.solve(a0, +1);
        expect(solver.getEigenvalues()).toEqual([2, 4, 6]);
    });
});
