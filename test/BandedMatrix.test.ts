import { describe, it, expect } from 'vitest';
import { BandedMatrix } from '../src/BandedMatrix';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The dense (row-major) form of a banded matrix.
function densify(A: BandedMatrix): number[] {
    const n = A.getSize();
    const dense = new Array<number>(n * n).fill(0);
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            dense[c + n * r] = A.get(r, c);
        }
    }
    return dense;
}

// A symmetric positive definite banded matrix with the given bandwidth. The
// diagonal is made dominant so that the Cholesky factorization succeeds.
function makeSPD(n: number, numBands: number, rand: () => number): BandedMatrix {
    const A = new BandedMatrix(n, numBands, numBands);
    for (let r = 0; r < n; ++r) {
        A.set(r, r, n + 2 * rand());
        for (let b = 1; b <= numBands && r + b < n; ++b) {
            const value = 2 * rand() - 1;
            A.set(r, r + b, value);
            A.set(r + b, r, value);
        }
    }
    return A;
}

function multiplyDense(n: number, A: readonly number[],
    B: readonly number[], k: number): number[] {
    const C = new Array<number>(n * k).fill(0);
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < k; ++c) {
            let sum = 0;
            for (let i = 0; i < n; ++i) {
                sum += A[i + n * r] * B[c + k * i];
            }
            C[c + k * r] = sum;
        }
    }
    return C;
}

describe('BandedMatrix construction and element access', () => {
    it('allocates bands of decreasing length', () => {
        const A = new BandedMatrix(5, 2, 3);
        expect(A.getSize()).toBe(5);
        expect(A.getDBand().length).toBe(5);
        expect(A.getLBands().map(b => b.length)).toEqual([4, 3]);
        expect(A.getUBands().map(b => b.length)).toEqual([4, 3, 2]);
    });

    it('starts zero-filled', () => {
        const A = new BandedMatrix(4, 1, 1);
        expect(densify(A)).toEqual(new Array<number>(16).fill(0));
    });

    it('rejects invalid sizes by producing an empty matrix', () => {
        for (const args of [[0, 0, 0], [-1, 0, 0], [3, 3, 0], [3, 0, 3],
            [3, -1, 0], [3, 0, -1]] as Array<[number, number, number]>) {
            const A = new BandedMatrix(args[0], args[1], args[2]);
            expect(A.getSize()).toBe(0);
            expect(A.getDBand().length).toBe(0);
        }
    });

    it('stores entries in the correct band', () => {
        const A = new BandedMatrix(4, 1, 2);
        A.set(0, 0, 10);
        A.set(1, 1, 11);
        A.set(1, 0, 20);  // first lower band
        A.set(0, 1, 30);  // first upper band
        A.set(0, 2, 40);  // second upper band

        expect(A.getDBand()[0]).toBe(10);
        expect(A.getDBand()[1]).toBe(11);
        expect(A.getLBands()[0][0]).toBe(20);
        expect(A.getUBands()[0][0]).toBe(30);
        expect(A.getUBands()[1][0]).toBe(40);

        expect(A.get(0, 0)).toBe(10);
        expect(A.get(1, 0)).toBe(20);
        expect(A.get(0, 1)).toBe(30);
        expect(A.get(0, 2)).toBe(40);
    });

    it('discards writes outside the bands and reads them as zero', () => {
        const A = new BandedMatrix(4, 1, 1);
        // (0, 2) is in the second upper band, which is not allocated.
        A.set(0, 2, 99);
        expect(A.get(0, 2)).toBe(0);
        // (3, 0) is in the third lower band, which is not allocated.
        A.set(3, 0, 99);
        expect(A.get(3, 0)).toBe(0);
        expect(densify(A)).toEqual(new Array<number>(16).fill(0));
    });

    it('discards writes outside the matrix and reads them as zero', () => {
        const A = new BandedMatrix(3, 1, 1);
        A.set(-1, 0, 5);
        A.set(0, 3, 5);
        A.set(3, 3, 5);
        expect(A.get(-1, 0)).toBe(0);
        expect(A.get(0, 3)).toBe(0);
        expect(A.get(3, 3)).toBe(0);
        expect(densify(A)).toEqual(new Array<number>(9).fill(0));
    });

    it('clones deeply', () => {
        const A = new BandedMatrix(3, 1, 1);
        A.set(0, 0, 1);
        A.set(1, 0, 2);
        A.set(0, 1, 3);
        const B = A.clone();
        B.set(0, 0, 99);
        expect(A.get(0, 0)).toBe(1);
        expect(B.get(1, 0)).toBe(2);
        expect(B.get(0, 1)).toBe(3);
    });
});

describe('BandedMatrix Cholesky factorization', () => {
    it('factors a tridiagonal SPD matrix so that L*L^T is the original', () => {
        const n = 6;
        const A = new BandedMatrix(n, 1, 1);
        for (let r = 0; r < n; ++r) {
            A.set(r, r, 4);
            if (r + 1 < n) {
                A.set(r, r + 1, 1);
                A.set(r + 1, r, 1);
            }
        }
        const original = densify(A);

        expect(A.choleskyFactor()).toBe(true);

        // The lower-triangular part is L; verify L*L^T within the band.
        const n2 = n;
        for (let r = 0; r < n2; ++r) {
            for (let c = 0; c < n2; ++c) {
                let sum = 0;
                for (let k = 0; k <= Math.min(r, c); ++k) {
                    const lrk = r >= k ? A.get(r, k) : 0;
                    const lck = c >= k ? A.get(c, k) : 0;
                    sum += lrk * lck;
                }
                expect(sum).toBeCloseTo(original[c + n2 * r], 10);
            }
        }
    });

    it('leaves the upper-triangular part equal to the transpose of L', () => {
        const rand = makeRandom(11);
        const A = makeSPD(5, 2, rand);
        expect(A.choleskyFactor()).toBe(true);
        for (let r = 0; r < 5; ++r) {
            for (let c = r + 1; c < 5; ++c) {
                expect(A.get(r, c)).toBe(A.get(c, r));
            }
        }
    });

    it('fails when the matrix is not positive definite', () => {
        const A = new BandedMatrix(3, 1, 1);
        for (let r = 0; r < 3; ++r) {
            A.set(r, r, -1);
        }
        expect(A.choleskyFactor()).toBe(false);
    });

    it('fails when the number of lower and upper bands differ', () => {
        const A = new BandedMatrix(4, 1, 2);
        for (let r = 0; r < 4; ++r) {
            A.set(r, r, 4);
        }
        expect(A.choleskyFactor()).toBe(false);
    });

    it('fails for the empty matrix', () => {
        expect(new BandedMatrix(0, 0, 0).choleskyFactor()).toBe(false);
    });
});

describe('BandedMatrix solvers', () => {
    it('solves a tridiagonal system with a known solution', () => {
        const n = 4;
        const A = new BandedMatrix(n, 1, 1);
        for (let r = 0; r < n; ++r) {
            A.set(r, r, 2);
            if (r + 1 < n) {
                A.set(r, r + 1, -1);
                A.set(r + 1, r, -1);
            }
        }
        // A is the 1D Laplacian. With x = (1, 2, 3, 4), b = A*x.
        const x = [1, 2, 3, 4];
        const dense = densify(A);
        const b = multiplyDense(n, dense, x, 1);

        expect(A.solveSystem(b)).toBe(true);
        for (let i = 0; i < n; ++i) {
            expect(b[i]).toBeCloseTo(x[i], 10);
        }
    });

    it('solves a system with a matrix right-hand side (row major)', () => {
        const rand = makeRandom(4321);
        const n = 6;
        const A = makeSPD(n, 2, rand);
        const dense = densify(A);
        const numColumns = 3;
        const X: number[] = [];
        for (let i = 0; i < n * numColumns; ++i) {
            X.push(2 * rand() - 1);
        }
        const B = multiplyDense(n, dense, X, numColumns);

        expect(A.solveSystemMatrix(B, numColumns)).toBe(true);
        for (let i = 0; i < n * numColumns; ++i) {
            expect(B[i]).toBeCloseTo(X[i], 9);
        }
    });

    it('solves a system with a matrix right-hand side (column major)', () => {
        const rand = makeRandom(4321);
        const n = 6;
        const A = makeSPD(n, 2, rand);
        const dense = densify(A);
        const numColumns = 3;
        const X: number[] = [];
        for (let i = 0; i < n * numColumns; ++i) {
            X.push(2 * rand() - 1);
        }
        const B = multiplyDense(n, dense, X, numColumns);

        // Transpose B and X into column-major storage.
        const Bcol = new Array<number>(n * numColumns).fill(0);
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < numColumns; ++c) {
                Bcol[r + n * c] = B[c + numColumns * r];
            }
        }

        expect(A.solveSystemMatrix(Bcol, numColumns, false)).toBe(true);
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < numColumns; ++c) {
                expect(Bcol[r + n * c])
                    .toBeCloseTo(X[c + numColumns * r], 9);
            }
        }
    });

    it('fails to solve a system that is not positive definite', () => {
        const A = new BandedMatrix(3, 1, 1);
        A.set(0, 0, 1);
        A.set(1, 1, -1);
        A.set(2, 2, 1);
        expect(A.solveSystem([1, 1, 1])).toBe(false);
    });

    it('solves random SPD banded systems (randomized cross-check)', () => {
        const rand = makeRandom(24680);
        for (let n = 2; n <= 8; ++n) {
            for (let numBands = 1; numBands < Math.min(n, 4); ++numBands) {
                const A = makeSPD(n, numBands, rand);
                const dense = densify(A);
                const x: number[] = [];
                for (let i = 0; i < n; ++i) {
                    x.push(2 * rand() - 1);
                }
                const b = multiplyDense(n, dense, x, 1);
                expect(A.solveSystem(b)).toBe(true);
                for (let i = 0; i < n; ++i) {
                    expect(b[i]).toBeCloseTo(x[i], 9);
                }
            }
        }
    });
});

describe('BandedMatrix inverse', () => {
    it('inverts a tridiagonal matrix without modifying it', () => {
        const n = 5;
        const A = new BandedMatrix(n, 1, 1);
        for (let r = 0; r < n; ++r) {
            A.set(r, r, 2);
            if (r + 1 < n) {
                A.set(r, r + 1, -1);
                A.set(r + 1, r, -1);
            }
        }
        const dense = densify(A);

        const inverse = new Array<number>(n * n).fill(0);
        expect(A.computeInverse(inverse)).toBe(true);
        expect(densify(A)).toEqual(dense);

        const identity = multiplyDense(n, dense, inverse, n);
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                expect(identity[c + n * r]).toBeCloseTo(r === c ? 1 : 0, 10);
            }
        }
    });

    it('inverts a nonsymmetric banded matrix', () => {
        const n = 5;
        const A = new BandedMatrix(n, 1, 2);
        for (let r = 0; r < n; ++r) {
            A.set(r, r, 5);
            if (r + 1 < n) {
                A.set(r, r + 1, 1);
                A.set(r + 1, r, 2);
            }
            if (r + 2 < n) {
                A.set(r, r + 2, -1);
            }
        }
        const dense = densify(A);
        const inverse = new Array<number>(n * n).fill(0);
        expect(A.computeInverse(inverse)).toBe(true);

        const identity = multiplyDense(n, dense, inverse, n);
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                expect(identity[c + n * r]).toBeCloseTo(r === c ? 1 : 0, 10);
            }
        }
    });

    it('produces the transposed inverse for column-major storage', () => {
        const n = 4;
        const A = new BandedMatrix(n, 1, 1);
        for (let r = 0; r < n; ++r) {
            A.set(r, r, 3);
            if (r + 1 < n) {
                A.set(r, r + 1, 1);
                A.set(r + 1, r, 1);
            }
        }
        const rowMajorInverse = new Array<number>(n * n).fill(0);
        const colMajorInverse = new Array<number>(n * n).fill(0);
        expect(A.computeInverse(rowMajorInverse)).toBe(true);
        expect(A.computeInverse(colMajorInverse, false)).toBe(true);
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                expect(colMajorInverse[r + n * c])
                    .toBeCloseTo(rowMajorInverse[c + n * r], 12);
            }
        }
    });

    it('reports failure for a singular matrix', () => {
        const n = 3;
        const A = new BandedMatrix(n, 1, 1);
        // A zero pivot on the first row.
        A.set(0, 0, 0);
        A.set(1, 1, 1);
        A.set(2, 2, 1);
        const inverse = new Array<number>(n * n).fill(0);
        expect(A.computeInverse(inverse)).toBe(false);
    });
});
