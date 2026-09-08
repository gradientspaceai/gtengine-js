import { describe, it, expect } from 'vitest';
import { BandedMatrix } from '../src/BandedMatrix.js';
import { GaussianElimination } from '../src/GaussianElimination.js';
import { check, expectClose, fc, scaled } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (group V38).
// ---------------------------------------------------------------------------
describe('BandedMatrix verification', () => {
    interface Spec {
        n: number;
        numL: number;
        numU: number;
        values: number[];
    }

    // n in [2,6], numL and numU in [0,n-1], plus a pool of entries.
    const specArb: fc.Arbitrary<Spec> = fc.integer({ min: 2, max: 6 })
        .chain(n => fc.record({
            n: fc.constant(n),
            numL: fc.integer({ min: 0, max: n - 1 }),
            numU: fc.integer({ min: 0, max: n - 1 }),
            values: fc.array(scaled(-4, 4), { minLength: 36, maxLength: 36 })
        }));

    // A banded matrix is stored iff the entry lies on the diagonal, on one of
    // the first numU superdiagonals or on one of the first numL subdiagonals.
    // (Each band array has exactly the length of its diagonal, so a band that
    // exists is stored in full.)
    function isStored(spec: Spec, r: number, c: number): boolean {
        const band = c - r;
        return band === 0 || (band > 0 && band <= spec.numU)
            || (band < 0 && -band <= spec.numL);
    }

    // Fill from the pool, then make the matrix strictly diagonally dominant so
    // that it is invertible (and, when symmetric, positive definite).
    function build(spec: Spec, symmetric: boolean): BandedMatrix {
        const m = new BandedMatrix(spec.n, spec.numL, spec.numU);
        let k = 0;
        for (let r = 0; r < spec.n; ++r) {
            for (let c = 0; c < spec.n; ++c) {
                if (r !== c && isStored(spec, r, c)) {
                    const v = spec.values[k++ % spec.values.length];
                    m.set(r, c, symmetric && c < r ? m.get(c, r) : v);
                }
            }
        }
        for (let r = 0; r < spec.n; ++r) {
            let sum = 0;
            for (let c = 0; c < spec.n; ++c) {
                if (c !== r) { sum += Math.abs(m.get(r, c)); }
            }
            m.set(r, r, sum + 1);
        }
        return m;
    }

    function dense(m: BandedMatrix): number[] {
        const n = m.getSize();
        const a = new Array<number>(n * n).fill(0);
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) { a[r * n + c] = m.get(r, c); }
        }
        return a;
    }

    it('stores exactly the diagonal and the requested bands', () => {
        check(specArb, spec => {
            const m = new BandedMatrix(spec.n, spec.numL, spec.numU);
            // The band arrays have the lengths upstream documents.
            expect(m.getDBand().length).toBe(spec.n);
            expect(m.getLBands().length).toBe(spec.numL);
            expect(m.getUBands().length).toBe(spec.numU);
            m.getLBands().forEach((band, i) =>
                expect(band.length).toBe(spec.n - 1 - i));
            m.getUBands().forEach((band, i) =>
                expect(band.length).toBe(spec.n - 1 - i));

            let k = 0;
            for (let r = 0; r < spec.n; ++r) {
                for (let c = 0; c < spec.n; ++c) {
                    const v = spec.values[k++ % spec.values.length] + 7;
                    m.set(r, c, v);
                    if (isStored(spec, r, c)) {
                        expect(m.get(r, c)).toBe(v);
                    } else {
                        // Upstream hands out a reference to a scratch member
                        // that is reset to zero on the next call, so the write
                        // is discarded and the read is zero.
                        expect(m.get(r, c)).toBe(0);
                    }
                }
            }
        });
    });

    it('discards reads and writes outside the matrix', () => {
        check(fc.tuple(specArb, fc.integer({ min: -3, max: 9 }),
            fc.integer({ min: -3, max: 9 })), ([spec, r, c]) => {
                if (0 <= r && r < spec.n && 0 <= c && c < spec.n) { return; }
                const m = build(spec, false);
                const before = dense(m);
                m.set(r, c, 12345);
                expect(m.get(r, c)).toBe(0);
                expect(dense(m)).toEqual(before);
            });
    });

    it('clone is a deep copy', () => {
        check(specArb, spec => {
            const m = build(spec, false);
            const copy = m.clone();
            expect(dense(copy)).toEqual(dense(m));
            copy.set(0, 0, 999);
            expect(m.get(0, 0)).not.toBe(999);
            expect(copy.getDBand()).not.toBe(m.getDBand());
        });
    });

    it('choleskyFactor reproduces A as L*L^T for symmetric dominant matrices',
        () => {
            check(specArb.filter(s => s.numL === s.numU), spec => {
                const m = build(spec, true);
                const original = dense(m);
                expect(m.choleskyFactor()).toBe(true);
                const n = spec.n;
                // L is the lower-triangular part of the factored matrix and
                // the upper-triangular part is L^T.
                for (let r = 0; r < n; ++r) {
                    for (let c = 0; c < n; ++c) {
                        let sum = 0;
                        for (let k = 0; k <= Math.min(r, c); ++k) {
                            sum += m.get(r, k) * m.get(c, k);
                        }
                        expectClose(sum, original[r * n + c], 1e-10, 1e-10);
                        if (r < c) {
                            expect(m.get(r, c)).toBe(m.get(c, r));
                        }
                    }
                }
            });
        });

    it('solveSystem agrees with dense Gaussian elimination', () => {
        check(fc.tuple(specArb.filter(s => s.numL === s.numU),
            fc.array(scaled(-5, 5), { minLength: 6, maxLength: 6 })),
            ([spec, rhs]) => {
                const m = build(spec, true);
                const a = dense(m);
                const b = rhs.slice(0, spec.n);
                const x = b.slice();
                expect(m.solveSystem(x)).toBe(true);
                const ref = new GaussianElimination().compute(spec.n, a,
                    { B: b });
                expect(ref.invertible).toBe(true);
                for (let i = 0; i < spec.n; ++i) {
                    expectClose(x[i], (ref.X as number[])[i], 1e-9, 1e-9);
                }
                // And A*x = b directly.
                for (let r = 0; r < spec.n; ++r) {
                    let sum = 0;
                    for (let c = 0; c < spec.n; ++c) {
                        sum += a[r * spec.n + c] * x[c];
                    }
                    expectClose(sum, b[r], 1e-9, 1e-9);
                }
            });
    });

    it('solveSystemMatrix solves each column, in either storage order', () => {
        check(fc.tuple(specArb.filter(s => s.numL === s.numU),
            fc.integer({ min: 1, max: 3 }),
            fc.array(scaled(-5, 5), { minLength: 18, maxLength: 18 })),
            ([spec, numCols, pool]) => {
                const n = spec.n;
                const columns: number[][] = [];
                for (let j = 0; j < numCols; ++j) {
                    columns.push(pool.slice(j * n, j * n + n));
                }

                // Per-column reference through the vector solver.
                const expected = columns.map(col => {
                    const x = col.slice();
                    expect(build(spec, true).solveSystem(x)).toBe(true);
                    return x;
                });

                const rowMajorB: number[] = [];
                for (let r = 0; r < n; ++r) {
                    for (let j = 0; j < numCols; ++j) {
                        rowMajorB.push(columns[j][r]);
                    }
                }
                const colMajorB: number[] = [];
                for (let j = 0; j < numCols; ++j) {
                    for (let r = 0; r < n; ++r) {
                        colMajorB.push(columns[j][r]);
                    }
                }

                expect(build(spec, true)
                    .solveSystemMatrix(rowMajorB, numCols, true)).toBe(true);
                expect(build(spec, true)
                    .solveSystemMatrix(colMajorB, numCols, false)).toBe(true);
                for (let r = 0; r < n; ++r) {
                    for (let j = 0; j < numCols; ++j) {
                        expectClose(rowMajorB[r * numCols + j], expected[j][r],
                            1e-9, 1e-9);
                        expectClose(colMajorB[j * n + r], expected[j][r],
                            1e-9, 1e-9);
                    }
                }
            });
    });

    it('computeInverse produces the two-sided inverse and leaves A unchanged',
        () => {
            check(specArb, spec => {
                const m = build(spec, false);
                const a = dense(m);
                const n = spec.n;
                const inv = new Array<number>(n * n).fill(0);
                expect(m.computeInverse(inv)).toBe(true);
                // computeInverse is const upstream.
                expect(dense(m)).toEqual(a);
                for (let r = 0; r < n; ++r) {
                    for (let c = 0; c < n; ++c) {
                        let left = 0;
                        let right = 0;
                        for (let k = 0; k < n; ++k) {
                            left += a[r * n + k] * inv[k * n + c];
                            right += inv[r * n + k] * a[k * n + c];
                        }
                        const target = (r === c ? 1 : 0);
                        expectClose(left, target, 1e-9, 1e-9);
                        expectClose(right, target, 1e-9, 1e-9);
                    }
                }
            });
        });

    it('computeInverse in column-major order transposes the row-major result',
        () => {
            check(specArb, spec => {
                const n = spec.n;
                const m = build(spec, false);
                const rowMajor = new Array<number>(n * n).fill(0);
                const colMajor = new Array<number>(n * n).fill(0);
                expect(m.computeInverse(rowMajor, true)).toBe(true);
                expect(m.computeInverse(colMajor, false)).toBe(true);
                for (let r = 0; r < n; ++r) {
                    for (let c = 0; c < n; ++c) {
                        expectClose(rowMajor[r * n + c], colMajor[c * n + r],
                            1e-12, 1e-12);
                    }
                }
            });
        });

    it('an invalid construction produces an inert matrix of size 0', () => {
        check(fc.tuple(fc.integer({ min: -2, max: 4 }),
            fc.integer({ min: -2, max: 6 }), fc.integer({ min: -2, max: 6 })),
            ([size, numL, numU]) => {
                const valid = size > 0 && 0 <= numL && numL < size
                    && 0 <= numU && numU < size;
                const m = new BandedMatrix(size, numL, numU);
                if (valid) {
                    expect(m.getSize()).toBe(size);
                    return;
                }
                expect(m.getSize()).toBe(0);
                expect(m.getDBand()).toEqual([]);
                m.set(0, 0, 5);
                expect(m.get(0, 0)).toBe(0);
                // Upstream returns false when mDBand is empty.
                expect(m.choleskyFactor()).toBe(false);
                expect(m.solveSystem([])).toBe(false);
            });
    });

    it('choleskyFactor fails when the band counts differ or A is not positive '
        + 'definite', () => {
            check(specArb.filter(s => s.numL !== s.numU), spec => {
                // Upstream requires L and U to have the same number of bands.
                expect(build(spec, false).choleskyFactor()).toBe(false);
            });
            check(specArb.filter(s => s.numL === s.numU && s.numL > 0), spec => {
                const m = build(spec, true);
                // Negate the diagonal: the matrix is now negative definite and
                // the very first pivot is nonpositive.
                for (let r = 0; r < spec.n; ++r) {
                    m.set(r, r, -m.get(r, r));
                }
                expect(m.choleskyFactor()).toBe(false);
            });
        });
});
