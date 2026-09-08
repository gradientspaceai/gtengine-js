import { describe, it, expect } from 'vitest';
import {
    BlockCholeskyDecomposition, CholeskyDecomposition
} from '../src/CholeskyDecomposition.js';
import type {
    CholeskyBlockMatrix, CholeskyBlockVector
} from '../src/CholeskyDecomposition.js';
import { LinearSystem } from '../src/LinearSystem.js';
import {
    Matrix, determinant, multiplyABT, mulMatrix
} from '../src/Matrix.js';
import { Vector, sub } from '../src/Vector.js';
import {
    check, expectClose, fc, scaled, wellScaledVector
} from './helpers/arbitraries.js';

// A deterministic pseudorandom generator so the randomized cross-checks are
// reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state + 0x6D2B79F5) >>> 0;
        let t = state;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

// A symmetric positive definite matrix A = M*M^T + n*I.
function randomSPD(n: number, rng: () => number): Matrix {
    const M = new Matrix(n, n);
    for (let i = 0; i < n * n; ++i) {
        M.values[i] = 2 * rng() - 1;
    }
    const A = multiplyABT(M, M);
    for (let i = 0; i < n; ++i) {
        A.set(i, i, A.get(i, i) + n);
    }
    return A;
}

// The strict lower triangle plus the diagonal of the factored matrix.
function lowerPart(A: Matrix): Matrix {
    const n = A.numRows;
    const L = new Matrix(n, n);
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c <= r; ++c) {
            L.set(r, c, A.get(r, c));
        }
    }
    return L;
}

function maxAbsDiff(M0: Matrix, M1: Matrix): number {
    let d = 0;
    for (let i = 0; i < M0.numElements; ++i) {
        d = Math.max(d, Math.abs(M0.values[i] - M1.values[i]));
    }
    return d;
}

function maxAbs(v: Vector): number {
    let d = 0;
    for (let i = 0; i < v.size; ++i) {
        d = Math.max(d, Math.abs(v.get(i)));
    }
    return d;
}

// Split a full matrix into the flat row-major block layout used by
// BlockCholeskyDecomposition.
function toBlockMatrix(A: Matrix, blockSize: number, numBlocks: number):
    CholeskyBlockMatrix {
    const blocks: CholeskyBlockMatrix = [];
    for (let r = 0; r < numBlocks; ++r) {
        for (let c = 0; c < numBlocks; ++c) {
            const block = new Matrix(blockSize, blockSize);
            for (let j = 0; j < blockSize; ++j) {
                for (let i = 0; i < blockSize; ++i) {
                    block.set(j, i,
                        A.get(r * blockSize + j, c * blockSize + i));
                }
            }
            blocks.push(block);
        }
    }
    return blocks;
}

function fromBlockMatrix(blocks: CholeskyBlockMatrix, blockSize: number,
    numBlocks: number): Matrix {
    const n = blockSize * numBlocks;
    const A = new Matrix(n, n);
    for (let r = 0; r < numBlocks; ++r) {
        for (let c = 0; c < numBlocks; ++c) {
            const block = blocks[c + r * numBlocks];
            for (let j = 0; j < blockSize; ++j) {
                for (let i = 0; i < blockSize; ++i) {
                    A.set(r * blockSize + j, c * blockSize + i,
                        block.get(j, i));
                }
            }
        }
    }
    return A;
}

function toBlockVector(v: Vector, blockSize: number, numBlocks: number):
    CholeskyBlockVector {
    const blocks: CholeskyBlockVector = [];
    for (let r = 0; r < numBlocks; ++r) {
        const block = new Vector(blockSize);
        for (let j = 0; j < blockSize; ++j) {
            block.set(j, v.get(r * blockSize + j));
        }
        blocks.push(block);
    }
    return blocks;
}

function fromBlockVector(blocks: CholeskyBlockVector, blockSize: number,
    numBlocks: number): Vector {
    const v = new Vector(blockSize * numBlocks);
    for (let r = 0; r < numBlocks; ++r) {
        for (let j = 0; j < blockSize; ++j) {
            v.set(r * blockSize + j, blocks[r].get(j));
        }
    }
    return v;
}

// The classic textbook example: A = L*L^T with an exactly representable L.
const A3 = () => Matrix.fromArray(3, 3, [
    4, 12, -16,
    12, 37, -43,
    -16, -43, 98
]);

describe('CholeskyDecomposition construction', () => {
    it('requires a positive size', () => {
        expect(() => new CholeskyDecomposition(0)).toThrow('Invalid size.');
        expect(() => new CholeskyDecomposition(-3)).toThrow('Invalid size.');
        expect(new CholeskyDecomposition(4).N).toBe(4);
    });
});

describe('CholeskyDecomposition.factor', () => {
    it('produces the known factor of the textbook 3x3 matrix', () => {
        const A = A3();
        const decomposer = new CholeskyDecomposition(3);
        expect(decomposer.factor(A)).toBe(true);

        // L = [[2,0,0],[6,1,0],[-8,5,3]].
        expect(A.get(0, 0)).toBeCloseTo(2, 12);
        expect(A.get(1, 0)).toBeCloseTo(6, 12);
        expect(A.get(1, 1)).toBeCloseTo(1, 12);
        expect(A.get(2, 0)).toBeCloseTo(-8, 12);
        expect(A.get(2, 1)).toBeCloseTo(5, 12);
        expect(A.get(2, 2)).toBeCloseTo(3, 12);
    });

    it('modifies only the lower-triangular portion', () => {
        const A = A3();
        new CholeskyDecomposition(3).factor(A);

        // The upper triangle keeps the original input values.
        expect(A.get(0, 1)).toBe(12);
        expect(A.get(0, 2)).toBe(-16);
        expect(A.get(1, 2)).toBe(-43);
    });

    it('reproduces the input as L*L^T', () => {
        const input = A3();
        const A = input.clone();
        new CholeskyDecomposition(3).factor(A);
        const L = lowerPart(A);
        expect(maxAbsDiff(multiplyABT(L, L), input)).toBeLessThan(1e-12);
    });

    it('handles the 1x1 case', () => {
        const A = Matrix.fromArray(1, 1, [9]);
        expect(new CholeskyDecomposition(1).factor(A)).toBe(true);
        expect(A.get(0, 0)).toBe(3);
    });

    it('returns false for a matrix with a nonpositive leading entry', () => {
        const A = Matrix.fromArray(2, 2, [0, 1, 1, 1]);
        expect(new CholeskyDecomposition(2).factor(A)).toBe(false);

        const B = Matrix.fromArray(2, 2, [-1, 0, 0, 1]);
        expect(new CholeskyDecomposition(2).factor(B)).toBe(false);
    });

    it('returns false for an indefinite symmetric matrix', () => {
        // A = [[1,2],[2,1]] has a negative second pivot (1 - 4 = -3).
        const A = Matrix.fromArray(2, 2, [1, 2, 2, 1]);
        expect(new CholeskyDecomposition(2).factor(A)).toBe(false);
    });

    it('returns false for a singular positive semidefinite matrix', () => {
        // A = [[1,1],[1,1]] has a zero second pivot.
        const A = Matrix.fromArray(2, 2, [1, 1, 1, 1]);
        expect(new CholeskyDecomposition(2).factor(A)).toBe(false);
    });

    it('throws when the matrix is not square', () => {
        const decomposer = new CholeskyDecomposition(3);
        expect(() => decomposer.factor(new Matrix(3, 4)))
            .toThrow('Matrix must be square.');
        expect(() => decomposer.factor(new Matrix(2, 2)))
            .toThrow('Matrix must be square.');
    });
});

describe('CholeskyDecomposition solvers', () => {
    it('solves A*x = b for the textbook 3x3 matrix', () => {
        const input = A3();
        const A = input.clone();
        const decomposer = new CholeskyDecomposition(3);
        expect(decomposer.factor(A)).toBe(true);

        const b = Vector.zero(3);
        b.set(0, 1); b.set(1, 2); b.set(2, 3);
        const x = b.clone();
        decomposer.solveLower(A, x);
        decomposer.solveUpper(A, x);

        expect(maxAbs(sub(mulMatrix(input, x), b))).toBeLessThan(1e-10);
    });

    it('solveLower solves the lower-triangular system exactly', () => {
        // L = [[2,0,0],[6,1,0],[-8,5,3]], y = (1,2,3), b = L*y.
        const L = Matrix.fromArray(3, 3, [2, 0, 0, 6, 1, 0, -8, 5, 3]);
        const y = Vector.zero(3);
        y.set(0, 1); y.set(1, 2); y.set(2, 3);
        const b = mulMatrix(L, y);
        const solved = b.clone();
        new CholeskyDecomposition(3).solveLower(L, solved);
        expect(maxAbs(sub(solved, y))).toBeLessThan(1e-12);
    });

    it('solveUpper solves the transposed system exactly', () => {
        const L = Matrix.fromArray(3, 3, [2, 0, 0, 6, 1, 0, -8, 5, 3]);
        const x = Vector.zero(3);
        x.set(0, -1); x.set(1, 4); x.set(2, 2);
        // b = L^T * x.
        const b = mulMatrix(x, L);
        const solved = b.clone();
        new CholeskyDecomposition(3).solveUpper(L, solved);
        expect(maxAbs(sub(solved, x))).toBeLessThan(1e-12);
    });

    it('throws on mismatched sizes', () => {
        const decomposer = new CholeskyDecomposition(3);
        expect(() => decomposer.solveLower(new Matrix(2, 2), new Vector(3)))
            .toThrow('Invalid size.');
        expect(() => decomposer.solveLower(new Matrix(3, 3), new Vector(2)))
            .toThrow('Invalid size.');
        expect(() => decomposer.solveUpper(new Matrix(3, 2), new Vector(3)))
            .toThrow('Invalid size.');
        expect(() => decomposer.solveUpper(new Matrix(3, 3), new Vector(4)))
            .toThrow('Invalid size.');
    });
});

describe('CholeskyDecomposition randomized cross-checks', () => {
    it('factors and solves random SPD matrices of sizes 1..8', () => {
        const rng = makeRng(20260901);
        let maxFactorError = 0;
        let maxResidual = 0;
        for (let n = 1; n <= 8; ++n) {
            for (let trial = 0; trial < 20; ++trial) {
                const input = randomSPD(n, rng);
                const A = input.clone();
                const decomposer = new CholeskyDecomposition(n);
                expect(decomposer.factor(A)).toBe(true);

                const L = lowerPart(A);
                maxFactorError = Math.max(maxFactorError,
                    maxAbsDiff(multiplyABT(L, L), input));

                const b = new Vector(n);
                for (let i = 0; i < n; ++i) {
                    b.set(i, 2 * rng() - 1);
                }
                const x = b.clone();
                decomposer.solveLower(A, x);
                decomposer.solveUpper(A, x);
                maxResidual = Math.max(maxResidual,
                    maxAbs(sub(mulMatrix(input, x), b)));
            }
        }
        expect(maxFactorError).toBeLessThan(1e-10);
        expect(maxResidual).toBeLessThan(1e-10);
    });
});

describe('BlockCholeskyDecomposition construction and element access', () => {
    it('requires positive sizes', () => {
        expect(() => new BlockCholeskyDecomposition(0, 2))
            .toThrow('Invalid input.');
        expect(() => new BlockCholeskyDecomposition(2, 0))
            .toThrow('Invalid input.');
        const decomposer = new BlockCholeskyDecomposition(2, 3);
        expect(decomposer.blockSize).toBe(2);
        expect(decomposer.numBlocks).toBe(3);
        expect(decomposer.numDimensions).toBe(6);
    });

    it('get/set address the full matrix as a 2D table of scalars', () => {
        const blockSize = 3, numBlocks = 2, n = blockSize * numBlocks;
        const decomposer = new BlockCholeskyDecomposition(blockSize, numBlocks);
        const blocks: CholeskyBlockMatrix = [];
        for (let i = 0; i < numBlocks * numBlocks; ++i) {
            blocks.push(new Matrix(blockSize, blockSize));
        }

        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                decomposer.set(blocks, r, c, r * n + c);
            }
        }
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                expect(decomposer.get(blocks, r, c)).toBe(r * n + c);
            }
        }
        // The blocked layout matches the flat matrix.
        const full = fromBlockMatrix(blocks, blockSize, numBlocks);
        expect(full.get(4, 1)).toBe(4 * n + 1);
    });
});

describe('BlockCholeskyDecomposition factor and solve', () => {
    // Exercise square and non-square (blockSize != numBlocks) partitions;
    // the latter is the case the upstream run-time specialization gets wrong.
    const cases: Array<[number, number]> = [
        [1, 4], [2, 2], [2, 3], [3, 2], [4, 1], [3, 3]
    ];

    it('reproduces the unblocked factorization', () => {
        const rng = makeRng(987654321);
        let maxError = 0;
        for (const [blockSize, numBlocks] of cases) {
            const n = blockSize * numBlocks;
            const input = randomSPD(n, rng);

            const unblocked = input.clone();
            expect(new CholeskyDecomposition(n).factor(unblocked)).toBe(true);
            const expectedL = lowerPart(unblocked);

            const blocks = toBlockMatrix(input, blockSize, numBlocks);
            const decomposer =
                new BlockCholeskyDecomposition(blockSize, numBlocks);
            expect(decomposer.factor(blocks)).toBe(true);
            const actualL =
                lowerPart(fromBlockMatrix(blocks, blockSize, numBlocks));

            maxError = Math.max(maxError, maxAbsDiff(actualL, expectedL));
        }
        expect(maxError).toBeLessThan(1e-10);
    });

    it('reproduces the input as L*L^T and solves A*x = b', () => {
        const rng = makeRng(13579);
        let maxFactorError = 0;
        let maxResidual = 0;
        for (const [blockSize, numBlocks] of cases) {
            const n = blockSize * numBlocks;
            const input = randomSPD(n, rng);
            const blocks = toBlockMatrix(input, blockSize, numBlocks);
            const decomposer =
                new BlockCholeskyDecomposition(blockSize, numBlocks);
            expect(decomposer.factor(blocks)).toBe(true);

            const L = lowerPart(fromBlockMatrix(blocks, blockSize, numBlocks));
            maxFactorError = Math.max(maxFactorError,
                maxAbsDiff(multiplyABT(L, L), input));

            const b = new Vector(n);
            for (let i = 0; i < n; ++i) {
                b.set(i, 2 * rng() - 1);
            }
            const x = toBlockVector(b, blockSize, numBlocks);
            decomposer.solveLower(blocks, x);
            decomposer.solveUpper(blocks, x);
            const solution = fromBlockVector(x, blockSize, numBlocks);
            maxResidual = Math.max(maxResidual,
                maxAbs(sub(mulMatrix(input, solution), b)));
        }
        expect(maxFactorError).toBeLessThan(1e-10);
        expect(maxResidual).toBeLessThan(1e-10);
    });

    it('matches the unblocked solution on the textbook 3x3 matrix with 3 1x1 blocks', () => {
        const input = A3();
        const b = Vector.zero(3);
        b.set(0, 1); b.set(1, 2); b.set(2, 3);

        const unblocked = input.clone();
        const scalar = new CholeskyDecomposition(3);
        scalar.factor(unblocked);
        const expected = b.clone();
        scalar.solveLower(unblocked, expected);
        scalar.solveUpper(unblocked, expected);

        const blocks = toBlockMatrix(input, 1, 3);
        const decomposer = new BlockCholeskyDecomposition(1, 3);
        expect(decomposer.factor(blocks)).toBe(true);
        const x = toBlockVector(b, 1, 3);
        decomposer.solveLower(blocks, x);
        decomposer.solveUpper(blocks, x);

        expect(maxAbs(sub(fromBlockVector(x, 1, 3), expected)))
            .toBeLessThan(1e-12);
    });

    it('a single block reduces to the unblocked decomposition', () => {
        const input = A3();
        const blocks = toBlockMatrix(input, 3, 1);
        const decomposer = new BlockCholeskyDecomposition(3, 1);
        expect(decomposer.factor(blocks)).toBe(true);
        const L = lowerPart(fromBlockMatrix(blocks, 3, 1));
        expect(maxAbsDiff(multiplyABT(L, L), input)).toBeLessThan(1e-12);
    });

    it('returns false when a diagonal block is not positive definite', () => {
        // The leading 2x2 block [[1,2],[2,1]] is indefinite.
        const A = Matrix.fromArray(4, 4, [
            1, 2, 0, 0,
            2, 1, 0, 0,
            0, 0, 5, 0,
            0, 0, 0, 5
        ]);
        const blocks = toBlockMatrix(A, 2, 2);
        expect(new BlockCholeskyDecomposition(2, 2).factor(blocks))
            .toBe(false);
    });

    it('returns false when a trailing Schur complement is not positive definite', () => {
        // A is symmetric but indefinite: the 2x2 trailing block of the Schur
        // complement fails.
        const A = Matrix.fromArray(4, 4, [
            4, 0, 4, 0,
            0, 4, 0, 4,
            4, 0, 4, 0,
            0, 4, 0, 4
        ]);
        const blocks = toBlockMatrix(A, 2, 2);
        expect(new BlockCholeskyDecomposition(2, 2).factor(blocks))
            .toBe(false);
    });
});

// ---------------------------------------------------------------------------
// Verification (V39): independent review against upstream
// CholeskyDecomposition.h.
// ---------------------------------------------------------------------------

describe('CholeskyDecomposition verification', () => {
    // A symmetric positive definite n-by-n matrix, well conditioned because
    // of the n*I shift.
    const spdMatrix = (n: number): fc.Arbitrary<Matrix> =>
        fc.array(scaled(-1, 1), { minLength: n * n, maxLength: n * n })
            .map(m => {
                const A = new Matrix(n, n);
                for (let r = 0; r < n; ++r) {
                    for (let c = 0; c < n; ++c) {
                        let sum = 0;
                        for (let k = 0; k < n; ++k) {
                            sum += m[r + n * k] * m[c + n * k];
                        }
                        A.set(r, c, sum + (r === c ? n : 0));
                    }
                }
                return A;
            });

    const sizedSpd = fc.integer({ min: 1, max: 7 })
        .chain(n => fc.tuple(fc.constant(n), spdMatrix(n),
            wellScaledVector(n, -5, 5)));

    it('produces a lower-triangular L with A = L*L^T', () => {
        check(sizedSpd, ([n, A, _b]) => {
            const original = A.clone();
            const decomposer = new CholeskyDecomposition(n);
            expect(decomposer.factor(A)).toBe(true);
            for (let r = 0; r < n; ++r) {
                for (let c = 0; c < n; ++c) {
                    if (c > r) {
                        // Only the lower-triangular portion is modified.
                        expect(A.get(r, c)).toBe(original.get(r, c));
                    } else {
                        let sum = 0;
                        for (let k = 0; k <= Math.min(r, c); ++k) {
                            sum += A.get(r, k) * A.get(c, k);
                        }
                        expectClose(sum, original.get(r, c), 1e-9, 1e-9);
                    }
                }
            }
            // The diagonal of L is positive, which is what makes the factor
            // unique.
            for (let r = 0; r < n; ++r) {
                expect(A.get(r, r)).toBeGreaterThan(0);
            }
        });
    });

    it('solves A*x = b the way the general linear solver does', () => {
        check(sizedSpd, ([n, A, b]) => {
            const original = A.clone();
            const decomposer = new CholeskyDecomposition(n);
            expect(decomposer.factor(A)).toBe(true);
            const x = b.clone();
            decomposer.solveLower(A, x);
            decomposer.solveUpper(A, x);
            const reference = LinearSystem.solve(n, original.values,
                b.values);
            expect(reference.invertible).toBe(true);
            for (let i = 0; i < n; ++i) {
                expectClose(x.get(i), reference.X[i], 1e-8, 1e-8);
            }
        });
    });

    it('succeeds exactly when every leading principal minor is positive',
        () => {
            // Sylvester's criterion. The pivot the factorization tests at
            // step c is minor(c+1)/minor(c), so the two conditions are the
            // same. Integer entries keep the minors far from zero, which is
            // what makes the sign comparison reliable in binary64.
            const symmetricLattice = fc.integer({ min: 1, max: 4 })
                .chain(n => fc.tuple(fc.constant(n),
                    fc.array(fc.integer({ min: -4, max: 4 }),
                        { minLength: n * n, maxLength: n * n })));
            check(symmetricLattice, ([n, values]) => {
                const A = new Matrix(n, n);
                for (let r = 0; r < n; ++r) {
                    for (let c = 0; c <= r; ++c) {
                        A.set(r, c, values[c + n * r]);
                        A.set(c, r, values[c + n * r]);
                    }
                }
                let allPositive = true;
                for (let k = 1; k <= n; ++k) {
                    const leading = new Matrix(k, k);
                    for (let r = 0; r < k; ++r) {
                        for (let c = 0; c < k; ++c) {
                            leading.set(r, c, A.get(r, c));
                        }
                    }
                    const minor = determinant(leading);
                    // A vanishing minor makes the criterion inconclusive for
                    // the floating-point pivot test; skip those draws.
                    fc.pre(Math.abs(minor) > 0.5);
                    allPositive = allPositive && minor > 0;
                }
                const decomposer = new CholeskyDecomposition(n);
                expect(decomposer.factor(A)).toBe(allPositive);
            });
        });

    it('block factorization agrees with the unblocked one', () => {
        // blockSize != numBlocks in every case below, which is what the
        // upstream run-time specialization gets wrong (#209): it indexes
        // inside a block with the block-level stride NumBlocks.
        const shapes: [number, number][] = [[2, 3], [3, 2], [1, 4], [4, 1],
        [2, 1], [1, 3]];
        for (const [blockSize, numBlocks] of shapes) {
            const n = blockSize * numBlocks;
            check(fc.tuple(spdMatrix(n), wellScaledVector(n, -5, 5)),
                ([A, b]) => {
                    const decomposer = new BlockCholeskyDecomposition(
                        blockSize, numBlocks);
                    const blocks: CholeskyBlockMatrix = [];
                    for (let i = 0; i < numBlocks * numBlocks; ++i) {
                        blocks.push(new Matrix(blockSize, blockSize));
                    }
                    for (let r = 0; r < n; ++r) {
                        for (let c = 0; c < n; ++c) {
                            decomposer.set(blocks, r, c, A.get(r, c));
                        }
                    }
                    expect(decomposer.factor(blocks)).toBe(true);

                    const Y: CholeskyBlockVector = [];
                    for (let r = 0; r < numBlocks; ++r) {
                        const v = new Vector(blockSize);
                        for (let j = 0; j < blockSize; ++j) {
                            v.set(j, b.get(j + blockSize * r));
                        }
                        Y.push(v);
                    }
                    decomposer.solveLower(blocks, Y);
                    decomposer.solveUpper(blocks, Y);

                    const reference = LinearSystem.solve(n, A.values,
                        b.values);
                    expect(reference.invertible).toBe(true);
                    for (let i = 0; i < n; ++i) {
                        const value = Y[Math.floor(i / blockSize)]
                            .get(i % blockSize);
                        expectClose(value, reference.X[i], 1e-8, 1e-8);
                    }
                }, 40);
        }
    });

    it('block get/set address the full matrix as a table of scalars', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 4 }),
            fc.integer({ min: 1, max: 4 })), ([blockSize, numBlocks]) => {
            const n = blockSize * numBlocks;
            const decomposer = new BlockCholeskyDecomposition(blockSize,
                numBlocks);
            const blocks: CholeskyBlockMatrix = [];
            for (let i = 0; i < numBlocks * numBlocks; ++i) {
                blocks.push(new Matrix(blockSize, blockSize));
            }
            for (let r = 0; r < n; ++r) {
                for (let c = 0; c < n; ++c) {
                    decomposer.set(blocks, r, c, c + n * r);
                }
            }
            for (let r = 0; r < n; ++r) {
                for (let c = 0; c < n; ++c) {
                    expect(decomposer.get(blocks, r, c)).toBe(c + n * r);
                }
            }
            // Every slot of every block was written exactly once, so no two
            // (row, col) pairs may collide.
            const seen = new Set<number>();
            for (const block of blocks) {
                for (const value of block.values) { seen.add(value); }
            }
            expect(seen.size).toBe(n * n);
        });
    });
});
