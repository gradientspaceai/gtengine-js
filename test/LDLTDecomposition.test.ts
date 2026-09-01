import { describe, it, expect } from 'vitest';
import {
    BlockLDLTDecomposition, LDLTDecomposition
} from '../src/LDLTDecomposition';
import type { LDLTBlockVector } from '../src/LDLTDecomposition';
import {
    Matrix, multiplyAB, multiplyABT, mulMatrix, transpose
} from '../src/Matrix';
import { Vector, sub } from '../src/Vector';

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

// L * D * L^T.
function reconstruct(L: Matrix, D: Matrix): Matrix {
    return multiplyAB(multiplyAB(L, D), transpose(L));
}

const A3 = () => Matrix.fromArray(3, 3, [
    4, 12, -16,
    12, 37, -43,
    -16, -43, 98
]);

describe('LDLTDecomposition construction', () => {
    it('requires a positive size', () => {
        expect(() => new LDLTDecomposition(0)).toThrow('Invalid size.');
        expect(() => new LDLTDecomposition(-2)).toThrow('Invalid size.');
        expect(new LDLTDecomposition(5).N).toBe(5);
    });
});

describe('LDLTDecomposition.factor', () => {
    it('produces the known factors of the textbook 3x3 matrix', () => {
        // A = L*L^T with L = [[2,0,0],[6,1,0],[-8,5,3]], so the LDLT factors
        // are D = diag(4,1,9) and unit lower L' = L * diag(1/2,1,1/3).
        const { success, L, D } = new LDLTDecomposition(3).factor(A3());
        expect(success).toBe(true);

        expect(D.get(0, 0)).toBeCloseTo(4, 12);
        expect(D.get(1, 1)).toBeCloseTo(1, 12);
        expect(D.get(2, 2)).toBeCloseTo(9, 12);
        // Off-diagonal entries of D are zero.
        expect(D.get(0, 1)).toBe(0);
        expect(D.get(2, 0)).toBe(0);

        expect(L.get(0, 0)).toBe(1);
        expect(L.get(1, 1)).toBe(1);
        expect(L.get(2, 2)).toBe(1);
        expect(L.get(1, 0)).toBeCloseTo(3, 12);
        expect(L.get(2, 0)).toBeCloseTo(-4, 12);
        expect(L.get(2, 1)).toBeCloseTo(5, 12);
        // L is lower triangular.
        expect(L.get(0, 1)).toBe(0);
        expect(L.get(0, 2)).toBe(0);
        expect(L.get(1, 2)).toBe(0);
    });

    it('reproduces the input as L*D*L^T', () => {
        const A = A3();
        const { L, D } = new LDLTDecomposition(3).factor(A);
        expect(maxAbsDiff(reconstruct(L, D), A)).toBeLessThan(1e-12);
    });

    it('does not modify the input matrix', () => {
        const A = A3();
        const copy = A.clone();
        new LDLTDecomposition(3).factor(A);
        expect(A.values).toEqual(copy.values);
    });

    it('handles the 1x1 case', () => {
        const { success, L, D } = new LDLTDecomposition(1)
            .factor(Matrix.fromArray(1, 1, [7]));
        expect(success).toBe(true);
        expect(L.get(0, 0)).toBe(1);
        expect(D.get(0, 0)).toBe(7);
    });

    it('factors an indefinite symmetric matrix, producing negative D entries', () => {
        // A = [[1,2],[2,1]] is indefinite: D = diag(1, -3), L(1,0) = 2.
        const A = Matrix.fromArray(2, 2, [1, 2, 2, 1]);
        const { success, L, D } = new LDLTDecomposition(2).factor(A);
        expect(success).toBe(true);
        expect(D.get(0, 0)).toBeCloseTo(1, 12);
        expect(D.get(1, 1)).toBeCloseTo(-3, 12);
        expect(L.get(1, 0)).toBeCloseTo(2, 12);
        expect(maxAbsDiff(reconstruct(L, D), A)).toBeLessThan(1e-12);
    });

    it('returns false when a pivot is exactly zero', () => {
        // A = [[0,1],[1,0]] has D(0,0) = 0.
        const zeroFirst = new LDLTDecomposition(2)
            .factor(Matrix.fromArray(2, 2, [0, 1, 1, 0]));
        expect(zeroFirst.success).toBe(false);

        // A = [[1,1],[1,1]] is singular: D(1,1) = 1 - 1 = 0.
        const zeroSecond = new LDLTDecomposition(2)
            .factor(Matrix.fromArray(2, 2, [1, 1, 1, 1]));
        expect(zeroSecond.success).toBe(false);
        expect(zeroSecond.D.get(0, 0)).toBe(1);
        expect(zeroSecond.D.get(1, 1)).toBe(0);
    });

    it('uses only the lower-triangular portion of the input', () => {
        const A = A3();
        const perturbed = A.clone();
        perturbed.set(0, 1, 1000);
        perturbed.set(0, 2, -2000);
        perturbed.set(1, 2, 3000);
        const expected = new LDLTDecomposition(3).factor(A);
        const actual = new LDLTDecomposition(3).factor(perturbed);
        expect(actual.L.values).toEqual(expected.L.values);
        expect(actual.D.values).toEqual(expected.D.values);
    });

    it('throws on a mismatched or non-square input', () => {
        const decomposer = new LDLTDecomposition(3);
        expect(() => decomposer.factor(new Matrix(3, 4)))
            .toThrow('Invalid size.');
        expect(() => decomposer.factor(new Matrix(2, 2)))
            .toThrow('Invalid size.');
    });
});

describe('LDLTDecomposition solvers', () => {
    it('solveFactored solves A*x = b', () => {
        const A = A3();
        const decomposer = new LDLTDecomposition(3);
        const { L, D } = decomposer.factor(A);
        const b = Vector.zero(3);
        b.set(0, 1); b.set(1, -2); b.set(2, 3);
        const x = decomposer.solveFactored(L, D, b);
        expect(maxAbs(sub(mulMatrix(A, x), b))).toBeLessThan(1e-10);
    });

    it('solveFactored does not modify B', () => {
        const A = A3();
        const decomposer = new LDLTDecomposition(3);
        const { L, D } = decomposer.factor(A);
        const b = Vector.zero(3);
        b.set(0, 1); b.set(1, -2); b.set(2, 3);
        const copy = b.clone();
        decomposer.solveFactored(L, D, b);
        expect(b.values).toEqual(copy.values);
    });

    it('solve factors during the call and reports success', () => {
        const A = A3();
        const b = Vector.zero(3);
        b.set(0, 5); b.set(1, 0); b.set(2, -1);
        const { success, X } = new LDLTDecomposition(3).solve(A, b);
        expect(success).toBe(true);
        expect(maxAbs(sub(mulMatrix(A, X), b))).toBeLessThan(1e-10);
    });

    it('solve reports failure for a singular matrix', () => {
        const A = Matrix.fromArray(2, 2, [1, 1, 1, 1]);
        const b = Vector.ones(2);
        const { success, X } = new LDLTDecomposition(2).solve(A, b);
        expect(success).toBe(false);
        expect(X.values).toEqual([0, 0]);
    });

    it('solves an indefinite symmetric system', () => {
        const A = Matrix.fromArray(2, 2, [1, 2, 2, 1]);
        const b = Vector.zero(2);
        b.set(0, 3); b.set(1, -1);
        const { success, X } = new LDLTDecomposition(2).solve(A, b);
        expect(success).toBe(true);
        expect(maxAbs(sub(mulMatrix(A, X), b))).toBeLessThan(1e-12);
    });

    it('throws on mismatched sizes', () => {
        const decomposer = new LDLTDecomposition(3);
        expect(() => decomposer.solveFactored(new Matrix(2, 2),
            new Matrix(3, 3), new Vector(3))).toThrow('Invalid size.');
        expect(() => decomposer.solveFactored(new Matrix(3, 3),
            new Matrix(3, 3), new Vector(2))).toThrow('Invalid size.');
        expect(() => decomposer.solve(new Matrix(3, 3), new Vector(4)))
            .toThrow('Invalid size.');
    });
});

describe('LDLTDecomposition randomized cross-checks', () => {
    it('factors and solves random SPD matrices of sizes 1..8', () => {
        const rng = makeRng(424242);
        let maxFactorError = 0;
        let maxResidual = 0;
        let maxPositive = -Infinity;
        for (let n = 1; n <= 8; ++n) {
            for (let trial = 0; trial < 20; ++trial) {
                const A = randomSPD(n, rng);
                const decomposer = new LDLTDecomposition(n);
                const { success, L, D } = decomposer.factor(A);
                expect(success).toBe(true);
                maxFactorError = Math.max(maxFactorError,
                    maxAbsDiff(reconstruct(L, D), A));

                // D is positive for an SPD input.
                for (let i = 0; i < n; ++i) {
                    maxPositive = Math.max(maxPositive, -D.get(i, i));
                }

                const b = new Vector(n);
                for (let i = 0; i < n; ++i) {
                    b.set(i, 2 * rng() - 1);
                }
                const x = decomposer.solveFactored(L, D, b);
                maxResidual = Math.max(maxResidual,
                    maxAbs(sub(mulMatrix(A, x), b)));
            }
        }
        expect(maxFactorError).toBeLessThan(1e-10);
        expect(maxResidual).toBeLessThan(1e-10);
        expect(maxPositive).toBeLessThan(0);
    });
});

describe('BlockLDLTDecomposition construction, conversion and access', () => {
    it('requires positive sizes', () => {
        expect(() => new BlockLDLTDecomposition(0, 2)).toThrow('Invalid size.');
        expect(() => new BlockLDLTDecomposition(2, 0)).toThrow('Invalid size.');
        const decomposer = new BlockLDLTDecomposition(3, 2);
        expect(decomposer.blockSize).toBe(3);
        expect(decomposer.numBlocks).toBe(2);
        expect(decomposer.numDimensions).toBe(6);
    });

    it('round-trips a matrix through the block layout', () => {
        const blockSize = 3, numBlocks = 2, n = 6;
        const decomposer = new BlockLDLTDecomposition(blockSize, numBlocks);
        const M = new Matrix(n, n);
        for (let i = 0; i < n * n; ++i) {
            M.values[i] = i + 1;
        }
        const blocks = decomposer.convertMatrixToBlock(M);
        expect(blocks.length).toBe(numBlocks * numBlocks);
        expect(blocks[0].numRows).toBe(blockSize);
        // Block (0,1) holds the upper-right 3x3 corner.
        expect(blocks[1].get(0, 0)).toBe(M.get(0, 3));
        expect(decomposer.convertBlockToMatrix(blocks).values)
            .toEqual(M.values);
    });

    it('round-trips a vector through the block layout', () => {
        const decomposer = new BlockLDLTDecomposition(3, 2);
        const V = new Vector(6);
        for (let i = 0; i < 6; ++i) {
            V.set(i, i * i);
        }
        const blocks = decomposer.convertVectorToBlock(V);
        expect(blocks.length).toBe(2);
        expect(blocks[1].get(0)).toBe(9);
        expect(decomposer.convertBlockToVector(blocks).values)
            .toEqual(V.values);
    });

    it('get/set address the full matrix as a 2D table of scalars', () => {
        const blockSize = 2, numBlocks = 3, n = 6;
        const decomposer = new BlockLDLTDecomposition(blockSize, numBlocks);
        const blocks = decomposer.convertMatrixToBlock(new Matrix(n, n));
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                decomposer.set(blocks, r, c, 10 * r + c);
            }
        }
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                expect(decomposer.get(blocks, r, c)).toBe(10 * r + c);
            }
        }
        expect(decomposer.convertBlockToMatrix(blocks).get(5, 2)).toBe(52);
    });

    it('validates sizes when verifySize is true', () => {
        const decomposer = new BlockLDLTDecomposition(2, 3);
        expect(() => decomposer.convertMatrixToBlock(new Matrix(4, 4)))
            .toThrow('Invalid size.');
        expect(() => decomposer.convertVectorToBlock(new Vector(4)))
            .toThrow('Invalid size.');
        expect(() => decomposer.convertBlockToMatrix([new Matrix(2, 2)]))
            .toThrow('Invalid size.');
        expect(() => decomposer.convertBlockToVector([new Vector(2)]))
            .toThrow('Invalid size.');
        expect(() => decomposer.get([new Matrix(2, 2)], 0, 0))
            .toThrow('Invalid size.');
        expect(() => decomposer.set([new Matrix(2, 2)], 0, 0, 1))
            .toThrow('Invalid size.');
        expect(() => decomposer.factor([new Matrix(2, 2)]))
            .toThrow('Invalid size.');
    });

    it('accepts blocks of the correct size when blockSize differs from numBlocks', () => {
        // Upstream's Convert(BlockVector, GVector) compares each block size
        // with NumBlocks instead of BlockSize; the port compares with
        // blockSize, so this 3x2 configuration is accepted.
        const decomposer = new BlockLDLTDecomposition(3, 2);
        const blocks: LDLTBlockVector = [new Vector(3), new Vector(3)];
        blocks[0].set(2, 7);
        blocks[1].set(0, 8);
        const V = decomposer.convertBlockToVector(blocks);
        expect(V.size).toBe(6);
        expect(V.get(2)).toBe(7);
        expect(V.get(3)).toBe(8);
    });
});

describe('BlockLDLTDecomposition factor and solve', () => {
    const cases: Array<[number, number]> = [
        [1, 4], [2, 2], [2, 3], [3, 2], [4, 1], [3, 3]
    ];

    it('reproduces the input as L*D*L^T for SPD block matrices', () => {
        const rng = makeRng(112233);
        let maxError = 0;
        for (const [blockSize, numBlocks] of cases) {
            const n = blockSize * numBlocks;
            const A = randomSPD(n, rng);
            const decomposer =
                new BlockLDLTDecomposition(blockSize, numBlocks);
            const ABlock = decomposer.convertMatrixToBlock(A);
            const { success, L, D } = decomposer.factor(ABlock);
            expect(success).toBe(true);

            const Lfull = decomposer.convertBlockToMatrix(L);
            const Dfull = decomposer.convertBlockToMatrix(D);
            maxError = Math.max(maxError,
                maxAbsDiff(reconstruct(Lfull, Dfull), A));
        }
        expect(maxError).toBeLessThan(1e-9);
    });

    it('produces identity diagonal blocks in L and zero off-diagonal blocks in L and D', () => {
        const rng = makeRng(556677);
        const decomposer = new BlockLDLTDecomposition(2, 3);
        const A = randomSPD(6, rng);
        const { success, L, D } = decomposer
            .factor(decomposer.convertMatrixToBlock(A));
        expect(success).toBe(true);
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 3; ++c) {
                const Lblock = L[c + 3 * r];
                const Dblock = D[c + 3 * r];
                if (r === c) {
                    expect(Lblock.values).toEqual([1, 0, 0, 1]);
                } else {
                    expect(Dblock.values).toEqual([0, 0, 0, 0]);
                    if (c > r) {
                        expect(Lblock.values).toEqual([0, 0, 0, 0]);
                    }
                }
            }
        }
    });

    it('solves A*x = b and matches the unblocked LDLT solution', () => {
        const rng = makeRng(778899);
        let maxResidual = 0;
        let maxDiff = 0;
        for (const [blockSize, numBlocks] of cases) {
            const n = blockSize * numBlocks;
            const A = randomSPD(n, rng);
            const b = new Vector(n);
            for (let i = 0; i < n; ++i) {
                b.set(i, 2 * rng() - 1);
            }

            const decomposer =
                new BlockLDLTDecomposition(blockSize, numBlocks);
            const { success, X } = decomposer.solve(
                decomposer.convertMatrixToBlock(A),
                decomposer.convertVectorToBlock(b));
            expect(success).toBe(true);
            const x = decomposer.convertBlockToVector(X);
            maxResidual = Math.max(maxResidual,
                maxAbs(sub(mulMatrix(A, x), b)));

            const expected = new LDLTDecomposition(n).solve(A, b);
            expect(expected.success).toBe(true);
            maxDiff = Math.max(maxDiff, maxAbs(sub(x, expected.X)));
        }
        expect(maxResidual).toBeLessThan(1e-9);
        expect(maxDiff).toBeLessThan(1e-9);
    });

    it('solveFactored matches solve', () => {
        const rng = makeRng(314159);
        const decomposer = new BlockLDLTDecomposition(3, 2);
        const A = randomSPD(6, rng);
        const b = new Vector(6);
        for (let i = 0; i < 6; ++i) {
            b.set(i, 2 * rng() - 1);
        }
        const ABlock = decomposer.convertMatrixToBlock(A);
        const BBlock = decomposer.convertVectorToBlock(b);

        const { success, L, D } = decomposer.factor(ABlock);
        expect(success).toBe(true);
        const X0 = decomposer.convertBlockToVector(
            decomposer.solveFactored(L, D, BBlock));
        const X1 = decomposer.convertBlockToVector(
            decomposer.solve(ABlock, BBlock).X);
        expect(X0.values).toEqual(X1.values);
    });

    it('a single 1x1 block reduces to the scalar solve', () => {
        const decomposer = new BlockLDLTDecomposition(1, 3);
        const A = A3();
        const b = Vector.zero(3);
        b.set(0, 1); b.set(1, 2); b.set(2, 3);
        const { success, X } = decomposer.solve(
            decomposer.convertMatrixToBlock(A),
            decomposer.convertVectorToBlock(b));
        expect(success).toBe(true);
        const expected = new LDLTDecomposition(3).solve(A, b);
        expect(maxAbs(sub(decomposer.convertBlockToVector(X), expected.X)))
            .toBeLessThan(1e-12);
    });

    it('returns false when a diagonal block of D is singular', () => {
        // The leading 2x2 block is singular, so its inverse does not exist.
        const A = Matrix.fromArray(4, 4, [
            1, 1, 0, 0,
            1, 1, 0, 0,
            0, 0, 5, 0,
            0, 0, 0, 5
        ]);
        const decomposer = new BlockLDLTDecomposition(2, 2);
        const ABlock = decomposer.convertMatrixToBlock(A);
        expect(decomposer.factor(ABlock).success).toBe(false);
        expect(decomposer.solve(ABlock,
            decomposer.convertVectorToBlock(Vector.ones(4))).success)
            .toBe(false);
    });

    it('factors an indefinite symmetric block matrix', () => {
        const A = Matrix.fromArray(4, 4, [
            1, 2, 0, 1,
            2, 1, 1, 0,
            0, 1, 3, 1,
            1, 0, 1, 3
        ]);
        const decomposer = new BlockLDLTDecomposition(2, 2);
        const { success, L, D } = decomposer
            .factor(decomposer.convertMatrixToBlock(A));
        expect(success).toBe(true);
        const Lfull = decomposer.convertBlockToMatrix(L);
        const Dfull = decomposer.convertBlockToMatrix(D);
        expect(maxAbsDiff(reconstruct(Lfull, Dfull), A)).toBeLessThan(1e-12);
    });
});
