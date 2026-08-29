// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) UnsymmetricEigenvalues.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// An implementation of the QR algorithm described in "Matrix Computations,
// 2nd edition" by G. H. Golub and C. F. Van Loan, The Johns Hopkins
// University Press, Baltimore MD, Fourth Printing 1993. In particular, the
// implementation is based on Chapter 7 (The Unsymmetric Eigenvalue Problem),
// Section 7.5 (The Practical QR Algorithm).

export class UnsymmetricEigenvalues {
    // The number N of rows and columns of the matrices to be processed.
    private mSize: number;
    private mSizeM1: number;

    // The maximum number of iterations for reducing the tridiagonal matrix
    // to a diagonal matrix.
    private mMaxIterations: number;

    // The internal copy of a matrix passed to the solver, NxN elements in
    // row-major order.
    private mMatrix: number[];

    // Temporary storage to compute Householder reflections, N elements.
    private mX: number[];
    private mV: number[];
    private mScaledV: number[];
    private mW: number[];

    // Flags about the zeroness of the subdiagonal entries. This is used to
    // detect uncoupled submatrices and apply the QR algorithm to the
    // corresponding subproblems. The storage is padded on both ends with
    // zeros to avoid additional code logic when packing the eigenvalues for
    // access by the caller. Upstream stores a pointer mSubdiagonalFlag to
    // mFlagStorage[1]; the port uses index accessors that add the offset.
    private mFlagStorage: Int32Array;

    private mNumEigenvalues: number;
    private mEigenvalues: number[];

    // The solver processes NxN matrices (not necessarily symmetric), where
    // N >= 3 ('size' is N) and the matrix is stored in row-major order. The
    // maximum number of iterations ('maxIterations') must be specified for
    // reducing an upper Hessenberg matrix to an upper quasi-triangular
    // matrix (upper triangular matrix of blocks where the diagonal blocks
    // are 1x1 or 2x2). The goal is to compute the real-valued eigenvalues.
    constructor(size: number, maxIterations: number) {
        this.mSize = 0;
        this.mSizeM1 = 0;
        this.mMaxIterations = 0;
        this.mMatrix = [];
        this.mX = [];
        this.mV = [];
        this.mScaledV = [];
        this.mW = [];
        this.mFlagStorage = new Int32Array(0);
        this.mNumEigenvalues = 0;
        this.mEigenvalues = [];

        if (size >= 3 && maxIterations > 0) {
            this.mSize = size;
            this.mSizeM1 = size - 1;
            this.mMaxIterations = maxIterations;
            this.mMatrix = new Array<number>(size * size).fill(0);
            this.mX = new Array<number>(size).fill(0);
            this.mV = new Array<number>(size).fill(0);
            this.mScaledV = new Array<number>(size).fill(0);
            this.mW = new Array<number>(size).fill(0);
            this.mFlagStorage = new Int32Array(size + 1);
            this.mEigenvalues = new Array<number>(size).fill(0);
        }
    }

    // A copy of the NxN input (row-major, N*N elements) is made internally.
    // The order of the eigenvalues is specified by sortType: -1 (decreasing),
    // 0 (no sorting), or +1 (increasing). The return value is the number of
    // iterations consumed when convergence occurred, or 0 when an invalid
    // size or maxIterations was passed to the constructor. If the return
    // value equals maxIterations, convergence did not occur at some block.
    solve(input: readonly number[], sortType: number): number {
        if (this.mSize > 0) {
            for (let k = 0; k < this.mSize * this.mSize; ++k) {
                this.mMatrix[k] = input[k];
            }
            this.reduceToUpperHessenberg();

            const block: [number, number] = [0, 0];
            let found = this.getBlock(block);
            let numIterations: number;
            for (numIterations = 0; numIterations < this.mMaxIterations; ++numIterations) {
                if (found) {
                    // Solve the current subproblem.
                    this.francisQRStep(block[0], block[1] + 1);

                    // Find another subproblem (if any).
                    found = this.getBlock(block);
                }
                else {
                    break;
                }
            }

            // The matrix is fully uncoupled, upper Hessenberg with 1x1 or
            // 2x2 diagonal blocks. Golub and Van Loan call this "upper
            // quasi-triangular".
            this.mNumEigenvalues = 0;
            this.mEigenvalues.fill(0);
            for (let i = 0; i < this.mSizeM1; ++i) {
                if (this.getFlag(i) === 0) {
                    if (this.getFlag(i - 1) === 0) {
                        // We have a 1x1 block with a real eigenvalue.
                        this.mEigenvalues[this.mNumEigenvalues++] = this.a(i, i);
                    }
                }
                else {
                    if (this.getFlag(i - 1) === 0 && this.getFlag(i + 1) === 0) {
                        // We have a 2x2 block that might have real
                        // eigenvalues.
                        const a00 = this.a(i, i);
                        const a01 = this.a(i, i + 1);
                        const a10 = this.a(i + 1, i);
                        const a11 = this.a(i + 1, i + 1);
                        const tr = a00 + a11;
                        const det = a00 * a11 - a01 * a10;
                        const halfTr = tr * 0.5;
                        const discr = halfTr * halfTr - det;
                        if (discr >= 0) {
                            const rootDiscr = Math.sqrt(discr);
                            this.mEigenvalues[this.mNumEigenvalues++] = halfTr - rootDiscr;
                            this.mEigenvalues[this.mNumEigenvalues++] = halfTr + rootDiscr;
                        }
                    }
                    // else:
                    // The QR iteration failed to converge at this block. It
                    // must also be the case that
                    // numIterations == mMaxIterations. The caller will be
                    // aware of this when testing the returned numIterations.
                }
            }

            // Port note: upstream's packing loop above runs i in
            // [0, mSizeM1), so a trailing 1x1 block at index mSizeM1 (row
            // N-1 decoupled from row N-2) is handled here: when the last
            // subdiagonal entry is effectively zero, A(N-1,N-1) is a real
            // eigenvalue of a 1x1 block. Without this, a matrix that
            // converges to a fully upper-triangular form would report only
            // N-1 of its N real eigenvalues. The i-loop cannot simply be
            // extended to include i = mSizeM1 because getFlag(i + 1) would
            // read out of bounds; this replicates the loop body's 1x1-block
            // case for the final index.
            if (this.getFlag(this.mSizeM1) === 0 && this.getFlag(this.mSizeM1 - 1) === 0) {
                this.mEigenvalues[this.mNumEigenvalues++] = this.a(this.mSizeM1, this.mSizeM1);
            }

            if (sortType !== 0 && this.mNumEigenvalues > 1) {
                const sorted = this.mEigenvalues.slice(0, this.mNumEigenvalues);
                if (sortType > 0) {
                    sorted.sort((x, y) => x - y);
                }
                else {
                    sorted.sort((x, y) => y - x);
                }
                for (let i = 0; i < this.mNumEigenvalues; ++i) {
                    this.mEigenvalues[i] = sorted[i];
                }
            }

            return numIterations;
        }
        return 0;
    }

    // Get the real-valued eigenvalues of the matrix passed to solve(...).
    // Upstream writes through output parameters; the port returns the count
    // and a copy of the first numEigenvalues values.
    getEigenvalues(): { numEigenvalues: number; eigenvalues: number[] } {
        if (this.mSize > 0) {
            return {
                numEigenvalues: this.mNumEigenvalues,
                eigenvalues: this.mEigenvalues.slice(0, this.mNumEigenvalues)
            };
        }
        else {
            return { numEigenvalues: 0, eigenvalues: [] };
        }
    }

    // 2D accessors to elements of mMatrix[] (row-major).
    private a(r: number, c: number): number {
        return this.mMatrix[c + r * this.mSize];
    }

    private setA(r: number, c: number, value: number): void {
        this.mMatrix[c + r * this.mSize] = value;
    }

    // The port of mSubdiagonalFlag[i] = mFlagStorage[i + 1], valid for i in
    // [-1, mSize - 1] with zero padding at both ends.
    private getFlag(i: number): number {
        return this.mFlagStorage[i + 1];
    }

    private setFlag(i: number, value: number): void {
        this.mFlagStorage[i + 1] = value;
    }

    // Compute the Householder vector for (X[rmin],...,X[rmax]). The input
    // vector is stored in mX in the index range [rmin,rmax]. The output
    // vector V is stored in mV in the index range [rmin,rmax]. The scaled
    // vector is S = (-2/Dot(V,V))*V and is stored in mScaledV in the index
    // range [rmin,rmax].
    private house(rmin: number, rmax: number): void {
        let length = 0;
        for (let r = rmin; r <= rmax; ++r) {
            length += this.mX[r] * this.mX[r];
        }
        length = Math.sqrt(length);
        if (length !== 0) {
            const sign = (this.mX[rmin] >= 0 ? 1 : -1);
            const invDenom = 1 / (this.mX[rmin] + sign * length);
            for (let r = rmin + 1; r <= rmax; ++r) {
                this.mV[r] = this.mX[r] * invDenom;
            }
        }
        this.mV[rmin] = 1;

        let dot = 1;
        for (let r = rmin + 1; r <= rmax; ++r) {
            dot += this.mV[r] * this.mV[r];
        }
        const scale = -2 / dot;
        for (let r = rmin; r <= rmax; ++r) {
            this.mScaledV[r] = scale * this.mV[r];
        }
    }

    // Support for replacing matrix A by P^T*A*P, where P is a Householder
    // reflection computed using house(...).
    private rowHouse(rmin: number, rmax: number, cmin: number, cmax: number): void {
        for (let c = cmin; c <= cmax; ++c) {
            this.mW[c] = 0;
            for (let r = rmin; r <= rmax; ++r) {
                this.mW[c] += this.mScaledV[r] * this.a(r, c);
            }
        }

        for (let r = rmin; r <= rmax; ++r) {
            for (let c = cmin; c <= cmax; ++c) {
                this.setA(r, c, this.a(r, c) + this.mV[r] * this.mW[c]);
            }
        }
    }

    private colHouse(rmin: number, rmax: number, cmin: number, cmax: number): void {
        for (let r = rmin; r <= rmax; ++r) {
            this.mW[r] = 0;
            for (let c = cmin; c <= cmax; ++c) {
                this.mW[r] += this.mScaledV[c] * this.a(r, c);
            }
        }

        for (let r = rmin; r <= rmax; ++r) {
            for (let c = cmin; c <= cmax; ++c) {
                this.setA(r, c, this.a(r, c) + this.mW[r] * this.mV[c]);
            }
        }
    }

    private reduceToUpperHessenberg(): void {
        for (let c = 0, cp1 = 1; c <= this.mSize - 3; ++c, ++cp1) {
            for (let r = cp1; r <= this.mSizeM1; ++r) {
                this.mX[r] = this.a(r, c);
            }

            this.house(cp1, this.mSizeM1);
            this.rowHouse(cp1, this.mSizeM1, c, this.mSizeM1);
            this.colHouse(0, this.mSizeM1, cp1, this.mSizeM1);
        }
    }

    private francisQRStep(rmin: number, rmax: number): void {
        // Apply the double implicit shift step.
        const i0 = rmax - 1, i1 = rmax;
        const a00 = this.a(i0, i0);
        const a01 = this.a(i0, i1);
        const a10 = this.a(i1, i0);
        const a11 = this.a(i1, i1);
        const tr = a00 + a11;
        const det = a00 * a11 - a01 * a10;

        const j0 = rmin, j1 = j0 + 1, j2 = j1 + 1;
        const b00 = this.a(j0, j0);
        const b01 = this.a(j0, j1);
        const b10 = this.a(j1, j0);
        const b11 = this.a(j1, j1);
        const b21 = this.a(j2, j1);
        this.mX[rmin] = b00 * (b00 - tr) + b01 * b10 + det;
        this.mX[rmin + 1] = b10 * (b00 + b11 - tr);
        this.mX[rmin + 2] = b10 * b21;

        this.house(rmin, rmin + 2);
        this.rowHouse(rmin, rmin + 2, rmin, rmax);
        this.colHouse(rmin, Math.min(rmax, rmin + 3), rmin, rmin + 2);

        // Apply Householder reflections to restore the matrix to upper
        // Hessenberg form.
        for (let c = 0, cp1 = 1; c <= this.mSize - 3; ++c, ++cp1) {
            const kmax = Math.min(cp1 + 2, this.mSizeM1);
            for (let r = cp1; r <= kmax; ++r) {
                this.mX[r] = this.a(r, c);
            }

            this.house(cp1, kmax);
            this.rowHouse(cp1, kmax, c, this.mSizeM1);
            this.colHouse(0, this.mSizeM1, cp1, kmax);
        }
    }

    private getBlock(block: [number, number]): boolean {
        for (let i = 0; i < this.mSizeM1; ++i) {
            const a00 = this.a(i, i);
            const a11 = this.a(i + 1, i + 1);
            const a21 = this.a(i + 1, i);
            const sum0 = a00 + a11;
            const sum1 = sum0 + a21;
            this.setFlag(i, sum1 !== sum0 ? 1 : 0);
        }

        for (let i = 0; i < this.mSizeM1; ++i) {
            if (this.getFlag(i) === 1) {
                block[0] = i;
                block[1] = -1;
                while (i < this.mSizeM1 && this.getFlag(i) === 1) {
                    block[1] = i++;
                }
                if (block[1] !== block[0]) {
                    return true;
                }
            }
        }
        return false;
    }
}
