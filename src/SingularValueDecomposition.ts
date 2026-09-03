// gtengine-js: TypeScript port of Geometric Tools Engine (GTE)
// SingularValueDecomposition.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The SingularValueDecomposition class is an implementation of Algorithm
// 8.3.2 (The SVD Algorithm) described in "Matrix Computations, 2nd edition"
// by G. H. Golub and Charles F. Van Loan, The Johns Hopkins Press,
// Baltimore MD, Fourth Printing 1993. Algorithm 5.4.2 (Householder
// Bidiagonalization) is used to reduce A to bidiagonal B. Algorithm 8.3.1
// (Golub-Kahan SVD Step) is used for the iterative reduction from bidiagonal
// to diagonal. If A is the original matrix, S is the matrix whose diagonal
// entries are the singular values, and U and V are corresponding matrices,
// then theoretically U^T*A*V = S; equivalently, A = U*S*V^T. Numerically, we
// have errors E = U^T*A*V-S. Algorithm 8.3.2 mentions that one expects |E|
// is approximately unitRoundoff*|A|, where |A| denotes the Frobenius norm of
// A and where unitRoundoff is 2^{-52} for double-precision arithmetic, which
// is Number.EPSILON = 2.2204460492503131e-16.
//
// During the iterations that process B, the bidiagonalized A, a superdiagonal
// entry is determined to be effectively zero when compared to its neighboring
// diagonal and superdiagonal elements,
//   |b(i,i+1)| <= e * (|b(i,i) + b(i+1,i+1)|)
// The suggestion by Golub and van Loan is to choose e to be a small positive
// multiple of the unit roundoff, e = multiplier * unitRoundoff. A diagonal
// entry is determined to be effectively zero relative to a norm of B,
//   |b(i,i)| <= e * |B|
// The implementation uses the L-infinity norm for |B|, which is the largest
// absolute value of the diagonal and superdiagonal elements of B.
//
// The authors suggest that once you have the bidiagonal matrix, a practical
// implementation will store the diagonal and superdiagonal entries in linear
// arrays, ignoring the theoretically zero values not in the 2-band. The
// implementation uses separate storage for the Householder u-vectors, so the
// essential parts of these vectors are not stored in the local copy of A (as
// suggested by Golub and van Loan) in order to make the implementation more
// readable.

import { logAssert } from './Logger.js';

// A Givens rotation is the identity with the following replacement entries:
// R(index0,index0) = cs, R(index0,index1) = sn, R(index1,index0) = -sn and
// R(index1,index1) = cs.
interface GivensRotation {
    index0: number;
    index1: number;
    cs: number;
    sn: number;
}

// Support for sorting singular values.
interface SingularInfo {
    singular: number;
    inversePermute: number;
}

// The port of std::numeric_limits<Real>::epsilon() for 'double'.
const unitRoundoff = Number.EPSILON;

export class SingularValueDecomposition {
    // The value returned by solve(...) when the reduction of the bidiagonal
    // matrix to a diagonal matrix does not converge within the maximum
    // number of iterations. Upstream returns
    // std::numeric_limits<size_t>::max().
    static readonly invalid: number = Number.MAX_SAFE_INTEGER;

    // The number of rows and columns of the matrices to be processed.
    private mNumRows: number;
    private mNumCols: number;

    // The maximum number of iterations for reducing the bidiagonal matrix
    // to a diagonal matrix.
    private mMaxIterations: number;

    // The internal copy of a matrix passed to the solver, MxN elements
    // stored in row-major order.
    private mMatrix: number[];

    // The U-matrix (MxM), V-matrix (NxN) and S-matrix (MxN) for which
    // U^T*A*V = S. These are stored in row-major order.
    private mUMatrix: number[];
    private mVMatrix: number[];
    private mSMatrix: number[];

    // The diagonal (N elements) and superdiagonal (N-1 elements) of the
    // bidiagonalized matrix.
    private mDiagonal: number[];
    private mSuperdiagonal: number[];

    // The Householder reflections used to reduce the input matrix to a
    // bidiagonal matrix.
    private mLHouseholder: number[][];
    private mRHouseholder: number[][];

    // The Givens rotations used to reduce the initial bidiagonal matrix to a
    // diagonal matrix.
    private mLGivens: GivensRotation[];
    private mRGivens: GivensRotation[];

    // The solver processes MxN matrices, where M >= N > 1 ('numRows' is M and
    // 'numCols' is N) and the matrix is stored in row-major order. The
    // maximum number of iterations ('maxIterations') must be specified for
    // the reduction of a bidiagonal matrix to a diagonal matrix. The goal is
    // to compute MxM orthogonal U, NxN orthogonal V and MxN matrix S for
    // which U^T*A*V = S. The only nonzero entries of S are on the diagonal;
    // the diagonal entries are the singular values of the original matrix.
    constructor(numRows: number, numCols: number, maxIterations: number) {
        this.mNumRows = numRows;
        this.mNumCols = numCols;
        this.mMaxIterations = maxIterations;

        logAssert(
            this.mNumCols >= 2 && this.mNumRows >= this.mNumCols && this.mMaxIterations > 0,
            'Invalid input.');

        this.mMatrix = new Array<number>(this.mNumRows * this.mNumCols).fill(0);
        this.mUMatrix = new Array<number>(this.mNumRows * this.mNumRows).fill(0);
        this.mVMatrix = new Array<number>(this.mNumCols * this.mNumCols).fill(0);
        this.mSMatrix = new Array<number>(this.mNumRows * this.mNumCols).fill(0);

        this.mDiagonal = new Array<number>(this.mNumCols).fill(0);
        this.mSuperdiagonal = new Array<number>(this.mNumCols - 1).fill(0);

        this.mLHouseholder = [];
        for (let col = 0; col < this.mNumCols; ++col) {
            this.mLHouseholder.push(new Array<number>(this.mNumRows).fill(0));
        }

        this.mRHouseholder = [];
        for (let row = 0; row < this.mNumCols - 2; ++row) {
            this.mRHouseholder.push(new Array<number>(this.mNumCols).fill(0));
        }

        this.mLGivens = [];
        this.mRGivens = [];
    }

    // A copy of the MxN input (row-major, M*N elements) is made internally.
    // The multiplier is a small positive number used to compute e that is
    // described in the preamble comments of this file. The default is 8, but
    // you can adjust this based on the needs of your application. The return
    // value is the number of iterations consumed when convergence occurs or
    // SingularValueDecomposition.invalid when convergence does not occur.
    solve(input: readonly number[], multiplier: number = 8): number {
        logAssert(
            input.length >= this.mNumRows * this.mNumCols && multiplier > 0,
            'Invalid input to Solve.');

        // Copy the input to mMatrix. The latter matrix is modified
        // internally by the solver.
        for (let i = 0; i < this.mMatrix.length; ++i) {
            this.mMatrix[i] = input[i];
        }

        // Reduce mMatrix to bidiagonal form, storing the diagonal
        // mMatrix(d,d) and superdiagonal mMatrix(d,d+1) in mDiagonal and
        // mSuperdiagonal, respectively.
        this.bidiagonalize();

        // The threshold is used to determine whether a diagonal entry of the
        // bidiagonal matrix B is sufficiently small so that it is considered
        // to be zero. It is defined by
        //   threshold = multiplier * unitRoundoff * |B|
        // where |B| is a matrix norm and the multiplier is a small number
        // [as suggested before Algorithm 8.3.2 (The SVD Algorithm) in Golub
        // and Van Loan]. The L-infinity norm is used for B.
        const { epsilon, threshold } = this.computeCutoffs(multiplier);

        this.mRGivens.length = 0;
        this.mLGivens.length = 0;
        for (let iteration = 0; iteration < this.mMaxIterations; ++iteration) {
            // Set superdiagonal entries to zero if they are effectively zero
            // compared to the neighboring diagonal entries.
            let numZero = 0;
            for (let i = 0; i <= this.mNumCols - 2; ++i) {
                const absSuper = Math.abs(this.mSuperdiagonal[i]);
                const absDiag0 = Math.abs(this.mDiagonal[i]);
                const absDiag1 = Math.abs(this.mDiagonal[i + 1]);
                if (absSuper <= epsilon * (absDiag0 + absDiag1)) {
                    this.mSuperdiagonal[i] = 0;
                    ++numZero;
                }
            }

            if (numZero === this.mNumCols - 1) {
                // The superdiagonal terms are all effectively zero, so the
                // algorithm has converged. Compute U, V and S.
                this.computeOrthogonalMatrices();
                return iteration;
            }

            // Find the largest sequence {iMin,...,iMax} for which the
            // superdiagonal entries are all not effectively zero. On loop
            // termination, the sequence is valid because if all
            // superdiagonal terms were zero, the previous if-statement
            // "numZero === mNumCols - 1" would ensure an exit from the
            // function. Upstream uses unsigned indices that wrap around to
            // the sentinel std::numeric_limits<size_t>::max(); the port uses
            // the signed value -1 for the same purpose.
            let iMax: number;
            for (iMax = this.mNumCols - 2; iMax !== 0; --iMax) {
                if (this.mSuperdiagonal[iMax] !== 0) {
                    break;
                }
            }
            ++iMax;

            let iMin: number;
            for (iMin = iMax - 1; iMin !== -1; --iMin) {
                if (this.mSuperdiagonal[iMin] === 0) {
                    break;
                }
            }
            ++iMin;

            // The subblock corresponding to {iMin,...,iMax} has all
            // superdiagonal entries not effectively zero. Determine whether
            // this subblock has a diagonal entry that is effectively zero.
            // If it does, use Givens rotations to zero-out the row
            // containing that entry.
            if (this.diagonalEntriesNonzero(iMin, iMax, threshold)) {
                this.doGolubKahanStep(iMin, iMax);
            }
        }
        return SingularValueDecomposition.invalid;
    }

    // Get the U-matrix, which is MxM and stored in row-major order. Upstream
    // copies into a caller-provided pointer; the port returns a copy.
    getU(): number[] {
        return this.mUMatrix.slice();
    }

    // Get the V-matrix, which is NxN and stored in row-major order.
    getV(): number[] {
        return this.mVMatrix.slice();
    }

    // Get the S-matrix, which is MxN and stored in row-major order.
    getS(): number[] {
        return this.mSMatrix.slice();
    }

    getSingularValue(index: number): number {
        logAssert(
            index < this.mNumCols,
            'Invalid index for singular value.');

        return this.mSMatrix[index + this.mNumCols * index];
    }

    getUColumn(index: number): number[] {
        logAssert(
            index < this.mNumRows,
            'Invalid index or null pointer for U-column.');

        const uColumn = new Array<number>(this.mNumRows).fill(0);
        for (let row = 0; row < this.mNumRows; ++row) {
            uColumn[row] = this.mUMatrix[index + this.mNumRows * row];
        }
        return uColumn;
    }

    getVColumn(index: number): number[] {
        logAssert(
            index < this.mNumCols,
            'Invalid index or null pointer for V-column.');

        const vColumn = new Array<number>(this.mNumCols).fill(0);
        for (let row = 0; row < this.mNumCols; ++row) {
            vColumn[row] = this.mVMatrix[index + this.mNumCols * row];
        }
        return vColumn;
    }

    // Get the singular values, which is an N-element array.
    getSingularValues(): number[] {
        const singularValues = new Array<number>(this.mNumCols).fill(0);
        for (let index = 0; index < this.mNumCols; ++index) {
            singularValues[index] = this.mSMatrix[index + this.mNumCols * index];
        }
        return singularValues;
    }

    // Algorithm 5.1.1 (Householder Vector). The matrix A has size
    // numRows-by-numCols with numRows >= numCols and the vector v has size
    // numRows.
    private static computeHouseholderU(numRows: number, numCols: number,
        A: readonly number[], selectCol: number, v: number[]): void {
        // Extract the column vector v[] where v[row] = A(row, selectCol).
        // The elements v[row] = 0 for 0 <= row < selectCol to avoid
        // conceptual uninitialized memory; the caller should not reference
        // these elements.
        let row: number;
        for (row = 0; row < selectCol; ++row) {
            v[row] = 0;
        }

        let mu = 0;
        for (; row < numRows; ++row) {
            const element = A[selectCol + numCols * row];
            mu += element * element;
            v[row] = element;
        }
        mu = Math.sqrt(mu);

        if (mu !== 0) {
            const beta = v[selectCol] + (v[selectCol] >= 0 ? mu : -mu);
            for (row = selectCol + 1; row < numRows; ++row) {
                v[row] /= beta;
            }
        }
        v[selectCol] = 1;
    }

    // Algorithm 5.1.1 (Householder Vector). The matrix A has size
    // numRows-by-numCols with numRows >= numCols and the vector v has size
    // numCols.
    private static computeHouseholderV(numCols: number, A: readonly number[],
        selectRow: number, v: number[]): void {
        // Extract the row vector v[] where v[col] = A(selectRow, col). The
        // elements v[col] = 0 for 0 <= col <= selectRow to avoid conceptual
        // uninitialized memory; the caller should not reference these
        // elements. (Upstream also takes an unused numRows parameter.)
        const selectRowP1 = selectRow + 1;
        let col: number;
        for (col = 0; col < selectRowP1; ++col) {
            v[col] = 0;
        }

        let mu = 0;
        for (; col < numCols; ++col) {
            const element = A[col + numCols * selectRow];
            mu += element * element;
            v[col] = element;
        }
        mu = Math.sqrt(mu);

        if (mu !== 0) {
            const beta = v[selectRowP1] + (v[selectRowP1] >= 0 ? mu : -mu);
            for (col = selectRowP1 + 1; col < numCols; ++col) {
                v[col] /= beta;
            }
        }
        v[selectRowP1] = 1;
    }

    // Algorithm 5.1.2 (Householder Pre-Multiplication)
    private static doHouseholderPremultiply(numRows: number, numCols: number,
        v: readonly number[], selectCol: number, A: number[]): void {
        let vSqrLength = 0;
        for (let row = selectCol; row < numRows; ++row) {
            vSqrLength += v[row] * v[row];
        }
        const beta = -2 / vSqrLength;

        const w = new Array<number>(numCols).fill(0);
        for (let col = selectCol; col < numCols; ++col) {
            w[col] = 0;
            for (let row = selectCol; row < numRows; ++row) {
                w[col] += v[row] * A[col + numCols * row];
            }
            w[col] *= beta;
        }

        for (let row = selectCol; row < numRows; ++row) {
            for (let col = selectCol; col < numCols; ++col) {
                A[col + numCols * row] += v[row] * w[col];
            }
        }
    }

    // Algorithm 5.1.3 (Householder Post-Multiplication)
    private static doHouseholderPostmultiply(numRows: number, numCols: number,
        v: readonly number[], selectRow: number, A: number[]): void {
        let vSqrLength = 0;
        for (let col = selectRow; col < numCols; ++col) {
            vSqrLength += v[col] * v[col];
        }
        const beta = -2 / vSqrLength;

        const w = new Array<number>(numRows).fill(0);
        for (let row = selectRow; row < numRows; ++row) {
            w[row] = 0;
            for (let col = selectRow; col < numCols; ++col) {
                w[row] += A[col + numCols * row] * v[col];
            }
            w[row] *= beta;
        }

        for (let row = selectRow; row < numRows; ++row) {
            for (let col = selectRow; col < numCols; ++col) {
                A[col + numCols * row] += w[row] * v[col];
            }
        }
    }

    // Bidiagonalize using Householder reflections. On input, mMatrix is a
    // copy of the input matrix passed to solve(...). On output, mDiagonal
    // and mSuperdiagonal contain the bidiagonalized results.
    private bidiagonalize(): void {
        for (let i = 0; i < this.mNumCols; ++i) {
            // Compute the u-Householder vector.
            SingularValueDecomposition.computeHouseholderU(this.mNumRows, this.mNumCols,
                this.mMatrix, i, this.mLHouseholder[i]);

            // Update A = (I - 2*u*u^T/u^T*u) * A.
            SingularValueDecomposition.doHouseholderPremultiply(this.mNumRows, this.mNumCols,
                this.mLHouseholder[i], i, this.mMatrix);

            if (i < this.mRHouseholder.length) {
                // Compute the v-Householder vectors.
                SingularValueDecomposition.computeHouseholderV(this.mNumCols,
                    this.mMatrix, i, this.mRHouseholder[i]);

                // Update A = A * (I - 2*v*v^T/v^T*v).
                SingularValueDecomposition.doHouseholderPostmultiply(this.mNumRows, this.mNumCols,
                    this.mRHouseholder[i], i, this.mMatrix);
            }
        }

        // Copy the diagonal and superdiagonal for cache coherence in the
        // Golub-Kahan iterations.
        for (let d = 0; d < this.mNumCols; ++d) {
            this.mDiagonal[d] = this.mMatrix[d + this.mNumCols * d];
        }
        for (let s = 0; s < this.mNumCols - 1; ++s) {
            this.mSuperdiagonal[s] = this.mMatrix[(s + 1) + this.mNumCols * s];
        }
    }

    private computeCutoffs(multiplier: number): { epsilon: number; threshold: number } {
        let norm = 0;
        for (let i = 0; i < this.mNumCols; ++i) {
            const abs = Math.abs(this.mDiagonal[i]);
            if (abs > norm) {
                norm = abs;
            }
        }

        for (let i = 0; i < this.mNumCols - 1; ++i) {
            const abs = Math.abs(this.mSuperdiagonal[i]);
            if (abs > norm) {
                norm = abs;
            }
        }

        const epsilon = multiplier * unitRoundoff;
        return { epsilon, threshold: epsilon * norm };
    }

    // A helper for generating Givens rotation sine and cosine robustly when
    // solving sn * x + cs * y = 0.
    private getSinCos(x: number, y: number): { cs: number; sn: number } {
        let tau: number;
        let cs: number;
        let sn: number;
        if (y !== 0) {
            if (Math.abs(y) > Math.abs(x)) {
                tau = -x / y;
                sn = 1 / Math.sqrt(1 + tau * tau);
                cs = sn * tau;
            }
            else {
                tau = -y / x;
                cs = 1 / Math.sqrt(1 + tau * tau);
                sn = cs * tau;
            }
        }
        else {
            cs = 1;
            sn = 0;
        }
        return { cs, sn };
    }

    // Test for diagonal entries that are effectively zero through all but
    // the last. For each such entry, the B matrix decouples. Perform that
    // decoupling. If there are no zero-valued entries, then the Golub-Kahan
    // step must be performed.
    private diagonalEntriesNonzero(iMin: number, iMax: number, threshold: number): boolean {
        for (let i = iMin; i < iMax; ++i) {
            if (Math.abs(this.mDiagonal[i]) <= threshold) {
                // Use planar rotations to chase the superdiagonal entry out
                // of the matrix, which produces a row of zeros.
                let y = this.mSuperdiagonal[i];
                this.mSuperdiagonal[i] = 0;
                for (let j = i + 1; j <= iMax; ++j) {
                    const x = this.mDiagonal[j];
                    const { cs, sn } = this.getSinCos(x, y);
                    // NOTE: The Givens parameters are (cs,-sn). The negative
                    // sign is not a coding error.
                    this.mLGivens.push({ index0: i, index1: j, cs: cs, sn: -sn });
                    this.mDiagonal[j] = cs * x - sn * y;
                    if (j < iMax) {
                        const z = this.mSuperdiagonal[j];
                        this.mSuperdiagonal[j] = cs * z;
                        y = sn * z;
                    }
                }
                return false;
            }
        }
        return true;
    }

    // Algorithm 8.3.1 (Golub-Kahan SVD Step).
    private doGolubKahanStep(iMin: number, iMax: number): void {
        // The implicit shift. Let A = {{a00,a01},{a01,a11}} be the lower
        // right 2x2 block of B^T*B. Compute the eigenvalue u of A that is
        // closer to a11 than to a00.
        let f0: number, f1: number, d1: number, d2: number;
        if (iMax > 1) {
            f0 = this.mSuperdiagonal[iMax - 2];
            f1 = this.mSuperdiagonal[iMax - 1];
            d1 = this.mDiagonal[iMax - 1];
            d2 = this.mDiagonal[iMax];
        }
        else {
            f0 = 0;
            f1 = this.mSuperdiagonal[0];
            d1 = this.mDiagonal[0];
            d2 = this.mDiagonal[1];
        }

        // Compute the lower right 2x2 block of B^T*B.
        const a00 = d1 * d1 + f0 * f0;
        const a01 = d1 * f1;
        let a11 = d2 * d2 + f1 * f1;

        // The eigenvalues are ((a00+a11) +/- sqrt((a00-a11)^2+a01^2))/2,
        // which are equidistant from (a00+a11)/2. If a11 >= a00, the
        // required eigenvalue uses the (+) sqrt term. If a11 <= a00, the
        // required eigenvalue uses the (-) sqrt term.
        const sum = a00 + a11;
        const dif = a00 - a11;
        const root = Math.sqrt(dif * dif + a01 * a01);
        const lambda = 0.5 * (a11 >= a00 ? sum + root : sum - root);
        let x = this.mDiagonal[iMin] * this.mDiagonal[iMin] - lambda;
        let y = this.mDiagonal[iMin] * this.mSuperdiagonal[iMin];

        let a12: number, a21: number, a22: number, a23: number;
        let a02 = 0;
        for (let i0 = iMin - 1, i1 = iMin, i2 = iMin + 1; i1 <= iMax - 1; ++i0, ++i1, ++i2) {
            // Compute the Givens rotation G and save it for use in computing
            // V in U^T*A*V = S.
            let { cs, sn } = this.getSinCos(x, y);
            this.mRGivens.push({ index0: i1, index1: i2, cs: cs, sn: sn });

            // Update B0 = B*G.
            if (i1 > iMin) {
                this.mSuperdiagonal[i0] = cs * this.mSuperdiagonal[i0] - sn * a02;
            }

            a11 = this.mDiagonal[i1];
            a12 = this.mSuperdiagonal[i1];
            a22 = this.mDiagonal[i2];
            this.mDiagonal[i1] = cs * a11 - sn * a12;
            this.mSuperdiagonal[i1] = sn * a11 + cs * a12;
            this.mDiagonal[i2] = cs * a22;
            a21 = -sn * a22;

            // Update the parameters for the next Givens rotations.
            x = this.mDiagonal[i1];
            y = a21;

            // Compute the Givens rotation G and save it for use in computing
            // U in U^T*A*V = S.
            ({ cs, sn } = this.getSinCos(x, y));
            this.mLGivens.push({ index0: i1, index1: i2, cs: cs, sn: sn });

            // Update B1 = G^T*B0.
            a11 = this.mDiagonal[i1];
            a12 = this.mSuperdiagonal[i1];
            a22 = this.mDiagonal[i2];
            this.mDiagonal[i1] = cs * a11 - sn * a21;
            this.mSuperdiagonal[i1] = cs * a12 - sn * a22;
            this.mDiagonal[i2] = sn * a12 + cs * a22;

            if (i1 < iMax - 1) {
                a23 = this.mSuperdiagonal[i2];
                a02 = -sn * a23;
                this.mSuperdiagonal[i2] = cs * a23;

                // Update the parameters for the next Givens rotations.
                x = this.mSuperdiagonal[i1];
                y = a02;
            }
        }
    }

    private computeOrthogonalMatrices(): void {
        // Compute U and V given the current signed singular values.
        this.computeUOrthogonal();
        this.computeVOrthogonal();

        // Ensure the singular values are nonnegative. The sign changes are
        // absorbed by the U-matrix. The nonnegative values are stored in the
        // S-matrix.
        this.ensureNonnegativeSingularValues();

        // Sort the singular values in descending order. The sort
        // permutations are absorbed by the U-matrix and V-matrix.
        this.sortSingularValues();
    }

    private computeUOrthogonal(): void {
        // Start with the identity matrix for U.
        this.mUMatrix.fill(0);
        for (let d = 0; d < this.mNumRows; ++d) {
            this.mUMatrix[d + this.mNumRows * d] = 1;
        }

        // Multiply the Householder reflections using backward accumulation.
        // This allows doHouseholderPremultiply. A forward accumulation using
        // doHouseholderPostmultiply does not work because the semantics of
        // doHouseholderPostmultiply are slightly different from those of
        // doHouseholderPremultiply.
        for (let k = 0, col = this.mNumCols - 1; k <= this.mNumCols - 1; ++k, --col) {
            SingularValueDecomposition.doHouseholderPremultiply(this.mNumRows, this.mNumRows,
                this.mLHouseholder[col], col, this.mUMatrix);
        }

        // Multiply the Givens rotations using forward accumulation.
        for (const givens of this.mLGivens) {
            let j0 = givens.index0;
            let j1 = givens.index1;
            for (let row = 0; row < this.mNumRows; ++row, j0 += this.mNumRows, j1 += this.mNumRows) {
                const q0 = this.mUMatrix[j0];
                const q1 = this.mUMatrix[j1];
                this.mUMatrix[j0] = givens.cs * q0 - givens.sn * q1;
                this.mUMatrix[j1] = givens.sn * q0 + givens.cs * q1;
            }
        }
    }

    private computeVOrthogonal(): void {
        // Start with the identity matrix for V.
        this.mVMatrix.fill(0);
        for (let d = 0; d < this.mNumCols; ++d) {
            this.mVMatrix[d + this.mNumCols * d] = 1;
        }

        // Multiply the Householder reflections using backward accumulation.
        if (this.mNumCols >= 3) {
            for (let k = 0, col = this.mNumCols - 3; k <= this.mNumCols - 3; ++k, --col) {
                SingularValueDecomposition.doHouseholderPremultiply(this.mNumCols, this.mNumCols,
                    this.mRHouseholder[col], col, this.mVMatrix);
            }
        }

        // Multiply the Givens rotations using forward accumulation.
        for (const givens of this.mRGivens) {
            let j0 = givens.index0;
            let j1 = givens.index1;
            for (let col = 0; col < this.mNumCols; ++col, j0 += this.mNumCols, j1 += this.mNumCols) {
                const q0 = this.mVMatrix[j0];
                const q1 = this.mVMatrix[j1];
                this.mVMatrix[j0] = givens.cs * q0 - givens.sn * q1;
                this.mVMatrix[j1] = givens.sn * q0 + givens.cs * q1;
            }
        }
    }

    private ensureNonnegativeSingularValues(): void {
        this.mSMatrix.fill(0);
        for (let i = 0; i < this.mNumCols; ++i) {
            if (this.mDiagonal[i] >= 0) {
                this.mSMatrix[i + this.mNumCols * i] = this.mDiagonal[i];
            }
            else {
                this.mSMatrix[i + this.mNumCols * i] = -this.mDiagonal[i];
                for (let row = 0; row < this.mNumRows; ++row) {
                    const index = i + this.mNumRows * row;
                    this.mUMatrix[index] = -this.mUMatrix[index];
                }
            }
        }
    }

    private sortSingularValues(): void {
        // Sort the nonnegative singular values in descending order.
        // Upstream uses std::sort with std::greater<SingularInfo>, which
        // compares only the singular value; the port uses the equivalent
        // comparison function.
        const sorted: SingularInfo[] = [];
        for (let i = 0; i < this.mNumCols; ++i) {
            sorted.push({
                singular: this.mSMatrix[i + this.mNumCols * i],
                inversePermute: i
            });
        }
        sorted.sort((info0, info1) => {
            if (info0.singular > info1.singular) {
                return -1;
            }
            if (info1.singular > info0.singular) {
                return 1;
            }
            return 0;
        });
        for (let i = 0; i < this.mNumCols; ++i) {
            this.mSMatrix[i + this.mNumCols * i] = sorted[i].singular;
        }

        // Compute the inverse permutation of the sorting.
        const permute = new Array<number>(this.mNumCols).fill(0);
        for (let i = 0; i < this.mNumCols; ++i) {
            permute[sorted[i].inversePermute] = i;
        }

        // Permute the columns of the U-matrix to be consistent with the
        // sorted singular values.
        const sortedUMatrix = new Array<number>(this.mNumRows * this.mNumRows).fill(0);
        let col: number;
        for (col = 0; col < this.mNumCols; ++col) {
            for (let row = 0; row < this.mNumRows; ++row) {
                sortedUMatrix[permute[col] + this.mNumRows * row] =
                    this.mUMatrix[col + this.mNumRows * row];
            }
        }
        for (; col < this.mNumRows; ++col) {
            for (let row = 0; row < this.mNumRows; ++row) {
                sortedUMatrix[col + this.mNumRows * row] =
                    this.mUMatrix[col + this.mNumRows * row];
            }
        }
        this.mUMatrix = sortedUMatrix;

        // Permute the columns of the V-matrix to be consistent with the
        // sorted singular values.
        const sortedVMatrix = new Array<number>(this.mNumCols * this.mNumCols).fill(0);
        for (col = 0; col < this.mNumCols; ++col) {
            for (let row = 0; row < this.mNumCols; ++row) {
                sortedVMatrix[permute[col] + this.mNumCols * row] =
                    this.mVMatrix[col + this.mNumCols * row];
            }
        }
        this.mVMatrix = sortedVMatrix;
    }
}
