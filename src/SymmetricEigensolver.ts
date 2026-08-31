// gtengine-js: TypeScript port of Geometric Tools Engine (GTE)
// SymmetricEigensolver.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The SymmetricEigensolver class is an implementation of Algorithm 8.2.3
// (Symmetric QR Algorithm) described in "Matrix Computations, 2nd edition"
// by G. H. Golub and C. F. Van Loan, The Johns Hopkins University Press,
// Baltimore MD, Fourth Printing 1993. Algorithm 8.2.1 (Householder
// Tridiagonalization) is used to reduce matrix A to tridiagonal T.
// Algorithm 8.2.2 (Implicit Symmetric QR Step with Wilkinson Shift) is used
// for the iterative reduction from tridiagonal to diagonal. If A is the
// original matrix, D is the diagonal matrix of eigenvalues, and Q is the
// orthogonal matrix of eigenvectors, then theoretically Q^T*A*Q = D.
// Numerically, we have errors E = Q^T*A*Q - D. Algorithm 8.2.3 mentions that
// one expects |E| is approximately u*|A|, where |M| denotes the Frobenius
// norm of M and where u is the unit roundoff for the floating-point
// arithmetic: 2^{-52} for double-precision arithmetic, which is
// Number.EPSILON = 2.2204460492503131e-16.
//
// The condition |a(i,i+1)| <= epsilon*(|a(i,i) + a(i+1,i+1)|) used to
// determine when the reduction decouples to smaller problems is implemented
// as: sum = |a(i,i)| + |a(i+1,i+1)|; sum + |a(i,i+1)| == sum. The idea is
// that the superdiagonal term is small relative to its diagonal neighbors,
// and so it is effectively zero.
//
// Once the tridiagonal matrix is available, the diagonal and superdiagonal
// entries are stored in linear arrays, ignoring the theoretically zero
// values not in the 3-band. The Householder vectors are stored in the
// lower-triangular portion of the matrix to save memory.
//
// Port notes. Matrices are row-major arrays of N*N numbers, as in the other
// ported linear-algebra solvers. Upstream fills caller-provided pointers for
// the eigenvalues and eigenvectors; per PORTING.md the port returns the
// arrays. Upstream uses the mutable scratch member mPVector inside
// GetEigenvector; the port uses local arrays there so that the returned
// array never aliases solver state (the arithmetic is unchanged).

// A Givens rotation is the identity with the following replacement entries:
// R(index,index) = cs, R(index,index+1) = sn, R(index+1,index) = -sn and
// R(index+1,index+1) = cs.
interface GivensRotation {
    index: number;
    cs: number;
    sn: number;
}

// Support for sorting eigenvalues.
interface SortItem {
    eigenvalue: number;
    index: number;
}

export class SymmetricEigensolver {
    // The value returned by solve(...) when the reduction of the tridiagonal
    // matrix to a diagonal matrix does not converge within the maximum
    // number of iterations.
    static readonly noConvergence: number = 0xFFFFFFFF;

    // The number N of rows and columns of the matrices to be processed.
    private mSize: number;

    // The maximum number of iterations for reducing the tridiagonal matrix
    // to a diagonal matrix.
    private mMaxIterations: number;

    // The internal copy of a matrix passed to the solver. See the comments
    // about tridiagonalize() for what is stored in the matrix.
    private mMatrix: number[];  // NxN elements

    // After the initial tridiagonalization by Householder reflections, the
    // full mMatrix is no longer needed. The diagonal and superdiagonal
    // entries are copied to linear arrays in order to be cache friendly.
    private mDiagonal: number[];  // N elements
    private mSuperdiagonal: number[];  // N-1 elements

    // The Givens rotations used to reduce the initial tridiagonal matrix to
    // a diagonal matrix.
    private mGivens: GivensRotation[];

    // When sorting is requested, the permutation associated with the sort is
    // stored in mPermutation. When sorting is not requested, mPermutation[0]
    // is set to -1. mVisited is used for finding cycles in the permutation.
    // mEigenvectorMatrixType is +1 if getEigenvectors returns a rotation
    // matrix, 0 if it returns a reflection matrix or -1 if an input to the
    // constructor is invalid.
    private mPermutation: number[];  // N elements
    private mVisited: number[];  // N elements
    private mEigenvectorMatrixType: number;

    // Temporary storage to compute Householder reflections and to support
    // sorting of eigenvectors.
    private mPVector: number[];  // N elements
    private mVVector: number[];  // N elements
    private mWVector: number[];  // N elements

    // The solver processes NxN symmetric matrices, where N > 1 ('size' is N)
    // and the matrix is stored in row-major order. The maximum number of
    // iterations ('maxIterations') must be specified for the reduction of a
    // tridiagonal matrix to a diagonal matrix. The goal is to compute NxN
    // orthogonal Q and NxN diagonal D for which Q^T*A*Q = D.
    constructor(size: number, maxIterations: number) {
        this.mSize = 0;
        this.mMaxIterations = 0;
        this.mEigenvectorMatrixType = -1;
        this.mMatrix = [];
        this.mDiagonal = [];
        this.mSuperdiagonal = [];
        this.mGivens = [];
        this.mPermutation = [];
        this.mVisited = [];
        this.mPVector = [];
        this.mVVector = [];
        this.mWVector = [];

        if (size > 1 && maxIterations > 0) {
            this.mSize = size;
            this.mMaxIterations = maxIterations;
            this.mMatrix = new Array<number>(size * size).fill(0);
            this.mDiagonal = new Array<number>(size).fill(0);
            this.mSuperdiagonal = new Array<number>(size - 1).fill(0);
            this.mPermutation = new Array<number>(size).fill(0);
            this.mVisited = new Array<number>(size).fill(0);
            this.mPVector = new Array<number>(size).fill(0);
            this.mVVector = new Array<number>(size).fill(0);
            this.mWVector = new Array<number>(size).fill(0);
        }
    }

    // A copy of the NxN symmetric input (row-major, N*N elements) is made
    // internally. The order of the eigenvalues is specified by sortType:
    // -1 (decreasing), 0 (no sorting), or +1 (increasing). When sorted, the
    // eigenvectors are ordered accordingly. The return value is the number
    // of iterations consumed when convergence occurred,
    // SymmetricEigensolver.noConvergence when convergence did not occur, or
    // 0 when N <= 1 was passed to the constructor.
    solve(input: readonly number[], sortType: number): number {
        this.mEigenvectorMatrixType = -1;

        if (this.mSize > 0) {
            for (let i = 0; i < this.mSize * this.mSize; ++i) {
                this.mMatrix[i] = input[i];
            }
            this.tridiagonalize();

            this.mGivens.length = 0;
            for (let j = 0; j < this.mMaxIterations; ++j) {
                let imin = -1, imax = -1;
                for (let i = this.mSize - 2; i >= 0; --i) {
                    // When a01 is much smaller than its diagonal neighbors,
                    // it is effectively zero.
                    const a00 = this.mDiagonal[i];
                    const a01 = this.mSuperdiagonal[i];
                    const a11 = this.mDiagonal[i + 1];
                    const sum = Math.abs(a00) + Math.abs(a11);
                    if (sum + Math.abs(a01) !== sum) {
                        if (imax === -1) {
                            imax = i;
                        }
                        imin = i;
                    }
                    else {
                        // The superdiagonal term is effectively zero
                        // compared to the neighboring diagonal terms.
                        if (imin >= 0) {
                            break;
                        }
                    }
                }

                if (imax === -1) {
                    // The algorithm has converged.
                    this.computePermutation(sortType);
                    return j;
                }

                // Process the lower-right-most unreduced tridiagonal block.
                this.doQRImplicitShift(imin, imax);
            }
            return SymmetricEigensolver.noConvergence;
        }
        else {
            return 0;
        }
    }

    // Get the eigenvalues of the matrix passed to solve(...). The returned
    // array has N elements (it is empty when N <= 1 was passed to the
    // constructor).
    getEigenvalues(): number[] {
        if (this.mSize > 0) {
            const eigenvalues = new Array<number>(this.mSize).fill(0);
            if (this.mPermutation[0] >= 0) {
                // Sorting was requested.
                for (let i = 0; i < this.mSize; ++i) {
                    eigenvalues[i] = this.mDiagonal[this.mPermutation[i]];
                }
            }
            else {
                // Sorting was not requested.
                for (let i = 0; i < this.mSize; ++i) {
                    eigenvalues[i] = this.mDiagonal[i];
                }
            }
            return eigenvalues;
        }
        return [];
    }

    // Accumulate the Householder reflections and Givens rotations to produce
    // the orthogonal matrix Q for which Q^T*A*Q = D. The returned array has
    // N*N elements, filled in as if the eigenvector matrix is stored in
    // row-major order. The i-th eigenvector is
    //   (eigenvectors[i+size*0], ... eigenvectors[i+size*(size - 1)])
    // which is the i-th column of 'eigenvectors' as an NxN matrix stored in
    // row-major order.
    getEigenvectors(): number[] {
        this.mEigenvectorMatrixType = -1;

        if (this.mSize === 0) {
            return [];
        }

        // Start with the identity matrix.
        const eigenvectors = new Array<number>(this.mSize * this.mSize).fill(0);
        for (let d = 0; d < this.mSize; ++d) {
            eigenvectors[d + this.mSize * d] = 1;
        }

        // Multiply the Householder reflections using backward accumulation.
        let r: number, c: number;
        for (let i = this.mSize - 3, rmin = this.mSize - 2; i >= 0; --i, --rmin) {
            // Copy the v vector and 2/Dot(v,v) from the matrix.
            const twoinvvdv = this.mMatrix[i + this.mSize * (i + 1)];
            for (r = 0; r < i + 1; ++r) {
                this.mVVector[r] = 0;
            }
            this.mVVector[r] = 1;
            for (++r; r < this.mSize; ++r) {
                this.mVVector[r] = this.mMatrix[i + this.mSize * r];
            }

            // Compute the w vector.
            for (r = 0; r < this.mSize; ++r) {
                this.mWVector[r] = 0;
                for (c = rmin; c < this.mSize; ++c) {
                    this.mWVector[r] += this.mVVector[c] * eigenvectors[r + this.mSize * c];
                }
                this.mWVector[r] *= twoinvvdv;
            }

            // Update the matrix, Q <- Q - v*w^T.
            for (r = rmin; r < this.mSize; ++r) {
                for (c = 0; c < this.mSize; ++c) {
                    eigenvectors[c + this.mSize * r] -= this.mVVector[r] * this.mWVector[c];
                }
            }
        }

        // Multiply the Givens rotations.
        for (const givens of this.mGivens) {
            for (r = 0; r < this.mSize; ++r) {
                const j = givens.index + this.mSize * r;
                const q0 = eigenvectors[j];
                const q1 = eigenvectors[j + 1];
                eigenvectors[j] = givens.cs * q0 - givens.sn * q1;
                eigenvectors[j + 1] = givens.sn * q0 + givens.cs * q1;
            }
        }

        // The number of Householder reflections is H = mSize - 2. If H is
        // even, the product of Householder reflections is a rotation;
        // otherwise, H is odd and the product is a reflection. The number of
        // Givens rotations does not influence the type of the product of
        // Householder reflections.
        this.mEigenvectorMatrixType = 1 - (this.mSize & 1);

        if (this.mPermutation[0] >= 0) {
            // Sorting was requested.
            this.mVisited.fill(0);
            for (let i = 0; i < this.mSize; ++i) {
                if (this.mVisited[i] === 0 && this.mPermutation[i] !== i) {
                    // The item starts a cycle with 2 or more elements.
                    const start = i;
                    let current = i;
                    let j: number, next: number;
                    for (j = 0; j < this.mSize; ++j) {
                        this.mPVector[j] = eigenvectors[i + this.mSize * j];
                    }
                    while ((next = this.mPermutation[current]) !== start) {
                        this.mEigenvectorMatrixType = 1 - this.mEigenvectorMatrixType;
                        this.mVisited[current] = 1;
                        for (j = 0; j < this.mSize; ++j) {
                            eigenvectors[current + this.mSize * j] =
                                eigenvectors[next + this.mSize * j];
                        }
                        current = next;
                    }
                    this.mVisited[current] = 1;
                    for (j = 0; j < this.mSize; ++j) {
                        eigenvectors[current + this.mSize * j] = this.mPVector[j];
                    }
                }
            }
        }

        return eigenvectors;
    }

    // The eigenvector matrix is a rotation (return +1) or a reflection
    // (return 0). If the input 'size' to the constructor is invalid, the
    // returned value is -1.
    getEigenvectorMatrixType(): number {
        return this.mEigenvectorMatrixType;
    }

    // Compute a single eigenvector, which amounts to computing column c of
    // matrix Q. The reflections and rotations are applied incrementally.
    // This is useful when you want only a small number of the eigenvectors.
    // The returned array has N elements; it is empty when c is invalid
    // (upstream leaves the caller's array unmodified in that case).
    getEigenvector(c: number): number[] {
        if (!(0 <= c && c < this.mSize)) {
            return [];
        }

        // y = H*x, then x and y are swapped for the next H.
        let x = new Array<number>(this.mSize).fill(0);
        let y = new Array<number>(this.mSize).fill(0);

        // Start with the Euclidean basis vector.
        if (this.mPermutation[0] >= 0) {
            // Sorting was requested.
            x[this.mPermutation[c]] = 1;
        }
        else {
            x[c] = 1;
        }

        // Apply the Givens rotations in reverse order.
        for (let k = this.mGivens.length - 1; k >= 0; --k) {
            const givens = this.mGivens[k];
            const xr = x[givens.index];
            const xrp1 = x[givens.index + 1];
            x[givens.index] = givens.cs * xr + givens.sn * xrp1;
            x[givens.index + 1] = -givens.sn * xr + givens.cs * xrp1;
        }

        // Apply the Householder reflections.
        for (let i = this.mSize - 3; i >= 0; --i) {
            // Get the Householder vector v.
            const twoinvvdv = this.mMatrix[i + this.mSize * (i + 1)];
            let r: number;
            for (r = 0; r <= i; ++r) {
                y[r] = x[r];
            }

            // Compute s = Dot(x,v) * 2/v^T*v. At this point r = i + 1 and
            // v[i+1] = 1.
            let s = x[r];
            for (let j = r + 1; j < this.mSize; ++j) {
                s += x[j] * this.mMatrix[i + this.mSize * j];
            }
            s *= twoinvvdv;

            y[r] = x[r] - s;  // v[i+1] = 1

            // Compute the remaining components of y.
            for (++r; r < this.mSize; ++r) {
                y[r] = x[r] - s * this.mMatrix[i + this.mSize * r];
            }

            const temp = x;
            x = y;
            y = temp;
        }

        // The final product is stored in x.
        return x;
    }

    getEigenvalue(c: number): number {
        if (this.mSize > 0) {
            if (this.mPermutation[0] >= 0) {
                // Sorting was requested.
                return this.mDiagonal[this.mPermutation[c]];
            }
            else {
                // Sorting was not requested.
                return this.mDiagonal[c];
            }
        }
        else {
            return Number.MAX_VALUE;
        }
    }

    // Tridiagonalize using Householder reflections. On input, mMatrix is a
    // copy of the input matrix. On output, the upper-triangular part of
    // mMatrix including the diagonal stores the tridiagonalization. The
    // lower-triangular part contains 2/Dot(v,v) that are used in computing
    // eigenvectors and the part below the subdiagonal stores the essential
    // parts of the Householder vectors v (the elements of v after the
    // leading 1-valued component).
    private tridiagonalize(): void {
        let r: number, c: number;
        for (let i = 0, ip1 = 1; i < this.mSize - 2; ++i, ++ip1) {
            // Compute the Householder vector. Read the initial vector from
            // the row of the matrix.
            let length = 0;
            for (r = 0; r < ip1; ++r) {
                this.mVVector[r] = 0;
            }
            for (r = ip1; r < this.mSize; ++r) {
                const vr = this.mMatrix[r + this.mSize * i];
                this.mVVector[r] = vr;
                length += vr * vr;
            }
            let vdv = 1;
            length = Math.sqrt(length);
            if (length > 0) {
                const v1 = this.mVVector[ip1];
                const sgn = (v1 >= 0 ? 1 : -1);
                const invDenom = 1 / (v1 + sgn * length);
                this.mVVector[ip1] = 1;
                for (r = ip1 + 1; r < this.mSize; ++r) {
                    const vr = this.mVVector[r] * invDenom;
                    this.mVVector[r] = vr;
                    vdv += vr * vr;
                }
            }

            // Compute the rank-1 offsets v*w^T and w*v^T.
            const invvdv = 1 / vdv;
            const twoinvvdv = invvdv * 2;

            // UPSTREAM BUG FIX (SymmetricEigensolver.h, Tridiagonalize).
            // When length is zero, the subcolumn below the subdiagonal is
            // already zero, mVVector is the zero vector and the reflection
            // applied to the matrix below is the identity. Upstream
            // nonetheless stores twoinvvdv = 2 as the reflection parameter
            // (see the assignment to mMatrix[i + mSize * ip1] at the end of
            // this iteration). GetEigenvectors/GetEigenvector rebuild the
            // reflection from that parameter with the implied v = e_{i+1},
            // producing H = I - 2*e_{i+1}*e_{i+1}^T rather than the identity
            // that was actually applied. The result is a spurious sign flip
            // of row i+1 of the eigenvector matrix, so Q no longer satisfies
            // A*Q = Q*D. This happens exactly when the tridiagonalization
            // decouples, for example for a block-diagonal input. The port
            // stores 0 instead, which makes the rebuilt reflection the
            // identity. The value used in the arithmetic below is unchanged
            // from upstream; since mVVector is identically zero in this
            // case, the two choices give the same mPVector and mWVector
            // (both zero), so the tridiagonal matrix is bit-for-bit the same
            // as upstream's and the eigenvalues are unaffected.
            const storedTwoInvVdV = (length > 0 ? twoinvvdv : 0);
            let pdvtvdv = 0;
            for (r = i; r < this.mSize; ++r) {
                this.mPVector[r] = 0;
                for (c = i; c < r; ++c) {
                    this.mPVector[r] += this.mMatrix[r + this.mSize * c] * this.mVVector[c];
                }
                for (/**/; c < this.mSize; ++c) {
                    this.mPVector[r] += this.mMatrix[c + this.mSize * r] * this.mVVector[c];
                }
                this.mPVector[r] *= twoinvvdv;
                pdvtvdv += this.mPVector[r] * this.mVVector[r];
            }

            pdvtvdv *= invvdv;
            for (r = i; r < this.mSize; ++r) {
                this.mWVector[r] = this.mPVector[r] - pdvtvdv * this.mVVector[r];
            }

            // Update the input matrix.
            for (r = i; r < this.mSize; ++r) {
                const vr = this.mVVector[r];
                const wr = this.mWVector[r];
                let offset = vr * wr * 2;
                this.mMatrix[r + this.mSize * r] -= offset;
                for (c = r + 1; c < this.mSize; ++c) {
                    offset = vr * this.mWVector[c] + wr * this.mVVector[c];
                    this.mMatrix[c + this.mSize * r] -= offset;
                }
            }

            // Copy the vector to column i of the matrix. The 0-valued
            // components at indices 0 through i are not stored. The 1-valued
            // component at index i+1 is also not stored; instead, the
            // quantity 2/Dot(v,v) is stored for use in eigenvector
            // construction. That construction must take into account the
            // implied components that are not stored.
            this.mMatrix[i + this.mSize * ip1] = storedTwoInvVdV;
            for (r = ip1 + 1; r < this.mSize; ++r) {
                this.mMatrix[i + this.mSize * r] = this.mVVector[r];
            }
        }

        // Copy the diagonal and subdiagonal entries for cache coherence in
        // the QR iterations.
        let k: number;
        const ksup = this.mSize - 1;
        let index = 0;
        const delta = this.mSize + 1;
        for (k = 0; k < ksup; ++k, index += delta) {
            this.mDiagonal[k] = this.mMatrix[index];
            this.mSuperdiagonal[k] = this.mMatrix[index + 1];
        }
        this.mDiagonal[k] = this.mMatrix[index];
    }

    // A helper for generating Givens rotation sine and cosine robustly.
    // Solves sn*x + cs*y = 0.
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

    // The QR step with implicit shift. Generally, the initial T is unreduced
    // tridiagonal (all subdiagonal entries are nonzero). If a QR step causes
    // a superdiagonal entry to become zero, the matrix decouples into a
    // block diagonal matrix with two tridiagonal blocks. These blocks can be
    // reduced independently of each other. The inputs imin and imax identify
    // the subblock of T to be processed. That block has upper-left element
    // T(imin,imin) and lower-right element T(imax,imax).
    private doQRImplicitShift(imin: number, imax: number): void {
        // The implicit shift. Compute the eigenvalue u of the lower-right
        // 2x2 block that is closer to a11.
        const a00 = this.mDiagonal[imax];
        const a01 = this.mSuperdiagonal[imax];
        let a11 = this.mDiagonal[imax + 1];
        const dif = (a00 - a11) * 0.5;
        const sgn = (dif >= 0 ? 1 : -1);
        const a01sqr = a01 * a01;
        const u = a11 - a01sqr / (dif + sgn * Math.sqrt(dif * dif + a01sqr));
        let x = this.mDiagonal[imin] - u;
        let y = this.mSuperdiagonal[imin];

        let a12: number, a22: number, a23: number;
        let tmp11: number, tmp12: number, tmp21: number, tmp22: number;
        let a02 = 0;
        for (let i0 = imin - 1, i1 = imin, i2 = imin + 1; i1 <= imax; ++i0, ++i1, ++i2) {
            // Compute the Givens rotation and save it for use in computing
            // the eigenvectors.
            const { cs, sn } = this.getSinCos(x, y);
            this.mGivens.push({ index: i1, cs: cs, sn: sn });

            // Update the tridiagonal matrix. This amounts to updating a 4x4
            // subblock,
            //   b00 b01 b02 b03
            //   b01 b11 b12 b13
            //   b02 b12 b22 b23
            //   b03 b13 b23 b33
            // The corners (b00, b03, b33) do not change values. The interior
            // block {{b11,b12},{b12,b22}} is updated on each pass. For the
            // first pass, the b0c values are out of range, so only the
            // values (b13, b23) change. For the last pass, the br3 values
            // are out of range, so only the values (b01, b02) change. For
            // passes between first and last, the values (b01, b02, b13, b23)
            // change.
            if (i1 > imin) {
                this.mSuperdiagonal[i0] = cs * this.mSuperdiagonal[i0] - sn * a02;
            }

            a11 = this.mDiagonal[i1];
            a12 = this.mSuperdiagonal[i1];
            a22 = this.mDiagonal[i2];
            tmp11 = cs * a11 - sn * a12;
            tmp12 = cs * a12 - sn * a22;
            tmp21 = sn * a11 + cs * a12;
            tmp22 = sn * a12 + cs * a22;
            this.mDiagonal[i1] = cs * tmp11 - sn * tmp12;
            this.mSuperdiagonal[i1] = sn * tmp11 + cs * tmp12;
            this.mDiagonal[i2] = sn * tmp21 + cs * tmp22;

            if (i1 < imax) {
                a23 = this.mSuperdiagonal[i2];
                a02 = -sn * a23;
                this.mSuperdiagonal[i2] = cs * a23;

                // Update the parameters for the next Givens rotation.
                x = this.mSuperdiagonal[i1];
                y = a02;
            }
        }
    }

    // Sort the eigenvalues and compute the corresponding permutation of the
    // indices of the array storing the eigenvalues. The permutation is used
    // for reordering the eigenvalues and eigenvectors in the calls to
    // getEigenvalues() and getEigenvectors().
    private computePermutation(sortType: number): void {
        // The number of Householder reflections is H = mSize - 2. If H is
        // even, the product of Householder reflections is a rotation;
        // otherwise, H is odd and the product is a reflection. The number of
        // Givens rotations does not influence the type of the product of
        // Householder reflections.
        this.mEigenvectorMatrixType = 1 - (this.mSize & 1);

        if (sortType === 0) {
            // Set a flag for getEigenvalues() and getEigenvectors() to know
            // that sorted output was not requested.
            this.mPermutation[0] = -1;
            return;
        }

        // Compute the permutation induced by sorting. Initially, we start
        // with the identity permutation I = (0,1,...,N-1).
        const items: SortItem[] = [];
        for (let i = 0; i < this.mSize; ++i) {
            items.push({ eigenvalue: this.mDiagonal[i], index: i });
        }

        if (sortType > 0) {
            items.sort((item0, item1) => {
                if (item0.eigenvalue < item1.eigenvalue) {
                    return -1;
                }
                if (item1.eigenvalue < item0.eigenvalue) {
                    return 1;
                }
                return 0;
            });
        }
        else {
            items.sort((item0, item1) => {
                if (item0.eigenvalue > item1.eigenvalue) {
                    return -1;
                }
                if (item1.eigenvalue > item0.eigenvalue) {
                    return 1;
                }
                return 0;
            });
        }

        for (let i = 0; i < this.mSize; ++i) {
            this.mPermutation[i] = items[i].index;
        }

        // getEigenvectors() has nontrivial code for computing the orthogonal
        // Q from the reflections and rotations. To avoid complicating the
        // code further when sorting is requested, Q is computed as in the
        // unsorted case. The columns of Q are then swapped to be consistent
        // with the sorting of the eigenvalues. To minimize copying due to
        // column swaps, the permutation P is used. The minimum number of
        // transpositions to obtain P from I is N minus the number of cycles
        // of P. Each cycle is reordered with a minimum number of
        // transpositions; that is, the eigenitems are cyclically swapped,
        // leading to a minimum amount of copying.
    }
}
