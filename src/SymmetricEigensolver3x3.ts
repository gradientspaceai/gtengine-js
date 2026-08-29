// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SymmetricEigensolver3x3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The document
// https://www.geometrictools.com/Documentation/RobustEigenSymmetric3x3.pdf
// describes algorithms for solving the eigensystem associated with a 3x3
// symmetric real-valued matrix. The iterative algorithm is implemented by
// class SymmetricEigensolver3x3. The noniterative algorithm is implemented
// by class NISymmetricEigensolver3x3. The code does not use GTE objects.
//
// Port notes: upstream 'eval' cannot be used as an identifier in strict-mode
// JavaScript, so the eigenvalue/eigenvector names are 'evals' and 'evecs'.
// Output reference parameters become returned object literals (PORTING.md),
// except SortEigenstuff whose operator() mutates its eigenvalue/eigenvector
// arguments in place exactly as upstream; it is ported as method 'sort'.

// A 3-tuple; evals[i] is the eigenvalue of the (row) eigenvector evecs[i].
export type EigenTriple = [number, number, number];
export type EigenBasis3 = [EigenTriple, EigenTriple, EigenTriple];

export interface SymmetricEigensolver3x3Result {
    // The number of iterations used by the algorithm.
    iterations: number;
    evals: EigenTriple;
    evecs: EigenBasis3;
}

export interface NISymmetricEigensolver3x3Result {
    evals: EigenTriple;
    evecs: EigenBasis3;
}

export class SortEigenstuff {
    // Sorts in place. The isRotation parameter is passed by value upstream;
    // the local copy is modified during sorting.
    sort(sortType: number, isRotation: boolean, evals: EigenTriple, evecs: EigenBasis3): void {
        if (sortType !== 0) {
            // Sort the eigenvalues to evals[0] <= evals[1] <= evals[2].
            const index: [number, number, number] = [0, 0, 0];
            if (evals[0] < evals[1]) {
                if (evals[2] < evals[0]) {
                    // even permutation
                    index[0] = 2;
                    index[1] = 0;
                    index[2] = 1;
                }
                else if (evals[2] < evals[1]) {
                    // odd permutation
                    index[0] = 0;
                    index[1] = 2;
                    index[2] = 1;
                    isRotation = !isRotation;
                }
                else {
                    // even permutation
                    index[0] = 0;
                    index[1] = 1;
                    index[2] = 2;
                }
            }
            else {
                if (evals[2] < evals[1]) {
                    // odd permutation
                    index[0] = 2;
                    index[1] = 1;
                    index[2] = 0;
                    isRotation = !isRotation;
                }
                else if (evals[2] < evals[0]) {
                    // even permutation
                    index[0] = 1;
                    index[1] = 2;
                    index[2] = 0;
                }
                else {
                    // odd permutation
                    index[0] = 1;
                    index[1] = 0;
                    index[2] = 2;
                    isRotation = !isRotation;
                }
            }

            if (sortType === -1) {
                // The request is for evals[0] >= evals[1] >= evals[2]. This
                // requires an odd permutation, (i0,i1,i2) -> (i2,i1,i0).
                const temp = index[0];
                index[0] = index[2];
                index[2] = temp;
                isRotation = !isRotation;
            }

            const unorderedEVal: EigenTriple = [evals[0], evals[1], evals[2]];
            const unorderedEVec: EigenBasis3 = [evecs[0], evecs[1], evecs[2]];
            for (let j = 0; j < 3; ++j) {
                const i = index[j];
                evals[j] = unorderedEVal[i];
                evecs[j] = [unorderedEVec[i][0], unorderedEVec[i][1], unorderedEVec[i][2]];
            }
        }

        // Ensure the ordered eigenvectors form a right-handed basis.
        if (!isRotation) {
            for (let j = 0; j < 3; ++j) {
                evecs[2][j] = -evecs[2][j];
            }
        }
    }
}

// Return the exponent of the frexp decomposition |value| = m * 2^exponent
// with m in [1/2, 1) (the port of the std::frexp exponent output). Returns
// 0 when value is 0, matching std::frexp.
function frexpExponent(value: number): number {
    if (value === 0) {
        return 0;
    }
    const dv = new DataView(new ArrayBuffer(8));
    dv.setFloat64(0, Math.abs(value));
    const biased = Number(dv.getBigUint64(0) >> 52n);
    if (biased === 0) {
        // Subnormal. Scale by 2^64 to normalize, then remove the scaling
        // from the reported exponent.
        dv.setFloat64(0, Math.abs(value) * 18446744073709551616);
        const biasedScaled = Number(dv.getBigUint64(0) >> 52n);
        return biasedScaled - 1022 - 64;
    }
    return biased - 1022;
}

export class SymmetricEigensolver3x3 {
    // The input matrix must be symmetric, so only the unique elements must
    // be specified: a00, a01, a02, a11, a12, and a22.
    //
    // If 'aggressive' is 'true', the iterations occur until a superdiagonal
    // entry is exactly zero. If 'aggressive' is 'false', the iterations
    // occur until a superdiagonal entry is effectively zero compared to the
    // sum of magnitudes of its diagonal neighbors. Generally, the
    // nonaggressive convergence is acceptable.
    //
    // The order of the eigenvalues is specified by sortType: -1 (decreasing),
    // 0 (no sorting) or +1 (increasing). When sorted, the eigenvectors are
    // ordered accordingly, and {evecs[0], evecs[1], evecs[2]} is guaranteed
    // to be a right-handed orthonormal set. The 'iterations' field of the
    // result is the number of iterations used by the algorithm.
    solve(a00: number, a01: number, a02: number, a11: number, a12: number,
        a22: number, aggressive: boolean, sortType: number): SymmetricEigensolver3x3Result {
        // Compute the Householder reflection H0 and B = H0*A*H0, where
        // b02 = 0. H0 = {{c,s,0},{s,-c,0},{0,0,1}} with each inner triple a
        // row of H0.
        const zero = 0, one = 1, half = 0.5;
        let isRotation = false;
        let { c, s } = SymmetricEigensolver3x3.getCosSin(a12, -a02);
        let term0 = c * a00 + s * a01;
        let term1 = c * a01 + s * a11;
        let term2 = s * a00 - c * a01;
        let term3 = s * a01 - c * a11;
        let b00 = c * term0 + s * term1;
        let b01 = s * term0 - c * term1;
        // b02 = c * a02 + s * a12;  // 0
        let b11 = s * term2 - c * term3;
        let b12 = s * a02 - c * a12;
        let b22 = a22;

        // Maintain Q as the product of the reflections. Initially, Q = H0.
        // Updates by Givens reflections G are Q <- Q * G. The columns of the
        // final Q are the estimates for the eigenvectors.
        const Q: EigenBasis3 = [
            [c, s, zero],
            [s, -c, zero],
            [zero, zero, one]
        ];

        // The smallest subnormal number is 2^{-alpha}. The value alpha is
        // 1074 for 'double' (std::numeric_limits<double>::digits (53) minus
        // std::numeric_limits<double>::min_exponent (-1021)).
        const alpha = 1074;
        let i = 0, imax = 0, power = 0;
        let c2 = zero, s2 = zero;

        if (Math.abs(b12) <= Math.abs(b01)) {
            // It is known that |currentB12| < 2^{-i/2} * |initialB12|.
            // Compute imax so that 0 is the closest floating-point number
            // to 2^{-imax/2} * |initialB12|.
            power = frexpExponent(b12);
            imax = 2 * (power + alpha + 1);

            for (i = 0; i < imax; ++i) {
                // Compute the Givens reflection
                // G = {{c,0,-s},{s,0,c},{0,1,0}} where each inner triple is
                // a row of G.
                ({ c: c2, s: s2 } = SymmetricEigensolver3x3.getCosSin(half * (b00 - b11), b01));
                s = Math.sqrt(half * (one - c2));
                c = half * s2 / s;

                // Update Q <- Q * G.
                for (let r = 0; r < 3; ++r) {
                    term0 = c * Q[r][0] + s * Q[r][1];
                    term1 = Q[r][2];
                    term2 = c * Q[r][1] - s * Q[r][0];
                    Q[r][0] = term0;
                    Q[r][1] = term1;
                    Q[r][2] = term2;
                }
                isRotation = !isRotation;

                // Update B <- Q^T * B * Q, ensuring that b02 is zero and
                // |b12| has strictly decreased.
                term0 = c * b00 + s * b01;
                term1 = c * b01 + s * b11;
                term2 = s * b00 - c * b01;
                term3 = s * b01 - c * b11;
                // b02 = s * c * (b11 - b00) + (c * c - s * s) * b01; // 0
                b00 = c * term0 + s * term1;
                b01 = s * b12;
                b11 = b22;
                b12 = c * b12;
                b22 = s * term2 - c * term3;

                if (SymmetricEigensolver3x3.converged(aggressive, b00, b11, b01)) {
                    // Compute the Householder reflection
                    // H1 = {{c,s,0},{s,-c,0},{0,0,1}} where each inner
                    // triple is a row of H1.
                    ({ c: c2, s: s2 } = SymmetricEigensolver3x3.getCosSin(half * (b00 - b11), b01));
                    s = Math.sqrt(half * (one - c2));
                    c = half * s2 / s;

                    // Update Q <- Q * H1.
                    for (let r = 0; r < 3; ++r) {
                        term0 = c * Q[r][0] + s * Q[r][1];
                        term1 = s * Q[r][0] - c * Q[r][1];
                        Q[r][0] = term0;
                        Q[r][1] = term1;
                    }
                    isRotation = !isRotation;

                    // Compute the diagonal estimate D = Q^T * B * Q.
                    term0 = c * b00 + s * b01;
                    term1 = c * b01 + s * b11;
                    term2 = s * b00 - c * b01;
                    term3 = s * b01 - c * b11;
                    b00 = c * term0 + s * term1;
                    b11 = s * term2 - c * term3;
                    break;
                }
            }
        }
        else {
            // It is known that |currentB01| < 2^{-i/2} * |initialB01|.
            // Compute imax so that 0 is the closest floating-point number
            // to 2^{-imax/2} * |initialB01|.
            power = frexpExponent(b01);
            imax = 2 * (power + alpha + 1);

            for (i = 0; i < imax; ++i) {
                // Compute the Givens reflection
                // G = {{0,1,0},{c,0,-s},{s,0,c}} where each inner triple is
                // a row of G.
                ({ c: c2, s: s2 } = SymmetricEigensolver3x3.getCosSin(half * (b11 - b22), b12));
                s = Math.sqrt(half * (one - c2));
                c = half * s2 / s;

                // Update Q <- Q * G.
                for (let r = 0; r < 3; ++r) {
                    term0 = c * Q[r][1] + s * Q[r][2];
                    term1 = Q[r][0];
                    term2 = c * Q[r][2] - s * Q[r][1];
                    Q[r][0] = term0;
                    Q[r][1] = term1;
                    Q[r][2] = term2;
                }
                isRotation = !isRotation;

                // Update B <- Q^T * B * Q, ensuring that b02 is zero and
                // |b01| has strictly decreased.
                term0 = c * b11 + s * b12;
                term1 = c * b12 + s * b22;
                term2 = s * b11 - c * b12;
                term3 = s * b12 - c * b22;
                // b02 = s * c * (b22 - b11) + (c * c - s * s) * b12;  // 0
                b22 = s * term2 - c * term3;
                b12 = -s * b01;
                b11 = b00;
                b01 = c * b01;
                b00 = c * term0 + s * term1;

                if (SymmetricEigensolver3x3.converged(aggressive, b11, b22, b12)) {
                    // Compute the Householder reflection
                    // H1 = {{1,0,0},{0,c,s},{0,s,-c}} where each inner
                    // triple is a row of H1.
                    ({ c: c2, s: s2 } = SymmetricEigensolver3x3.getCosSin(half * (b11 - b22), b12));
                    s = Math.sqrt(half * (one - c2));
                    c = half * s2 / s;

                    // Update Q <- Q * H1.
                    for (let r = 0; r < 3; ++r) {
                        term0 = c * Q[r][1] + s * Q[r][2];
                        term1 = s * Q[r][1] - c * Q[r][2];
                        Q[r][1] = term0;
                        Q[r][2] = term1;
                    }
                    isRotation = !isRotation;

                    // Compute the diagonal estimate D = Q^T * B * Q.
                    term0 = c * b11 + s * b12;
                    term1 = c * b12 + s * b22;
                    term2 = s * b11 - c * b12;
                    term3 = s * b12 - c * b22;
                    b11 = c * term0 + s * term1;
                    b22 = s * term2 - c * term3;
                    break;
                }
            }
        }

        const evals: EigenTriple = [b00, b11, b22];
        const evecs: EigenBasis3 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        for (let row = 0; row < 3; ++row) {
            for (let col = 0; col < 3; ++col) {
                evecs[row][col] = Q[col][row];
            }
        }

        new SortEigenstuff().sort(sortType, isRotation, evals, evecs);
        return { iterations: i, evals, evecs };
    }

    // Normalize (u,v) to (c,s) with c <= 0 when (u,v) is not (0,0). If
    // (u,v) = (0,0), the function returns (c,s) = (-1,0). When used to
    // generate a Householder reflection, it does not matter whether (c,s)
    // or (-c,-s) is returned. When generating a Givens reflection,
    // c = cos(2*theta) and s = sin(2*theta). Having a negative cosine for
    // the double-angle term ensures that the single-angle terms
    // c = cos(theta) and s = sin(theta) satisfy |c| < 1/sqrt(2) < |s|.
    private static getCosSin(u: number, v: number): { c: number; s: number } {
        const length = Math.sqrt(u * u + v * v);
        if (length > 0) {
            let c = u / length;
            let s = v / length;
            if (c > 0) {
                c = -c;
                s = -s;
            }
            return { c, s };
        }
        else {
            return { c: -1, s: 0 };
        }
    }

    private static converged(aggressive: boolean, diagonal0: number,
        diagonal1: number, superdiagonal: number): boolean {
        if (aggressive) {
            // Test whether the superdiagonal term is zero.
            return superdiagonal === 0;
        }
        else {
            // Test whether the superdiagonal term is effectively zero
            // compared to its diagonal neighbors.
            const sum = Math.abs(diagonal0) + Math.abs(diagonal1);
            return sum + Math.abs(superdiagonal) === sum;
        }
    }
}

// Module-local 3-tuple arithmetic used by the noniterative solver (the port
// of its private static member functions).
function niMultiply(s: number, U: EigenTriple): EigenTriple {
    return [s * U[0], s * U[1], s * U[2]];
}

function niSubtract(U: EigenTriple, V: EigenTriple): EigenTriple {
    return [U[0] - V[0], U[1] - V[1], U[2] - V[2]];
}

function niDivide(U: EigenTriple, s: number): EigenTriple {
    const invS = 1 / s;
    return [U[0] * invS, U[1] * invS, U[2] * invS];
}

function niDot(U: EigenTriple, V: EigenTriple): number {
    return U[0] * V[0] + U[1] * V[1] + U[2] * V[2];
}

function niCross(U: EigenTriple, V: EigenTriple): EigenTriple {
    return [
        U[1] * V[2] - U[2] * V[1],
        U[2] * V[0] - U[0] * V[2],
        U[0] * V[1] - U[1] * V[0]
    ];
}

export class NISymmetricEigensolver3x3 {
    // The input matrix must be symmetric, so only the unique elements must
    // be specified: a00, a01, a02, a11, a12, and a22. The eigenvalues are
    // sorted in ascending order: evals[0] <= evals[1] <= evals[2] before the
    // sortType-based reordering is applied.
    solve(a00: number, a01: number, a02: number, a11: number, a12: number,
        a22: number, sortType: number): NISymmetricEigensolver3x3Result {
        // Precondition the matrix by factoring out the maximum absolute
        // value of the components. This guards against floating-point
        // overflow when computing the eigenvalues.
        const max0 = Math.max(Math.abs(a00), Math.abs(a01));
        const max1 = Math.max(Math.abs(a02), Math.abs(a11));
        const max2 = Math.max(Math.abs(a12), Math.abs(a22));
        const maxAbsElement = Math.max(Math.max(max0, max1), max2);
        if (maxAbsElement === 0) {
            // A is the zero matrix.
            return {
                evals: [0, 0, 0],
                evecs: [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
            };
        }

        const invMaxAbsElement = 1 / maxAbsElement;
        a00 *= invMaxAbsElement;
        a01 *= invMaxAbsElement;
        a02 *= invMaxAbsElement;
        a11 *= invMaxAbsElement;
        a12 *= invMaxAbsElement;
        a22 *= invMaxAbsElement;

        const evals: EigenTriple = [0, 0, 0];
        const evecs: EigenBasis3 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];

        const norm = a01 * a01 + a02 * a02 + a12 * a12;
        if (norm > 0) {
            // Compute the eigenvalues of A.

            // In the PDF mentioned previously, B = (A - q*I)/p, where
            // q = tr(A)/3 with tr(A) the trace of A (sum of the diagonal
            // entries of A) and where p = sqrt(tr((A - q*I)^2)/6).
            const q = (a00 + a11 + a22) / 3;

            // The matrix A - q*I is represented by the following, where b00,
            // b11 and b22 are computed after these comments,
            //   +-           -+
            //   | b00 a01 a02 |
            //   | a01 b11 a12 |
            //   | a02 a12 b22 |
            //   +-           -+
            const b00 = a00 - q;
            const b11 = a11 - q;
            const b22 = a22 - q;

            // This is the variable p mentioned in the PDF.
            const p = Math.sqrt((b00 * b00 + b11 * b11 + b22 * b22 + norm * 2) / 6);

            // We need det(B) = det((A - q*I)/p) = det(A - q*I)/p^3. The
            // value det(A - q*I) is computed using a cofactor expansion by
            // the first row of A - q*I. The cofactors are c00, c01 and c02
            // and the determinant is b00*c00 - a01*c01 + a02*c02. The det(B)
            // is then computed finally by the division with p^3.
            const c00 = b11 * b22 - a12 * a12;
            const c01 = a01 * b22 - a12 * a02;
            const c02 = a01 * a12 - b11 * a02;
            const det = (b00 * c00 - a01 * c01 + a02 * c02) / (p * p * p);

            // The halfDet value is cos(3*theta) mentioned in the PDF. The
            // acos(z) function requires |z| <= 1, but will fail silently and
            // return NaN if the input is larger than 1 in magnitude. To
            // avoid this problem due to rounding errors, the halfDet value
            // is clamped to [-1,1].
            let halfDet = det * 0.5;
            halfDet = Math.min(Math.max(halfDet, -1), 1);

            // The eigenvalues of B are ordered as beta0 <= beta1 <= beta2.
            // The number of digits in twoThirdsPi is chosen so that, whether
            // float or double, the floating-point number is the closest to
            // theoretical 2*pi/3.
            const angle = Math.acos(halfDet) / 3;
            const twoThirdsPi = 2.09439510239319549;
            const beta2 = Math.cos(angle) * 2;
            const beta0 = Math.cos(angle + twoThirdsPi) * 2;
            const beta1 = -(beta0 + beta2);

            // The eigenvalues of A are ordered as alpha0 <= alpha1 <= alpha2.
            evals[0] = q + p * beta0;
            evals[1] = q + p * beta1;
            evals[2] = q + p * beta2;

            // Compute the eigenvectors so that the set
            // {evecs[0], evecs[1], evecs[2]} is right handed and orthonormal.
            if (halfDet >= 0) {
                evecs[2] = NISymmetricEigensolver3x3.computeEigenvector0(a00, a01, a02, a11, a12, a22, evals[2]);
                evecs[1] = NISymmetricEigensolver3x3.computeEigenvector1(a00, a01, a02, a11, a12, a22, evecs[2], evals[1]);
                evecs[0] = niCross(evecs[1], evecs[2]);
            }
            else {
                evecs[0] = NISymmetricEigensolver3x3.computeEigenvector0(a00, a01, a02, a11, a12, a22, evals[0]);
                evecs[1] = NISymmetricEigensolver3x3.computeEigenvector1(a00, a01, a02, a11, a12, a22, evecs[0], evals[1]);
                evecs[2] = niCross(evecs[0], evecs[1]);
            }
        }
        else {
            // The matrix is diagonal.
            evals[0] = a00;
            evals[1] = a11;
            evals[2] = a22;
            evecs[0] = [1, 0, 0];
            evecs[1] = [0, 1, 0];
            evecs[2] = [0, 0, 1];
        }

        // The preconditioning scaled the matrix A, which scales the
        // eigenvalues. Revert the scaling.
        evals[0] *= maxAbsElement;
        evals[1] *= maxAbsElement;
        evals[2] *= maxAbsElement;

        new SortEigenstuff().sort(sortType, true, evals, evecs);
        return { evals, evecs };
    }

    // Robustly compute a right-handed orthonormal set { U, V, W }. The
    // vector W is guaranteed to be unit-length, in which case there is no
    // need to worry about a division by zero when computing invLength.
    private static computeOrthogonalComplement(W: EigenTriple): { U: EigenTriple; V: EigenTriple } {
        let U: EigenTriple;
        if (Math.abs(W[0]) > Math.abs(W[1])) {
            // The component of maximum absolute value is either W[0] or
            // W[2].
            const invLength = 1 / Math.sqrt(W[0] * W[0] + W[2] * W[2]);
            U = [-W[2] * invLength, 0, +W[0] * invLength];
        }
        else {
            // The component of maximum absolute value is either W[1] or
            // W[2].
            const invLength = 1 / Math.sqrt(W[1] * W[1] + W[2] * W[2]);
            U = [0, +W[2] * invLength, -W[1] * invLength];
        }
        const V = niCross(W, U);
        return { U, V };
    }

    private static computeEigenvector0(a00: number, a01: number, a02: number,
        a11: number, a12: number, a22: number, eval0: number): EigenTriple {
        // Compute a unit-length eigenvector for eigenvalue[i0]. The matrix
        // is rank 2, so two of the rows are linearly independent. For a
        // robust computation of the eigenvector, select the two rows whose
        // cross product has largest length of all pairs of rows.
        const row0: EigenTriple = [a00 - eval0, a01, a02];
        const row1: EigenTriple = [a01, a11 - eval0, a12];
        const row2: EigenTriple = [a02, a12, a22 - eval0];
        const r0xr1 = niCross(row0, row1);
        const r0xr2 = niCross(row0, row2);
        const r1xr2 = niCross(row1, row2);
        const d0 = niDot(r0xr1, r0xr1);
        const d1 = niDot(r0xr2, r0xr2);
        const d2 = niDot(r1xr2, r1xr2);

        let dmax = d0;
        let imax = 0;
        if (d1 > dmax) {
            dmax = d1;
            imax = 1;
        }
        if (d2 > dmax) {
            imax = 2;
        }

        if (imax === 0) {
            return niDivide(r0xr1, Math.sqrt(d0));
        }
        else if (imax === 1) {
            return niDivide(r0xr2, Math.sqrt(d1));
        }
        else {
            return niDivide(r1xr2, Math.sqrt(d2));
        }
    }

    private static computeEigenvector1(a00: number, a01: number, a02: number,
        a11: number, a12: number, a22: number, evec0: EigenTriple,
        eval1: number): EigenTriple {
        // Robustly compute a right-handed orthonormal set { U, V, evec0 }.
        const { U, V } = NISymmetricEigensolver3x3.computeOrthogonalComplement(evec0);

        // Let e be eval1 and let E be a corresponding eigenvector which is
        // a solution to the linear system (A - e*I)*E = 0. The matrix
        // (A - e*I) is 3x3, not invertible (so infinitely many solutions),
        // and has rank 2 when eval1 and eval2 are different. It has rank 1
        // when eval1 and eval2 are equal. Numerically, it is difficult to
        // compute robustly the rank of a matrix. Instead, the 3x3 linear
        // system is reduced to a 2x2 system as follows. Define the 3x2
        // matrix J = [U V] whose columns are the U and V computed
        // previously. Define the 2x1 vector X = J*E. The 2x2 system is
        // 0 = M * X = (J^T * (A - e*I) * J) * X where J^T is the transpose
        // of J and M = J^T * (A - e*I) * J is a 2x2 matrix. The system may
        // be written as
        //     +-                        -++-  -+       +-  -+
        //     | U^T*A*U - e  U^T*A*V     || x0 | = e * | x0 |
        //     | V^T*A*U      V^T*A*V - e || x1 |       | x1 |
        //     +-                        -++   -+       +-  -+
        // where X has row entries x0 and x1.

        const AU: EigenTriple = [
            a00 * U[0] + a01 * U[1] + a02 * U[2],
            a01 * U[0] + a11 * U[1] + a12 * U[2],
            a02 * U[0] + a12 * U[1] + a22 * U[2]
        ];

        const AV: EigenTriple = [
            a00 * V[0] + a01 * V[1] + a02 * V[2],
            a01 * V[0] + a11 * V[1] + a12 * V[2],
            a02 * V[0] + a12 * V[1] + a22 * V[2]
        ];

        let m00 = U[0] * AU[0] + U[1] * AU[1] + U[2] * AU[2] - eval1;
        let m01 = U[0] * AV[0] + U[1] * AV[1] + U[2] * AV[2];
        let m11 = V[0] * AV[0] + V[1] * AV[1] + V[2] * AV[2] - eval1;

        // For robustness, choose the largest-length row of M to compute the
        // eigenvector. The 2-tuple of coefficients of U and V in the
        // assignments to eigenvector[1] lies on a circle, and U and V are
        // unit length and perpendicular, so eigenvector[1] is unit length
        // (within numerical tolerance).
        const absM00 = Math.abs(m00);
        const absM01 = Math.abs(m01);
        const absM11 = Math.abs(m11);
        let maxAbsComp: number;
        if (absM00 >= absM11) {
            maxAbsComp = Math.max(absM00, absM01);
            if (maxAbsComp > 0) {
                if (absM00 >= absM01) {
                    m01 /= m00;
                    m00 = 1 / Math.sqrt(1 + m01 * m01);
                    m01 *= m00;
                }
                else {
                    m00 /= m01;
                    m01 = 1 / Math.sqrt(1 + m00 * m00);
                    m00 *= m01;
                }
                return niSubtract(niMultiply(m01, U), niMultiply(m00, V));
            }
            else {
                return U;
            }
        }
        else {
            maxAbsComp = Math.max(absM11, absM01);
            if (maxAbsComp > 0) {
                if (absM11 >= absM01) {
                    m01 /= m11;
                    m11 = 1 / Math.sqrt(1 + m01 * m01);
                    m01 *= m11;
                }
                else {
                    m11 /= m01;
                    m01 = 1 / Math.sqrt(1 + m11 * m11);
                    m11 *= m01;
                }
                return niSubtract(niMultiply(m11, U), niMultiply(m01, V));
            }
            else {
                return U;
            }
        }
    }
}
