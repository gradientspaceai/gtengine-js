// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) QuarticRootsQR.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// An implementation of the QR algorithm described in "Matrix Computations,
// 2nd edition" by G. H. Golub and C. F. Van Loan, The Johns Hopkins
// University Press, Baltimore MD, Fourth Printing 1993. In particular, the
// implementation is based on Chapter 7 (The Unsymmetric Eigenvalue Problem),
// Section 7.5 (The Practical QR Algorithm). The algorithm is specialized for
// the companion matrix associated with a quartic polynomial.
//
// Port notes (mirroring the CubicRootsQR port): upstream exposes two
// operator() overloads; the port names them solve (from the polynomial
// coefficients) and solveMatrix (from an upper Hessenberg matrix). Upstream
// returns the iteration count and writes numRoots/roots to reference
// parameters; the port returns { iterations, numRoots, roots }.

import { CubicRootsQR, type CubicRootsQRMatrix } from './CubicRootsQR';

// A 4x4 matrix stored as rows: A[r][c] is row r, column c, matching the
// upstream std::array<std::array<Real, 4>, 4> layout.
export type QuarticRootsQRMatrix = number[][];

export interface QuarticRootsQRResult {
    // The number of Francis QR iterations applied. This equals the requested
    // maxIterations when the matrix did not uncouple (in which case numRoots
    // is 0). When the deflation path runs the cubic solver, the count
    // includes the cubic solver's iterations, as upstream.
    iterations: number;

    // The number of real roots found (0 to 4). Only the first numRoots
    // elements of 'roots' are valid; the rest are 0.
    numRoots: number;
    roots: number[];
}

export class QuarticRootsQR {
    // Solve p(x) = c0 + c1 * x + c2 * x^2 + c3 * x^3 + x^4 = 0.
    solve(maxIterations: number, c0: number, c1: number, c2: number, c3: number):
        QuarticRootsQRResult {
        // Create the companion matrix for the polynomial. The matrix is in
        // upper Hessenberg form.
        const A: QuarticRootsQRMatrix = [
            [0, 0, 0, -c0],
            [1, 0, 0, -c1],
            [0, 1, 0, -c2],
            [0, 0, 1, -c3]
        ];

        // Avoid the QR-cycle when c1 = c2 = 0 and avoid the slow convergence
        // when c1 and c2 are nearly zero.
        const V = [1, 0.36602540378443865, 0.36602540378443865];
        this.doIteration(V, A);

        return this.solveMatrix(maxIterations, A);
    }

    // Compute the real eigenvalues of the upper Hessenberg matrix A. The
    // matrix is modified by in-place operations, so if you need to remember
    // A, you must make your own copy before calling this function.
    solveMatrix(maxIterations: number, A: QuarticRootsQRMatrix): QuarticRootsQRResult {
        const roots = [0, 0, 0, 0];

        for (let numIterations = 0; numIterations < maxIterations; ++numIterations) {
            // Apply a Francis QR iteration.
            const tr = A[2][2] + A[3][3];
            const det = A[2][2] * A[3][3] - A[2][3] * A[3][2];
            const X = [
                A[0][0] * A[0][0] + A[0][1] * A[1][0] - tr * A[0][0] + det,
                A[1][0] * (A[0][0] + A[1][1] - tr),
                A[1][0] * A[2][1]
            ];
            const V = this.house(X);
            this.doIteration(V, A);

            // Test for uncoupling of A.
            const tr12 = A[1][1] + A[2][2];
            if (tr12 + A[2][1] === tr12) {
                let numRoots = this.getQuadraticRoots(0, 1, A, roots, 0);
                numRoots = this.getQuadraticRoots(2, 3, A, roots, numRoots);
                return { iterations: numIterations, numRoots, roots };
            }

            const tr01 = A[0][0] + A[1][1];
            if (tr01 + A[1][0] === tr01) {
                let numRoots = 1;
                roots[0] = A[0][0];

                // The cubic solver is not designed to process 3x3 submatrices
                // of an NxN matrix, so the copy of a submatrix of A to B is a
                // simple workaround for running the solver.
                const subMaxIterations = maxIterations - numIterations;
                const B: CubicRootsQRMatrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
                for (let r = 0, rp1 = 1; r < 3; ++r, ++rp1) {
                    for (let c = 0, cp1 = 1; c < 3; ++c, ++cp1) {
                        B[r][c] = A[rp1][cp1];
                    }
                }

                const sub = new CubicRootsQR().solveMatrix(subMaxIterations, B);
                for (let i = 0; i < sub.numRoots; ++i) {
                    roots[numRoots] = sub.roots[i];
                    ++numRoots;
                }
                return { iterations: numIterations + sub.iterations, numRoots, roots };
            }

            const tr23 = A[2][2] + A[3][3];
            if (tr23 + A[3][2] === tr23) {
                let numRoots = 1;
                roots[0] = A[3][3];

                // See the comment about the submatrix copy above.
                const subMaxIterations = maxIterations - numIterations;
                const B: CubicRootsQRMatrix = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 3; ++c) {
                        B[r][c] = A[r][c];
                    }
                }

                const sub = new CubicRootsQR().solveMatrix(subMaxIterations, B);
                for (let i = 0; i < sub.numRoots; ++i) {
                    roots[numRoots] = sub.roots[i];
                    ++numRoots;
                }
                return { iterations: numIterations + sub.iterations, numRoots, roots };
            }
        }
        return { iterations: maxIterations, numRoots: 0, roots };
    }

    private doIteration(V: number[], A: QuarticRootsQRMatrix): void {
        let multV = -2 / (V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
        let MV = [multV * V[0], multV * V[1], multV * V[2]];
        this.rowHouse(0, 2, 0, 3, V, MV, A);
        this.colHouse(0, 3, 0, 2, V, MV, A);

        const X = [A[1][0], A[2][0], A[3][0]];
        const locV = this.house(X);
        multV = -2 / (locV[0] * locV[0] + locV[1] * locV[1] + locV[2] * locV[2]);
        MV = [multV * locV[0], multV * locV[1], multV * locV[2]];
        this.rowHouse(1, 3, 0, 3, locV, MV, A);
        this.colHouse(0, 3, 1, 3, locV, MV, A);

        const Y = [A[2][1], A[3][1]];
        const W = this.house(Y);
        const multW = -2 / (W[0] * W[0] + W[1] * W[1]);
        const MW = [multW * W[0], multW * W[1]];
        this.rowHouse(2, 3, 0, 3, W, MW, A);
        this.colHouse(0, 3, 2, 3, W, MW, A);
    }

    // The upstream House<N> template; N is X.length (2 or 3).
    private house(X: number[]): number[] {
        const N = X.length;
        const V = new Array<number>(N).fill(0);
        let length = 0;
        for (let i = 0; i < N; ++i) {
            length += X[i] * X[i];
        }
        length = Math.sqrt(length);
        if (length !== 0) {
            const sign = (X[0] >= 0 ? 1 : -1);
            const denom = X[0] + sign * length;
            for (let i = 1; i < N; ++i) {
                V[i] = X[i] / denom;
            }
        }
        V[0] = 1;
        return V;
    }

    private rowHouse(rmin: number, rmax: number, cmin: number, cmax: number,
        V: number[], MV: number[], A: QuarticRootsQRMatrix): void {
        // Only the elements cmin through cmax are used.
        const W = [0, 0, 0, 0];

        for (let c = cmin; c <= cmax; ++c) {
            W[c] = 0;
            for (let r = rmin, k = 0; r <= rmax; ++r, ++k) {
                W[c] += V[k] * A[r][c];
            }
        }

        for (let r = rmin, k = 0; r <= rmax; ++r, ++k) {
            for (let c = cmin; c <= cmax; ++c) {
                A[r][c] += MV[k] * W[c];
            }
        }
    }

    private colHouse(rmin: number, rmax: number, cmin: number, cmax: number,
        V: number[], MV: number[], A: QuarticRootsQRMatrix): void {
        // Only the elements rmin through rmax are used.
        const W = [0, 0, 0, 0];

        for (let r = rmin; r <= rmax; ++r) {
            W[r] = 0;
            for (let c = cmin, k = 0; c <= cmax; ++c, ++k) {
                W[r] += V[k] * A[r][c];
            }
        }

        for (let r = rmin; r <= rmax; ++r) {
            for (let c = cmin, k = 0; c <= cmax; ++c, ++k) {
                A[r][c] += W[r] * MV[k];
            }
        }
    }

    // Appends the real roots (if any) of the quadratic associated with the
    // 2x2 block of A defined by indices i0 and i1, starting at index
    // numRoots of 'roots'. Returns the updated number of roots.
    private getQuadraticRoots(i0: number, i1: number, A: QuarticRootsQRMatrix,
        roots: number[], numRoots: number): number {
        // Solve x^2 - t * x + d = 0, where t is the trace and d is the
        // determinant of the 2x2 matrix defined by indices i0 and i1. The
        // discriminant is D = (t/2)^2 - d. When D >= 0, the roots are real
        // values t/2 - sqrt(D) and t/2 + sqrt(D). To avoid potential
        // numerical issues with subtractive cancellation, the roots are
        // computed as
        //   r0 = t/2 + sign(t/2)*sqrt(D), r1 = trace - r0.
        const trace = A[i0][i0] + A[i1][i1];
        const halfTrace = trace * 0.5;
        const determinant = A[i0][i0] * A[i1][i1] - A[i0][i1] * A[i1][i0];
        const discriminant = halfTrace * halfTrace - determinant;
        if (discriminant >= 0) {
            const sign = (trace >= 0 ? 1 : -1);
            const root = halfTrace + sign * Math.sqrt(discriminant);
            roots[numRoots++] = root;
            roots[numRoots++] = trace - root;
        }
        return numRoots;
    }
}
