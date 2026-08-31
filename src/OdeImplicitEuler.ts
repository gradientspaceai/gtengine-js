// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) OdeImplicitEuler.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The function F(t,x) has input t, a scalar, and input x, an N-vector. The
// first derivative matrix with respect to x is DF(t,x), an N-by-N matrix.
// Entry DF(r,c) is the derivative of F[r] with respect to x[c].
//
// Port notes. Upstream is templated on TVector (Vector<N,Real> or
// GVector<Real>) and on the corresponding TMatrix (Matrix<N,N,Real> or
// GMatrix<Real>). The port has a single Vector class whose dimension is a
// run-time value, so the solver is a concrete class derived from
// OdeSolver<Vector>. Matrix.h/GMatrix.h are not yet ported, so DF returns
// the N-by-N derivative matrix as a row-major array of N*N numbers, the same
// storage convention used by the other ported linear-algebra solvers. The
// upstream call Inverse(dgMatrix), which uses the full-pivoting Gaussian
// elimination of GaussianElimination.h, is replaced by the module-private
// inverse below (Gauss-Jordan with partial pivoting); as upstream, a
// noninvertible matrix produces the zero matrix rather than an exception.

import { OdeSolver, type OdeFunction } from './OdeSolver';
import { Vector, add, mul } from './Vector';

// The port of 'std::function<TMatrix(Real, TVector const&)>'. The returned
// array stores the N-by-N matrix in row-major order: DF(r,c) is element
// [c + N * r].
export type OdeDerivativeFunction = (t: number, x: Vector) => number[];

export class OdeImplicitEuler extends OdeSolver<Vector> {
    private mDerivativeFunction: OdeDerivativeFunction;

    constructor(tDelta: number, F: OdeFunction<Vector>, DF: OdeDerivativeFunction) {
        super(tDelta, F);
        this.mDerivativeFunction = DF;
    }

    // Estimate x(t + tDelta) from x(t) using dx/dt = F(t,x). The estimate is
    // a new Vector, so xIn is never modified.
    update(tIn: number, xIn: Vector): { tOut: number; xOut: Vector } {
        let fVector = this.mFunction(tIn, xIn);
        const dfMatrix = this.mDerivativeFunction(tIn, xIn);
        const n = xIn.size;

        // dgMatrix = I - tDelta * DF.
        const dgMatrix = new Array<number>(n * n).fill(0);
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                const index = c + n * r;
                const identity = (r === c ? 1 : 0);
                dgMatrix[index] = identity - this.mTDelta * dfMatrix[index];
            }
        }

        const dgInverse = inverse(dgMatrix, n);
        fVector = multiply(dgInverse, fVector);
        const tOut = tIn + this.mTDelta;
        const xOut = add(xIn, mul(this.mTDelta, fVector));
        return { tOut, xOut };
    }
}

// Compute the inverse of the N-by-N row-major matrix m using Gauss-Jordan
// elimination with partial pivoting. The zero matrix is returned when m is
// not invertible, which is the behavior of the upstream Inverse(GMatrix).
function inverse(m: readonly number[], n: number): number[] {
    // The augmented working copy: a[.] is reduced to the identity while
    // inv[.] is transformed from the identity to the inverse.
    const a = m.slice();
    const inv = new Array<number>(n * n).fill(0);
    for (let d = 0; d < n; ++d) {
        inv[d + n * d] = 1;
    }

    for (let col = 0; col < n; ++col) {
        // Select the pivot row, the one with the largest absolute value in
        // the current column.
        let pivotRow = col;
        let maxValue = Math.abs(a[col + n * col]);
        for (let row = col + 1; row < n; ++row) {
            const absValue = Math.abs(a[col + n * row]);
            if (absValue > maxValue) {
                maxValue = absValue;
                pivotRow = row;
            }
        }

        if (maxValue === 0) {
            // The matrix is not invertible.
            return new Array<number>(n * n).fill(0);
        }

        if (pivotRow !== col) {
            for (let c = 0; c < n; ++c) {
                let temp = a[c + n * col];
                a[c + n * col] = a[c + n * pivotRow];
                a[c + n * pivotRow] = temp;
                temp = inv[c + n * col];
                inv[c + n * col] = inv[c + n * pivotRow];
                inv[c + n * pivotRow] = temp;
            }
        }

        // Scale the pivot row so that the pivot element is 1.
        const invPivot = 1 / a[col + n * col];
        for (let c = 0; c < n; ++c) {
            a[c + n * col] *= invPivot;
            inv[c + n * col] *= invPivot;
        }

        // Eliminate the current column from the other rows.
        for (let row = 0; row < n; ++row) {
            if (row !== col) {
                const factor = a[col + n * row];
                if (factor !== 0) {
                    for (let c = 0; c < n; ++c) {
                        a[c + n * row] -= factor * a[c + n * col];
                        inv[c + n * row] -= factor * inv[c + n * col];
                    }
                }
            }
        }
    }

    return inv;
}

// Compute the product of the N-by-N row-major matrix m and the N-vector v.
function multiply(m: readonly number[], v: Vector): Vector {
    const n = v.size;
    const result = new Vector(n);
    for (let r = 0; r < n; ++r) {
        let sum = 0;
        for (let c = 0; c < n; ++c) {
            sum += m[c + n * r] * v.values[c];
        }
        result.values[r] = sum;
    }
    return result;
}
