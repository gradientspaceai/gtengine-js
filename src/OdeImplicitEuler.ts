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
// GMatrix<Real>). The port has a single Vector class and a single Matrix
// class whose dimensions are run-time values, so the solver is a concrete
// class derived from OdeSolver<Vector> and DF returns a Matrix. The upstream
// call Inverse(dgMatrix) is the port's inverse() from Matrix.ts, which - as
// upstream - is the full-pivoting Gaussian elimination of
// GaussianElimination.h and produces the zero matrix for a noninvertible
// input rather than throwing.

import { Matrix, inverse, mulMatrix } from './Matrix.js';
import { OdeSolver, type OdeFunction } from './OdeSolver.js';
import { Vector, add, mul } from './Vector.js';

// The port of 'std::function<TMatrix(Real, TVector const&)>'. DF(t,x) is the
// N-by-N matrix whose entry (r,c) is the derivative of F[r] with respect to
// x[c].
export type OdeDerivativeFunction = (t: number, x: Vector) => Matrix;

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
        const dgMatrix = new Matrix(n, n);
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                const identity = (r === c ? 1 : 0);
                dgMatrix.set(r, c,
                    identity - this.mTDelta * dfMatrix.get(r, c));
            }
        }

        const dgInverse = inverse(dgMatrix).inverse;
        fVector = mulMatrix(dgInverse, fVector);
        const tOut = tIn + this.mTDelta;
        const xOut = add(xIn, mul(this.mTDelta, fVector));
        return { tOut, xOut };
    }
}
