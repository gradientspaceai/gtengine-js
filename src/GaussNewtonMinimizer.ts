// gtengine-js: TypeScript port of Geometric Tools Engine (GTE)
// GaussNewtonMinimizer.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Let F(p) = (F_{0}(p), F_{1}(p), ..., F_{n-1}(p)) be a vector-valued
// function of the parameters p = (p_{0}, p_{1}, ..., p_{m-1}). The nonlinear
// least-squares problem is to minimize the real-valued error function
// E(p) = |F(p)|^2, which is the squared length of F(p).
//
// Let J = dF/dp = [dF_{r}/dp_{c}] denote the Jacobian matrix, which is the
// matrix of first-order partial derivatives of F. The matrix has n rows and
// m columns, and the indexing (r,c) refers to row r and column c. A
// first-order approximation is F(p + d) = F(p) + J(p)d, where d is an m-by-1
// vector with small length. Consequently, an approximation to E is
// E(p + d) = |F(p + d)|^2 = |F(p) + J(p)d|^2. The goal is to choose d to
// minimize |F(p) + J(p)d|^2 and, hopefully, with E(p + d) < E(p). Choosing
// an initial p_{0}, the hope is that the algorithm generates a sequence
// p_{i} for which E(p_{i+1}) < E(p_{i}) and, in the limit, E(p_{j})
// approaches the global minimum of E. The algorithm is referred to as
// Gauss-Newton iteration. If E does not decrease for a step of the
// algorithm, one can modify the algorithm to the Levenberg-Marquardt
// iteration. See LevenbergMarquardtMinimizer.ts.
//
// For a single Gauss-Newton iteration, we need to choose d to minimize
// |F(p) + J(p)d|^2 where p is fixed. This is a linear least-squares problem
// which can be formulated using the normal equations
// (J^T(p)*J(p))*d = -J^T(p)*F(p). The matrix J^T*J is positive semidefinite.
// If it is invertible, then d = -(J^T(p)*J(p))^{-1}*F(p). If it is not
// invertible, some other algorithm must be used to choose d; one option is
// to use gradient descent for the step. A Cholesky decomposition is used to
// solve the linear system.
//
// Although an implementation can allow the caller to pass an array of
// functions F_{i}(p) and an array of derivatives dF_{r}/dp_{c}, some
// applications might involve a very large n that precludes storing all the
// computed Jacobian matrix entries because of excessive memory requirements.
// In such an application, it is better to compute instead the entries of the
// m-by-m matrix J^T*J and the m-by-1 vector J^T*F. Typically, m is small, so
// the memory requirements are not excessive. Both approaches are supported.
//
// Port notes:
// - Upstream's GVector<T>/GMatrix<T> (run-time sizes) are the port's Vector
//   and Matrix, which already carry run-time dimensions.
// - The two constructors (one taking a J function, the other taking a
//   "J plus" function that produces J^T*J and -J^T*F directly) become the
//   static factories fromJFunction/fromJPlusFunction, since the two callback
//   types are not distinguishable overloads in TypeScript.
// - operator() becomes minimize(); the nested Result struct becomes the
//   exported GaussNewtonMinimizerResult interface.
// - The callbacks write into the output objects handed to them, matching the
//   upstream reference-parameter convention.

import { CholeskyDecomposition } from './CholeskyDecomposition.js';
import { logAssert } from './Logger.js';
import { Matrix, multiplyATB, mulMatrix } from './Matrix.js';
import { Vector, add, dot, length, negate } from './Vector.js';

// F: given p (numPDimensions), fill f (numFDimensions) with F(p).
export type GaussNewtonFFunction = (p: Vector, f: Vector) => void;

// J: given p, fill j (numFDimensions-by-numPDimensions) with the Jacobian.
export type GaussNewtonJFunction = (p: Vector, j: Matrix) => void;

// J plus: given p, fill jtj (numPDimensions-by-numPDimensions) with
// J^T(p)*J(p) and negJTF (numPDimensions) with -J^T(p)*F(p).
export type GaussNewtonJPlusFunction =
    (p: Vector, jtj: Matrix, negJTF: Vector) => void;

export interface GaussNewtonMinimizerResult {
    minLocation: Vector;
    minError: number;
    minErrorDifference: number;
    minUpdateLength: number;
    numIterations: number;
    converged: boolean;
}

export class GaussNewtonMinimizer {
    private readonly mNumPDimensions: number;
    private readonly mNumFDimensions: number;
    private readonly mFFunction: GaussNewtonFFunction;
    private readonly mJFunction: GaussNewtonJFunction | null;
    private readonly mJPlusFunction: GaussNewtonJPlusFunction | null;

    // Storage for F(p), J(p), J^T(p)*J(p) and -J^T(p)*F(p) during the
    // iterations.
    private mF: Vector;
    private readonly mJ: Matrix;
    private mJTJ: Matrix;
    private mNegJTF: Vector;

    private readonly mDecomposer: CholeskyDecomposition;
    private readonly mUseJFunction: boolean;

    private constructor(numPDimensions: number, numFDimensions: number,
        fFunction: GaussNewtonFFunction,
        jFunction: GaussNewtonJFunction | null,
        jPlusFunction: GaussNewtonJPlusFunction | null) {
        logAssert(numPDimensions > 0 && numFDimensions > 0,
            'Invalid dimensions.');
        this.mNumPDimensions = numPDimensions;
        this.mNumFDimensions = numFDimensions;
        this.mFFunction = fFunction;
        this.mJFunction = jFunction;
        this.mJPlusFunction = jPlusFunction;
        this.mF = new Vector(numFDimensions);
        this.mJ = new Matrix(numFDimensions, numPDimensions);
        this.mJTJ = new Matrix(numPDimensions, numPDimensions);
        this.mNegJTF = new Vector(numPDimensions);
        this.mDecomposer = new CholeskyDecomposition(numPDimensions);
        this.mUseJFunction = (jFunction !== null);
    }

    // Create the minimizer that computes F(p) and J(p) directly.
    static fromJFunction(numPDimensions: number, numFDimensions: number,
        fFunction: GaussNewtonFFunction,
        jFunction: GaussNewtonJFunction): GaussNewtonMinimizer {
        return new GaussNewtonMinimizer(numPDimensions, numFDimensions,
            fFunction, jFunction, null);
    }

    // Create the minimizer that computes J^T(p)*J(p) and -J^T(p)*F(p).
    static fromJPlusFunction(numPDimensions: number, numFDimensions: number,
        fFunction: GaussNewtonFFunction,
        jPlusFunction: GaussNewtonJPlusFunction): GaussNewtonMinimizer {
        return new GaussNewtonMinimizer(numPDimensions, numFDimensions,
            fFunction, null, jPlusFunction);
    }

    getNumPDimensions(): number {
        return this.mNumPDimensions;
    }

    getNumFDimensions(): number {
        return this.mNumFDimensions;
    }

    minimize(p0: Vector, maxIterations: number,
        updateLengthTolerance: number,
        errorDifferenceTolerance: number): GaussNewtonMinimizerResult {
        const result: GaussNewtonMinimizerResult = {
            minLocation: p0.clone(),
            minError: Number.MAX_VALUE,
            minErrorDifference: Number.MAX_VALUE,
            minUpdateLength: 0,
            numIterations: 0,
            converged: false
        };

        // As a simple precaution, ensure the tolerances are nonnegative.
        updateLengthTolerance = Math.max(updateLengthTolerance, 0);
        errorDifferenceTolerance = Math.max(errorDifferenceTolerance, 0);

        // Compute the initial error.
        this.mFFunction(p0, this.mF);
        result.minError = dot(this.mF, this.mF);

        // Do the Gauss-Newton iterations.
        let pCurrent = p0.clone();
        for (result.numIterations = 1;
            result.numIterations <= maxIterations;
            ++result.numIterations) {
            this.computeLinearSystemInputs(pCurrent);
            if (!this.mDecomposer.factor(this.mJTJ)) {
                // TODO (upstream): The matrix J^T*J is positive
                // semi-definite, so the failure can occur when it has a zero
                // eigenvalue, in which case it is not invertible. Generate an
                // iterate anyway, perhaps using gradient descent?
                return result;
            }
            this.mDecomposer.solveLower(this.mJTJ, this.mNegJTF);
            this.mDecomposer.solveUpper(this.mJTJ, this.mNegJTF);

            const pNext = add(pCurrent, this.mNegJTF);
            this.mFFunction(pNext, this.mF);
            const error = dot(this.mF, this.mF);
            if (error < result.minError) {
                result.minErrorDifference = result.minError - error;
                result.minUpdateLength = length(this.mNegJTF);
                result.minLocation = pNext.clone();
                result.minError = error;
                if (result.minErrorDifference <= errorDifferenceTolerance
                    || result.minUpdateLength <= updateLengthTolerance) {
                    result.converged = true;
                    return result;
                }
            }

            pCurrent = pNext;
        }

        return result;
    }

    private computeLinearSystemInputs(pCurrent: Vector): void {
        if (this.mUseJFunction) {
            (this.mJFunction as GaussNewtonJFunction)(pCurrent, this.mJ);
            this.mJTJ = multiplyATB(this.mJ, this.mJ);
            this.mNegJTF = negate(mulMatrix(this.mF, this.mJ));
        } else {
            (this.mJPlusFunction as GaussNewtonJPlusFunction)(
                pCurrent, this.mJTJ, this.mNegJTF);
        }
    }
}
