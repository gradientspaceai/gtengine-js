// gtengine-js: TypeScript port of Geometric Tools Engine (GTE)
// LevenbergMarquardtMinimizer.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// See GaussNewtonMinimizer.ts for a formulation of the minimization problem
// and how Levenberg-Marquardt relates to Gauss-Newton. The difference is
// that the normal-equations matrix J^T*J is modified by adding
// lambda*average(diagonal(J^T*J)) to its diagonal, which biases the step
// toward gradient descent when lambda is large and toward Gauss-Newton when
// lambda is small.
//
// Port notes:
// - Upstream's GVector<T>/GMatrix<T> (run-time sizes) are the port's Vector
//   and Matrix.
// - The two constructors become the static factories fromJFunction and
//   fromJPlusFunction (see GaussNewtonMinimizer.ts for the rationale).
// - operator() becomes minimize(); the nested Result struct becomes the
//   exported LevenbergMarquardtMinimizerResult interface. The
//   std::pair<bool,bool> returned by DoIteration becomes an object with the
//   named fields 'stop' and 'errorReduced'.
// - Upstream bug (fixed here): DoIteration recomputes the Jacobian at
//   pCurrent but reuses the member mF, which at that point holds F evaluated
//   at the *previous* candidate pNext (the tail of the previous DoIteration
//   call), not F(pCurrent). Whenever the inner lambda-adjustment loop calls
//   DoIteration more than once for the same pCurrent, the right-hand side
//   -J^T(pCurrent)*F is therefore built from a stale residual vector, which
//   corrupts the step. The port tracks the point at which mF was last
//   evaluated and re-evaluates F(pCurrent) when it is stale, so no extra
//   function evaluations occur on the common path. The JPlus path is
//   unaffected because that callback computes both J^T*J and -J^T*F from
//   pCurrent itself.

import { CholeskyDecomposition } from './CholeskyDecomposition.js';
import { logAssert } from './Logger.js';
import { Matrix, multiplyATB, mulMatrix } from './Matrix.js';
import { Vector, add, dot, length, negate } from './Vector.js';

// F: given p (numPDimensions), fill f (numFDimensions) with F(p).
export type LevenbergMarquardtFFunction = (p: Vector, f: Vector) => void;

// J: given p, fill j (numFDimensions-by-numPDimensions) with the Jacobian.
export type LevenbergMarquardtJFunction = (p: Vector, j: Matrix) => void;

// J plus: given p, fill jtj (numPDimensions-by-numPDimensions) with
// J^T(p)*J(p) and negJTF (numPDimensions) with -J^T(p)*F(p).
export type LevenbergMarquardtJPlusFunction =
    (p: Vector, jtj: Matrix, negJTF: Vector) => void;

export interface LevenbergMarquardtMinimizerResult {
    minLocation: Vector;
    minError: number;
    minErrorDifference: number;
    minUpdateLength: number;
    numIterations: number;
    numAdjustments: number;
    converged: boolean;
}

// The 'stop' flag is true when the linear system cannot be solved
// (result.converged is false in this case) or when the error is reduced to
// within the tolerances specified by the caller (result.converged is true in
// this case). The 'errorReduced' flag is true when the error is reduced.
interface DoIterationStatus {
    stop: boolean;
    errorReduced: boolean;
}

export class LevenbergMarquardtMinimizer {
    private readonly mNumPDimensions: number;
    private readonly mNumFDimensions: number;
    private readonly mFFunction: LevenbergMarquardtFFunction;
    private readonly mJFunction: LevenbergMarquardtJFunction | null;
    private readonly mJPlusFunction: LevenbergMarquardtJPlusFunction | null;

    // Storage for F(p), J(p), J^T(p)*J(p) and -J^T(p)*F(p) during the
    // iterations.
    private mF: Vector;
    private readonly mJ: Matrix;
    private mJTJ: Matrix;
    private mNegJTF: Vector;

    // The domain point at which mF was last evaluated. See the upstream-bug
    // port note at the top of the file.
    private mFPoint: Vector | null;

    private readonly mDecomposer: CholeskyDecomposition;
    private readonly mUseJFunction: boolean;

    private constructor(numPDimensions: number, numFDimensions: number,
        fFunction: LevenbergMarquardtFFunction,
        jFunction: LevenbergMarquardtJFunction | null,
        jPlusFunction: LevenbergMarquardtJPlusFunction | null) {
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
        this.mFPoint = null;
        this.mDecomposer = new CholeskyDecomposition(numPDimensions);
        this.mUseJFunction = (jFunction !== null);
    }

    // Create the minimizer that computes F(p) and J(p) directly.
    static fromJFunction(numPDimensions: number, numFDimensions: number,
        fFunction: LevenbergMarquardtFFunction,
        jFunction: LevenbergMarquardtJFunction): LevenbergMarquardtMinimizer {
        return new LevenbergMarquardtMinimizer(numPDimensions, numFDimensions,
            fFunction, jFunction, null);
    }

    // Create the minimizer that computes J^T(p)*J(p) and -J^T(p)*F(p).
    static fromJPlusFunction(numPDimensions: number, numFDimensions: number,
        fFunction: LevenbergMarquardtFFunction,
        jPlusFunction: LevenbergMarquardtJPlusFunction):
        LevenbergMarquardtMinimizer {
        return new LevenbergMarquardtMinimizer(numPDimensions, numFDimensions,
            fFunction, null, jPlusFunction);
    }

    getNumPDimensions(): number {
        return this.mNumPDimensions;
    }

    getNumFDimensions(): number {
        return this.mNumFDimensions;
    }

    // The lambda is positive, the multiplier is positive, and the initial
    // guess for the p-parameter is p0. Typical choices are
    // lambdaFactor = 0.001 and lambdaAdjust = 10.
    minimize(p0: Vector, maxIterations: number,
        updateLengthTolerance: number, errorDifferenceTolerance: number,
        lambdaFactor: number, lambdaAdjust: number,
        maxAdjustments: number): LevenbergMarquardtMinimizerResult {
        const result: LevenbergMarquardtMinimizerResult = {
            minLocation: p0.clone(),
            minError: Number.MAX_VALUE,
            minErrorDifference: Number.MAX_VALUE,
            minUpdateLength: 0,
            numIterations: 0,
            numAdjustments: 0,
            converged: false
        };

        // As a simple precaution, ensure that the lambda inputs are valid.
        // If invalid, fall back to Gauss-Newton iteration.
        if (lambdaFactor <= 0 || lambdaAdjust <= 0) {
            maxAdjustments = 1;
            lambdaFactor = 0;
            lambdaAdjust = 1;
        }

        // As a simple precaution, ensure the tolerances are nonnegative.
        updateLengthTolerance = Math.max(updateLengthTolerance, 0);
        errorDifferenceTolerance = Math.max(errorDifferenceTolerance, 0);

        // Compute the initial error.
        this.mFFunction(p0, this.mF);
        this.mFPoint = p0;
        result.minError = dot(this.mF, this.mF);

        // Do the Levenberg-Marquardt iterations.
        let pCurrent = p0;
        for (result.numIterations = 1;
            result.numIterations <= maxIterations;
            ++result.numIterations) {
            let status: DoIterationStatus = { stop: false, errorReduced: false };
            let pNext = new Vector(this.mNumPDimensions);
            for (result.numAdjustments = 0;
                result.numAdjustments < maxAdjustments;
                ++result.numAdjustments) {
                const iterate = this.doIteration(pCurrent, lambdaFactor,
                    updateLengthTolerance, errorDifferenceTolerance, result);
                status = iterate.status;
                pNext = iterate.pNext;
                if (status.stop) {
                    // Either the Cholesky decomposition failed or the
                    // iterates converged within tolerance.
                    return result;
                }

                if (status.errorReduced) {
                    // The error has been reduced but we have not yet
                    // converged within tolerance.
                    break;
                }

                lambdaFactor *= lambdaAdjust;
            }

            if (result.numAdjustments < maxAdjustments) {
                // The current value of lambda led us to an update that
                // reduced the error, but the error is not yet small enough to
                // conclude we converged. Reduce lambda for the next
                // outer-loop iteration.
                lambdaFactor /= lambdaAdjust;
            } else {
                // All lambdas tried during the inner-loop iteration did not
                // lead to a reduced error. If we do nothing here, the next
                // inner-loop iteration will continue to multiply lambda,
                // risking eventual floating-point overflow. To avoid this,
                // fall back to a Gauss-Newton iterate.
                const iterate = this.doIteration(pCurrent, lambdaFactor,
                    updateLengthTolerance, errorDifferenceTolerance, result);
                status = iterate.status;
                pNext = iterate.pNext;
                if (status.stop) {
                    // Either the Cholesky decomposition failed or the
                    // iterates converged within tolerance.
                    return result;
                }
            }

            pCurrent = pNext;
        }

        return result;
    }

    private computeLinearSystemInputs(pCurrent: Vector, lambda: number): void {
        if (this.mUseJFunction) {
            // Ensure mF stores F(pCurrent). See the upstream-bug port note.
            if (this.mFPoint !== pCurrent) {
                this.mFFunction(pCurrent, this.mF);
                this.mFPoint = pCurrent;
            }
            (this.mJFunction as LevenbergMarquardtJFunction)(pCurrent, this.mJ);
            this.mJTJ = multiplyATB(this.mJ, this.mJ);
            this.mNegJTF = negate(mulMatrix(this.mF, this.mJ));
        } else {
            (this.mJPlusFunction as LevenbergMarquardtJPlusFunction)(
                pCurrent, this.mJTJ, this.mNegJTF);
        }

        let diagonalSum = 0;
        for (let i = 0; i < this.mNumPDimensions; ++i) {
            diagonalSum += this.mJTJ.get(i, i);
        }

        const diagonalAdjust = lambda * diagonalSum / this.mNumPDimensions;
        for (let i = 0; i < this.mNumPDimensions; ++i) {
            this.mJTJ.set(i, i, this.mJTJ.get(i, i) + diagonalAdjust);
        }
    }

    private doIteration(pCurrent: Vector, lambdaFactor: number,
        updateLengthTolerance: number, errorDifferenceTolerance: number,
        result: LevenbergMarquardtMinimizerResult):
        { status: DoIterationStatus, pNext: Vector } {
        this.computeLinearSystemInputs(pCurrent, lambdaFactor);
        if (!this.mDecomposer.factor(this.mJTJ)) {
            // TODO (upstream): The matrix J^T*J is positive semi-definite, so
            // the failure can occur when it has a zero eigenvalue, in which
            // case it is not invertible. Generate an iterate anyway, perhaps
            // using gradient descent?
            return {
                status: { stop: true, errorReduced: false },
                pNext: pCurrent
            };
        }
        this.mDecomposer.solveLower(this.mJTJ, this.mNegJTF);
        this.mDecomposer.solveUpper(this.mJTJ, this.mNegJTF);

        const pNext = add(pCurrent, this.mNegJTF);
        this.mFFunction(pNext, this.mF);
        this.mFPoint = pNext;
        const error = dot(this.mF, this.mF);
        if (error < result.minError) {
            result.minErrorDifference = result.minError - error;
            result.minUpdateLength = length(this.mNegJTF);
            result.minLocation = pNext.clone();
            result.minError = error;
            if (result.minErrorDifference <= errorDifferenceTolerance
                || result.minUpdateLength <= updateLengthTolerance) {
                result.converged = true;
                return { status: { stop: true, errorReduced: true }, pNext };
            } else {
                return { status: { stop: false, errorReduced: true }, pNext };
            }
        } else {
            return { status: { stop: false, errorReduced: false }, pNext };
        }
    }
}
