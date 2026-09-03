// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) LCPSolver.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

import { logError } from './Logger.js';

// A class for solving the Linear Complementarity Problem (LCP)
//   w = q + M * z,  w^T * z = 0,  w >= 0,  z >= 0
// The vectors q, w and z are n-tuples and the matrix M is n-by-n. The inputs
// to solve(...) are q and M. The outputs are w and z, which are valid when
// the returned 'success' is true but are invalid when it is false.
//
// The algorithm is Lemke's complementary pivoting method with a symbolic
// perturbation of q (each q[r] is replaced by the polynomial
// q[r] + t^(r+1)) to avoid the degeneracies that arise when a q-term becomes
// zero during the iterations.
//
// Port notes:
// - Upstream has an abstract LCPSolverShared<T> plus two derived classes: one
//   with the dimension as a compile-time template parameter and one with the
//   dimension known only at run time. TypeScript has no compile-time
//   dimensions, so the port is a single class with a run-time dimension,
//   equivalent to upstream's LCPSolver<T>. The shared base is folded in.
// - Upstream's Variable::tuple is a T* pointing at &w[0] or &z[0]. The port
//   stores a reference to the w or z array directly, which is equivalent
//   because JavaScript arrays are references.
// - Upstream's mPoly is an array of T* pointing into the columns
//   [n+1, 2n+1] of each row of mAugmented. The port stores the corresponding
//   flat offsets into mAugmented, and the polynomial helpers take an array
//   plus an offset.
// - The nested enum LCPSolverShared<T>::Result is exported as the top-level
//   enum LCPSolverResult (export names must be globally unique).
// - GTE_THROW_ON_LCPSOLVER_ERRORS is ported as the static flag
//   LCPSolver.throwOnErrors (default false, matching upstream's
//   commented-out #define).
// - solve(...) returns { success, result, w, z } instead of writing through
//   output parameters.
// - Upstream is templated on T so that the solver can be instantiated with
//   an exact rational type (the header suggests BSRational<UIntegerAP32> for
//   diagnosing round-off failures) or with QFNumber, which is why there is a
//   second constructor taking the representations of zero and one. Following
//   PORTING.md, only the floating-point instantiation is ported: T is
//   'number', zero and one are the literals 0 and 1, and the extra
//   constructor has no purpose here and is omitted.
// - Upstream's dynamic-size LCPSolver<T> allows construction with n <= 0
//   (mDimension becomes 0 and no arrays are allocated), and Solve then
//   dereferences the null mPoly and mQMin pointers. The port reports
//   INVALID_INPUT for that case; see the guard in solve().

export enum LCPSolverResult {
    HAS_TRIVIAL_SOLUTION,
    HAS_NONTRIVIAL_SOLUTION,
    NO_SOLUTION,
    FAILED_TO_CONVERGE,
    INVALID_INPUT
}

export interface LCPSolverOutput {
    // Upstream's bool return value from Solve.
    success: boolean;
    result: LCPSolverResult;
    w: number[];
    z: number[];
}

// Bookkeeping of variables during the iterations of the solver. The name is
// either 'w' or 'z' and is used for human-readable debugging help. The
// 'index' is that for the original variables w[index] or z[index]. The
// 'complementary' index is the location of the complementary variable in
// mVarBasic[] or in mVarNonbasic[]. The 'tuple' is the w or z array, and is
// used to fill in the solution values (the variables are permuted during the
// pivoting algorithm).
interface LCPVariable {
    name: string;
    index: number;
    complementary: number;
    tuple: number[];
}

export class LCPSolver {
    // Port of GTE_THROW_ON_LCPSOLVER_ERRORS. When true, a failure to
    // converge throws instead of returning FAILED_TO_CONVERGE.
    static throwOnErrors = false;

    private mDimension: number;
    private mMaxIterations: number;
    private mNumIterations: number;

    private mVarBasic: LCPVariable[];
    private mVarNonbasic: LCPVariable[];
    private mNumCols: number;
    private mAugmented: number[];
    private mQMin: number[];
    private mMinRatio: number[];
    private mRatio: number[];
    // Offsets into mAugmented at which each row's perturbation polynomial
    // begins (upstream's T* mPoly[r]).
    private mPolyOffset: number[];

    // Construction. The member mMaxIterations is set by this call to the
    // default value n*n.
    constructor(n: number) {
        if (n > 0) {
            this.mDimension = n;
            this.mMaxIterations = n * n;
        } else {
            this.mDimension = 0;
            this.mMaxIterations = 0;
        }
        this.mNumIterations = 0;

        this.mVarBasic = [];
        this.mVarNonbasic = [];
        this.mNumCols = 0;
        this.mAugmented = [];
        this.mQMin = [];
        this.mMinRatio = [];
        this.mRatio = [];
        this.mPolyOffset = [];

        if (n > 0) {
            const np1 = n + 1;
            for (let i = 0; i < np1; ++i) {
                this.mVarBasic.push({ name: 'w', index: i, complementary: i, tuple: [] });
                this.mVarNonbasic.push({ name: 'z', index: i, complementary: i, tuple: [] });
            }
            this.mNumCols = 2 * np1;
            this.mAugmented = new Array<number>(2 * np1 * n).fill(0);
            this.mQMin = new Array<number>(np1).fill(0);
            this.mMinRatio = new Array<number>(np1).fill(0);
            this.mRatio = new Array<number>(np1).fill(0);
            this.mPolyOffset = new Array<number>(n).fill(0);
        }
    }

    getDimension(): number {
        return this.mDimension;
    }

    // Theoretically, when there is a solution the algorithm must converge in
    // a finite number of iterations. The number of iterations depends on the
    // problem at hand, but we need to guard against an infinite loop by
    // limiting the number. The implementation uses a maximum number of n*n
    // (chosen arbitrarily). You can set the number yourself, perhaps when a
    // call to solve fails -- increase the number of iterations and solve
    // again.
    setMaxIterations(maxIterations: number): void {
        this.mMaxIterations = (maxIterations > 0
            ? maxIterations : this.mDimension * this.mDimension);
    }

    getMaxIterations(): number {
        return this.mMaxIterations;
    }

    // Access the actual number of iterations used in a call to solve.
    getNumIterations(): number {
        return this.mNumIterations;
    }

    // The input q must have n elements and the input M must be an n-by-n
    // matrix stored in row-major order. The outputs w and z have n elements.
    solve(q: readonly number[], M: readonly number[]): LCPSolverOutput {
        const n = this.mDimension;
        const w = new Array<number>(Math.max(n, 0)).fill(0);
        const z = new Array<number>(Math.max(n, 0)).fill(0);

        // The 'n <= 0' guard is a port addition. Upstream's dynamic-size
        // LCPSolver<T> allows construction with n <= 0 (mDimension is set to
        // 0 and no arrays are allocated), but LCPSolverShared<T>::Solve then
        // evaluates Copy(mPoly[0], mQMin) on an empty vector, which
        // dereferences a null pointer. The port reports INVALID_INPUT
        // instead.
        if (n <= 0 || n > q.length || n * n > M.length) {
            return {
                success: false,
                result: LCPSolverResult.INVALID_INPUT,
                w,
                z
            };
        }

        // Perturb the q[r] constants to be polynomials of degree r+1
        // represented as an array of n+1 coefficients. The coefficient with
        // index r+1 is 1 and the coefficients with indices larger than r+1
        // are 0.
        for (let r = 0; r < n; ++r) {
            this.mPolyOffset[r] = this.augmentedIndex(r, n + 1);
            this.makeZero(this.mAugmented, this.mPolyOffset[r]);
            this.mAugmented[this.mPolyOffset[r]] = q[r];
            this.mAugmented[this.mPolyOffset[r] + r + 1] = 1;
        }

        // Determine whether there is the trivial solution w = z = 0.
        this.copyPoly(this.mAugmented, this.mPolyOffset[0], this.mQMin, 0);
        let basic = 0;
        for (let r = 1; r < n; ++r) {
            if (this.lessThan(this.mAugmented, this.mPolyOffset[r], this.mQMin, 0)) {
                this.copyPoly(this.mAugmented, this.mPolyOffset[r], this.mQMin, 0);
                basic = r;
            }
        }

        if (!this.lessThanZero(this.mQMin, 0)) {
            for (let r = 0; r < n; ++r) {
                w[r] = q[r];
                z[r] = 0;
            }
            return {
                success: true,
                result: LCPSolverResult.HAS_TRIVIAL_SOLUTION,
                w,
                z
            };
        }

        // Initialize the remainder of the augmented matrix with M and U.
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                this.mAugmented[this.augmentedIndex(r, c)] = M[c + n * r];
            }
            this.mAugmented[this.augmentedIndex(r, n)] = 1;
        }

        // Keep track of when the variables enter and exit the dictionary,
        // including where complementary variables are relocated.
        for (let i = 0; i <= n; ++i) {
            this.mVarBasic[i] = { name: 'w', index: i, complementary: i, tuple: w };
            this.mVarNonbasic[i] = { name: 'z', index: i, complementary: i, tuple: z };
        }

        // The augmented variable z[n] is the initial driving variable for
        // pivoting. The equation 'basic' is the one to solve for z[n] and
        // pivoting with w[basic]. The last column of M remains all 1-values
        // for this initial step, so no algebraic computations occur for
        // M[r][n].
        let driving = n;
        for (let r = 0; r < n; ++r) {
            if (r !== basic) {
                for (let c = 0; c < this.mNumCols; ++c) {
                    if (c !== n) {
                        this.mAugmented[this.augmentedIndex(r, c)] -=
                            this.mAugmented[this.augmentedIndex(basic, c)];
                    }
                }
            }
        }

        for (let c = 0; c < this.mNumCols; ++c) {
            if (c !== n) {
                this.mAugmented[this.augmentedIndex(basic, c)] =
                    -this.mAugmented[this.augmentedIndex(basic, c)];
            }
        }

        this.mNumIterations = 0;
        for (let i = 0; i < this.mMaxIterations; ++i, ++this.mNumIterations) {
            // The basic variable of equation 'basic' exited the dictionary,
            // so its complementary (nonbasic) variable must become the next
            // driving variable in order for it to enter the dictionary.
            const nextDriving = this.mVarBasic[basic].complementary;
            this.mVarNonbasic[nextDriving].complementary = driving;
            const swapTmp = this.mVarBasic[basic];
            this.mVarBasic[basic] = this.mVarNonbasic[driving];
            this.mVarNonbasic[driving] = swapTmp;
            if (this.mVarNonbasic[driving].index === n) {
                // The algorithm has converged.
                for (let r = 0; r < n; ++r) {
                    this.mVarBasic[r].tuple[this.mVarBasic[r].index] =
                        this.mAugmented[this.mPolyOffset[r]];
                }
                for (let c = 0; c <= n; ++c) {
                    const index = this.mVarNonbasic[c].index;
                    if (index < n) {
                        this.mVarNonbasic[c].tuple[index] = 0;
                    }
                }
                return {
                    success: true,
                    result: LCPSolverResult.HAS_NONTRIVIAL_SOLUTION,
                    w,
                    z
                };
            }

            // Determine the 'basic' equation for which the ratio
            // -q[r]/M(r,driving) is minimized among all equations r with
            // M(r,driving) < 0.
            driving = nextDriving;
            basic = -1;
            for (let r = 0; r < n; ++r) {
                if (this.mAugmented[this.augmentedIndex(r, driving)] < 0) {
                    const factor = -1 / this.mAugmented[this.augmentedIndex(r, driving)];
                    this.multiplyPoly(this.mAugmented, this.mPolyOffset[r], factor,
                        this.mRatio, 0);
                    if (basic === -1 || this.lessThan(this.mRatio, 0, this.mMinRatio, 0)) {
                        this.copyPoly(this.mRatio, 0, this.mMinRatio, 0);
                        basic = r;
                    }
                }
            }

            if (basic === -1) {
                // The coefficients of z[driving] in all the equations are
                // nonnegative, so the z[driving] variable cannot leave the
                // dictionary. There is no solution to the LCP.
                for (let r = 0; r < n; ++r) {
                    w[r] = 0;
                    z[r] = 0;
                }
                return {
                    success: false,
                    result: LCPSolverResult.NO_SOLUTION,
                    w,
                    z
                };
            }

            // Solve the basic equation so that z[driving] enters the
            // dictionary and w[basic] exits the dictionary.
            const invDenom = 1 / this.mAugmented[this.augmentedIndex(basic, driving)];
            for (let r = 0; r < n; ++r) {
                const arDriving = this.mAugmented[this.augmentedIndex(r, driving)];
                if (r !== basic && arDriving !== 0) {
                    const multiplier = arDriving * invDenom;
                    for (let c = 0; c < this.mNumCols; ++c) {
                        if (c !== driving) {
                            this.mAugmented[this.augmentedIndex(r, c)] -=
                                this.mAugmented[this.augmentedIndex(basic, c)] * multiplier;
                        } else {
                            this.mAugmented[this.augmentedIndex(r, driving)] = multiplier;
                        }
                    }
                }
            }

            for (let c = 0; c < this.mNumCols; ++c) {
                if (c !== driving) {
                    this.mAugmented[this.augmentedIndex(basic, c)] =
                        -this.mAugmented[this.augmentedIndex(basic, c)] * invDenom;
                } else {
                    this.mAugmented[this.augmentedIndex(basic, driving)] = invDenom;
                }
            }
        }

        // Numerical round-off errors can cause the Lemke algorithm not to
        // converge. In particular, the code above has a test
        //   if (mAugmented[r][driving] < 0) { ... }
        // to determine the 'basic' equation with which to pivot. It is
        // possible that theoretically mAugmented[r][driving] is zero but
        // rounding errors cause it to be slightly negative. If theoretically
        // all mAugmented[r][driving] >= 0, there is no solution to the LCP.
        // With the rounding errors, if the algorithm fails to converge within
        // the specified number of iterations, NO_SOLUTION is hopefully the
        // correct result. It is also possible that the rounding errors lead
        // to a NO_SOLUTION (returned from inside the loop) when in fact there
        // is a solution. When the LCP solver is used by intersection testing
        // algorithms, the hope is that misclassifications occur only when the
        // two objects are nearly in tangential contact.
        if (LCPSolver.throwOnErrors) {
            logError('LCPSolver.solve failed to converge.');
        }
        return {
            success: false,
            result: LCPSolverResult.FAILED_TO_CONVERGE,
            w,
            z
        };
    }

    // Access mAugmented as a 2-dimensional array.
    private augmentedIndex(row: number, col: number): number {
        return col + this.mNumCols * row;
    }

    // Support for polynomials with n+1 coefficients and degree no larger
    // than n. Each polynomial is a run of n+1 entries of 'data' starting at
    // 'offset'.
    private makeZero(data: number[], offset: number): void {
        for (let i = 0; i <= this.mDimension; ++i) {
            data[offset + i] = 0;
        }
    }

    private copyPoly(src: readonly number[], srcOffset: number,
        dst: number[], dstOffset: number): void {
        for (let i = 0; i <= this.mDimension; ++i) {
            dst[dstOffset + i] = src[srcOffset + i];
        }
    }

    private lessThan(poly0: readonly number[], offset0: number,
        poly1: readonly number[], offset1: number): boolean {
        for (let i = 0; i <= this.mDimension; ++i) {
            if (poly0[offset0 + i] < poly1[offset1 + i]) {
                return true;
            }
            if (poly0[offset0 + i] > poly1[offset1 + i]) {
                return false;
            }
        }
        return false;
    }

    private lessThanZero(poly: readonly number[], offset: number): boolean {
        for (let i = 0; i <= this.mDimension; ++i) {
            if (poly[offset + i] < 0) {
                return true;
            }
            if (poly[offset + i] > 0) {
                return false;
            }
        }
        return false;
    }

    private multiplyPoly(poly: readonly number[], offset: number, scalar: number,
        product: number[], productOffset: number): void {
        for (let i = 0; i <= this.mDimension; ++i) {
            product[productOffset + i] = poly[offset + i] * scalar;
        }
    }
}
