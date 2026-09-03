// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Hyperellipsoid.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A hyperellipsoid has center K; axis directions U[0] through U[N-1], all
// unit-length vectors; and extents e[0] through e[N-1], all positive numbers.
// A point X = K + sum_{d=0}^{N-1} y[d]*U[d] is on the hyperellipsoid whenever
// sum_{d=0}^{N-1} (y[d]/e[d])^2 = 1. An algebraic representation for the
// hyperellipsoid is (X-K)^T * M * (X-K) = 1, where M is the NxN symmetric
// matrix M = sum_{d=0}^{N-1} U[d]*U[d]^T/e[d]^2, where the superscript T
// denotes transpose. Observe that U[i]*U[i]^T is a matrix, not a scalar dot
// product. The hyperellipsoid is also represented by a quadratic equation
// 0 = C + B^T*X + X^T*A*X, where C is a scalar, B is an Nx1 vector, and A is
// an NxN symmetric matrix with positive eigenvalues. The coefficients can be
// stored from lowest degree to highest degree,
//   C = k[0]
//   B = k[1], ..., k[N]
//   A = k[N+1], ..., k[(N+1)(N+2)/2 - 1]
// where the A-coefficients are the upper-triangular elements of A listed in
// row-major order. For N = 2, X = (x[0],x[1]) and
//   0 = k[0] +
//       k[1]*x[0] + k[2]*x[1] +
//       k[3]*x[0]*x[0] + k[4]*x[0]*x[1]
//                      + k[5]*x[1]*x[1]
// For N = 3, X = (x[0],x[1],x[2]) and
//   0 = k[0] +
//       k[1]*x[0] + k[2]*x[1] + k[3]*x[2] +
//       k[4]*x[0]*x[0] + k[5]*x[0]*x[1] + k[6]*x[0]*x[2] +
//                      + k[7]*x[1]*x[1] + k[8]*x[1]*x[2] +
//                                       + k[9]*x[2]*x[2]
// This equation can be factored to the form (X-K)^T * M * (X-K) = 1, where
// K = -A^{-1}*B/2, M = A/(B^T*A^{-1}*B/4-C).
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (runtime dimension, 'new Hyperellipsoid(n)' for the default
// constructor, named static factories that copy their arguments, comparison
// methods). The upstream overloads of ToCoefficients/FromCoefficients are
// distinguished by name: toCoefficients()/fromCoefficients() use the packed
// coefficient array, toCoefficientsABC()/fromCoefficientsABC() use the
// (A, B, C) triple. Out-parameters become return values.

import { logAssert } from './Logger.js';
import {
    Matrix, addMatrix, inverse, mulMatrix, outerProduct
} from './Matrix.js';
import { SymmetricEigensolver } from './SymmetricEigensolver.js';
import { Vector, div, dot, mul } from './Vector.js';

// The number of coefficients (N+1)*(N+2)/2 of the quadratic equation for a
// hyperellipsoid of dimension N.
export function hyperellipsoidNumCoefficients(n: number): number {
    return ((n + 1) * (n + 2)) / 2;
}

// The port of the C++ lexicographic comparison of std::array<Vector,N>.
function axisLess(a0: readonly Vector[], a1: readonly Vector[]): boolean {
    const n = Math.min(a0.length, a1.length);
    for (let d = 0; d < n; ++d) {
        if (a0[d].lessThan(a1[d])) {
            return true;
        }
        if (a1[d].lessThan(a0[d])) {
            return false;
        }
    }
    return a0.length < a1.length;
}

function axisEquals(a0: readonly Vector[], a1: readonly Vector[]): boolean {
    if (a0.length !== a1.length) {
        return false;
    }
    for (let d = 0; d < a0.length; ++d) {
        if (!a0[d].equals(a1[d])) {
            return false;
        }
    }
    return true;
}

export class Hyperellipsoid {
    // Public member access.
    center: Vector;
    axis: Vector[];
    extent: Vector;

    // The port of the default constructor, which sets the center to
    // (0,...,0), the axes to Unit(d) and all extents to 1. The dimension N
    // of the C++ template is a constructor argument here.
    constructor(n: number) {
        logAssert(n >= 2, 'Invalid dimension.');

        this.center = new Vector(n);
        this.axis = new Array<Vector>(n);
        this.extent = new Vector(n);
        for (let d = 0; d < n; ++d) {
            const u = new Vector(n);
            u.makeUnit(d);
            this.axis[d] = u;
            this.extent.values[d] = 1;
        }
    }

    // The port of 'Hyperellipsoid(inCenter, inAxis, inExtent)'. The vectors
    // are copied, matching C++ value semantics.
    static fromCenterAxisExtent(inCenter: Vector, inAxis: readonly Vector[],
        inExtent: Vector): Hyperellipsoid {
        const n = inCenter.size;
        logAssert(inAxis.length === n && inExtent.size === n,
            'Hyperellipsoid: mismatched sizes.');
        const hyperellipsoid = new Hyperellipsoid(n);
        hyperellipsoid.center = inCenter.clone();
        hyperellipsoid.axis = inAxis.map(u => u.clone());
        hyperellipsoid.extent = inExtent.clone();
        return hyperellipsoid;
    }

    // The dimension N of the hyperellipsoid.
    get dimension(): number {
        return this.center.size;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Hyperellipsoid {
        return Hyperellipsoid.fromCenterAxisExtent(this.center, this.axis,
            this.extent);
    }

    // Compute M = sum_{d=0}^{N-1} U[d]*U[d]^T/e[d]^2.
    getM(): Matrix {
        const n = this.dimension;
        let m = Matrix.zero(n, n);
        for (let d = 0; d < n; ++d) {
            const ratio = div(this.axis[d], this.extent.values[d]);
            m = addMatrix(m, outerProduct(ratio, ratio));
        }
        return m;
    }

    // Compute M^{-1} = sum_{d=0}^{N-1} U[d]*U[d]^T*e[d]^2.
    getMInverse(): Matrix {
        const n = this.dimension;
        let mInverse = Matrix.zero(n, n);
        for (let d = 0; d < n; ++d) {
            const product = mul(this.axis[d], this.extent.values[d]);
            mInverse = addMatrix(mInverse, outerProduct(product, product));
        }
        return mInverse;
    }

    // Construct the coefficients in the quadratic equation that represents
    // the hyperellipsoid.
    toCoefficients(): number[] {
        const n = this.dimension;
        const numCoefficients = hyperellipsoidNumCoefficients(n);
        const { A, B, C } = this.toCoefficientsABC();
        const coeff = convertABCToCoeff(A, B, C);

        // For numerical robustness, divide the coefficients by the quadratic
        // coefficient of largest magnitude. The resulting coefficients are in
        // [-1,1]. The number of coefficients is (N+1)*(N+2)/2. Of these, 1 is
        // the constant term, N are the linear terms, and
        // (N+1)*(N+2)/2 - N - 1 = N(N+1)/2 are the quadratic terms. The
        // i-values in coeff[i] for the quadratic terms satisfy
        // N + 1 <= i < (N+1)*(N+2)/2.
        //
        // Upstream also tracks the index 'maxIndex' of the largest quadratic
        // coefficient, but never uses it; the dead store is dropped here.
        let maxValue = 0;
        for (let i = n + 1; i < numCoefficients; ++i) {
            const absValue = Math.abs(coeff[i]);
            if (absValue > maxValue) {
                maxValue = absValue;
            }
        }

        for (let i = 0; i < numCoefficients; ++i) {
            coeff[i] /= maxValue;
        }

        return coeff;
    }

    // The port of 'ToCoefficients(Matrix& A, Vector& B, Real& C)'.
    toCoefficientsABC(): { A: Matrix, B: Vector, C: number } {
        const A = this.getM();
        const product = mulMatrix(A, this.center) as Vector;
        const B = mul(product, -2);
        const C = dot(this.center, product) - 1;
        return { A, B, C };
    }

    // Construct the center, U[i] and e[i] from the equation. The return
    // value is 'true' if and only if the input coefficients represent a
    // hyperellipsoid. If the function returns 'false', the hyperellipsoid
    // data members are undefined.
    fromCoefficients(coeff: readonly number[]): boolean {
        const n = this.dimension;
        logAssert(coeff.length >= hyperellipsoidNumCoefficients(n),
            'Hyperellipsoid: invalid number of coefficients.');
        const { A, B, C } = convertCoeffToABC(n, coeff);
        return this.fromCoefficientsABC(A, B, C);
    }

    // The port of 'FromCoefficients(Matrix const& A, Vector const& B, Real C)'.
    fromCoefficientsABC(A: Matrix, B: Vector, C: number): boolean {
        const n = this.dimension;
        logAssert(A.numRows === n && A.numCols === n && B.size === n,
            'Hyperellipsoid: mismatched sizes.');

        // Compute the center K = -A^{-1}*B/2.
        const { inverse: invA, invertible } = inverse(A);
        if (!invertible) {
            return false;
        }

        this.center = mul(mulMatrix(invA, B) as Vector, -0.5);

        // Compute B^T*A^{-1}*B/4 - C = K^T*A*K - C = -K^T*B/2 - C.
        const rightSide = -0.5 * dot(this.center, B) - C;
        if (rightSide === 0) {
            return false;
        }

        // Compute M = A/(K^T*A*K - C).
        const invRightSide = 1 / rightSide;
        const M = mulMatrix(A, invRightSide) as Matrix;

        // Factor into M = R*D*R^T. M is symmetric, so it does not matter
        // whether the matrix is stored in row-major or column-major order;
        // they are equivalent. The output R, however, is in row-major order.
        const es = new SymmetricEigensolver(n, 32);
        es.solve(M.values, +1);  // diagonal[i] are nondecreasing
        const diagonal = es.getEigenvalues();
        const rotation = Matrix.fromArray(n, n, es.getEigenvectors());
        if (es.getEigenvectorMatrixType() === 0) {
            const negLast = mul(rotation.getCol(n - 1), -1);
            rotation.setCol(n - 1, negLast);
        }

        for (let d = 0; d < n; ++d) {
            if (diagonal[d] <= 0) {
                return false;
            }

            this.extent.values[d] = 1 / Math.sqrt(diagonal[d]);
            this.axis[d] = rotation.getCol(d);
        }

        return true;
    }

    // Comparisons to support sorted containers.
    equals(hyperellipsoid: Hyperellipsoid): boolean {
        return this.center.equals(hyperellipsoid.center)
            && axisEquals(this.axis, hyperellipsoid.axis)
            && this.extent.equals(hyperellipsoid.extent);
    }

    notEquals(hyperellipsoid: Hyperellipsoid): boolean {
        return !this.equals(hyperellipsoid);
    }

    lessThan(hyperellipsoid: Hyperellipsoid): boolean {
        if (this.center.lessThan(hyperellipsoid.center)) {
            return true;
        }

        if (this.center.greaterThan(hyperellipsoid.center)) {
            return false;
        }

        if (axisLess(this.axis, hyperellipsoid.axis)) {
            return true;
        }

        if (axisLess(hyperellipsoid.axis, this.axis)) {
            return false;
        }

        return this.extent.lessThan(hyperellipsoid.extent);
    }

    lessThanOrEqual(hyperellipsoid: Hyperellipsoid): boolean {
        return !hyperellipsoid.lessThan(this);
    }

    greaterThan(hyperellipsoid: Hyperellipsoid): boolean {
        return hyperellipsoid.lessThan(this);
    }

    greaterThanOrEqual(hyperellipsoid: Hyperellipsoid): boolean {
        return !this.lessThan(hyperellipsoid);
    }
}

// The port of the private static 'Convert(coeff, A, B, C)'.
function convertCoeffToABC(n: number, coeff: readonly number[]):
    { A: Matrix, B: Vector, C: number } {
    const A = Matrix.zero(n, n);
    const B = new Vector(n);
    let i = 0;
    const C = coeff[i++];

    for (let j = 0; j < n; ++j, ++i) {
        B.values[j] = coeff[i];
    }

    i = n + 1;
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < r; ++c) {
            A.set(r, c, A.get(c, r));
        }

        // When r = N-1, i = (N+1)*(N+2)/2 - 1 which corresponds to the last
        // element of coeff[]. The assignment is valid. After the assignment,
        // i is incremented and now out of range for coeff[]. However, the
        // loop after the assignment starts at c = N and the loop body is not
        // executed, after which the r-loop terminates.
        A.set(r, r, coeff[i]);
        ++i;
        for (let c = r + 1; c < n; ++c, ++i) {
            A.set(r, c, coeff[i] * 0.5);
        }
    }

    return { A, B, C };
}

// The port of the private static 'Convert(A, B, C, coeff)'.
function convertABCToCoeff(A: Matrix, B: Vector, C: number): number[] {
    const n = B.size;
    const numCoefficients = hyperellipsoidNumCoefficients(n);
    const coeff = new Array<number>(numCoefficients).fill(0);

    let i = 0;
    coeff[i++] = C;

    for (let j = 0; j < n; ++j, ++i) {
        coeff[i] = B.values[j];
    }

    i = n + 1;
    for (let r = 0; r < n; ++r) {
        coeff[i] = A.get(r, r);
        ++i;
        for (let c = r + 1; c < n && i < numCoefficients; ++c, ++i) {
            coeff[i] = A.get(r, c) * 2;
        }
    }

    return coeff;
}

// Aliases for convenience (the ports of the upstream template aliases).
export type Ellipse2 = Hyperellipsoid;
export type Ellipsoid3 = Hyperellipsoid;
