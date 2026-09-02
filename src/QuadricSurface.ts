// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) QuadricSurface.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A quadric surface is defined implicitly by
//
//   0 = q0 + q1 * x[0] + q2 * x[1] + q3 * x[2] +
//       q4 * x[0]^2 + q5 * x[0] * x[1] + q6 * x[0] * x[2] +
//       q7 * x[1]^2 + q8 * x[1] * x[2] + q9 * x[2]^2
//
//                                   +-              -+
//                                   | q4   q5/2 q6/2 |
//     = q0 + [q1 q2 q3] * X + X^T * | q5/2 q7   q8/2 | * X
//                                   | q6/2 q8/2 q9   |
//                                   +-              -+
//
//     = C + B^T*X + X^T*A*X
//
// A document describing the classification of the solution set is
// https://www.geometrictools.com/Documentation/ClassifyingQuadrics.pdf
//
// Port notes:
// - The two nondefault C++ constructors become the static factories
//   fromMatrix(A, b, c) and fromCoefficients(q), per the PORTING.md
//   precedent for ambiguous constructor overloads. The default constructor
//   zero-fills, as upstream does.
// - The member functions F, FX, ..., FZZ become camelCase f, fx, ..., fzz.
// - The nested enum class Classification becomes the file-qualified exported
//   enum QuadricSurfaceClassification (the bare name Classification is
//   already used privately by BSPPolygon2 and would not be unique under the
//   library-wide flat export).
// - Upstream's classification arithmetic uses BSRational<UIntegerAP32> for
//   exactness. The port uses the ported BSRational the same way. The port's
//   Matrix3x3/Matrix2x2 are aliases of the number-valued Matrix, so they
//   cannot carry rational entries; the small amount of rational linear
//   algebra needed (3x3 and 2x2 inverse, matrix-vector product, dot, cross)
//   is implemented by module-private helpers on flat BSRational arrays.
// - The rational inverses are computed as adjugate/determinant, which is the
//   value the upstream Inverse() returns for an invertible matrix. Both call
//   sites are reached only when the matrix is invertible.

import { logAssert } from './Logger';
import { BSRational } from './BSRational';
import { Matrix, mulMatrix } from './Matrix';
import { Vector, add, dot } from './Vector';

// The classification of the solution set of the quadric equation.
export enum QuadricSurfaceClassification {
    NO_SOLUTION,
    POINT,
    LINE,
    PLANE,
    TWO_PLANES,
    PARABOLIC_CYLINDER,
    ELLIPTIC_CYLINDER,
    HYPERBOLIC_CYLINDER,
    ELLIPTIC_PARABOLOID,
    HYPERBOLIC_PARABOLOID,
    ELLIPTIC_CONE,
    HYPERBOLOID_ONE_SHEET,
    HYPERBOLOID_TWO_SHEETS,
    ELLIPSOID,
    ENTIRE_SPACE,
    UNKNOWN
}

export class QuadricSurface {
    // Floating-point quadratic coefficients. The classification is cached
    // (upstream's 'mutable' member) and computed on demand.
    private mClassification: QuadricSurfaceClassification;
    private mA: Matrix;
    private mB: Vector;
    private mC: number;

    // The default constructor creates the zero quadric, whose solution set
    // is all of space.
    constructor() {
        this.mClassification = QuadricSurfaceClassification.UNKNOWN;
        this.mA = new Matrix(3, 3);
        this.mB = new Vector(3);
        this.mC = 0;
    }

    // 0 = c + b^T * X + X^T * A * X, where A is symmetric.
    static fromMatrix(A: Matrix, b: Vector, c: number): QuadricSurface {
        logAssert(A.numRows === 3 && A.numCols === 3 && b.size === 3,
            'Invalid size.');
        const surface = new QuadricSurface();
        surface.mA = A.clone();
        surface.mB = b.clone();
        surface.mC = c;
        return surface;
    }

    // The 10 coefficients q[0] through q[9] of the implicit equation.
    static fromCoefficients(q: readonly number[]): QuadricSurface {
        logAssert(q.length === 10, 'Invalid size.');
        const surface = new QuadricSurface();
        surface.mA = Matrix.fromArray(3, 3, [
            q[4], 0.5 * q[5], 0.5 * q[6],
            0.5 * q[5], q[7], 0.5 * q[8],
            0.5 * q[6], 0.5 * q[8], q[9]
        ]);
        surface.mB = Vector.fromArray([q[1], q[2], q[3]]);
        surface.mC = q[0];
        return surface;
    }

    // Member access.
    getA(): Matrix {
        return this.mA;
    }

    getB(): Vector {
        return this.mB;
    }

    getC(): number {
        return this.mC;
    }

    getQ(): number[] {
        const q = new Array<number>(10).fill(0);
        q[0] = this.mC;
        q[1] = this.mB.values[0];
        q[2] = this.mB.values[1];
        q[3] = this.mB.values[2];
        q[4] = this.mA.get(0, 0);
        q[5] = 2 * this.mA.get(0, 1);
        q[6] = 2 * this.mA.get(0, 2);
        q[7] = this.mA.get(1, 1);
        q[8] = 2 * this.mA.get(1, 2);
        q[9] = this.mA.get(2, 2);
        return q;
    }

    // Evaluate the function.
    f(position: Vector): number {
        return dot(position, add(mulMatrix(this.mA, position), this.mB)) + this.mC;
    }

    // Evaluate the first-order partial derivatives (gradient).
    fx(position: Vector): number {
        const sum = this.mA.get(0, 0) * position.values[0] +
            this.mA.get(0, 1) * position.values[1] +
            this.mA.get(0, 2) * position.values[2];
        return 2 * sum + this.mB.values[0];
    }

    fy(position: Vector): number {
        const sum = this.mA.get(1, 0) * position.values[0] +
            this.mA.get(1, 1) * position.values[1] +
            this.mA.get(1, 2) * position.values[2];
        return 2 * sum + this.mB.values[1];
    }

    fz(position: Vector): number {
        const sum = this.mA.get(2, 0) * position.values[0] +
            this.mA.get(2, 1) * position.values[1] +
            this.mA.get(2, 2) * position.values[2];
        return 2 * sum + this.mB.values[2];
    }

    // Evaluate the second-order partial derivatives (Hessian). The position
    // is unused; it is kept for signature symmetry with upstream.
    fxx(): number {
        return 2 * this.mA.get(0, 0);
    }

    fxy(): number {
        return 2 * this.mA.get(0, 1);
    }

    fxz(): number {
        return 2 * this.mA.get(0, 2);
    }

    fyy(): number {
        return 2 * this.mA.get(1, 1);
    }

    fyz(): number {
        return 2 * this.mA.get(1, 2);
    }

    fzz(): number {
        return 2 * this.mA.get(2, 2);
    }

    // Classification of the quadric. The implementation uses exact rational
    // arithmetic to avoid misclassification due to floating-point rounding
    // errors.
    getClassification(): QuadricSurfaceClassification {
        if (this.mClassification !== QuadricSurfaceClassification.UNKNOWN) {
            return this.mClassification;
        }

        // Convert the coefficients to their rational representations and
        // compute various derived quantities. The upper triangle of the
        // input is used, so an asymmetric input matrix is symmetrized the
        // way upstream symmetrizes it.
        const rA = rZeroMatrix(3);
        rSet(rA, 3, 0, 0, num(this.mA.get(0, 0)));
        rSet(rA, 3, 0, 1, num(this.mA.get(0, 1)));
        rSet(rA, 3, 0, 2, num(this.mA.get(0, 2)));
        rSet(rA, 3, 1, 0, rGet(rA, 3, 0, 1));
        rSet(rA, 3, 1, 1, num(this.mA.get(1, 1)));
        rSet(rA, 3, 1, 2, num(this.mA.get(1, 2)));
        rSet(rA, 3, 2, 0, rGet(rA, 3, 0, 2));
        rSet(rA, 3, 2, 1, rGet(rA, 3, 1, 2));
        rSet(rA, 3, 2, 2, num(this.mA.get(2, 2)));
        const rB: BSRational[] = [
            num(this.mB.values[0]),
            num(this.mB.values[1]),
            num(this.mB.values[2])
        ];
        const rC = num(this.mC);

        // Compute the polynomial det(lambda * I - A) with rational
        // coefficients.
        const a = (r: number, c: number) => rGet(rA, 3, r, c);
        const rS00 = a(1, 1).mul(a(2, 2)).sub(a(1, 2).mul(a(1, 2)));
        const rS01 = a(0, 1).mul(a(2, 2)).sub(a(1, 2).mul(a(0, 2)));
        const rS02 = a(0, 1).mul(a(1, 2)).sub(a(0, 2).mul(a(1, 1)));
        const rS11 = a(0, 0).mul(a(2, 2)).sub(a(0, 2).mul(a(0, 2)));
        // Upstream also computes rS12 = A00*A12 - A02*A01 here, but never
        // uses it; the port drops the dead computation.
        const rS22 = a(0, 0).mul(a(1, 1)).sub(a(0, 1).mul(a(0, 1)));
        const rP: BSRational[] = [
            a(0, 0).mul(rS00).sub(a(0, 1).mul(rS01)).add(a(0, 2).mul(rS02)).negated(),
            rS00.add(rS11).add(rS22),
            a(0, 0).add(a(1, 1)).add(a(2, 2)).negated(),
            num(1)
        ];

        // Determine the signs of the roots.
        const { numPositive, numNegative, numZero } = computeRootSigns(rP);

        // Classify the solution set to the equation.
        if (numZero === 0) {
            this.mClassification = allNonzero(rA, rB, rC, numPositive);
        } else if (numZero === 1) {
            this.mClassification = twoNonzero(rA, rB, rC, numPositive, numNegative);
        } else if (numZero === 2) {
            this.mClassification = oneNonzero(rA, rB, rC, numPositive);
        } else {  // numZero = 3
            this.mClassification = allZero(rB, rC);
        }
        return this.mClassification;
    }
}

// Rational helpers. A rational NxN matrix is a flat row-major BSRational[].

function num(value: number): BSRational {
    return BSRational.fromNumber(value);
}

function rZeroMatrix(n: number): BSRational[] {
    const m = new Array<BSRational>(n * n);
    for (let i = 0; i < n * n; ++i) {
        m[i] = new BSRational();
    }
    return m;
}

function rGet(m: readonly BSRational[], n: number, r: number, c: number): BSRational {
    return m[c + n * r];
}

function rSet(m: BSRational[], n: number, r: number, c: number,
    value: BSRational): void {
    m[c + n * r] = value;
}

function rIsZeroVector(v: readonly BSRational[]): boolean {
    for (const component of v) {
        if (component.getSign() !== 0) {
            return false;
        }
    }
    return true;
}

function rDot(u: readonly BSRational[], v: readonly BSRational[]): BSRational {
    let sum = u[0].mul(v[0]);
    for (let i = 1; i < u.length; ++i) {
        sum = sum.add(u[i].mul(v[i]));
    }
    return sum;
}

function rCross(u: readonly BSRational[], v: readonly BSRational[]): BSRational[] {
    return [
        u[1].mul(v[2]).sub(u[2].mul(v[1])),
        u[2].mul(v[0]).sub(u[0].mul(v[2])),
        u[0].mul(v[1]).sub(u[1].mul(v[0]))
    ];
}

function rMulMatVec(m: readonly BSRational[], n: number,
    v: readonly BSRational[]): BSRational[] {
    const result = new Array<BSRational>(n);
    for (let r = 0; r < n; ++r) {
        let sum = rGet(m, n, r, 0).mul(v[0]);
        for (let c = 1; c < n; ++c) {
            sum = sum.add(rGet(m, n, r, c).mul(v[c]));
        }
        result[r] = sum;
    }
    return result;
}

function rGetRow(m: readonly BSRational[], n: number, r: number): BSRational[] {
    const row = new Array<BSRational>(n);
    for (let c = 0; c < n; ++c) {
        row[c] = rGet(m, n, r, c);
    }
    return row;
}

// The inverse of an invertible rational 3x3 matrix, computed as the adjugate
// divided by the determinant. Both call sites guarantee invertibility.
function rInverse3x3(m: readonly BSRational[]): BSRational[] {
    const g = (r: number, c: number) => rGet(m, 3, r, c);
    const c00 = g(1, 1).mul(g(2, 2)).sub(g(1, 2).mul(g(2, 1)));
    const c01 = g(1, 2).mul(g(2, 0)).sub(g(1, 0).mul(g(2, 2)));
    const c02 = g(1, 0).mul(g(2, 1)).sub(g(1, 1).mul(g(2, 0)));
    const det = g(0, 0).mul(c00).add(g(0, 1).mul(c01)).add(g(0, 2).mul(c02));
    logAssert(det.getSign() !== 0, 'The matrix must be invertible.');

    const adj: BSRational[] = [
        c00,
        g(0, 2).mul(g(2, 1)).sub(g(0, 1).mul(g(2, 2))),
        g(0, 1).mul(g(1, 2)).sub(g(0, 2).mul(g(1, 1))),
        c01,
        g(0, 0).mul(g(2, 2)).sub(g(0, 2).mul(g(2, 0))),
        g(0, 2).mul(g(1, 0)).sub(g(0, 0).mul(g(1, 2))),
        c02,
        g(0, 1).mul(g(2, 0)).sub(g(0, 0).mul(g(2, 1))),
        g(0, 0).mul(g(1, 1)).sub(g(0, 1).mul(g(1, 0)))
    ];
    return adj.map(value => value.div(det));
}

// The inverse of an invertible rational 2x2 matrix.
function rInverse2x2(m: readonly BSRational[]): BSRational[] {
    const g = (r: number, c: number) => rGet(m, 2, r, c);
    const det = g(0, 0).mul(g(1, 1)).sub(g(0, 1).mul(g(1, 0)));
    logAssert(det.getSign() !== 0, 'The matrix must be invertible.');
    return [
        g(1, 1).div(det), g(0, 1).negated().div(det),
        g(1, 0).negated().div(det), g(0, 0).div(det)
    ];
}

// Use Descartes' rule of signs to determine the root signs. Because A is
// symmetric, all three roots of the characteristic polynomial are real, so
// the rule reports the exact counts rather than only upper bounds.
function computeRootSigns(rP: readonly BSRational[]):
    { numPositive: number, numNegative: number, numZero: number } {
    // Collect the nonzero signs of rP[0], rP[1] and rP[2].
    const degree = 3;
    const signs = new Array<number>(degree + 1).fill(0);
    for (let i = 0; i <= degree; ++i) {
        signs[i] = rP[i].getSign();
    }

    // Compute the number of positive roots of p(lambda).
    let currentSign = signs[degree];
    let numPositive = 0;
    for (let i = degree - 1; i >= 0; --i) {
        if (signs[i] === -currentSign) {
            currentSign = signs[i];
            ++numPositive;
        }
    }

    // Compute the signs of the coefficients of p(-lambda).
    for (let i = 1; i <= degree; i += 2) {
        signs[i] = -signs[i];
    }

    // Compute the number of positive roots of p(-lambda).
    currentSign = signs[degree];
    let numNegative = 0;
    for (let i = degree - 1; i >= 0; --i) {
        if (signs[i] === -currentSign) {
            currentSign = signs[i];
            ++numNegative;
        }
    }

    // Compute the number of zero roots of p(lambda).
    const numZero = 3 - numPositive - numNegative;
    return { numPositive, numNegative, numZero };
}

function allNonzero(A: readonly BSRational[], b: readonly BSRational[],
    c: BSRational, numPositiveRoots: number): QuadricSurfaceClassification {
    const r = rDot(b, rMulMatVec(rInverse3x3(A), 3, b)).div(num(4)).sub(c);
    const sign = r.getSign();

    if (sign > 0) {
        if (numPositiveRoots === 3) {
            return QuadricSurfaceClassification.ELLIPSOID;
        } else if (numPositiveRoots === 2) {
            return QuadricSurfaceClassification.HYPERBOLOID_ONE_SHEET;
        } else if (numPositiveRoots === 1) {
            return QuadricSurfaceClassification.HYPERBOLOID_TWO_SHEETS;
        } else {
            return QuadricSurfaceClassification.NO_SOLUTION;
        }
    } else if (sign < 0) {
        if (numPositiveRoots === 3) {
            return QuadricSurfaceClassification.NO_SOLUTION;
        } else if (numPositiveRoots === 2) {
            return QuadricSurfaceClassification.HYPERBOLOID_TWO_SHEETS;
        } else if (numPositiveRoots === 1) {
            return QuadricSurfaceClassification.HYPERBOLOID_ONE_SHEET;
        } else {
            return QuadricSurfaceClassification.ELLIPSOID;
        }
    } else {  // sign == 0
        if (numPositiveRoots === 3 || numPositiveRoots === 0) {
            return QuadricSurfaceClassification.POINT;
        } else {
            return QuadricSurfaceClassification.ELLIPTIC_CONE;
        }
    }
}

// A has rank 2. The vector w0 spans the null space of A and {w1,w2} is an
// orthogonal basis for the row space of A.
function computeOrthogonalSetTwoNonzero(A: readonly BSRational[]):
    { w0: BSRational[], w1: BSRational[], w2: BSRational[] } {
    let w0: BSRational[], w1: BSRational[], w2: BSRational[];
    w1 = rGetRow(A, 3, 0);
    if (!rIsZeroVector(w1)) {
        w2 = rGetRow(A, 3, 1);
        w0 = rCross(w1, w2);
        if (rIsZeroVector(w0)) {
            w2 = rGetRow(A, 3, 2);
            w0 = rCross(w1, w2);
        }
    } else {
        w1 = rGetRow(A, 3, 1);
        w2 = rGetRow(A, 3, 2);
        w0 = rCross(w1, w2);
    }
    w2 = rCross(w0, w1);
    return { w0, w1, w2 };
}

function twoNonzero(A: readonly BSRational[], b: readonly BSRational[],
    c: BSRational, numPositive: number,
    numNegative: number): QuadricSurfaceClassification {
    const { w0, w1, w2 } = computeOrthogonalSetTwoNonzero(A);
    const d0 = rDot(w0, b);
    if (d0.getSign() !== 0) {
        if (numPositive === numNegative) {
            return QuadricSurfaceClassification.HYPERBOLIC_PARABOLOID;
        } else {
            return QuadricSurfaceClassification.ELLIPTIC_PARABOLOID;
        }
    }

    const Aw1 = rMulMatVec(A, 3, w1);
    const Aw2 = rMulMatVec(A, 3, w2);
    const E01 = rDot(w1, Aw2);
    const E: BSRational[] = [
        rDot(w1, Aw1), E01,
        E01, rDot(w2, Aw2)
    ];
    const f: BSRational[] = [rDot(w1, b), rDot(w2, b)];
    const r = rDot(f, rMulMatVec(rInverse2x2(E), 2, f)).div(num(4)).sub(c);
    const sign = r.getSign();

    if (numPositive === 2) {
        if (sign > 0) {
            return QuadricSurfaceClassification.ELLIPTIC_CYLINDER;
        } else if (sign < 0) {
            return QuadricSurfaceClassification.NO_SOLUTION;
        } else {
            return QuadricSurfaceClassification.LINE;
        }
    } else if (numNegative === 2) {
        if (sign < 0) {
            return QuadricSurfaceClassification.ELLIPTIC_CYLINDER;
        } else if (sign > 0) {
            return QuadricSurfaceClassification.NO_SOLUTION;
        } else {
            return QuadricSurfaceClassification.LINE;
        }
    } else {  // numPositive = numNegative = 1
        if (sign !== 0) {
            return QuadricSurfaceClassification.HYPERBOLIC_CYLINDER;
        } else {
            return QuadricSurfaceClassification.TWO_PLANES;
        }
    }
}

// A has rank 1. The vector w2 spans the row space of A and {w0,w1} is an
// orthogonal basis for the null space of A.
function computeOrthogonalSetOneNonzero(A: readonly BSRational[]):
    { w0: BSRational[], w1: BSRational[], w2: BSRational[] } {
    let w2 = rGetRow(A, 3, 0);
    if (rIsZeroVector(w2)) {
        w2 = rGetRow(A, 3, 1);
        if (rIsZeroVector(w2)) {
            w2 = rGetRow(A, 3, 2);
        }
    }

    let w0: BSRational[];
    if (BSRational.fabs(w2[0]).greaterThan(BSRational.fabs(w2[1]))) {
        w0 = [w2[2].negated(), new BSRational(), w2[0]];
    } else {
        w0 = [new BSRational(), w2[2], w2[1].negated()];
    }
    const w1 = rCross(w2, w0);
    return { w0, w1, w2 };
}

function oneNonzero(A: readonly BSRational[], b: readonly BSRational[],
    c: BSRational, numPositive: number): QuadricSurfaceClassification {
    const { w0, w1, w2 } = computeOrthogonalSetOneNonzero(A);
    const d0 = rDot(w0, b);
    const d1 = rDot(w1, b);
    if (d0.getSign() !== 0 || d1.getSign() !== 0) {
        return QuadricSurfaceClassification.PARABOLIC_CYLINDER;
    }

    const E = rDot(w2, rMulMatVec(A, 3, w2));
    const f = rDot(w2, b);
    const r = f.mul(f).div(num(4).mul(E)).sub(c);
    const sign = r.getSign();

    if (numPositive === 1) {  // numNegative = 0
        if (sign > 0) {
            return QuadricSurfaceClassification.TWO_PLANES;
        } else if (sign < 0) {
            return QuadricSurfaceClassification.NO_SOLUTION;
        } else {
            return QuadricSurfaceClassification.PLANE;
        }
    } else {  // numPositive = 0, numNegative = 1
        if (sign < 0) {
            return QuadricSurfaceClassification.TWO_PLANES;
        } else if (sign > 0) {
            return QuadricSurfaceClassification.NO_SOLUTION;
        } else {
            return QuadricSurfaceClassification.PLANE;
        }
    }
}

function allZero(b: readonly BSRational[],
    c: BSRational): QuadricSurfaceClassification {
    if (!rIsZeroVector(b)) {
        return QuadricSurfaceClassification.PLANE;
    } else {
        if (c.getSign() === 0) {
            return QuadricSurfaceClassification.ENTIRE_SPACE;
        } else {
            return QuadricSurfaceClassification.NO_SOLUTION;
        }
    }
}
