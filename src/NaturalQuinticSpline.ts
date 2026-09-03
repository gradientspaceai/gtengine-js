// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) NaturalQuinticSpline.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Documentation for natural splines is found in
// https://www.geometrictools.com/Documentation/NaturalSplines.pdf
// The number of points must be 3 or larger. The points[] and times[] arrays
// must have the same number of elements. The times[] values must be strictly
// increasing.
//
// The spline interpolates both the values f0[] and the first derivatives
// f1[] at the knots. Each segment is a quintic with 6 coefficients, 4 of
// which are consumed by the Hermite interpolation conditions, so the spline
// is C3 at the interior knots (the fourth derivative jumps in general). A
// free spline has zero third derivative at the endpoints, a clamped spline
// has the specified second derivatives there, and a closed spline matches
// the second and third derivatives at the two endpoints.
//
// Port notes (mirroring NaturalCubicSpline):
// - Upstream 'template <int32_t N, typename T>' becomes a runtime dimension,
//   taken from the interpolation points f0[0]; all points and derivatives
//   are asserted to have that dimension.
// - The four upstream constructors become three static factories:
//   createFree, createClosed (upstream's isFree = true/false) and
//   createClamped.
// - evaluate supports order <= 5, but ParametricCurve.SUP_ORDER is 4 (the
//   maximum for the curve family). createJet() is therefore overridden to
//   allocate SUP_ORDER_QUINTIC = 6 entries so callers can request the
//   fourth- and fifth-order derivatives.

import { logAssert } from './Logger.js';
import { Matrix } from './Matrix.js';
import { inverse4x4 } from './Matrix4x4.js';
import { ParametricCurve } from './ParametricCurve.js';
import { Vector, add, mul, sub } from './Vector.js';

// A quintic polynomial segment: the six Vector coefficients of
// p(u) = c[0] + u*(c[1] + u*(c[2] + u*(c[3] + u*(c[4] + u*c[5])))) for u in
// [0,1].
export type NaturalQuinticSplinePolynomial = Vector[];

type NaturalQuinticSplineKind = 'free' | 'closed' | 'clamped';

export class NaturalQuinticSpline extends ParametricCurve {
    // The number of entries a 'jet' array must have to hold the position and
    // the derivatives of orders 1 through 5.
    static readonly SUP_ORDER_QUINTIC = 6;

    private mPolynomials: NaturalQuinticSplinePolynomial[];
    private mDelta: number[];

    // Construct a free spline (upstream's isFree = true). The function values
    // are f0[] and the first derivative values are f1[].
    static createFree(f0: readonly Vector[], f1: readonly Vector[],
        times: readonly number[]): NaturalQuinticSpline {
        return new NaturalQuinticSpline('free', f0, f1, times, null, null);
    }

    // Construct a closed spline (upstream's isFree = false).
    static createClosed(f0: readonly Vector[], f1: readonly Vector[],
        times: readonly number[]): NaturalQuinticSpline {
        return new NaturalQuinticSpline('closed', f0, f1, times, null, null);
    }

    // Construct a clamped spline: the second derivatives at the endpoints
    // are the specified values.
    static createClamped(f0: readonly Vector[], f1: readonly Vector[],
        times: readonly number[], derivative0: Vector,
        derivative1: Vector): NaturalQuinticSpline {
        return new NaturalQuinticSpline('clamped', f0, f1, times, derivative0,
            derivative1);
    }

    private constructor(kind: NaturalQuinticSplineKind, f0: readonly Vector[],
        f1: readonly Vector[], times: readonly number[],
        derivative0: Vector | null, derivative1: Vector | null) {
        super(f0.length > 0 ? f0[0].size : 0, f0.length - 1, times);

        const numPoints = f0.length;
        logAssert(numPoints >= 3 && f1.length === numPoints
            && times.length === numPoints, 'Invalid input.');
        for (let i = 0; i < numPoints; ++i) {
            logAssert(f0[i].size === this.mDimension
                && f1[i].size === this.mDimension,
                'Mismatched point dimensions.');
        }

        const numPm1 = numPoints - 1;
        this.mPolynomials = new Array<NaturalQuinticSplinePolynomial>(numPm1);
        for (let i = 0; i < numPm1; ++i) {
            const poly = new Array<Vector>(6);
            for (let j = 0; j < 6; ++j) {
                poly[j] = new Vector(this.mDimension);
            }
            this.mPolynomials[i] = poly;
        }
        this.mDelta = new Array<number>(numPm1).fill(0);
        for (let i0 = 0, i1 = 1; i1 < numPoints; i0 = i1++) {
            this.mDelta[i0] = times[i1] - times[i0];
        }

        const numPm2 = numPoints - 2;
        let boundary0: Vector;
        let boundary1: Vector;
        if (kind === 'clamped') {
            const coeff0 = 0.5 * this.mDelta[0] * this.mDelta[0];
            const coeff1 = 0.5 * this.mDelta[numPm2] * this.mDelta[numPm2];
            boundary0 = mul(coeff0, derivative0 as Vector);
            boundary1 = mul(coeff1, derivative1 as Vector);
        }
        else {
            // Free splines and closed splines have the last two B-entries
            // set to the zero vector.
            boundary0 = new Vector(this.mDimension);
            boundary1 = new Vector(this.mDimension);
        }

        const R = new Matrix(4, 4);
        const numBElements = 4 * numPm1;
        const B = new Array<Vector>(numBElements);
        for (let i = 0; i < numBElements; ++i) {
            B[i] = new Vector(this.mDimension);
        }
        this.onPresolve(numPoints, f0, f1, boundary0, boundary1, R, B);

        if (kind === 'free') {
            R.set(2, 1, 1);
            R.set(2, 2, 4);
            R.set(2, 3, 10);
            this.solve(0, 1, numPoints, f0, f1, R, B);
        }
        else if (kind === 'closed') {
            const lambda = this.mDelta[0] / this.mDelta[numPm2];
            const lambdasqr = lambda * lambda;
            const lambdacub = lambdasqr * lambda;
            R.set(2, 0, -lambdasqr);
            R.set(2, 1, -3 * lambdasqr);
            R.set(2, 2, -6 * lambdasqr);
            R.set(2, 3, -10 * lambdasqr);
            R.set(3, 1, -lambdacub);
            R.set(3, 2, -4 * lambdacub);
            R.set(3, 3, -10 * lambdacub);
            this.solve(1, 1, numPoints, f0, f1, R, B);
        }
        else {
            R.set(3, 0, 1);
            R.set(3, 1, 3);
            R.set(3, 2, 6);
            R.set(3, 3, 10);
            this.solve(1, 0, numPoints, f0, f1, R, B);
        }

        this.mConstructed = true;
    }

    // The live array of polynomial segments; do not modify it.
    getPolynomials(): readonly NaturalQuinticSplinePolynomial[] {
        return this.mPolynomials;
    }

    // Allocate a jet with SUP_ORDER_QUINTIC entries so that evaluate can be
    // called with order up to 5 (the base class allocates only SUP_ORDER = 4
    // entries, enough for order 3).
    override createJet(): Vector[] {
        const n = NaturalQuinticSpline.SUP_ORDER_QUINTIC;
        const jet = new Array<Vector>(n);
        for (let i = 0; i < n; ++i) {
            jet[i] = new Vector(this.mDimension);
        }
        return jet;
    }

    // Evaluation of the function and its derivatives through order 5. If you
    // want only the position, pass in order 0. If you want the position and
    // first derivative, pass in order of 1 and so on. The output array 'jet'
    // must have 'order + 1' elements (use createJet()). The values are
    // ordered as position, first derivative, second derivative and so on.
    override evaluate(t: number, order: number, jet: Vector[]): void {
        if (!this.mConstructed) {
            // Return a zero-valued jet for invalid state.
            for (let i = 0; i <= order; ++i) {
                jet[i] = new Vector(this.mDimension);
            }
            return;
        }

        const info = this.getKeyInfo(t);
        const key = info.key;
        const u = info.u;
        const poly = this.mPolynomials[key];

        // Compute position.
        jet[0] = add(poly[0], mul(u, add(poly[1], mul(u, add(poly[2],
            mul(u, add(poly[3], mul(u, add(poly[4], mul(u, poly[5]))))))))));
        if (order >= 1) {
            // Compute first derivative.
            let denom = this.mDelta[key];
            jet[1] = mul(add(poly[1], mul(u, add(mul(2, poly[2]),
                mul(u, add(mul(3, poly[3]), mul(u, add(mul(4, poly[4]),
                    mul(u, mul(5, poly[5]))))))))), 1 / denom);
            if (order >= 2) {
                // Compute second derivative.
                denom *= this.mDelta[key];
                jet[2] = mul(add(mul(2, poly[2]), mul(u, add(mul(6, poly[3]),
                    mul(u, add(mul(12, poly[4]), mul(u, mul(20, poly[5]))))))),
                    1 / denom);
                if (order >= 3) {
                    // Compute third derivative.
                    denom *= this.mDelta[key];
                    jet[3] = mul(add(mul(6, poly[3]), mul(u,
                        add(mul(24, poly[4]), mul(u, mul(60, poly[5]))))),
                        1 / denom);
                    if (order >= 4) {
                        // Compute fourth derivative.
                        denom *= this.mDelta[key];
                        jet[4] = mul(add(mul(24, poly[4]),
                            mul(u, mul(120, poly[5]))), 1 / denom);

                        if (order >= 5) {
                            // Compute fifth derivative.
                            denom *= this.mDelta[key];
                            jet[5] = mul(mul(120, poly[5]), 1 / denom);

                            for (let i = 6; i <= order; ++i) {
                                // Derivatives of order 6 and higher are zero.
                                jet[i] = new Vector(this.mDimension);
                            }
                        }
                    }
                }
            }
        }
    }

    private onPresolve(numPoints: number, f0: readonly Vector[],
        f1: readonly Vector[], boundary0: Vector, boundary1: Vector,
        R: Matrix, B: Vector[]): void {
        const numPm1 = numPoints - 1;
        const numPm2 = numPoints - 2;
        const numPm3 = numPoints - 3;

        const coeff0 = [10, -20, 15, -4];
        const coeff1 = [-6, 14, -11, 3];
        for (let i0 = 0, i1 = 1; i0 <= numPm3; i0 = i1++) {
            const diff0 = sub(sub(f0[i1], f0[i0]), mul(this.mDelta[i0], f1[i0]));
            const diff1 = mul(this.mDelta[i0], sub(f1[i1], f1[i0]));
            for (let j = 0, k = 4 * i0; j < 4; ++j, ++k) {
                B[k] = add(mul(coeff0[j], diff0), mul(coeff1[j], diff1));
            }
        }

        B[B.length - 4] = sub(sub(f0[numPm1], f0[numPm2]),
            mul(this.mDelta[numPm2], f1[numPm2]));
        B[B.length - 3] = mul(this.mDelta[numPm2], sub(f1[numPm1], f1[numPm2]));
        B[B.length - 2] = boundary0;
        B[B.length - 1] = boundary1;

        R.set(0, 0, 1);
        R.set(0, 1, 1);
        R.set(0, 2, 1);
        R.set(0, 3, 1);
        R.set(1, 0, 2);
        R.set(1, 1, 3);
        R.set(1, 2, 4);
        R.set(1, 3, 5);
    }

    private solve(ell20: number, ell31: number, numPoints: number,
        f0: readonly Vector[], f1: readonly Vector[], R: Matrix,
        B: Vector[]): void {
        this.rowReduce(ell20, ell31, numPoints, R, B);
        this.backSubstitute(f0, f1, R, B);
    }

    // Apply the row reductions that convert the matrix system to an
    // upper-triangular block-matrix system.
    private rowReduce(ell20: number, ell31: number, numPoints: number,
        R: Matrix, B: Vector[]): void {
        const numPm3 = numPoints - 3;

        if (ell20 === 1) {
            const trg = B.length - 2;
            B[trg] = sub(B[trg], B[0]);
            let sigma = this.mDelta[0] / this.mDelta[1];
            let sigmasqr = sigma * sigma;
            let sigmacub = sigmasqr * sigma;
            let luProd0 = -3 * sigmasqr;
            let luProd1 = sigmacub;
            let sign = -1;

            for (let i = 1; i <= numPm3; ++i) {
                B[trg] = sub(B[trg], mul(sign,
                    add(mul(luProd0, B[4 * i]), mul(luProd1, B[4 * i + 1]))));
                sigma = this.mDelta[i] / this.mDelta[i + 1];
                sigmasqr = sigma * sigma;
                sigmacub = sigmasqr * sigma;
                const temp0 = sigmasqr * (-3 * luProd0 + 8 * luProd1);
                const temp1 = sigmacub * (luProd0 - 3 * luProd1);
                luProd0 = temp0;
                luProd1 = temp1;
                sign = -sign;
            }

            R.set(2, 0, R.get(2, 0) + sign * luProd0);
            R.set(2, 1, R.get(2, 1) + sign * luProd1);
        }

        if (ell31 === 1) {
            const trg = B.length - 1;
            B[trg] = sub(B[trg], B[1]);
            let sigma = this.mDelta[0] / this.mDelta[1];
            let sigmasqr = sigma * sigma;
            let sigmacub = sigmasqr * sigma;
            let luProd0 = 8 * sigmasqr;
            let luProd1 = -3 * sigmacub;
            let sign = -1;

            for (let i = 1; i <= numPm3; ++i) {
                B[trg] = sub(B[trg], mul(sign,
                    add(mul(luProd0, B[4 * i]), mul(luProd1, B[4 * i + 1]))));
                sigma = this.mDelta[i] / this.mDelta[i + 1];
                sigmasqr = sigma * sigma;
                sigmacub = sigmasqr * sigma;
                const temp0 = sigmasqr * (-3 * luProd0 + 8 * luProd1);
                const temp1 = sigmacub * (luProd0 - 3 * luProd1);
                luProd0 = temp0;
                luProd1 = temp1;
                sign = -sign;
            }

            R.set(3, 0, R.get(3, 0) + sign * luProd0);
            R.set(3, 1, R.get(3, 1) + sign * luProd1);
        }
    }

    private backSubstitute(f0: readonly Vector[], f1: readonly Vector[],
        R: Matrix, B: readonly Vector[]): void {
        const inv = inverse4x4(R);
        logAssert(inv.invertible, 'R matrix is not invertible.');
        const invR = inv.inverse;

        const numPolynomials = this.mPolynomials.length;
        const poly = this.mPolynomials[numPolynomials - 1];
        let j0 = B.length - 4;
        let j1 = j0 + 1;
        let j2 = j0 + 2;
        let j3 = j0 + 3;

        poly[0] = f0[numPolynomials - 1].clone();
        poly[1] = mul(f1[numPolynomials - 1], this.mDelta[numPolynomials - 1]);
        for (let r = 0; r < 4; ++r) {
            poly[r + 2] = add(add(mul(invR.get(r, 0), B[j0]),
                mul(invR.get(r, 1), B[j1])),
                add(mul(invR.get(r, 2), B[j2]), mul(invR.get(r, 3), B[j3])));
        }

        for (let i1 = numPolynomials - 2, i0 = i1 + 1; i1 >= 0; i0 = i1--) {
            const prev = this.mPolynomials[i0];
            const curr = this.mPolynomials[i1];
            const sigma = this.mDelta[i1] / this.mDelta[i0];
            const sigmasqr = sigma * sigma;
            const sigmacub = sigmasqr * sigma;
            const u00 = -3 * sigmasqr;
            const u01 = sigmacub;
            const u10 = 8 * sigmasqr;
            const u11 = -3 * sigmacub;
            const u20 = -7 * sigmasqr;
            const u21 = 3 * sigmacub;
            const u30 = 2 * sigmasqr;
            const u31 = -sigmacub;

            j0 -= 4;
            j1 -= 4;
            j2 -= 4;
            j3 -= 4;

            curr[0] = f0[i1].clone();
            curr[1] = mul(f1[i1], this.mDelta[i1]);
            curr[2] = sub(B[j0], add(mul(u00, prev[2]), mul(u01, prev[3])));
            curr[3] = sub(B[j1], add(mul(u10, prev[2]), mul(u11, prev[3])));
            curr[4] = sub(B[j2], add(mul(u20, prev[2]), mul(u21, prev[3])));
            curr[5] = sub(B[j3], add(mul(u30, prev[2]), mul(u31, prev[3])));
        }
    }

    // Determine the index key for which times[key] <= t < times[key+1].
    // Return u = (t - times[key]) / delta[key], which is in [0,1].
    private getKeyInfo(t: number): { key: number, u: number } {
        const numSegments = this.getNumSegments();
        if (t > this.mTime[0]) {
            if (t < this.mTime[numSegments]) {
                for (let i = 0, ip1 = 1; i < numSegments; i = ip1++) {
                    if (t < this.mTime[ip1]) {
                        return {
                            key: i,
                            u: (t - this.mTime[i]) / this.mDelta[i]
                        };
                    }
                }
                // Unreachable: t < mTime[numSegments] guarantees a match.
                return { key: numSegments - 1, u: 1 };
            }
            else {
                return { key: numSegments - 1, u: 1 };
            }
        }
        else {
            return { key: 0, u: 0 };
        }
    }
}
