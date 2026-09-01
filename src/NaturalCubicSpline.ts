// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) NaturalCubicSpline.h
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
// The spline interpolates the values f0[] at the knots and is C2 at the
// interior knots. A free spline has zero second derivative at the
// endpoints, a clamped spline has the specified first derivatives there,
// and a closed spline matches the first and second derivatives at the two
// endpoints.
//
// Port notes:
// - Upstream 'template <int32_t N, typename T>' becomes a runtime dimension.
//   It is taken from the interpolation points f0[0] rather than being passed
//   in, because f0 is required to be nonempty; all points are asserted to
//   have the same dimension.
// - Upstream has four constructors: (isFree, numPoints, f0, times) and its
//   std::vector overload build a free (isFree = true) or a closed spline,
//   and the remaining pair build a clamped spline. The (pointer, count)
//   forms have no TypeScript analogue and the remaining two overloads are
//   ambiguous, so the port exposes three static factories instead:
//   createFree, createClosed and createClamped (the "ambiguous C++
//   constructor overloads become static factories" precedent).
// - getPolynomials returns the live array of polynomial segments; each
//   segment is an array of 4 Vector coefficients, the port of upstream's
//   'std::array<Vector<N,T>, 4>'.
// - evaluate supports order <= 3, so ParametricCurve.createJet() (which
//   allocates SUP_ORDER = 4 entries) is the right allocation for callers.

import { logAssert } from './Logger';
import { Matrix } from './Matrix';
import { inverse3x3 } from './Matrix3x3';
import { ParametricCurve } from './ParametricCurve';
import { Vector, add, mul, sub } from './Vector';

// A cubic polynomial segment: the four Vector coefficients of
// p(u) = c[0] + u*(c[1] + u*(c[2] + u*c[3])) for u in [0,1].
export type NaturalCubicSplinePolynomial = Vector[];

type NaturalCubicSplineKind = 'free' | 'closed' | 'clamped';

export class NaturalCubicSpline extends ParametricCurve {
    private mPolynomials: NaturalCubicSplinePolynomial[];
    private mDelta: number[];

    // Construct a free spline (upstream's isFree = true).
    static createFree(f0: readonly Vector[],
        times: readonly number[]): NaturalCubicSpline {
        return new NaturalCubicSpline('free', f0, times, null, null);
    }

    // Construct a closed spline (upstream's isFree = false).
    static createClosed(f0: readonly Vector[],
        times: readonly number[]): NaturalCubicSpline {
        return new NaturalCubicSpline('closed', f0, times, null, null);
    }

    // Construct a clamped spline: the first derivatives at the endpoints are
    // the specified values.
    static createClamped(f0: readonly Vector[], times: readonly number[],
        derivative0: Vector, derivative1: Vector): NaturalCubicSpline {
        return new NaturalCubicSpline('clamped', f0, times, derivative0,
            derivative1);
    }

    private constructor(kind: NaturalCubicSplineKind, f0: readonly Vector[],
        times: readonly number[], derivative0: Vector | null,
        derivative1: Vector | null) {
        super(f0.length > 0 ? f0[0].size : 0, f0.length - 1, times);

        const numPoints = f0.length;
        logAssert(numPoints >= 3 && times.length === numPoints,
            'Invalid input.');
        for (let i = 0; i < numPoints; ++i) {
            logAssert(f0[i].size === this.mDimension,
                'Mismatched point dimensions.');
        }

        const numPm1 = numPoints - 1;
        this.mPolynomials = new Array<NaturalCubicSplinePolynomial>(numPm1);
        for (let i = 0; i < numPm1; ++i) {
            this.mPolynomials[i] = [
                new Vector(this.mDimension), new Vector(this.mDimension),
                new Vector(this.mDimension), new Vector(this.mDimension)
            ];
        }
        this.mDelta = new Array<number>(numPm1).fill(0);
        for (let i0 = 0, i1 = 1; i1 < numPoints; i0 = i1++) {
            this.mDelta[i0] = times[i1] - times[i0];
        }

        const numPm2 = numPoints - 2;
        let boundary0: Vector;
        let boundary1: Vector;
        if (kind === 'clamped') {
            boundary0 = mul(this.mDelta[0], derivative0 as Vector);
            boundary1 = mul(this.mDelta[numPm2], derivative1 as Vector);
        }
        else {
            // Free splines and closed splines have the last two B-entries
            // set to the zero vector.
            boundary0 = new Vector(this.mDimension);
            boundary1 = new Vector(this.mDimension);
        }

        const R = new Matrix(3, 3);
        const numBElements = 3 * numPm1;
        const B = new Array<Vector>(numBElements);
        for (let i = 0; i < numBElements; ++i) {
            B[i] = new Vector(this.mDimension);
        }
        this.onPresolve(numPoints, f0, boundary0, boundary1, R, B);

        if (kind === 'free') {
            R.set(1, 1, 1);
            R.set(1, 2, 3);
            this.solve(0, 1, numPoints, f0, R, B);
        }
        else if (kind === 'closed') {
            const lambda = this.mDelta[0] / this.mDelta[numPm2];
            const lambdasqr = lambda * lambda;
            R.set(1, 0, -lambda);
            R.set(1, 1, -2 * lambda);
            R.set(1, 2, -3 * lambda);
            R.set(2, 1, -lambdasqr);
            R.set(2, 2, -3 * lambdasqr);
            this.solve(1, 1, numPoints, f0, R, B);
        }
        else {
            R.set(2, 0, 1);
            R.set(2, 1, 2);
            R.set(2, 2, 3);
            this.solve(1, 0, numPoints, f0, R, B);
        }

        this.mConstructed = true;
    }

    // The live array of polynomial segments; do not modify it.
    getPolynomials(): readonly NaturalCubicSplinePolynomial[] {
        return this.mPolynomials;
    }

    // Evaluation of the function and its derivatives through order 3. If you
    // want only the position, pass in order 0. If you want the position and
    // first derivative, pass in order of 1 and so on. The output array 'jet'
    // must have 'order + 1' elements. The values are ordered as position,
    // first derivative, second derivative and so on.
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
        jet[0] = add(poly[0],
            mul(u, add(poly[1], mul(u, add(poly[2], mul(u, poly[3]))))));
        if (order >= 1) {
            // Compute first derivative.
            let denom = this.mDelta[key];
            jet[1] = mul(add(poly[1],
                mul(u, add(mul(2, poly[2]), mul(u, mul(3, poly[3]))))),
                1 / denom);
            if (order >= 2) {
                // Compute second derivative.
                denom *= this.mDelta[key];
                jet[2] = mul(add(mul(2, poly[2]), mul(u, mul(6, poly[3]))),
                    1 / denom);
                if (order >= 3) {
                    // Compute third derivative.
                    denom *= this.mDelta[key];
                    jet[3] = mul(mul(6, poly[3]), 1 / denom);

                    for (let i = 4; i <= order; ++i) {
                        // Derivatives of order 4 and higher are zero.
                        jet[i] = new Vector(this.mDimension);
                    }
                }
            }
        }
    }

    private onPresolve(numPoints: number, f0: readonly Vector[],
        boundary0: Vector, boundary1: Vector, R: Matrix, B: Vector[]): void {
        const numPm1 = numPoints - 1;
        const numPm2 = numPoints - 2;
        const numPm3 = numPoints - 3;

        const coeff = [3, -3, 1];
        for (let i0 = 0, i1 = 1; i0 <= numPm3; i0 = i1++) {
            const diff = sub(f0[i1], f0[i0]);
            for (let j = 0, k = 3 * i0; j < 3; ++j, ++k) {
                B[k] = mul(coeff[j], diff);
            }
        }

        B[B.length - 3] = sub(f0[numPm1], f0[numPm2]);
        B[B.length - 2] = boundary0;
        B[B.length - 1] = boundary1;

        R.set(0, 0, 1);
        R.set(0, 1, 1);
        R.set(0, 2, 1);
    }

    private solve(ell10: number, ell21: number, numPoints: number,
        f0: readonly Vector[], R: Matrix, B: Vector[]): void {
        this.rowReduce(ell10, ell21, numPoints, R, B);
        this.backSubstitute(f0, R, B);
    }

    // Apply the row reductions that convert the matrix system to an
    // upper-triangular block-matrix system.
    private rowReduce(ell10: number, ell21: number, numPoints: number,
        R: Matrix, B: Vector[]): void {
        const numPm3 = numPoints - 3;

        if (ell10 === 1) {
            const trg = B.length - 2;
            B[trg] = sub(B[trg], B[0]);
            let sigma = this.mDelta[0] / this.mDelta[1];
            let sigmasqr = sigma * sigma;
            let luProd0 = 2 * sigma;
            let luProd1 = -sigmasqr;
            let sign = -1;

            for (let i = 1; i <= numPm3; ++i) {
                B[trg] = sub(B[trg], mul(sign,
                    add(mul(luProd0, B[3 * i]), mul(luProd1, B[3 * i + 1]))));
                sigma = this.mDelta[i] / this.mDelta[i + 1];
                sigmasqr = sigma * sigma;
                const temp0 = sigma * (2 * luProd0 - 3 * luProd1);
                const temp1 = sigmasqr * (-luProd0 + 2 * luProd1);
                luProd0 = temp0;
                luProd1 = temp1;
                sign = -sign;
            }

            R.set(1, 0, R.get(1, 0) + sign * luProd0);
            R.set(1, 1, R.get(1, 1) + sign * luProd1);
        }

        if (ell21 === 1) {
            const trg = B.length - 1;
            B[trg] = sub(B[trg], B[1]);
            let sigma = this.mDelta[0] / this.mDelta[1];
            let sigmasqr = sigma * sigma;
            let luProd0 = -3 * sigma;
            let luProd1 = 2 * sigmasqr;
            let sign = -1;

            for (let i = 1; i <= numPm3; ++i) {
                B[trg] = sub(B[trg], mul(sign,
                    add(mul(luProd0, B[3 * i]), mul(luProd1, B[3 * i + 1]))));
                sigma = this.mDelta[i] / this.mDelta[i + 1];
                sigmasqr = sigma * sigma;
                const temp0 = sigma * (2 * luProd0 - 3 * luProd1);
                const temp1 = sigmasqr * (-luProd0 + 2 * luProd1);
                luProd0 = temp0;
                luProd1 = temp1;
                sign = -sign;
            }

            R.set(2, 0, R.get(2, 0) + sign * luProd0);
            R.set(2, 1, R.get(2, 1) + sign * luProd1);
        }
    }

    private backSubstitute(f0: readonly Vector[], R: Matrix,
        B: readonly Vector[]): void {
        const inv = inverse3x3(R);
        logAssert(inv.invertible, 'R matrix is not invertible.');
        const invR = inv.inverse;

        const numPolynomials = this.mPolynomials.length;
        const poly = this.mPolynomials[numPolynomials - 1];
        let j0 = B.length - 3;
        let j1 = j0 + 1;
        let j2 = j0 + 2;

        poly[0] = f0[numPolynomials - 1].clone();
        poly[1] = add(add(mul(invR.get(0, 0), B[j0]),
            mul(invR.get(0, 1), B[j1])), mul(invR.get(0, 2), B[j2]));
        poly[2] = add(add(mul(invR.get(1, 0), B[j0]),
            mul(invR.get(1, 1), B[j1])), mul(invR.get(1, 2), B[j2]));
        poly[3] = add(add(mul(invR.get(2, 0), B[j0]),
            mul(invR.get(2, 1), B[j1])), mul(invR.get(2, 2), B[j2]));

        for (let i1 = numPolynomials - 2, i0 = i1 + 1; i1 >= 0; i0 = i1--) {
            const prev = this.mPolynomials[i0];
            const curr = this.mPolynomials[i1];
            const sigma = this.mDelta[i1] / this.mDelta[i0];
            const sigmasqr = sigma * sigma;
            const u00 = 2 * sigma;
            const u01 = -sigmasqr;
            const u10 = -3 * sigma;
            const u11 = 2 * sigmasqr;
            const u20 = sigma;
            const u21 = -sigmasqr;

            j0 -= 3;
            j1 -= 3;
            j2 -= 3;

            curr[0] = f0[i1].clone();
            curr[1] = sub(B[j0], add(mul(u00, prev[1]), mul(u01, prev[2])));
            curr[2] = sub(B[j1], add(mul(u10, prev[1]), mul(u11, prev[2])));
            curr[3] = sub(B[j2], add(mul(u20, prev[1]), mul(u21, prev[2])));
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
