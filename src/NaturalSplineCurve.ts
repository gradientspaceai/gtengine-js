// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) NaturalSplineCurve.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// NOTE (upstream): this class is deprecated and will not be ported to GTL.
// Use instead NaturalCubicSpline; there is also an extension of the idea in
// NaturalQuinticSpline. It is ported here for completeness of the GTE
// surface area.
//
// The curve is a C2 cubic spline that interpolates the input points at the
// input times. On segment i the curve is
//   X(t) = A[i] + dt*(B[i] + dt*(C[i] + dt*D[i])),  dt = t - times[i]
// A free spline has zero second derivative at the endpoints, a closed spline
// wraps around (the first and last points are treated as adjacent), and a
// clamped spline has the specified first derivatives at the endpoints.
//
// Port notes:
// - Upstream 'template <int32_t N, typename Real>' becomes a runtime
//   dimension (the ParametricCurve/BezierCurve precedent from B40). The
//   dimension is taken from points[0]; all points must have that dimension.
// - Upstream has two constructors, (isFree, numPoints, points, times) and
//   (numPoints, points, times, derivative0, derivative1). Following the
//   NaturalCubicSpline precedent (B85), the port exposes three static
//   factories instead: createFree, createClosed and createClamped.
// - Upstream packs A, B, C and D into one std::vector<Vector<N,Real>> of
//   4*numPoints-2 entries and aliases four raw pointers into it. The port
//   uses four separate arrays with the same lengths: A has numPoints
//   entries, B and D have numSegments entries, and C has numSegments+1
//   entries.
// - The upstream 'Evaluate' becomes 'evaluate'; 'GetKeyInfo' becomes the
//   protected 'getKeyInfo', which returns { key, dt } instead of using
//   output reference parameters.

import { LinearSystem } from './LinearSystem';
import { logAssert } from './Logger';
import { ParametricCurve } from './ParametricCurve';
import { Vector, add, div, mul, sub } from './Vector';

// Which boundary conditions the spline satisfies.
export enum NaturalSplineCurveKind {
    FREE,
    CLOSED,
    CLAMPED
}

export class NaturalSplineCurve extends ParametricCurve {
    // Polynomial coefficients. mA are the points (constant coefficients of
    // the polynomials), mB are the degree 1 coefficients, mC are the degree 2
    // coefficients and mD are the degree 3 coefficients.
    protected mNumPoints: number;
    protected mNumSegments: number;
    protected mA: Vector[];
    protected mB: Vector[];
    protected mC: Vector[];
    protected mD: Vector[];

    // A spline with second derivatives zero at the endpoints.
    static createFree(points: readonly Vector[],
        times: readonly number[]): NaturalSplineCurve {
        return new NaturalSplineCurve(NaturalSplineCurveKind.FREE, points,
            times, null, null);
    }

    // A closed spline.
    static createClosed(points: readonly Vector[],
        times: readonly number[]): NaturalSplineCurve {
        return new NaturalSplineCurve(NaturalSplineCurveKind.CLOSED, points,
            times, null, null);
    }

    // A clamped spline with the specified first derivatives at the
    // endpoints. Usually, derivative0 = points[1] - points[0] at the first
    // point and derivative1 = points[M-1] - points[M-2] at the last point.
    static createClamped(points: readonly Vector[], times: readonly number[],
        derivative0: Vector, derivative1: Vector): NaturalSplineCurve {
        return new NaturalSplineCurve(NaturalSplineCurveKind.CLAMPED, points,
            times, derivative0, derivative1);
    }

    private constructor(kind: NaturalSplineCurveKind,
        points: readonly Vector[], times: readonly number[],
        derivative0: Vector | null, derivative1: Vector | null) {
        logAssert(points !== null && points !== undefined &&
            times !== null && times !== undefined &&
            points.length >= 2 && times.length >= points.length,
            'Invalid input.');

        super(points[0].size, points.length - 1, times);

        const dimension = points[0].size;
        this.mNumPoints = points.length;
        this.mNumSegments = this.mNumPoints - 1;
        this.mA = new Array<Vector>(this.mNumPoints);
        this.mB = new Array<Vector>(this.mNumSegments);
        this.mC = new Array<Vector>(this.mNumSegments + 1);
        this.mD = new Array<Vector>(this.mNumSegments);
        for (let i = 0; i < this.mNumPoints; ++i) {
            logAssert(points[i].size === dimension, 'Invalid input.');
            this.mA[i] = points[i].clone();
        }
        for (let i = 0; i < this.mNumSegments; ++i) {
            this.mB[i] = new Vector(dimension);
            this.mD[i] = new Vector(dimension);
        }
        for (let i = 0; i <= this.mNumSegments; ++i) {
            this.mC[i] = new Vector(dimension);
        }

        if (kind === NaturalSplineCurveKind.FREE) {
            this.createFreeInternal();
        }
        else if (kind === NaturalSplineCurveKind.CLOSED) {
            this.createClosedInternal();
        }
        else {
            logAssert(derivative0 !== null && derivative1 !== null &&
                derivative0.size === dimension &&
                derivative1.size === dimension, 'Invalid input.');
            this.createClampedInternal(derivative0, derivative1);
        }

        this.mConstructed = true;
    }

    // Member access.
    getNumPoints(): number {
        return this.mNumPoints;
    }

    // The live array of interpolation points (upstream returns the raw mA
    // pointer); it must not be modified by the caller.
    getPoints(): readonly Vector[] {
        return this.mA;
    }

    // Evaluation of the function and its derivatives through order 3. If you
    // want only the position, pass in order 0. If you want the position and
    // first derivative, pass in an order of 1 and so on. The output array
    // 'jet' must have 'order + 1' elements. The values are ordered as
    // position, first derivative, second derivative and so on.
    override evaluate(t: number, order: number, jet: Vector[]): void {
        if (!this.mConstructed) {
            // Return a zero-valued jet for invalid state.
            for (let i = 0; i <= order; ++i) {
                jet[i].makeZero();
            }
            return;
        }

        const info = this.getKeyInfo(t);
        const key = info.key;
        const dt = info.dt;

        // Compute position.
        jet[0] = add(this.mA[key], mul(dt, add(this.mB[key],
            mul(dt, add(this.mC[key], mul(dt, this.mD[key]))))));
        if (order >= 1) {
            // Compute first derivative.
            jet[1] = add(this.mB[key], mul(dt, add(mul(2, this.mC[key]),
                mul(3 * dt, this.mD[key]))));
            if (order >= 2) {
                // Compute second derivative.
                jet[2] = add(mul(2, this.mC[key]), mul(6 * dt, this.mD[key]));
                if (order >= 3) {
                    jet[3] = mul(6, this.mD[key]);

                    for (let i = 4; i <= order; ++i) {
                        jet[i].makeZero();
                    }
                }
            }
        }
    }

    protected createFreeInternal(): void {
        const numS = this.mNumSegments;
        const numSm1 = numS - 1;
        const dimension = this.mDimension;

        const dt = new Array<number>(numS).fill(0);
        const d2t = new Array<number>(numS).fill(0);
        const alpha = new Array<Vector>(numS);
        const ell = new Array<number>(numS + 1).fill(0);
        const mu = new Array<number>(numS).fill(0);
        const z = new Array<Vector>(numS + 1);
        for (let i = 0; i < numS; ++i) {
            alpha[i] = new Vector(dimension);
        }
        for (let i = 0; i <= numS; ++i) {
            z[i] = new Vector(dimension);
        }

        for (let i = 0, ip1 = 1; i < numS; ++i, ++ip1) {
            dt[i] = this.mTime[ip1] - this.mTime[i];
        }

        // d2t[0] is unused.
        for (let im1 = 0, i = 1, ip1 = 2; i < numS; im1 = i, i = ip1++) {
            d2t[i] = this.mTime[ip1] - this.mTime[im1];
        }

        // alpha[0] is unused.
        for (let im1 = 0, i = 1, ip1 = 2; i < numS; im1 = i, i = ip1++) {
            const numer = mul(3, add(sub(mul(dt[im1], this.mA[ip1]),
                mul(d2t[i], this.mA[i])), mul(dt[i], this.mA[im1])));
            const denom = dt[im1] * dt[i];
            alpha[i] = div(numer, denom);
        }

        ell[0] = 1;
        mu[0] = 0;
        // z[0] is already zero.
        for (let im1 = 0, i = 1; i < numS; im1 = i++) {
            ell[i] = 2 * d2t[i] - dt[im1] * mu[im1];
            mu[i] = dt[i] / ell[i];
            z[i] = div(sub(alpha[i], mul(dt[im1], z[im1])), ell[i]);
        }
        ell[numS] = 1;
        // z[numS] is already zero.

        this.mC[numS].makeZero();
        for (let j = 0, i = numSm1; j < numS; ++j, --i) {
            this.mC[i] = sub(z[i], mul(mu[i], this.mC[i + 1]));
            this.mB[i] = sub(div(sub(this.mA[i + 1], this.mA[i]), dt[i]),
                div(mul(dt[i], add(this.mC[i + 1], mul(2, this.mC[i]))), 3));
            this.mD[i] = div(sub(this.mC[i + 1], this.mC[i]), 3 * dt[i]);
        }
    }

    protected createClosedInternal(): void {
        const numP = this.mNumPoints;
        const numS = this.mNumSegments;
        const numSm1 = numS - 1;
        const dimension = this.mDimension;

        const dt = new Array<number>(numS).fill(0);
        // The matrix of the system, in row-major order.
        const mat = new Array<number>(numP * numP).fill(0);

        for (let i = 0, ip1 = 1; i < numS; ++i, ++ip1) {
            dt[i] = this.mTime[ip1] - this.mTime[i];
        }

        // Construct the matrix of the system.
        mat[0 + numP * 0] = 1;
        mat[numS + numP * 0] = -1;
        for (let im1 = 0, i = 1, ip1 = 2; i <= numSm1; im1 = i, i = ip1++) {
            mat[im1 + numP * i] = dt[im1];
            mat[i + numP * i] = 2 * (dt[im1] + dt[i]);
            mat[ip1 + numP * i] = dt[i];
        }
        // Upstream bug fix. Upstream assigns these three entries of the
        // wrap-around row:
        //     mat(numS,numSm1) = dt[numSm1];
        //     mat(numS,0)      = 2*(dt[numSm1] + dt[0]);
        //     mat(numS,1)      = dt[0];
        // The three columns numSm1, 0 and 1 are distinct only when
        // numSegments >= 3. For numSegments == 2 the columns numSm1 and 1
        // coincide, so the assignment of dt[0] silently discards the
        // dt[numSm1] term and the solved system is wrong (and similarly for
        // numSegments == 1, where columns numSm1 and 0 coincide). Because the
        // wrap-around equation genuinely involves the same unknown twice in
        // those cases, the terms must accumulate. 'mat' is zero-filled, so
        // for numSegments >= 3 this is identical to upstream.
        mat[numSm1 + numP * numS] += dt[numSm1];
        mat[0 + numP * numS] += 2 * (dt[numSm1] + dt[0]);
        mat[1 + numP * numS] += dt[0];

        // Construct the right-hand side of the system.
        this.mC[0].makeZero();
        for (let im1 = 0, i = 1, ip1 = 2; ip1 <= numS; im1 = i, i = ip1++) {
            this.mC[i] = mul(3, sub(
                div(sub(this.mA[ip1], this.mA[i]), dt[i]),
                div(sub(this.mA[i], this.mA[im1]), dt[im1])));
        }
        this.mC[numS] = mul(3, sub(
            div(sub(this.mA[1], this.mA[0]), dt[0]),
            div(sub(this.mA[0], this.mA[numSm1]), dt[numSm1])));

        // Solve the linear systems. The right-hand side is the numP-by-N
        // row-major array of the C values.
        const rhs = new Array<number>(numP * dimension).fill(0);
        for (let i = 0, k = 0; i <= numS; ++i) {
            for (let j = 0; j < dimension; ++j, ++k) {
                rhs[k] = this.mC[i].values[j];
            }
        }
        const result = LinearSystem.solveMultiple(numP, dimension, mat, rhs);
        logAssert(result.invertible, 'Failed to solve linear system.');

        for (let i = 0, k = 0; i <= numS; ++i) {
            for (let j = 0; j < dimension; ++j, ++k) {
                this.mC[i].values[j] = result.X[k];
            }
        }

        for (let i = 0; i < numS; ++i) {
            this.mB[i] = sub(div(sub(this.mA[i + 1], this.mA[i]), dt[i]),
                div(mul(add(this.mC[i + 1], mul(2, this.mC[i])), dt[i]), 3));
            this.mD[i] = div(sub(this.mC[i + 1], this.mC[i]), 3 * dt[i]);
        }
    }

    protected createClampedInternal(derivative0: Vector,
        derivative1: Vector): void {
        const numS = this.mNumSegments;
        const numSm1 = numS - 1;
        const dimension = this.mDimension;

        const dt = new Array<number>(numS).fill(0);
        const d2t = new Array<number>(numS).fill(0);
        const alpha = new Array<Vector>(numS + 1);
        const ell = new Array<number>(numS + 1).fill(0);
        const mu = new Array<number>(numS).fill(0);
        const z = new Array<Vector>(numS + 1);
        for (let i = 0; i <= numS; ++i) {
            alpha[i] = new Vector(dimension);
            z[i] = new Vector(dimension);
        }

        for (let i = 0, ip1 = 1; i < numS; i = ip1++) {
            dt[i] = this.mTime[ip1] - this.mTime[i];
        }

        for (let im1 = 0, i = 1, ip1 = 2; i < numS; im1 = i, i = ip1++) {
            d2t[i] = this.mTime[ip1] - this.mTime[im1];
        }

        alpha[0] = mul(3, sub(div(sub(this.mA[1], this.mA[0]), dt[0]),
            derivative0));
        alpha[numS] = mul(3, sub(derivative1,
            div(sub(this.mA[numS], this.mA[numS - 1]), dt[numS - 1])));
        for (let im1 = 0, i = 1, ip1 = 2; i < numS; im1 = i, i = ip1++) {
            const numer = mul(3, add(sub(mul(dt[im1], this.mA[ip1]),
                mul(d2t[i], this.mA[i])), mul(dt[i], this.mA[im1])));
            const denom = dt[im1] * dt[i];
            alpha[i] = div(numer, denom);
        }

        ell[0] = 2 * dt[0];
        mu[0] = 0.5;
        z[0] = div(alpha[0], ell[0]);
        for (let im1 = 0, i = 1; i < numS; im1 = i++) {
            ell[i] = 2 * d2t[i] - dt[im1] * mu[im1];
            mu[i] = dt[i] / ell[i];
            z[i] = div(sub(alpha[i], mul(dt[im1], z[im1])), ell[i]);
        }
        ell[numS] = dt[numSm1] * (2 - mu[numSm1]);
        z[numS] = div(sub(alpha[numS], mul(dt[numSm1], z[numSm1])), ell[numS]);

        this.mC[numS] = z[numS].clone();
        for (let j = 0, i = numSm1; j < numS; ++j, --i) {
            this.mC[i] = sub(z[i], mul(mu[i], this.mC[i + 1]));
            this.mB[i] = sub(div(sub(this.mA[i + 1], this.mA[i]), dt[i]),
                div(mul(dt[i], add(this.mC[i + 1], mul(2, this.mC[i]))), 3));
            this.mD[i] = div(sub(this.mC[i + 1], this.mC[i]), 3 * dt[i]);
        }
    }

    // Determine the index i for which times[i] <= t < times[i+1].
    protected getKeyInfo(t: number): { key: number, dt: number } {
        if (t <= this.mTime[0]) {
            return { key: 0, dt: 0 };
        }

        if (t >= this.mTime[this.mNumSegments]) {
            return {
                key: this.mNumSegments - 1,
                dt: this.mTime[this.mNumSegments] -
                    this.mTime[this.mNumSegments - 1]
            };
        }

        for (let i = 0, ip1 = 1; i < this.mNumSegments; i = ip1++) {
            if (t < this.mTime[ip1]) {
                return { key: i, dt: t - this.mTime[i] };
            }
        }

        // Unreachable given the branches above, but a return value is
        // required (as upstream's loop-with-break requires 'key' and 'dt' to
        // have been assigned).
        return { key: 0, dt: 0 };
    }
}
