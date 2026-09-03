// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ConvexHull2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the convex hull of 2D points using a divide-and-conquer algorithm.
// This is an O(N log N) algorithm for N input points. The only way to ensure
// a correct result for the input vertices is to use an exact predicate for
// computing signs of various expressions. The implementation uses interval
// arithmetic and rational arithmetic for the predicate.
//
// Port notes:
// * The upstream operator() functions become compute(points). The
//   (numPoints, points) overload is subsumed by the array argument.
// * Upstream selects Rational = BSNumber<UIntegerFP32<NumWords>> where
//   NumWords is a compile-time worst-case bound (18 for float, 132 for
//   double). The port's BSNumber is bigint-backed and grows as needed, so
//   the word-count bound is unnecessary and is dropped. The exact-arithmetic
//   behavior is identical (no divisions are performed, so BSNumber suffices).
// * The private enum class Order stays module-private; it is not part of the
//   public API upstream.
// * The mConverted/mRationalPoints memoization pair becomes a single array of
//   nullable rational points.

import { logAssert } from './Logger.js';
import { Line } from './Line.js';
import type { Line2 } from './Line.js';
import { Vector, sub, normalize } from './Vector.js';
import { BSNumber } from './BSNumber.js';
import { SWInterval } from './SWInterval.js';

// A rational 2D point, the port of Vector2<Rational>.
type RationalPoint2 = [BSNumber, BSNumber];

// An extended classification of the relationship of a point to a line
// segment. For noncollinear points, the return value is
//   POSITIVE when <P,Q0,Q1> is a counterclockwise triangle
//   NEGATIVE when <P,Q0,Q1> is a clockwise triangle
// For collinear points, the line direction is Q1-Q0. The return value is
//   COLLINEAR_LEFT when the line ordering is <P,Q0,Q1>
//   COLLINEAR_RIGHT when the line ordering is <Q0,Q1,P>
//   COLLINEAR_CONTAIN when the line ordering is <Q0,P,Q1>
enum Order {
    Q0_EQUALS_Q1,
    P_EQUALS_Q0,
    P_EQUALS_Q1,
    POSITIVE,
    NEGATIVE,
    COLLINEAR_LEFT,
    COLLINEAR_RIGHT,
    COLLINEAR_CONTAIN
}

export class ConvexHull2 {
    // If the dimension is 0 or 1, compute() returns false. The caller is
    // responsible for retrieving the dimension and taking an alternate path
    // should the dimension be smaller than 2. If the dimension is 0, the
    // points[] are all the same point. If the dimension is 1, the caller can
    // query for the approximating line and project points[] onto it for
    // further processing.
    private mDimension: number;
    private mLine: Line2;

    // The array of rational points used for the exact predicate. An entry is
    // null until the corresponding floating-point point is first needed by a
    // predicate computation, after which the conversion is memoized to avoid
    // converting again.
    private mRationalPoints: (RationalPoint2 | null)[];

    private mNumPoints: number;
    private mNumUniquePoints: number;
    private mPoints: readonly Vector[];
    private mMerged: number[];
    private mHull: number[];

    constructor() {
        this.mDimension = 0;
        this.mLine = new Line(2);
        this.mLine.direction.makeZero();
        this.mRationalPoints = [];
        this.mNumPoints = 0;
        this.mNumUniquePoints = 0;
        this.mPoints = [];
        this.mMerged = [];
        this.mHull = [];
    }

    // The input is the array of points whose convex hull is required. The
    // return value is 'true' if the hull is 2-dimensional. It is 'false' if
    // the hull is 0-dimensional or 1-dimensional. The input must have at
    // least 1 element.
    compute(points: readonly Vector[]): boolean {
        logAssert(points.length > 0, 'Invalid input to ConvexHull2 compute().');
        for (const point of points) {
            logAssert(point.size === 2, 'ConvexHull2 requires 2D points.');
        }

        this.mDimension = 0;
        this.mLine.origin.makeZero();
        this.mLine.direction.makeZero();
        this.mNumPoints = points.length;
        this.mNumUniquePoints = 0;
        this.mPoints = points;
        this.mMerged = [];
        this.mHull = [];

        // Allocate storage for any rational points that must be computed in
        // the exact predicate.
        this.mRationalPoints = new Array<RationalPoint2 | null>(this.mNumPoints).fill(null);

        // Sort the points indirectly. The mHull array is used to store the
        // unique indices.
        const hull: number[] = new Array<number>(this.mNumPoints);
        for (let i = 0; i < this.mNumPoints; ++i) {
            hull[i] = i;
        }
        hull.sort((s0, s1) => {
            const p0 = points[s0].values;
            const p1 = points[s1].values;
            if (p0[0] !== p1[0]) {
                return p0[0] < p1[0] ? -1 : 1;
            }
            if (p0[1] !== p1[1]) {
                return p0[1] < p1[1] ? -1 : 1;
            }
            return 0;
        });

        // The port of std::unique with the equal-points predicate.
        this.mHull = [];
        for (let i = 0; i < hull.length; ++i) {
            if (i === 0 || !points[hull[i]].equals(points[this.mHull[this.mHull.length - 1]])) {
                this.mHull.push(hull[i]);
            }
        }
        this.mNumUniquePoints = this.mHull.length;

        // Use a divide-and-conquer algorithm. The merge step computes the
        // convex hull of two convex polygons. The merge storage is allocated
        // once to avoid reallocations during the recursive chain of the
        // getHull and merge member functions.
        this.mMerged = new Array<number>(this.mNumUniquePoints).fill(0);
        const range = { i0: 0, i1: this.mNumUniquePoints - 1 };
        this.getHullRecursive(range);
        const hullSize = range.i1 - range.i0 + 1;
        this.mHull.length = hullSize;
        if (hullSize === 1) {
            // The input points are all the same point.
            this.mDimension = 0;
            return false;
        }
        else if (hullSize === 2) {
            // The input points are collinear.
            this.mDimension = 1;
            this.mLine.origin = points[this.mHull[0]].clone();
            this.mLine.direction = sub(points[this.mHull[1]], points[this.mHull[0]]);
            normalize(this.mLine.direction);
            return false;
        }
        else {  // hullSize > 2
            this.mDimension = 2;
            return true;
        }
    }

    // The dimension is 0 (hull is a single point), 1 (hull is a line
    // segment) or 2 (hull is a convex polygon).
    getDimension(): number {
        return this.mDimension;
    }

    // When dimension is 1, the line is a floating-point approximation to the
    // line containing the hull points.
    getLine(): Line2 {
        return this.mLine;
    }

    // Member access. getNumPoints() returns the number of elements of the
    // points[] array passed to compute(). getPoints() returns the points.
    // getNumUniquePoints() returns the number of unique points in points[].
    getNumPoints(): number {
        return this.mNumPoints;
    }

    getNumUniquePoints(): number {
        return this.mNumUniquePoints;
    }

    getPoints(): readonly Vector[] {
        return this.mPoints;
    }

    // Get the indices into the input 'points[]' that correspond to hull
    // vertices. The returned array is organized according to the hull
    // dimension.
    //   0: The hull is a single point. The returned array has size 1 with
    //      index corresponding to that point.
    //   1: The hull is a line segment. The returned array has size 2 with
    //      indices corresponding to the segment endpoints.
    //   2: The hull is a convex polygon. The returned array has size N with
    //      indices corresponding to the polygon vertices. The vertices are
    //      counterclockwise ordered.
    getHull(): readonly number[] {
        return this.mHull;
    }

    // Support for divide-and-conquer. The upstream (int32_t& i0, int32_t& i1)
    // in/out reference pair becomes a mutable range object.
    private getHullRecursive(range: { i0: number; i1: number }): void {
        const numVertices = range.i1 - range.i0 + 1;
        if (numVertices > 1) {
            // Compute the middle index of input range.
            const mid = Math.floor((range.i0 + range.i1) / 2);

            // Compute the hull of subsets (mid-i0+1 >= i1-mid).
            const left = { i0: range.i0, i1: mid };
            const right = { i0: mid + 1, i1: range.i1 };
            this.getHullRecursive(left);
            this.getHullRecursive(right);

            // Merge the convex hulls into a single convex hull.
            this.merge(left.i0, left.i1, right.i0, right.i1, range);
        }
        // else: The convex hull is a single point.
    }

    private merge(j0: number, j1: number, j2: number, j3: number,
        range: { i0: number; i1: number }): void {
        // Subhull0 is to the left of subhull1 because of the initial sorting
        // of the points by x-components. We need to find two mutually visible
        // points, one on the left subhull and one on the right subhull.
        const size0 = j1 - j0 + 1;
        const size1 = j3 - j2 + 1;

        let i: number;
        let p: Vector;

        // Find the right-most point of the left subhull.
        let pmax0 = this.mPoints[this.mHull[j0]];
        let imax0 = j0;
        for (i = j0 + 1; i <= j1; ++i) {
            p = this.mPoints[this.mHull[i]];
            if (pmax0.lessThan(p)) {
                pmax0 = p;
                imax0 = i;
            }
        }

        // Find the left-most point of the right subhull.
        let pmin1 = this.mPoints[this.mHull[j2]];
        let imin1 = j2;
        for (i = j2 + 1; i <= j3; ++i) {
            p = this.mPoints[this.mHull[i]];
            if (p.lessThan(pmin1)) {
                pmin1 = p;
                imin1 = i;
            }
        }

        // Get the lower tangent to hulls (LL = lower-left, LR = lower-right).
        const lower = { i0: imax0, i1: imin1 };
        this.getTangent(j0, j1, j2, j3, lower);
        const iLL = lower.i0, iLR = lower.i1;

        // Get the upper tangent to hulls (UL = upper-left, UR = upper-right).
        const upper = { i0: imin1, i1: imax0 };
        this.getTangent(j2, j3, j0, j1, upper);
        const iUR = upper.i0, iUL = upper.i1;

        // Construct the counterclockwise-ordered merged-hull vertices.
        let numMerged = 0;
        let k: number;

        i = iUL;
        for (k = 0; k < size0; ++k) {
            this.mMerged[numMerged++] = this.mHull[i];
            if (i === iLL) {
                break;
            }
            i = (i < j1 ? i + 1 : j0);
        }
        logAssert(k < size0, 'Unexpected condition.');

        i = iLR;
        for (k = 0; k < size1; ++k) {
            this.mMerged[numMerged++] = this.mHull[i];
            if (i === iUR) {
                break;
            }
            i = (i < j3 ? i + 1 : j2);
        }
        logAssert(k < size1, 'Unexpected condition.');

        let next = j0;
        for (k = 0; k < numMerged; ++k) {
            this.mHull[next] = this.mMerged[k];
            ++next;
        }

        range.i0 = j0;
        range.i1 = next - 1;
    }

    private getTangent(j0: number, j1: number, j2: number, j3: number,
        range: { i0: number; i1: number }): void {
        // The loop terminates in a finite number of steps, but the upper
        // bound for the loop variable is used as a guard against an infinite
        // loop. The infinite loop should not occur because rational
        // arithmetic is used in toLineExtended.
        const size0 = j1 - j0 + 1;
        const size1 = j3 - j2 + 1;
        const imax = size0 + size1;

        for (let i = 0; i < imax; ++i) {
            // Get the endpoints of the potential tangent.
            const L1index = this.mHull[range.i0];
            const R0index = this.mHull[range.i1];

            // Walk along the left hull to find the point of tangency.
            if (size0 > 1) {
                const iLm1 = (range.i0 > j0 ? range.i0 - 1 : j1);
                const L0index = this.mHull[iLm1];
                const order = this.toLineExtended(R0index, L0index, L1index);
                if (order === Order.NEGATIVE || order === Order.COLLINEAR_RIGHT) {
                    range.i0 = iLm1;
                    continue;
                }
            }

            // Walk along right hull to find the point of tangency.
            if (size1 > 1) {
                const iRp1 = (range.i1 < j3 ? range.i1 + 1 : j2);
                const R1index = this.mHull[iRp1];
                const order = this.toLineExtended(L1index, R0index, R1index);
                if (order === Order.NEGATIVE || order === Order.COLLINEAR_LEFT) {
                    range.i1 = iRp1;
                    continue;
                }
            }

            // The tangent segment has been found.
            break;
        }
    }

    // Memoized access to the rational representation of the points.
    private getRationalPoint(index: number): RationalPoint2 {
        let rPoint = this.mRationalPoints[index];
        if (rPoint === null) {
            const point = this.mPoints[index].values;
            rPoint = [BSNumber.fromNumber(point[0]), BSNumber.fromNumber(point[1])];
            this.mRationalPoints[index] = rPoint;
        }
        return rPoint;
    }

    private toLineExtended(pIndex: number, q0Index: number, q1Index: number): Order {
        const P = this.mPoints[pIndex].values;
        const Q0 = this.mPoints[q0Index].values;
        const Q1 = this.mPoints[q1Index].values;

        if (Q1[0] === Q0[0] && Q1[1] === Q0[1]) {
            return Order.Q0_EQUALS_Q1;
        }

        if (P[0] === Q0[0] && P[1] === Q0[1]) {
            return Order.P_EQUALS_Q0;
        }

        if (P[0] === Q1[0] && P[1] === Q1[1]) {
            return Order.P_EQUALS_Q1;
        }

        // The theoretical classification relies on computing exactly the sign
        // of the determinant. Numerical roundoff errors can cause
        // misclassification.
        const zero = 0;
        const ip0 = new SWInterval(P[0]), ip1 = new SWInterval(P[1]);
        const iq00 = new SWInterval(Q0[0]), iq01 = new SWInterval(Q0[1]);
        const iq10 = new SWInterval(Q1[0]), iq11 = new SWInterval(Q1[1]);
        const ix0 = iq10.sub(iq00), iy0 = iq11.sub(iq01);
        const ix1 = ip0.sub(iq00), iy1 = ip1.sub(iq01);
        const ix0y1 = ix0.mul(iy1);
        const ix1y0 = ix1.mul(iy0);
        const iDet = ix0y1.sub(ix1y0);
        let sign: number;

        let rDiff0: RationalPoint2 | null = null;
        let rDiff1: RationalPoint2 | null = null;
        let rDot: BSNumber | null = null;

        if (iDet.get(0) > zero) {
            sign = +1;
        }
        else if (iDet.get(1) < zero) {
            sign = -1;
        }
        else {
            // The exact sign of the determinant is not known, so compute the
            // determinant using rational arithmetic.
            const rP = this.getRationalPoint(pIndex);
            const rQ0 = this.getRationalPoint(q0Index);
            const rQ1 = this.getRationalPoint(q1Index);
            rDiff0 = [rQ1[0].sub(rQ0[0]), rQ1[1].sub(rQ0[1])];
            rDiff1 = [rP[0].sub(rQ0[0]), rP[1].sub(rQ0[1])];
            // DotPerp(rDiff0, rDiff1)
            const rDet = rDiff0[0].mul(rDiff1[1]).sub(rDiff0[1].mul(rDiff1[0]));
            sign = rDet.getSign();
        }

        if (sign > 0) {
            // The points form a counterclockwise triangle <P,Q0,Q1>.
            return Order.POSITIVE;
        }
        else if (sign < 0) {
            // The points form a clockwise triangle <P,Q1,Q0>.
            return Order.NEGATIVE;
        }
        else {
            // The points are collinear. P is on the line through Q0 and Q1.
            const iDot = ix0.mul(ix1).add(iy0.mul(iy1));
            if (iDot.get(0) > zero) {
                sign = +1;
            }
            else if (iDot.get(1) < zero) {
                sign = -1;
            }
            else {
                // The exact sign of the dot product is not known, so compute
                // the dot product using rational arithmetic.
                const rP = this.getRationalPoint(pIndex);
                const rQ0 = this.getRationalPoint(q0Index);
                const rQ1 = this.getRationalPoint(q1Index);
                if (rDiff0 === null) {
                    rDiff0 = [rQ1[0].sub(rQ0[0]), rQ1[1].sub(rQ0[1])];
                }
                if (rDiff1 === null) {
                    rDiff1 = [rP[0].sub(rQ0[0]), rP[1].sub(rQ0[1])];
                }
                rDot = rDiff0[0].mul(rDiff1[0]).add(rDiff0[1].mul(rDiff1[1]));
                sign = rDot.getSign();
            }

            if (sign < 0) {
                // The line ordering is <P,Q0,Q1>.
                return Order.COLLINEAR_LEFT;
            }

            const iSqrLength = ix0.mul(ix0).add(iy0.mul(iy0));
            const iTest = iDot.sub(iSqrLength);
            if (iTest.get(0) > zero) {
                sign = +1;
            }
            else if (iTest.get(1) < zero) {
                sign = -1;
            }
            else {
                // The exact sign of the test is not known, so compute the
                // test using rational arithmetic.
                const rP = this.getRationalPoint(pIndex);
                const rQ0 = this.getRationalPoint(q0Index);
                const rQ1 = this.getRationalPoint(q1Index);
                if (rDiff0 === null) {
                    rDiff0 = [rQ1[0].sub(rQ0[0]), rQ1[1].sub(rQ0[1])];
                }
                if (rDiff1 === null) {
                    rDiff1 = [rP[0].sub(rQ0[0]), rP[1].sub(rQ0[1])];
                }
                if (rDot === null) {
                    rDot = rDiff0[0].mul(rDiff1[0]).add(rDiff0[1].mul(rDiff1[1]));
                }
                const rSqrLength = rDiff0[0].mul(rDiff0[0]).add(rDiff0[1].mul(rDiff0[1]));
                const rTest = rDot.sub(rSqrLength);
                sign = rTest.getSign();
            }

            if (sign > 0) {
                // The line ordering is <Q0,Q1,P>.
                return Order.COLLINEAR_RIGHT;
            }

            // The line ordering is <Q0,P,Q1> with P strictly between Q0 and
            // Q1.
            return Order.COLLINEAR_CONTAIN;
        }
    }
}
