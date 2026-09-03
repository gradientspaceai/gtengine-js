// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ConvexHullSimplePolygon.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the convex hull of a simple polygon. The implementation is for
// the algorithm published in
//   On-line construction of the convex hull of a simple polyline
//   Avraham A. Melkman
//   Information Processing Letters 25 (1987), pages 11-12
//   North Holland Publishing Co.
//
// A freely downloadable copy can be obtained by
//   https://www.ime.usp.br/~walterfm/cursos/mac0331/2006/melkman.pdf
//
// A related webpage with a description of algorithm details is
//   https://cgm.cs.mcgill.ca/~athens/cs601/Melkman.html
//
// Port notes: the upstream operator() becomes compute(...) and the upstream
// output parameter 'std::vector<std::size_t>& hull' becomes the return value.
// The C++ std::deque becomes a plain array used with push/pop/unshift/shift;
// the polygon sizes for which this class is used are small enough that the
// O(n) cost of unshift/shift is not a concern relative to the clarity of
// matching the upstream code.

import { logAssert } from './Logger.js';
import { Vector, sub } from './Vector.js';
import { dotPerp } from './Vector2.js';

export class ConvexHullSimplePolygon {
    // The polygon must be counterclockwise ordered, because the Minkowski
    // sum of convex polygon and disk assumes counterclockwise ordering.
    // The returned hull[] is an ordered list of indices into polygon[].
    // The hull vertices are
    //   { polygon[hull[0]], ..., polygon[hull[hull.length - 1]] }
    // and the hull is counterclockwise ordered.
    compute(polygon: readonly Vector[]): number[] {
        const n = polygon.length;
        logAssert(n >= 3, 'The input polygon must have at least 3 vertices.');

        // Melkman's algorithm converted to use an array as the double-ended
        // queue rather than the double-ended queue in the pseudocode of the
        // paper. Notice that the array does not require the 'b' and 't'
        // indices used in the pseudocode of the paper.
        const dq: number[] = [];
        if (this.whichSide(polygon, 0, 1, 2) > 0) {
            dq.push(0);
            dq.push(1);
        }
        else {
            dq.push(1);
            dq.push(0);
        }
        dq.push(2);
        dq.unshift(2);

        let i = 2;
        for (;;) {
            if (++i >= n) {
                break;
            }

            // The incrementing of i in "++i < n" and the following block of
            // code are Step 2 of the PDF pseudocode. The author's comment
            // before the pseudocode is: "The algorithm halts when its input
            // is exhausted." It is unclear whether the exhaustion occurs in
            // "while (++i < n)" or in the incrementing of i inside the
            // while-loop below (or perhaps it can occur in either based on
            // the input polygon). Just to be safe, range checking is
            // performed after the while-loop below terminates.
            while (
                i < n &&
                this.whichSide(polygon, i, dq[0], dq[1]) >= 0 &&
                this.whichSide(polygon, dq[dq.length - 2], dq[dq.length - 1], i) >= 0) {
                ++i;
            }
            if (i === n) {
                break;
            }

            // This block of code is Step 3 of the PDF pseudocode.
            while (this.whichSide(polygon, dq[dq.length - 2], dq[dq.length - 1], i) <= 0) {
                dq.pop();
            }
            dq.push(i);

            // This block of code is Step 4 of the PDF pseudocode.
            while (this.whichSide(polygon, i, dq[0], dq[1]) <= 0) {
                dq.shift();
            }
            dq.unshift(i);
        }

        return this.hullFromDoubleEndedQueue(dq);
    }

    // Given directed edge <p0,p1>, determine which side of the line of the
    // directed edge contains the point p2. The function returns
    //   +1: p2 is on the right of the line
    //    0: p2 is on the line (p0, p1, and p2 are colinear)
    //   -1: p2 is on the left of the line
    private whichSide(polygon: readonly Vector[], i0: number, i1: number,
        i2: number): number {
        const zero = 0;
        const diff10 = sub(polygon[i1], polygon[i0]);
        const diff20 = sub(polygon[i2], polygon[i0]);
        const test = dotPerp(diff20, diff10);
        return (test > zero ? +1 : (test < zero ? -1 : 0));
    }

    // Copy the double-ended queue into an array, reversing the order.
    private hullFromDoubleEndedQueue(dq: readonly number[]): number[] {
        // Guard against a negative hull size.
        logAssert(dq.length >= 2, 'Invalid double-ended queue size.');

        const lastIndex = dq.length - 1;
        const hull = new Array<number>(lastIndex);
        for (let i0 = 0, i1 = lastIndex; i0 < hull.length; ++i0, --i1) {
            hull[i0] = dq[i1];
        }
        return hull;
    }
}
