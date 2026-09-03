// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContEllipse2MinCR.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the minimum-area ellipse, (X-C)^T R D R^T (X-C) = 1, given the
// center C and the orientation matrix R. The columns of R are the axes of
// the ellipse. The algorithm computes the diagonal matrix D. The minimum
// area is pi/sqrt(D[0]*D[1]), where D = diag(D[0],D[1]). The problem is
// equivalent to maximizing the product D[0]*D[1] for a given C and R, and
// subject to the constraints
//   (P[i]-C)^T R D R^T (P[i]-C) <= 1
// for all input points P[i] with 0 <= i < N. Each constraint has the form
//   A[0]*D[0] + A[1]*D[1] <= 1
// where A[0] >= 0 and A[1] >= 0.
//
// Port notes: upstream's functor class ContEllipse2MinCR becomes the free
// function getContainerEllipse2MinCR, following the Cont* naming precedent
// (container type as the suffix so the exported name is globally unique; see
// ContAlignedBox.ts and ContOrientedBox2.ts). Upstream writes the result into
// a caller-provided 'Real D[2]'; the port returns the two diagonal entries as
// a number[] of length 2. The private static MaxProduct helper becomes a
// module-private function. An empty point set is rejected with logAssert,
// matching the Cont* empty-input precedent (upstream's LogAssert on
// 'iXMin != -1 && iYMin != -1' fires in that case, but only after the sorts).
// Upstream's 'x1 = xMax;' at the bottom of the hull-walk loop is a dead store
// (x1 is reassigned at the top of the next iteration) and has no port.
//
// One upstream bug is fixed: when the hull walk steps onto a vertical
// constraint line (one whose D[1] coefficient is 0, which happens whenever an
// input point lies on the second ellipse axis), upstream evaluates 0/0 and
// returns NaN for D[1]. See the comment at the 'r < x0' test in maxProduct.

import { logAssert } from './Logger.js';
import { type Matrix2x2 } from './Matrix2x2.js';
import { mulMatrix } from './Matrix.js';
import { Vector, compMul, sub } from './Vector.js';
import { dotPerp } from './Vector2.js';

// Compute the diagonal matrix D = diag(D[0],D[1]) of the minimum-area
// ellipse with the specified center and orientation that contains the input
// points. The returned array has length 2.
export function getContainerEllipse2MinCR(points: readonly Vector[],
    C: Vector, R: Matrix2x2): number[] {
    logAssert(points.length > 0, 'getContainerEllipse2MinCR: no points.');
    logAssert(C.size === 2, 'getContainerEllipse2MinCR: center must be 2D.');
    logAssert(R.numRows === 2 && R.numCols === 2,
        'getContainerEllipse2MinCR: rotation must be 2x2.');
    for (const point of points) {
        logAssert(point.size === 2,
            'getContainerEllipse2MinCR: points must be 2D.');
    }

    // Compute the constraint coefficients, of the form (A[0],A[1]) for each i.
    let A: Vector[] = new Array<Vector>(points.length);
    for (let i = 0; i < points.length; ++i) {
        const diff = sub(points[i], C);  // P[i] - C
        const prod = mulMatrix(diff, R);  // R^T*(P[i] - C) = (u,v)
        A[i] = compMul(prod, prod);  // (u^2, v^2)
    }

    // Use a lexicographical sort to eliminate redundant constraints. Remove
    // all but the first entry in blocks with x0 = x1 because the
    // corresponding constraint line for the first entry hides all the others
    // from the origin.
    A.sort((P0, P1) => {
        if (P0.values[0] > P1.values[0]) { return -1; }
        if (P0.values[0] < P1.values[0]) { return 1; }
        if (P0.values[1] > P1.values[1]) { return -1; }
        if (P0.values[1] < P1.values[1]) { return 1; }
        return 0;
    });
    A = uniqueBy(A, 0);

    // Use a lexicographical sort to eliminate redundant constraints. Remove
    // all but the first entry in blocks with y0 = y1 because the
    // corresponding constraint line for the first entry hides all the others
    // from the origin.
    A.sort((P0, P1) => {
        if (P0.values[1] > P1.values[1]) { return -1; }
        if (P0.values[1] < P1.values[1]) { return 1; }
        if (P0.values[0] > P1.values[0]) { return -1; }
        if (P0.values[0] < P1.values[0]) { return 1; }
        return 0;
    });
    A = uniqueBy(A, 1);

    const D = [0, 0];
    maxProduct(A, D);
    return D;
}

// The port of std::unique with a predicate that compares only component
// 'index' of consecutive (sorted) entries: keep the first entry of each block
// of equal component values.
function uniqueBy(A: readonly Vector[], index: number): Vector[] {
    const result: Vector[] = [];
    for (const a of A) {
        if (result.length === 0 ||
            result[result.length - 1].values[index] !== a.values[index]) {
            result.push(a);
        }
    }
    return result;
}

// Maximize x*y subject to x >= 0, y >= 0 and A[i][0]*x + A[i][1]*y <= 1 for
// 0 <= i < N, where A[i][0] >= 0 and A[i][1] >= 0.
function maxProduct(A: readonly Vector[], D: number[]): void {
    // Keep track of which constraint lines have already been used in the
    // search.
    const numConstraints = A.length;
    const used = new Array<boolean>(numConstraints).fill(false);

    // Find the constraint line whose y-intercept (0,ymin) is closest to the
    // origin. This line contributes to the convex hull of the constraints and
    // the search for the maximum starts here. Also find the constraint line
    // whose x-intercept (xmin,0) is closest to the origin. This line
    // contributes to the convex hull of the constraints and the search for
    // the maximum terminates before or at this line.
    let i: number;
    let iYMin = -1;
    let iXMin = -1;
    let axMax = 0, ayMax = 0;  // A[i] >= (0,0) by design
    for (i = 0; i < numConstraints; ++i) {
        // The minimum x-intercept is 1/A[iXMin][0] for A[iXMin][0] the
        // maximum of the A[i][0].
        if (A[i].values[0] > axMax) {
            axMax = A[i].values[0];
            iXMin = i;
        }

        // The minimum y-intercept is 1/A[iYMin][1] for A[iYMin][1] the
        // maximum of the A[i][1].
        if (A[i].values[1] > ayMax) {
            ayMax = A[i].values[1];
            iYMin = i;
        }
    }
    logAssert(iXMin !== -1 && iYMin !== -1, 'Unexpected condition.');
    used[iYMin] = true;

    // The convex hull is searched in a clockwise manner starting with the
    // constraint line constructed above. The next vertex of the hull occurs
    // as the closest point to the first vertex on the current constraint
    // line. The following loop finds each consecutive vertex.
    let x0 = 0;
    const xMax = 1 / axMax;
    // The line the walk was on before the current one. It is needed only for
    // the fix to upstream's division by zero described at the 'r < x0' test
    // below.
    let iPrev = -1;
    let j: number;
    for (j = 0; j < numConstraints; ++j) {
        // Find the line whose intersection with the current line is closest
        // to the last hull vertex. The last vertex is at (x0,y0) on the
        // current line.
        let x1 = xMax;
        let line = -1;
        for (i = 0; i < numConstraints; ++i) {
            if (!used[i]) {
                // This line not yet visited, process it. Given current
                // constraint line a0*x+b0*y = 1 and candidate line
                // a1*x+b1*y = 1, find the point of intersection. The
                // determinant of the system is d = a0*b1-a1*b0. We care only
                // about lines that have more negative slope than the previous
                // one, that is, -a1/b1 < -a0/b0, in which case we process
                // only lines for which d < 0.
                const det = dotPerp(A[iYMin], A[i]);
                if (det < 0) {
                    // Compute the x-value for the point of intersection,
                    // (x1,y1). There may be floating point error issues in
                    // the comparison 'D[0] <= x1'. Consider modifying to
                    // 'D[0] <= x1 + epsilon'.
                    D[0] = (A[i].values[1] - A[iYMin].values[1]) / det;
                    if (x0 < D[0] && D[0] <= x1) {
                        line = i;
                        x1 = D[0];
                    }
                }
            }
        }

        // Next vertex is at (x1,y1) whose x-value was computed above. First
        // check for the maximum of x*y on the current line for x in [x0,x1].
        // On this interval the function is f(x) = x*(1-a0*x)/b0. The
        // derivative is f'(x) = (1-2*a0*x)/b0 and f'(r) = 0 when
        // r = 1/(2*a0). The three candidates for the maximum are f(x0), f(r)
        // and f(x1). Comparisons are made between r and the endpoints x0 and
        // x1. Because a0 = 0 is possible (constraint line is horizontal and f
        // is increasing on line), the division in r is not performed and the
        // comparisons are made between 1/2 = a0*r and a0*x0 or a0*x1.

        // Compare r < x0.
        if (0.5 < A[iYMin].values[0] * x0) {
            // The maximum is f(x0) since the quadratic f decreases for x > r.
            // The value D[1] is f(x0).
            D[0] = x0;
            if (A[iYMin].values[1] !== 0) {
                D[1] = (1 - A[iYMin].values[0] * D[0]) / A[iYMin].values[1];
            } else {
                // Upstream bug: when the walk steps onto a vertical
                // constraint line (b0 = 0, which happens whenever a point is
                // on the second ellipse axis, so that its constraint is
                // a0*x <= 1), this test is entered immediately with
                // a0*x0 = 1 > 1/2 and upstream evaluates 0/0, returning NaN
                // for D[1]. The vertex (x0,y0) is shared by the current line
                // and the previous line of the walk, so evaluate y0 on the
                // previous line, which has b != 0. (The first line of the
                // walk has b > 0 by construction, and this test cannot be
                // entered on the first iteration because x0 = 0 there, so
                // iPrev is a valid index with a nonzero b whenever the walk
                // reaches a vertical line.)
                D[1] = (1 - A[iPrev].values[0] * D[0]) / A[iPrev].values[1];
            }
            break;
        }

        // Compare r < x1.
        if (0.5 < A[iYMin].values[0] * x1) {
            // The maximum is f(r). The search ends here because the current
            // line is tangent to the level curve of f(x) = f(r) and x*y can
            // therefore only decrease as we traverse farther around the hull
            // in the clockwise direction. The value D[1] is f(r).
            D[0] = 0.5 / A[iYMin].values[0];
            D[1] = 0.5 / A[iYMin].values[1];
            break;
        }

        // The maximum is f(x1). The function x*y is potentially larger on the
        // next line, so continue the search.
        logAssert(line !== -1, 'Unexpected condition.');
        x0 = x1;
        used[line] = true;
        iPrev = iYMin;
        iYMin = line;
    }

    logAssert(j < numConstraints, 'Unexpected condition.');
}
