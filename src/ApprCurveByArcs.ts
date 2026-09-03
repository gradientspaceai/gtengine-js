// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprCurveByArcs.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Approximate a 2-dimensional parametric curve X(t) for t in [tmin,tmax] by
// a collection of circular arcs. Some of the arcs can be degenerate in that
// the arc center is a point at infinity. In this case, the arc represents a
// line segment connecting its endpoints, and the arc radius is set to
// Number.MAX_VALUE (the port of std::numeric_limits<T>::max()) to let the
// caller know the object is actually a line segment. The algorithm is
// described in
//   https://www.geometrictools.com/Documentation/ApproximateCurveByArcs.pdf
// The collection of arcs form a C0-continuous curve. Generally, the
// derivatives at a curve point shared by two arcs are not equal.
//
// Port notes: upstream takes a std::shared_ptr<ParametricCurve<2,T>> and
// three output vectors that it resizes. The port takes the curve object
// directly (the port's ParametricCurve carries its dimension at runtime,
// checked here) and returns the three containers in an object literal, per
// PORTING.md. The 'static_assert(std::is_floating_point<T>)' has no port
// (every numeric type is 'number').

import { Arc2 } from './Arc2.js';
import { logAssert } from './Logger.js';
import { ParametricCurve } from './ParametricCurve.js';
import { Vector, add, dot, length, mul, sub } from './Vector.js';
import { dotPerp } from './Vector2.js';

export interface ApproximateCurveByArcsResult {
    // The 2*numArcs+1 curve parameters. The even-indexed entries are the
    // arc-length subdivision times; the odd-indexed entries are the times of
    // the arc midpoints.
    times: number[];

    // The curve points X(times[i]).
    points: Vector[];

    // The numArcs arcs approximating the curve.
    arcs: Arc2[];
}

// The number of arcs N (input numArcs) must be positive. The returned
// times[] and points[] have 2*N+1 elements. The parametric curve times and
// samples are stored in the even-indexed locations of these containers. The
// odd-indexed locations store the times and midpoints that are used for
// fitting arcs to subcurves.
//
// An arc has endpoints {P0,P1} = {points[2*i],points[2*i+2]} corresponding
// to parameters {t0,t1} = {times[2*i],times[2*i+2]}. The midpoint of the arc
// is at M = points[2*i+1] corresponding to parameter tmid = times[2*i+1].
//
// The arc has a center point C and a radius r. If {P0,M,P1} are not
// colinear, the radius is finite and the arc is truly an arc. If the point
// triple is colinear (or nearly colinear), the center components and radius
// are set to Number.MAX_VALUE to let the caller know that the arc represents
// a line segment; the segment endpoints are the arc endpoints.
//
// The arc center is computed as the solution to a linear system of equations
// in the components of C. If the determinant of this system is nearly 0, the
// triple {P0,M,P1} is nearly collinear. The epsilon input to the function
// must be nonnegative and is a lower threshold for the determinant; that is,
// if the determinant is greater or equal to epsilon, the linear system is
// solved for C. If the determinant is smaller than epsilon, the arc center
// and radius are set as mentioned in the previous paragraph to let the
// caller know the arc really represents a line segment. You can set epsilon
// to 0, but nearly colinear {P0,M,P1} can cause floating-point rounding
// errors to produce inaccurate center and radius.
export function approximateCurveByArcs(curve: ParametricCurve,
    numArcs: number, epsilon: number = 0): ApproximateCurveByArcsResult {
    logAssert(curve !== null && curve !== undefined && numArcs >= 1,
        'Invalid input.');
    logAssert(curve.getDimension() === 2, 'The curve must be 2-dimensional.');

    const numTimes = 2 * numArcs + 1;
    const times = new Array<number>(numTimes).fill(0);
    const points = new Array<Vector>(numTimes);
    const arcs = new Array<Arc2>(numArcs);
    for (let i = 0; i < numTimes; ++i) {
        points[i] = new Vector(2);
    }

    // Subdivide the curve by arc length. The arc length between any pair of
    // consecutive points is constant. The consecutive points are stored in
    // the even-indexed locations. The odd-indexed locations are assigned in
    // the block of code after this one.
    const totalLength = curve.getTotalLength();
    const deltaLength = totalLength / (numTimes - 1);
    for (let i = 0; i < numTimes; i += 2) {
        const arcLength = deltaLength * i;
        times[i] = curve.getTime(arcLength);
        points[i] = curve.getPosition(times[i]);
    }

    const half = 0.5;
    for (let i = 0, j0 = 0, j1 = 1, j2 = 2; i < numArcs;
        ++i, j0 += 2, j1 += 2, j2 += 2) {
        const arc = new Arc2();
        arcs[i] = arc;
        arc.end[0] = points[j0].clone();
        arc.end[1] = points[j2].clone();
        const P0 = arc.end[0];
        const P1 = arc.end[1];

        // Let P0 = arc.end[0] and P1 = arc.end[1]. Compute a point of
        // intersection between the bisector of segment <P0, P1> and the
        // curve X(t). This is accomplished using bisection for
        // F(t) = Dot(D, X(t) - A) on [t0,t1] with P0 = X(t0), P1 = X(t1),
        // D = P1 - P0, and A = (P0 + P1) / 2. Observe that
        //   F(t0) = Dot(D, P0 - A) = -|D|^2/2 < 0
        //   F(t1) = Dot(D, P1 - A) = +|D|^2/2 > 0
        // There must be a tRoot in [t0,t1] for which F(tRoot) = 0.
        //
        // The loop is guaranteed to terminate because a sufficient number of
        // iterations will either find a tRoot where F(tRoot) = 0 using
        // floating-point computations, or the interval to bisect has
        // consecutive floating-point endpoints and the interval midpoint
        // rounds to one of those endpoints.
        const D = sub(P1, P0);
        const A = mul(half, add(P0, P1));
        let t0 = times[j0], t1 = times[j2], tRoot = 0;
        for (; ;) {
            tRoot = half * (t0 + t1);
            const fAtRoot = dot(D, sub(curve.getPosition(tRoot), A));
            const signRoot = (fAtRoot > 0 ? +1 : (fAtRoot < 0 ? -1 : 0));
            if (signRoot === 0 || tRoot === t0 || tRoot === t1) {
                break;
            }

            if (signRoot === -1) {
                t0 = tRoot;
            }
            else {  // signRoot = +1
                t1 = tRoot;
            }
        }

        // Fill in the odd-indexed values.
        times[j1] = tRoot;
        points[j1] = curve.getPosition(tRoot);
        const M = points[j1];

        // The points P0, X(tRoot), and P1 are circumscribed to determine the
        // arc. If the three points are colinear, the center and radius of
        // the arc are set to Number.MAX_VALUE as a signal to the caller that
        // the arc represents a line segment.
        const diffP0M = sub(P0, M);
        const diffP1M = sub(P1, M);
        const avrgP0M = mul(half, add(P0, M));
        const avrgP1M = mul(half, add(P1, M));
        const dot0 = dot(diffP0M, avrgP0M);
        const dot1 = dot(diffP1M, avrgP1M);
        const det = dotPerp(diffP0M, diffP1M);
        if (Math.abs(det) >= epsilon) {
            arc.center.values[0] =
                (diffP1M.values[1] * dot0 - diffP0M.values[1] * dot1) / det;
            arc.center.values[1] =
                (diffP0M.values[0] * dot1 - diffP1M.values[0] * dot0) / det;
            arc.radius = length(sub(M, arc.center));
        }
        else {
            const tmax = Number.MAX_VALUE;
            arc.center.values[0] = tmax;
            arc.center.values[1] = tmax;
            arc.radius = tmax;
        }
    }

    return { times: times, points: points, arcs: arcs };
}
