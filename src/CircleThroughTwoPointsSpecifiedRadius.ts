// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) CircleThroughTwoPointsSpecifiedRadius.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// This file provides an implementation of the algorithm in Section 8.6 of the
// book
//    Geometric Tools for Computer Graphics,
//    Philip J. Schneider and David H. Eberly,
//    Morgan Kaufmann, San Francisco CA, 2002
//
// Given two distinct points P and Q and given a radius r, compute the centers
// of circles, each containing the points and having the specified radius.
//
// The book states that the circle centers are the points of intersection of
// circles |X-P|^2 = r^2 and |X-Q|^2 = r^2. The pseudocode simply calls a
// function to compute these intersections.
//
// A simpler approach uses the fact that the bisector of the line segment with
// endpoints P and Q is a line that contains the centers. The bisector is
// parameterized by X(t) = t*Perp(P-Q)+(P+Q)/2, where Perp(P-Q) is
// perpendicular to P-Q and has the same length as that of P-Q. We need values
// of t for which X(t)-P has length r,
//   X(t)-P = t*Perp(P-Q)-(P-Q)/2
//   r^2 = |X(t)-P|^2
//       = |t*Perp(P-Q)-(P-Q)/2|^2
//       = |Perp(P-Q)|^2 * t^2 - 2*t*Dot(Perp(P-Q),P-Q) + |P-Q|^2/4
//       = |P-Q|^2 * t^2 + |P-Q|^2/4
//       = |P-Q|^2 * (t^2 + 1/4)
// Observe that t^2+1/4 >= 1/4, which implies that r >= |P-Q|/2. This
// condition is clear geometrically. The radius must be at least half the
// length of the segment connecting P and Q.
//
// If r = |P-Q|/2, there is a single circle with center (P+Q)/2. If
// r > |P-Q|/2, there are two circles whose centers occur when
// t^2 = r^2/|P-Q|^2 - 1/4, which implies t = +/- sqrt(r^2/|P-Q|^2-1/4).
//
// Port notes: the upstream output parameter 'std::array<Circle2<T>, 2>&
// circle' becomes the 'circle' field of the returned object, and the upstream
// return value (the number of circles) becomes the 'numCircles' field. The
// 'circle' array always has two elements; those with index >= numCircles have
// a zero center and a zero radius, exactly as upstream.

import { Hypersphere } from './Hypersphere.js';
import { logAssert } from './Logger.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import { perp } from './Vector2.js';

export interface CircleThroughTwoPointsSpecifiedRadiusResult {
    // The number of circles satisfying the constraints. The number is 2 when
    // r > |P-Q|/2, 1 when r = |P-Q|/2, or 0 when P = Q or r < |P-Q|/2.
    numCircles: number;

    // The two candidate circles. Any circle[i] with i >= numCircles has
    // members set to zero.
    circle: Hypersphere[];
}

// A circle with a zero center and a zero radius, the port of the upstream
// assignments 'circle[i].center = { zero, zero }; circle[i].radius = zero;'.
function zeroCircle(): Hypersphere {
    const circle = new Hypersphere(2);
    circle.radius = 0;
    return circle;
}

function circleFromCenter(center: Vector, radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(center, radius);
}

export function circleThroughTwoPointsSpecifiedRadius(P: Vector, Q: Vector,
    r: number): CircleThroughTwoPointsSpecifiedRadiusResult {
    logAssert(P.size === 2 && Q.size === 2,
        'CircleThroughTwoPointsSpecifiedRadius: the points must have size 2.');

    const zero = 0;

    const PmQ = sub(P, Q);
    const sqrLengthPmQ = dot(PmQ, PmQ);
    if (sqrLengthPmQ !== zero) {
        const argument = r * r / sqrLengthPmQ - 0.25;
        if (argument > zero) {
            const root = Math.sqrt(argument);
            const bisectorOrigin = mul(add(P, Q), 0.5);
            const bisectorDirection = perp(PmQ);
            return {
                numCircles: 2,
                circle: [
                    circleFromCenter(sub(bisectorOrigin, mul(bisectorDirection, root)), r),
                    circleFromCenter(add(bisectorOrigin, mul(bisectorDirection, root)), r)
                ]
            };
        }

        if (argument === zero) {
            return {
                numCircles: 1,
                circle: [
                    circleFromCenter(mul(add(P, Q), 0.5), r),
                    zeroCircle()
                ]
            };
        }
    }

    return {
        numCircles: 0,
        circle: [zeroCircle(), zeroCircle()]
    };
}
