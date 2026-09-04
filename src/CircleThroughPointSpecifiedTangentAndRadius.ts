// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) CircleThroughPointSpecifiedTangentAndRadius.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// This file provides an implementation of the algorithm in Section 8.7 of the
// book
//    Geometric Tools for Computer Graphics,
//    Philip J. Schneider and David H. Eberly,
//    Morgan Kaufmann, San Francisco CA, 2002
//
// Given a point P, a radius r and a line Dot(N,X-A) = 0, where A is a point
// on the line and N is a unit-length normal to the line. Compute the centers
// of circles, each containing the point, having the specified radius and
// having the line as a tangent line. The book describes one algebraic
// approach to solving the problem. The implementation here is another
// approach, a portion using the algorithm of Section 8.6
// (CircleThroughTwoPointsSpecifiedRadius).
//
// Let N = (n0,n1) and define the unit-length perpendicular D = Perp(N) =
// (n1,-n0). (The upstream comment says (-n1,n0), but GTE's Perp is (n1,-n0);
// the code uses Perp, so the port's comment is corrected here. Either
// perpendicular works, the choice only swaps circle[0] and circle[1].)
// Represent P = A+u*D+s*N with parameters u = Dot(D,P-A) and
// s = Dot(N,P-A). The parameter s is the signed distance from P to the line.
// To simplify the logic of the implementation, if s < 0, the values of s and
// N are negated (and therefore D is derived from the negated N). The
// discussion below assumes s >= 0.
//
// The cases are
//
//   (1) s = 0: P is on the line. There are two circles containing P and
//       tangent to the line at P. The circle centers are C0 = P-r*N and
//       C1 = P+r*N.
//
//   (2) s = r: The two circles have a single point of intersection, which is
//       P. The circle centers are C0 = P-r*D and C1 = P+r*D.
//
//   (3) s = 2*r: P is the farthest point on a circle of radius r which has
//       the line as the tangent line. The circle center is C0 = P-r*N.
//
//   (4) s > 2*r: The distance from P to the tangent line is larger than the
//       desired circle diameter, so there is no circle that satisfies the
//       constraints.
//
//   (5a) 0 < s < r: The two circles intersect in P. The other point of
//        intersection is Q = A+u*D+(2*r-s)*N. The bisector of segment <P,Q>
//        has origin B = (P+Q)/2 = A+u*D+r*N and direction D. If a circle
//        center is C, the triangle <P,B,C> is a right triangle at B. Using
//        the Pythagorean theorem, the length of segment <B,C> is
//        h = |B-C| = sqrt(r^2 - (r-s)^2). The circle centers are C0 = B-h*D
//        and C1 = B+h*D.
//
//   (5b) r < s < 2*r: This is analogous to (5a) with the roles of the two
//        intersection points swapped; Q is the intersection point closest to
//        the tangent line and P is the one farthest from it. The construction
//        of the centers is the same as that of (5a).
//
// Port notes: the upstream output parameter 'std::array<Circle2<T>, 2>&
// circle' becomes the 'circle' field of the returned object, and the upstream
// return value (the number of circles) becomes the 'numCircles' field. The
// 'circle' array always has two elements; those with index >= numCircles have
// a zero center and a zero radius, exactly as upstream. Upstream takes the
// normal N by value and negates it in place when s < 0; the port copies N
// before negating so the caller's vector is not modified.

import { Hypersphere } from './Hypersphere.js';
import { logAssert } from './Logger.js';
import { Vector, add, dot, mul, negate, sub } from './Vector.js';
import { perp } from './Vector2.js';

export interface CircleThroughPointSpecifiedTangentAndRadiusResult {
    // The number of circles satisfying the constraints.
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

export function circleThroughPointSpecifiedTangentAndRadius(P: Vector,
    A: Vector, N: Vector, r: number): CircleThroughPointSpecifiedTangentAndRadiusResult {
    logAssert(P.size === 2 && A.size === 2 && N.size === 2,
        'CircleThroughPointSpecifiedTangentAndRadius: the inputs must have size 2.');

    const zero = 0;

    // Upstream passes N by value; the port copies so the caller's vector is
    // not modified by the negation below.
    let normal = N.clone();

    const PmA = sub(P, A);
    let s = dot(normal, PmA);
    if (s === zero) {
        // Case (1).
        return {
            numCircles: 2,
            circle: [
                Hypersphere.fromCenterRadius(sub(P, mul(normal, r)), r),
                Hypersphere.fromCenterRadius(add(P, mul(normal, r)), r)
            ]
        };
    }

    if (s < zero) {
        normal = negate(normal);
        s = -s;
    }

    if (s === r) {
        // Case (2).
        const D = perp(normal);
        return {
            numCircles: 2,
            circle: [
                Hypersphere.fromCenterRadius(sub(P, mul(D, r)), r),
                Hypersphere.fromCenterRadius(add(P, mul(D, r)), r)
            ]
        };
    }

    const twoR = 2 * r;
    if (s === twoR) {
        // Case (3).
        return {
            numCircles: 1,
            circle: [
                Hypersphere.fromCenterRadius(sub(P, mul(normal, r)), r),
                zeroCircle()
            ]
        };
    }

    if (s > twoR) {
        // Case (4).
        return {
            numCircles: 0,
            circle: [zeroCircle(), zeroCircle()]
        };
    }

    // The bisector direction is D = Perp(N) and the bisector origin is
    // B = (P + Q) / 2 = A + t * D + r * N with t = Dot(D, P - A).
    const bisectorDirection = perp(normal);
    const t = dot(bisectorDirection, PmA);
    const bisectorOrigin = add(add(A, mul(bisectorDirection, t)), mul(normal, r));

    const diffRS = r - s;
    const argument = r * r - diffRS * diffRS;
    if (argument > zero) {
        const h = Math.sqrt(argument);
        return {
            numCircles: 2,
            circle: [
                Hypersphere.fromCenterRadius(sub(bisectorOrigin, mul(bisectorDirection, h)), r),
                Hypersphere.fromCenterRadius(add(bisectorOrigin, mul(bisectorDirection, h)), r)
            ]
        };
    }
    else {
        // Theoretically this code cannot be reached, but floating-point
        // rounding errors might trigger it. This corresponds to Case (3)
        // where r = s.
        return {
            numCircles: 1,
            circle: [
                Hypersphere.fromCenterRadius(bisectorOrigin, r),
                zeroCircle()
            ]
        };
    }
}
