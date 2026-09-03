// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrPlane3Circle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the circle to be a 1-dimensional object (the circular
// curve, not the disk).
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has both
// a TIQuery and an FIQuery specialization, which become IntrPlane3Circle3TI
// and IntrPlane3Circle3FI. 'std::numeric_limits<size_t>::max()' becomes
// Number.MAX_SAFE_INTEGER (the SIZE_MAX precedent of BVTree.ts), exported as
// 'intrPlane3Circle3InfinitePoints'.

import { Circle3 } from './Circle3.js';
import type { FIQuery } from './FIQuery.js';
import { Hyperplane } from './Hyperplane.js';
import type { Plane3 } from './Hyperplane.js';
import { IntrPlane3Plane3FI } from './IntrPlane3Plane3.js';
import { logAssert } from './Logger.js';
import { Vector, dot, mul, sub } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The 'numIntersections' value that reports the entire circle as the set of
// intersection (the port of std::numeric_limits<size_t>::max()).
export const intrPlane3Circle3InfinitePoints = Number.MAX_SAFE_INTEGER;

// The result of IntrPlane3Circle3TI.test.
export interface IntrPlane3Circle3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrPlane3Circle3TIResult():
    IntrPlane3Circle3TIResult {
    return { intersect: false };
}

// The result of IntrPlane3Circle3FI.find.
//
// If 'intersect' is false, the set of intersection is empty.
// 'numIntersections' is 0 and 'point[]' and 'circle' have members all set
// to 0.
//
// If 'intersect' is true, the set of intersection contains either 1 or 2
// points or the entire circle.
//
// (1) When the set of intersection has 1 point, the circle is just touching
//     the plane. 'numIntersections' is 1 and 'point[0]' and 'point[1]' are
//     the same point. The 'circle' is set to invalid (center at the origin,
//     normal is the zero vector, radius is 0).
//
// (2) When the set of intersection has 2 points, the plane cuts the circle
//     into 2 arcs. 'numIntersections' is 2 and 'point[0]' and 'point[1]' are
//     the distinct intersection points. The 'circle' is set to invalid.
//
// (3) When the set of intersection contains the entire circle, the plane of
//     the circle and the input plane are the same. 'numIntersections' is
//     intrPlane3Circle3InfinitePoints. 'point[0]' and 'point[1]' are set to
//     the zero vector. 'circle' is set to the input circle.
export interface IntrPlane3Circle3FIResult {
    intersect: boolean;
    numIntersections: number;
    point: [Vector, Vector];
    circle: Circle3;
}

// The port of the upstream FIQuery::Result default constructor. Upstream
// builds the circle from an explicit (zero center, zero normal, zero radius);
// that is not the Circle3 default constructor, which has a unit normal and
// radius 1.
export function defaultIntrPlane3Circle3FIResult():
    IntrPlane3Circle3FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        point: [Vector.zero(3), Vector.zero(3)],
        circle: Circle3.fromCenterNormalRadius(Vector.zero(3),
            Vector.zero(3), 0)
    };
}

// Test-intersection query for a plane and a circle in 3D.
export class IntrPlane3Circle3TI implements
    TIQuery<Plane3, Circle3, IntrPlane3Circle3TIResult> {

    test(plane: Plane3, circle: Circle3): IntrPlane3Circle3TIResult {
        logAssert(plane.dimension === 3, 'IntrPlane3Circle3TI: mismatched sizes.');
        const result = defaultIntrPlane3Circle3TIResult();

        // Construct the plane of the circle.
        const cPlane = Hyperplane.fromNormalOrigin(circle.normal, circle.center);

        // Compute the intersection of this plane with the input plane.
        const ppQuery = new IntrPlane3Plane3FI();
        const ppResult = ppQuery.find(plane, cPlane);
        if (!ppResult.intersect) {
            // The planes are parallel and nonintersecting.
            result.intersect = false;
            return result;
        }

        if (!ppResult.isLine) {
            // The planes are the same, so the circle is the set of
            // intersection.
            result.intersect = true;
            return result;
        }

        // The planes intersect in a line. Locate one or two points that are
        // on the circle and line. If the line is t*D+P, the circle center is
        // C and the circle radius is r, then
        //   r^2 = |t*D+P-C|^2 = |D|^2*t^2 + 2*Dot(D,P-C)*t + |P-C|^2
        // This is a quadratic equation of the form
        // a2*t^2 + 2*a1*t + a0 = 0.
        const diff = sub(ppResult.line.origin, circle.center);
        const a2 = dot(ppResult.line.direction, ppResult.line.direction);
        const a1 = dot(diff, ppResult.line.direction);
        const a0 = dot(diff, diff) - circle.radius * circle.radius;

        // Real-valued roots imply an intersection.
        const discr = a1 * a1 - a0 * a2;
        result.intersect = (discr >= 0);
        return result;
    }
}

// Find-intersection query for a plane and a circle in 3D.
export class IntrPlane3Circle3FI implements
    FIQuery<Plane3, Circle3, IntrPlane3Circle3FIResult> {

    find(plane: Plane3, circle: Circle3): IntrPlane3Circle3FIResult {
        logAssert(plane.dimension === 3, 'IntrPlane3Circle3FI: mismatched sizes.');

        // The 'result' members have initial values set by the default
        // factory. In each return block, only the relevant members are
        // modified.
        const result = defaultIntrPlane3Circle3FIResult();

        // Construct the plane of the circle.
        const cPlane = Hyperplane.fromNormalOrigin(circle.normal, circle.center);

        // Compute the intersection of this plane with the input plane.
        const ppQuery = new IntrPlane3Plane3FI();
        const ppResult = ppQuery.find(plane, cPlane);
        if (!ppResult.intersect) {
            // The planes are parallel and nonintersecting.
            return result;
        }

        if (!ppResult.isLine) {
            // The planes are the same, so the circle is the set of
            // intersection.
            result.intersect = true;
            result.numIntersections = intrPlane3Circle3InfinitePoints;
            result.circle = circle.clone();
            return result;
        }

        // The planes intersect in a line. Locate one or two points that are
        // on the circle and line. If the line is t*D+P, the circle center is
        // C, and the circle radius is r, then
        //   r^2 = |t*D+P-C|^2 = |D|^2*t^2 + 2*Dot(D,P-C)*t + |P-C|^2
        // This is a quadratic equation of the form
        // a2*t^2 + 2*a1*t + a0 = 0.
        const diff = sub(ppResult.line.origin, circle.center);
        const a2 = dot(ppResult.line.direction, ppResult.line.direction);
        const a1 = dot(diff, ppResult.line.direction);
        const a0 = dot(diff, diff) - circle.radius * circle.radius;

        const discr = a1 * a1 - a0 * a2;
        if (discr < 0) {
            // No real roots, the circle does not intersect the plane.
            return result;
        }

        if (discr === 0) {
            // The quadratic polynomial has 1 real-valued repeated root. The
            // circle just touches the plane.
            result.intersect = true;
            result.numIntersections = 1;
            result.point[0] = sub(ppResult.line.origin,
                mul(a1 / a2, ppResult.line.direction));
            result.point[1] = result.point[0].clone();
            return result;
        }

        // The quadratic polynomial has 2 distinct, real-valued roots. The
        // circle intersects the plane in two points.
        const root = Math.sqrt(discr);
        result.intersect = true;
        result.numIntersections = 2;
        result.point[0] = sub(ppResult.line.origin,
            mul((a1 + root) / a2, ppResult.line.direction));
        result.point[1] = sub(ppResult.line.origin,
            mul((a1 - root) / a2, ppResult.line.direction));
        return result;
    }
}
