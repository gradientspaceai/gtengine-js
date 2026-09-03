// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPoint2Circle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a point and a circle in 2D. The circle is
// considered to be a curve, not a solid disk.
//
// The input point is stored in the member closest[0]. If a single point on
// the circle is closest to the input point, the member closest[1] is set to
// that point and the equidistant member is set to false. If the entire circle
// is equidistant to the point, the member closest[1] is set to C+r*(1,0),
// where C is the circle center and r is the circle radius. Moreover, the
// equidistant member is set to true.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector2<T>, Circle2<T>>' becomes the
// class DistPoint2Circle2 with the result type DistPoint2Circle2Result.

import type { DCPQuery } from './DCPQuery.js';
import type { Circle2 } from './Hypersphere.js';
import { Vector, add, div, dot, mul, sub } from './Vector.js';

export interface DistPoint2Circle2Result {
    distance: number;
    sqrDistance: number;

    // closest[0] is the input point, closest[1] is the closest circle point.
    closest: [Vector, Vector];

    // True when the input point is the circle center, in which case every
    // circle point is equidistant from it.
    equidistant: boolean;
}

export class DistPoint2Circle2
    implements DCPQuery<Vector, Circle2, DistPoint2Circle2Result> {
    compute(point: Vector, circle: Circle2): DistPoint2Circle2Result {
        let diff = sub(point, circle.center);
        const sqrLength = dot(diff, diff);
        const length = Math.sqrt(sqrLength);
        if (length > 0) {
            diff = div(diff, length);
            const distance = Math.abs(length - circle.radius);
            return {
                distance,
                sqrDistance: distance * distance,
                closest: [
                    point.clone(),
                    add(circle.center, mul(circle.radius, diff))
                ],
                equidistant: false
            };
        }
        else {
            return {
                distance: circle.radius,
                sqrDistance: circle.radius * circle.radius,
                closest: [
                    point.clone(),
                    add(circle.center,
                        mul(circle.radius, Vector.unit(2, 0)))
                ],
                equidistant: true
            };
        }
    }
}
