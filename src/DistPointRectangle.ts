// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPointRectangle.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a point and a rectangle in nD.
//
// The rectangle has center C, unit-length axis directions W[0] and W[1], and
// extents e[0] and e[1]. A rectangle point is X = C + sum_{i=0}^1 s[i] * W[i]
// where |s[i]| <= e[i] for all i.
//
// The input point is stored in closest[0]. The closest point on the
// rectangle is stored in closest[1] with W-coordinates (s[0],s[1]). When
// there are infinitely many choices for the pair of closest points, only one
// of them is returned.
//
// Upstream TODO: modify to support non-unit-length W[].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector<N,T>, Rectangle<N,T>>' becomes
// the class DistPointRectangle with the result type
// DistPointRectangleResult.

import type { DCPQuery } from './DCPQuery';
import type { Rectangle } from './Rectangle';
import { Vector, add, dot, mul, sub } from './Vector';

export interface DistPointRectangleResult {
    distance: number;
    sqrDistance: number;

    // The W-coordinates (s[0],s[1]) of the closest rectangle point.
    cartesian: [number, number];

    // closest[0] is the input point, closest[1] is the closest rectangle
    // point.
    closest: [Vector, Vector];
}

export class DistPointRectangle
    implements DCPQuery<Vector, Rectangle, DistPointRectangleResult> {
    compute(point: Vector, rectangle: Rectangle): DistPointRectangleResult {
        let diff = sub(point, rectangle.center);
        const closest0 = point.clone();
        let closest1 = rectangle.center.clone();
        const cartesian: [number, number] = [0, 0];
        for (let i = 0; i < 2; ++i) {
            cartesian[i] = dot(rectangle.axis[i], diff);
            if (cartesian[i] < -rectangle.extent.values[i]) {
                cartesian[i] = -rectangle.extent.values[i];
            }
            else if (cartesian[i] > rectangle.extent.values[i]) {
                cartesian[i] = rectangle.extent.values[i];
            }
            closest1 = add(closest1, mul(cartesian[i], rectangle.axis[i]));
        }

        diff = sub(closest0, closest1);
        const sqrDistance = dot(diff, diff);

        return {
            distance: Math.sqrt(sqrDistance),
            sqrDistance,
            cartesian,
            closest: [closest0, closest1]
        };
    }
}
