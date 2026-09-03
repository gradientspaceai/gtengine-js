// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRay2Circle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a ray and a circle in 2D. The circle is
// considered to be a curve, not a solid disk.
//
// The ray is P + t * D, where P is a point on the ray and D is not required
// to be unit length. The t-value satisfies t >= 0.
//
// The circle is C + r * U(s), where C is the center, r > 0 is the radius, and
// U(s) = (cos(s), sin(s)) for s in [0,2*pi).
//
// The number of pairs of closest points is result.numClosestPairs which is
// 1 or 2. If result.numClosestPairs is 1, result.parameter[0] is the ray
// t-value for its closest point result.closest[0][0]. The circle closest
// point is result.closest[0][1]. If result.numClosestPairs is 2,
// result.parameter[0] and result.parameter[1] are the ray t-values for its
// closest points result.closest[0][0] and result.closest[1][0]. The circle
// closest points are result.closest[0][1] and result.closest[1][1].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Ray2<T>, Circle2<T>>' becomes the
// class DistRay2Circle2 with the result type DistRay2Circle2Result, which is
// structurally identical to the line-circle result that upstream aliases. The
// two private static 'Update' overloads become the module-private functions
// updateForRay and updateForOrigin.

import type { DCPQuery } from './DCPQuery.js';
import { DistLine2Circle2 } from './DistLine2Circle2.js';
import { DistPoint2Circle2 } from './DistPoint2Circle2.js';
import type { Circle2 } from './Hypersphere.js';
import { Line } from './Line.js';
import type { Ray2 } from './Ray.js';
import { Vector } from './Vector.js';

export interface DistRay2Circle2Result {
    distance: number;
    sqrDistance: number;

    // The number of pairs of closest points, 1 or 2.
    numClosestPairs: number;

    // parameter[j] is the ray t-value of closest[j][0].
    parameter: [number, number];

    // closest[j][0] is on the ray, closest[j][1] is on the circle.
    closest: [[Vector, Vector], [Vector, Vector]];
}

// Compute the closest circle point to the ray origin and overwrite the
// line-circle result with it. The ray does not intersect the circle even
// though the line does.
function updateForOrigin(origin: Vector, circle: Circle2,
    lcResult: DistRay2Circle2Result): void {
    const pcResult = new DistPoint2Circle2().compute(origin, circle);

    lcResult.distance = pcResult.distance;
    lcResult.sqrDistance = pcResult.sqrDistance;
    lcResult.numClosestPairs = 1;
    lcResult.parameter[0] = 0;
    lcResult.parameter[1] = 0;
    lcResult.closest[0][0] = pcResult.closest[0];
    lcResult.closest[0][1] = pcResult.closest[1];
    lcResult.closest[1][0] = new Vector(2);
    lcResult.closest[1][1] = new Vector(2);
}

function updateForRay(ray: Ray2, circle: Circle2,
    lcResult: DistRay2Circle2Result): void {
    const t0 = lcResult.parameter[0];
    const t1 = lcResult.parameter[1];

    if (t1 <= 0) {
        // The ray origin is the closest point to the circle.
        updateForOrigin(ray.origin, circle, lcResult);
    }
    else if (t0 < 0) {
        // The ray origin is strictly inside the circle. Remove the t0-point.
        lcResult.numClosestPairs = 1;
        lcResult.parameter[0] = lcResult.parameter[1];
        lcResult.parameter[1] = 0;
        lcResult.closest[0][0] = lcResult.closest[1][0];
        lcResult.closest[0][1] = lcResult.closest[1][1];
        lcResult.closest[1][0] = new Vector(2);
        lcResult.closest[1][1] = new Vector(2);
    }
    // Otherwise 0 <= t0 < t1 and the line-circle intersection points are
    // contained by the ray.
}

export class DistRay2Circle2
    implements DCPQuery<Ray2, Circle2, DistRay2Circle2Result> {
    compute(ray: Ray2, circle: Circle2): DistRay2Circle2Result {
        // Execute the query for line-circle.
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const lcResult: DistRay2Circle2Result =
            new DistLine2Circle2().compute(line, circle);

        // Restrict the analysis to ray-circle.
        if (lcResult.numClosestPairs === 2) {
            // The segment connecting the line-circle intersection points has
            // parameter interval [t0,t1]. Determine how this intersects with
            // the ray interval [0,+infinity) and modify lcResult accordingly.
            updateForRay(ray, circle, lcResult);
        }
        else {
            // The line does not intersect the circle or is tangent to the
            // circle. If the closest line point to the circle has a negative
            // parameter, then the ray is outside the circle and the ray
            // origin is the closest ray point to the circle.
            if (lcResult.parameter[0] < 0) {
                updateForOrigin(ray.origin, circle, lcResult);
            }
        }

        return lcResult;
    }
}
