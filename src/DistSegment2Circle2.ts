// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistSegment2Circle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a segment and a circle in 2D. The circle is
// considered to be a curve, not a solid disk.
//
// The segment has endpoints P0 and P1 and is parameterized by
// P0 + t * (P1 - P0). The t-value satisfies 0 <= t <= 1.
//
// The circle is C + r * U(s), where C is the center, r > 0 is the radius, and
// U(s) = (cos(s), sin(s)) for s in [0,2*pi).
//
// The number of pairs of closest points is result.numClosestPairs which is
// 1 or 2. If result.numClosestPairs is 1, result.parameter[0] is the segment
// t-value for its closest point result.closest[0][0]. The circle closest
// point is result.closest[0][1]. If result.numClosestPairs is 2,
// result.parameter[0] and result.parameter[1] are the segment t-values for
// its closest points result.closest[0][0] and result.closest[1][0]. The
// circle closest points are result.closest[0][1] and result.closest[1][1].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Segment2<T>, Circle2<T>>' becomes the
// class DistSegment2Circle2 with the result type DistSegment2Circle2Result,
// which is structurally identical to the line-circle result that upstream
// aliases. The two private static 'Update' overloads become the
// module-private functions updateForSegment and updateForEndpoint.

import type { DCPQuery } from './DCPQuery.js';
import { DistLine2Circle2 } from './DistLine2Circle2.js';
import { DistPoint2Circle2 } from './DistPoint2Circle2.js';
import type { Circle2 } from './Hypersphere.js';
import { Line } from './Line.js';
import type { Segment2 } from './Segment.js';
import { Vector, sub } from './Vector.js';

export interface DistSegment2Circle2Result {
    distance: number;
    sqrDistance: number;

    // The number of pairs of closest points, 1 or 2.
    numClosestPairs: number;

    // parameter[j] is the segment t-value of closest[j][0].
    parameter: [number, number];

    // closest[j][0] is on the segment, closest[j][1] is on the circle.
    closest: [[Vector, Vector], [Vector, Vector]];
}

// Compute the closest circle point to a segment endpoint and overwrite the
// line-circle result with it. The segment does not intersect the circle even
// though the line does.
function updateForEndpoint(endpoint: Vector, parameter: number,
    circle: Circle2, lcResult: DistSegment2Circle2Result): void {
    const pcResult = new DistPoint2Circle2().compute(endpoint, circle);

    lcResult.distance = pcResult.distance;
    lcResult.sqrDistance = pcResult.sqrDistance;
    lcResult.numClosestPairs = 1;
    lcResult.parameter[0] = parameter;
    lcResult.parameter[1] = 0;
    lcResult.closest[0][0] = pcResult.closest[0];
    lcResult.closest[0][1] = pcResult.closest[1];
    lcResult.closest[1][0] = new Vector(2);
    lcResult.closest[1][1] = new Vector(2);
}

function updateForSegment(segment: Segment2, circle: Circle2,
    lcResult: DistSegment2Circle2Result): void {
    const t0 = lcResult.parameter[0];
    const t1 = lcResult.parameter[1];

    if (t0 > 1) {
        // The segment endpoint p[1] is the closest point to the circle.
        updateForEndpoint(segment.p[1], 1, circle, lcResult);
    }
    else if (t1 < 0) {
        // The segment endpoint p[0] is the closest point to the circle.
        updateForEndpoint(segment.p[0], 0, circle, lcResult);
    }
    else if (t0 < 0 && t1 < 1) {
        // The segment overlaps the t1-point. Remove the t0-point.
        lcResult.numClosestPairs = 1;
        lcResult.parameter[0] = lcResult.parameter[1];
        lcResult.parameter[1] = 0;
        lcResult.closest[0][0] = lcResult.closest[1][0];
        lcResult.closest[0][1] = lcResult.closest[1][1];
        lcResult.closest[1][0] = new Vector(2);
        lcResult.closest[1][1] = new Vector(2);
    }
    else if (t0 > 0 && t1 > 1) {
        // The segment overlaps the t0-point. Remove the t1-point.
        lcResult.numClosestPairs = 1;
        lcResult.parameter[1] = 0;
        lcResult.closest[1][0] = new Vector(2);
        lcResult.closest[1][1] = new Vector(2);
    }
    else if (t0 < 0 && t1 > 1) {
        // The segment is strictly inside the circle. Remove both the
        // t0-point and the t1-point.
        lcResult.distance = 0;
        lcResult.sqrDistance = 0;
        lcResult.numClosestPairs = 0;
        lcResult.parameter[0] = 0;
        lcResult.parameter[1] = 0;
        lcResult.closest[0][0] = new Vector(2);
        lcResult.closest[0][1] = new Vector(2);
        lcResult.closest[1][0] = new Vector(2);
        lcResult.closest[1][1] = new Vector(2);
    }
    // Otherwise 0 <= t0 < t1 <= 1 and the line-circle intersection points are
    // contained by the segment.
}

export class DistSegment2Circle2
    implements DCPQuery<Segment2, Circle2, DistSegment2Circle2Result> {
    compute(segment: Segment2, circle: Circle2): DistSegment2Circle2Result {
        // Execute the query for line-circle.
        const line = Line.fromOriginDirection(segment.p[0],
            sub(segment.p[1], segment.p[0]));
        const lcResult: DistSegment2Circle2Result =
            new DistLine2Circle2().compute(line, circle);

        // Restrict the analysis to segment-circle.
        if (lcResult.numClosestPairs === 2) {
            // The segment connecting the line-circle intersection points has
            // parameter interval [t0,t1]. Determine how this intersects with
            // the segment interval [0,1] and modify lcResult accordingly.
            updateForSegment(segment, circle, lcResult);
        }
        else {
            // The line does not intersect the circle or is tangent to the
            // circle. If the closest point on the line has a parameter not in
            // [0,1], then a segment endpoint is the closest segment point to
            // the circle.
            if (lcResult.parameter[0] < 0) {
                updateForEndpoint(segment.p[0], 0, circle, lcResult);
            }
            else if (lcResult.parameter[0] > 1) {
                updateForEndpoint(segment.p[1], 1, circle, lcResult);
            }
        }

        return lcResult;
    }
}
