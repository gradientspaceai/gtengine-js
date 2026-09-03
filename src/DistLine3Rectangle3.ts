// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLine3Rectangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a line and a solid rectangle in 3D.
//
// The line is P + t * D, where D is not required to be unit length.
//
// The rectangle has center C, unit-length axis directions W[0] and W[1], and
// extents e[0] and e[1]. A rectangle point is X = C + sum_{i=0}^1 s[i] * W[i]
// where |s[i]| <= e[i] for all i.
//
// The closest point on the line is stored in closest[0] with parameter t. The
// closest point on the rectangle is stored in closest[1] with W-coordinates
// (s[0],s[1]). When there are infinitely many choices for the pair of closest
// points, only one of them is returned.
//
// Upstream TODO: modify to support non-unit-length W[].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Line3<T>, Rectangle3<T>>' becomes the
// class DistLine3Rectangle3 with the result type DistLine3Rectangle3Result.
// The convenience alias DCPLine3Rectangle3 is dropped.

import type { DCPQuery } from './DCPQuery.js';
import { DistLineSegment } from './DistLineSegment.js';
import type { Line3 } from './Line.js';
import type { Rectangle3 } from './Rectangle.js';
import { Segment } from './Segment.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import { cross } from './Vector3.js';

export interface DistLine3Rectangle3Result {
    distance: number;
    sqrDistance: number;

    // The line parameter t of the closest line point.
    parameter: number;

    // The W-coordinates (s[0],s[1]) of closest[1].
    cartesian: [number, number];

    // closest[0] is on the line, closest[1] is on the rectangle.
    closest: [Vector, Vector];
}

export class DistLine3Rectangle3
    implements DCPQuery<Line3, Rectangle3, DistLine3Rectangle3Result> {
    compute(line: Line3, rectangle: Rectangle3): DistLine3Rectangle3Result {
        const result: DistLine3Rectangle3Result = {
            distance: 0,
            sqrDistance: 0,
            parameter: 0,
            cartesian: [0, 0],
            closest: [new Vector(3), new Vector(3)]
        };

        // Test whether the line intersects the rectangle. If so, the squared
        // distance is zero. The normal of the plane of the rectangle does not
        // have to be normalized to unit length.
        const N = cross(rectangle.axis[0], rectangle.axis[1]);
        const NdD = dot(N, line.direction);
        if (Math.abs(NdD) > 0) {
            // The line and rectangle are not parallel, so the line intersects
            // the plane of the rectangle at a point Y. Determine whether Y is
            // contained by the rectangle.
            const PmC = sub(line.origin, rectangle.center);
            const NdDiff = dot(N, PmC);
            const tIntersect = -NdDiff / NdD;
            const Y = add(line.origin, mul(tIntersect, line.direction));
            const YmC = sub(Y, rectangle.center);

            // Compute the rectangle coordinates of the intersection.
            const s0 = dot(rectangle.axis[0], YmC);
            const s1 = dot(rectangle.axis[1], YmC);

            if (Math.abs(s0) <= rectangle.extent.values[0] &&
                Math.abs(s1) <= rectangle.extent.values[1]) {
                // The point Y is contained by the rectangle.
                result.sqrDistance = 0;
                result.distance = 0;
                result.parameter = tIntersect;
                result.cartesian[0] = s0;
                result.cartesian[1] = s1;
                result.closest[0] = Y;
                result.closest[1] = Y.clone();
                return result;
            }
        }

        // Either (1) the line is not parallel to the rectangle and the point
        // of intersection of the line and the plane of the rectangle is
        // outside the rectangle or (2) the line and rectangle are parallel.
        // Regardless, the closest point on the rectangle is on an edge of the
        // rectangle. Compare the line to all four edges of the rectangle. To
        // allow for arbitrary precision arithmetic, the initial distance and
        // sqrDistance are initialized to a negative number rather than a
        // floating-point maximum value. Tracking the minimum requires a small
        // amount of extra logic.
        const lsQuery = new DistLineSegment();
        const segment = new Segment(3);

        const invalid = -1;
        result.distance = invalid;
        result.sqrDistance = invalid;

        const sign = [-1, 1, -1, 1];
        const j0 = [0, 0, 1, 1];
        const j1 = [1, 1, 0, 0];
        const edges: readonly (readonly [number, number])[] = [
            // horizontal edges (y = -e1 or +e1)
            [0, 1], [2, 3],
            // vertical edges (x = -e0 or +e0)
            [0, 2], [1, 3]
        ];
        const vertices = rectangle.getVertices();

        for (let i = 0; i < 4; ++i) {
            const edge = edges[i];
            segment.p[0] = vertices[edge[0]].clone();
            segment.p[1] = vertices[edge[1]].clone();

            const lsResult = lsQuery.compute(line, segment);
            if (result.sqrDistance === invalid ||
                lsResult.sqrDistance < result.sqrDistance) {
                result.sqrDistance = lsResult.sqrDistance;
                result.distance = lsResult.distance;
                result.parameter = lsResult.parameter[0];
                result.closest = lsResult.closest;

                const scale = 2 * lsResult.parameter[1] - 1;
                result.cartesian[j0[i]] =
                    scale * rectangle.extent.values[j0[i]];
                result.cartesian[j1[i]] =
                    sign[i] * rectangle.extent.values[j1[i]];
            }
        }

        return result;
    }
}
