// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLine3Triangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a line and a solid triangle in 3D.
//
// The line is P + t * D, where D is not required to be unit length.
//
// The triangle has vertices <V[0],V[1],V[2]>. A triangle point is
// X = sum_{i=0}^2 b[i] * V[i], where 0 <= b[i] <= 1 for all i and
// sum_{i=0}^2 b[i] = 1.
//
// The closest point on the line is stored in closest[0] with parameter t. The
// closest point on the triangle is closest[1] with barycentric coordinates
// (b[0],b[1],b[2]). When there are infinitely many choices for the pair of
// closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Line3<T>, Triangle3<T>>' becomes the
// class DistLine3Triangle3 with the result type DistLine3Triangle3Result.

import type { DCPQuery } from './DCPQuery.js';
import { DistLineSegment } from './DistLineSegment.js';
import type { Line3 } from './Line.js';
import { Segment } from './Segment.js';
import type { Triangle3 } from './Triangle.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import { cross } from './Vector3.js';

export interface DistLine3Triangle3Result {
    distance: number;
    sqrDistance: number;

    // The line parameter t of the closest line point.
    parameter: number;

    // The barycentric coordinates of closest[1].
    barycentric: [number, number, number];

    // closest[0] is on the line, closest[1] is on the triangle.
    closest: [Vector, Vector];
}

export class DistLine3Triangle3
    implements DCPQuery<Line3, Triangle3, DistLine3Triangle3Result> {
    compute(line: Line3, triangle: Triangle3): DistLine3Triangle3Result {
        // The line points are X = P + t * D and the triangle points are
        // Y = b[0] * V[0] + b[1] * V[1] + b[2] * V[2], where the barycentric
        // coordinates satisfy b[i] in [0,1] and b[0] + b[1] + b[2] = 1.
        // Define the triangle edge directions by E[1] = V[1] - V[0] and
        // E[2] = V[2] - V[0]; then Y = V[0] + b1 * E[1] + b2 * E[2]. If Y is
        // specified the barycentric coordinates are the solution to
        //
        // +-                        -+ +-    -+   +-                 -+
        // | Dot(E1, E1)  Dot(E1, E2) | | b[1] | = | Dot(E1, Y - V[0]) |
        // | Dot(E1, E2)  Dot(E2, E2) | | b[2] |   | Dot(E2, Y - V[0]) |
        // +-                        -+ +-    -+   +-                 -+
        //
        // and b[0] = 1 - b[1] - b[2].
        const result: DistLine3Triangle3Result = {
            distance: 0,
            sqrDistance: 0,
            parameter: 0,
            barycentric: [0, 0, 0],
            closest: [new Vector(3), new Vector(3)]
        };

        // Test whether the line intersects the triangle. If so, the squared
        // distance is zero. The normal of the plane of the triangle does not
        // have to be normalized to unit length.
        const E1 = sub(triangle.v[1], triangle.v[0]);
        const E2 = sub(triangle.v[2], triangle.v[0]);
        const N = cross(E1, E2);
        const NdD = dot(N, line.direction);
        if (Math.abs(NdD) > 0) {
            // The line and triangle are not parallel, so the line intersects
            // the plane of the triangle at a point Y. Determine whether Y is
            // contained by the triangle.
            const PmV0 = sub(line.origin, triangle.v[0]);
            const NdDiff = dot(N, PmV0);
            const tIntersect = -NdDiff / NdD;
            const Y = add(line.origin, mul(tIntersect, line.direction));
            const YmV0 = sub(Y, triangle.v[0]);

            // Compute the barycentric coordinates of the intersection.
            const E1dE1 = dot(E1, E1);
            const E1dE2 = dot(E1, E2);
            const E2dE2 = dot(E2, E2);
            const E1dYmV0 = dot(E1, YmV0);
            const E2dYmV0 = dot(E2, YmV0);
            const det = E1dE1 * E2dE2 - E1dE2 * E1dE2;
            const b1 = (E2dE2 * E1dYmV0 - E1dE2 * E2dYmV0) / det;
            const b2 = (E1dE1 * E2dYmV0 - E1dE2 * E1dYmV0) / det;
            const b0 = 1 - b1 - b2;

            if (b0 >= 0 && b1 >= 0 && b2 >= 0) {
                // The point Y is contained by the triangle.
                result.sqrDistance = 0;
                result.distance = 0;
                result.parameter = tIntersect;
                result.barycentric[0] = b0;
                result.barycentric[1] = b1;
                result.barycentric[2] = b2;
                result.closest[0] = Y;
                result.closest[1] = Y.clone();
                return result;
            }
        }

        // Either (1) the line is not parallel to the triangle and the point
        // of intersection of the line and the plane of the triangle is
        // outside the triangle or (2) the line and triangle are parallel.
        // Regardless, the closest point on the triangle is on an edge of the
        // triangle. Compare the line to all three edges of the triangle. To
        // allow for arbitrary precision arithmetic, the initial distance and
        // sqrDistance are initialized to a negative number rather than a
        // floating-point maximum value. Tracking the minimum requires a small
        // amount of extra logic.
        const invalid = -1;
        result.distance = invalid;
        result.sqrDistance = invalid;

        const lsQuery = new DistLineSegment();
        const segment = new Segment(3);

        for (let i0 = 2, i1 = 0, i2 = 1; i1 < 3;) {
            segment.p[0] = triangle.v[i0].clone();
            segment.p[1] = triangle.v[i1].clone();

            const lsResult = lsQuery.compute(line, segment);
            if (result.sqrDistance === invalid ||
                lsResult.sqrDistance < result.sqrDistance) {
                result.sqrDistance = lsResult.sqrDistance;
                result.distance = lsResult.distance;
                result.parameter = lsResult.parameter[0];
                result.barycentric[i0] = 1 - lsResult.parameter[1];
                result.barycentric[i1] = lsResult.parameter[1];
                result.barycentric[i2] = 0;
                result.closest = lsResult.closest;
            }

            // The upstream loop increment 'i2 = i0, i0 = i1++'.
            i2 = i0;
            i0 = i1;
            ++i1;
        }

        return result;
    }
}
