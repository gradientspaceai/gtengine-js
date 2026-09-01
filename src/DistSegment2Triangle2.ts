// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistSegment2Triangle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a segment and a solid triangle in 2D. (The
// upstream file comment says "in 3D", which is stale; the specialization is
// for Segment2 and Triangle2.)
//
// The segment is P0 + t * (P1 - P0) for 0 <= t <= 1. The direction D = P1-P0
// is generally not unit length.
//
// The triangle has vertices <V[0],V[1],V[2]>. A triangle point is
// X = sum_{i=0}^2 b[i] * V[i], where 0 <= b[i] <= 1 for all i and
// sum_{i=0}^2 b[i] = 1.
//
// The closest point on the segment is stored in closest[0] with parameter t.
// The closest point on the triangle is closest[1] with barycentric
// coordinates (b[0],b[1],b[2]). When there are infinitely many choices for
// the pair of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Segment2<T>, Triangle2<T>>' becomes
// the class DistSegment2Triangle2. As upstream does, the result type is the
// line-triangle result type, re-exported here as the alias
// DistSegment2Triangle2Result.

import type { DCPQuery } from './DCPQuery';
import { DistLine2Triangle2 } from './DistLine2Triangle2';
import type { DistLine2Triangle2Result } from './DistLine2Triangle2';
import { DistPointTriangle } from './DistPointTriangle';
import { Line } from './Line';
import type { Segment2 } from './Segment';
import type { Triangle2 } from './Triangle';
import { Vector, sub } from './Vector';

// Upstream reuses the line-triangle result type ('using Result = typename
// LTQuery::Result').
export type DistSegment2Triangle2Result = DistLine2Triangle2Result;

export class DistSegment2Triangle2
    implements DCPQuery<Segment2, Triangle2, DistSegment2Triangle2Result> {
    compute(segment: Segment2, triangle: Triangle2):
        DistSegment2Triangle2Result {
        const result: DistSegment2Triangle2Result = {
            distance: 0,
            sqrDistance: 0,
            parameter: 0,
            barycentric: [0, 0, 0],
            closest: [new Vector(2), new Vector(2)]
        };

        // Compute the distance between the triangle and the line containing
        // the segment. If the closest line point has a parameter in [0,1], it
        // is on the segment and the line result is the segment result.
        // Otherwise the squared distance is a convex function of the line
        // parameter, so clamping the parameter to the nearest endpoint and
        // re-solving as a point-triangle query gives the segment result.
        const segDirection = sub(segment.p[1], segment.p[0]);
        const line = Line.fromOriginDirection(segment.p[0], segDirection);
        const ltQuery = new DistLine2Triangle2();
        const ltResult = ltQuery.compute(line, triangle);
        if (ltResult.parameter >= 0) {
            if (ltResult.parameter <= 1) {
                return ltResult;
            }

            const ptQuery = new DistPointTriangle();
            const ptResult = ptQuery.compute(segment.p[1], triangle);
            result.distance = ptResult.distance;
            result.sqrDistance = ptResult.sqrDistance;
            result.parameter = 1;
            result.barycentric = [
                ptResult.barycentric[0],
                ptResult.barycentric[1],
                ptResult.barycentric[2]
            ];
            result.closest[0] = segment.p[1].clone();
            result.closest[1] = ptResult.closest[1];
        }
        else {
            const ptQuery = new DistPointTriangle();
            const ptResult = ptQuery.compute(segment.p[0], triangle);
            result.distance = ptResult.distance;
            result.sqrDistance = ptResult.sqrDistance;
            result.parameter = 0;
            result.barycentric = [
                ptResult.barycentric[0],
                ptResult.barycentric[1],
                ptResult.barycentric[2]
            ];
            result.closest[0] = segment.p[0].clone();
            result.closest[1] = ptResult.closest[1];
        }
        return result;
    }
}
