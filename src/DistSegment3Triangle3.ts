// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistSegment3Triangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a segment and a solid triangle in 3D.
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
// upstream specialization 'DCPQuery<T, Segment3<T>, Triangle3<T>>' becomes
// the class DistSegment3Triangle3. As upstream does, the result type is the
// line-triangle result type, re-exported here as the alias
// DistSegment3Triangle3Result.

import type { DCPQuery } from './DCPQuery';
import { DistLine3Triangle3 } from './DistLine3Triangle3';
import type { DistLine3Triangle3Result } from './DistLine3Triangle3';
import { DistPointTriangle } from './DistPointTriangle';
import { Line } from './Line';
import type { Segment3 } from './Segment';
import type { Triangle3 } from './Triangle';
import { Vector, sub } from './Vector';

// Upstream reuses the line-triangle result type ('using Result = typename
// LTQuery::Result').
export type DistSegment3Triangle3Result = DistLine3Triangle3Result;

export class DistSegment3Triangle3
    implements DCPQuery<Segment3, Triangle3, DistSegment3Triangle3Result> {
    compute(segment: Segment3, triangle: Triangle3):
        DistSegment3Triangle3Result {
        const result: DistSegment3Triangle3Result = {
            distance: 0,
            sqrDistance: 0,
            parameter: 0,
            barycentric: [0, 0, 0],
            closest: [new Vector(3), new Vector(3)]
        };

        // Compute the distance between the triangle and the line containing
        // the segment. If the closest line point has a parameter in [0,1], it
        // is on the segment and the line result is the segment result.
        // Otherwise the squared distance is a convex function of the line
        // parameter, so clamping the parameter to the nearest endpoint and
        // re-solving as a point-triangle query gives the segment result.
        const segDirection = sub(segment.p[1], segment.p[0]);
        const line = Line.fromOriginDirection(segment.p[0], segDirection);
        const ltQuery = new DistLine3Triangle3();
        const ltOutput = ltQuery.compute(line, triangle);
        if (ltOutput.parameter >= 0) {
            if (ltOutput.parameter <= 1) {
                return ltOutput;
            }

            const ptQuery = new DistPointTriangle();
            const ptOutput = ptQuery.compute(segment.p[1], triangle);
            result.sqrDistance = ptOutput.sqrDistance;
            result.distance = ptOutput.distance;
            result.parameter = 1;
            result.barycentric = [ptOutput.barycentric[0],
                ptOutput.barycentric[1], ptOutput.barycentric[2]];
            result.closest[0] = segment.p[1].clone();
            result.closest[1] = ptOutput.closest[1];
        }
        else {
            const ptQuery = new DistPointTriangle();
            const ptOutput = ptQuery.compute(segment.p[0], triangle);
            result.sqrDistance = ptOutput.sqrDistance;
            result.distance = ptOutput.distance;
            result.parameter = 0;
            result.barycentric = [ptOutput.barycentric[0],
                ptOutput.barycentric[1], ptOutput.barycentric[2]];
            result.closest[0] = segment.p[0].clone();
            result.closest[1] = ptOutput.closest[1];
        }
        return result;
    }
}
