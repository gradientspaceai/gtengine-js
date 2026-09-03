// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistSegment3Rectangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a segment and a solid rectangle in 3D.
//
// The segment is P0 + t * (P1 - P0) for 0 <= t <= 1. The direction D = P1-P0
// is generally not unit length.
//
// The rectangle has center C, unit-length axis directions W[0] and W[1], and
// extents e[0] and e[1]. A rectangle point is X = C + sum_{i=0}^1 s[i] * W[i]
// where |s[i]| <= e[i] for all i.
//
// The closest point on the segment is stored in closest[0] with parameter t.
// The closest point on the rectangle is closest[1] with W-coordinates
// (s[0],s[1]). When there are infinitely many choices for the pair of closest
// points, only one of them is returned.
//
// Upstream TODO: modify to support non-unit-length W[].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Segment3<T>, Rectangle3<T>>' becomes
// the class DistSegment3Rectangle3. As upstream does, the result type is the
// line-rectangle result type, re-exported here as the alias
// DistSegment3Rectangle3Result.

import type { DCPQuery } from './DCPQuery.js';
import { DistLine3Rectangle3 } from './DistLine3Rectangle3.js';
import type { DistLine3Rectangle3Result } from './DistLine3Rectangle3.js';
import { DistPointRectangle } from './DistPointRectangle.js';
import { Line } from './Line.js';
import type { Rectangle3 } from './Rectangle.js';
import type { Segment3 } from './Segment.js';
import { Vector, sub } from './Vector.js';

// Upstream reuses the line-rectangle result type ('using Result = typename
// LRQuery::Result').
export type DistSegment3Rectangle3Result = DistLine3Rectangle3Result;

export class DistSegment3Rectangle3
    implements DCPQuery<Segment3, Rectangle3,
    DistSegment3Rectangle3Result> {
    compute(segment: Segment3, rectangle: Rectangle3):
        DistSegment3Rectangle3Result {
        const result: DistSegment3Rectangle3Result = {
            distance: 0,
            sqrDistance: 0,
            parameter: 0,
            cartesian: [0, 0],
            closest: [new Vector(3), new Vector(3)]
        };

        // Compute the distance between the rectangle and the line containing
        // the segment. If the closest line point has a parameter in [0,1], it
        // is on the segment and the line result is the segment result.
        // Otherwise the squared distance is a convex function of the line
        // parameter, so clamping the parameter to the nearest endpoint and
        // re-solving as a point-rectangle query gives the segment result.
        const segDirection = sub(segment.p[1], segment.p[0]);
        const line = Line.fromOriginDirection(segment.p[0], segDirection);
        const lrQuery = new DistLine3Rectangle3();
        const lrResult = lrQuery.compute(line, rectangle);
        if (lrResult.parameter >= 0) {
            if (lrResult.parameter <= 1) {
                return lrResult;
            }

            const prQuery = new DistPointRectangle();
            const prResult = prQuery.compute(segment.p[1], rectangle);
            result.distance = prResult.distance;
            result.sqrDistance = prResult.sqrDistance;
            result.parameter = 1;
            result.cartesian = [prResult.cartesian[0], prResult.cartesian[1]];
            result.closest[0] = segment.p[1].clone();
            result.closest[1] = prResult.closest[1];
        }
        else {
            const prQuery = new DistPointRectangle();
            const prResult = prQuery.compute(segment.p[0], rectangle);
            result.distance = prResult.distance;
            result.sqrDistance = prResult.sqrDistance;
            result.parameter = 0;
            result.cartesian = [prResult.cartesian[0], prResult.cartesian[1]];
            result.closest[0] = segment.p[0].clone();
            result.closest[1] = prResult.closest[1];
        }
        return result;
    }
}
