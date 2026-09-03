// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistSegment2AlignedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a segment and a solid aligned box in 2D.
//
// The segment is P0 + t * (P1 - P0) for 0 <= t <= 1. The direction D = P1-P0
// is generally not unit length.
//
// The aligned box has minimum corner A and maximum corner B. A box point is X
// where A <= X <= B; the comparisons are componentwise.
//
// The closest point on the segment is stored in closest[0] with parameter t.
// The closest point on the box is stored in closest[1]. When there are
// infinitely many choices for the pair of closest points, only one of them is
// returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Segment2<T>, AlignedBox2<T>>' becomes
// the class DistSegment2AlignedBox2. As upstream does, the result type is the
// line-box result type, re-exported here as the alias
// DistSegment2AlignedBox2Result.

import type { AlignedBox2 } from './AlignedBox.js';
import type { DCPQuery } from './DCPQuery.js';
import { DistLine2AlignedBox2 } from './DistLine2AlignedBox2.js';
import type { DistLine2AlignedBox2Result } from './DistLine2AlignedBox2.js';
import { DistPointAlignedBox } from './DistPointAlignedBox.js';
import { Line } from './Line.js';
import type { Segment2 } from './Segment.js';
import { Vector, sub } from './Vector.js';

// Upstream reuses the line-box result type ('using Result = typename
// AlignedQuery::Result').
export type DistSegment2AlignedBox2Result = DistLine2AlignedBox2Result;

export class DistSegment2AlignedBox2
    implements DCPQuery<Segment2, AlignedBox2, DistSegment2AlignedBox2Result> {
    compute(segment: Segment2, box: AlignedBox2): DistSegment2AlignedBox2Result {
        const result: DistSegment2AlignedBox2Result = {
            distance: 0,
            sqrDistance: 0,
            parameter: 0,
            closest: [new Vector(2), new Vector(2)]
        };

        // Compute the distance between the box and the line containing the
        // segment. If the closest line point has a parameter in [0,1], it is
        // on the segment and the line result is the segment result. Otherwise
        // the squared distance is a convex function of the line parameter, so
        // clamping the parameter to the nearest endpoint and re-solving as a
        // point-box query gives the segment result.
        const direction = sub(segment.p[1], segment.p[0]);
        const line = Line.fromOriginDirection(segment.p[0], direction);
        const lbResult = new DistLine2AlignedBox2().compute(line, box);
        if (lbResult.parameter >= 0) {
            if (lbResult.parameter <= 1) {
                return lbResult;
            }

            const pbResult = new DistPointAlignedBox().compute(segment.p[1], box);
            result.sqrDistance = pbResult.sqrDistance;
            result.distance = pbResult.distance;
            result.parameter = 1;
            result.closest[0] = segment.p[1].clone();
            result.closest[1] = pbResult.closest[1];
        }
        else {
            const pbResult = new DistPointAlignedBox().compute(segment.p[0], box);
            result.sqrDistance = pbResult.sqrDistance;
            result.distance = pbResult.distance;
            result.parameter = 0;
            result.closest[0] = segment.p[0].clone();
            result.closest[1] = pbResult.closest[1];
        }

        return result;
    }
}
