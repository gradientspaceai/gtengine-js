// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistSegment2OrientedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a segment and a solid oriented box in 2D.
//
// The segment is P0 + t * (P1 - P0) for 0 <= t <= 1. The direction D = P1-P0
// is generally not unit length.
//
// The oriented box has center C, unit-length axis directions U[i] and extents
// e[i] for 0 <= i < 2. A box point is X = C + sum_i y[i] * U[i], where
// |y[i]| <= e[i] for all i.
//
// The closest point on the segment is stored in closest[0] with parameter t.
// The closest point on the box is stored in closest[1]. When there are
// infinitely many choices for the pair of closest points, only one of them is
// returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Segment2<T>, OrientedBox2<T>>' becomes
// the class DistSegment2OrientedBox2. As upstream does, the result type is
// the line-box result type, re-exported here as the alias
// DistSegment2OrientedBox2Result.

import type { DCPQuery } from './DCPQuery';
import { DistLine2OrientedBox2 } from './DistLine2OrientedBox2';
import type { DistLine2OrientedBox2Result } from './DistLine2OrientedBox2';
import { DistPointOrientedBox } from './DistPointOrientedBox';
import { Line } from './Line';
import type { OrientedBox2 } from './OrientedBox';
import type { Segment2 } from './Segment';
import { Vector, sub } from './Vector';

// Upstream reuses the line-box result type ('using Result = typename
// OrientedQuery::Result').
export type DistSegment2OrientedBox2Result = DistLine2OrientedBox2Result;

export class DistSegment2OrientedBox2
    implements DCPQuery<Segment2, OrientedBox2, DistSegment2OrientedBox2Result> {
    compute(segment: Segment2, box: OrientedBox2): DistSegment2OrientedBox2Result {
        const result: DistSegment2OrientedBox2Result = {
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
        const lbResult = new DistLine2OrientedBox2().compute(line, box);
        if (lbResult.parameter >= 0) {
            if (lbResult.parameter <= 1) {
                return lbResult;
            }

            const pbResult = new DistPointOrientedBox().compute(segment.p[1], box);
            result.sqrDistance = pbResult.sqrDistance;
            result.distance = pbResult.distance;
            result.parameter = 1;
            result.closest[0] = segment.p[1].clone();
            result.closest[1] = pbResult.closest[1];
        }
        else {
            const pbResult = new DistPointOrientedBox().compute(segment.p[0], box);
            result.sqrDistance = pbResult.sqrDistance;
            result.distance = pbResult.distance;
            result.parameter = 0;
            result.closest[0] = segment.p[0].clone();
            result.closest[1] = pbResult.closest[1];
        }

        return result;
    }
}
