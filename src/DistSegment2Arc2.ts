// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistSegment2Arc2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a segment and an arc in 2D.
//
// The segment has endpoints P0 and P1 and is parameterized by
// P0 + t * (P1 - P0). The t-value satisfies 0 <= t <= 1.
//
// The circle containing the arc has center C and radius r. The arc has two
// endpoints E0 and E1 on the circle so that E1 is obtained from E0 by
// traversing counterclockwise. The application is responsible for ensuring
// that E0 and E1 are on the circle and that they are properly ordered.
//
// The number of pairs of closest points is result.numClosestPairs which is
// 1 or 2. If result.numClosestPairs is 1, result.parameter[0] is the segment
// t-value for its closest point result.closest[0][0]. The arc closest point
// is result.closest[0][1]. If result.numClosestPairs is 2,
// result.parameter[0] and result.parameter[1] are the segment t-values for
// its closest points result.closest[0][0] and result.closest[1][0]. The arc
// closest points are result.closest[0][1] and result.closest[1][1].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Segment2<T>, Arc2<T>>' becomes the
// class DistSegment2Arc2. As upstream does, the result type is the
// line-circle result type, re-exported here as the alias
// DistSegment2Arc2Result. The private 'SortItem' struct becomes the
// module-private interface SortItem; upstream's 'std::sort' with 'operator<'
// on sqrDistance becomes Array.sort with the equivalent comparator
// (JavaScript's sort is stable, so equal-distance items keep their
// construction order; std::sort leaves that order unspecified). The upstream
// single-argument 'Arc2::Contains' (which assumes the point is on the circle)
// is the port's 'containsOnCircle'.

import type { Arc2 } from './Arc2';
import type { DCPQuery } from './DCPQuery';
import type { DistLine2Circle2Result } from './DistLine2Circle2';
import { DistPoint2Arc2 } from './DistPoint2Arc2';
import { DistPointSegment } from './DistPointSegment';
import { DistSegment2Circle2 } from './DistSegment2Circle2';
import { Hypersphere } from './Hypersphere';
import type { Segment2 } from './Segment';
import { Vector } from './Vector';

// Upstream reuses the line-circle result type ('using Result = typename
// LCQuery::Result').
export type DistSegment2Arc2Result = DistLine2Circle2Result;

interface SortItem {
    distance: number;
    sqrDistance: number;
    parameter: number;
    closest: [Vector, Vector];
}

export class DistSegment2Arc2
    implements DCPQuery<Segment2, Arc2, DistSegment2Arc2Result> {
    compute(segment: Segment2, arc: Arc2): DistSegment2Arc2Result {
        const result: DistSegment2Arc2Result = {
            distance: 0,
            sqrDistance: 0,
            numClosestPairs: 0,
            parameter: [0, 0],
            closest: [
                [new Vector(2), new Vector(2)],
                [new Vector(2), new Vector(2)]
            ]
        };

        // Execute the query for segment-circle. Test whether the circle
        // closest points are on or off the arc. If any closest point is on
        // the arc, there is no need to test arc endpoints for closeness.
        const circle = Hypersphere.fromCenterRadius(arc.center, arc.radius);
        const scResult = new DistSegment2Circle2().compute(segment, circle);
        for (let i = 0; i < scResult.numClosestPairs; ++i) {
            if (arc.containsOnCircle(scResult.closest[i][1])) {
                const j = result.numClosestPairs++;
                result.distance = scResult.distance;
                result.sqrDistance = scResult.sqrDistance;
                result.parameter[j] = scResult.parameter[i];
                result.closest[j][0] = scResult.closest[i][0];
                result.closest[j][1] = scResult.closest[i][1];
            }
        }

        if (result.numClosestPairs > 0) {
            // At least one circle closest point is on the arc. There is no
            // need to test arc endpoints.
            return result;
        }

        // No circle closest points are on the arc. Compute distances to the
        // arc endpoints and from segment endpoints to the arc and then select
        // the minima.
        const psQuery = new DistPointSegment();
        const paQuery = new DistPoint2Arc2();
        const psResult0 = psQuery.compute(arc.end[0], segment);
        const psResult1 = psQuery.compute(arc.end[1], segment);
        const paResult2 = paQuery.compute(segment.p[0], arc);
        const paResult3 = paQuery.compute(segment.p[1], arc);

        const items: SortItem[] = [
            {
                distance: Math.sqrt(psResult0.sqrDistance),
                sqrDistance: psResult0.sqrDistance,
                parameter: psResult0.parameter,
                closest: [psResult0.closest[1], arc.end[0].clone()]
            },
            {
                distance: Math.sqrt(psResult1.sqrDistance),
                sqrDistance: psResult1.sqrDistance,
                parameter: psResult1.parameter,
                closest: [psResult1.closest[1], arc.end[1].clone()]
            },
            {
                distance: paResult2.distance,
                sqrDistance: paResult2.sqrDistance,
                parameter: 0,
                closest: [paResult2.closest[0], paResult2.closest[1]]
            },
            {
                distance: paResult3.distance,
                sqrDistance: paResult3.sqrDistance,
                parameter: 1,
                closest: [paResult3.closest[0], paResult3.closest[1]]
            }
        ];
        items.sort((a, b) => a.sqrDistance - b.sqrDistance);

        const item0 = items[0];
        const item1 = items[1];
        if (item0.sqrDistance < item1.sqrDistance
            || item0.closest[1].equals(item1.closest[1])) {
            // The arc point closest to the segment is unique.
            result.distance = item0.distance;
            result.sqrDistance = item0.sqrDistance;
            result.numClosestPairs = 1;
            result.parameter[0] = item0.parameter;
            result.closest[0][0] = item0.closest[0];
            result.closest[0][1] = item0.closest[1];
        }
        else {
            // Two arc points are equidistant from the segment.
            result.distance = item0.distance;
            result.sqrDistance = item0.sqrDistance;
            result.numClosestPairs = 2;
            result.parameter[0] = item0.parameter;
            result.parameter[1] = item1.parameter;
            result.closest[0][0] = item0.closest[0];
            result.closest[0][1] = item0.closest[1];
            result.closest[1][0] = item1.closest[0];
            result.closest[1][1] = item1.closest[1];
        }
        return result;
    }
}
