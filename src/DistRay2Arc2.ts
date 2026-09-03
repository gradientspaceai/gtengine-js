// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRay2Arc2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a ray and an arc in 2D.
//
// The ray is P + t * D, where P is a point on the ray and D is not required
// to be unit length. The t-value satisfies t >= 0.
//
// The circle containing the arc has center C and radius r. The arc has two
// endpoints E0 and E1 on the circle so that E1 is obtained from E0 by
// traversing counterclockwise. The application is responsible for ensuring
// that E0 and E1 are on the circle and that they are properly ordered.
//
// The number of pairs of closest points is result.numClosestPairs which is
// 1 or 2. If result.numClosestPairs is 1, result.parameter[0] is the ray
// t-value for its closest point result.closest[0][0]. The arc closest point
// is result.closest[0][1]. If result.numClosestPairs is 2,
// result.parameter[0] and result.parameter[1] are the ray t-values for its
// closest points result.closest[0][0] and result.closest[1][0]. The arc
// closest points are result.closest[0][1] and result.closest[1][1].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Ray2<T>, Arc2<T>>' becomes the class
// DistRay2Arc2. As upstream does, the result type is the line-circle result
// type, re-exported here as the alias DistRay2Arc2Result. The private
// 'SortItem' struct becomes the module-private interface SortItem; upstream's
// 'std::sort' with 'operator<' on sqrDistance becomes Array.sort with the
// equivalent comparator (JavaScript's sort is stable, so equal-distance items
// keep their construction order; std::sort leaves that order unspecified).
// The upstream single-argument 'Arc2::Contains' (which assumes the point is
// on the circle) is the port's 'containsOnCircle'.

import type { Arc2 } from './Arc2.js';
import type { DCPQuery } from './DCPQuery.js';
import type { DistLine2Circle2Result } from './DistLine2Circle2.js';
import { DistPoint2Arc2 } from './DistPoint2Arc2.js';
import { DistPointRay } from './DistPointRay.js';
import { DistRay2Circle2 } from './DistRay2Circle2.js';
import { Hypersphere } from './Hypersphere.js';
import type { Ray2 } from './Ray.js';
import { Vector } from './Vector.js';

// Upstream reuses the line-circle result type ('using Result = typename
// LCQuery::Result').
export type DistRay2Arc2Result = DistLine2Circle2Result;

interface SortItem {
    distance: number;
    sqrDistance: number;
    parameter: number;
    closest: [Vector, Vector];
}

export class DistRay2Arc2 implements DCPQuery<Ray2, Arc2, DistRay2Arc2Result> {
    compute(ray: Ray2, arc: Arc2): DistRay2Arc2Result {
        const result: DistRay2Arc2Result = {
            distance: 0,
            sqrDistance: 0,
            numClosestPairs: 0,
            parameter: [0, 0],
            closest: [
                [new Vector(2), new Vector(2)],
                [new Vector(2), new Vector(2)]
            ]
        };

        // Execute the query for ray-circle. Test whether the circle closest
        // points are on or off the arc. If any closest point is on the arc,
        // there is no need to test arc endpoints for closeness.
        const circle = Hypersphere.fromCenterRadius(arc.center, arc.radius);
        const rcResult = new DistRay2Circle2().compute(ray, circle);
        for (let i = 0; i < rcResult.numClosestPairs; ++i) {
            if (arc.containsOnCircle(rcResult.closest[i][1])) {
                const j = result.numClosestPairs++;
                result.distance = rcResult.distance;
                result.sqrDistance = rcResult.sqrDistance;
                result.parameter[j] = rcResult.parameter[i];
                result.closest[j][0] = rcResult.closest[i][0];
                result.closest[j][1] = rcResult.closest[i][1];
            }
        }

        if (result.numClosestPairs > 0) {
            // At least one circle closest point is on the arc. There is no
            // need to test arc endpoints.
            return result;
        }

        // No circle closest points are on the arc. Compute distances to the
        // arc endpoints and from ray origin to the arc and then select the
        // minima.
        const prQuery = new DistPointRay();
        const paQuery = new DistPoint2Arc2();
        const prResult0 = prQuery.compute(arc.end[0], ray);
        const prResult1 = prQuery.compute(arc.end[1], ray);
        const paResult2 = paQuery.compute(ray.origin, arc);

        const items: SortItem[] = [
            {
                distance: Math.sqrt(prResult0.sqrDistance),
                sqrDistance: prResult0.sqrDistance,
                parameter: prResult0.parameter,
                closest: [prResult0.closest[1], arc.end[0].clone()]
            },
            {
                distance: Math.sqrt(prResult1.sqrDistance),
                sqrDistance: prResult1.sqrDistance,
                parameter: prResult1.parameter,
                closest: [prResult1.closest[1], arc.end[1].clone()]
            },
            {
                distance: paResult2.distance,
                sqrDistance: paResult2.sqrDistance,
                parameter: 0,
                closest: [paResult2.closest[0], paResult2.closest[1]]
            }
        ];
        items.sort((a, b) => a.sqrDistance - b.sqrDistance);

        const item0 = items[0];
        const item1 = items[1];
        if (item0.sqrDistance < item1.sqrDistance
            || item0.closest[1].equals(item1.closest[1])) {
            // The arc point closest to the ray is unique.
            result.distance = item0.distance;
            result.sqrDistance = item0.sqrDistance;
            result.numClosestPairs = 1;
            result.parameter[0] = item0.parameter;
            result.closest[0][0] = item0.closest[0];
            result.closest[0][1] = item0.closest[1];
        }
        else {
            // Two arc points are equidistant from the ray.
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
