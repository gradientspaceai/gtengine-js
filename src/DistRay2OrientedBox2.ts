// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRay2OrientedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a ray and a solid oriented box in 2D.
//
// The ray is P + t * D for t >= 0, where D is not required to be unit length.
//
// The oriented box has center C, unit-length axis directions U[i] and extents
// e[i] for 0 <= i < 2. A box point is X = C + sum_i y[i] * U[i], where
// |y[i]| <= e[i] for all i.
//
// The closest point on the ray is stored in closest[0] with parameter t. The
// closest point on the box is stored in closest[1]. When there are infinitely
// many choices for the pair of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Ray2<T>, OrientedBox2<T>>' becomes the
// class DistRay2OrientedBox2. As upstream does, the result type is the
// line-box result type, re-exported here as the alias
// DistRay2OrientedBox2Result.

import type { DCPQuery } from './DCPQuery';
import { DistLine2OrientedBox2 } from './DistLine2OrientedBox2';
import type { DistLine2OrientedBox2Result } from './DistLine2OrientedBox2';
import { DistPointOrientedBox } from './DistPointOrientedBox';
import { Line } from './Line';
import type { OrientedBox2 } from './OrientedBox';
import type { Ray2 } from './Ray';

// Upstream reuses the line-box result type ('using Result = typename
// OrientedQuery::Result').
export type DistRay2OrientedBox2Result = DistLine2OrientedBox2Result;

export class DistRay2OrientedBox2
    implements DCPQuery<Ray2, OrientedBox2, DistRay2OrientedBox2Result> {
    compute(ray: Ray2, box: OrientedBox2): DistRay2OrientedBox2Result {
        // Compute the distance between the box and the line containing the
        // ray. If the closest line point has a nonnegative parameter, it is
        // on the ray and the line result is the ray result. Otherwise the
        // squared distance is a convex function of the line parameter, so
        // clamping the parameter to 0 and re-solving as a point-box query
        // gives the ray result.
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const lbResult = new DistLine2OrientedBox2().compute(line, box);
        if (lbResult.parameter >= 0) {
            return lbResult;
        }

        const pbResult = new DistPointOrientedBox().compute(ray.origin, box);
        return {
            distance: pbResult.distance,
            sqrDistance: pbResult.sqrDistance,
            parameter: 0,
            closest: [ray.origin.clone(), pbResult.closest[1]]
        };
    }
}
