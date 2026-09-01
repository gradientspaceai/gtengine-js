// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRay2AlignedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a ray and a solid aligned box in 2D.
//
// The ray is P + t * D for t >= 0, where D is not required to be unit length.
//
// The aligned box has minimum corner A and maximum corner B. A box point is X
// where A <= X <= B; the comparisons are componentwise.
//
// The closest point on the ray is stored in closest[0] with parameter t. The
// closest point on the box is stored in closest[1]. When there are infinitely
// many choices for the pair of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Ray2<T>, AlignedBox2<T>>' becomes the
// class DistRay2AlignedBox2. As upstream does, the result type is the
// line-box result type, re-exported here as the alias
// DistRay2AlignedBox2Result.

import type { AlignedBox2 } from './AlignedBox';
import type { DCPQuery } from './DCPQuery';
import { DistLine2AlignedBox2 } from './DistLine2AlignedBox2';
import type { DistLine2AlignedBox2Result } from './DistLine2AlignedBox2';
import { DistPointAlignedBox } from './DistPointAlignedBox';
import { Line } from './Line';
import type { Ray2 } from './Ray';

// Upstream reuses the line-box result type ('using Result = typename
// AlignedQuery::Result').
export type DistRay2AlignedBox2Result = DistLine2AlignedBox2Result;

export class DistRay2AlignedBox2
    implements DCPQuery<Ray2, AlignedBox2, DistRay2AlignedBox2Result> {
    compute(ray: Ray2, box: AlignedBox2): DistRay2AlignedBox2Result {
        // Compute the distance between the box and the line containing the
        // ray. If the closest line point has a nonnegative parameter, it is
        // on the ray and the line result is the ray result. Otherwise the
        // squared distance is a convex function of the line parameter, so
        // clamping the parameter to 0 and re-solving as a point-box query
        // gives the ray result.
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const lbResult = new DistLine2AlignedBox2().compute(line, box);
        if (lbResult.parameter >= 0) {
            return lbResult;
        }

        const pbResult = new DistPointAlignedBox().compute(ray.origin, box);
        return {
            distance: pbResult.distance,
            sqrDistance: pbResult.sqrDistance,
            parameter: 0,
            closest: [ray.origin.clone(), pbResult.closest[1]]
        };
    }
}
