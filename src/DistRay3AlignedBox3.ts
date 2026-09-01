// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRay3AlignedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a ray and a solid aligned box in 3D.
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
// upstream specialization 'DCPQuery<T, Ray3<T>, AlignedBox3<T>>' becomes the
// class DistRay3AlignedBox3. As upstream does, the result type is the
// line-box result type, re-exported here as the alias
// DistRay3AlignedBox3Result.

import type { AlignedBox3 } from './AlignedBox';
import type { DCPQuery } from './DCPQuery';
import { DistLine3AlignedBox3 } from './DistLine3AlignedBox3';
import type { DistLine3AlignedBox3Result } from './DistLine3AlignedBox3';
import { DistPointAlignedBox } from './DistPointAlignedBox';
import { Line } from './Line';
import type { Ray3 } from './Ray';

// Upstream reuses the line-box result type ('using Result = typename
// AlignedQuery::Result').
export type DistRay3AlignedBox3Result = DistLine3AlignedBox3Result;

export class DistRay3AlignedBox3
    implements DCPQuery<Ray3, AlignedBox3, DistRay3AlignedBox3Result> {
    compute(ray: Ray3, box: AlignedBox3): DistRay3AlignedBox3Result {
        // Compute the distance between the box and the line containing
        // the ray. If the closest line point has a nonnegative parameter, it
        // is on the ray and the line result is the ray result. Otherwise the
        // squared distance is a convex function of the line parameter, so
        // clamping the parameter to 0 and re-solving as a point-box
        // query gives the ray result.
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const lrResult = new DistLine3AlignedBox3().compute(line, box);
        if (lrResult.parameter >= 0) {
            return lrResult;
        }

        const prResult = new DistPointAlignedBox().compute(ray.origin, box);
        return {
            distance: prResult.distance,
            sqrDistance: prResult.sqrDistance,
            parameter: 0,
            closest: [ray.origin.clone(), prResult.closest[1]]
        };
    }
}
