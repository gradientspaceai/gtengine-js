// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRay3OrientedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a ray and a solid oriented box in 3D.
//
// The ray is P + t * D for t >= 0, where D is not required to be unit length.
//
// The oriented box has center C, unit-length axis directions U[i] and extents
// e[i] for all i. A box point is X = C + sum_i y[i] * U[i], where
// |y[i]| <= e[i] for all i.
//
// The closest point on the ray is stored in closest[0] with parameter t. The
// closest point on the box is stored in closest[1]. When there are infinitely
// many choices for the pair of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Ray3<T>, OrientedBox3<T>>' becomes the
// class DistRay3OrientedBox3. As upstream does, the result type is the
// line-box result type, re-exported here as the alias
// DistRay3OrientedBox3Result.

import type { DCPQuery } from './DCPQuery';
import { DistLine3OrientedBox3 } from './DistLine3OrientedBox3';
import type { DistLine3OrientedBox3Result } from './DistLine3OrientedBox3';
import { DistPointOrientedBox } from './DistPointOrientedBox';
import { Line } from './Line';
import type { OrientedBox3 } from './OrientedBox';
import type { Ray3 } from './Ray';

// Upstream reuses the line-box result type ('using Result = typename
// OrientedQuery::Result').
export type DistRay3OrientedBox3Result = DistLine3OrientedBox3Result;

export class DistRay3OrientedBox3
    implements DCPQuery<Ray3, OrientedBox3, DistRay3OrientedBox3Result> {
    compute(ray: Ray3, box: OrientedBox3): DistRay3OrientedBox3Result {
        // Compute the distance between the box and the line containing
        // the ray. If the closest line point has a nonnegative parameter, it
        // is on the ray and the line result is the ray result. Otherwise the
        // squared distance is a convex function of the line parameter, so
        // clamping the parameter to 0 and re-solving as a point-box
        // query gives the ray result.
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const lrResult = new DistLine3OrientedBox3().compute(line, box);
        if (lrResult.parameter >= 0) {
            return lrResult;
        }

        const prResult = new DistPointOrientedBox().compute(ray.origin, box);
        return {
            distance: prResult.distance,
            sqrDistance: prResult.sqrDistance,
            parameter: 0,
            closest: [ray.origin.clone(), prResult.closest[1]]
        };
    }
}
