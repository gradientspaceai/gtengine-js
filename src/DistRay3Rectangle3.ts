// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRay3Rectangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a ray and a solid rectangle in 3D.
//
// The ray is P + t * D for t >= 0, where D is not required to be unit length.
//
// The rectangle has center C, unit-length axis directions W[0] and W[1], and
// extents e[0] and e[1]. A rectangle point is X = C + s[0] * W[0] +
// s[1] * W[1] where |s[i]| <= e[i] for all i.
//
// The closest point on the ray is stored in closest[0] with parameter t. The
// closest point on the rectangle is stored in closest[1] with W-coordinates
// (s[0],s[1]). When there are infinitely many choices for the pair of closest
// points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Ray3<T>, Rectangle3<T>>' becomes the
// class DistRay3Rectangle3. As upstream does, the result type is the
// line-rectangle result type, re-exported here as the alias
// DistRay3Rectangle3Result.

import type { DCPQuery } from './DCPQuery';
import { DistLine3Rectangle3 } from './DistLine3Rectangle3';
import type { DistLine3Rectangle3Result } from './DistLine3Rectangle3';
import { DistPointRectangle } from './DistPointRectangle';
import { Line } from './Line';
import type { Ray3 } from './Ray';
import type { Rectangle3 } from './Rectangle';

// Upstream reuses the line-rectangle result type ('using Result = typename
// LRQuery::Result').
export type DistRay3Rectangle3Result = DistLine3Rectangle3Result;

export class DistRay3Rectangle3
    implements DCPQuery<Ray3, Rectangle3, DistRay3Rectangle3Result> {
    compute(ray: Ray3, rectangle: Rectangle3): DistRay3Rectangle3Result {
        // Compute the distance between the rectangle and the line containing
        // the ray. If the closest line point has a nonnegative parameter, it
        // is on the ray and the line result is the ray result. Otherwise the
        // squared distance is a convex function of the line parameter, so
        // clamping the parameter to 0 and re-solving as a point-rectangle
        // query gives the ray result.
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const lrResult = new DistLine3Rectangle3().compute(line, rectangle);
        if (lrResult.parameter >= 0) {
            return lrResult;
        }

        const prResult = new DistPointRectangle().compute(ray.origin, rectangle);
        return {
            distance: prResult.distance,
            sqrDistance: prResult.sqrDistance,
            parameter: 0,
            cartesian: [prResult.cartesian[0], prResult.cartesian[1]],
            closest: [ray.origin.clone(), prResult.closest[1]]
        };
    }
}
