// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRay3Triangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a ray and a triangle in 3D.
//
// The ray is P + t * D for t >= 0, where D is not required to be unit length.
//
// The triangle has vertices <V[0],V[1],V[2]>. A triangle point is
// X = b[0] * V[0] + b[1] * V[1] + b[2] * V[2], where 0 <= b[i] <= 1 for all i
// and b[0] + b[1] + b[2] = 1.
//
// The closest point on the ray is stored in closest[0] with parameter t. The
// closest point on the triangle is closest[1] with barycentric coordinates
// (b[0],b[1],b[2]). When there are infinitely many choices for the pair of
// closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Ray3<T>, Triangle3<T>>' becomes the
// class DistRay3Triangle3. As upstream does, the result type is the
// line-triangle result type, re-exported here as the alias
// DistRay3Triangle3Result.

import type { DCPQuery } from './DCPQuery.js';
import { DistLine3Triangle3 } from './DistLine3Triangle3.js';
import type { DistLine3Triangle3Result } from './DistLine3Triangle3.js';
import { DistPointTriangle } from './DistPointTriangle.js';
import { Line } from './Line.js';
import type { Ray3 } from './Ray.js';
import type { Triangle3 } from './Triangle.js';

// Upstream reuses the line-triangle result type ('using Result = typename
// LTQuery::Result').
export type DistRay3Triangle3Result = DistLine3Triangle3Result;

export class DistRay3Triangle3
    implements DCPQuery<Ray3, Triangle3, DistRay3Triangle3Result> {
    compute(ray: Ray3, triangle: Triangle3): DistRay3Triangle3Result {
        // Compute the distance between the triangle and the line containing
        // the ray. If the closest line point has a nonnegative parameter, it
        // is on the ray and the line result is the ray result. Otherwise the
        // squared distance is a convex function of the line parameter, so
        // clamping the parameter to 0 and re-solving as a point-triangle
        // query gives the ray result.
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const lrResult = new DistLine3Triangle3().compute(line, triangle);
        if (lrResult.parameter >= 0) {
            return lrResult;
        }

        const prResult = new DistPointTriangle().compute(ray.origin, triangle);
        return {
            distance: prResult.distance,
            sqrDistance: prResult.sqrDistance,
            parameter: 0,
            barycentric: [prResult.barycentric[0],
                prResult.barycentric[1], prResult.barycentric[2]],
            closest: [ray.origin.clone(), prResult.closest[1]]
        };
    }
}
