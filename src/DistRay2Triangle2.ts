// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRay2Triangle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a ray and a triangle in 2D.
//
// The ray is P + t * D for t >= 0, where D is not required to be unit length.
//
// The triangle has vertices <V[0],V[1],V[2]>. A triangle point is
// X = sum_{i=0}^2 b[i] * V[i], where 0 <= b[i] <= 1 for all i and
// sum_{i=0}^2 b[i] = 1.
//
// The closest point on the ray is stored in closest[0] with parameter t. The
// closest point on the triangle is closest[1] with barycentric coordinates
// (b[0],b[1],b[2]). When there are infinitely many choices for the pair of
// closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Ray2<T>, Triangle2<T>>' becomes the
// class DistRay2Triangle2 with the result type DistRay2Triangle2Result, which
// is structurally identical to the line-triangle result that upstream
// aliases.

import type { DCPQuery } from './DCPQuery';
import { DistLine2Triangle2 } from './DistLine2Triangle2';
import { DistPointTriangle } from './DistPointTriangle';
import { Line } from './Line';
import type { Ray2 } from './Ray';
import type { Triangle2 } from './Triangle';
import type { Vector } from './Vector';

export interface DistRay2Triangle2Result {
    distance: number;
    sqrDistance: number;

    // The ray parameter t of the closest ray point.
    parameter: number;

    // The barycentric coordinates of closest[1].
    barycentric: [number, number, number];

    // closest[0] is on the ray, closest[1] is on the triangle.
    closest: [Vector, Vector];
}

export class DistRay2Triangle2
    implements DCPQuery<Ray2, Triangle2, DistRay2Triangle2Result> {
    compute(ray: Ray2, triangle: Triangle2): DistRay2Triangle2Result {
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const ltResult = new DistLine2Triangle2().compute(line, triangle);
        if (ltResult.parameter >= 0) {
            return ltResult;
        }

        // The closest line point has a negative parameter, so the closest ray
        // point is the ray origin.
        const ptResult = new DistPointTriangle().compute(ray.origin, triangle);
        return {
            distance: ptResult.distance,
            sqrDistance: ptResult.sqrDistance,
            parameter: 0,
            barycentric: [ptResult.barycentric[0], ptResult.barycentric[1],
                ptResult.barycentric[2]],
            closest: [ray.origin.clone(), ptResult.closest[1]]
        };
    }
}
