// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPointRay.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a point and a ray in nD.
//
// The ray is P + t * D for t >= 0, where D is not required to be unit length.
//
// The input point is stored in closest[0]. The closest point on the ray is
// stored in closest[1].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector<N,T>, Ray<N,T>>' becomes the
// class DistPointRay with the result type DistPointRayResult.

import type { DCPQuery } from './DCPQuery.js';
import type { Ray } from './Ray.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

export interface DistPointRayResult {
    distance: number;
    sqrDistance: number;

    // The ray parameter t >= 0 of the closest ray point.
    parameter: number;

    // closest[0] is the input point, closest[1] is the closest ray point.
    closest: [Vector, Vector];
}

export class DistPointRay
    implements DCPQuery<Vector, Ray, DistPointRayResult> {
    compute(point: Vector, ray: Ray): DistPointRayResult {
        let diff = sub(point, ray.origin);
        let parameter = dot(ray.direction, diff)
            / dot(ray.direction, ray.direction);
        const closest0 = point.clone();
        let closest1: Vector;
        if (parameter > 0) {
            closest1 = add(ray.origin, mul(parameter, ray.direction));
        }
        else {
            parameter = 0;
            closest1 = ray.origin.clone();
        }
        diff = sub(closest0, closest1);
        const sqrDistance = dot(diff, diff);

        return {
            distance: Math.sqrt(sqrDistance),
            sqrDistance,
            parameter,
            closest: [closest0, closest1]
        };
    }
}
