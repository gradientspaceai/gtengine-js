// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLineRay.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a line and a ray in nD.
//
// The line is P[0] + s[0] * D[0] and the ray is P[1] + s[1] * D[1] for
// s[1] >= 0. The D[i] are not required to be unit length.
//
// The closest point on the line is stored in closest[0] with parameter[0]
// storing s[0]. The closest point on the ray is stored in closest[1] with
// parameter[1] storing s[1]. When there are infinitely many choices for the
// pair of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Line<N,T>, Ray<N,T>>' becomes the
// class DistLineRay with the result type DistLineRayResult.

import type { DCPQuery } from './DCPQuery';
import type { Line } from './Line';
import type { Ray } from './Ray';
import { Vector, add, dot, mul, sub } from './Vector';

export interface DistLineRayResult {
    distance: number;
    sqrDistance: number;

    // parameter[0] is the line parameter, parameter[1] is the ray parameter.
    parameter: [number, number];

    // closest[0] is on the line, closest[1] is on the ray.
    closest: [Vector, Vector];
}

export class DistLineRay implements DCPQuery<Line, Ray, DistLineRayResult> {
    compute(line: Line, ray: Ray): DistLineRayResult {
        let diff = sub(line.origin, ray.origin);
        const a00 = dot(line.direction, line.direction);
        const a01 = -dot(line.direction, ray.direction);
        const a11 = dot(ray.direction, ray.direction);
        const b0 = dot(line.direction, diff);
        const det = Math.max(a00 * a11 - a01 * a01, 0);
        let s0: number;
        let s1: number;

        if (det > 0) {
            // The line and ray are not parallel.
            const b1 = -dot(ray.direction, diff);
            s1 = a01 * b0 - a00 * b1;

            if (s1 >= 0) {
                // Two interior points are closest, one on the line and one on
                // the ray.
                s0 = (a01 * b1 - a11 * b0) / det;
                s1 /= det;
            }
            else {
                // The origin of the ray is the closest ray point.
                s0 = -b0 / a00;
                s1 = 0;
            }
        }
        else {
            // The line and ray are parallel. Select the pair of closest
            // points where the closest ray point is the ray origin.
            s0 = -b0 / a00;
            s1 = 0;
        }

        const closest0 = add(line.origin, mul(s0, line.direction));
        const closest1 = add(ray.origin, mul(s1, ray.direction));
        diff = sub(closest0, closest1);
        const sqrDistance = dot(diff, diff);

        return {
            distance: Math.sqrt(sqrDistance),
            sqrDistance,
            parameter: [s0, s1],
            closest: [closest0, closest1]
        };
    }
}
