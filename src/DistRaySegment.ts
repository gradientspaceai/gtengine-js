// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRaySegment.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a ray and a segment in nD.
//
// The ray is P[0] + s[0] * D[0] for s[0] >= 0. D[0] is not required to be
// unit length.
//
// The segment is Q[0] + s[1] * (Q[1] - Q[0]) for 0 <= s[1] <= 1. The
// direction D = Q[1] - Q[0] is generally not unit length.
//
// The closest point on the ray is stored in closest[0] with parameter[0]
// storing s[0]. The closest point on the segment is stored in closest[1]
// with parameter[1] storing s[1]. When there are infinitely many choices for
// the pair of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Ray<N,T>, Segment<N,T>>' becomes the
// class DistRaySegment with the result type DistRaySegmentResult.

import type { DCPQuery } from './DCPQuery.js';
import type { Ray } from './Ray.js';
import type { Segment } from './Segment.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

export interface DistRaySegmentResult {
    distance: number;
    sqrDistance: number;

    // parameter[0] is the nonnegative ray parameter, parameter[1] is the
    // segment parameter in [0,1].
    parameter: [number, number];

    // closest[0] is on the ray, closest[1] is on the segment.
    closest: [Vector, Vector];
}

export class DistRaySegment
    implements DCPQuery<Ray, Segment, DistRaySegmentResult> {
    compute(ray: Ray, segment: Segment): DistRaySegmentResult {
        const segDirection = sub(segment.p[1], segment.p[0]);
        let diff = sub(ray.origin, segment.p[0]);
        const a00 = dot(ray.direction, ray.direction);
        const a01 = -dot(ray.direction, segDirection);
        const a11 = dot(segDirection, segDirection);
        const b0 = dot(ray.direction, diff);
        const b1 = -dot(segDirection, diff);
        const det = Math.max(a00 * a11 - a01 * a01, 0);
        let s0: number;
        let s1: number;

        if (det > 0) {
            // The ray and segment are not parallel.
            s0 = a01 * b1 - a11 * b0;
            s1 = a01 * b0 - a00 * b1;

            if (s0 >= 0) {
                if (s1 >= 0) {
                    if (s1 <= det) {
                        // region 0: the minimum occurs at interior points of
                        // the ray and the segment.
                        s0 /= det;
                        s1 /= det;
                    }
                    else {
                        // region 1: the endpoint Q1 of the segment and an
                        // interior point of the ray are closest.
                        s0 = -(a01 + b0) / a00;
                        s1 = 1;
                    }
                }
                else {
                    // region 5: the endpoint Q0 of the segment and an
                    // interior point of the ray are closest.
                    s0 = -b0 / a00;
                    s1 = 0;
                }
            }
            else {
                // s0 < 0
                if (s1 <= 0) {
                    // region 4
                    s0 = -b0;
                    if (s0 > 0) {
                        s0 /= a00;
                        s1 = 0;
                    }
                    else {
                        s0 = 0;
                        s1 = -b1;
                        if (s1 < 0) {
                            s1 = 0;
                        }
                        else if (s1 > a11) {
                            s1 = 1;
                        }
                        else {
                            s1 /= a11;
                        }
                    }
                }
                else if (s1 <= det) {
                    // region 3
                    s0 = 0;
                    s1 = -b1;
                    if (s1 < 0) {
                        s1 = 0;
                    }
                    else if (s1 > a11) {
                        s1 = 1;
                    }
                    else {
                        s1 /= a11;
                    }
                }
                else {
                    // region 2
                    s0 = -(a01 + b0);
                    if (s0 > 0) {
                        s0 /= a00;
                        s1 = 1;
                    }
                    else {
                        s0 = 0;
                        s1 = -b1;
                        if (s1 < 0) {
                            s1 = 0;
                        }
                        else if (s1 > a11) {
                            s1 = 1;
                        }
                        else {
                            s1 /= a11;
                        }
                    }
                }
            }
        }
        else {
            // The ray and segment are parallel. Choose the segment endpoint
            // whose projection onto the ray line is farthest along the ray
            // direction; the ray parameter is that projection.
            if (a01 > 0) {
                // Opposite direction vectors; Q0 is the farthest endpoint.
                s0 = -b0 / a00;
                s1 = 0;
            }
            else {
                // Same direction vectors; Q1 is the farthest endpoint. This
                // also covers a degenerate segment (a01 = a11 = 0), for
                // which Q0 and Q1 are the same point.
                s0 = -(a01 + b0) / a00;
                s1 = 1;
            }
        }

        // UPSTREAM BUG FIX (DistRaySegment.h). Several branches minimize the
        // squared distance on the segment edge s1 = 0 or s1 = 1 with s0
        // unconstrained, and none of them clamps the resulting s0 to the ray
        // domain s0 >= 0:
        //
        //   * the parallel branch (both cases), when the whole segment
        //     projects behind the ray origin;
        //   * region 1 (s1 clamped to 1) and region 5 (s1 clamped to 0) of
        //     the nonparallel branch, where s0 = -(a01+b0)/a00 and
        //     s0 = -b0/a00 can be negative even though the unconstrained
        //     line/line parameter is nonnegative. For example, ray origin
        //     (-1.985,-2.485,-1.567) with direction
        //     (-2.705,-0.448,-1.664) and segment endpoints
        //     (1.403,-2.981,2.256) and (-1.674,-2.963,0.385) reaches
        //     region 1 with s0 = -0.377.
        //
        // Upstream then reports a closest "ray" point that is not on the ray
        // and a distance smaller than the true one. The sibling query
        // DistRayRay.h does clamp its analogous branches. Since the squared
        // distance is convex on the domain {s0 >= 0, 0 <= s1 <= 1}, a
        // minimizer with s0 < 0 means the true minimum lies on the face
        // s0 = 0, where the optimal segment parameter is -b1/a11 clamped to
        // [0,1]. (For a degenerate segment, a11 = 0 and the segment point is
        // Q0 = Q1.)
        if (s0 < 0) {
            s0 = 0;
            if (a11 > 0) {
                s1 = -b1 / a11;
                if (s1 < 0) {
                    s1 = 0;
                }
                else if (s1 > 1) {
                    s1 = 1;
                }
            }
            else {
                s1 = 0;
            }
        }

        const closest0 = add(ray.origin, mul(s0, ray.direction));
        const closest1 = add(segment.p[0], mul(s1, segDirection));
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
