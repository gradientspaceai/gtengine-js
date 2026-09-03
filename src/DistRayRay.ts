// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRayRay.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between two rays in nD.
//
// The rays are P[i] + s[i] * D[i] for s[i] >= 0, where D[i] is not required
// to be unit length.
//
// The closest point on ray[i] is stored in closest[i] with parameter[i]
// storing s[i]. When there are infinitely many choices for the pair of
// closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Ray<N,T>, Ray<N,T>>' becomes the
// class DistRayRay with the result type DistRayRayResult.

import type { DCPQuery } from './DCPQuery.js';
import type { Ray } from './Ray.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

export interface DistRayRayResult {
    distance: number;
    sqrDistance: number;

    // The nonnegative ray parameters s[0] and s[1] of the closest points.
    parameter: [number, number];

    // closest[0] is on ray0, closest[1] is on ray1.
    closest: [Vector, Vector];
}

export class DistRayRay implements DCPQuery<Ray, Ray, DistRayRayResult> {
    compute(ray0: Ray, ray1: Ray): DistRayRayResult {
        let diff = sub(ray0.origin, ray1.origin);
        const a00 = dot(ray0.direction, ray0.direction);
        const a01 = -dot(ray0.direction, ray1.direction);
        const a11 = dot(ray1.direction, ray1.direction);
        const b0 = dot(ray0.direction, diff);
        const det = Math.max(a00 * a11 - a01 * a01, 0);
        let s0: number;
        let s1: number;

        if (det > 0) {
            // The rays are not parallel.
            const b1 = -dot(ray1.direction, diff);
            s0 = a01 * b1 - a11 * b0;
            s1 = a01 * b0 - a00 * b1;

            if (s0 >= 0) {
                if (s1 >= 0) {
                    // region 0 (interior): the minimum occurs at two
                    // interior points of the rays.
                    s0 /= det;
                    s1 /= det;
                }
                else {
                    // region 3 (side)
                    if (b0 >= 0) {
                        s0 = 0;
                    }
                    else {
                        s0 = -b0 / a00;
                    }
                    s1 = 0;
                }
            }
            else {
                if (s1 >= 0) {
                    // region 1 (side)
                    s0 = 0;
                    if (b1 >= 0) {
                        s1 = 0;
                    }
                    else {
                        s1 = -b1 / a11;
                    }
                }
                else {
                    // region 2 (corner)
                    if (b0 < 0) {
                        s0 = -b0 / a00;
                        s1 = 0;
                    }
                    else {
                        s0 = 0;
                        if (b1 >= 0) {
                            s1 = 0;
                        }
                        else {
                            s1 = -b1 / a11;
                        }
                    }
                }
            }
        }
        else {
            // The rays are parallel.
            if (a01 > 0) {
                // Opposite direction vectors.
                s1 = 0;
                if (b0 >= 0) {
                    s0 = 0;
                }
                else {
                    s0 = -b0 / a00;
                }
            }
            else {
                // Same direction vectors.
                if (b0 >= 0) {
                    const b1 = -dot(ray1.direction, diff);
                    s0 = 0;
                    s1 = -b1 / a11;
                }
                else {
                    s0 = -b0 / a00;
                    s1 = 0;
                }
            }
        }

        const closest0 = add(ray0.origin, mul(s0, ray0.direction));
        const closest1 = add(ray1.origin, mul(s1, ray1.direction));
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
