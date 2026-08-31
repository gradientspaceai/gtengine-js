// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPointSegment.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a point and a segment in nD.
//
// The segment is P0 + t * (P1 - P0) for 0 <= t <= 1. The direction
// D = P1 - P0 is generally not unit length.
//
// The input point is stored in closest[0]. The closest point on the segment
// is stored in closest[1]. When there are infinitely many choices for the
// pair of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector<N,T>, Segment<N,T>>' becomes
// the class DistPointSegment with the result type DistPointSegmentResult.

import type { DCPQuery } from './DCPQuery';
import type { Segment } from './Segment';
import { Vector, add, dot, mul, sub } from './Vector';

export interface DistPointSegmentResult {
    distance: number;
    sqrDistance: number;

    // The segment parameter t in [0,1] of the closest segment point.
    parameter: number;

    // closest[0] is the input point, closest[1] is the closest segment point.
    closest: [Vector, Vector];
}

export class DistPointSegment
    implements DCPQuery<Vector, Segment, DistPointSegmentResult> {
    compute(point: Vector, segment: Segment): DistPointSegmentResult {
        // The direction vector is not unit length. The normalization is
        // deferred until it is needed.
        const direction = sub(segment.p[1], segment.p[0]);
        let diff = sub(point, segment.p[1]);
        let t = dot(direction, diff);
        let parameter: number;
        let closest1: Vector;
        if (t >= 0) {
            parameter = 1;
            closest1 = segment.p[1].clone();
        }
        else {
            diff = sub(point, segment.p[0]);
            t = dot(direction, diff);
            if (t <= 0) {
                parameter = 0;
                closest1 = segment.p[0].clone();
            }
            else {
                const sqrLength = dot(direction, direction);
                if (sqrLength > 0) {
                    t /= sqrLength;
                    parameter = t;
                    closest1 = add(segment.p[0], mul(t, direction));
                }
                else {
                    parameter = 0;
                    closest1 = segment.p[0].clone();
                }
            }
        }
        const closest0 = point.clone();

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
