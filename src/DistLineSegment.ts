// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLineSegment.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a line and a segment in nD.
//
// The line is P + s[0] * D, where D is not required to be unit length.
//
// The segment is Q[0] + s[1] * (Q[1] - Q[0]) for 0 <= s[1] <= 1. The
// direction Q[1] - Q[0] is generally not unit length.
//
// The closest point on the line is stored in closest[0] with parameter[0]
// storing s[0]. The closest point on the segment is stored in closest[1] with
// parameter[1] storing s[1]. When there are infinitely many choices for the
// pair of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Line<N,T>, Segment<N,T>>' becomes the
// class DistLineSegment with the result type DistLineSegmentResult. The
// upstream file comment swaps the roles of closest[0] and closest[1]
// (it describes closest[0] as the segment point); the code stores the line
// point in closest[0], which is what the port documents.

import type { DCPQuery } from './DCPQuery';
import type { Line } from './Line';
import type { Segment } from './Segment';
import { Vector, add, dot, mul, sub } from './Vector';

export interface DistLineSegmentResult {
    distance: number;
    sqrDistance: number;

    // parameter[0] is the line parameter, parameter[1] is the segment
    // parameter in [0,1].
    parameter: [number, number];

    // closest[0] is on the line, closest[1] is on the segment.
    closest: [Vector, Vector];
}

export class DistLineSegment
    implements DCPQuery<Line, Segment, DistLineSegmentResult> {
    compute(line: Line, segment: Segment): DistLineSegmentResult {
        const segDirection = sub(segment.p[1], segment.p[0]);
        let diff = sub(line.origin, segment.p[0]);
        const a00 = dot(line.direction, line.direction);
        const a01 = -dot(line.direction, segDirection);
        const a11 = dot(segDirection, segDirection);
        const b0 = dot(line.direction, diff);
        const det = Math.max(a00 * a11 - a01 * a01, 0);
        let s0: number;
        let s1: number;

        if (det > 0) {
            // The line and segment are not parallel.
            const b1 = -dot(segDirection, diff);
            s1 = a01 * b0 - a00 * b1;

            if (s1 >= 0) {
                if (s1 <= det) {
                    // Two interior points are closest, one on the line and
                    // one on the segment.
                    s0 = (a01 * b1 - a11 * b0) / det;
                    s1 /= det;
                }
                else {
                    // The endpoint Q1 of the segment and an interior point of
                    // the line are closest.
                    s0 = -(a01 + b0) / a00;
                    s1 = 1;
                }
            }
            else {
                // The endpoint Q0 of the segment and an interior point of the
                // line are closest.
                s0 = -b0 / a00;
                s1 = 0;
            }
        }
        else {
            // The line and segment are parallel. Select the pair of closest
            // points where the closest segment point is the endpoint Q0.
            s0 = -b0 / a00;
            s1 = 0;
        }

        const closest0 = add(line.origin, mul(s0, line.direction));
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
