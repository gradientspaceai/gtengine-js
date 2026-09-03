// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistSegment3Circle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The 3D segment-circle distance algorithm is described in
// https://www.geometrictools.com/Documentation/DistanceToCircle3.pdf
// The notation used in the code matches that of the document. The circle has
// center C and the plane of the circle has unit-length normal N. The segment
// has endpoints P0 and P1. The parameterization is P(t) = P0 + t * (P1 - P0)
// = P0 + t * M, where M is generally not a unit-length vector.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Segment3<T>, Circle3<T>>' becomes the
// class DistSegment3Circle3. As upstream does, the result type is the
// line-circle result type (with 'linearClosest' naming the segment points),
// re-exported here as the alias DistSegment3Circle3Result. The upstream
// friend access to 'DCPQuery<T,Line3,Circle3>::Execute' and its private
// 'Critical' struct are the exported 'distLine3Circle3Execute' and
// 'DistLine3Circle3Critical' from DistLine3Circle3.ts. The private helpers
// Execute, HasOneCriticalPoint, HasTwoCriticalPoints,
// SegmentEndpointClosest and SelectClosestPoint become module-private
// functions.

import type { Circle3 } from './Circle3.js';
import type { DCPQuery } from './DCPQuery.js';
import { distLine3Circle3Execute } from './DistLine3Circle3.js';
import type { DistLine3Circle3Critical, DistLine3Circle3Result }
    from './DistLine3Circle3.js';
import { DistPoint3Circle3 } from './DistPoint3Circle3.js';
import { Line } from './Line.js';
import type { Segment3 } from './Segment.js';
import { Vector, sub } from './Vector.js';

// Upstream reuses the line-circle result type ('using Result = typename
// LCQuery::Result').
export type DistSegment3Circle3Result = DistLine3Circle3Result;

function segmentEndpointClosest(segmentEndpoint: Vector, circle: Circle3,
    result: DistSegment3Circle3Result): void {
    const pcResult = new DistPoint3Circle3().compute(segmentEndpoint, circle);
    result.numClosestPairs = 1;
    result.linearClosest[0] = segmentEndpoint.clone();
    result.linearClosest[1] = new Vector(3);
    result.circularClosest[0] = pcResult.closest[1];
    result.circularClosest[1] = new Vector(3);
    result.distance = pcResult.distance;
    result.sqrDistance = result.distance * result.distance;
}

function selectClosestPoint(point0: Vector, point1: Vector, circle: Circle3,
    result: DistSegment3Circle3Result): void {
    const pcQuery = new DistPoint3Circle3();
    const pcResult0 = pcQuery.compute(point0, circle);
    const pcResult1 = pcQuery.compute(point1, circle);
    if (pcResult0.distance < pcResult1.distance) {
        result.numClosestPairs = 1;
        result.linearClosest[0] = point0.clone();
        result.linearClosest[1] = new Vector(3);
        result.circularClosest[0] = pcResult0.closest[1];
        result.circularClosest[1] = new Vector(3);
        result.distance = pcResult0.distance;
        result.sqrDistance = result.distance * result.distance;
    }
    else if (pcResult0.distance > pcResult1.distance) {
        result.numClosestPairs = 1;
        result.linearClosest[0] = point1.clone();
        result.linearClosest[1] = new Vector(3);
        result.circularClosest[0] = pcResult1.closest[1];
        result.circularClosest[1] = new Vector(3);
        result.distance = pcResult1.distance;
        result.sqrDistance = result.distance * result.distance;
    }
    else {
        // pcResult0.distance = pcResult1.distance
        result.numClosestPairs = 2;
        result.linearClosest[0] = point0.clone();
        result.linearClosest[1] = point1.clone();
        result.circularClosest[0] = pcResult0.closest[1];
        result.circularClosest[1] = pcResult1.closest[1];
        result.distance = pcResult0.distance;
        result.sqrDistance = result.distance * result.distance;
    }
}

function hasOneCriticalPoint(segment: Segment3, circle: Circle3,
    critical: DistLine3Circle3Critical,
    result: DistSegment3Circle3Result): void {
    const t0 = critical.parameter[0];

    if (t0 >= 1) {
        // The critical point is not on the segment except possibly the first
        // critical point being the right endpoint of the segment. The right
        // endpoint is the segment point closest to the circle. See the left
        // red segment of the one-critical-point graph of figure 8 in the PDF.
        segmentEndpointClosest(segment.p[1], circle, result);
        return;
    }

    if (t0 <= 0) {
        // The critical points are not on the segment except possibly the
        // critical point being the left endpoint of the segment. The left
        // endpoint is the segment point closest to the circle. See the right
        // red segment of the one-critical-point graph of figure 8 in the PDF.
        segmentEndpointClosest(segment.p[0], circle, result);
        return;
    }

    // At this time, 0 < t0 < 1. The closest line-circle pair is the closest
    // segment-circle pair. The output does not need to be modified. See the
    // green segment of the one-critical-point graph of figure 8 in the PDF.
}

function hasTwoCriticalPoints(segment: Segment3, circle: Circle3,
    critical: DistLine3Circle3Critical,
    result: DistSegment3Circle3Result): void {
    const t0 = critical.parameter[0];
    const t1 = critical.parameter[1];

    if (t0 >= 1) {
        // The critical points are not on the segment except possibly the
        // first critical point being the right endpoint of the segment. The
        // right endpoint is the segment point closest to the circle. See the
        // left red segment of the two-point-critical graphs of figure 8 in
        // the PDF.
        segmentEndpointClosest(segment.p[1], circle, result);
        return;
    }

    if (t1 <= 0) {
        // The critical points are not on the segment except possibly the
        // second critical point being the left endpoint of the segment. The
        // left endpoint is the segment point closest to the circle. See the
        // right red segment of the two-point-critical graphs of figure 8 in
        // the PDF.
        segmentEndpointClosest(segment.p[0], circle, result);
        return;
    }

    // At this time, t0 < 1 and t1 > 0.
    if (0 <= t0 && t1 <= 1) {
        // At this time, 0 <= t0 < t1 <= 1. The critical points are on the
        // segment, so the closest segment-circle pairs are the closest
        // line-circle pairs. The output does not need to be modified. See the
        // green segment of the two-critical-point graphs of figure 8 in the
        // PDF.
        return;
    }

    // At this time, t0 < 0 or t1 > 1. At most one critical point is on the
    // segment.
    if (t0 < 0) {
        if (t1 >= 1) {
            // At this time, t0 < 0 < 1 <= t1. The critical points are not on
            // the segment except possibly the second critical point is the
            // right endpoint. See the orange segment of the
            // two-critical-point graphs of figure 8 in the PDF.
            selectClosestPoint(segment.p[0], segment.p[1], circle, result);
        }
        else {
            // At this time, t0 < 0 < t1 < 1. The critical point at t1 is on
            // the segment but is not an endpoint. See the purple segment of
            // the two-critical-point graphs of figure 8 in the PDF.
            selectClosestPoint(segment.p[0], critical.linearPoint[1], circle,
                result);
        }
    }
    else {
        // t1 > 1
        if (t0 <= 0) {
            // At this time, t0 <= 0 < 1 < t1. The critical points are not on
            // the segment except possibly the first critical point is the
            // left endpoint. See the orange segment of the
            // two-critical-point graphs of figure 8 in the PDF.
            selectClosestPoint(segment.p[0], segment.p[1], circle, result);
        }
        else {
            // At this time, 0 < t0 < 1 < t1. The critical point at t0 is on
            // the segment but is not an endpoint. See the gold segment of the
            // two-critical-point graphs of figure 8 in the PDF.
            selectClosestPoint(segment.p[1], critical.linearPoint[0], circle,
                result);
        }
    }
}

export class DistSegment3Circle3
    implements DCPQuery<Segment3, Circle3, DistSegment3Circle3Result> {
    compute(segment: Segment3, circle: Circle3): DistSegment3Circle3Result {
        // Compute the line points closest to the circle. The line is
        // L(t) = P + t * D for any real-valued t. The segment restricts
        // 0 <= t <= 1 and has endpoints P0 = L(0) and P1 = L(1) with
        // D = P1 - P0.
        const line = Line.fromOriginDirection(segment.p[0],
            sub(segment.p[1], segment.p[0]));
        const { result, critical } = distLine3Circle3Execute(line, circle);

        // Clamp the query output to the segment domain.
        if (critical.numPoints === 1) {
            hasOneCriticalPoint(segment, circle, critical, result);
        }
        else {
            hasTwoCriticalPoints(segment, circle, critical, result);
        }

        return result;
    }
}
