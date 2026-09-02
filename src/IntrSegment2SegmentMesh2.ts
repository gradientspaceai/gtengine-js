// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment2SegmentMesh2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query uses a line-segment intersection test with the segments of the
// mesh. The resulting set of intersection points is trimmed by discarding
// those for which the line parameters are outside [0,1].
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// an FIQuery specialization, which becomes IntrSegment2SegmentMesh2FI. The
// nested Intersection struct becomes the exported interface
// IntrSegment2SegmentMesh2Intersection (as in IntrLine2SegmentMesh2.ts). The
// 'std::remove_if' followed by a copy of the surviving range is an
// order-preserving filter. Note that the line used by the query is
// (segment.p[0], segment.p[1] - segment.p[0]), whose direction is not unit
// length, so the reported segmentParameter is in [0,1] rather than in
// [-extent, extent].

import type { FIQuery } from './FIQuery';
import { IntrLine2SegmentMesh2FI } from './IntrLine2SegmentMesh2';
import { Line } from './Line';
import type { Segment2 } from './Segment';
import type { SegmentMesh2 } from './SegmentMesh';
import { sub } from './Vector';
import type { Vector } from './Vector';

// One point at which the segment intersects a mesh segment.
export interface IntrSegment2SegmentMesh2Intersection {
    // The indices of the mesh vertices that are the endpoints of the mesh
    // segment containing the intersection point.
    indexPair: [number, number];

    // The intersection point is (1-t) * segment.p[0] + t * segment.p[1],
    // where t is the segmentParameter.
    segmentParameter: number;

    // The intersection point is (1-t) * V0 + t * V1, where t is the
    // meshSegmentParameter and V0, V1 are the mesh segment endpoints.
    meshSegmentParameter: number;

    // The point of intersection.
    point: Vector;
}

// The result of IntrSegment2SegmentMesh2FI.find.
export interface IntrSegment2SegmentMesh2FIResult {
    intersections: IntrSegment2SegmentMesh2Intersection[];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment2SegmentMesh2FIResult():
    IntrSegment2SegmentMesh2FIResult {
    return { intersections: [] };
}

// Find-intersection query for a segment and a segment mesh in 2D.
export class IntrSegment2SegmentMesh2FI implements
    FIQuery<Segment2, SegmentMesh2, IntrSegment2SegmentMesh2FIResult> {

    find(segment: Segment2, mesh: SegmentMesh2):
        IntrSegment2SegmentMesh2FIResult {
        const result = defaultIntrSegment2SegmentMesh2FIResult();

        const lsQuery = new IntrLine2SegmentMesh2FI();
        const line = Line.fromOriginDirection(segment.p[0],
            sub(segment.p[1], segment.p[0]));
        const lsResult = lsQuery.find(line, mesh);
        for (const object of lsResult.intersections) {
            if (object.lineParameter >= 0 && object.lineParameter <= 1) {
                result.intersections.push({
                    indexPair: [object.indexPair[0], object.indexPair[1]],
                    segmentParameter: object.lineParameter,
                    meshSegmentParameter: object.meshSegmentParameter,
                    point: object.point.clone()
                });
            }
        }

        return result;
    }
}
