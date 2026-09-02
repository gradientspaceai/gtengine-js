// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay2SegmentMesh2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query uses a line-segment intersection test with the segments of the
// mesh. The resulting set of intersection points is trimmed by discarding
// those for which the line parameters are negative.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// an FIQuery specialization, which becomes IntrRay2SegmentMesh2FI. The nested
// Intersection struct becomes the exported interface
// IntrRay2SegmentMesh2Intersection (as in IntrLine2SegmentMesh2.ts). The
// 'std::remove_if' followed by a copy of the surviving range is an
// order-preserving filter.

import type { FIQuery } from './FIQuery';
import { IntrLine2SegmentMesh2FI } from './IntrLine2SegmentMesh2';
import { Line } from './Line';
import type { Ray2 } from './Ray';
import type { SegmentMesh2 } from './SegmentMesh';
import type { Vector } from './Vector';

// One point at which the ray intersects a mesh segment.
export interface IntrRay2SegmentMesh2Intersection {
    // The indices of the mesh vertices that are the endpoints of the mesh
    // segment containing the intersection point.
    indexPair: [number, number];

    // The intersection point is ray.origin + rayParameter * ray.direction.
    rayParameter: number;

    // The intersection point is (1-t) * V0 + t * V1, where t is the
    // meshSegmentParameter and V0, V1 are the mesh segment endpoints.
    meshSegmentParameter: number;

    // The point of intersection.
    point: Vector;
}

// The result of IntrRay2SegmentMesh2FI.find.
export interface IntrRay2SegmentMesh2FIResult {
    intersections: IntrRay2SegmentMesh2Intersection[];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrRay2SegmentMesh2FIResult():
    IntrRay2SegmentMesh2FIResult {
    return { intersections: [] };
}

// Find-intersection query for a ray and a segment mesh in 2D.
export class IntrRay2SegmentMesh2FI implements
    FIQuery<Ray2, SegmentMesh2, IntrRay2SegmentMesh2FIResult> {

    find(ray: Ray2, mesh: SegmentMesh2): IntrRay2SegmentMesh2FIResult {
        const result = defaultIntrRay2SegmentMesh2FIResult();

        // Execute the line-mesh query and then remove intersections for which
        // the line parameter is negative. The remaining intersections are on
        // the ray.
        const lsQuery = new IntrLine2SegmentMesh2FI();
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const lsResult = lsQuery.find(line, mesh);
        for (const object of lsResult.intersections) {
            if (object.lineParameter >= 0) {
                result.intersections.push({
                    indexPair: [object.indexPair[0], object.indexPair[1]],
                    rayParameter: object.lineParameter,
                    meshSegmentParameter: object.meshSegmentParameter,
                    point: object.point.clone()
                });
            }
        }

        return result;
    }
}
