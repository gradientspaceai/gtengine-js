// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine2SegmentMesh2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query performs an exhaustive search of the segments and finds
// line-segment intersections.
//
// Upstream TODOs (not ported): multithreading for large meshes, and a
// preprocessed bounding-region tree to reduce the O(n) line-segment
// intersection tests to O(log n).
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// an FIQuery specialization, so the port has only IntrLine2SegmentMesh2FI.
// The nested 'Intersection' struct becomes the exported interface
// IntrLine2SegmentMesh2Intersection.

import type { FIQuery } from './FIQuery.js';
import { IntrLine2Segment2FI } from './IntrLine2Segment2.js';
import type { Line2 } from './Line.js';
import { Segment } from './Segment.js';
import type { SegmentMesh2 } from './SegmentMesh.js';
import type { Vector } from './Vector.js';

// One line-mesh intersection record.
export interface IntrLine2SegmentMesh2Intersection {
    // The indices of the mesh vertices of the intersected mesh segment.
    indexPair: [number, number];

    // The line parameter t for the intersection point P = origin + t*direction.
    lineParameter: number;

    // The mesh-segment parameter t for P = (1-t)*p[0] + t*p[1].
    meshSegmentParameter: number;

    // The intersection point.
    point: Vector;
}

// The result of IntrLine2SegmentMesh2FI.find.
export interface IntrLine2SegmentMesh2FIResult {
    intersections: IntrLine2SegmentMesh2Intersection[];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrLine2SegmentMesh2FIResult():
    IntrLine2SegmentMesh2FIResult {
    return { intersections: [] };
}

// Find-intersection query for a line and a segment mesh in 2D.
export class IntrLine2SegmentMesh2FI implements
    FIQuery<Line2, SegmentMesh2, IntrLine2SegmentMesh2FIResult> {

    find(line: Line2, mesh: SegmentMesh2): IntrLine2SegmentMesh2FIResult {
        const result = defaultIntrLine2SegmentMesh2FIResult();

        const lsQuery = new IntrLine2Segment2FI();

        const vertices = mesh.getVertices();
        const indices = mesh.getIndices();
        for (let i = 0; i < indices.length; ++i) {
            const segment = Segment.fromEndpoints(
                vertices[indices[i][0]], vertices[indices[i][1]]);
            const lsResult = lsQuery.find(line, segment);
            if (lsResult.intersect) {
                if (lsResult.numIntersections === 1) {
                    // The line and segment intersect in a unique point.
                    result.intersections.push({
                        indexPair: [indices[i][0], indices[i][1]],
                        lineParameter: lsResult.lineParameter[0],
                        meshSegmentParameter: lsResult.segmentParameter[0],
                        point: lsResult.point.clone()
                    });
                }
                else {
                    // The line and segment are coincident. Report both
                    // segment endpoints as intersections.
                    for (let j = 0; j < 2; ++j) {
                        result.intersections.push({
                            indexPair: [indices[i][0], indices[i][1]],
                            lineParameter: lsResult.lineParameter[j],
                            meshSegmentParameter: lsResult.segmentParameter[j],
                            point: segment.p[j].clone()
                        });
                    }
                }
            }
        }

        // Sort the intersection points by line parameter. This makes it
        // easier to implement the ray-mesh and segment-mesh queries than by
        // using the lower-level Ray2-Segment2 and Segment2-Segment2
        // intersection queries.
        result.intersections.sort(
            (object0, object1) => object0.lineParameter - object1.lineParameter);

        return result;
    }
}
