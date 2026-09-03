// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment3Plane3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the plane to be a 2-dimensional object embedded in 3D
// and the segment to be a 1-dimensional object.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has both
// a TIQuery and an FIQuery specialization, which become IntrSegment3Plane3TI
// and IntrSegment3Plane3FI. The upstream FIQuery derives from
// FIQuery<T, Line3, Plane3> and adds no result members, so the port aliases
// IntrLine3Plane3FIResult (the precedent set by IntrSegment2Triangle2.ts) and
// reuses the exported 'intrLine3Plane3FIDoQuery'. The protected
// 'FIQuery::DoQuery' becomes the exported module function
// 'intrSegment3Plane3FIDoQuery'.
//
// Note that the FIQuery reports the intersection parameter relative to the
// centered form of the segment, that is, the segment is
// center + t * direction for |t| <= extent with unit-length direction.

import { DistPointHyperplane } from './DistPointHyperplane';
import type { FIQuery } from './FIQuery';
import {
    defaultIntrLine3Plane3FIResult, intrLine3Plane3FIDoQuery
} from './IntrLine3Plane3';
import type { IntrLine3Plane3FIResult } from './IntrLine3Plane3';
import type { Plane3 } from './Hyperplane';
import type { Segment3 } from './Segment';
import type { TIQuery } from './TIQuery';
import { Vector, add, mul } from './Vector';

// The result of IntrSegment3Plane3TI.test.
export interface IntrSegment3Plane3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSegment3Plane3TIResult():
    IntrSegment3Plane3TIResult {
    return { intersect: false };
}

// The result of IntrSegment3Plane3FI.find. Upstream derives its Result from
// the line-plane Result and adds no members.
export type IntrSegment3Plane3FIResult = IntrLine3Plane3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment3Plane3FIResult():
    IntrSegment3Plane3FIResult {
    return defaultIntrLine3Plane3FIResult();
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
export function intrSegment3Plane3FIDoQuery(segOrigin: Vector,
    segDirection: Vector, segExtent: number, plane: Plane3,
    result: IntrSegment3Plane3FIResult): void {
    intrLine3Plane3FIDoQuery(segOrigin, segDirection, plane, result);

    if (result.intersect) {
        // The line intersects the plane in a point that might not be on the
        // segment.
        if (Math.abs(result.parameter) > segExtent) {
            result.intersect = false;
            result.numIntersections = 0;
        }
    }
}

// Test-intersection query for a segment and a plane in 3D.
export class IntrSegment3Plane3TI implements
    TIQuery<Segment3, Plane3, IntrSegment3Plane3TIResult> {

    test(segment: Segment3, plane: Plane3): IntrSegment3Plane3TIResult {
        const result = defaultIntrSegment3Plane3TIResult();

        // Compute the (signed) distance from the segment endpoints to the
        // plane.
        const vpQuery = new DistPointHyperplane();
        const sdistance0 = vpQuery.compute(segment.p[0], plane).signedDistance;
        if (sdistance0 === 0) {
            // Endpoint p[0] is on the plane.
            result.intersect = true;
            return result;
        }

        const sdistance1 = vpQuery.compute(segment.p[1], plane).signedDistance;
        if (sdistance1 === 0) {
            // Endpoint p[1] is on the plane.
            result.intersect = true;
            return result;
        }

        // Test whether the segment transversely intersects the plane.
        result.intersect = (sdistance0 * sdistance1 < 0);
        return result;
    }
}

// Find-intersection query for a segment and a plane in 3D.
export class IntrSegment3Plane3FI implements
    FIQuery<Segment3, Plane3, IntrSegment3Plane3FIResult> {

    find(segment: Segment3, plane: Plane3): IntrSegment3Plane3FIResult {
        const { center: segOrigin, direction: segDirection, extent: segExtent } =
            segment.getCenteredForm();

        const result = defaultIntrSegment3Plane3FIResult();
        intrSegment3Plane3FIDoQuery(segOrigin, segDirection, segExtent, plane,
            result);
        if (result.intersect) {
            result.point = add(segOrigin, mul(result.parameter, segDirection));
        }
        return result;
    }
}
