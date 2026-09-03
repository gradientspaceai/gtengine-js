// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment3Capsule3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the capsule to be a solid.
//
// The test-intersection queries are based on distance computations.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// the FIQuery from the Line3-vs-Capsule3 query to reuse the protected DoQuery
// helper, which the port exports as the module function
// 'intrLine3Capsule3FIDoQuery'. The segment-specific DoQuery helper is
// exported as the module function 'intrSegment3Capsule3FIDoQuery'. The
// reported parameters are relative to the centered form of the segment,
// C + t * D with |t| <= e, as upstream reports them.

import type { Capsule3 } from './Capsule.js';
import { DistSegmentSegment } from './DistSegmentSegment.js';
import type { FIQuery } from './FIQuery.js';
import { IntrIntervalsFI } from './IntrIntervals.js';
import {
    intrLine3Capsule3FIDoQuery,
    defaultIntrLine3Capsule3FIResult
} from './IntrLine3Capsule3.js';
import type { IntrLine3Capsule3FIResult } from './IntrLine3Capsule3.js';
import type { Segment3 } from './Segment.js';
import { Vector, add, mul } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrSegment3Capsule3TI.test.
export interface IntrSegment3Capsule3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrSegment3Capsule3TIResult {
    return { intersect: false };
}

// The upstream derived FIQuery::Result adds no members.
export type IntrSegment3Capsule3FIResult = IntrLine3Capsule3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment3Capsule3FIResult():
    IntrSegment3Capsule3FIResult {
    return defaultIntrLine3Capsule3FIResult();
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrSegment3Capsule3FIDoQuery(segOrigin: Vector,
    segDirection: Vector, segExtent: number, capsule: Capsule3,
    result: IntrSegment3Capsule3FIResult): void {
    intrLine3Capsule3FIDoQuery(segOrigin, segDirection, capsule, result);

    if (result.intersect) {
        // The line containing the segment intersects the capsule; the
        // t-interval is [t0,t1]. The segment intersects the capsule as long
        // as [t0,t1] overlaps the segment t-interval [-segExtent,+segExtent].
        const iiQuery = new IntrIntervalsFI();
        const segInterval: [number, number] = [-segExtent, segExtent];
        const iiResult = iiQuery.find(result.parameter, segInterval);
        if (iiResult.intersect) {
            result.numIntersections = iiResult.numIntersections;
            result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
        }
        else {
            // The segment does not intersect the capsule.
            const empty = defaultIntrSegment3Capsule3FIResult();
            result.intersect = empty.intersect;
            result.numIntersections = empty.numIntersections;
            result.parameter = empty.parameter;
            result.point = empty.point;
        }
    }
}

// Test-intersection query for a segment and a solid capsule in 3D.
export class IntrSegment3Capsule3TI implements
    TIQuery<Segment3, Capsule3, IntrSegment3Capsule3TIResult> {

    test(segment: Segment3, capsule: Capsule3):
        IntrSegment3Capsule3TIResult {
        const result = defaultTIResult();
        const ssQuery = new DistSegmentSegment();
        const ssResult = ssQuery.compute(segment, capsule.segment);
        result.intersect = (ssResult.distance <= capsule.radius);
        return result;
    }
}

// Find-intersection query for a segment and a solid capsule in 3D.
export class IntrSegment3Capsule3FI implements
    FIQuery<Segment3, Capsule3, IntrSegment3Capsule3FIResult> {

    find(segment: Segment3, capsule: Capsule3):
        IntrSegment3Capsule3FIResult {
        const { center: segOrigin, direction: segDirection,
            extent: segExtent } = segment.getCenteredForm();

        const result = defaultIntrSegment3Capsule3FIResult();
        intrSegment3Capsule3FIDoQuery(segOrigin, segDirection, segExtent,
            capsule, result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(segOrigin,
                    mul(result.parameter[i], segDirection));
            }
        }
        return result;
    }
}
