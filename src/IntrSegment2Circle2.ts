// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment2Circle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the circle to be a solid (disk).
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// the FIQuery from the Line2-vs-Circle2 query to reuse the protected DoQuery
// helper; as in IntrSegment3Cylinder3, the port reaches that helper through a
// module-private accessor subclass. The segment-specific DoQuery helper is
// exported as the module function 'intrSegment2Circle2FIDoQuery'. The
// reported parameters are relative to the centered form of the segment,
// C + t * D with |t| <= e, as upstream reports them.

import type { FIQuery } from './FIQuery';
import type { Hypersphere } from './Hypersphere';
import { IntrIntervalsFI } from './IntrIntervals';
import {
    IntrLine2Circle2FI,
    defaultIntrLine2Circle2FIResult
} from './IntrLine2Circle2';
import type { IntrLine2Circle2FIResult } from './IntrLine2Circle2';
import type { Segment2 } from './Segment';
import { Vector, add, mul } from './Vector';
import type { TIQuery } from './TIQuery';

// The result of IntrSegment2Circle2TI.test.
export interface IntrSegment2Circle2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrSegment2Circle2TIResult {
    return { intersect: false };
}

// The upstream derived FIQuery::Result adds no members.
export type IntrSegment2Circle2FIResult = IntrLine2Circle2FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment2Circle2FIResult():
    IntrSegment2Circle2FIResult {
    return defaultIntrLine2Circle2FIResult();
}

// Expose the protected line-circle helper to this module.
class LineCircleFIAccess extends IntrLine2Circle2FI {
    run(lineOrigin: Vector, lineDirection: Vector, circle: Hypersphere,
        result: IntrLine2Circle2FIResult): void {
        this.doQuery(lineOrigin, lineDirection, circle, result);
    }
}

// The port of the protected 'FIQuery::DoQuery'.
export function intrSegment2Circle2FIDoQuery(segOrigin: Vector,
    segDirection: Vector, segExtent: number, circle: Hypersphere,
    result: IntrSegment2Circle2FIResult): void {
    new LineCircleFIAccess().run(segOrigin, segDirection, circle, result);

    if (result.intersect) {
        // The line containing the segment intersects the disk; the t-interval
        // is [t0,t1]. The segment intersects the disk as long as [t0,t1]
        // overlaps the segment t-interval [-segExtent,+segExtent].
        const segInterval: [number, number] = [-segExtent, segExtent];
        const iiQuery = new IntrIntervalsFI();
        const iiResult = iiQuery.find(result.parameter, segInterval);
        result.intersect = iiResult.intersect;
        result.numIntersections = iiResult.numIntersections;
        result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
    }
}

// Test-intersection query for a segment and a solid circle (disk) in 2D.
export class IntrSegment2Circle2TI implements
    TIQuery<Segment2, Hypersphere, IntrSegment2Circle2TIResult> {

    test(segment: Segment2, circle: Hypersphere):
        IntrSegment2Circle2TIResult {
        const result = defaultTIResult();
        const scQuery = new IntrSegment2Circle2FI();
        result.intersect = scQuery.find(segment, circle).intersect;
        return result;
    }
}

// Find-intersection query for a segment and a solid circle (disk) in 2D.
export class IntrSegment2Circle2FI implements
    FIQuery<Segment2, Hypersphere, IntrSegment2Circle2FIResult> {

    find(segment: Segment2, circle: Hypersphere):
        IntrSegment2Circle2FIResult {
        const { center: segOrigin, direction: segDirection,
            extent: segExtent } = segment.getCenteredForm();

        const result = defaultIntrSegment2Circle2FIResult();
        intrSegment2Circle2FIDoQuery(segOrigin, segDirection, segExtent,
            circle, result);
        for (let i = 0; i < result.numIntersections; ++i) {
            result.point[i] = add(segOrigin,
                mul(result.parameter[i], segDirection));
        }
        return result;
    }
}
