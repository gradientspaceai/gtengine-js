// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment3Cylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query considers the cylinder to be a solid.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream provides
// only an FIQuery specialization for this pair of primitives, which becomes
// IntrSegment3Cylinder3FI, and derives it from the Line3-vs-Cylinder3 query
// to reuse the protected DoQuery helper; as in IntrSegment3AlignedBox3, the
// port reaches that helper through a module-private accessor subclass because
// a TypeScript subclass cannot change an inherited method signature. The
// segment-specific DoQuery helper is exported as the module function
// 'intrSegment3Cylinder3FIDoQuery'. The reported parameters are relative to
// the centered form of the segment, C + t * D with |t| <= e, as upstream
// reports them.
//
// As in IntrLine3Cylinder3, upstream has no infinite-cylinder branch: it
// reads cylinder.height directly, so an infinite cylinder (height = -1) would
// produce meaningless results. Following the precedent set by
// IntrHalfspace3Cylinder3 and IntrTriangle3Cylinder3, the port asserts that
// the cylinder is finite instead of silently computing nonsense.
//
// NOTE (upstream behavior, preserved): a zero-length segment has no
// well-defined centered-form direction (Normalize of the zero vector leaves
// it zero), so the line-cylinder DoQuery divides by a vanishing quadratic
// leading coefficient. For a zero-length segment whose point is inside the
// cylinder, the result reports an intersection with NaN parameters. Upstream
// behaves identically; callers must pass a segment with distinct endpoints.

import type { Cylinder3 } from './Cylinder3';
import type { FIQuery } from './FIQuery';
import {
    IntrLine3Cylinder3FI,
    defaultIntrLine3Cylinder3FIResult
} from './IntrLine3Cylinder3';
import type { IntrLine3Cylinder3FIResult } from './IntrLine3Cylinder3';
import { IntrIntervalsFI } from './IntrIntervals';
import { logAssert } from './Logger';
import type { Segment3 } from './Segment';
import { Vector, add, mul } from './Vector';

// The upstream derived FIQuery::Result adds no members.
export type IntrSegment3Cylinder3FIResult = IntrLine3Cylinder3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment3Cylinder3FIResult():
    IntrSegment3Cylinder3FIResult {
    return defaultIntrLine3Cylinder3FIResult();
}

// Expose the protected line-cylinder helper to this module.
class FIHelper extends IntrLine3Cylinder3FI {
    runDoQuery(lineOrigin: Vector, lineDirection: Vector,
        cylinder: Cylinder3, result: IntrLine3Cylinder3FIResult): void {
        this.doQuery(lineOrigin, lineDirection, cylinder, result);
    }
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrSegment3Cylinder3FIDoQuery(segOrigin: Vector,
    segDirection: Vector, segExtent: number, cylinder: Cylinder3,
    result: IntrSegment3Cylinder3FIResult): void {
    new FIHelper().runDoQuery(segOrigin, segDirection, cylinder, result);

    if (result.intersect) {
        // The line containing the segment intersects the cylinder; the
        // t-interval is [t0,t1]. The segment intersects the cylinder as long
        // as [t0,t1] overlaps the segment t-interval [-segExtent,+segExtent].
        const iiQuery = new IntrIntervalsFI();
        const segInterval: [number, number] = [-segExtent, segExtent];
        const iiResult = iiQuery.find(result.parameter, segInterval);
        if (iiResult.intersect) {
            result.numIntersections = iiResult.numIntersections;
            result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
        }
        else {
            // The segment does not intersect the cylinder.
            const empty = defaultIntrSegment3Cylinder3FIResult();
            result.intersect = empty.intersect;
            result.numIntersections = empty.numIntersections;
            result.parameter = empty.parameter;
            result.point = empty.point;
        }
    }
}

// Find-intersection query for a segment and a solid cylinder in 3D.
export class IntrSegment3Cylinder3FI implements
    FIQuery<Segment3, Cylinder3, IntrSegment3Cylinder3FIResult> {

    find(segment: Segment3, cylinder: Cylinder3):
        IntrSegment3Cylinder3FIResult {
        logAssert(cylinder.isFinite(),
            'Infinite cylinders are not yet supported.');

        const { center: segOrigin, direction: segDirection,
            extent: segExtent } = segment.getCenteredForm();

        const result = defaultIntrSegment3Cylinder3FIResult();
        intrSegment3Cylinder3FIDoQuery(segOrigin, segDirection, segExtent,
            cylinder, result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(segOrigin,
                    mul(result.parameter[i], segDirection));
            }
        }
        return result;
    }
}
