// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment3AlignedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The test-intersection queries use the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection queries use Liang-Barsky parametric clipping against
// the six faces of the box. The queries consider the box to be a solid. The
// algorithms are described in
// https://www.geometrictools.com/Documentation/IntersectionLineBox.pdf
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// these queries from the Line3-vs-AlignedBox3 queries to reuse the protected
// DoQuery helpers; as in IntrLine2OrientedBox2, the port reaches those
// helpers through module-private subclasses because a TypeScript subclass
// cannot change an inherited method signature. The segment-specific DoQuery
// helpers are exported as the module functions
// 'intrSegment3AlignedBox3TIDoQuery' and 'intrSegment3AlignedBox3FIDoQuery'.
// The reported parameters are relative to the centered form of the segment,
// C + t * D with |t| <= e, as upstream reports them.

import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';
import type { AlignedBox } from './AlignedBox';
import { Segment } from './Segment';
import { Vector, add, mul, sub } from './Vector';
import {
    IntrLine3AlignedBox3TI,
    IntrLine3AlignedBox3FI,
    defaultIntrLine3AlignedBox3TIResult,
    defaultIntrLine3AlignedBox3FIResult
} from './IntrLine3AlignedBox3';
import type {
    IntrLine3AlignedBox3TIResult,
    IntrLine3AlignedBox3FIResult
} from './IntrLine3AlignedBox3';
import { IntrIntervalsFI } from './IntrIntervals';

// The upstream derived TIQuery::Result adds no members.
export type IntrSegment3AlignedBox3TIResult = IntrLine3AlignedBox3TIResult;

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSegment3AlignedBox3TIResult(): IntrSegment3AlignedBox3TIResult {
    return defaultIntrLine3AlignedBox3TIResult();
}

// The upstream derived FIQuery::Result adds no members.
export type IntrSegment3AlignedBox3FIResult = IntrLine3AlignedBox3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment3AlignedBox3FIResult(): IntrSegment3AlignedBox3FIResult {
    return defaultIntrLine3AlignedBox3FIResult();
}

// Expose the protected line-box helpers to this module.
class TIHelper extends IntrLine3AlignedBox3TI {
    runDoQuery(lineOrigin: Vector, lineDirection: Vector, boxExtent: Vector,
        result: IntrLine3AlignedBox3TIResult): void {
        this.doQuery(lineOrigin, lineDirection, boxExtent, result);
    }
}

class FIHelper extends IntrLine3AlignedBox3FI {
    runDoQuery(lineOrigin: Vector, lineDirection: Vector, boxExtent: Vector,
        result: IntrLine3AlignedBox3FIResult): void {
        this.doQuery(lineOrigin, lineDirection, boxExtent, result);
    }
}

// Transform the segment to a centered form in the aligned-box coordinate
// system.
function transformSegment(segment: Segment, boxCenter: Vector):
    { segOrigin: Vector, segDirection: Vector, segExtent: number } {
    const transformedSegment = Segment.fromEndpoints(
        sub(segment.p[0], boxCenter), sub(segment.p[1], boxCenter));
    const { center, direction, extent } = transformedSegment.getCenteredForm();
    return { segOrigin: center, segDirection: direction, segExtent: extent };
}

// The port of the protected 'TIQuery::DoQuery'.
export function intrSegment3AlignedBox3TIDoQuery(segOrigin: Vector,
    segDirection: Vector, segExtent: number, boxExtent: Vector,
    result: IntrSegment3AlignedBox3TIResult): void {
    for (let i = 0; i < 3; ++i) {
        if (Math.abs(segOrigin.values[i]) > boxExtent.values[i] +
            segExtent * Math.abs(segDirection.values[i])) {
            result.intersect = false;
            return;
        }
    }

    new TIHelper().runDoQuery(segOrigin, segDirection, boxExtent, result);
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrSegment3AlignedBox3FIDoQuery(segOrigin: Vector,
    segDirection: Vector, segExtent: number, boxExtent: Vector,
    result: IntrSegment3AlignedBox3FIResult): void {
    new FIHelper().runDoQuery(segOrigin, segDirection, boxExtent, result);

    if (result.intersect) {
        // The line containing the segment intersects the box; the t-interval
        // is [t0,t1]. The segment intersects the box as long as [t0,t1]
        // overlaps the segment t-interval [-segExtent,+segExtent].
        const iiQuery = new IntrIntervalsFI();
        const segInterval: [number, number] = [-segExtent, segExtent];
        const iiResult = iiQuery.find(result.parameter, segInterval);
        if (iiResult.intersect) {
            result.numIntersections = iiResult.numIntersections;
            result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
        }
        else {
            // The line containing the segment does not intersect the box.
            const empty = defaultIntrSegment3AlignedBox3FIResult();
            result.intersect = empty.intersect;
            result.numIntersections = empty.numIntersections;
            result.parameter = empty.parameter;
            result.point = empty.point;
        }
    }
}

// Test-intersection query for a segment and a solid aligned box in 3D.
export class IntrSegment3AlignedBox3TI implements
    TIQuery<Segment, AlignedBox, IntrSegment3AlignedBox3TIResult> {

    test(segment: Segment, box: AlignedBox): IntrSegment3AlignedBox3TIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector3::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        const { segOrigin, segDirection, segExtent } =
            transformSegment(segment, boxCenter);

        const result = defaultIntrSegment3AlignedBox3TIResult();
        intrSegment3AlignedBox3TIDoQuery(segOrigin, segDirection, segExtent,
            boxExtent, result);
        return result;
    }
}

// Find-intersection query for a segment and a solid aligned box in 3D.
export class IntrSegment3AlignedBox3FI implements
    FIQuery<Segment, AlignedBox, IntrSegment3AlignedBox3FIResult> {

    find(segment: Segment, box: AlignedBox): IntrSegment3AlignedBox3FIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector3::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        const { segOrigin, segDirection, segExtent } =
            transformSegment(segment, boxCenter);

        const result = defaultIntrSegment3AlignedBox3FIResult();
        intrSegment3AlignedBox3FIDoQuery(segOrigin, segDirection, segExtent,
            boxExtent, result);
        if (result.intersect) {
            // The segment origin is in aligned-box coordinates. Transform it
            // back to the original space.
            const worldSegOrigin = add(segOrigin, boxCenter);
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(worldSegOrigin,
                    mul(result.parameter[i], segDirection));
            }
        }
        return result;
    }
}
