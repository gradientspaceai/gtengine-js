// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment2AlignedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the box to be a solid.
//
// The test-intersection queries use the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection queries use parametric clipping against the four
// edges of the box.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// these queries from the Line2-vs-AlignedBox2 queries to reuse the protected
// DoQuery helpers; the port calls the exported module functions
// 'intrLine2AlignedBox2TIDoQuery' and 'intrLine2AlignedBox2FIDoQuery'
// instead. The segment-specific DoQuery helpers are exported as the module
// functions 'intrSegment2AlignedBox2TIDoQuery' and
// 'intrSegment2AlignedBox2FIDoQuery'.

import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import type { AlignedBox } from './AlignedBox.js';
import { Segment } from './Segment.js';
import { Vector, add, mul, sub } from './Vector.js';
import {
    intrLine2AlignedBox2TIDoQuery,
    intrLine2AlignedBox2FIDoQuery,
    defaultIntrLine2AlignedBox2TIResult,
    defaultIntrLine2AlignedBox2FIResult
} from './IntrLine2AlignedBox2.js';
import type {
    IntrLine2AlignedBox2TIResult,
    IntrLine2AlignedBox2FIResult
} from './IntrLine2AlignedBox2.js';
import { IntrIntervalsFI } from './IntrIntervals.js';
import { inContainerAlignedBox } from './ContAlignedBox.js';

// The upstream derived TIQuery::Result adds no members.
export type IntrSegment2AlignedBox2TIResult = IntrLine2AlignedBox2TIResult;

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSegment2AlignedBox2TIResult(): IntrSegment2AlignedBox2TIResult {
    return defaultIntrLine2AlignedBox2TIResult();
}

// The result of IntrSegment2AlignedBox2FI queries.
export interface IntrSegment2AlignedBox2FIResult extends IntrLine2AlignedBox2FIResult {
    // The base 'parameter' values are t-values for the segment
    // parameterization (1-t)*p[0] + t*p[1], where t is in [0,1]. The
    // 'cdeParameter' values are s-values for the centered form C + s * D,
    // where s is in [-e,e] and e is the extent of the segment.
    cdeParameter: [number, number];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment2AlignedBox2FIResult(): IntrSegment2AlignedBox2FIResult {
    return { ...defaultIntrLine2AlignedBox2FIResult(), cdeParameter: [0, 0] };
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
export function intrSegment2AlignedBox2TIDoQuery(segOrigin: Vector,
    segDirection: Vector, segExtent: number, boxExtent: Vector,
    result: IntrSegment2AlignedBox2TIResult): void {
    for (let i = 0; i < 2; ++i) {
        const lhs = Math.abs(segOrigin.values[i]);
        const rhs = boxExtent.values[i] +
            segExtent * Math.abs(segDirection.values[i]);
        if (lhs > rhs) {
            result.intersect = false;
            return;
        }
    }

    intrLine2AlignedBox2TIDoQuery(segOrigin, segDirection, boxExtent, result);
}

// The port of the protected 'FIQuery::DoQuery'.
export function intrSegment2AlignedBox2FIDoQuery(segOrigin: Vector,
    segDirection: Vector, segExtent: number, boxExtent: Vector,
    result: IntrSegment2AlignedBox2FIResult): void {
    intrLine2AlignedBox2FIDoQuery(segOrigin, segDirection, boxExtent, result);

    if (result.intersect) {
        // The line containing the segment intersects the box; the t-interval
        // is [t0,t1]. The segment intersects the box as long as [t0,t1]
        // overlaps the segment t-interval [-segExtent,+segExtent].
        const segInterval: [number, number] = [-segExtent, segExtent];
        const iiQuery = new IntrIntervalsFI();
        const iiResult = iiQuery.find(result.parameter, segInterval);
        result.intersect = iiResult.intersect;
        result.numIntersections = iiResult.numIntersections;
        result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];

        // If a segment intersects a box at an endpoint, and if that endpoint
        // is the only point of intersection, ensure the caller computes 2
        // points of intersection for a degenerate line segment representing a
        // single point.
        if (result.numIntersections === 1) {
            result.numIntersections = 2;
        }
    }
}

// Test-intersection query for a segment and a solid aligned box in 2D.
export class IntrSegment2AlignedBox2TI implements
    TIQuery<Segment, AlignedBox, IntrSegment2AlignedBox2TIResult> {

    test(segment: Segment, box: AlignedBox): IntrSegment2AlignedBox2TIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector2::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        const { segOrigin, segDirection, segExtent } =
            transformSegment(segment, boxCenter);

        const result = defaultIntrSegment2AlignedBox2TIResult();
        intrSegment2AlignedBox2TIDoQuery(segOrigin, segDirection, segExtent,
            boxExtent, result);
        return result;
    }
}

// Find-intersection query for a segment and a solid aligned box in 2D.
export class IntrSegment2AlignedBox2FI implements
    FIQuery<Segment, AlignedBox, IntrSegment2AlignedBox2FIResult> {

    find(segment: Segment, box: AlignedBox): IntrSegment2AlignedBox2FIResult {
        // The default result indicates no intersection.
        const result = defaultIntrSegment2AlignedBox2FIResult();

        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector2::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        const { segOrigin, segDirection, segExtent } =
            transformSegment(segment, boxCenter);

        if (segExtent > 0) {
            intrSegment2AlignedBox2FIDoQuery(segOrigin, segDirection,
                segExtent, boxExtent, result);
            for (let i = 0; i < result.numIntersections; ++i) {
                // Compute the point in the aligned-box coordinate system and
                // then translate it back to the original coordinates using
                // the box center.
                result.point[i] = add(boxCenter,
                    add(segOrigin, mul(result.parameter[i], segDirection)));
                result.cdeParameter[i] = result.parameter[i];

                // Convert the parameter from the centered form to the
                // endpoint form.
                result.parameter[i] =
                    (result.parameter[i] / segExtent + 1) * 0.5;
            }
        }
        else {
            // The segment is degenerate, representing a single point. Report
            // an intersection when this point is contained by the box.
            if (inContainerAlignedBox(segment.p[0], box)) {
                result.intersect = true;
                result.numIntersections = 1;
                result.parameter[0] = 0;
                result.parameter[1] = 0;
                // Upstream has the self-assignment
                // 'result.cdeParameter = result.cdeParameter;' here, which is
                // a no-op; cdeParameter is already [0,0] from the default
                // constructor, so the port simply omits it.
                result.point[0] = segment.p[0].clone();
                result.point[1] = segment.p[1].clone();
            }
        }

        return result;
    }
}
