// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment2OrientedBox2.h
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
// these queries from the Segment2-vs-AlignedBox2 queries to reuse their
// protected DoQuery helpers, which the port already exposes as the module
// functions intrSegment2AlignedBox2TIDoQuery / intrSegment2AlignedBox2FIDoQuery.
//
// Upstream bug (fixed here): the FIQuery computes the intersection point as
//   result.point[i] = box.center + (segOrigin + parameter[i] * segDirection)
// where segOrigin and segDirection are expressed in the box coordinate
// system. Adding a box-frame vector directly to the world-space box center is
// valid only when the box axes are the standard coordinate axes (which is why
// the identical line in IntrSegment2AlignedBox2.h is correct). For a rotated
// box the reported points are wrong. The port computes the points from the
// untransformed (world-space) centered form of the segment, which is the
// mathematically equivalent correction.
//
// NOTE (upstream quirk, preserved): for a degenerate (zero-length) segment
// contained by the box, this query reports numIntersections = 2, whereas
// IntrSegment2AlignedBox2 reports 1 for the same configuration. Upstream also
// has the no-op self-assignment 'result.cdeParameter = result.cdeParameter;'
// in that branch; cdeParameter is already [0,0] from the default constructor,
// so the port simply omits it.

import { inContainerOrientedBox2 } from './ContOrientedBox2.js';
import type { FIQuery } from './FIQuery.js';
import {
    intrSegment2AlignedBox2TIDoQuery,
    intrSegment2AlignedBox2FIDoQuery,
    defaultIntrSegment2AlignedBox2TIResult,
    defaultIntrSegment2AlignedBox2FIResult
} from './IntrSegment2AlignedBox2.js';
import type {
    IntrSegment2AlignedBox2TIResult,
    IntrSegment2AlignedBox2FIResult
} from './IntrSegment2AlignedBox2.js';
import { logAssert } from './Logger.js';
import type { OrientedBox2 } from './OrientedBox.js';
import type { Segment2 } from './Segment.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The upstream derived TIQuery::Result adds no members.
export type IntrSegment2OrientedBox2TIResult = IntrSegment2AlignedBox2TIResult;

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSegment2OrientedBox2TIResult():
    IntrSegment2OrientedBox2TIResult {
    return defaultIntrSegment2AlignedBox2TIResult();
}

// The upstream derived FIQuery::Result adds the same 'cdeParameter' member
// that the segment-vs-aligned-box result already has.
export type IntrSegment2OrientedBox2FIResult = IntrSegment2AlignedBox2FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment2OrientedBox2FIResult():
    IntrSegment2OrientedBox2FIResult {
    return defaultIntrSegment2AlignedBox2FIResult();
}

// Transform the segment to a centered form in the oriented-box coordinate
// system. The world-space centered form is returned as well, because the
// intersection points are computed from it.
function transformSegment(segment: Segment2, box: OrientedBox2): {
    tmpOrigin: Vector, tmpDirection: Vector,
    segOrigin: Vector, segDirection: Vector, segExtent: number
} {
    const { center: tmpOrigin, direction: tmpDirection, extent: segExtent } =
        segment.getCenteredForm();
    const diff = sub(tmpOrigin, box.center);
    const segOrigin = Vector.zero(2);
    const segDirection = Vector.zero(2);
    for (let i = 0; i < 2; ++i) {
        segOrigin.values[i] = dot(diff, box.axis[i]);
        segDirection.values[i] = dot(tmpDirection, box.axis[i]);
    }
    return { tmpOrigin, tmpDirection, segOrigin, segDirection, segExtent };
}

// Test-intersection query for a segment and a solid oriented box in 2D.
export class IntrSegment2OrientedBox2TI implements
    TIQuery<Segment2, OrientedBox2, IntrSegment2OrientedBox2TIResult> {

    test(segment: Segment2, box: OrientedBox2):
        IntrSegment2OrientedBox2TIResult {
        logAssert(box.dimension === 2,
            'IntrSegment2OrientedBox2TI: mismatched sizes.');

        const { segOrigin, segDirection, segExtent } =
            transformSegment(segment, box);

        const result = defaultIntrSegment2OrientedBox2TIResult();
        intrSegment2AlignedBox2TIDoQuery(segOrigin, segDirection, segExtent,
            box.extent, result);
        return result;
    }
}

// Find-intersection query for a segment and a solid oriented box in 2D.
export class IntrSegment2OrientedBox2FI implements
    FIQuery<Segment2, OrientedBox2, IntrSegment2OrientedBox2FIResult> {

    find(segment: Segment2, box: OrientedBox2):
        IntrSegment2OrientedBox2FIResult {
        logAssert(box.dimension === 2,
            'IntrSegment2OrientedBox2FI: mismatched sizes.');

        // The default result indicates no intersection.
        const result = defaultIntrSegment2OrientedBox2FIResult();

        const { tmpOrigin, tmpDirection, segOrigin, segDirection, segExtent } =
            transformSegment(segment, box);

        if (segExtent > 0) {
            intrSegment2AlignedBox2FIDoQuery(segOrigin, segDirection,
                segExtent, box.extent, result);
            for (let i = 0; i < result.numIntersections; ++i) {
                // Evaluate the world-space centered form of the segment. See
                // the upstream bug note at the top of this file.
                result.point[i] = add(tmpOrigin,
                    mul(result.parameter[i], tmpDirection));
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
            if (inContainerOrientedBox2(segment.p[0], box)) {
                result.intersect = true;
                result.numIntersections = 2;
                result.parameter[0] = 0;
                result.parameter[1] = 0;
                result.point[0] = segment.p[0].clone();
                result.point[1] = segment.p[1].clone();
            }
        }

        return result;
    }
}
