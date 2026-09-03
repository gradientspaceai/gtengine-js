// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment3OrientedBox3.h
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
// these queries from the Segment3-vs-AlignedBox3 queries to reuse the
// protected DoQuery helpers; the port calls the module functions that
// IntrSegment3AlignedBox3.ts exports for them. The derived Result structs add
// no members, so they are type aliases here. The reported parameters are
// relative to the centered form of the segment, C + t * D with |t| <= e, as
// upstream reports them.

import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import {
    intrSegment3AlignedBox3TIDoQuery,
    intrSegment3AlignedBox3FIDoQuery,
    defaultIntrSegment3AlignedBox3TIResult,
    defaultIntrSegment3AlignedBox3FIResult
} from './IntrSegment3AlignedBox3.js';
import type {
    IntrSegment3AlignedBox3TIResult,
    IntrSegment3AlignedBox3FIResult
} from './IntrSegment3AlignedBox3.js';
import type { OrientedBox3 } from './OrientedBox.js';
import type { Segment3 } from './Segment.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

// The upstream derived TIQuery::Result adds no members.
export type IntrSegment3OrientedBox3TIResult =
    IntrSegment3AlignedBox3TIResult;

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSegment3OrientedBox3TIResult():
    IntrSegment3OrientedBox3TIResult {
    return defaultIntrSegment3AlignedBox3TIResult();
}

// The upstream derived FIQuery::Result adds no members.
export type IntrSegment3OrientedBox3FIResult =
    IntrSegment3AlignedBox3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment3OrientedBox3FIResult():
    IntrSegment3OrientedBox3FIResult {
    return defaultIntrSegment3AlignedBox3FIResult();
}

// Transform the segment to the oriented-box coordinate system, returning the
// centered form of the transformed segment.
function transformSegment(segment: Segment3, box: OrientedBox3):
    { segOrigin: Vector, segDirection: Vector, segExtent: number } {
    const { center: tmpOrigin, direction: tmpDirection, extent: segExtent } =
        segment.getCenteredForm();
    const diff = sub(tmpOrigin, box.center);
    const segOrigin = Vector.fromArray([
        dot(diff, box.axis[0]),
        dot(diff, box.axis[1]),
        dot(diff, box.axis[2])
    ]);
    const segDirection = Vector.fromArray([
        dot(tmpDirection, box.axis[0]),
        dot(tmpDirection, box.axis[1]),
        dot(tmpDirection, box.axis[2])
    ]);
    return { segOrigin, segDirection, segExtent };
}

// Test-intersection query for a segment and a solid oriented box in 3D.
export class IntrSegment3OrientedBox3TI implements
    TIQuery<Segment3, OrientedBox3, IntrSegment3OrientedBox3TIResult> {

    test(segment: Segment3, box: OrientedBox3):
        IntrSegment3OrientedBox3TIResult {
        const { segOrigin, segDirection, segExtent } =
            transformSegment(segment, box);

        const result = defaultIntrSegment3OrientedBox3TIResult();
        intrSegment3AlignedBox3TIDoQuery(segOrigin, segDirection, segExtent,
            box.extent, result);
        return result;
    }
}

// Find-intersection query for a segment and a solid oriented box in 3D.
export class IntrSegment3OrientedBox3FI implements
    FIQuery<Segment3, OrientedBox3, IntrSegment3OrientedBox3FIResult> {

    find(segment: Segment3, box: OrientedBox3):
        IntrSegment3OrientedBox3FIResult {
        const { segOrigin, segDirection, segExtent } =
            transformSegment(segment, box);

        const result = defaultIntrSegment3OrientedBox3FIResult();
        intrSegment3AlignedBox3FIDoQuery(segOrigin, segDirection, segExtent,
            box.extent, result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                // Compute the intersection point in the oriented-box
                // coordinate system.
                const y = add(segOrigin,
                    mul(result.parameter[i], segDirection));

                // Transform the intersection point to the original
                // coordinate system.
                let point = box.center.clone();
                for (let j = 0; j < 3; ++j) {
                    point = add(point, mul(y.values[j], box.axis[j]));
                }
                result.point[i] = point;
            }
        }
        return result;
    }
}
