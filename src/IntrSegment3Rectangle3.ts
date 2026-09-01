// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment3Rectangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the intersection between a segment and a solid rectangle in 3D.
//
// The segment is P0 + t * (P1 - P0) for 0 <= t <= 1. The direction D = P1-P0
// is generally not unit length.
//
// The rectangle has center C, unit-length axis directions W[0] and W[1], and
// extents e[0] and e[1]. A rectangle point is X = C + sum_i s[i] * W[i] where
// |s[i]| <= e[i] for all i.
//
// The intersection point, if any, is stored in result.point. The
// corresponding segment parameter t is stored in result.parameter. The
// corresponding rectangle parameters s[] are stored in result.rectCoord[].
// When the segment is in the plane of the rectangle and intersects the
// rectangle, the queries state that there are no intersections.
//
// Upstream TODO: modify to support non-unit-length W[]; return the point or
// segment of intersection when the segment is in the plane of the rectangle
// and intersects the rectangle.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The 'rectCoord'
// field has three components even though the rectangle has two parameters;
// that quirk is inherited from IntrLine3Rectangle3 (issue #141), where the
// third component is always left at zero.

import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';
import type { Rectangle } from './Rectangle';
import type { Segment } from './Segment';
import { Line } from './Line';
import { Vector, sub } from './Vector';
import { IntrLine3Rectangle3FI } from './IntrLine3Rectangle3';

// The result of IntrSegment3Rectangle3TI queries.
export interface IntrSegment3Rectangle3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSegment3Rectangle3TIResult(): IntrSegment3Rectangle3TIResult {
    return { intersect: false };
}

// The result of IntrSegment3Rectangle3FI queries.
export interface IntrSegment3Rectangle3FIResult {
    intersect: boolean;

    // The segment parameter t at the intersection point.
    parameter: number;

    // The rectangle parameters s[] at the intersection point.
    rectCoord: [number, number, number];

    // The point of intersection.
    point: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment3Rectangle3FIResult(): IntrSegment3Rectangle3FIResult {
    return {
        intersect: false,
        parameter: 0,
        rectCoord: [0, 0, 0],
        point: Vector.zero(3)
    };
}

// Test-intersection query for a segment and a solid rectangle in 3D.
export class IntrSegment3Rectangle3TI implements
    TIQuery<Segment, Rectangle, IntrSegment3Rectangle3TIResult> {

    test(segment: Segment, rectangle: Rectangle): IntrSegment3Rectangle3TIResult {
        const result = defaultIntrSegment3Rectangle3TIResult();

        const lrQuery = new IntrLine3Rectangle3FI();
        const line = Line.fromOriginDirection(segment.p[0],
            sub(segment.p[1], segment.p[0]));
        const lrResult = lrQuery.find(line, rectangle);
        if (lrResult.intersect) {
            if (0 <= lrResult.parameter && lrResult.parameter <= 1) {
                // The line-rectangle intersection is on the segment.
                result.intersect = true;
                return result;
            }
        }

        result.intersect = false;
        return result;
    }
}

// Find-intersection query for a segment and a solid rectangle in 3D.
export class IntrSegment3Rectangle3FI implements
    FIQuery<Segment, Rectangle, IntrSegment3Rectangle3FIResult> {

    find(segment: Segment, rectangle: Rectangle): IntrSegment3Rectangle3FIResult {
        const result = defaultIntrSegment3Rectangle3FIResult();

        const lrQuery = new IntrLine3Rectangle3FI();
        const line = Line.fromOriginDirection(segment.p[0],
            sub(segment.p[1], segment.p[0]));
        const lrResult = lrQuery.find(line, rectangle);
        if (lrResult.intersect) {
            if (0 <= lrResult.parameter && lrResult.parameter <= 1) {
                // The line-rectangle intersection is on the segment.
                result.intersect = true;
                result.parameter = lrResult.parameter;
                result.rectCoord = [
                    lrResult.rectCoord[0],
                    lrResult.rectCoord[1],
                    lrResult.rectCoord[2]
                ];
                result.point = lrResult.point.clone();
                return result;
            }
        }

        result.intersect = false;
        return result;
    }
}
