// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment2Arc2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the arc to be a 1-dimensional object.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has both
// a TIQuery and an FIQuery specialization, which become IntrSegment2Arc2TI
// and IntrSegment2Arc2FI. Upstream calls the single-argument Arc2::Contains,
// which the port names 'containsOnCircle' (see IntrLine2Arc2.ts).
//
// Port fix for an upstream bug. Upstream runs FIQuery<Segment2,Circle2>,
// which treats the circle as a SOLID disk: it clips the line t-interval to
// [-extent,+extent], so when a segment endpoint is inside the disk the
// reported "intersection point" is that endpoint, which is not on the circle.
// Upstream then hands the interior point to Arc2::Contains, whose
// single-argument form assumes the point is on the circle, and the query can
// report an intersection at a point that is not on the arc. The port instead
// intersects the segment with the circular CURVE (the line-circle query in
// the segment's centered form, filtered by |t| <= extent, the same technique
// upstream itself uses in IntrSegment2SegmentMesh2.h), which yields the same
// answers whenever both segment endpoints are outside the disk and the
// correct answers otherwise. As upstream, the reported parameters are for the
// centered form segOrigin + t * segDirection with t in [-extent, extent].
//
// Upstream's FIQuery::operator() re-initializes the result members that the
// Result default constructor already set, and in doing so assigns
// 'result.parameter[0] = 0' twice rather than setting 'result.parameter[1]'.
// The typo is harmless (both parameters are already zero) and the redundant
// re-initialization is dropped here; the default result factory covers it.

import type { Arc2 } from './Arc2';
import type { FIQuery } from './FIQuery';
import { Hypersphere } from './Hypersphere';
import { IntrLine2Circle2FI } from './IntrLine2Circle2';
import { Line } from './Line';
import type { Segment2 } from './Segment';
import { Vector } from './Vector';
import type { TIQuery } from './TIQuery';

// The result of IntrSegment2Arc2TI.test.
export interface IntrSegment2Arc2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSegment2Arc2TIResult(): IntrSegment2Arc2TIResult {
    return { intersect: false };
}

// The result of IntrSegment2Arc2FI.find.
export interface IntrSegment2Arc2FIResult {
    intersect: boolean;
    numIntersections: number;
    parameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment2Arc2FIResult(): IntrSegment2Arc2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        parameter: [0, 0],
        point: [Vector.zero(2), Vector.zero(2)]
    };
}

// Test-intersection query for a segment and an arc in 2D.
export class IntrSegment2Arc2TI implements
    TIQuery<Segment2, Arc2, IntrSegment2Arc2TIResult> {

    test(segment: Segment2, arc: Arc2): IntrSegment2Arc2TIResult {
        const result = defaultIntrSegment2Arc2TIResult();
        const saQuery = new IntrSegment2Arc2FI();
        const saResult = saQuery.find(segment, arc);
        result.intersect = saResult.intersect;
        return result;
    }
}

// Find-intersection query for a segment and an arc in 2D.
export class IntrSegment2Arc2FI implements
    FIQuery<Segment2, Arc2, IntrSegment2Arc2FIResult> {

    find(segment: Segment2, arc: Arc2): IntrSegment2Arc2FIResult {
        const result = defaultIntrSegment2Arc2FIResult();

        const { center: segOrigin, direction: segDirection,
            extent: segExtent } = segment.getCenteredForm();

        const lcQuery = new IntrLine2Circle2FI();
        const circle = Hypersphere.fromCenterRadius(arc.center, arc.radius);
        const line = Line.fromOriginDirection(segOrigin, segDirection);
        const scResult = lcQuery.find(line, circle);
        if (scResult.intersect) {
            // Test whether the line-circle intersections are on the segment
            // and on the arc.
            for (let i = 0; i < scResult.numIntersections; ++i) {
                if (Math.abs(scResult.parameter[i]) <= segExtent
                    && arc.containsOnCircle(scResult.point[i])) {
                    result.intersect = true;
                    result.parameter[result.numIntersections] =
                        scResult.parameter[i];
                    result.point[result.numIntersections] =
                        scResult.point[i].clone();
                    ++result.numIntersections;
                }
            }
        }

        return result;
    }
}
