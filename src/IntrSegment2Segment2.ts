// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment2Segment2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test-intersection and find-intersection queries for two segments in 2D.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Each upstream
// class has two query functions: 'operator()', which uses the centered form
// of the segments (a Normalize call, so floating-point only), and 'Exact',
// which uses the endpoint form and supports rational arithmetic. Per
// PORTING.md the canonical two-argument query keeps 'test'/'find'; the second
// becomes 'testExact'/'findExact'. Both are ported, since 'Exact' is also the
// numerically cleaner query for float/double input (no normalization) and its
// parameters have a different meaning (endpoint form, t in [0,1]).

import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import type { Segment } from './Segment.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import { Line } from './Line.js';
import { IntrLine2Line2FI } from './IntrLine2Line2.js';
import { IntrIntervalsFI } from './IntrIntervals.js';

// The port of std::numeric_limits<int32_t>::max(), the IntrLine2Line2
// 'numIntersections' value meaning "the lines are the same".
const INT32_MAX = 2147483647;

// The result of IntrSegment2Segment2TI queries. The 'numIntersections' is 0
// (no intersection), 1 (segments intersect in a single point) or 2 (segments
// are collinear and intersect in a segment).
export interface IntrSegment2Segment2TIResult {
    intersect: boolean;
    numIntersections: number;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSegment2Segment2TIResult(): IntrSegment2Segment2TIResult {
    return { intersect: false, numIntersections: 0 };
}

// The result of IntrSegment2Segment2FI queries.
//
// The 'numIntersections' is 0 (no intersection), 1 (segments intersect in a
// single point) or 2 (segments are collinear and intersect in a segment).
//
// If numIntersections is 1, the intersection is
//   point[0]
//   = segment0.origin + segment0Parameter[0] * segment0.direction
//   = segment1.origin + segment1Parameter[0] * segment1.direction
// If numIntersections is 2, the endpoints of the segment of intersection are
//   point[i]
//   = segment0.origin + segment0Parameter[i] * segment0.direction
//   = segment1.origin + segment1Parameter[i] * segment1.direction
// with segment0Parameter[0] <= segment0Parameter[1] and
// segment1Parameter[0] <= segment1Parameter[1].
export interface IntrSegment2Segment2FIResult {
    intersect: boolean;
    numIntersections: number;
    segment0Parameter: [number, number];
    segment1Parameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment2Segment2FIResult(): IntrSegment2Segment2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        segment0Parameter: [0, 0],
        segment1Parameter: [0, 0],
        point: [Vector.zero(2), Vector.zero(2)]
    };
}

// Test-intersection queries for two segments in 2D.
export class IntrSegment2Segment2TI implements
    TIQuery<Segment, Segment, IntrSegment2Segment2TIResult> {

    // This version of the query uses Segment.getCenteredForm, which has a
    // normalize call. This generates rounding errors, so the query should be
    // used only with floating-point arithmetic.
    test(segment0: Segment, segment1: Segment): IntrSegment2Segment2TIResult {
        const result = defaultIntrSegment2Segment2TIResult();
        const { center: seg0Origin, direction: seg0Direction, extent: seg0Extent } =
            segment0.getCenteredForm();
        const { center: seg1Origin, direction: seg1Direction, extent: seg1Extent } =
            segment1.getCenteredForm();

        const llQuery = new IntrLine2Line2FI();
        const line0 = Line.fromOriginDirection(seg0Origin, seg0Direction);
        const line1 = Line.fromOriginDirection(seg1Origin, seg1Direction);
        const llResult = llQuery.find(line0, line1);
        if (llResult.numIntersections === 1) {
            // Test whether the line-line intersection is on the segments.
            if (Math.abs(llResult.line0Parameter[0]) <= seg0Extent &&
                Math.abs(llResult.line1Parameter[0]) <= seg1Extent) {
                result.intersect = true;
                result.numIntersections = 1;
            }
            else {
                result.intersect = false;
                result.numIntersections = 0;
            }
        }
        else if (llResult.numIntersections === INT32_MAX) {
            // Compute the location of segment1 endpoints relative to
            // segment0.
            const diff = sub(seg1Origin, seg0Origin);
            const t = dot(seg0Direction, diff);

            // Get the parameter intervals of the segments relative to
            // segment0.
            const interval0: [number, number] = [-seg0Extent, seg0Extent];
            const interval1: [number, number] = [t - seg1Extent, t + seg1Extent];

            // Compute the intersection of the intervals.
            const iiQuery = new IntrIntervalsFI();
            const iiResult = iiQuery.find(interval0, interval1);
            result.intersect = iiResult.intersect;
            result.numIntersections = iiResult.numIntersections;
        }
        else {
            result.intersect = false;
            result.numIntersections = 0;
        }

        return result;
    }

    // This version of the query supports rational arithmetic; it uses the
    // endpoint form of the segments and does not normalize.
    testExact(segment0: Segment, segment1: Segment): IntrSegment2Segment2TIResult {
        const result = defaultIntrSegment2Segment2TIResult();

        const llQuery = new IntrLine2Line2FI();
        const seg0Direction = sub(segment0.p[1], segment0.p[0]);
        const seg1Direction = sub(segment1.p[1], segment1.p[0]);
        const line0 = Line.fromOriginDirection(segment0.p[0], seg0Direction);
        const line1 = Line.fromOriginDirection(segment1.p[0], seg1Direction);
        const llResult = llQuery.find(line0, line1);
        if (llResult.numIntersections === 1) {
            // The lines are not parallel, so they intersect in a single
            // point. Test whether the line-line intersection is on the
            // segments. NOTE: upstream compares line0Parameter[1] and
            // line1Parameter[1] to 1 rather than the [0] values; this is
            // harmless because IntrLine2Line2FI sets both components of each
            // parameter pair to the same value when the lines intersect in a
            // single point. The comparison is preserved as upstream writes it.
            if (0 <= llResult.line0Parameter[0] && llResult.line0Parameter[1] <= 1 &&
                0 <= llResult.line1Parameter[0] && llResult.line1Parameter[1] <= 1) {
                result.intersect = true;
                result.numIntersections = 1;
            }
            else {
                result.intersect = false;
                result.numIntersections = 0;
            }
        }
        else if (llResult.numIntersections === INT32_MAX) {
            // The lines are the same. Compute the location of segment1
            // endpoints relative to segment0.
            const dotD0D0 = dot(seg0Direction, seg0Direction);
            let diff = sub(segment1.p[0], segment0.p[0]);
            const t0 = dot(seg0Direction, diff) / dotD0D0;
            diff = sub(segment1.p[1], segment0.p[0]);
            const t1 = dot(seg0Direction, diff) / dotD0D0;

            // Get the parameter intervals of the segments relative to
            // segment0.
            const interval0: [number, number] = [0, 1];
            const interval1: [number, number] =
                (t1 >= t0 ? [t0, t1] : [t1, t0]);

            // Compute the intersection of the intervals.
            const iiQuery = new IntrIntervalsFI();
            const iiResult = iiQuery.find(interval0, interval1);
            result.intersect = iiResult.intersect;
            result.numIntersections = iiResult.numIntersections;
        }
        else {
            // The lines are parallel but not the same, so the segments cannot
            // intersect.
            result.intersect = false;
            result.numIntersections = 0;
        }

        return result;
    }
}

// Find-intersection queries for two segments in 2D.
export class IntrSegment2Segment2FI implements
    FIQuery<Segment, Segment, IntrSegment2Segment2FIResult> {

    // This version of the query uses Segment.getCenteredForm, which has a
    // normalize call. This generates rounding errors, so the query should be
    // used only with floating-point arithmetic. NOTE: the parameters are
    // relative to the centered form of the segment. Each segment has a center
    // C, a unit-length direction D and an extent e > 0. A segment point is
    // C+t*D where |t| <= e.
    find(segment0: Segment, segment1: Segment): IntrSegment2Segment2FIResult {
        const result = defaultIntrSegment2Segment2FIResult();
        const { center: seg0Origin, direction: seg0Direction, extent: seg0Extent } =
            segment0.getCenteredForm();
        const { center: seg1Origin, direction: seg1Direction, extent: seg1Extent } =
            segment1.getCenteredForm();

        const llQuery = new IntrLine2Line2FI();
        const line0 = Line.fromOriginDirection(seg0Origin, seg0Direction);
        const line1 = Line.fromOriginDirection(seg1Origin, seg1Direction);
        const llResult = llQuery.find(line0, line1);
        if (llResult.numIntersections === 1) {
            // Test whether the line-line intersection is on the segments.
            if (Math.abs(llResult.line0Parameter[0]) <= seg0Extent &&
                Math.abs(llResult.line1Parameter[0]) <= seg1Extent) {
                result.intersect = true;
                result.numIntersections = 1;
                result.segment0Parameter[0] = llResult.line0Parameter[0];
                result.segment0Parameter[1] = result.segment0Parameter[0];
                result.segment1Parameter[0] = llResult.line1Parameter[0];
                result.segment1Parameter[1] = result.segment1Parameter[0];
                result.point[0] = llResult.point;
                result.point[1] = result.point[0];
            }
            else {
                result.intersect = false;
                result.numIntersections = 0;
            }
        }
        else if (llResult.numIntersections === INT32_MAX) {
            // Compute the location of segment1 endpoints relative to
            // segment0.
            const diff = sub(seg1Origin, seg0Origin);
            const t = dot(seg0Direction, diff);

            // Get the parameter intervals of the segments relative to
            // segment0.
            const interval0: [number, number] = [-seg0Extent, seg0Extent];
            const interval1: [number, number] = [t - seg1Extent, t + seg1Extent];

            // Compute the intersection of the intervals.
            const iiQuery = new IntrIntervalsFI();
            const iiResult = iiQuery.find(interval0, interval1);
            if (iiResult.intersect) {
                result.intersect = true;
                result.numIntersections = iiResult.numIntersections;
                for (let i = 0; i < iiResult.numIntersections; ++i) {
                    result.segment0Parameter[i] = iiResult.overlap[i];
                    result.segment1Parameter[i] = iiResult.overlap[i] - t;
                    result.point[i] = add(seg0Origin,
                        mul(result.segment0Parameter[i], seg0Direction));
                }
                if (iiResult.numIntersections === 1) {
                    result.segment0Parameter[1] = result.segment0Parameter[0];
                    result.segment1Parameter[1] = result.segment1Parameter[0];
                    result.point[1] = result.point[0];
                }
            }
            else {
                result.intersect = false;
                result.numIntersections = 0;
            }
        }
        else {
            result.intersect = false;
            result.numIntersections = 0;
        }

        return result;
    }

    // This version of the query supports rational arithmetic. NOTE: the
    // parameters are relative to the endpoint form of the segment. Each
    // segment has endpoints P0 and P1. A segment point is P0+t*(P1-P0) where
    // 0 <= t <= 1.
    findExact(segment0: Segment, segment1: Segment): IntrSegment2Segment2FIResult {
        const result = defaultIntrSegment2Segment2FIResult();

        const llQuery = new IntrLine2Line2FI();
        const seg0Direction = sub(segment0.p[1], segment0.p[0]);
        const seg1Direction = sub(segment1.p[1], segment1.p[0]);
        const line0 = Line.fromOriginDirection(segment0.p[0], seg0Direction);
        const line1 = Line.fromOriginDirection(segment1.p[0], seg1Direction);
        const llResult = llQuery.find(line0, line1);
        if (llResult.numIntersections === 1) {
            // The lines are not parallel, so they intersect in a single
            // point. Test whether the line-line intersection is on the
            // segments. (See the note in testExact about the [1] indices.)
            if (0 <= llResult.line0Parameter[0] && llResult.line0Parameter[1] <= 1 &&
                0 <= llResult.line1Parameter[0] && llResult.line1Parameter[1] <= 1) {
                result.intersect = true;
                result.numIntersections = 1;
                result.segment0Parameter[0] = llResult.line0Parameter[0];
                result.segment0Parameter[1] = result.segment0Parameter[0];
                result.segment1Parameter[0] = llResult.line1Parameter[0];
                result.segment1Parameter[1] = result.segment1Parameter[0];
                result.point[0] = llResult.point;
                result.point[1] = result.point[0];
            }
            else {
                result.intersect = false;
                result.numIntersections = 0;
            }
        }
        else if (llResult.numIntersections === INT32_MAX) {
            // The lines are the same. Compute the location of segment1
            // endpoints relative to segment0.
            const dotD0D0 = dot(seg0Direction, seg0Direction);
            let diff = sub(segment1.p[0], segment0.p[0]);
            const t0 = dot(seg0Direction, diff) / dotD0D0;
            diff = sub(segment1.p[1], segment0.p[0]);
            const t1 = dot(seg0Direction, diff) / dotD0D0;

            // Get the parameter intervals of the segments relative to
            // segment0.
            const interval0: [number, number] = [0, 1];
            const interval1: [number, number] =
                (t1 >= t0 ? [t0, t1] : [t1, t0]);

            // Compute the intersection of the intervals.
            const iiQuery = new IntrIntervalsFI();
            const iiResult = iiQuery.find(interval0, interval1);
            if (iiResult.intersect) {
                result.intersect = true;
                result.numIntersections = iiResult.numIntersections;

                // Compute the results for segment0.
                for (let i = 0; i < iiResult.numIntersections; ++i) {
                    result.segment0Parameter[i] = iiResult.overlap[i];
                    result.point[i] = add(segment0.p[0],
                        mul(result.segment0Parameter[i], seg0Direction));
                }

                // Compute the results for segment1. The interval1 was
                // computed relative to segment0, so we have to reverse the
                // process to obtain the parameters.
                const dotD1D1 = dot(seg1Direction, seg1Direction);
                for (let i = 0; i < iiResult.numIntersections; ++i) {
                    diff = sub(result.point[i], segment1.p[0]);
                    result.segment1Parameter[i] =
                        dot(seg1Direction, diff) / dotD1D1;
                }

                if (iiResult.numIntersections === 1) {
                    result.segment0Parameter[1] = result.segment0Parameter[0];
                    result.segment1Parameter[1] = result.segment1Parameter[0];
                    result.point[1] = result.point[0];
                }
                else {
                    if (t1 < t0) {
                        const tmp = result.segment1Parameter[0];
                        result.segment1Parameter[0] = result.segment1Parameter[1];
                        result.segment1Parameter[1] = tmp;
                    }
                }
            }
            else {
                result.intersect = false;
                result.numIntersections = 0;
            }
        }
        else {
            // The lines are parallel but not the same, so the segments cannot
            // intersect.
            result.intersect = false;
            result.numIntersections = 0;
        }

        return result;
    }
}
