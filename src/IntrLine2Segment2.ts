// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine2Segment2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test-intersection and find-intersection queries for a line and a segment in
// 2D.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent.

import { Line } from './Line.js';
import type { Segment } from './Segment.js';
import { Vector, sub } from './Vector.js';
import { IntrLine2Line2FI } from './IntrLine2Line2.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';

// The port of std::numeric_limits<int32_t>::max(), the 'numIntersections'
// value meaning "the line and segment are collinear".
const INT32_MAX = 2147483647;

// The port of std::numeric_limits<T>::max() for T = double.
const MAX_T = Number.MAX_VALUE;

// The result of IntrLine2Segment2TI.test.
//
// If the line and segment do not intersect,
//   intersect = false
//   numIntersections = 0
//
// If the line and segment intersect in a single point,
//   intersect = true
//   numIntersections = 1
//
// If the line and segment are collinear,
//   intersect = true
//   numIntersections = 2147483647
export interface IntrLine2Segment2TIResult {
    intersect: boolean;
    numIntersections: number;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrLine2Segment2TIResult {
    return { intersect: false, numIntersections: 0 };
}

// The result of IntrLine2Segment2FI.find.
//
// If the line and segment do not intersect,
//   intersect = false
//   numIntersections = 0
//   lineParameter = [0, 0]  // invalid
//   segmentParameter = [0, 0]  // invalid
//   point = (0, 0)  // invalid
//
// If the line and segment intersect in a single point, the parameter for the
// line is s0 and the parameter for the segment is s1 in [0,1],
//   intersect = true
//   numIntersections = 1
//   lineParameter = [s0, s0]
//   segmentParameter = [s1, s1]
//   point = line.origin + s0 * line.direction
//         = segment.p[0] + s1 * (segment.p[1] - segment.p[0])
//
// If the line and segment are collinear, let maxT = Number.MAX_VALUE,
//   intersect = true
//   numIntersections = 2147483647
//   lineParameter = [-maxT, +maxT]
//   segmentParameter = [0, 1]
//   point = (0, 0)  // invalid
export interface IntrLine2Segment2FIResult {
    intersect: boolean;
    numIntersections: number;
    lineParameter: [number, number];
    segmentParameter: [number, number];
    point: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrLine2Segment2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        lineParameter: [0, 0],
        segmentParameter: [0, 0],
        point: Vector.zero(2)
    };
}

export class IntrLine2Segment2TI implements
    TIQuery<Line, Segment, IntrLine2Segment2TIResult> {

    test(line: Line, segment: Segment): IntrLine2Segment2TIResult {
        const result = defaultTIResult();

        const llQuery = new IntrLine2Line2FI();
        const segLine = Line.fromOriginDirection(segment.p[0],
            sub(segment.p[1], segment.p[0]));
        const llResult = llQuery.find(line, segLine);
        if (llResult.numIntersections === 1) {
            // Test whether the line-line intersection is on the segment.
            if (llResult.line1Parameter[0] >= 0 &&
                llResult.line1Parameter[1] <= 1) {
                result.intersect = true;
                result.numIntersections = 1;
            } else {
                result.intersect = false;
                result.numIntersections = 0;
            }
        } else {
            result.intersect = llResult.intersect;
            result.numIntersections = llResult.numIntersections;
        }

        return result;
    }
}

export class IntrLine2Segment2FI implements
    FIQuery<Line, Segment, IntrLine2Segment2FIResult> {

    find(line: Line, segment: Segment): IntrLine2Segment2FIResult {
        const result = defaultFIResult();

        const llQuery = new IntrLine2Line2FI();
        const segLine = Line.fromOriginDirection(segment.p[0],
            sub(segment.p[1], segment.p[0]));
        const llResult = llQuery.find(line, segLine);
        if (llResult.numIntersections === 1) {
            // Test whether the line-line intersection is on the segment.
            if (llResult.line1Parameter[0] >= 0 &&
                llResult.line1Parameter[1] <= 1) {
                result.intersect = true;
                result.numIntersections = 1;
                result.lineParameter[0] = llResult.line0Parameter[0];
                result.lineParameter[1] = result.lineParameter[0];
                result.segmentParameter[0] = llResult.line1Parameter[0];
                result.segmentParameter[1] = result.segmentParameter[0];
                result.point = llResult.point;
            } else {
                result.intersect = false;
                result.numIntersections = 0;
            }
        } else if (llResult.numIntersections === INT32_MAX) {
            result.intersect = true;
            result.numIntersections = INT32_MAX;
            result.lineParameter[0] = -MAX_T;
            result.lineParameter[1] = +MAX_T;
            result.segmentParameter[0] = 0;
            result.segmentParameter[1] = 1;
        } else {
            result.intersect = false;
            result.numIntersections = 0;
        }

        return result;
    }
}
