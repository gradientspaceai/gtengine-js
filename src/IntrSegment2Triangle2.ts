// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment2Triangle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the triangle to be a solid.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// FIQuery<Segment2,Triangle2> from FIQuery<Line2,Triangle2> only to reuse the
// protected DoQuery; the derived Result adds no members, so the result type
// is an alias of the line-triangle result type. The protected line-triangle
// helper is the exported module function 'intrLine2Triangle2DoQuery'.

import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import type { Segment } from './Segment.js';
import type { Triangle } from './Triangle.js';
import { Vector, add, mul, sub } from './Vector.js';
import {
    intrLine2Triangle2DoQuery,
    defaultIntrLine2Triangle2FIResult
} from './IntrLine2Triangle2.js';
import type { IntrLine2Triangle2FIResult } from './IntrLine2Triangle2.js';
import { IntrIntervalsFI } from './IntrIntervals.js';

// The result of IntrSegment2Triangle2TI queries.
export interface IntrSegment2Triangle2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSegment2Triangle2TIResult(): IntrSegment2Triangle2TIResult {
    return { intersect: false };
}

// The upstream derived FIQuery::Result adds no members.
export type IntrSegment2Triangle2FIResult = IntrLine2Triangle2FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment2Triangle2FIResult(): IntrSegment2Triangle2FIResult {
    return defaultIntrLine2Triangle2FIResult();
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrSegment2Triangle2DoQuery(origin: Vector,
    direction: Vector, triangle: Triangle,
    result: IntrSegment2Triangle2FIResult): void {
    intrLine2Triangle2DoQuery(origin, direction, triangle, result);

    if (result.intersect) {
        // The line containing the segment intersects the triangle; the
        // t-interval is [t0,t1]. The segment intersects the triangle as long
        // as [t0,t1] overlaps the segment t-interval [0,1].
        const iiQuery = new IntrIntervalsFI();
        const segInterval: [number, number] = [0, 1];
        const iiResult = iiQuery.find(result.parameter, segInterval);
        if (iiResult.intersect) {
            result.numIntersections = iiResult.numIntersections;
            result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
        }
        else {
            // The line containing the segment does not intersect the
            // triangle.
            const empty = defaultIntrSegment2Triangle2FIResult();
            result.intersect = empty.intersect;
            result.numIntersections = empty.numIntersections;
            result.parameter = empty.parameter;
            result.point = empty.point;
        }
    }
}

// Test-intersection query for a segment and a solid triangle in 2D.
export class IntrSegment2Triangle2TI implements
    TIQuery<Segment, Triangle, IntrSegment2Triangle2TIResult> {

    // The segment is P0 + t * (P1 - P0) for t in [0,1].
    test(segment: Segment, triangle: Triangle): IntrSegment2Triangle2TIResult {
        const result = defaultIntrSegment2Triangle2TIResult();
        const stQuery = new IntrSegment2Triangle2FI();
        result.intersect = stQuery.find(segment, triangle).intersect;
        return result;
    }
}

// Find-intersection query for a segment and a solid triangle in 2D.
export class IntrSegment2Triangle2FI implements
    FIQuery<Segment, Triangle, IntrSegment2Triangle2FIResult> {

    // The segment is P0 + t * (P1 - P0) for t in [0,1].
    find(segment: Segment, triangle: Triangle): IntrSegment2Triangle2FIResult {
        const result = defaultIntrSegment2Triangle2FIResult();
        const segOrigin = segment.p[0];
        const segDirection = sub(segment.p[1], segment.p[0]);
        intrSegment2Triangle2DoQuery(segOrigin, segDirection, triangle, result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(segOrigin,
                    mul(result.parameter[i], segDirection));
            }
        }
        return result;
    }
}
