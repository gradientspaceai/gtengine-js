// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay2Segment2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test-intersection and find-intersection queries for a ray and a segment in
// 2D. The queries are built on the line-line queries, clamping the line
// parameters to the ray and segment intervals.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The two upstream
// specializations become IntrRay2Segment2TI and IntrRay2Segment2FI.

import { Line } from './Line';
import type { Ray } from './Ray';
import type { Segment } from './Segment';
import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';
import { Vector, add, dot, mul, sub } from './Vector';
import { IntrLine2Line2FI } from './IntrLine2Line2';
import { IntrIntervalsFI } from './IntrIntervals';

// The port of std::numeric_limits<int32_t>::max().
const INT32_MAX = 2147483647;

// The port of std::numeric_limits<T>::max() for T = double.
const MAX_T = Number.MAX_VALUE;

// The result of IntrRay2Segment2TI.test. The 'numIntersections' value is 0
// (no intersection), 1 (ray and segment intersect in a single point) or 2
// (ray and segment are collinear and intersect in a segment).
export interface IntrRay2Segment2TIResult {
    intersect: boolean;
    numIntersections: number;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrRay2Segment2TIResult {
    return { intersect: false, numIntersections: 0 };
}

// The result of IntrRay2Segment2FI.find.
//
// The 'numIntersections' value is 0 (no intersection), 1 (ray and segment
// intersect in a single point) or 2 (ray and segment are collinear and
// intersect in a segment).
//
// If numIntersections is 1, the intersection is
//   point[0] = ray.origin + rayParameter[0] * ray.direction
//            = segment.center + segmentParameter[0] * segment.direction
// If numIntersections is 2, the endpoints of the segment of intersection are
//   point[i] = ray.origin + rayParameter[i] * ray.direction
//            = segment.center + segmentParameter[i] * segment.direction
// with rayParameter[0] <= rayParameter[1] and
// segmentParameter[0] <= segmentParameter[1].
export interface IntrRay2Segment2FIResult {
    intersect: boolean;
    numIntersections: number;
    rayParameter: [number, number];
    segmentParameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrRay2Segment2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        rayParameter: [0, 0],
        segmentParameter: [0, 0],
        point: [Vector.zero(2), Vector.zero(2)]
    };
}

// Test-intersection query for a ray and a segment in 2D.
export class IntrRay2Segment2TI implements
    TIQuery<Ray, Segment, IntrRay2Segment2TIResult> {

    test(ray: Ray, segment: Segment): IntrRay2Segment2TIResult {
        const result = defaultTIResult();

        const { center: segOrigin, direction: segDirection, extent: segExtent }
            = segment.getCenteredForm();

        const llQuery = new IntrLine2Line2FI();
        const line0 = Line.fromOriginDirection(ray.origin, ray.direction);
        const line1 = Line.fromOriginDirection(segOrigin, segDirection);
        const llResult = llQuery.find(line0, line1);
        if (llResult.numIntersections === 1) {
            // Test whether the line-line intersection is on the ray and
            // segment.
            if (llResult.line0Parameter[0] >= 0 &&
                Math.abs(llResult.line1Parameter[0]) <= segExtent) {
                result.intersect = true;
                result.numIntersections = 1;
            }
            else {
                result.intersect = false;
                result.numIntersections = 0;
            }
        }
        else if (llResult.numIntersections === INT32_MAX) {
            // Compute the location of the right-most point of the segment
            // relative to the ray direction.
            const diff = sub(segOrigin, ray.origin);
            const t = dot(ray.direction, diff) + segExtent;
            if (t > 0) {
                result.intersect = true;
                result.numIntersections = 2;
            }
            else if (t < 0) {
                result.intersect = false;
                result.numIntersections = 0;
            }
            else {
                // t == 0
                result.intersect = true;
                result.numIntersections = 1;
            }
        }
        else {
            result.intersect = false;
            result.numIntersections = 0;
        }

        return result;
    }
}

// Find-intersection query for a ray and a segment in 2D.
export class IntrRay2Segment2FI implements
    FIQuery<Ray, Segment, IntrRay2Segment2FIResult> {

    find(ray: Ray, segment: Segment): IntrRay2Segment2FIResult {
        const result = defaultFIResult();

        const { center: segOrigin, direction: segDirection, extent: segExtent }
            = segment.getCenteredForm();

        const llQuery = new IntrLine2Line2FI();
        const line0 = Line.fromOriginDirection(ray.origin, ray.direction);
        const line1 = Line.fromOriginDirection(segOrigin, segDirection);
        const llResult = llQuery.find(line0, line1);
        if (llResult.numIntersections === 1) {
            // Test whether the line-line intersection is on the ray and
            // segment.
            if (llResult.line0Parameter[0] >= 0 &&
                Math.abs(llResult.line1Parameter[0]) <= segExtent) {
                result.intersect = true;
                result.numIntersections = 1;
                result.rayParameter[0] = llResult.line0Parameter[0];
                result.segmentParameter[0] = llResult.line1Parameter[0];
                result.point[0] = llResult.point;
            }
            else {
                result.intersect = false;
                result.numIntersections = 0;
            }
        }
        else if (llResult.numIntersections === INT32_MAX) {
            // Compute t for which segment.origin =
            // ray.origin + t * ray.direction.
            const diff = sub(segOrigin, ray.origin);
            const t = dot(ray.direction, diff);

            // Get the ray interval.
            const interval0: [number, number] = [0, MAX_T];

            // Compute the location of the segment endpoints relative to the
            // ray.
            const interval1: [number, number] = [t - segExtent, t + segExtent];

            // Compute the intersection of [0,+infinity) and [tmin,tmax].
            const iiQuery = new IntrIntervalsFI();
            const iiResult = iiQuery.find(interval0, interval1);
            if (iiResult.intersect) {
                result.intersect = true;
                result.numIntersections = iiResult.numIntersections;
                for (let i = 0; i < iiResult.numIntersections; ++i) {
                    result.rayParameter[i] = iiResult.overlap[i];
                    result.segmentParameter[i] = iiResult.overlap[i] - t;
                    result.point[i] = add(ray.origin,
                        mul(result.rayParameter[i], ray.direction));
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
}
