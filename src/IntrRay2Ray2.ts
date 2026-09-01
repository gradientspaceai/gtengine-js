// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay2Ray2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test-intersection and find-intersection queries for two rays in 2D. The
// queries are built on the line-line queries, clamping the line parameters to
// the ray intervals.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The two upstream
// specializations become IntrRay2Ray2TI and IntrRay2Ray2FI.

import { Line } from './Line';
import type { Ray } from './Ray';
import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';
import { Vector, dot, sub } from './Vector';
import { IntrLine2Line2FI } from './IntrLine2Line2';

// The port of std::numeric_limits<int32_t>::max().
const INT32_MAX = 2147483647;

// The port of std::numeric_limits<T>::max() for T = double.
const MAX_T = Number.MAX_VALUE;

// The result of IntrRay2Ray2TI.test. The 'numIntersections' value is 0 (no
// intersection), 1 (rays intersect in a single point), 2 (rays are collinear
// and intersect in a segment; ray directions are opposite of each other) or
// 2147483647 (intersection is a ray; ray directions are the same).
export interface IntrRay2Ray2TIResult {
    intersect: boolean;
    numIntersections: number;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrRay2Ray2TIResult {
    return { intersect: false, numIntersections: 0 };
}

// The result of IntrRay2Ray2FI.find.
//
// The 'numIntersections' value is 0 (no intersection), 1 (rays intersect in a
// single point), 2 (rays are collinear and intersect in a segment; ray
// directions are opposite of each other) or 2147483647 (intersection is a
// ray; ray directions are the same).
//
// If numIntersections is 1, the intersection is
//   point[0] = ray0.origin + ray0Parameter[0] * ray0.direction
//            = ray1.origin + ray1Parameter[0] * ray1.direction
//
// If numIntersections is 2, the segment of intersection is formed by the ray
// origins,
//   ray0Parameter[0] = ray1Parameter[0] = 0
//   point[0] = ray0.origin
//            = ray1.origin + ray1Parameter[1] * ray1.direction
//   point[1] = ray1.origin
//            = ray0.origin + ray0Parameter[1] * ray0.direction
// where ray0Parameter[1] >= 0 and ray1Parameter[1] >= 0.
//
// If numIntersections is 2147483647, let
//   ray1.origin = ray0.origin + t * ray0.direction
// then
//   ray0Parameter = [max(t,0), +maxReal]
//   ray1Parameter = [-min(t,0), +maxReal]
//   point[0] = ray0.origin + ray0Parameter[0] * ray0.direction
export interface IntrRay2Ray2FIResult {
    intersect: boolean;
    numIntersections: number;
    ray0Parameter: [number, number];
    ray1Parameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrRay2Ray2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        ray0Parameter: [0, 0],
        ray1Parameter: [0, 0],
        point: [Vector.zero(2), Vector.zero(2)]
    };
}

// Test-intersection query for two rays in 2D.
export class IntrRay2Ray2TI implements
    TIQuery<Ray, Ray, IntrRay2Ray2TIResult> {

    test(ray0: Ray, ray1: Ray): IntrRay2Ray2TIResult {
        const result = defaultTIResult();

        const llQuery = new IntrLine2Line2FI();
        const line0 = Line.fromOriginDirection(ray0.origin, ray0.direction);
        const line1 = Line.fromOriginDirection(ray1.origin, ray1.direction);
        const llResult = llQuery.find(line0, line1);
        if (llResult.numIntersections === 1) {
            // Test whether the line-line intersection is on the rays.
            if (llResult.line0Parameter[0] >= 0 &&
                llResult.line1Parameter[0] >= 0) {
                result.intersect = true;
                result.numIntersections = 1;
            }
            else {
                result.intersect = false;
                result.numIntersections = 0;
            }
        }
        else if (llResult.numIntersections === INT32_MAX) {
            if (dot(ray0.direction, ray1.direction) > 0) {
                // The rays are collinear and in the same direction, so they
                // must overlap.
                result.intersect = true;
                result.numIntersections = INT32_MAX;
            }
            else {
                // The rays are collinear but in opposite directions. Test
                // whether they overlap. Ray0 has interval [0,+infinity) and
                // ray1 has interval (-infinity,t] relative to ray0.direction.
                const diff = sub(ray1.origin, ray0.origin);
                const t = dot(ray0.direction, diff);
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
        }
        else {
            result.intersect = false;
            result.numIntersections = 0;
        }

        return result;
    }
}

// Find-intersection query for two rays in 2D.
export class IntrRay2Ray2FI implements
    FIQuery<Ray, Ray, IntrRay2Ray2FIResult> {

    find(ray0: Ray, ray1: Ray): IntrRay2Ray2FIResult {
        const result = defaultFIResult();

        const llQuery = new IntrLine2Line2FI();
        const line0 = Line.fromOriginDirection(ray0.origin, ray0.direction);
        const line1 = Line.fromOriginDirection(ray1.origin, ray1.direction);
        const llResult = llQuery.find(line0, line1);
        if (llResult.numIntersections === 1) {
            // Test whether the line-line intersection is on the rays.
            if (llResult.line0Parameter[0] >= 0 &&
                llResult.line1Parameter[0] >= 0) {
                result.intersect = true;
                result.numIntersections = 1;
                result.ray0Parameter[0] = llResult.line0Parameter[0];
                result.ray1Parameter[0] = llResult.line1Parameter[0];
                result.point[0] = llResult.point;
            }
            else {
                result.intersect = false;
                result.numIntersections = 0;
            }
        }
        else if (llResult.numIntersections === INT32_MAX) {
            // Compute t for which ray1.origin =
            // ray0.origin + t * ray0.direction.
            const diff = sub(ray1.origin, ray0.origin);
            const t = dot(ray0.direction, diff);
            if (dot(ray0.direction, ray1.direction) > 0) {
                // The rays are collinear and in the same direction, so they
                // must overlap.
                result.intersect = true;
                result.numIntersections = INT32_MAX;
                if (t >= 0) {
                    result.ray0Parameter = [t, MAX_T];
                    result.ray1Parameter = [0, MAX_T];
                    result.point[0] = ray1.origin.clone();
                }
                else {
                    result.ray0Parameter = [0, MAX_T];
                    result.ray1Parameter = [-t, MAX_T];
                    result.point[0] = ray0.origin.clone();
                }
            }
            else {
                // The rays are collinear but in opposite directions. Test
                // whether they overlap. Ray0 has interval [0,+infinity) and
                // ray1 has interval (-infinity,t] relative to ray0.direction.
                if (t >= 0) {
                    result.intersect = true;
                    result.numIntersections = 2;
                    result.ray0Parameter = [0, t];
                    result.ray1Parameter = [0, t];
                    result.point[0] = ray0.origin.clone();
                    result.point[1] = ray1.origin.clone();
                }
                else {
                    result.intersect = false;
                    result.numIntersections = 0;
                }
            }
        }
        else {
            result.intersect = false;
            result.numIntersections = 0;
        }

        return result;
    }
}
