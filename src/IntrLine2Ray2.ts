// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine2Ray2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test-intersection and find-intersection queries for a line and a ray in 2D.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent.

import { Line } from './Line.js';
import type { Ray } from './Ray.js';
import { Vector } from './Vector.js';
import { IntrLine2Line2FI } from './IntrLine2Line2.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';

// The port of std::numeric_limits<int32_t>::max(), the 'numIntersections'
// value meaning "the line and ray are collinear".
const INT32_MAX = 2147483647;

// The port of std::numeric_limits<T>::max() for T = double.
const MAX_T = Number.MAX_VALUE;

// The result of IntrLine2Ray2TI.test.
//
// If the line and ray do not intersect,
//   intersect = false
//   numIntersections = 0
//
// If the line and ray intersect in a single point,
//   intersect = true
//   numIntersections = 1
//
// If the line and ray are collinear,
//   intersect = true
//   numIntersections = 2147483647
export interface IntrLine2Ray2TIResult {
    intersect: boolean;
    numIntersections: number;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrLine2Ray2TIResult {
    return { intersect: false, numIntersections: 0 };
}

// The result of IntrLine2Ray2FI.find.
//
// If the line and ray do not intersect,
//   intersect = false
//   numIntersections = 0
//   lineParameter = [0, 0]  // invalid
//   rayParameter = [0, 0]  // invalid
//   point = (0, 0)  // invalid
//
// If the line and ray intersect in a single point, the parameter for the line
// is s0 and the parameter for the ray is s1 >= 0,
//   intersect = true
//   numIntersections = 1
//   lineParameter = [s0, s0]
//   rayParameter = [s1, s1]
//   point = line.origin + s0 * line.direction
//         = ray.origin + s1 * ray.direction
//
// If the line and ray are collinear, let maxT = Number.MAX_VALUE,
//   intersect = true
//   numIntersections = 2147483647
//   lineParameter = [-maxT, +maxT]
//   rayParameter = [0, +maxT]
//   point = (0, 0)  // invalid
export interface IntrLine2Ray2FIResult {
    intersect: boolean;
    numIntersections: number;
    lineParameter: [number, number];
    rayParameter: [number, number];
    point: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrLine2Ray2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        lineParameter: [0, 0],
        rayParameter: [0, 0],
        point: Vector.zero(2)
    };
}

export class IntrLine2Ray2TI implements
    TIQuery<Line, Ray, IntrLine2Ray2TIResult> {

    test(line: Line, ray: Ray): IntrLine2Ray2TIResult {
        const result = defaultTIResult();

        const llQuery = new IntrLine2Line2FI();
        const llResult = llQuery.find(line,
            Line.fromOriginDirection(ray.origin, ray.direction));
        if (llResult.numIntersections === 1) {
            // Test whether the line-line intersection is on the ray.
            if (llResult.line1Parameter[0] >= 0) {
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

export class IntrLine2Ray2FI implements
    FIQuery<Line, Ray, IntrLine2Ray2FIResult> {

    find(line: Line, ray: Ray): IntrLine2Ray2FIResult {
        const result = defaultFIResult();

        const llQuery = new IntrLine2Line2FI();
        const llResult = llQuery.find(line,
            Line.fromOriginDirection(ray.origin, ray.direction));
        if (llResult.numIntersections === 1) {
            // Test whether the line-line intersection is on the ray.
            if (llResult.line1Parameter[0] >= 0) {
                result.intersect = true;
                result.numIntersections = 1;
                result.lineParameter[0] = llResult.line0Parameter[0];
                result.lineParameter[1] = result.lineParameter[0];
                result.rayParameter[0] = llResult.line1Parameter[0];
                result.rayParameter[1] = result.rayParameter[0];
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
            result.rayParameter[0] = 0;
            result.rayParameter[1] = +MAX_T;
        } else {
            result.intersect = false;
            result.numIntersections = 0;
        }

        return result;
    }
}
