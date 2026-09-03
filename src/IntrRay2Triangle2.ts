// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay2Triangle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the triangle to be a solid.
//
// Port notes (see IntrIntervals.ts for the Intr* precedent): upstream derives
// the ray-triangle FIQuery from the line-triangle FIQuery only to reuse the
// protected DoQuery member, which the port exports as the module function
// 'intrLine2Triangle2DoQuery'. The ray-specific DoQuery is exported here as
// 'intrRay2Triangle2DoQuery' for the same reason. The upstream FIQuery Result
// adds no members to the line-triangle result, so the port exports a type
// alias.

import type { Ray } from './Ray.js';
import type { Triangle } from './Triangle.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import { Vector, add, mul } from './Vector.js';
import {
    intrLine2Triangle2DoQuery,
    defaultIntrLine2Triangle2FIResult
} from './IntrLine2Triangle2.js';
import type { IntrLine2Triangle2FIResult } from './IntrLine2Triangle2.js';
import { IntrIntervalsFI } from './IntrIntervals.js';

// The result of IntrRay2Triangle2TI queries.
export interface IntrRay2Triangle2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrRay2Triangle2TIResult {
    return { intersect: false };
}

// The result of IntrRay2Triangle2FI.find. Upstream adds no members to the
// line-triangle result.
export type IntrRay2Triangle2FIResult = IntrLine2Triangle2FIResult;

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrRay2Triangle2DoQuery(origin: Vector, direction: Vector,
    triangle: Triangle, result: IntrRay2Triangle2FIResult): void {
    intrLine2Triangle2DoQuery(origin, direction, triangle, result);

    if (result.intersect) {
        // The line containing the ray intersects the triangle; the
        // t-interval is [t0,t1]. The ray intersects the triangle as long
        // as [t0,t1] overlaps the ray t-interval [0,+infinity).
        const iiQuery = new IntrIntervalsFI();
        const iiResult = iiQuery.findFiniteSemiInfinite(result.parameter,
            0, true);
        if (iiResult.intersect) {
            result.numIntersections = iiResult.numIntersections;
            result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
        }
        else {
            // The ray does not intersect the triangle.
            result.intersect = false;
            result.numIntersections = 0;
            result.parameter = [0, 0];
            result.point = [Vector.zero(2), Vector.zero(2)];
        }
    }
}

// Test-intersection query for a ray and a solid triangle in 2D.
export class IntrRay2Triangle2TI implements
    TIQuery<Ray, Triangle, IntrRay2Triangle2TIResult> {

    // The ray is P + t * D, where P is a point on the line and D is a
    // direction vector that does not have to be unit length. This is useful
    // when using a 2-point representation P0 + t * (P1 - P0). The t-parameter
    // is constrained by t >= 0.
    test(ray: Ray, triangle: Triangle): IntrRay2Triangle2TIResult {
        const result = defaultTIResult();
        const rtQuery = new IntrRay2Triangle2FI();
        result.intersect = rtQuery.find(ray, triangle).intersect;
        return result;
    }
}

// Find-intersection query for a ray and a solid triangle in 2D.
export class IntrRay2Triangle2FI implements
    FIQuery<Ray, Triangle, IntrRay2Triangle2FIResult> {

    // The ray is P + t * D, where P is a point on the line and D is a
    // direction vector that does not have to be unit length. This is useful
    // when using a 2-point representation P0 + t * (P1 - P0).
    find(ray: Ray, triangle: Triangle): IntrRay2Triangle2FIResult {
        const result = defaultIntrLine2Triangle2FIResult();
        intrRay2Triangle2DoQuery(ray.origin, ray.direction, triangle,
            result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(ray.origin,
                    mul(result.parameter[i], ray.direction));
            }
        }
        return result;
    }
}
