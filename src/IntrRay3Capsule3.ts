// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay3Capsule3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the capsule to be a solid.
//
// The test-intersection queries are based on distance computations.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// the FIQuery from the Line3-vs-Capsule3 query to reuse the protected DoQuery
// helper, which the port exports as the module function
// 'intrLine3Capsule3FIDoQuery'. The ray-specific DoQuery helper is exported
// as the module function 'intrRay3Capsule3FIDoQuery'.

import type { Capsule3 } from './Capsule.js';
import { DistRaySegment } from './DistRaySegment.js';
import type { FIQuery } from './FIQuery.js';
import { IntrIntervalsFI } from './IntrIntervals.js';
import {
    intrLine3Capsule3FIDoQuery,
    defaultIntrLine3Capsule3FIResult
} from './IntrLine3Capsule3.js';
import type { IntrLine3Capsule3FIResult } from './IntrLine3Capsule3.js';
import type { Ray3 } from './Ray.js';
import { Vector, add, mul } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrRay3Capsule3TI.test.
export interface IntrRay3Capsule3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrRay3Capsule3TIResult {
    return { intersect: false };
}

// The upstream derived FIQuery::Result adds no members.
export type IntrRay3Capsule3FIResult = IntrLine3Capsule3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrRay3Capsule3FIResult(): IntrRay3Capsule3FIResult {
    return defaultIntrLine3Capsule3FIResult();
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrRay3Capsule3FIDoQuery(rayOrigin: Vector,
    rayDirection: Vector, capsule: Capsule3,
    result: IntrRay3Capsule3FIResult): void {
    intrLine3Capsule3FIDoQuery(rayOrigin, rayDirection, capsule, result);

    if (result.intersect) {
        // The line containing the ray intersects the capsule; the t-interval
        // is [t0,t1]. The ray intersects the capsule as long as [t0,t1]
        // overlaps the ray t-interval [0,+infinity).
        const iiQuery = new IntrIntervalsFI();
        const iiResult = iiQuery.findFiniteSemiInfinite(result.parameter, 0,
            true);
        if (iiResult.intersect) {
            result.numIntersections = iiResult.numIntersections;
            result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
        }
        else {
            // The ray does not intersect the capsule.
            const empty = defaultIntrRay3Capsule3FIResult();
            result.intersect = empty.intersect;
            result.numIntersections = empty.numIntersections;
            result.parameter = empty.parameter;
            result.point = empty.point;
        }
    }
}

// Test-intersection query for a ray and a solid capsule in 3D.
export class IntrRay3Capsule3TI implements
    TIQuery<Ray3, Capsule3, IntrRay3Capsule3TIResult> {

    test(ray: Ray3, capsule: Capsule3): IntrRay3Capsule3TIResult {
        const result = defaultTIResult();
        const rsQuery = new DistRaySegment();
        const rsResult = rsQuery.compute(ray, capsule.segment);
        result.intersect = (rsResult.distance <= capsule.radius);
        return result;
    }
}

// Find-intersection query for a ray and a solid capsule in 3D.
export class IntrRay3Capsule3FI implements
    FIQuery<Ray3, Capsule3, IntrRay3Capsule3FIResult> {

    find(ray: Ray3, capsule: Capsule3): IntrRay3Capsule3FIResult {
        const result = defaultIntrRay3Capsule3FIResult();
        intrRay3Capsule3FIDoQuery(ray.origin, ray.direction, capsule, result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(ray.origin,
                    mul(result.parameter[i], ray.direction));
            }
        }
        return result;
    }
}
