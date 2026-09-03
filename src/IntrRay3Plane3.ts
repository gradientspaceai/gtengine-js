// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay3Plane3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the plane to be a 2-dimensional object embedded in 3D
// and the ray to be a 1-dimensional object.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has both
// a TIQuery and an FIQuery specialization, which become IntrRay3Plane3TI and
// IntrRay3Plane3FI. The upstream FIQuery derives from
// FIQuery<Real, Line3, Plane3> and adds no result members, so the port
// aliases IntrLine3Plane3FIResult (the precedent set by
// IntrSegment2Triangle2.ts) and reuses the exported
// 'intrLine3Plane3FIDoQuery'. The protected 'FIQuery::DoQuery' becomes the
// exported module function 'intrRay3Plane3FIDoQuery'.

import { DistPointHyperplane } from './DistPointHyperplane';
import type { FIQuery } from './FIQuery';
import {
    defaultIntrLine3Plane3FIResult, intrLine3Plane3FIDoQuery
} from './IntrLine3Plane3';
import type { IntrLine3Plane3FIResult } from './IntrLine3Plane3';
import type { Plane3 } from './Hyperplane';
import type { Ray3 } from './Ray';
import type { TIQuery } from './TIQuery';
import { Vector, add, dot, mul } from './Vector';

// The result of IntrRay3Plane3TI.test.
export interface IntrRay3Plane3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrRay3Plane3TIResult(): IntrRay3Plane3TIResult {
    return { intersect: false };
}

// The result of IntrRay3Plane3FI.find. Upstream derives its Result from the
// line-plane Result and adds no members.
export type IntrRay3Plane3FIResult = IntrLine3Plane3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrRay3Plane3FIResult(): IntrRay3Plane3FIResult {
    return defaultIntrLine3Plane3FIResult();
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
export function intrRay3Plane3FIDoQuery(rayOrigin: Vector,
    rayDirection: Vector, plane: Plane3,
    result: IntrRay3Plane3FIResult): void {
    intrLine3Plane3FIDoQuery(rayOrigin, rayDirection, plane, result);

    if (result.intersect) {
        // The line intersects the plane in a point that might not be on the
        // ray.
        if (result.parameter < 0) {
            result.intersect = false;
            result.numIntersections = 0;
        }
    }
}

// Test-intersection query for a ray and a plane in 3D.
export class IntrRay3Plane3TI implements
    TIQuery<Ray3, Plane3, IntrRay3Plane3TIResult> {

    test(ray: Ray3, plane: Plane3): IntrRay3Plane3TIResult {
        const result = defaultIntrRay3Plane3TIResult();

        // Compute the (signed) distance from the ray origin to the plane.
        const vpQuery = new DistPointHyperplane();
        const vpResult = vpQuery.compute(ray.origin, plane);

        const DdN = dot(ray.direction, plane.normal);
        if (DdN > 0) {
            // The ray is not parallel to the plane and is directed toward the
            // +normal side of the plane.
            result.intersect = (vpResult.signedDistance <= 0);
        }
        else if (DdN < 0) {
            // The ray is not parallel to the plane and is directed toward the
            // -normal side of the plane.
            result.intersect = (vpResult.signedDistance >= 0);
        }
        else {
            // The ray and plane are parallel.
            result.intersect = (vpResult.distance === 0);
        }

        return result;
    }
}

// Find-intersection query for a ray and a plane in 3D.
export class IntrRay3Plane3FI implements
    FIQuery<Ray3, Plane3, IntrRay3Plane3FIResult> {

    find(ray: Ray3, plane: Plane3): IntrRay3Plane3FIResult {
        const result = defaultIntrRay3Plane3FIResult();
        intrRay3Plane3FIDoQuery(ray.origin, ray.direction, plane, result);
        if (result.intersect) {
            result.point = add(ray.origin, mul(result.parameter, ray.direction));
        }
        return result;
    }
}
