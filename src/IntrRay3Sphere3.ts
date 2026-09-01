// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay3Sphere3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the sphere to be a solid.
//
// The sphere is (X-C)^T*(X-C)-r^2 = 0 and the ray is X = P+t*D for t >= 0.
// Substitute the ray equation into the sphere equation to obtain a quadratic
// equation Q(t) = t^2 + 2*a1*t + a0 = 0, where a1 = D^T*(P-C) and
// a0 = (P-C)^T*(P-C)-r^2. The algorithm involves an analysis of the
// real-valued roots of Q(t) for t >= 0.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// FIQuery<Ray3,Sphere3> from FIQuery<Line3,Sphere3> only to reuse the
// protected DoQuery; the derived Result adds no members, so the result type
// is an alias of the line-sphere result type. The line-sphere DoQuery is the
// exported module function 'intrLine3Sphere3DoQuery', and the ray-sphere one
// is exported here as 'intrRay3Sphere3DoQuery' for the same reason.

import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';
import type { Hypersphere } from './Hypersphere';
import type { Ray } from './Ray';
import { Vector, add, dot, mul, sub } from './Vector';
import {
    intrLine3Sphere3DoQuery,
    defaultIntrLine3Sphere3FIResult
} from './IntrLine3Sphere3';
import type { IntrLine3Sphere3FIResult } from './IntrLine3Sphere3';
import { IntrIntervalsFI } from './IntrIntervals';

// The result of IntrRay3Sphere3TI queries.
export interface IntrRay3Sphere3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrRay3Sphere3TIResult(): IntrRay3Sphere3TIResult {
    return { intersect: false };
}

// The upstream derived FIQuery::Result adds no members.
export type IntrRay3Sphere3FIResult = IntrLine3Sphere3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrRay3Sphere3FIResult(): IntrRay3Sphere3FIResult {
    return defaultIntrLine3Sphere3FIResult();
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrRay3Sphere3DoQuery(rayOrigin: Vector,
    rayDirection: Vector, sphere: Hypersphere,
    result: IntrRay3Sphere3FIResult): void {
    intrLine3Sphere3DoQuery(rayOrigin, rayDirection, sphere, result);

    if (result.intersect) {
        // The line containing the ray intersects the sphere; the t-interval
        // is [t0,t1]. The ray intersects the sphere as long as [t0,t1]
        // overlaps the ray t-interval [0,+infinity).
        const iiQuery = new IntrIntervalsFI();
        const iiResult = iiQuery.findFiniteSemiInfinite(result.parameter, 0, true);
        if (iiResult.intersect) {
            result.numIntersections = iiResult.numIntersections;
            result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
        }
        else {
            // The line containing the ray does not intersect the sphere.
            const empty = defaultIntrRay3Sphere3FIResult();
            result.intersect = empty.intersect;
            result.numIntersections = empty.numIntersections;
            result.parameter = empty.parameter;
            result.point = empty.point;
        }
    }
}

// Test-intersection query for a ray and a solid sphere in 3D.
export class IntrRay3Sphere3TI implements
    TIQuery<Ray, Hypersphere, IntrRay3Sphere3TIResult> {

    test(ray: Ray, sphere: Hypersphere): IntrRay3Sphere3TIResult {
        const result = defaultIntrRay3Sphere3TIResult();

        const diff = sub(ray.origin, sphere.center);
        const a0 = dot(diff, diff) - sphere.radius * sphere.radius;
        if (a0 <= 0) {
            // P is inside the sphere.
            result.intersect = true;
            return result;
        }
        // else: P is outside the sphere

        const a1 = dot(ray.direction, diff);
        if (a1 >= 0) {
            result.intersect = false;
            return result;
        }

        // Intersection occurs when Q(t) has real roots.
        const discr = a1 * a1 - a0;
        result.intersect = (discr >= 0);
        return result;
    }
}

// Find-intersection query for a ray and a solid sphere in 3D.
export class IntrRay3Sphere3FI implements
    FIQuery<Ray, Hypersphere, IntrRay3Sphere3FIResult> {

    find(ray: Ray, sphere: Hypersphere): IntrRay3Sphere3FIResult {
        const result = defaultIntrRay3Sphere3FIResult();
        intrRay3Sphere3DoQuery(ray.origin, ray.direction, sphere, result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(ray.origin,
                    mul(result.parameter[i], ray.direction));
            }
        }
        return result;
    }
}
