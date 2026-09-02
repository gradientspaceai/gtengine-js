// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay3Ellipsoid3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the ellipsoid to be a solid.
//
// The ellipsoid is (X-C)^T*M*(X-C)-1 = 0 and the ray is X = P+t*D for t >= 0.
// Substitute the ray equation into the ellipsoid equation to obtain a
// quadratic equation Q(t) = a2*t^2 + 2*a1*t + a0 = 0, where a2 = D^T*M*D,
// a1 = D^T*M*(P-C) and a0 = (P-C)^T*M*(P-C)-1. The algorithm involves an
// analysis of the real-valued roots of Q(t) for t >= 0.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// FIQuery<Ray3,Ellipsoid3> from FIQuery<Line3,Ellipsoid3> only to reuse the
// protected DoQuery; the derived Result adds no members, so the result type
// is an alias of the line-ellipsoid result type. The line-ellipsoid DoQuery
// is the exported module function 'intrLine3Ellipsoid3FIDoQuery', and the
// ray-ellipsoid one is exported here as 'intrRay3Ellipsoid3FIDoQuery' for the
// same reason (the precedent set by IntrRay3Sphere3.ts).

import type { Ellipsoid3 } from './Hyperellipsoid';
import type { FIQuery } from './FIQuery';
import { IntrIntervalsFI } from './IntrIntervals';
import {
    intrLine3Ellipsoid3FIDoQuery,
    defaultIntrLine3Ellipsoid3FIResult
} from './IntrLine3Ellipsoid3';
import type { IntrLine3Ellipsoid3FIResult } from './IntrLine3Ellipsoid3';
import { logAssert } from './Logger';
import { mulMatrix } from './Matrix';
import type { Ray3 } from './Ray';
import { Vector, add, dot, mul, sub } from './Vector';
import type { TIQuery } from './TIQuery';

// The result of IntrRay3Ellipsoid3TI.test.
export interface IntrRay3Ellipsoid3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrRay3Ellipsoid3TIResult():
    IntrRay3Ellipsoid3TIResult {
    return { intersect: false };
}

// The upstream derived FIQuery::Result adds no members.
export type IntrRay3Ellipsoid3FIResult = IntrLine3Ellipsoid3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrRay3Ellipsoid3FIResult():
    IntrRay3Ellipsoid3FIResult {
    return defaultIntrLine3Ellipsoid3FIResult();
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrRay3Ellipsoid3FIDoQuery(rayOrigin: Vector,
    rayDirection: Vector, ellipsoid: Ellipsoid3,
    result: IntrRay3Ellipsoid3FIResult): void {
    intrLine3Ellipsoid3FIDoQuery(rayOrigin, rayDirection, ellipsoid, result);

    if (result.intersect) {
        // The line containing the ray intersects the ellipsoid; the
        // t-interval is [t0,t1]. The ray intersects the ellipsoid as long as
        // [t0,t1] overlaps the ray t-interval [0,+infinity).
        const iiQuery = new IntrIntervalsFI();
        const iiResult = iiQuery.findFiniteSemiInfinite(result.parameter, 0, true);
        if (iiResult.intersect) {
            result.numIntersections = iiResult.numIntersections;
            result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
        }
        else {
            // The line containing the ray does not intersect the ellipsoid.
            const empty = defaultIntrRay3Ellipsoid3FIResult();
            result.intersect = empty.intersect;
            result.numIntersections = empty.numIntersections;
            result.parameter = empty.parameter;
            result.point = empty.point;
        }
    }
}

// Test-intersection query for a ray and a solid ellipsoid in 3D.
export class IntrRay3Ellipsoid3TI implements
    TIQuery<Ray3, Ellipsoid3, IntrRay3Ellipsoid3TIResult> {

    test(ray: Ray3, ellipsoid: Ellipsoid3): IntrRay3Ellipsoid3TIResult {
        logAssert(ellipsoid.dimension === 3 && ray.origin.size === 3,
            'IntrRay3Ellipsoid3TI: mismatched sizes.');
        const result = defaultIntrRay3Ellipsoid3TIResult();

        const M = ellipsoid.getM();
        const diff = sub(ray.origin, ellipsoid.center);
        const matDir = mulMatrix(M, ray.direction) as Vector;
        const matDiff = mulMatrix(M, diff) as Vector;
        const a0 = dot(diff, matDiff) - 1;
        if (a0 <= 0) {
            // P is inside the ellipsoid.
            result.intersect = true;
            return result;
        }
        // else: P is outside the ellipsoid

        const a1 = dot(ray.direction, matDiff);
        if (a1 >= 0) {
            // Q(t) >= a0 > 0 for t >= 0, so Q(t) cannot be zero for t in
            // [0,+infinity) and the ray does not intersect the ellipsoid.
            result.intersect = false;
            return result;
        }

        // The minimum of Q(t) occurs for some t in (0,+infinity). An
        // intersection occurs when Q(t) has real roots.
        const a2 = dot(ray.direction, matDir);
        const discr = a1 * a1 - a0 * a2;
        result.intersect = (discr >= 0);
        return result;
    }
}

// Find-intersection query for a ray and a solid ellipsoid in 3D.
export class IntrRay3Ellipsoid3FI implements
    FIQuery<Ray3, Ellipsoid3, IntrRay3Ellipsoid3FIResult> {

    find(ray: Ray3, ellipsoid: Ellipsoid3): IntrRay3Ellipsoid3FIResult {
        const result = defaultIntrRay3Ellipsoid3FIResult();
        intrRay3Ellipsoid3FIDoQuery(ray.origin, ray.direction, ellipsoid,
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
