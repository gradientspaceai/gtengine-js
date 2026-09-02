// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay2Circle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the circle to be a solid (disk).
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// the FIQuery from the Line2-vs-Circle2 query to reuse the protected DoQuery
// helper; as in IntrSegment3Cylinder3, the port reaches that helper through a
// module-private accessor subclass because a TypeScript subclass cannot
// change an inherited method signature. The ray-specific DoQuery helper is
// exported as the module function 'intrRay2Circle2FIDoQuery'.

import type { FIQuery } from './FIQuery';
import type { Hypersphere } from './Hypersphere';
import { IntrIntervalsFI } from './IntrIntervals';
import {
    IntrLine2Circle2FI,
    defaultIntrLine2Circle2FIResult
} from './IntrLine2Circle2';
import type { IntrLine2Circle2FIResult } from './IntrLine2Circle2';
import type { Ray2 } from './Ray';
import { Vector, add, mul } from './Vector';
import type { TIQuery } from './TIQuery';

// The result of IntrRay2Circle2TI.test.
export interface IntrRay2Circle2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrRay2Circle2TIResult {
    return { intersect: false };
}

// The upstream derived FIQuery::Result adds no members.
export type IntrRay2Circle2FIResult = IntrLine2Circle2FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrRay2Circle2FIResult(): IntrRay2Circle2FIResult {
    return defaultIntrLine2Circle2FIResult();
}

// Expose the protected line-circle helper to this module.
class LineCircleFIAccess extends IntrLine2Circle2FI {
    run(lineOrigin: Vector, lineDirection: Vector, circle: Hypersphere,
        result: IntrLine2Circle2FIResult): void {
        this.doQuery(lineOrigin, lineDirection, circle, result);
    }
}

// The port of the protected 'FIQuery::DoQuery'.
export function intrRay2Circle2FIDoQuery(rayOrigin: Vector,
    rayDirection: Vector, circle: Hypersphere,
    result: IntrRay2Circle2FIResult): void {
    new LineCircleFIAccess().run(rayOrigin, rayDirection, circle, result);

    if (result.intersect) {
        // The line containing the ray intersects the disk; the t-interval is
        // [t0,t1]. The ray intersects the disk as long as [t0,t1] overlaps
        // the ray t-interval [0,+infinity).
        const rayInterval: [number, number] = [0, Number.MAX_VALUE];
        const iiQuery = new IntrIntervalsFI();
        const iiResult = iiQuery.find(result.parameter, rayInterval);
        result.intersect = iiResult.intersect;
        result.numIntersections = iiResult.numIntersections;
        result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
    }
}

// Test-intersection query for a ray and a solid circle (disk) in 2D.
export class IntrRay2Circle2TI implements
    TIQuery<Ray2, Hypersphere, IntrRay2Circle2TIResult> {

    test(ray: Ray2, circle: Hypersphere): IntrRay2Circle2TIResult {
        const result = defaultTIResult();
        const rcQuery = new IntrRay2Circle2FI();
        result.intersect = rcQuery.find(ray, circle).intersect;
        return result;
    }
}

// Find-intersection query for a ray and a solid circle (disk) in 2D.
export class IntrRay2Circle2FI implements
    FIQuery<Ray2, Hypersphere, IntrRay2Circle2FIResult> {

    find(ray: Ray2, circle: Hypersphere): IntrRay2Circle2FIResult {
        const result = defaultIntrRay2Circle2FIResult();
        intrRay2Circle2FIDoQuery(ray.origin, ray.direction, circle, result);
        for (let i = 0; i < result.numIntersections; ++i) {
            result.point[i] = add(ray.origin,
                mul(result.parameter[i], ray.direction));
        }
        return result;
    }
}
