// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay3Cylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query considers the cylinder to be a solid.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream provides
// only an FIQuery specialization for this pair of primitives, which becomes
// IntrRay3Cylinder3FI, and derives it from the Line3-vs-Cylinder3 query to
// reuse the protected DoQuery helper; as in IntrSegment3Cylinder3, the port
// reaches that helper through a module-private accessor subclass. The
// ray-specific DoQuery helper is exported as the module function
// 'intrRay3Cylinder3FIDoQuery'.
//
// As in IntrLine3Cylinder3 and IntrSegment3Cylinder3, upstream has no
// infinite-cylinder branch: it reads cylinder.height directly, so an infinite
// cylinder (height = -1) would produce meaningless results. Following the
// precedent set by IntrHalfspace3Cylinder3 and IntrTriangle3Cylinder3, the
// port asserts that the cylinder is finite instead of silently computing
// nonsense.

import type { Cylinder3 } from './Cylinder3.js';
import type { FIQuery } from './FIQuery.js';
import { IntrIntervalsFI } from './IntrIntervals.js';
import {
    IntrLine3Cylinder3FI,
    defaultIntrLine3Cylinder3FIResult
} from './IntrLine3Cylinder3.js';
import type { IntrLine3Cylinder3FIResult } from './IntrLine3Cylinder3.js';
import { logAssert } from './Logger.js';
import type { Ray3 } from './Ray.js';
import { Vector, add, mul } from './Vector.js';

// The upstream derived FIQuery::Result adds no members.
export type IntrRay3Cylinder3FIResult = IntrLine3Cylinder3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrRay3Cylinder3FIResult():
    IntrRay3Cylinder3FIResult {
    return defaultIntrLine3Cylinder3FIResult();
}

// Expose the protected line-cylinder helper to this module.
class LineCylinderFIAccess extends IntrLine3Cylinder3FI {
    run(lineOrigin: Vector, lineDirection: Vector, cylinder: Cylinder3,
        result: IntrLine3Cylinder3FIResult): void {
        this.doQuery(lineOrigin, lineDirection, cylinder, result);
    }
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrRay3Cylinder3FIDoQuery(rayOrigin: Vector,
    rayDirection: Vector, cylinder: Cylinder3,
    result: IntrRay3Cylinder3FIResult): void {
    new LineCylinderFIAccess().run(rayOrigin, rayDirection, cylinder, result);

    if (result.intersect) {
        // The line containing the ray intersects the cylinder; the t-interval
        // is [t0,t1]. The ray intersects the cylinder as long as [t0,t1]
        // overlaps the ray t-interval [0,+infinity).
        const iiQuery = new IntrIntervalsFI();
        const iiResult = iiQuery.findFiniteSemiInfinite(result.parameter, 0,
            true);
        if (iiResult.intersect) {
            result.numIntersections = iiResult.numIntersections;
            result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
        }
        else {
            // The ray does not intersect the cylinder.
            const empty = defaultIntrRay3Cylinder3FIResult();
            result.intersect = empty.intersect;
            result.numIntersections = empty.numIntersections;
            result.parameter = empty.parameter;
            result.point = empty.point;
        }
    }
}

// Find-intersection query for a ray and a solid cylinder in 3D.
export class IntrRay3Cylinder3FI implements
    FIQuery<Ray3, Cylinder3, IntrRay3Cylinder3FIResult> {

    find(ray: Ray3, cylinder: Cylinder3): IntrRay3Cylinder3FIResult {
        logAssert(cylinder.isFinite(),
            'Infinite cylinders are not yet supported.');

        const result = defaultIntrRay3Cylinder3FIResult();
        intrRay3Cylinder3FIDoQuery(ray.origin, ray.direction, cylinder,
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
