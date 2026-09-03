// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay3AlignedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The test-intersection queries use the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection queries use parametric clipping against the six
// faces of the box (Liang-Barsky clipping). The queries consider the box to
// be a solid. The algorithms are described in
// https://www.geometrictools.com/Documentation/IntersectionLineBox.pdf
//
// Port notes (see IntrIntervals.ts for the Intr* precedent): upstream derives
// the ray queries from the line queries only to reuse the protected DoQuery
// members. The line-box helpers are the exported module functions
// 'intrLine3AlignedBox3TIDoQuery' and 'intrLine3AlignedBox3FIDoQuery', and
// the ray-specific ones are exported here as 'intrRay3AlignedBox3TIDoQuery'
// and 'intrRay3AlignedBox3FIDoQuery' for the same reason. The upstream Result
// structs add no members to the line results, so the port exports type
// aliases.

import type { AlignedBox } from './AlignedBox.js';
import type { Ray } from './Ray.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import { Vector, add, mul, sub } from './Vector.js';
import {
    intrLine3AlignedBox3TIDoQuery,
    intrLine3AlignedBox3FIDoQuery,
    defaultIntrLine3AlignedBox3TIResult,
    defaultIntrLine3AlignedBox3FIResult
} from './IntrLine3AlignedBox3.js';
import type {
    IntrLine3AlignedBox3TIResult,
    IntrLine3AlignedBox3FIResult
} from './IntrLine3AlignedBox3.js';
import { IntrIntervalsFI } from './IntrIntervals.js';

// The result of IntrRay3AlignedBox3TI.test. Upstream adds no members to the
// line-box result.
export type IntrRay3AlignedBox3TIResult = IntrLine3AlignedBox3TIResult;

// The result of IntrRay3AlignedBox3FI.find. Upstream adds no members to the
// line-box result.
export type IntrRay3AlignedBox3FIResult = IntrLine3AlignedBox3FIResult;

// The port of the protected 'TIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrRay3AlignedBox3TIDoQuery(rayOrigin: Vector,
    rayDirection: Vector, boxExtent: Vector,
    result: IntrRay3AlignedBox3TIResult): void {
    const o = rayOrigin.values;
    const d = rayDirection.values;
    const e = boxExtent.values;
    for (let i = 0; i < 3; ++i) {
        if (Math.abs(o[i]) > e[i] && o[i] * d[i] >= 0) {
            result.intersect = false;
            return;
        }
    }

    intrLine3AlignedBox3TIDoQuery(rayOrigin, rayDirection, boxExtent, result);
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrRay3AlignedBox3FIDoQuery(rayOrigin: Vector,
    rayDirection: Vector, boxExtent: Vector,
    result: IntrRay3AlignedBox3FIResult): void {
    intrLine3AlignedBox3FIDoQuery(rayOrigin, rayDirection, boxExtent, result);

    if (result.intersect) {
        // The line containing the ray intersects the box; the t-interval
        // is [t0,t1]. The ray intersects the box as long as [t0,t1]
        // overlaps the ray t-interval [0,+infinity).
        const iiQuery = new IntrIntervalsFI();
        const iiResult = iiQuery.findFiniteSemiInfinite(result.parameter,
            0, true);
        if (iiResult.intersect) {
            result.numIntersections = iiResult.numIntersections;
            result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
        }
        else {
            // The ray does not intersect the box.
            result.intersect = false;
            result.numIntersections = 0;
            result.parameter = [0, 0];
            result.point = [Vector.zero(3), Vector.zero(3)];
        }
    }
}

// Test-intersection query for a ray and a solid aligned box in 3D.
export class IntrRay3AlignedBox3TI implements
    TIQuery<Ray, AlignedBox, IntrRay3AlignedBox3TIResult> {

    test(ray: Ray, box: AlignedBox): IntrRay3AlignedBox3TIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector3::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        // Transform the ray to the aligned-box coordinate system.
        const rayOrigin = sub(ray.origin, boxCenter);

        const result = defaultIntrLine3AlignedBox3TIResult();
        intrRay3AlignedBox3TIDoQuery(rayOrigin, ray.direction, boxExtent,
            result);
        return result;
    }
}

// Find-intersection query for a ray and a solid aligned box in 3D.
export class IntrRay3AlignedBox3FI implements
    FIQuery<Ray, AlignedBox, IntrRay3AlignedBox3FIResult> {

    find(ray: Ray, box: AlignedBox): IntrRay3AlignedBox3FIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector3::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        // Transform the ray to the aligned-box coordinate system.
        const rayOrigin = sub(ray.origin, boxCenter);

        const result = defaultIntrLine3AlignedBox3FIResult();
        intrRay3AlignedBox3FIDoQuery(rayOrigin, ray.direction, boxExtent,
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
