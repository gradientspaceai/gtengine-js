// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay2AlignedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the box to be a solid.
//
// The test-intersection queries use the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection queries use parametric clipping against the four
// edges of the box.
//
// Port notes (see IntrIntervals.ts for the Intr* precedent): upstream derives
// the ray queries from the line queries only to reuse the protected DoQuery
// members. The line-box helpers are the exported module functions
// 'intrLine2AlignedBox2TIDoQuery' and 'intrLine2AlignedBox2FIDoQuery', and
// the ray-specific ones are exported here as 'intrRay2AlignedBox2TIDoQuery'
// and 'intrRay2AlignedBox2FIDoQuery' for the same reason. The upstream Result
// structs add no members to the line results, so the port exports type
// aliases.

import type { AlignedBox } from './AlignedBox.js';
import type { Ray } from './Ray.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import { Vector, add, mul, sub } from './Vector.js';
import {
    intrLine2AlignedBox2TIDoQuery,
    intrLine2AlignedBox2FIDoQuery,
    defaultIntrLine2AlignedBox2TIResult,
    defaultIntrLine2AlignedBox2FIResult
} from './IntrLine2AlignedBox2.js';
import type {
    IntrLine2AlignedBox2TIResult,
    IntrLine2AlignedBox2FIResult
} from './IntrLine2AlignedBox2.js';
import { IntrIntervalsFI } from './IntrIntervals.js';

// The port of std::numeric_limits<T>::max() for T = double.
const MAX_T = Number.MAX_VALUE;

// The result of IntrRay2AlignedBox2TI.test. Upstream adds no members to the
// line-box result.
export type IntrRay2AlignedBox2TIResult = IntrLine2AlignedBox2TIResult;

// The result of IntrRay2AlignedBox2FI.find. Upstream adds no members to the
// line-box result.
export type IntrRay2AlignedBox2FIResult = IntrLine2AlignedBox2FIResult;

// The port of the protected 'TIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed. The 'result' values are modified
// in place.
export function intrRay2AlignedBox2TIDoQuery(rayOrigin: Vector,
    rayDirection: Vector, boxExtent: Vector,
    result: IntrRay2AlignedBox2TIResult): void {
    const o = rayOrigin.values;
    const d = rayDirection.values;
    const e = boxExtent.values;
    for (let i = 0; i < 2; ++i) {
        if (Math.abs(o[i]) > e[i] && o[i] * d[i] >= 0) {
            result.intersect = false;
            return;
        }
    }

    intrLine2AlignedBox2TIDoQuery(rayOrigin, rayDirection, boxExtent, result);
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed. The 'result' values are modified
// in place.
export function intrRay2AlignedBox2FIDoQuery(rayOrigin: Vector,
    rayDirection: Vector, boxExtent: Vector,
    result: IntrRay2AlignedBox2FIResult): void {
    intrLine2AlignedBox2FIDoQuery(rayOrigin, rayDirection, boxExtent, result);

    if (result.intersect) {
        // The line containing the ray intersects the box; the t-interval
        // is [t0,t1]. The ray intersects the box as long as [t0,t1]
        // overlaps the ray t-interval [0,+infinity).
        const rayInterval: [number, number] = [0, MAX_T];
        const iiQuery = new IntrIntervalsFI();
        const iiResult = iiQuery.find(result.parameter, rayInterval);
        result.intersect = iiResult.intersect;
        result.numIntersections = iiResult.numIntersections;
        result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
    }
}

// Test-intersection query for a ray and a solid aligned box in 2D.
export class IntrRay2AlignedBox2TI implements
    TIQuery<Ray, AlignedBox, IntrRay2AlignedBox2TIResult> {

    test(ray: Ray, box: AlignedBox): IntrRay2AlignedBox2TIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector2::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        // Transform the ray to the aligned-box coordinate system.
        const rayOrigin = sub(ray.origin, boxCenter);

        const result = defaultIntrLine2AlignedBox2TIResult();
        intrRay2AlignedBox2TIDoQuery(rayOrigin, ray.direction, boxExtent,
            result);
        return result;
    }
}

// Find-intersection query for a ray and a solid aligned box in 2D.
export class IntrRay2AlignedBox2FI implements
    FIQuery<Ray, AlignedBox, IntrRay2AlignedBox2FIResult> {

    find(ray: Ray, box: AlignedBox): IntrRay2AlignedBox2FIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector2::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        // Transform the ray to the aligned-box coordinate system.
        const rayOrigin = sub(ray.origin, boxCenter);

        const result = defaultIntrLine2AlignedBox2FIResult();
        intrRay2AlignedBox2FIDoQuery(rayOrigin, ray.direction, boxExtent,
            result);
        for (let i = 0; i < result.numIntersections; ++i) {
            result.point[i] = add(ray.origin,
                mul(result.parameter[i], ray.direction));
        }
        return result;
    }
}
