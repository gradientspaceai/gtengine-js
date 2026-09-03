// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay2OrientedBox2.h
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
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// these queries from the Ray2-vs-AlignedBox2 queries to reuse their protected
// DoQuery helpers, which the port exports as the module functions
// 'intrRay2AlignedBox2TIDoQuery' and 'intrRay2AlignedBox2FIDoQuery'.

import {
    intrRay2AlignedBox2TIDoQuery,
    intrRay2AlignedBox2FIDoQuery
} from './IntrRay2AlignedBox2.js';
import type {
    IntrRay2AlignedBox2TIResult,
    IntrRay2AlignedBox2FIResult
} from './IntrRay2AlignedBox2.js';
import {
    defaultIntrLine2AlignedBox2TIResult,
    defaultIntrLine2AlignedBox2FIResult
} from './IntrLine2AlignedBox2.js';
import { logAssert } from './Logger.js';
import type { OrientedBox2 } from './OrientedBox.js';
import type { Ray2 } from './Ray.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';

// The upstream derived Result structs add no members.
export type IntrRay2OrientedBox2TIResult = IntrRay2AlignedBox2TIResult;
export type IntrRay2OrientedBox2FIResult = IntrRay2AlignedBox2FIResult;

// The ports of the upstream Result default constructors.
export function defaultIntrRay2OrientedBox2TIResult():
    IntrRay2OrientedBox2TIResult {
    return defaultIntrLine2AlignedBox2TIResult();
}

export function defaultIntrRay2OrientedBox2FIResult():
    IntrRay2OrientedBox2FIResult {
    return defaultIntrLine2AlignedBox2FIResult();
}

// Transform the ray to the oriented-box coordinate system.
function transformRay(ray: Ray2, box: OrientedBox2):
    { rayOrigin: Vector, rayDirection: Vector } {
    const diff = sub(ray.origin, box.center);
    const rayOrigin = Vector.zero(2);
    const rayDirection = Vector.zero(2);
    for (let i = 0; i < 2; ++i) {
        rayOrigin.values[i] = dot(diff, box.axis[i]);
        rayDirection.values[i] = dot(ray.direction, box.axis[i]);
    }
    return { rayOrigin, rayDirection };
}

// Test-intersection query for a ray and a solid oriented box in 2D.
export class IntrRay2OrientedBox2TI implements
    TIQuery<Ray2, OrientedBox2, IntrRay2OrientedBox2TIResult> {

    test(ray: Ray2, box: OrientedBox2): IntrRay2OrientedBox2TIResult {
        logAssert(box.dimension === 2,
            'IntrRay2OrientedBox2TI: mismatched sizes.');

        const { rayOrigin, rayDirection } = transformRay(ray, box);

        const result = defaultIntrRay2OrientedBox2TIResult();
        intrRay2AlignedBox2TIDoQuery(rayOrigin, rayDirection, box.extent,
            result);
        return result;
    }
}

// Find-intersection query for a ray and a solid oriented box in 2D.
export class IntrRay2OrientedBox2FI implements
    FIQuery<Ray2, OrientedBox2, IntrRay2OrientedBox2FIResult> {

    find(ray: Ray2, box: OrientedBox2): IntrRay2OrientedBox2FIResult {
        logAssert(box.dimension === 2,
            'IntrRay2OrientedBox2FI: mismatched sizes.');

        const { rayOrigin, rayDirection } = transformRay(ray, box);

        const result = defaultIntrRay2OrientedBox2FIResult();
        intrRay2AlignedBox2FIDoQuery(rayOrigin, rayDirection, box.extent,
            result);
        for (let i = 0; i < result.numIntersections; ++i) {
            result.point[i] = add(ray.origin,
                mul(result.parameter[i], ray.direction));
        }
        return result;
    }
}
