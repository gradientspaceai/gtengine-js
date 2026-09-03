// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay3OrientedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The test-intersection queries use the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection queries use parametric clipping against the six faces
// of the box; they use Liang-Barsky clipping. The queries consider the box to
// be a solid. The algorithms are described in
// https://www.geometrictools.com/Documentation/IntersectionLineBox.pdf
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// these queries from the Ray3-vs-AlignedBox3 queries to reuse their protected
// DoQuery helpers, which the port exports as the module functions
// 'intrRay3AlignedBox3TIDoQuery' and 'intrRay3AlignedBox3FIDoQuery'.

import {
    intrRay3AlignedBox3TIDoQuery,
    intrRay3AlignedBox3FIDoQuery
} from './IntrRay3AlignedBox3.js';
import type {
    IntrRay3AlignedBox3TIResult,
    IntrRay3AlignedBox3FIResult
} from './IntrRay3AlignedBox3.js';
import {
    defaultIntrLine3AlignedBox3TIResult,
    defaultIntrLine3AlignedBox3FIResult
} from './IntrLine3AlignedBox3.js';
import { logAssert } from './Logger.js';
import type { OrientedBox3 } from './OrientedBox.js';
import type { Ray3 } from './Ray.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';

// The upstream derived Result structs add no members.
export type IntrRay3OrientedBox3TIResult = IntrRay3AlignedBox3TIResult;
export type IntrRay3OrientedBox3FIResult = IntrRay3AlignedBox3FIResult;

// The ports of the upstream Result default constructors.
export function defaultIntrRay3OrientedBox3TIResult():
    IntrRay3OrientedBox3TIResult {
    return defaultIntrLine3AlignedBox3TIResult();
}

export function defaultIntrRay3OrientedBox3FIResult():
    IntrRay3OrientedBox3FIResult {
    return defaultIntrLine3AlignedBox3FIResult();
}

// Transform the ray to the oriented-box coordinate system.
function transformRay(ray: Ray3, box: OrientedBox3):
    { rayOrigin: Vector, rayDirection: Vector } {
    const diff = sub(ray.origin, box.center);
    const rayOrigin = Vector.zero(3);
    const rayDirection = Vector.zero(3);
    for (let i = 0; i < 3; ++i) {
        rayOrigin.values[i] = dot(diff, box.axis[i]);
        rayDirection.values[i] = dot(ray.direction, box.axis[i]);
    }
    return { rayOrigin, rayDirection };
}

// Test-intersection query for a ray and a solid oriented box in 3D.
export class IntrRay3OrientedBox3TI implements
    TIQuery<Ray3, OrientedBox3, IntrRay3OrientedBox3TIResult> {

    test(ray: Ray3, box: OrientedBox3): IntrRay3OrientedBox3TIResult {
        logAssert(box.dimension === 3,
            'IntrRay3OrientedBox3TI: mismatched sizes.');

        const { rayOrigin, rayDirection } = transformRay(ray, box);

        const result = defaultIntrRay3OrientedBox3TIResult();
        intrRay3AlignedBox3TIDoQuery(rayOrigin, rayDirection, box.extent,
            result);
        return result;
    }
}

// Find-intersection query for a ray and a solid oriented box in 3D.
export class IntrRay3OrientedBox3FI implements
    FIQuery<Ray3, OrientedBox3, IntrRay3OrientedBox3FIResult> {

    find(ray: Ray3, box: OrientedBox3): IntrRay3OrientedBox3FIResult {
        logAssert(box.dimension === 3,
            'IntrRay3OrientedBox3FI: mismatched sizes.');

        const { rayOrigin, rayDirection } = transformRay(ray, box);

        const result = defaultIntrRay3OrientedBox3FIResult();
        intrRay3AlignedBox3FIDoQuery(rayOrigin, rayDirection, box.extent,
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
