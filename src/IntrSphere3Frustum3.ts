// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSphere3Frustum3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test whether a solid sphere and a solid orthogonal frustum intersect. The
// sphere intersects the frustum if and only if the distance from the sphere
// center to the frustum is at most the sphere radius.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream provides
// only a TIQuery specialization for this pair of primitives, which becomes
// IntrSphere3Frustum3TI with the result type IntrSphere3Frustum3TIResult.

import type { TIQuery } from './TIQuery.js';
import { DistPoint3Frustum3 } from './DistPoint3Frustum3.js';
import type { Frustum3 } from './Frustum3.js';
import type { Sphere3 } from './Hypersphere.js';

export interface IntrSphere3Frustum3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSphere3Frustum3TIResult():
    IntrSphere3Frustum3TIResult {
    return { intersect: false };
}

// Test-intersection query for a solid sphere and a solid frustum in 3D.
export class IntrSphere3Frustum3TI implements
    TIQuery<Sphere3, Frustum3, IntrSphere3Frustum3TIResult> {

    test(sphere: Sphere3, frustum: Frustum3): IntrSphere3Frustum3TIResult {
        const result = defaultIntrSphere3Frustum3TIResult();
        const vfQuery = new DistPoint3Frustum3();
        const distance = vfQuery.compute(sphere.center, frustum).distance;
        result.intersect = (distance <= sphere.radius);
        return result;
    }
}
