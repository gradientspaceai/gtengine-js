// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrHalfspace3Sphere3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Queries for intersection of objects with halfspaces. These are useful for
// containment testing, object culling, and clipping.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream provides
// only a TIQuery specialization for this pair of primitives.

import { Halfspace } from './Halfspace.js';
import { Hypersphere } from './Hypersphere.js';
import { dot } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrHalfspace3Sphere3TI.test.
export interface IntrHalfspace3Sphere3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrHalfspace3Sphere3TIResult {
    return { intersect: false };
}

export class IntrHalfspace3Sphere3TI implements
    TIQuery<Halfspace, Hypersphere, IntrHalfspace3Sphere3TIResult> {

    test(halfspace: Halfspace, sphere: Hypersphere): IntrHalfspace3Sphere3TIResult {
        const result = defaultTIResult();

        // Project the sphere center onto the normal line. The plane of the
        // halfspace occurs at the origin (zero) of the normal line.
        const center = dot(halfspace.normal, sphere.center) - halfspace.constant;

        // The sphere and halfspace intersect when the projection interval
        // maximum is nonnegative.
        result.intersect = (center + sphere.radius >= 0);
        return result;
    }
}
