// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrHalfspace3OrientedBox3.h
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
import { OrientedBox } from './OrientedBox.js';
import { dot } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrHalfspace3OrientedBox3TI.test.
export interface IntrHalfspace3OrientedBox3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrHalfspace3OrientedBox3TIResult {
    return { intersect: false };
}

export class IntrHalfspace3OrientedBox3TI implements
    TIQuery<Halfspace, OrientedBox, IntrHalfspace3OrientedBox3TIResult> {

    test(halfspace: Halfspace, box: OrientedBox): IntrHalfspace3OrientedBox3TIResult {
        const result = defaultTIResult();

        // Project the box center onto the normal line. The plane of the
        // halfspace occurs at the origin (zero) of the normal line.
        const center = dot(halfspace.normal, box.center) - halfspace.constant;

        // Compute the radius of the interval of projection.
        const radius =
            Math.abs(box.extent.values[0] * dot(halfspace.normal, box.axis[0])) +
            Math.abs(box.extent.values[1] * dot(halfspace.normal, box.axis[1])) +
            Math.abs(box.extent.values[2] * dot(halfspace.normal, box.axis[2]));

        // The box and halfspace intersect when the projection interval
        // maximum is nonnegative.
        result.intersect = (center + radius >= 0);
        return result;
    }
}
