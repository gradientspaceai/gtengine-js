// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrHalfspace3Capsule3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Queries for intersection of objects with halfspaces. These are useful for
// containment testing, object culling, and clipping.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, so the port has only IntrHalfspace3Capsule3TI.

import type { Capsule } from './Capsule.js';
import type { Halfspace } from './Halfspace.js';
import { logAssert } from './Logger.js';
import { dot } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrHalfspace3Capsule3TI.test.
export interface IntrHalfspace3Capsule3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrHalfspace3Capsule3TIResult {
    return { intersect: false };
}

export class IntrHalfspace3Capsule3TI implements
    TIQuery<Halfspace, Capsule, IntrHalfspace3Capsule3TIResult> {

    test(halfspace: Halfspace, capsule: Capsule): IntrHalfspace3Capsule3TIResult {
        logAssert(halfspace.dimension === 3 && capsule.dimension === 3,
            'IntrHalfspace3Capsule3TI: mismatched sizes.');

        const result = defaultTIResult();

        // Project the capsule onto the normal line. The plane of the
        // halfspace occurs at the origin (zero) of the normal line.
        const e0 = dot(halfspace.normal, capsule.segment.p[0])
            - halfspace.constant;
        const e1 = dot(halfspace.normal, capsule.segment.p[1])
            - halfspace.constant;

        // The capsule and halfspace intersect when the projection interval
        // maximum is nonnegative.
        result.intersect = (Math.max(e0, e1) + capsule.radius >= 0);
        return result;
    }
}
