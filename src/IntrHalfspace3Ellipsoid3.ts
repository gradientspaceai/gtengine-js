// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrHalfspace3Ellipsoid3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Queries for intersection of objects with halfspaces. These are useful for
// containment testing, object culling, and clipping.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, so the port has only IntrHalfspace3Ellipsoid3TI.

import type { Ellipsoid3 } from './Hyperellipsoid';
import type { Halfspace3 } from './Halfspace';
import { logAssert } from './Logger';
import { mulMatrix } from './Matrix';
import type { Vector } from './Vector';
import { dot } from './Vector';
import type { TIQuery } from './TIQuery';

// The result of IntrHalfspace3Ellipsoid3TI.test.
export interface IntrHalfspace3Ellipsoid3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrHalfspace3Ellipsoid3TIResult {
    return { intersect: false };
}

// Test-intersection query for a halfspace and a solid ellipsoid in 3D.
export class IntrHalfspace3Ellipsoid3TI implements
    TIQuery<Halfspace3, Ellipsoid3, IntrHalfspace3Ellipsoid3TIResult> {

    test(halfspace: Halfspace3, ellipsoid: Ellipsoid3):
        IntrHalfspace3Ellipsoid3TIResult {
        logAssert(halfspace.dimension === 3 && ellipsoid.dimension === 3,
            'IntrHalfspace3Ellipsoid3TI: mismatched sizes.');

        // Project the ellipsoid onto the normal line. The plane of the
        // halfspace occurs at the origin (zero) of the normal line.
        const result = defaultTIResult();
        const mInverse = ellipsoid.getMInverse();
        const discr = dot(halfspace.normal,
            mulMatrix(mInverse, halfspace.normal) as Vector);
        const extent = Math.sqrt(Math.max(discr, 0));
        const center = dot(halfspace.normal, ellipsoid.center)
            - halfspace.constant;
        const tmax = center + extent;

        // The ellipsoid and halfspace intersect when the projection interval
        // maximum is nonnegative.
        result.intersect = (tmax >= 0);
        return result;
    }
}
