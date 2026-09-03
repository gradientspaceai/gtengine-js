// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrHalfspace3Cylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Queries for intersection of objects with halfspaces. These are useful for
// containment testing, object culling, and clipping.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, so the port has only IntrHalfspace3Cylinder3TI.
//
// Upstream bug (fixed here): the C++ computes
//   root = std::sqrt(std::max((T)1, (T)1 - absNdW * absNdW))
// but 1 - absNdW^2 is never larger than 1, so the max() always selects 1 and
// the code computes root = 1 regardless of the axis direction. The comment
// immediately above it states the intended formula r*sqrt(1-Dot(N,W)^2), so
// the clamp is meant to guard the square root against a slightly negative
// argument, i.e. std::max((T)0, ...). The upstream expression overestimates
// the projection interval and therefore reports intersections that do not
// exist. The port uses Math.max(0, 1 - absNdW * absNdW).
//
// Upstream bug (guarded here): the C++ uses cylinder.height unconditionally,
// but Cylinder3.h represents an infinite cylinder with the sentinel
// height = -1. Feeding that sentinel into the finite formula silently
// produces a wrong (too small) projection interval. The port asserts that the
// cylinder is finite, using the same message as the upstream assertion in
// IntrCanonicalBox3Cylinder3.h. See also upstream issue #187.

import type { Cylinder3 } from './Cylinder3.js';
import type { Halfspace } from './Halfspace.js';
import { logAssert } from './Logger.js';
import { dot } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrHalfspace3Cylinder3TI.test.
export interface IntrHalfspace3Cylinder3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrHalfspace3Cylinder3TIResult {
    return { intersect: false };
}

export class IntrHalfspace3Cylinder3TI implements
    TIQuery<Halfspace, Cylinder3, IntrHalfspace3Cylinder3TIResult> {

    test(halfspace: Halfspace, cylinder: Cylinder3):
        IntrHalfspace3Cylinder3TIResult {
        logAssert(halfspace.dimension === 3,
            'IntrHalfspace3Cylinder3TI: mismatched sizes.');
        logAssert(cylinder.isFinite(),
            'Infinite cylinders are not yet supported.');

        const result = defaultTIResult();

        // Compute extremes of signed distance Dot(N,X)-d for points on the
        // cylinder. These are
        //   min = (Dot(N,C)-d) - r*sqrt(1-Dot(N,W)^2) - (h/2)*|Dot(N,W)|
        //   max = (Dot(N,C)-d) + r*sqrt(1-Dot(N,W)^2) + (h/2)*|Dot(N,W)|
        const center = dot(halfspace.normal, cylinder.axis.origin)
            - halfspace.constant;
        const absNdW = Math.abs(dot(halfspace.normal,
            cylinder.axis.direction));
        const root = Math.sqrt(Math.max(0, 1 - absNdW * absNdW));
        const tmax = center + cylinder.radius * root
            + 0.5 * cylinder.height * absNdW;

        // The cylinder and halfspace intersect when the projection interval
        // maximum is nonnegative.
        result.intersect = (tmax >= 0);
        return result;
    }
}
