// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrPlane3Capsule3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query considers the capsule to be a solid. The signed distances of the
// capsule segment endpoints from the plane determine the answer; the
// comparisons with the capsule radius assume the plane normal is unit length.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, which becomes IntrPlane3Capsule3TI.

import type { Capsule3 } from './Capsule';
import { DistPointHyperplane } from './DistPointHyperplane';
import { logAssert } from './Logger';
import type { Plane3 } from './Hyperplane';
import type { TIQuery } from './TIQuery';

// The result of IntrPlane3Capsule3TI.test.
export interface IntrPlane3Capsule3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrPlane3Capsule3TIResult():
    IntrPlane3Capsule3TIResult {
    return { intersect: false };
}

// Test-intersection query for a plane and a solid capsule in 3D.
export class IntrPlane3Capsule3TI implements
    TIQuery<Plane3, Capsule3, IntrPlane3Capsule3TIResult> {

    test(plane: Plane3, capsule: Capsule3): IntrPlane3Capsule3TIResult {
        logAssert(plane.dimension === 3 && capsule.dimension === 3,
            'IntrPlane3Capsule3TI: mismatched sizes.');
        const result = defaultIntrPlane3Capsule3TIResult();

        const vpQuery = new DistPointHyperplane();
        const sdistance0 = vpQuery.compute(capsule.segment.p[0], plane).signedDistance;
        const sdistance1 = vpQuery.compute(capsule.segment.p[1], plane).signedDistance;
        if (sdistance0 * sdistance1 <= 0) {
            // A capsule segment endpoint is on the plane or the two endpoints
            // are on opposite sides of the plane.
            result.intersect = true;
            return result;
        }

        // The endpoints are on the same side of the plane, but the endpoint
        // spheres might intersect the plane.
        result.intersect =
            Math.abs(sdistance0) <= capsule.radius ||
            Math.abs(sdistance1) <= capsule.radius;
        return result;
    }
}
