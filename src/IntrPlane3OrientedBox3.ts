// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrPlane3OrientedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query considers the box to be a solid. The box projects onto the plane
// normal in an interval of radius
//   r = sum_i |extent[i] * Dot(N, axis[i])|
// centered at the projection of the box center, so the plane and box
// intersect when the distance from the box center to the plane is at most r.
// The comparison assumes the plane normal is unit length.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, which becomes IntrPlane3OrientedBox3TI.

import { DistPointHyperplane } from './DistPointHyperplane';
import { logAssert } from './Logger';
import type { OrientedBox3 } from './OrientedBox';
import type { Plane3 } from './Hyperplane';
import { dot } from './Vector';
import type { TIQuery } from './TIQuery';

// The result of IntrPlane3OrientedBox3TI.test.
export interface IntrPlane3OrientedBox3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrPlane3OrientedBox3TIResult():
    IntrPlane3OrientedBox3TIResult {
    return { intersect: false };
}

// Test-intersection query for a plane and a solid oriented box in 3D.
export class IntrPlane3OrientedBox3TI implements
    TIQuery<Plane3, OrientedBox3, IntrPlane3OrientedBox3TIResult> {

    test(plane: Plane3, box: OrientedBox3): IntrPlane3OrientedBox3TIResult {
        logAssert(plane.dimension === 3 && box.dimension === 3,
            'IntrPlane3OrientedBox3TI: mismatched sizes.');
        const result = defaultIntrPlane3OrientedBox3TIResult();

        const radius =
            Math.abs(box.extent.values[0] * dot(plane.normal, box.axis[0])) +
            Math.abs(box.extent.values[1] * dot(plane.normal, box.axis[1])) +
            Math.abs(box.extent.values[2] * dot(plane.normal, box.axis[2]));

        const ppQuery = new DistPointHyperplane();
        const ppResult = ppQuery.compute(box.center, plane);
        result.intersect = (ppResult.distance <= radius);
        return result;
    }
}
