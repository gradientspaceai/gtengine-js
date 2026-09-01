// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRay3CanonicalBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a ray and a solid canonical box in 3D.
//
// The ray is P + t * D for t >= 0, where D is not required to be unit length.
//
// The canonical box has center at the origin and is aligned with the
// coordinate axes. The extents are E = (e[0],e[1],e[2]). A box point is
// Y = (y[0],y[1],y[2]) with |y[i]| <= e[i] for all i.
//
// The closest point on the ray is stored in closest[0] with parameter t. The
// closest point on the box is stored in closest[1]. When there are infinitely
// many choices for the pair of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Ray3<T>, CanonicalBox3<T>>' becomes
// the class DistRay3CanonicalBox3 with the result type
// DistRay3CanonicalBox3Result, which is structurally identical to the
// line-canonical-box result that upstream aliases.

import type { CanonicalBox3 } from './CanonicalBox';
import type { DCPQuery } from './DCPQuery';
import { DistLine3CanonicalBox3 } from './DistLine3CanonicalBox3';
import { DistPointCanonicalBox } from './DistPointCanonicalBox';
import { Line } from './Line';
import type { Ray3 } from './Ray';
import type { Vector } from './Vector';

export interface DistRay3CanonicalBox3Result {
    distance: number;
    sqrDistance: number;

    // The ray parameter t of the closest ray point.
    parameter: number;

    // closest[0] is on the ray, closest[1] is on the box.
    closest: [Vector, Vector];
}

export class DistRay3CanonicalBox3
    implements DCPQuery<Ray3, CanonicalBox3, DistRay3CanonicalBox3Result> {
    compute(ray: Ray3, box: CanonicalBox3): DistRay3CanonicalBox3Result {
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const lbResult = new DistLine3CanonicalBox3().compute(line, box);
        if (lbResult.parameter >= 0) {
            return lbResult;
        }

        // The closest line point has a negative parameter, so the closest ray
        // point is the ray origin.
        const pbResult = new DistPointCanonicalBox().compute(ray.origin, box);
        return {
            distance: pbResult.distance,
            sqrDistance: pbResult.sqrDistance,
            parameter: 0,
            closest: [ray.origin.clone(), pbResult.closest[1]]
        };
    }
}
