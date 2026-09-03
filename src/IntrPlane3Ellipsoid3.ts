// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrPlane3Ellipsoid3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query considers the ellipsoid to be a solid.
//
// The ellipsoid is (X-K)^T*M*(X-K) = 1 and the plane is Dot(N,X) = c. The
// plane intersects the ellipsoid when the distance from the ellipsoid center
// to the plane is no larger than sqrt(N^T*M^{-1}*N)/|N|. Upstream compares
// the (already |N|-normalized) point-plane distance to sqrt(N^T*M^{-1}*N),
// which requires the plane normal to be unit length.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, which becomes IntrPlane3Ellipsoid3TI.
// 'Ellipsoid3::GetMInverse(M)' returns its matrix in the port.

import { DistPointHyperplane } from './DistPointHyperplane.js';
import type { Ellipsoid3 } from './Hyperellipsoid.js';
import { logAssert } from './Logger.js';
import { mulMatrix } from './Matrix.js';
import type { Plane3 } from './Hyperplane.js';
import { Vector, dot } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrPlane3Ellipsoid3TI.test.
export interface IntrPlane3Ellipsoid3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrPlane3Ellipsoid3TIResult():
    IntrPlane3Ellipsoid3TIResult {
    return { intersect: false };
}

// Test-intersection query for a plane and a solid ellipsoid in 3D.
export class IntrPlane3Ellipsoid3TI implements
    TIQuery<Plane3, Ellipsoid3, IntrPlane3Ellipsoid3TIResult> {

    test(plane: Plane3, ellipsoid: Ellipsoid3): IntrPlane3Ellipsoid3TIResult {
        logAssert(plane.dimension === 3 && ellipsoid.dimension === 3,
            'IntrPlane3Ellipsoid3TI: mismatched sizes.');
        const result = defaultIntrPlane3Ellipsoid3TIResult();
        const MInverse = ellipsoid.getMInverse();
        const discr = dot(plane.normal,
            mulMatrix(MInverse, plane.normal) as Vector);
        const root = Math.sqrt(Math.max(discr, 0));
        const vpQuery = new DistPointHyperplane();
        const distance = vpQuery.compute(ellipsoid.center, plane).distance;
        result.intersect = (distance <= root);
        return result;
    }
}
