// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPlane3OrientedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a plane and a solid oriented box in 3D.
//
// The plane is defined by Dot(N, X - P) = 0, where P is the plane origin and
// N is a unit-length normal for the plane.
//
// The oriented box has center C, unit-length axis directions U[i] and extents
// e[i] for all i. A box point is X = C + sum_i y[i] * U[i], where
// |y[i]| <= e[i] for all i.
//
// The closest point on the plane is stored in closest[0]. The closest point
// on the box is stored in closest[1]. When there are infinitely many choices
// for the pair of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Plane3<T>, OrientedBox3<T>>' becomes
// the class DistPlane3OrientedBox3. As upstream does, the result type is the
// plane-canonical-box result type, re-exported here as the alias
// DistPlane3OrientedBox3Result.

import { CanonicalBox } from './CanonicalBox';
import type { DCPQuery } from './DCPQuery';
import { DistPlane3CanonicalBox3 } from './DistPlane3CanonicalBox3';
import type { DistPlane3CanonicalBox3Result } from './DistPlane3CanonicalBox3';
import { Hyperplane } from './Hyperplane';
import type { Plane3 } from './Hyperplane';
import type { OrientedBox3 } from './OrientedBox';
import { Vector, add, dot, mul, sub } from './Vector';

// Upstream reuses the plane-box result type ('using Result = typename
// PCQuery::Result').
export type DistPlane3OrientedBox3Result = DistPlane3CanonicalBox3Result;

export class DistPlane3OrientedBox3
    implements DCPQuery<Plane3, OrientedBox3, DistPlane3OrientedBox3Result> {
    compute(plane: Plane3, box: OrientedBox3): DistPlane3OrientedBox3Result {
        // Rotate and translate the plane and box so that the box is aligned
        // and has center at the origin.
        const cbox = CanonicalBox.fromExtent(box.extent);
        const delta = sub(plane.origin, box.center);
        const xfrmOrigin = new Vector(3);
        const xfrmNormal = new Vector(3);
        for (let i = 0; i < 3; ++i) {
            xfrmOrigin.values[i] = dot(box.axis[i], delta);
            xfrmNormal.values[i] = dot(box.axis[i], plane.normal);
        }

        // The query computes 'result' relative to the box with center at the
        // origin.
        const xfrmPlane = Hyperplane.fromNormalOrigin(xfrmNormal, xfrmOrigin);
        const pcQuery = new DistPlane3CanonicalBox3();
        const result = pcQuery.compute(xfrmPlane, cbox);

        // Rotate and translate the closest points to the original
        // coordinates.
        const closest: [Vector, Vector] = [box.center.clone(),
            box.center.clone()];
        for (let i = 0; i < 2; ++i) {
            for (let j = 0; j < 3; ++j) {
                closest[i] = add(closest[i],
                    mul(result.closest[i].values[j], box.axis[j]));
            }
        }
        result.closest = closest;

        return result;
    }
}
