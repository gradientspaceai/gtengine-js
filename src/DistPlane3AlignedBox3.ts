// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPlane3AlignedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a plane and a solid aligned box in 3D.
//
// The plane is defined by Dot(N, X - P) = 0, where P is the plane origin and
// N is a unit-length normal for the plane.
//
// The aligned box has minimum corner A and maximum corner B. A box point is X
// where A <= X <= B; the comparisons are componentwise.
//
// The closest point on the plane is stored in closest[0]. The closest point
// on the box is stored in closest[1]. When there are infinitely many choices
// for the pair of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Plane3<T>, AlignedBox3<T>>' becomes
// the class DistPlane3AlignedBox3. As upstream does, the result type is the
// plane-canonical-box result type, re-exported here as the alias
// DistPlane3AlignedBox3Result.

import type { AlignedBox3 } from './AlignedBox.js';
import { CanonicalBox } from './CanonicalBox.js';
import type { DCPQuery } from './DCPQuery.js';
import { DistPlane3CanonicalBox3 } from './DistPlane3CanonicalBox3.js';
import type { DistPlane3CanonicalBox3Result } from './DistPlane3CanonicalBox3.js';
import { Hyperplane } from './Hyperplane.js';
import type { Plane3 } from './Hyperplane.js';
import { add, sub } from './Vector.js';

// Upstream reuses the plane-box result type ('using Result = typename
// PCQuery::Result').
export type DistPlane3AlignedBox3Result = DistPlane3CanonicalBox3Result;

export class DistPlane3AlignedBox3
    implements DCPQuery<Plane3, AlignedBox3, DistPlane3AlignedBox3Result> {
    compute(plane: Plane3, box: AlignedBox3): DistPlane3AlignedBox3Result {
        // Translate the plane and box so that the box has center at the
        // origin.
        const centeredForm = box.getCenteredForm();
        const boxCenter = centeredForm.center;
        const cbox = CanonicalBox.fromExtent(centeredForm.extent);
        const xfrmOrigin = sub(plane.origin, boxCenter);

        // The query computes 'result' relative to the box with center at the
        // origin.
        const xfrmPlane = Hyperplane.fromNormalOrigin(plane.normal, xfrmOrigin);
        const pcQuery = new DistPlane3CanonicalBox3();
        const result = pcQuery.compute(xfrmPlane, cbox);

        // Translate the closest points to the original coordinates.
        result.closest[0] = add(result.closest[0], boxCenter);
        result.closest[1] = add(result.closest[1], boxCenter);

        return result;
    }
}
