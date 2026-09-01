// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPointAlignedBox.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance from a point to a solid aligned box in nD.
//
// The aligned box has minimum corner A and maximum corner B. A box point is X
// where A <= X <= B; the comparisons are componentwise.
//
// The input point is stored in closest[0]. The closest point on the box is
// stored in closest[1]. When there are infinitely many choices for the pair
// of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector<N,T>, AlignedBox<N,T>>' becomes
// the class DistPointAlignedBox. Upstream aliases its Result to the
// point-canonical-box Result; the port declares the structurally identical
// DistPointAlignedBoxResult so every Dist* file owns its result type. The
// dimension aliases DCPPoint2AlignedBox2/DCPPoint3AlignedBox3 are dropped
// since the runtime-dimension Vector serves every N.

import type { AlignedBox } from './AlignedBox';
import { CanonicalBox } from './CanonicalBox';
import type { DCPQuery } from './DCPQuery';
import { DistPointCanonicalBox } from './DistPointCanonicalBox';
import { Vector, add, sub } from './Vector';

export interface DistPointAlignedBoxResult {
    distance: number;
    sqrDistance: number;

    // closest[0] is the input point, closest[1] is the closest box point.
    closest: [Vector, Vector];
}

export class DistPointAlignedBox
    implements DCPQuery<Vector, AlignedBox, DistPointAlignedBoxResult> {
    compute(point: Vector, box: AlignedBox): DistPointAlignedBoxResult {
        // Translate the point and box so that the box has center at the
        // origin.
        const { center: boxCenter, extent } = box.getCenteredForm();
        const cbox = CanonicalBox.fromExtent(extent);
        const xfrmPoint = sub(point, boxCenter);

        // The query computes the result relative to the box with center at
        // the origin.
        const pcResult = new DistPointCanonicalBox().compute(xfrmPoint, cbox);

        return {
            distance: pcResult.distance,
            sqrDistance: pcResult.sqrDistance,
            closest: [
                // Store the input point.
                point.clone(),
                // Translate the closest box point to the original
                // coordinates.
                add(pcResult.closest[1], boxCenter)
            ]
        };
    }
}
