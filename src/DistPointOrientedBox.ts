// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPointOrientedBox.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance from a point to a solid oriented box in nD.
//
// The oriented box has center C, unit-length axis directions U[i] and extents
// e[i] for all i. A box point is X = C + sum_i y[i] * U[i], where
// |y[i]| <= e[i] for all i.
//
// The input point is stored in closest[0]. The closest point on the box is
// stored in closest[1].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector<N,T>, OrientedBox<N,T>>'
// becomes the class DistPointOrientedBox with the result type
// DistPointOrientedBoxResult. The dimension aliases
// DCPPoint2OrientedBox2/DCPPoint3OrientedBox3 are dropped since the
// runtime-dimension Vector serves every N.

import { CanonicalBox } from './CanonicalBox';
import type { DCPQuery } from './DCPQuery';
import { DistPointCanonicalBox } from './DistPointCanonicalBox';
import type { OrientedBox } from './OrientedBox';
import { Vector, add, dot, mul, sub } from './Vector';

export interface DistPointOrientedBoxResult {
    distance: number;
    sqrDistance: number;

    // closest[0] is the input point, closest[1] is the closest box point.
    closest: [Vector, Vector];
}

export class DistPointOrientedBox
    implements DCPQuery<Vector, OrientedBox, DistPointOrientedBoxResult> {
    compute(point: Vector, box: OrientedBox): DistPointOrientedBoxResult {
        // Rotate and translate the point and box so that the box is aligned
        // and has center at the origin.
        const n = box.extent.size;
        const cbox = CanonicalBox.fromExtent(box.extent);
        const delta = sub(point, box.center);
        const xfrmPoint = new Vector(n);
        for (let i = 0; i < n; ++i) {
            xfrmPoint.values[i] = dot(box.axis[i], delta);
        }

        // The query computes the result relative to the box with center at
        // the origin.
        const pcResult = new DistPointCanonicalBox().compute(xfrmPoint, cbox);

        // Rotate and translate the closest box point to the original
        // coordinates.
        let closest1 = box.center.clone();
        for (let i = 0; i < n; ++i) {
            closest1 = add(closest1,
                mul(pcResult.closest[1].values[i], box.axis[i]));
        }

        return {
            distance: pcResult.distance,
            sqrDistance: pcResult.sqrDistance,
            closest: [point.clone(), closest1]
        };
    }
}
