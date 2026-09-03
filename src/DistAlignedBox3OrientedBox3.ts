// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistAlignedBox3OrientedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between solid aligned and oriented boxes in 3D.
//
// The aligned box has minimum corner A and maximum corner B. A box point is X
// where A <= X <= B; the comparisons are componentwise.
//
// The oriented box has center C, unit-length axis directions U[i] and extents
// e[i] for all i. A box point is X = C + sum_i y[i] * U[i], where
// |y[i]| <= e[i] for all i.
//
// The closest point of the aligned box is stored in closest[0]. The closest
// point of the oriented box is stored in closest[1].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, AlignedBox3<T>, OrientedBox3<T>>'
// becomes the class DistAlignedBox3OrientedBox3. As upstream does, the
// result type is the box-box result type, re-exported here as the alias
// DistAlignedBox3OrientedBox3Result.

import type { AlignedBox3 } from './AlignedBox';
import type { DCPQuery } from './DCPQuery';
import { DistOrientedBox3OrientedBox3 } from './DistOrientedBox3OrientedBox3';
import type { DistOrientedBox3OrientedBox3Result } from './DistOrientedBox3OrientedBox3';
import { OrientedBox } from './OrientedBox';
import type { OrientedBox3 } from './OrientedBox';
import { Vector, add, mul, sub } from './Vector';

// Upstream reuses the box-box result type ('using Result = typename
// BBQuery::Result').
export type DistAlignedBox3OrientedBox3Result =
    DistOrientedBox3OrientedBox3Result;

export class DistAlignedBox3OrientedBox3
    implements DCPQuery<AlignedBox3, OrientedBox3,
    DistAlignedBox3OrientedBox3Result> {
    compute(box0: AlignedBox3, box1: OrientedBox3):
        DistAlignedBox3OrientedBox3Result {
        // Convert the aligned box to an oriented box.
        const half = 0.5;
        const obox0 = new OrientedBox(3);
        obox0.center = mul(add(box0.max, box0.min), half);
        obox0.extent = mul(sub(box0.max, box0.min), half);
        obox0.axis = [Vector.unit(3, 0), Vector.unit(3, 1), Vector.unit(3, 2)];

        // Execute the query for two oriented boxes.
        const bbQuery = new DistOrientedBox3OrientedBox3();
        return bbQuery.compute(obox0, box1);
    }
}
