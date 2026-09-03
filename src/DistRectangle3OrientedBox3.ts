// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRectangle3OrientedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a rectangle and a solid oriented box in 3D.
//
// The rectangle has center C, unit-length axis directions W[0] and W[1], and
// extents e[0] and e[1]. A rectangle point is X = C + sum_{i=0}^1 s[i] * W[i]
// where |s[i]| <= e[i] for all i.
//
// The oriented box has center C, unit-length axis directions U[i] and extents
// e[i] for all i. A box point is X = C + sum_i y[i] * U[i], where
// |y[i]| <= e[i] for all i.
//
// The closest point on the rectangle is stored in closest[0] with
// W-coordinates (s[0],s[1]). The closest point on the box is stored in
// closest[1]. When there are infinitely many choices for the pair of closest
// points, only one of them is returned.
//
// TODO (upstream): Modify to support non-unit-length W[].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Rectangle3<T>, OrientedBox3<T>>'
// becomes the class DistRectangle3OrientedBox3. As upstream does, the result
// type is the rectangle-canonical-box result type, re-exported here as the
// alias DistRectangle3OrientedBox3Result.

import { CanonicalBox } from './CanonicalBox.js';
import type { DCPQuery } from './DCPQuery.js';
import { DistRectangle3CanonicalBox3 } from './DistRectangle3CanonicalBox3.js';
import type { DistRectangle3CanonicalBox3Result } from './DistRectangle3CanonicalBox3.js';
import type { OrientedBox3 } from './OrientedBox.js';
import { Rectangle } from './Rectangle.js';
import type { Rectangle3 } from './Rectangle.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

// Upstream reuses the rectangle-canonical-box result type ('using Result =
// typename RBQuery::Result').
export type DistRectangle3OrientedBox3Result =
    DistRectangle3CanonicalBox3Result;

export class DistRectangle3OrientedBox3
    implements DCPQuery<Rectangle3, OrientedBox3,
    DistRectangle3OrientedBox3Result> {
    compute(rectangle: Rectangle3, box: OrientedBox3):
        DistRectangle3OrientedBox3Result {
        // Rotate and translate the rectangle and box so that the box is
        // aligned and has center at the origin.
        const cbox = CanonicalBox.fromExtent(box.extent);
        const delta = sub(rectangle.center, box.center);
        const xfrmCenter = new Vector(3);
        const xfrmAxis: Vector[] = [new Vector(3), new Vector(3)];
        for (let i = 0; i < 3; ++i) {
            xfrmCenter.values[i] = dot(box.axis[i], delta);
            for (let j = 0; j < 2; ++j) {
                xfrmAxis[j].values[i] = dot(box.axis[i], rectangle.axis[j]);
            }
        }

        // The query computes 'result' relative to the box with center at the
        // origin.
        const xfrmRectangle = Rectangle.fromCenterAxisExtent(xfrmCenter,
            xfrmAxis, rectangle.extent);
        const rbQuery = new DistRectangle3CanonicalBox3();
        const result = rbQuery.compute(xfrmRectangle, cbox);

        // Rotate and translate the closest points to the original
        // coordinates.
        const closest: [Vector, Vector] =
            [box.center.clone(), box.center.clone()];
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
