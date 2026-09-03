// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRectangle3AlignedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a rectangle and a solid aligned box in 3D.
//
// The rectangle has center C, unit-length axis directions W[0] and W[1], and
// extents e[0] and e[1]. A rectangle point is X = C + sum_{i=0}^1 s[i] * W[i]
// where |s[i]| <= e[i] for all i.
//
// The aligned box has minimum corner A and maximum corner B. A box point is X
// where A <= X <= B; the comparisons are componentwise.
//
// The closest point on the rectangle is stored in closest[0] with
// W-coordinates (s[0],s[1]). The closest point on the box is stored in
// closest[1]. When there are infinitely many choices for the pair of closest
// points, only one of them is returned.
//
// TODO (upstream): Modify to support non-unit-length W[].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Rectangle3<T>, AlignedBox3<T>>'
// becomes the class DistRectangle3AlignedBox3. As upstream does, the result
// type is the rectangle-canonical-box result type, re-exported here as the
// alias DistRectangle3AlignedBox3Result.

import type { AlignedBox3 } from './AlignedBox';
import { CanonicalBox } from './CanonicalBox';
import type { DCPQuery } from './DCPQuery';
import { DistRectangle3CanonicalBox3 } from './DistRectangle3CanonicalBox3';
import type { DistRectangle3CanonicalBox3Result } from './DistRectangle3CanonicalBox3';
import { Rectangle } from './Rectangle';
import type { Rectangle3 } from './Rectangle';
import { add, sub } from './Vector';

// Upstream reuses the rectangle-canonical-box result type ('using Result =
// typename RBQuery::Result').
export type DistRectangle3AlignedBox3Result =
    DistRectangle3CanonicalBox3Result;

export class DistRectangle3AlignedBox3
    implements DCPQuery<Rectangle3, AlignedBox3,
    DistRectangle3AlignedBox3Result> {
    compute(rectangle: Rectangle3, box: AlignedBox3):
        DistRectangle3AlignedBox3Result {
        // Translate the rectangle and box so that the box has center at the
        // origin.
        const centeredForm = box.getCenteredForm();
        const boxCenter = centeredForm.center;
        const cbox = CanonicalBox.fromExtent(centeredForm.extent);
        const xfrmCenter = sub(rectangle.center, boxCenter);

        // The query computes 'result' relative to the box with center at the
        // origin.
        const xfrmRectangle = Rectangle.fromCenterAxisExtent(xfrmCenter,
            rectangle.axis, rectangle.extent);
        const rbQuery = new DistRectangle3CanonicalBox3();
        const result = rbQuery.compute(xfrmRectangle, cbox);

        // Translate the closest points to the original coordinates.
        result.closest[0] = add(result.closest[0], boxCenter);
        result.closest[1] = add(result.closest[1], boxCenter);

        return result;
    }
}
