// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistTriangle3AlignedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a solid triangle and a solid aligned box
// in 3D.
//
// The triangle has vertices <V[0],V[1],V[2]>. A triangle point is
// X = sum_{i=0}^2 b[i] * V[i], where 0 <= b[i] <= 1 for all i and
// sum_{i=0}^2 b[i] = 1.
//
// The aligned box has minimum corner A and maximum corner B. A box point is X
// where A <= X <= B; the comparisons are componentwise.
//
// The closest point on the triangle is stored in closest[0] with barycentric
// coordinates (b[0],b[1],b[2]). The closest point on the box is stored in
// closest[1]. When there are infinitely many choices for the pair of closest
// points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Triangle3<T>, AlignedBox3<T>>'
// becomes the class DistTriangle3AlignedBox3. As upstream does, the result
// type is the triangle-canonical-box result type, re-exported here as the
// alias DistTriangle3AlignedBox3Result.

import type { AlignedBox3 } from './AlignedBox.js';
import { CanonicalBox } from './CanonicalBox.js';
import type { DCPQuery } from './DCPQuery.js';
import { DistTriangle3CanonicalBox3 } from './DistTriangle3CanonicalBox3.js';
import type { DistTriangle3CanonicalBox3Result } from './DistTriangle3CanonicalBox3.js';
import { Triangle } from './Triangle.js';
import type { Triangle3 } from './Triangle.js';
import { add, sub } from './Vector.js';

// Upstream reuses the triangle-canonical-box result type ('using Result =
// typename TBQuery::Result').
export type DistTriangle3AlignedBox3Result =
    DistTriangle3CanonicalBox3Result;

export class DistTriangle3AlignedBox3
    implements DCPQuery<Triangle3, AlignedBox3,
    DistTriangle3AlignedBox3Result> {
    compute(triangle: Triangle3, box: AlignedBox3):
        DistTriangle3AlignedBox3Result {
        // Translate the triangle and box so that the box has center at the
        // origin.
        const centeredForm = box.getCenteredForm();
        const boxCenter = centeredForm.center;
        const cbox = CanonicalBox.fromExtent(centeredForm.extent);
        const xfrmTriangle = Triangle.fromVertices(
            sub(triangle.v[0], boxCenter),
            sub(triangle.v[1], boxCenter),
            sub(triangle.v[2], boxCenter));

        // The query computes 'result' relative to the box with center at the
        // origin.
        const tbQuery = new DistTriangle3CanonicalBox3();
        const result = tbQuery.compute(xfrmTriangle, cbox);

        // Translate the closest points to the original coordinates.
        result.closest[0] = add(result.closest[0], boxCenter);
        result.closest[1] = add(result.closest[1], boxCenter);

        return result;
    }
}
