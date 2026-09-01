// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLine3AlignedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a line and a solid aligned box in 3D.
//
// The line is P + t * D, where D is not required to be unit length.
//
// The aligned box has minimum corner A and maximum corner B. A box point is X
// where A <= X <= B; the comparisons are componentwise.
//
// The closest point on the line is stored in closest[0] with parameter t. The
// closest point on the box is stored in closest[1]. When there are infinitely
// many choices for the pair of closest points, only one of them is returned.
//
// The doQueryND functions are described in Section 10.9.4 Linear Component to
// Oriented Bounding Box of
//    Geometric Tools for Computer Graphics,
//    Philip J. Schneider and David H. Eberly,
//    Morgan Kaufmann, San Francisco CA, 2002
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Line3<T>, AlignedBox3<T>>' becomes the
// class DistLine3AlignedBox3 with the result type DistLine3AlignedBox3Result,
// which is structurally identical to the canonical-box result that upstream
// aliases.

import type { AlignedBox3 } from './AlignedBox';
import { CanonicalBox } from './CanonicalBox';
import type { DCPQuery } from './DCPQuery';
import { DistLine3CanonicalBox3 } from './DistLine3CanonicalBox3';
import { Line, type Line3 } from './Line';
import { Vector, add, mul, sub } from './Vector';

export interface DistLine3AlignedBox3Result {
    distance: number;
    sqrDistance: number;

    // The line parameter t of the closest line point.
    parameter: number;

    // closest[0] is on the line, closest[1] is on the box.
    closest: [Vector, Vector];
}

export class DistLine3AlignedBox3
    implements DCPQuery<Line3, AlignedBox3, DistLine3AlignedBox3Result> {
    compute(line: Line3, box: AlignedBox3): DistLine3AlignedBox3Result {
        // Translate the line and box so that the box has center at the
        // origin.
        const { center: boxCenter, extent } = box.getCenteredForm();
        const cbox = CanonicalBox.fromExtent(extent);
        const xfrmOrigin = sub(line.origin, boxCenter);

        // The query computes the result relative to the box with center at
        // the origin.
        const xfrmLine = Line.fromOriginDirection(xfrmOrigin, line.direction);
        const result: DistLine3AlignedBox3Result =
            new DistLine3CanonicalBox3().compute(xfrmLine, cbox);

        // Compute the closest point on the line.
        result.closest[0] = add(line.origin,
            mul(result.parameter, line.direction));

        // Translate the closest box point to the original coordinates.
        result.closest[1] = add(result.closest[1], boxCenter);
        return result;
    }
}
