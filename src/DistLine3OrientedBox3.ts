// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLine3OrientedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a line and a solid oriented box in 3D.
//
// The line is P + t * D, where D is not required to be unit length.
//
// The oriented box has center C, unit-length axis directions U[i] and extents
// e[i] for all i. A box point is X = C + sum_i y[i] * U[i], where
// |y[i]| <= e[i] for all i.
//
// The closest point on the line is stored in closest[0] with parameter t. The
// closest point on the box is stored in closest[1]. When there are infinitely
// many choices for the pair of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Line3<T>, OrientedBox3<T>>' becomes
// the class DistLine3OrientedBox3 with the result type
// DistLine3OrientedBox3Result, which is structurally identical to the
// canonical-box result that upstream aliases.

import { CanonicalBox } from './CanonicalBox.js';
import type { DCPQuery } from './DCPQuery.js';
import { DistLine3CanonicalBox3 } from './DistLine3CanonicalBox3.js';
import { Line, type Line3 } from './Line.js';
import type { OrientedBox3 } from './OrientedBox.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

export interface DistLine3OrientedBox3Result {
    distance: number;
    sqrDistance: number;

    // The line parameter t of the closest line point.
    parameter: number;

    // closest[0] is on the line, closest[1] is on the box.
    closest: [Vector, Vector];
}

export class DistLine3OrientedBox3
    implements DCPQuery<Line3, OrientedBox3, DistLine3OrientedBox3Result> {
    compute(line: Line3, box: OrientedBox3): DistLine3OrientedBox3Result {
        // Rotate and translate the line and box so that the box is aligned
        // and has center at the origin.
        const cbox = CanonicalBox.fromExtent(box.extent);
        const delta = sub(line.origin, box.center);
        const xfrmOrigin = new Vector(3);
        const xfrmDirection = new Vector(3);
        for (let i = 0; i < 3; ++i) {
            xfrmOrigin.values[i] = dot(box.axis[i], delta);
            xfrmDirection.values[i] = dot(box.axis[i], line.direction);
        }

        // The query computes the result relative to the box with center at
        // the origin.
        const xfrmLine = Line.fromOriginDirection(xfrmOrigin, xfrmDirection);
        const result: DistLine3OrientedBox3Result =
            new DistLine3CanonicalBox3().compute(xfrmLine, cbox);

        // Rotate and translate the closest points to the original
        // coordinates. Upstream first assigns closest[0] to
        // line.origin + t * line.direction and then overwrites both closest
        // points with the rotated ones; the two expressions produce the same
        // line point, so the port omits the redundant assignment.
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
