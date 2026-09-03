// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLine2OrientedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a line and a solid oriented box in 2D.
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
// upstream specialization 'DCPQuery<T, Line2<T>, OrientedBox2<T>>' becomes
// the class DistLine2OrientedBox2 with the result type
// DistLine2OrientedBox2Result, which is structurally identical to the
// aligned-box result that upstream aliases. The aligned-box 'DoQuery' helper
// is the exported distLine2AlignedBox2DoQuery.

import type { DCPQuery } from './DCPQuery.js';
import { distLine2AlignedBox2DoQuery } from './DistLine2AlignedBox2.js';
import type { Line2 } from './Line.js';
import type { OrientedBox2 } from './OrientedBox.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

export interface DistLine2OrientedBox2Result {
    distance: number;
    sqrDistance: number;

    // The line parameter t of the closest line point.
    parameter: number;

    // closest[0] is on the line, closest[1] is on the box.
    closest: [Vector, Vector];
}

export class DistLine2OrientedBox2
    implements DCPQuery<Line2, OrientedBox2, DistLine2OrientedBox2Result> {
    compute(line: Line2, box: OrientedBox2): DistLine2OrientedBox2Result {
        const result: DistLine2OrientedBox2Result = {
            distance: 0,
            sqrDistance: 0,
            parameter: 0,
            closest: [new Vector(2), new Vector(2)]
        };

        // Rotate and translate the line and box so that the box is aligned
        // and has center at the origin.
        const delta = sub(line.origin, box.center);
        const origin = new Vector(2);
        const direction = new Vector(2);
        for (let i = 0; i < 2; ++i) {
            origin.values[i] = dot(box.axis[i], delta);
            direction.values[i] = dot(box.axis[i], line.direction);
        }

        // The query computes 'result' relative to the box with center at the
        // origin.
        distLine2AlignedBox2DoQuery(origin, direction, box.extent, result);

        // Rotate and translate the closest points to the original
        // coordinates.
        const temp = [result.closest[0], result.closest[1]];
        for (let i = 0; i < 2; ++i) {
            result.closest[i] = add(box.center,
                add(mul(temp[i].values[0], box.axis[0]),
                    mul(temp[i].values[1], box.axis[1])));
        }

        // Compute the distance and squared distance.
        const diff = sub(result.closest[0], result.closest[1]);
        result.sqrDistance = dot(diff, diff);
        result.distance = Math.sqrt(result.sqrDistance);
        return result;
    }
}
