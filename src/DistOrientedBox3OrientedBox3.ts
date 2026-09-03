// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistOrientedBox3OrientedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between two solid oriented boxes in 3D.
//
// Each oriented box has center C, unit-length axis directions U[i], and
// extents e[i] for all i. A box point is X = C + sum_i y[i] * U[i], where
// |y[i]| <= e[i] for all i.
//
// The closest point of the first oriented box is stored in closest[0]. The
// closest point of the second oriented box is stored in closest[1]. When
// there are infinitely many choices for the pair of closest points, only one
// of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, OrientedBox3<T>, OrientedBox3<T>>'
// becomes the class DistOrientedBox3OrientedBox3 with the result type
// DistOrientedBox3OrientedBox3Result.

import type { DCPQuery } from './DCPQuery.js';
import { DistRectangle3OrientedBox3 } from './DistRectangle3OrientedBox3.js';
import type { OrientedBox3 } from './OrientedBox.js';
import { Rectangle } from './Rectangle.js';
import { Vector, add, mul, sub } from './Vector.js';

export interface DistOrientedBox3OrientedBox3Result {
    distance: number;
    sqrDistance: number;

    // closest[0] is on box0, closest[1] is on box1.
    closest: [Vector, Vector];
}

export class DistOrientedBox3OrientedBox3
    implements DCPQuery<OrientedBox3, OrientedBox3,
    DistOrientedBox3OrientedBox3Result> {
    compute(box0: OrientedBox3, box1: OrientedBox3):
        DistOrientedBox3OrientedBox3Result {
        const result: DistOrientedBox3OrientedBox3Result = {
            distance: 0,
            sqrDistance: 0,
            closest: [new Vector(3), new Vector(3)]
        };

        const rbQuery = new DistRectangle3OrientedBox3();
        const rectangle = new Rectangle(3);

        // The sentinel is negative and a squared distance is nonnegative, so
        // the first candidate always replaces it.
        const invalid = -1;
        result.distance = invalid;
        result.sqrDistance = invalid;

        // Compare faces of box0 to box1. The closest points of the two solids
        // (when the solids are separated) must occur on their boundaries, and
        // each boundary is the union of the six faces.
        for (let i0 = 2, i1 = 0, i2 = 1; i1 < 3; i2 = i0, i0 = i1++) {
            rectangle.axis[0] = box0.axis[i0].clone();
            rectangle.axis[1] = box0.axis[i1].clone();
            rectangle.extent.values[0] = box0.extent.values[i0];
            rectangle.extent.values[1] = box0.extent.values[i1];

            const scaledAxis = mul(box0.extent.values[i2], box0.axis[i2]);
            for (const sign of [+1, -1]) {
                rectangle.center = sign > 0
                    ? add(box0.center, scaledAxis)
                    : sub(box0.center, scaledAxis);
                const rbOutput = rbQuery.compute(rectangle, box1);
                if (result.sqrDistance === invalid ||
                    rbOutput.sqrDistance < result.sqrDistance) {
                    result.distance = rbOutput.distance;
                    result.sqrDistance = rbOutput.sqrDistance;
                    result.closest = [rbOutput.closest[0],
                        rbOutput.closest[1]];
                }
            }
        }

        // Compare faces of box1 to box0. The rectangle-box query reports the
        // rectangle point in closest[0] and the box point in closest[1], so
        // the pair must be swapped for the result.
        for (let i0 = 2, i1 = 0, i2 = 1; i1 < 3; i2 = i0, i0 = i1++) {
            rectangle.axis[0] = box1.axis[i0].clone();
            rectangle.axis[1] = box1.axis[i1].clone();
            rectangle.extent.values[0] = box1.extent.values[i0];
            rectangle.extent.values[1] = box1.extent.values[i1];

            const scaledAxis = mul(box1.extent.values[i2], box1.axis[i2]);
            for (const sign of [+1, -1]) {
                rectangle.center = sign > 0
                    ? add(box1.center, scaledAxis)
                    : sub(box1.center, scaledAxis);
                const rbOutput = rbQuery.compute(rectangle, box0);
                if (result.sqrDistance === invalid ||
                    rbOutput.sqrDistance < result.sqrDistance) {
                    result.distance = rbOutput.distance;
                    result.sqrDistance = rbOutput.sqrDistance;
                    result.closest = [rbOutput.closest[1],
                        rbOutput.closest[0]];
                }
            }
        }

        return result;
    }
}
