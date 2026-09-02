// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistTriangle3OrientedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a solid triangle and a solid oriented box
// in 3D.
//
// The triangle has vertices <V[0],V[1],V[2]>. A triangle point is
// X = sum_{i=0}^2 b[i] * V[i], where 0 <= b[i] <= 1 for all i and
// sum_{i=0}^2 b[i] = 1.
//
// The oriented box has center C, unit-length axis directions U[i] and extents
// e[i] for all i. A box point is X = C + sum_i y[i] * U[i], where
// |y[i]| <= e[i] for all i.
//
// The closest point on the triangle is stored in closest[0] with barycentric
// coordinates (b[0],b[1],b[2]). The closest point on the box is stored in
// closest[1]. When there are infinitely many choices for the pair of closest
// points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Triangle3<T>, OrientedBox3<T>>'
// becomes the class DistTriangle3OrientedBox3. As upstream does, the result
// type is the triangle-canonical-box result type, re-exported here as the
// alias DistTriangle3OrientedBox3Result.

import { CanonicalBox } from './CanonicalBox';
import type { DCPQuery } from './DCPQuery';
import { DistTriangle3CanonicalBox3 } from './DistTriangle3CanonicalBox3';
import type { DistTriangle3CanonicalBox3Result } from './DistTriangle3CanonicalBox3';
import type { OrientedBox3 } from './OrientedBox';
import { Triangle } from './Triangle';
import type { Triangle3 } from './Triangle';
import { Vector, add, dot, mul, sub } from './Vector';

// Upstream reuses the triangle-canonical-box result type ('using Result =
// typename TBQuery::Result').
export type DistTriangle3OrientedBox3Result =
    DistTriangle3CanonicalBox3Result;

export class DistTriangle3OrientedBox3
    implements DCPQuery<Triangle3, OrientedBox3,
    DistTriangle3OrientedBox3Result> {
    compute(triangle: Triangle3, box: OrientedBox3):
        DistTriangle3OrientedBox3Result {
        // Rotate and translate the triangle and box so that the box is
        // aligned and has center at the origin.
        const cbox = CanonicalBox.fromExtent(box.extent);
        const xfrmTriangle = new Triangle(3);
        for (let j = 0; j < 3; ++j) {
            const delta = sub(triangle.v[j], box.center);
            for (let i = 0; i < 3; ++i) {
                xfrmTriangle.v[j].values[i] = dot(box.axis[i], delta);
            }
        }

        // The query computes 'result' relative to the box with center at the
        // origin.
        const tbQuery = new DistTriangle3CanonicalBox3();
        const result = tbQuery.compute(xfrmTriangle, cbox);

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
