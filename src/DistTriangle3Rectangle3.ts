// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistTriangle3Rectangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a solid triangle and a solid rectangle in 3D.
//
// The triangle has vertices <V[0],V[1],V[2]>. A triangle point is
// X = sum_{i=0}^2 b[i] * V[i], where 0 <= b[i] <= 1 for all i and
// sum_{i=0}^2 b[i] = 1.
//
// The rectangle has center C, unit-length axis directions W[0] and W[1], and
// extents e[0] and e[1]. A rectangle point is X = C + sum_{i=0}^1 s[i] * W[i]
// where |s[i]| <= e[i] for all i.
//
// The closest point on the triangle is stored in closest[0] with barycentric
// coordinates (b[0],b[1],b[2]). The closest point on the rectangle is stored
// in closest[1] with cartesian coordinates (s[0],s[1]). When there are
// infinitely many choices for the pair of closest points, only one of them is
// returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Triangle3<T>, Rectangle3<T>>' becomes
// the class DistTriangle3Rectangle3 with the result type
// DistTriangle3Rectangle3Result.

import type { DCPQuery } from './DCPQuery';
import { DistSegment3Rectangle3 } from './DistSegment3Rectangle3';
import { DistSegment3Triangle3 } from './DistSegment3Triangle3';
import type { Rectangle3 } from './Rectangle';
import { Segment } from './Segment';
import type { Triangle3 } from './Triangle';
import { Vector } from './Vector';

export interface DistTriangle3Rectangle3Result {
    distance: number;
    sqrDistance: number;

    // The barycentric coordinates of closest[0] relative to the triangle.
    barycentric: [number, number, number];

    // The W-coordinates (s[0],s[1]) of closest[1].
    cartesian: [number, number];

    // closest[0] is on the triangle, closest[1] is on the rectangle.
    closest: [Vector, Vector];
}

// The rectangle-edge index pairs, listed as {horizontal, horizontal,
// vertical, vertical}. Rectangle.getVertices() returns the vertices in
// bit-pattern order, so edges {0,1} and {2,3} have y = -e[1] and y = +e[1],
// and edges {0,2} and {1,3} have x = -e[0] and x = +e[0].
const edges: readonly (readonly number[])[] =
    [[0, 1], [2, 3], [0, 2], [1, 3]];
const sign: readonly number[] = [-1, +1, -1, +1];
const j0: readonly number[] = [0, 0, 1, 1];
const j1: readonly number[] = [1, 1, 0, 0];

export class DistTriangle3Rectangle3
    implements DCPQuery<Triangle3, Rectangle3, DistTriangle3Rectangle3Result> {
    compute(triangle: Triangle3, rectangle: Rectangle3):
        DistTriangle3Rectangle3Result {
        const result: DistTriangle3Rectangle3Result = {
            distance: 0,
            sqrDistance: 0,
            barycentric: [0, 0, 0],
            cartesian: [0, 0],
            closest: [new Vector(3), new Vector(3)]
        };

        const stQuery = new DistSegment3Triangle3();
        const srQuery = new DistSegment3Rectangle3();

        const invalid = -1;
        result.distance = invalid;
        result.sqrDistance = invalid;

        // Compare edges of the triangle to the interior of the rectangle.
        for (let i0 = 2, i1 = 0, i2 = 1; i1 < 3; i2 = i0, i0 = i1++) {
            const segment = Segment.fromEndpoints(triangle.v[i0],
                triangle.v[i1]);

            const srResult = srQuery.compute(segment, rectangle);
            if (result.sqrDistance === invalid ||
                srResult.sqrDistance < result.sqrDistance) {
                result.distance = srResult.distance;
                result.sqrDistance = srResult.sqrDistance;
                result.barycentric[i0] = 1 - srResult.parameter;
                result.barycentric[i1] = srResult.parameter;
                result.barycentric[i2] = 0;
                result.cartesian = [srResult.cartesian[0],
                    srResult.cartesian[1]];
                result.closest = [srResult.closest[0], srResult.closest[1]];
            }
        }

        // Compare edges of the rectangle to the interior of the triangle.
        const vertices = rectangle.getVertices();
        for (let i = 0; i < 4; ++i) {
            const edge = edges[i];
            const segment = Segment.fromEndpoints(vertices[edge[0]],
                vertices[edge[1]]);

            const stResult = stQuery.compute(segment, triangle);
            if (result.sqrDistance === invalid ||
                stResult.sqrDistance < result.sqrDistance) {
                result.distance = stResult.distance;
                result.sqrDistance = stResult.sqrDistance;
                result.barycentric = [stResult.barycentric[0],
                    stResult.barycentric[1], stResult.barycentric[2]];
                const scale = 2 * stResult.parameter - 1;
                result.cartesian[j0[i]] =
                    scale * rectangle.extent.values[j0[i]];
                result.cartesian[j1[i]] =
                    sign[i] * rectangle.extent.values[j1[i]];
                result.closest = [stResult.closest[1], stResult.closest[0]];
            }
        }

        return result;
    }
}
