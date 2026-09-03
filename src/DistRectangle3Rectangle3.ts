// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRectangle3Rectangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between two rectangles in 3D.
//
// Each rectangle has center C, unit-length axis directions W[0] and W[1], and
// extents e[0] and e[1]. A rectangle point is X = C + sum_{i=0}^1 s[i] * W[i]
// where |s[i]| <= e[i] for all i.
//
// The closest point on rectangle0 is stored in closest[0] with W-coordinates
// (s[0],s[1]) corresponding to its W-axes. The closest point on rectangle1 is
// stored in closest[1] with W-coordinates (s[0],s[1]) corresponding to its
// W-axes. When there are infinitely many choices for the pair of closest
// points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Rectangle3<T>, Rectangle3<T>>' becomes
// the class DistRectangle3Rectangle3 with the result type
// DistRectangle3Rectangle3Result.

import type { DCPQuery } from './DCPQuery.js';
import { DistSegment3Rectangle3 } from './DistSegment3Rectangle3.js';
import type { Rectangle3 } from './Rectangle.js';
import { Segment } from './Segment.js';
import { Vector } from './Vector.js';

export interface DistRectangle3Rectangle3Result {
    distance: number;
    sqrDistance: number;

    // The W-coordinates of closest[0] relative to rectangle0.
    cartesian0: [number, number];

    // The W-coordinates of closest[1] relative to rectangle1.
    cartesian1: [number, number];

    // closest[0] is on rectangle0, closest[1] is on rectangle1.
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

export class DistRectangle3Rectangle3
    implements DCPQuery<Rectangle3, Rectangle3,
    DistRectangle3Rectangle3Result> {
    compute(rectangle0: Rectangle3, rectangle1: Rectangle3):
        DistRectangle3Rectangle3Result {
        const result: DistRectangle3Rectangle3Result = {
            distance: 0,
            sqrDistance: 0,
            cartesian0: [0, 0],
            cartesian1: [0, 0],
            closest: [new Vector(3), new Vector(3)]
        };

        const srQuery = new DistSegment3Rectangle3();

        const invalid = -1;
        result.distance = invalid;
        result.sqrDistance = invalid;

        // Compare edges of rectangle0 to the interior of rectangle1.
        const vertices0 = rectangle0.getVertices();
        for (let i = 0; i < 4; ++i) {
            const edge = edges[i];
            const segment = Segment.fromEndpoints(vertices0[edge[0]],
                vertices0[edge[1]]);

            const srResult = srQuery.compute(segment, rectangle1);
            if (result.sqrDistance === invalid ||
                srResult.sqrDistance < result.sqrDistance) {
                result.distance = srResult.distance;
                result.sqrDistance = srResult.sqrDistance;
                const scale = 2 * srResult.parameter - 1;
                result.cartesian0[j0[i]] =
                    scale * rectangle0.extent.values[j0[i]];
                result.cartesian0[j1[i]] =
                    sign[i] * rectangle0.extent.values[j1[i]];
                result.cartesian1 = [srResult.cartesian[0],
                    srResult.cartesian[1]];
                result.closest = [srResult.closest[0], srResult.closest[1]];
            }
        }

        // Compare edges of rectangle1 to the interior of rectangle0.
        const vertices1 = rectangle1.getVertices();
        for (let i = 0; i < 4; ++i) {
            const edge = edges[i];
            const segment = Segment.fromEndpoints(vertices1[edge[0]],
                vertices1[edge[1]]);

            const srResult = srQuery.compute(segment, rectangle0);
            if (result.sqrDistance === invalid ||
                srResult.sqrDistance < result.sqrDistance) {
                result.distance = srResult.distance;
                result.sqrDistance = srResult.sqrDistance;
                const scale = 2 * srResult.parameter - 1;
                result.cartesian0 = [srResult.cartesian[0],
                    srResult.cartesian[1]];
                result.cartesian1[j0[i]] =
                    scale * rectangle1.extent.values[j0[i]];
                result.cartesian1[j1[i]] =
                    sign[i] * rectangle1.extent.values[j1[i]];
                result.closest = [srResult.closest[1], srResult.closest[0]];
            }
        }

        return result;
    }
}
