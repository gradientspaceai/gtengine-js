// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRectangle3CanonicalBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a rectangle and a solid canonical box in 3D.
//
// The rectangle has center C, unit-length axis directions W[0] and W[1], and
// extents e[0] and e[1]. A rectangle point is X = C + sum_{i=0}^1 s[i] * W[i]
// where |s[i]| <= e[i] for all i.
//
// The canonical box has center at the origin and is aligned with the
// coordinate axes. The extents are E = (e[0],e[1],e[2]). A box point is
// Y = (y[0],y[1],y[2]) with |y[i]| <= e[i] for all i.
//
// The closest point on the rectangle is stored in closest[0] with
// W-coordinates (s[0],s[1]). The closest point on the box is stored in
// closest[1]. When there are infinitely many choices for the pair of closest
// points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Rectangle3<T>, CanonicalBox3<T>>'
// becomes the class DistRectangle3CanonicalBox3 with the result type
// DistRectangle3CanonicalBox3Result.

import type { CanonicalBox3 } from './CanonicalBox.js';
import type { DCPQuery } from './DCPQuery.js';
import { DistPlane3CanonicalBox3 } from './DistPlane3CanonicalBox3.js';
import { DistSegment3CanonicalBox3 } from './DistSegment3CanonicalBox3.js';
import { Hyperplane } from './Hyperplane.js';
import type { Rectangle3 } from './Rectangle.js';
import { Segment } from './Segment.js';
import { Vector, dot, sub } from './Vector.js';
import { cross } from './Vector3.js';

export interface DistRectangle3CanonicalBox3Result {
    distance: number;
    sqrDistance: number;

    // The W-coordinates (s[0],s[1]) of closest[0].
    cartesian: [number, number];

    // closest[0] is on the rectangle, closest[1] is on the box.
    closest: [Vector, Vector];
}

// The rectangle-edge indices, listed as {horizontal, horizontal, vertical,
// vertical}. Rectangle.getVertices() returns the vertices in bit-pattern
// order, so edges {0,1} and {2,3} have y = -e[1] and y = +e[1] and edges
// {0,2} and {1,3} have x = -e[0] and x = +e[0].
const edges: readonly (readonly [number, number])[] =
    [[0, 1], [2, 3], [0, 2], [1, 3]];
const sign: readonly number[] = [-1, +1, -1, +1];
const j0: readonly number[] = [0, 0, 1, 1];
const j1: readonly number[] = [1, 1, 0, 0];

export class DistRectangle3CanonicalBox3
    implements DCPQuery<Rectangle3, CanonicalBox3,
    DistRectangle3CanonicalBox3Result> {
    compute(rectangle: Rectangle3, box: CanonicalBox3):
        DistRectangle3CanonicalBox3Result {
        const result: DistRectangle3CanonicalBox3Result = {
            distance: 0,
            sqrDistance: 0,
            cartesian: [0, 0],
            closest: [new Vector(3), new Vector(3)]
        };

        const pbQuery = new DistPlane3CanonicalBox3();
        const normal = cross(rectangle.axis[0], rectangle.axis[1]);
        const plane = Hyperplane.fromNormalOrigin(normal, rectangle.center);
        const pbOutput = pbQuery.compute(plane, box);
        const delta = sub(pbOutput.closest[0], rectangle.center);
        result.cartesian[0] = dot(rectangle.axis[0], delta);
        result.cartesian[1] = dot(rectangle.axis[1], delta);

        if (Math.abs(result.cartesian[0]) <= rectangle.extent.values[0] &&
            Math.abs(result.cartesian[1]) <= rectangle.extent.values[1]) {
            result.distance = pbOutput.distance;
            result.sqrDistance = pbOutput.sqrDistance;
            result.closest = [pbOutput.closest[0], pbOutput.closest[1]];
        }
        else {
            // The closest plane point is outside the rectangle, although it
            // is possible there are points inside the rectangle that also are
            // closest points to the box. Regardless, locate a point on an
            // edge of the rectangle that is closest to the box.
            const sbQuery = new DistSegment3CanonicalBox3();

            const invalid = -1;
            result.distance = invalid;
            result.sqrDistance = invalid;

            const vertices = rectangle.getVertices();

            for (let i = 0; i < 4; ++i) {
                const edge = edges[i];
                const segment = Segment.fromEndpoints(vertices[edge[0]],
                    vertices[edge[1]]);

                const sbOutput = sbQuery.compute(segment, box);
                if (result.sqrDistance === invalid ||
                    sbOutput.sqrDistance < result.sqrDistance) {
                    result.distance = sbOutput.distance;
                    result.sqrDistance = sbOutput.sqrDistance;
                    result.closest = [sbOutput.closest[0],
                        sbOutput.closest[1]];

                    const scale = 2 * sbOutput.parameter - 1;
                    result.cartesian[j0[i]] =
                        scale * rectangle.extent.values[j0[i]];
                    result.cartesian[j1[i]] =
                        sign[i] * rectangle.extent.values[j1[i]];
                }
            }
        }

        return result;
    }
}
