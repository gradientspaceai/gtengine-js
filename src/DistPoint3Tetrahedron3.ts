// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPoint3Tetrahedron3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a point and a solid tetrahedron in 3D.
//
// The tetrahedron is represented as an array of four vertices, V[i] for
// 0 <= i <= 3. The vertices are ordered so that the triangular faces are
// counterclockwise-ordered triangles when viewed by an observer outside the
// tetrahedron: face 0 = <V[0],V[2],V[1]>, face 1 = <V[0],V[1],V[3]>,
// face 2 = <V[0],V[3],V[2]> and face 3 = <V[1],V[2],V[3]>. The canonical
// tetrahedron has V[0] = (0,0,0), V[1] = (1,0,0), V[2] = (0,1,0) and
// V[3] = (0,0,1). A tetrahedron point is X = sum_{i=0}^3 b[i] * V[i],
// where 0 <= b[i] <= 1 for all i and sum_{i=0}^3 b[i] = 1.
//
// The input P is stored in closest[0]. The closest point on the tetrahedron
// is stored in closest[1] with barycentric coordinates
// (b[0],b[1],b[2],b[3]).
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector3<T>, Tetrahedron3<T>>' becomes
// the class DistPoint3Tetrahedron3 with the result type
// DistPoint3Tetrahedron3Result. Tetrahedron3.getPlanes() returns the four
// planes and Tetrahedron3.getFaceIndices(i) is a static method in the port.

import type { DCPQuery } from './DCPQuery.js';
import { DistPointTriangle } from './DistPointTriangle.js';
import { Tetrahedron3 } from './Tetrahedron3.js';
import { Triangle } from './Triangle.js';
import { Vector, dot } from './Vector.js';
import { computeBarycentrics3 } from './Vector3.js';

export interface DistPoint3Tetrahedron3Result {
    distance: number;
    sqrDistance: number;

    // The barycentric coordinates (b[0],b[1],b[2],b[3]) of closest[1]
    // relative to the tetrahedron vertices; they are nonnegative and sum
    // to 1.
    barycentric: [number, number, number, number];

    // closest[0] is the input point, closest[1] is the closest tetrahedron
    // point.
    closest: [Vector, Vector];
}

export class DistPoint3Tetrahedron3
    implements DCPQuery<Vector, Tetrahedron3, DistPoint3Tetrahedron3Result> {
    compute(point: Vector, tetrahedron: Tetrahedron3):
        DistPoint3Tetrahedron3Result {
        const result: DistPoint3Tetrahedron3Result = {
            distance: 0,
            sqrDistance: 0,
            barycentric: [0, 0, 0, 0],
            closest: [new Vector(3), new Vector(3)]
        };

        // Construct the planes for the faces of the tetrahedron. The normals
        // are outer pointing. We only need to know the sidedness of the query
        // point.
        const planes = tetrahedron.getPlanes();

        // Determine which faces are visible to the query point. Only these
        // need to be processed by point-to-triangle distance queries. To
        // allow for arbitrary precision arithmetic, the initial sqrDistance
        // is initialized to an invalid value rather than a floating-point
        // maximum value. Tracking the minimum requires a small amount of
        // extra logic.
        const invalid = -1;
        result.sqrDistance = invalid;
        for (let i = 0; i < 4; ++i) {
            if (dot(planes[i].normal, point) >= planes[i].constant) {
                const indices = Tetrahedron3.getFaceIndices(i);
                const triangle = Triangle.fromVertices(
                    tetrahedron.v[indices[0]],
                    tetrahedron.v[indices[1]],
                    tetrahedron.v[indices[2]]);

                const ptQuery = new DistPointTriangle();
                const ptResult = ptQuery.compute(point, triangle);
                if (result.sqrDistance === invalid ||
                    ptResult.sqrDistance < result.sqrDistance) {
                    result.sqrDistance = ptResult.sqrDistance;
                    result.closest = [ptResult.closest[0], ptResult.closest[1]];
                }
            }
        }

        if (result.sqrDistance === invalid) {
            // The query point is inside the solid tetrahedron. Report a zero
            // distance. The closest points are identical.
            result.sqrDistance = 0;
            result.closest[0] = point.clone();
            result.closest[1] = point.clone();
        }
        result.distance = Math.sqrt(result.sqrDistance);

        const bary = computeBarycentrics3(result.closest[1], tetrahedron.v[0],
            tetrahedron.v[1], tetrahedron.v[2], tetrahedron.v[3]);
        result.barycentric = bary.bary;

        return result;
    }
}
