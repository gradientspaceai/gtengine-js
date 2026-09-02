// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistTriangle3CanonicalBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a solid triangle and a solid canonical box
// in 3D.
//
// The triangle has vertices <V[0],V[1],V[2]>. A triangle point is
// X = sum_{i=0}^2 b[i] * V[i], where 0 <= b[i] <= 1 for all i and
// sum_{i=0}^2 b[i] = 1.
//
// The canonical box has center at the origin and is aligned with the
// coordinate axes. The extents are E = (e[0],e[1],e[2]). A box point is
// Y = (y[0],y[1],y[2]) with |y[i]| <= e[i] for all i.
//
// The closest point on the triangle is stored in closest[0] with barycentric
// coordinates (b[0],b[1],b[2]). The closest point on the box is stored in
// closest[1]. When there are infinitely many choices for the pair of closest
// points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Triangle3<T>, CanonicalBox3<T>>'
// becomes the class DistTriangle3CanonicalBox3 with the result type
// DistTriangle3CanonicalBox3Result.

import type { CanonicalBox3 } from './CanonicalBox';
import type { DCPQuery } from './DCPQuery';
import { DistPlane3CanonicalBox3 } from './DistPlane3CanonicalBox3';
import { DistSegment3CanonicalBox3 } from './DistSegment3CanonicalBox3';
import { Hyperplane } from './Hyperplane';
import { Segment } from './Segment';
import type { Triangle3 } from './Triangle';
import { Vector, dot, normalize, sub } from './Vector';
import { cross } from './Vector3';

export interface DistTriangle3CanonicalBox3Result {
    distance: number;
    sqrDistance: number;

    // The barycentric coordinates (b[0],b[1],b[2]) of closest[0] relative to
    // the triangle vertices.
    barycentric: [number, number, number];

    // closest[0] is on the triangle, closest[1] is on the box.
    closest: [Vector, Vector];
}

export class DistTriangle3CanonicalBox3
    implements DCPQuery<Triangle3, CanonicalBox3,
    DistTriangle3CanonicalBox3Result> {
    compute(triangle: Triangle3, box: CanonicalBox3):
        DistTriangle3CanonicalBox3Result {
        const result: DistTriangle3CanonicalBox3Result = {
            distance: 0,
            sqrDistance: 0,
            barycentric: [0, 0, 0],
            closest: [new Vector(3), new Vector(3)]
        };

        const E10 = sub(triangle.v[1], triangle.v[0]);
        const E20 = sub(triangle.v[2], triangle.v[0]);
        const K = cross(E10, E20);
        const sqrLength = dot(K, K);
        const N = K.clone();
        normalize(N);

        const pbQuery = new DistPlane3CanonicalBox3();
        const plane = Hyperplane.fromNormalOrigin(N, triangle.v[0]);
        const pbOutput = pbQuery.compute(plane, box);

        // closest[0] = b[0] * V[0] + b[1] * V[1] + b[2] * V[2]
        // = V[0] + b[1] * (V[1] - V[0]) + b[2] * (V[2] - V[0]);
        // delta = closest[0] - V[0] = b[1] * E10 + b[2] * E20
        const delta = sub(pbOutput.closest[0], triangle.v[0]);
        const KxDelta = cross(K, delta);
        result.barycentric[1] = dot(E20, KxDelta) / sqrLength;
        result.barycentric[2] = -dot(E10, KxDelta) / sqrLength;
        result.barycentric[0] = 1 - result.barycentric[1]
            - result.barycentric[2];

        if (0 <= result.barycentric[0] && result.barycentric[0] <= 1 &&
            0 <= result.barycentric[1] && result.barycentric[1] <= 1 &&
            0 <= result.barycentric[2] && result.barycentric[2] <= 1) {
            result.distance = pbOutput.distance;
            result.sqrDistance = pbOutput.sqrDistance;
            result.closest = [pbOutput.closest[0], pbOutput.closest[1]];
        }
        else {
            // The closest plane point is outside the triangle, although it is
            // possible there are points inside the triangle that also are
            // closest points to the box. Regardless, locate a point on an
            // edge of the triangle that is closest to the box.
            const sbQuery = new DistSegment3CanonicalBox3();

            const invalid = -1;
            result.distance = invalid;
            result.sqrDistance = invalid;

            // Compare edges of the triangle to the box.
            for (let i0 = 2, i1 = 0, i2 = 1; i1 < 3; i2 = i0, i0 = i1++) {
                const segment = Segment.fromEndpoints(triangle.v[i0],
                    triangle.v[i1]);

                const sbOutput = sbQuery.compute(segment, box);
                if (result.sqrDistance === invalid ||
                    sbOutput.sqrDistance < result.sqrDistance) {
                    result.distance = sbOutput.distance;
                    result.sqrDistance = sbOutput.sqrDistance;
                    result.barycentric[i0] = 1 - sbOutput.parameter;
                    result.barycentric[i1] = sbOutput.parameter;
                    result.barycentric[i2] = 0;
                    result.closest = [sbOutput.closest[0], sbOutput.closest[1]];
                }
            }
        }

        return result;
    }
}
