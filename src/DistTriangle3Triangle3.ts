// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistTriangle3Triangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between two solid triangles in 3D.
//
// Each triangle has vertices <V[0],V[1],V[2]>. A triangle point is
// X = sum_{i=0}^2 b[i] * V[i], where 0 <= b[i] <= 1 for all i and
// sum_{i=0}^2 b[i] = 1.
//
// The closest point on triangle0 is stored in closest[0] with barycentric
// coordinates relative to its vertices. The closest point on triangle1 is
// stored in closest[1] with barycentric coordinates relative to its
// vertices. When there are infinitely many choices for the pair of closest
// points, only one pair is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Triangle3<T>, Triangle3<T>>' becomes
// the class DistTriangle3Triangle3 with the result type
// DistTriangle3Triangle3Result.

import type { DCPQuery } from './DCPQuery';
import { DistSegment3Triangle3 } from './DistSegment3Triangle3';
import { Segment } from './Segment';
import type { Triangle3 } from './Triangle';
import { Vector } from './Vector';

export interface DistTriangle3Triangle3Result {
    distance: number;
    sqrDistance: number;

    // The barycentric coordinates of closest[0] relative to triangle0.
    barycentric0: [number, number, number];

    // The barycentric coordinates of closest[1] relative to triangle1.
    barycentric1: [number, number, number];

    // closest[0] is on triangle0, closest[1] is on triangle1.
    closest: [Vector, Vector];
}

export class DistTriangle3Triangle3
    implements DCPQuery<Triangle3, Triangle3, DistTriangle3Triangle3Result> {
    compute(triangle0: Triangle3, triangle1: Triangle3):
        DistTriangle3Triangle3Result {
        const result: DistTriangle3Triangle3Result = {
            distance: 0,
            sqrDistance: 0,
            barycentric0: [0, 0, 0],
            barycentric1: [0, 0, 0],
            closest: [new Vector(3), new Vector(3)]
        };

        const stQuery = new DistSegment3Triangle3();

        const invalid = -1;
        result.distance = invalid;
        result.sqrDistance = invalid;

        // Compare edges of triangle0 to the interior of triangle1.
        for (let i0 = 2, i1 = 0, i2 = 1; i1 < 3; i2 = i0, i0 = i1++) {
            const segment = Segment.fromEndpoints(triangle0.v[i0],
                triangle0.v[i1]);

            const stResult = stQuery.compute(segment, triangle1);
            if (result.sqrDistance === invalid ||
                stResult.sqrDistance < result.sqrDistance) {
                result.distance = stResult.distance;
                result.sqrDistance = stResult.sqrDistance;
                result.barycentric0[i0] = 1 - stResult.parameter;
                result.barycentric0[i1] = stResult.parameter;
                result.barycentric0[i2] = 0;
                result.barycentric1 = [stResult.barycentric[0],
                    stResult.barycentric[1], stResult.barycentric[2]];
                result.closest = [stResult.closest[0], stResult.closest[1]];
            }
        }

        // Compare edges of triangle1 to the interior of triangle0.
        for (let i0 = 2, i1 = 0, i2 = 1; i1 < 3; i2 = i0, i0 = i1++) {
            const segment = Segment.fromEndpoints(triangle1.v[i0],
                triangle1.v[i1]);

            const stResult = stQuery.compute(segment, triangle0);
            if (result.sqrDistance === invalid ||
                stResult.sqrDistance < result.sqrDistance) {
                result.distance = stResult.distance;
                result.sqrDistance = stResult.sqrDistance;
                result.barycentric0 = [stResult.barycentric[0],
                    stResult.barycentric[1], stResult.barycentric[2]];
                result.barycentric1[i0] = 1 - stResult.parameter;
                result.barycentric1[i1] = stResult.parameter;
                result.barycentric1[i2] = 0;
                result.closest = [stResult.closest[1], stResult.closest[0]];
            }
        }

        return result;
    }
}
