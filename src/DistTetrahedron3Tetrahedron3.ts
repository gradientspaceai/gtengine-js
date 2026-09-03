// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistTetrahedron3Tetrahedron3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between two solid tetrahedra in 3D.
//
// Each tetrahedron has vertices <V[0],V[1],V[2],V[3]>. A tetrahedron point
// is X = sum_{i=0}^3 b[i] * V[i], where 0 <= b[i] <= 1 for all i and
// sum_{i=0}^3 b[i] = 1.
//
// The closest point on tetra0 is stored in closest[0] with barycentric
// coordinates relative to its vertices. The closest point on tetra1 is
// stored in closest[1] with barycentric coordinates relative to its
// vertices. When there are infinitely many choices for the pair of closest
// points, only one pair is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Tetrahedron3<T>, Tetrahedron3<T>>'
// becomes the class DistTetrahedron3Tetrahedron3 with the result type
// DistTetrahedron3Tetrahedron3Result. Upstream discards the 'valid' flag of
// ComputeBarycentrics with a (void) cast; the port does the same, so the
// barycentric coordinates are all zero for a degenerate (zero-volume)
// tetrahedron.

import { inContainerTetrahedron3 } from './ContTetrahedron3.js';
import type { DCPQuery } from './DCPQuery.js';
import { DistTriangle3Triangle3 } from './DistTriangle3Triangle3.js';
import { Triangle } from './Triangle.js';
import { Tetrahedron3 } from './Tetrahedron3.js';
import { Vector } from './Vector.js';
import { computeBarycentrics3 } from './Vector3.js';

export interface DistTetrahedron3Tetrahedron3Result {
    distance: number;
    sqrDistance: number;

    // The barycentric coordinates of closest[0] relative to tetra0.
    barycentric0: [number, number, number, number];

    // The barycentric coordinates of closest[1] relative to tetra1.
    barycentric1: [number, number, number, number];

    // closest[0] is on tetra0, closest[1] is on tetra1.
    closest: [Vector, Vector];
}

export class DistTetrahedron3Tetrahedron3
    implements DCPQuery<Tetrahedron3, Tetrahedron3,
    DistTetrahedron3Tetrahedron3Result> {
    compute(tetra0: Tetrahedron3, tetra1: Tetrahedron3):
        DistTetrahedron3Tetrahedron3Result {
        const result: DistTetrahedron3Tetrahedron3Result = {
            distance: 0,
            sqrDistance: 0,
            barycentric0: [0, 0, 0, 0],
            barycentric1: [0, 0, 0, 0],
            closest: [new Vector(3), new Vector(3)]
        };

        const ttQuery = new DistTriangle3Triangle3();

        // The sentinel is negative and a squared distance is nonnegative, so
        // the first candidate always replaces it.
        const invalid = -1;
        result.distance = invalid;
        result.sqrDistance = invalid;

        // Compute the distances between pairs of faces, each pair having a
        // face from tetra0 and a face from tetra1.
        let foundZeroDistance = false;
        for (let face0 = 0; face0 < 4 && !foundZeroDistance; ++face0) {
            const triangle0 = new Triangle(3);
            const indices0 = Tetrahedron3.getFaceIndices(face0);
            for (let j = 0; j < 3; ++j) {
                triangle0.v[j] = tetra0.v[indices0[j]].clone();
            }

            for (let face1 = 0; face1 < 4; ++face1) {
                const triangle1 = new Triangle(3);
                const indices1 = Tetrahedron3.getFaceIndices(face1);
                for (let j = 0; j < 3; ++j) {
                    triangle1.v[j] = tetra1.v[indices1[j]].clone();
                }

                const ttResult = ttQuery.compute(triangle0, triangle1);
                if (ttResult.sqrDistance === 0) {
                    result.distance = 0;
                    result.sqrDistance = 0;
                    result.closest[0] = ttResult.closest[0];
                    result.closest[1] = ttResult.closest[1];
                    foundZeroDistance = true;
                    break;
                }

                if (result.sqrDistance === invalid ||
                    ttResult.sqrDistance < result.sqrDistance) {
                    result.distance = ttResult.distance;
                    result.sqrDistance = ttResult.sqrDistance;
                    result.closest[0] = ttResult.closest[0];
                    result.closest[1] = ttResult.closest[1];
                }
            }
        }

        if (!foundZeroDistance) {
            // The tetrahedra are either nested or separated. Test for
            // containment of the centroids to decide which case.
            const centroid0 = tetra0.computeCentroid();
            if (inContainerTetrahedron3(centroid0, tetra1)) {
                // Tetra0 is nested inside tetra1. Choose the centroid of
                // tetra0 as the closest point for both tetrahedra.
                result.distance = 0;
                result.sqrDistance = 0;
                result.closest[0] = centroid0;
                result.closest[1] = centroid0.clone();
            }

            const centroid1 = tetra1.computeCentroid();
            if (inContainerTetrahedron3(centroid1, tetra0)) {
                // Tetra1 is nested inside tetra0. Choose the centroid of
                // tetra1 as the closest point for both tetrahedra.
                result.distance = 0;
                result.sqrDistance = 0;
                result.closest[0] = centroid1;
                result.closest[1] = centroid1.clone();
            }

            // With exact arithmetic, at this point the tetrahedra are
            // separated. The result object already contains the distance
            // information. However, with floating-point arithmetic, it is
            // possible that a tetrahedron with volume nearly zero is close
            // enough to the other tetrahedron yet separated, but rounding
            // errors make it appear that the nearly-zero-volume tetrahedron
            // has centroid inside the other tetrahedron. This situation is
            // trapped by the previous two if-blocks.
        }

        // Compute the barycentric coordinates of the closest points.
        result.barycentric0 = computeBarycentrics3(result.closest[0],
            tetra0.v[0], tetra0.v[1], tetra0.v[2], tetra0.v[3]).bary;

        result.barycentric1 = computeBarycentrics3(result.closest[1],
            tetra1.v[0], tetra1.v[1], tetra1.v[2], tetra1.v[3]).bary;

        return result;
    }
}
