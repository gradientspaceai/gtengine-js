// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrTetrahedron3Tetrahedron3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the tetrahedron to be a solid.
//
// The test-intersection query uses the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The set of potential separating directions includes the 4 face normals of
// tetra0, the 4 face normals of tetra1, and 36 directions, each of which is
// the cross product of an edge of tetra0 and an edge of tetra1.
//
// The separating axes involving cross products of edges has numerical
// robustness problems when the two edges are nearly parallel. The cross
// product of the edges is nearly the zero vector, so normalization of the
// cross product may produce unit-length directions that are not close to the
// true direction. Such a pair of edges occurs when an object0 face normal N0
// and an object1 face normal N1 are nearly parallel. In this case, you may
// skip the edge-edge directions, which is equivalent to projecting the
// objects onto the plane with normal N0 and applying a 2D separating axis
// test. The ability to do so involves choosing a small nonnegative epsilon.
// It is used to determine whether two edge directions, one from each object,
// are nearly parallel: |Dot(U0,U1)| >= 1 - epsilon for unit-length U0 and U1,
// where 0 <= epsilon <= 1. The epsilon input to test(...) is clamped to
// [0,1]. An epsilon of 0 skips only exactly-parallel edge pairs, for which
// the cross product is the zero vector and carries no separating direction.
//
// The pair of integers 'separating', say, (i0,i1), identifies the axes that
// reported separation; there may be more than one but only one is reported.
// If the separating axis is a face normal N[i0] of object0, then
// (i0,invalid) is returned. If the axis is a face normal N[i1] of object1,
// then (invalid,i1) is returned. If the separating axis is a cross product
// of the edges with indices i0 and i1, then (i0,i1) is returned. If
// 'intersect' is true, the separating[] values are invalid because there is
// no separation.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The upstream
// TIQuery specialization becomes IntrTetrahedron3Tetrahedron3TI; the private
// static WhichSide helper becomes a module-private function. The SIZE_MAX
// invalid index becomes Number.MAX_SAFE_INTEGER (the BVTree.ts precedent),
// exported as intrTetrahedron3Tetrahedron3InvalidIndex.
//
// The port fixes two upstream defects in the edge-edge phase (both flagged
// in the code below): the parallelism test compares unnormalized edge dot
// products against a cosine cutoff, and the separation test uses a plane
// through an edge endpoint rather than testing the projection intervals for
// overlap. Both cause separated tetrahedra to be reported as intersecting;
// the second also makes the upstream query asymmetric in its arguments.

import { Tetrahedron3 } from './Tetrahedron3.js';
import type { TIQuery } from './TIQuery.js';
import { Vector, dot, sub } from './Vector.js';
import { cross } from './Vector3.js';

// The port of std::numeric_limits<size_t>::max() used to mark a 'separating'
// entry that does not participate in the reported separating axis.
export const intrTetrahedron3Tetrahedron3InvalidIndex = Number.MAX_SAFE_INTEGER;

// The result of IntrTetrahedron3Tetrahedron3TI.test.
export interface IntrTetrahedron3Tetrahedron3TIResult {
    intersect: boolean;

    // The indices of the axes that reported separation. See the file comments
    // for the interpretation.
    separating: [number, number];
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrTetrahedron3Tetrahedron3TIResult():
    IntrTetrahedron3Tetrahedron3TIResult {
    return {
        intersect: false,
        separating: [intrTetrahedron3Tetrahedron3InvalidIndex,
            intrTetrahedron3Tetrahedron3InvalidIndex]
    };
}

// The vertices of 'tetra' are projected to the form P + t * N. The return
// value is +1 if all t > 0, -1 if all t < 0, but 0 otherwise, in which case
// tetra has points on both sides of the plane Dot(N,X-P) = 0.
function whichSide(tetra: Tetrahedron3, P: Vector, N: Vector): number {
    let positive = 0, negative = 0;
    for (let i = 0; i < 4; ++i) {
        // Project a vertex onto the normal line.
        const t = dot(N, sub(tetra.v[i], P));
        if (t > 0) {
            ++positive;
        } else if (t < 0) {
            ++negative;
        }

        if (positive > 0 && negative > 0) {
            // Tetra has vertices on both sides of the line, so the line is
            // not a separating axis.
            return 0;
        }
    }

    // Either positive > 0 or negative > 0 but not both are positive.
    return positive > 0 ? +1 : -1;
}

// The interval of projection of the tetrahedron vertices onto the direction
// N; the direction does not have to be unit length.
function getProjection(tetra: Tetrahedron3, N: Vector):
    { min: number, max: number } {
    let min = dot(N, tetra.v[0]);
    let max = min;
    for (let i = 1; i < 4; ++i) {
        const d = dot(N, tetra.v[i]);
        if (d < min) {
            min = d;
        } else if (d > max) {
            max = d;
        }
    }
    return { min, max };
}

// Test-intersection query for two solid tetrahedra.
export class IntrTetrahedron3Tetrahedron3TI implements
    TIQuery<Tetrahedron3, Tetrahedron3, IntrTetrahedron3Tetrahedron3TIResult> {

    test(tetra0: Tetrahedron3, tetra1: Tetrahedron3, epsilon: number = 0):
        IntrTetrahedron3Tetrahedron3TIResult {
        const result = defaultIntrTetrahedron3Tetrahedron3TIResult();
        const invalid = intrTetrahedron3Tetrahedron3InvalidIndex;

        // Test the face normals of tetra0 for separation. Because of the
        // counterclockwise ordering of the face vertices relative to an
        // observer outside the tetrahedron, the projection interval for the
        // face is [T,0], where T < 0. Determine whether tetra1 is on the
        // positive side of the face-normal line.
        for (let i = 0; i < 4; ++i) {
            const faceIndices = Tetrahedron3.getFaceIndices(i);
            const P = tetra0.v[faceIndices[0]];
            const N = tetra0.computeFaceNormal(i);
            if (whichSide(tetra1, P, N) > 0) {
                // Tetra1 is entirely on the positive side of the normal line
                // P + t * N.
                result.intersect = false;
                result.separating[0] = i;
                result.separating[1] = invalid;
                return result;
            }
        }

        // Test the face normals of tetra1 for separation.
        for (let i = 0; i < 4; ++i) {
            const faceIndices = Tetrahedron3.getFaceIndices(i);
            const P = tetra1.v[faceIndices[0]];
            const N = tetra1.computeFaceNormal(i);
            if (whichSide(tetra0, P, N) > 0) {
                // Tetra0 is entirely on the positive side of the normal line
                // P + t * N.
                result.intersect = false;
                result.separating[0] = invalid;
                result.separating[1] = i;
                return result;
            }
        }

        // Test cross products of pairs of edge directions, one edge from each
        // tetrahedron.
        const cutoff = Math.min(Math.max(1 - epsilon, 0), 1);
        for (let i0 = 0; i0 < 6; ++i0) {
            const edge0Indices = Tetrahedron3.getEdgeIndices(i0);
            const P0 = tetra0.v[edge0Indices[0]];
            const E0 = sub(tetra0.v[edge0Indices[1]], P0);
            const lengthE0 = Math.sqrt(dot(E0, E0));
            for (let i1 = 0; i1 < 6; ++i1) {
                const edge1Indices = Tetrahedron3.getEdgeIndices(i1);
                const E1 = sub(tetra1.v[edge1Indices[1]],
                    tetra1.v[edge1Indices[0]]);
                const lengthE1 = Math.sqrt(dot(E1, E1));

                // Upstream bug (fixed here): upstream tests
                // |Dot(E0,E1)| < cutoff with UNNORMALIZED edge vectors, but
                // the cutoff is 1 - epsilon in [0,1], which is a threshold
                // for the cosine of the angle between the directions. For
                // tetrahedra whose edges are longer than one unit, almost
                // every |Dot(E0,E1)| exceeds the cutoff, so the edge-edge
                // axes are skipped and separated tetrahedra are reported as
                // intersecting. The port divides by the edge lengths so that
                // the test really is on the cosine, as the file comments
                // describe. Zero-length edges give no direction, so they are
                // skipped.
                const denom = lengthE0 * lengthE1;
                if (denom === 0) {
                    continue;
                }
                const cosAngle = Math.abs(dot(E0, E1)) / denom;

                if (cosAngle < cutoff) {
                    const N = cross(E0, E1);
                    if (dot(N, N) === 0) {
                        // The edges are parallel, so the cross product is the
                        // zero vector and carries no separating direction.
                        // The cosine test above can round to slightly less
                        // than 1 for parallel edges, so guard explicitly; a
                        // zero direction would otherwise report the two
                        // degenerate projection intervals as "touching" and
                        // hence separated.
                        continue;
                    }

                    // Upstream bug (fixed here): upstream tests the sides of
                    // the plane Dot(N,X-P0) = 0, where P0 is an endpoint of
                    // the tetra0 edge, and reports separation only when
                    // WhichSide(tetra0,P0,N) and WhichSide(tetra1,P0,N) are
                    // nonzero with opposite signs. That plane separates only
                    // when the tetra0 edge is extreme in the direction that
                    // faces tetra1. A separating direction N can have the
                    // tetra0 edge extreme on the far side instead, in which
                    // case upstream misses the separation and reports a
                    // false intersection; the outcome even depends on the
                    // argument order, so the upstream query is not symmetric.
                    // The method of separating axes asks whether the
                    // projection intervals of the two objects onto N are
                    // disjoint, which is what the port tests. Touching
                    // intervals count as separated, matching the strict-sign
                    // convention of the face-normal tests above.
                    const proj0 = getProjection(tetra0, N);
                    const proj1 = getProjection(tetra1, N);
                    if (proj0.max <= proj1.min || proj1.max <= proj0.min) {
                        // The projections of tetra0 and tetra1 onto the line
                        // P + t * N are disjoint intervals.
                        result.intersect = false;
                        result.separating[0] = i0;
                        result.separating[1] = i1;
                        return result;
                    }
                }
            }
        }

        result.intersect = true;
        result.separating[0] = invalid;
        result.separating[1] = invalid;
        return result;
    }
}
