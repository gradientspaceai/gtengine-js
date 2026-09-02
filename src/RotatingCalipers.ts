// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RotatingCalipers.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The rotating calipers algorithm finds all antipodal vertex-edge pairs for a
// convex polygon. The algorithm is O(n) in time for n polygon edges. The
// brute-force method that finds extreme points for a perpendicular direction
// for each edge and searching all polygon vertices is O(n^2). The search for
// extreme points can use a form of bisection, which reduces the algorithm to
// O(n log n). A description can be found at
// http://www-cgrl.cs.mcgill.ca/~godfried/research/calipers.html
// https://web.archive.org/web/20150330010154/http://cgm.cs.mcgill.ca/~orm/rotcal.html
//
// Port notes:
// * Upstream selects Rational = BSNumber<UIntegerFP32<NumWords>> where
//   NumWords is a compile-time worst-case bound (54 for float, 394 for
//   double). The port's BSNumber is bigint-backed and grows as needed, so
//   the word-count bound is unnecessary and is dropped. The exact-arithmetic
//   behavior is identical (no divisions are performed). This follows the
//   ConvexHull2 precedent.
// * The static member function ComputeAntipodes(vertices, antipodes) becomes
//   the static computeAntipodes(vertices) that returns the antipode array.
// * The private Rational vertex array is a plain array of [BSNumber,
//   BSNumber] tuples (the port of Vector2<Rational>). Upstream converts only
//   the retained polygon vertices; the port converts every vertex to keep
//   the array type non-nullable, which is behavior-preserving because the
//   entries for removed vertices are never read.

import { BSNumber } from './BSNumber';
import { logAssert } from './Logger';
import { Vector } from './Vector';

// The Antipode members are lookups into the input vertices[] to
// computeAntipodes(...). 'vertex' is the antipodal vertex and 'edge' holds
// the two endpoints of the antipodal edge.
export interface RotatingCalipersAntipode {
    vertex: number;
    edge: [number, number];
}

// The port of Vector2<Rational>.
type RationalPoint2 = [BSNumber, BSNumber];

function rSub(v0: RationalPoint2, v1: RationalPoint2): RationalPoint2 {
    return [v0[0].sub(v1[0]), v0[1].sub(v1[1])];
}

function rNegate(v: RationalPoint2): RationalPoint2 {
    return [v[0].negated(), v[1].negated()];
}

function rDot(v0: RationalPoint2, v1: RationalPoint2): BSNumber {
    return v0[0].mul(v1[0]).add(v0[1].mul(v1[1]));
}

function rDotPerp(v0: RationalPoint2, v1: RationalPoint2): BSNumber {
    return v0[0].mul(v1[1]).sub(v0[1].mul(v1[0]));
}

export class RotatingCalipers {
    // Compute all antipodal vertex-edge pairs of the convex polygon whose
    // vertices are given in counterclockwise order. Duplicate points and
    // collinear points are removed first; the polygon must have at least
    // three noncollinear vertices.
    static computeAntipodes(vertices: readonly Vector[]):
        RotatingCalipersAntipode[] {
        // Upstream reads vertices.back() and vertices[i0 + 1] without a
        // guard, which is undefined behavior for fewer than 2 vertices. The
        // final LogAssert on the retained-index count fires only after the
        // out-of-range reads, so the port checks the input size up front.
        logAssert(vertices.length >= 3,
            'The convex polygon must have at least 3 noncollinear vertices.');
        for (const vertex of vertices) {
            logAssert(vertex.size === 2,
                'RotatingCalipers: vertices must be 2D.');
        }

        // Internally, the antipode members are lookups into indices[]. The
        // members are re-mapped to lookups into vertices[] after all
        // antipodes are created.
        const { rVertices, indices } = RotatingCalipers.createPolygon(vertices);
        logAssert(indices.length >= 3,
            'The convex polygon must have at least 3 noncollinear vertices.');

        const antipode = RotatingCalipers.computeInitialAntipode(
            rVertices, indices);
        const antipodes: RotatingCalipersAntipode[] = [
            { vertex: antipode.vertex, edge: [antipode.edge[0], antipode.edge[1]] }
        ];

        for (let i = 1; i < indices.length; ++i) {
            RotatingCalipers.computeNextAntipode(rVertices, indices, antipode);
            antipodes.push({
                vertex: antipode.vertex,
                edge: [antipode.edge[0], antipode.edge[1]]
            });
        }

        // Re-map the antipode members to be lookups into vertices[].
        for (const element of antipodes) {
            element.vertex = indices[element.vertex];
            element.edge[0] = indices[element.edge[0]];
            element.edge[1] = indices[element.edge[1]];
        }

        return antipodes;
    }

    // The rotating calipers algorithm requires the convex polygon to have no
    // duplicate points and no collinear points. Such points must be removed
    // first. To ensure correctness, rational arithmetic is used. This
    // requires converting the floating-point vertices to rational vertices.
    //
    // A vertex is retained when the edge arriving at it and the edge leaving
    // it are not parallel, so that duplicate points (zero-length edges) and
    // collinear points are dropped.
    //
    // Upstream bug (RotatingCalipers.h, CreatePolygon): the "arriving" edge
    // is the immediately preceding edge of the input array, which may be the
    // zero-length edge produced by a duplicate point. DotPerp is then zero
    // and a genuine corner of the polygon is discarded along with the
    // duplicate. For example, the square (0,0),(0,0),(1,0),(1,1),(1,1),(0,1)
    // retains only 2 of its 4 corners and the LogAssert on the retained
    // count then fires. The port fixes this by comparing against the most
    // recent *nonzero* edge, which is identical to upstream whenever the
    // input has no duplicate points.
    private static createPolygon(vertices: readonly Vector[]):
        { rVertices: RationalPoint2[], indices: number[] } {
        const numVertices = vertices.length;
        const rzero = BSNumber.fromNumber(0);

        const rVertices: RationalPoint2[] = vertices.map(v =>
            [BSNumber.fromNumber(v.get(0)),
                BSNumber.fromNumber(v.get(1))] as RationalPoint2);

        const edges: RationalPoint2[] = [];
        const isZeroEdge: boolean[] = [];
        for (let i = 0; i < numVertices; ++i) {
            const edge = rSub(rVertices[(i + 1) % numVertices], rVertices[i]);
            edges.push(edge);
            isZeroEdge.push(edge[0].equals(rzero) && edge[1].equals(rzero));
        }

        // The most recent nonzero edge before index 0, cyclically.
        let prevNonzero = -1;
        for (let i = numVertices - 1; i >= 0; --i) {
            if (!isZeroEdge[i]) {
                prevNonzero = i;
                break;
            }
        }

        const indices: number[] = [];
        if (prevNonzero < 0) {
            // All input points are identical; there are no edges at all.
            return { rVertices, indices };
        }

        for (let i0 = 0; i0 < numVertices; ++i0) {
            if (isZeroEdge[i0]) {
                continue;
            }
            const rDP = rDotPerp(edges[prevNonzero], edges[i0]);
            if (rDP.notEquals(rzero)) {
                indices.push(i0);
            }
            prevNonzero = i0;
        }

        return { rVertices, indices };
    }

    private static computeInitialAntipode(vertices: readonly RationalPoint2[],
        indices: readonly number[]): RotatingCalipersAntipode {
        const numIndices = indices.length;
        const antipode: RotatingCalipersAntipode = {
            vertex: 0,
            edge: [numIndices - 1, 0]
        };

        const origin = vertices[indices[antipode.edge[0]]];
        const U = rSub(vertices[indices[antipode.edge[1]]], origin);

        const zero = BSNumber.fromNumber(0);
        let extreme: RationalPoint2 = [zero, zero];
        antipode.vertex = 0;
        for (let i = 0; i < numIndices; ++i) {
            const diff = rSub(vertices[indices[i]], origin);
            const c: RationalPoint2 = [
                U[0].mul(diff[0]).add(U[1].mul(diff[1])),
                U[0].mul(diff[1]).sub(U[1].mul(diff[0]))
            ];

            if (c[1].greaterThan(extreme[1])
                || (c[1].equals(extreme[1]) && c[0].lessThan(extreme[0]))) {
                antipode.vertex = i;
                extreme = c;
            }
        }

        return antipode;
    }

    private static computeNextAntipode(vertices: readonly RationalPoint2[],
        indices: readonly number[],
        antipode: RotatingCalipersAntipode): void {
        // Given edges E0 and E1 we know that the angle between them is
        // determined by Dot(E0,E1)/(|E0|*|E1|) = cos(angle). The angle is in
        // (0,pi/2] when Dot(E0,E1) >= 0 or in (pi/2,pi) when Dot(E0,E1) < 0.
        // To allow for exact arithmetic, observe that
        //   sin^2(angle) = 1 - cos^2(angle)
        //                = 1 - Dot(E0,E1)^2/(|E0|^2*|E1|^2)
        // The comparator function for angles in (0,pi) compares the squared
        // sine values and the signs of the dot product of edges.
        const numIndices = indices.length;

        // Compute the edges associated with the current antipodal edge.
        const i0 = indices[antipode.edge[0]];
        const i1 = indices[antipode.edge[1]];
        let enext = antipode.edge[1] + 1;
        if (enext === numIndices) {
            enext = 0;
        }
        const i2 = indices[enext];

        // Compute the edges associated with the current antipodal vertex.
        const j0 = indices[antipode.vertex];
        let vnext = antipode.vertex + 1;
        if (vnext === numIndices) {
            vnext = 0;
        }
        const j1 = indices[vnext];

        const D0: [RationalPoint2, RationalPoint2] = [
            rSub(vertices[j1], vertices[j0]),
            rSub(vertices[i0], vertices[i1])
        ];

        const D1: [RationalPoint2, RationalPoint2] = [
            rNegate(D0[1]),
            rSub(vertices[i2], vertices[i1])
        ];

        if (RotatingCalipers.angleLessThan(D0, D1)) {
            // The angle at the antipodal vertex is minimum.
            const temp = antipode.vertex;
            antipode.vertex = antipode.edge[1];
            antipode.edge[1] = temp;
            antipode.edge[0] = antipode.edge[1];
            antipode.edge[1] = vnext;
        } else {
            // The angle at the antipodal edge is minimum. The antipodal
            // vertex does not change.
            antipode.edge[0] = antipode.edge[1];
            antipode.edge[1] = enext;
        }
    }

    // Test Angle(D0[0],D0[1]) < Angle(D1[0],D1[1]). It is known that
    // D1[0] = -D0[1].
    private static angleLessThan(D0: readonly [RationalPoint2, RationalPoint2],
        D1: readonly [RationalPoint2, RationalPoint2]): boolean {
        const zero = BSNumber.fromNumber(0);
        const dot0 = rDot(D0[0], D0[1]);
        const dot1 = rDot(D1[0], D1[1]);

        if (dot0.greaterThanOrEqual(zero)) {
            // angle0 in (0,pi/2]
            if (dot1.lessThan(zero)) {
                // angle1 in (pi/2,pi), so angle0 < angle1
                return true;
            }

            // angle1 in (0,pi/2], sin^2(angle) is an increasing function
            const sqrLen00 = rDot(D0[0], D0[0]);
            const sqrLen11 = rDot(D1[1], D1[1]);
            return dot0.mul(dot0).mul(sqrLen11)
                .greaterThan(dot1.mul(dot1).mul(sqrLen00));
        } else {
            // angle0 in (pi/2,pi)
            if (dot1.greaterThanOrEqual(zero)) {
                // angle1 in (0,pi/2], so angle1 < angle0
                return false;
            }

            // angle1 in (pi/2,pi), sin^2(angle) is a decreasing function
            const sqrLen00 = rDot(D0[0], D0[0]);
            const sqrLen11 = rDot(D1[1], D1[1]);
            return dot0.mul(dot0).mul(sqrLen11)
                .lessThan(dot1.mul(dot1).mul(sqrLen00));
        }
    }
}
