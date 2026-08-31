// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLine2Triangle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a line and a solid triangle in 2D.
//
// The line is P + t * D, where D is not required to be unit length.
//
// The triangle has vertices <V[0],V[1],V[2]>. A triangle point is
// X = sum_{i=0}^2 b[i] * V[i], where 0 <= b[i] <= 1 for all i and
// sum_{i=0}^2 b[i] = 1.
//
// The closest point on the line is stored in closest[0] with parameter t. The
// closest point on the triangle is closest[1] with barycentric coordinates
// (b[0],b[1],b[2]). When there are infinitely many choices for the pair of
// closest points, only one of them is returned. Upstream TODO: compute the
// entire set of intersection when it is a line segment.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Line2<T>, Triangle2<T>>' becomes the
// class DistLine2Triangle2 with the result type DistLine2Triangle2Result. The
// private helpers LineContainsVertex, LineIntersectsTwoEdges and
// NoCommonPoints become module-private functions.

import type { DCPQuery } from './DCPQuery';
import type { Line2 } from './Line';
import type { Triangle2 } from './Triangle';
import { Vector, add, dot, mul, sub } from './Vector';
import { dotPerp, perp } from './Vector2';

export interface DistLine2Triangle2Result {
    distance: number;
    sqrDistance: number;

    // The line parameter t of the closest line point.
    parameter: number;

    // The barycentric coordinates of the closest triangle point.
    barycentric: [number, number, number];

    // closest[0] is on the line, closest[1] is on the triangle.
    closest: [Vector, Vector];
}

function defaultResult(): DistLine2Triangle2Result {
    return {
        distance: 0,
        sqrDistance: 0,
        parameter: 0,
        barycentric: [0, 0, 0],
        closest: [new Vector(2), new Vector(2)]
    };
}

// The line contains the triangle vertex V[i0], so the distance is 0.
function lineContainsVertex(P: Vector, D: Vector, V: readonly Vector[],
    i0: number, i1: number, i2: number,
    result: DistLine2Triangle2Result): void {
    result.distance = 0;
    result.sqrDistance = 0;
    result.parameter = dot(D, sub(V[i0], P)) / dot(D, D);
    result.barycentric[i0] = 1;
    result.barycentric[i1] = 0;
    result.barycentric[i2] = 0;
    result.closest[0] = V[i0].clone();
    result.closest[1] = V[i0].clone();
}

// At V[i0] and V[i1] the signs satisfy sign[i0] * sign[i1] < 0.
function lineIntersectsTwoEdges(P: Vector, D: Vector, V: readonly Vector[],
    i0: number, i1: number, i2: number,
    result: DistLine2Triangle2Result): void {
    const s = dotPerp(D, sub(P, V[i0])) / dotPerp(D, sub(V[i1], V[i0]));
    const oms = 1 - s;
    const Q = add(mul(oms, V[i0]), mul(s, V[i1]));
    result.distance = 0;
    result.sqrDistance = 0;
    result.parameter = dot(D, sub(Q, P)) / dot(D, D);
    result.barycentric[i0] = oms;
    result.barycentric[i1] = s;
    result.barycentric[i2] = 0;
    result.closest[0] = Q.clone();
    result.closest[1] = Q.clone();
}

// The triangle is strictly on one side of the line, so the closest triangle
// point is the vertex with the smallest absolute line-normal component.
function noCommonPoints(P: Vector, D: Vector, V: readonly Vector[],
    ncomp: readonly number[], result: DistLine2Triangle2Result): void {
    let minDistance = Math.abs(ncomp[0]);
    let minIndex = 0;
    let distance = Math.abs(ncomp[1]);
    if (distance < minDistance) {
        minDistance = distance;
        minIndex = 1;
    }
    distance = Math.abs(ncomp[2]);
    if (distance < minDistance) {
        minDistance = distance;
        minIndex = 2;
    }

    result.distance = minDistance;
    result.sqrDistance = minDistance * minDistance;
    result.parameter = dot(D, sub(V[minIndex], P)) / dot(D, D);
    result.barycentric[0] = (minIndex === 0 ? 1 : 0);
    result.barycentric[1] = (minIndex === 1 ? 1 : 0);
    result.barycentric[2] = (minIndex === 2 ? 1 : 0);
    result.closest[0] = add(P, mul(result.parameter, D));
    result.closest[1] = V[minIndex].clone();
}

export class DistLine2Triangle2
    implements DCPQuery<Line2, Triangle2, DistLine2Triangle2Result> {
    compute(line: Line2, triangle: Triangle2): DistLine2Triangle2Result {
        const result = defaultResult();

        // Test whether the triangle is strictly on one side of the line, in
        // which case the distance is the smallest absolute line-normal
        // component of the vertices. If at least one vertex is on the line,
        // then the distance is 0. If the triangle has at least one vertex
        // strictly on one side of the line and at least one vertex strictly
        // on the other side of the line, then the distance is 0.
        const P = line.origin;
        const D = line.direction;
        const V = triangle.v;
        const N = perp(D);
        const ncomp: number[] = [0, 0, 0];
        const sign: number[] = [0, 0, 0];
        for (let i = 0; i < 3; ++i) {
            ncomp[i] = dot(N, sub(V[i], P));
            if (ncomp[i] > 0) {
                sign[i] = +1;
            }
            else if (ncomp[i] < 0) {
                sign[i] = -1;
            }
            else {  // ncomp[i] = 0
                sign[i] = 0;
            }
        }

        // In the ensuing blocks of code, s0s1s2 represents the signs of the
        // normal component (ncomp) at the vertices. The term s? is in
        // {+,0,-}.
        if (sign[0] > 0) {
            if (sign[1] > 0) {
                if (sign[2] > 0) {
                    // +++ The triangle is strictly on the positive side of
                    // the line.
                    noCommonPoints(P, D, V, ncomp, result);
                }
                else if (sign[2] < 0) {
                    // ++- The line intersects triangle edges <V[2],V[0]> and
                    // <V[2],V[1]> at interior edge points.
                    lineIntersectsTwoEdges(P, D, V, 2, 0, 1, result);
                }
                else {
                    // ++0 The line intersects the triangle vertex V[2].
                    lineContainsVertex(P, D, V, 2, 0, 1, result);
                }
            }
            else if (sign[1] < 0) {
                if (sign[2] > 0) {
                    // +-+ The line intersects triangle edges <V[0],V[1]> and
                    // <V[2],V[1]> at interior edge points.
                    lineIntersectsTwoEdges(P, D, V, 0, 1, 2, result);
                }
                else if (sign[2] < 0) {
                    // +-- The line intersects triangle edges <V[0],V[1]> and
                    // <V[0],V[2]> at interior edge points.
                    lineIntersectsTwoEdges(P, D, V, 0, 1, 2, result);
                }
                else {
                    // +-0 The line intersects triangle edge <V[0],V[1]> at an
                    // interior point and at triangle vertex V[2].
                    lineContainsVertex(P, D, V, 2, 0, 1, result);
                }
            }
            else {
                if (sign[2] > 0) {
                    // +0+ The line intersects the triangle vertex V[1].
                    lineContainsVertex(P, D, V, 1, 2, 0, result);
                }
                else if (sign[2] < 0) {
                    // +0- The line intersects triangle edge <V[0],V[2]> at an
                    // interior point and at triangle vertex V[1].
                    lineContainsVertex(P, D, V, 1, 2, 0, result);
                }
                else {
                    // +00 The line contains triangle edge <V[1],V[2]>.
                    lineContainsVertex(P, D, V, 1, 2, 0, result);
                }
            }
        }
        else if (sign[0] < 0) {
            if (sign[1] > 0) {
                if (sign[2] > 0) {
                    // -++ The line intersects triangle edges <V[1],V[0]> and
                    // <V[2],V[0]> at interior edge points.
                    lineIntersectsTwoEdges(P, D, V, 0, 1, 2, result);
                }
                else if (sign[2] < 0) {
                    // -+- The line intersects triangle edges <V[1],V[0]> and
                    // <V[1],V[2]> at interior edge points.
                    lineIntersectsTwoEdges(P, D, V, 0, 1, 2, result);
                }
                else {
                    // -+0 The line intersects triangle edge <V[0],V[1]> at an
                    // interior point and at triangle vertex V[2].
                    lineContainsVertex(P, D, V, 2, 0, 1, result);
                }
            }
            else if (sign[1] < 0) {
                if (sign[2] > 0) {
                    // --+ The line intersects triangle edges <V[2],V[0]> and
                    // <V[2],V[1]> at interior edge points.
                    lineIntersectsTwoEdges(P, D, V, 1, 2, 0, result);
                }
                else if (sign[2] < 0) {
                    // --- The triangle is strictly on the negative side of
                    // the line.
                    noCommonPoints(P, D, V, ncomp, result);
                }
                else {
                    // --0 The line intersects the triangle vertex V[2].
                    lineContainsVertex(P, D, V, 2, 0, 1, result);
                }
            }
            else {
                if (sign[2] > 0) {
                    // -0+ The line intersects triangle edge <V[0],V[2]> at an
                    // interior point and at triangle vertex V[1].
                    lineContainsVertex(P, D, V, 1, 2, 0, result);
                }
                else if (sign[2] < 0) {
                    // -0- The line intersects the triangle vertex V[1].
                    lineContainsVertex(P, D, V, 1, 2, 0, result);
                }
                else {
                    // -00 The line contains triangle edge <V[1],V[2]>.
                    lineContainsVertex(P, D, V, 1, 2, 0, result);
                }
            }
        }
        else {
            if (sign[1] > 0) {
                if (sign[2] > 0) {
                    // 0++ The line intersects the triangle vertex V[0].
                    lineContainsVertex(P, D, V, 0, 1, 2, result);
                }
                else if (sign[2] < 0) {
                    // 0+- The line intersects triangle edge <V[1],V[2]> at an
                    // interior point and at triangle vertex V[0].
                    lineContainsVertex(P, D, V, 0, 1, 2, result);
                }
                else {
                    // 0+0 The line contains triangle edge <V[2],V[0]>.
                    lineContainsVertex(P, D, V, 2, 0, 1, result);
                }
            }
            else if (sign[1] < 0) {
                if (sign[2] > 0) {
                    // 0-+ The line intersects triangle edge <V[1],V[2]> at an
                    // interior point and at triangle vertex V[0].
                    lineContainsVertex(P, D, V, 0, 1, 2, result);
                }
                else if (sign[2] < 0) {
                    // 0-- The line intersects the triangle vertex V[0].
                    lineContainsVertex(P, D, V, 0, 1, 2, result);
                }
                else {
                    // 0-0 The line contains the triangle edge <V[2],V[0]>.
                    lineContainsVertex(P, D, V, 2, 0, 1, result);
                }
            }
            else {
                if (sign[2] > 0) {
                    // 00+ The line contains triangle edge <V[0],V[1]>.
                    lineContainsVertex(P, D, V, 0, 1, 2, result);
                }
                else if (sign[2] < 0) {
                    // 00- The line contains triangle edge <V[0],V[1]>.
                    lineContainsVertex(P, D, V, 0, 1, 2, result);
                }
                else {
                    // 000 The triangle is degenerate, a single point.
                    lineContainsVertex(P, D, V, 0, 1, 2, result);
                }
            }
        }

        const diff = sub(result.closest[0], result.closest[1]);
        result.sqrDistance = dot(diff, diff);
        result.distance = Math.sqrt(result.sqrDistance);
        return result;
    }
}
