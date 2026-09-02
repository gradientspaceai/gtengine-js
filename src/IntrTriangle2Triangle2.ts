// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrTriangle2Triangle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The test-intersection query is based on the document
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection query for stationary triangles is based on clipping
// one triangle against the edges of the other to compute the intersection
// set (if it exists).
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The upstream
// TIQuery and FIQuery specializations become IntrTriangle2Triangle2TI and
// IntrTriangle2Triangle2FI. The protected WhichSide/Separated helpers become
// module-private functions; upstream declares them protected so that
// IntrTriangle3Triangle3 can reuse them, but the port has no such subclass
// requirement. Upstream also documents a moving-triangle find-intersection
// query in the file comments, but the header contains no such code; the
// remark is dropped.

import type { Triangle } from './Triangle';
import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';
import { Hyperplane } from './Hyperplane';
import { IntrConvexPolygonHyperplaneFI } from './IntrConvexPolygonHyperplane';
import { Vector, dot, sub } from './Vector';
import { perp } from './Vector2';

// The result of IntrTriangle2Triangle2TI.test.
export interface IntrTriangle2Triangle2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrTriangle2Triangle2TIResult():
    IntrTriangle2Triangle2TIResult {
    return { intersect: false };
}

// The result of IntrTriangle2Triangle2FI.find.
export interface IntrTriangle2Triangle2FIResult {
    // An intersection exists iff intersection.length > 0. The vertices are
    // counterclockwise ordered.
    intersection: Vector[];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrTriangle2Triangle2FIResult():
    IntrTriangle2Triangle2FIResult {
    return { intersection: [] };
}

// The triangle vertices are projected to t-values for the line P + t * D. The
// D-vector is nonzero but does not have to be unit length. The return value
// is +1 if all t >= 0, -1 if all t <= 0, but 0 otherwise, in which case the
// line splits the triangle into two subtriangles, each of positive area.
function whichSide(triangle: Triangle, P: Vector, D: Vector): number {
    // Upstream decrements 'negative' rather than incrementing it and relies
    // on the C++ truthiness of a nonzero (negative) int in the 'positive &&
    // negative' test. The port counts both up; the behavior is identical.
    let positive = 0, negative = 0;
    for (let i = 0; i < 3; ++i) {
        const t = dot(D, sub(triangle.v[i], P));
        if (t > 0) {
            ++positive;
        } else if (t < 0) {
            ++negative;
        }

        if (positive > 0 && negative > 0) {
            // The triangle has vertices strictly on both sides of the line,
            // so the line splits the triangle into two subtriangles each of
            // positive area.
            return 0;
        }
    }

    // Either positive > 0 or negative > 0 but not both are positive.
    return positive > 0 ? +1 : -1;
}

// Test the edges of triangle0 for separation. Because of the
// counterclockwise ordering, the projection interval for triangle0 is [T,0]
// for some T < 0. Determine whether triangle1 is on the positive side of the
// line; if it is, the triangles are separated.
function separated(triangle0: Triangle, triangle1: Triangle): boolean {
    for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
        // The potential separating axis is P + t * D.
        const P = triangle0.v[i0];
        const D = perp(sub(triangle0.v[i1], triangle0.v[i0]));
        if (whichSide(triangle1, P, D) > 0) {
            // The triangle1 projection interval is [a,b] where a > 0, so the
            // triangles are separated.
            return true;
        }
    }
    return false;
}

// Test whether two triangles intersect using the method of separating axes.
// The set of intersection, if it exists, is not computed. The input
// triangles' vertices must be counterclockwise ordered.
export class IntrTriangle2Triangle2TI implements
    TIQuery<Triangle, Triangle, IntrTriangle2Triangle2TIResult> {

    test(triangle0: Triangle, triangle1: Triangle):
        IntrTriangle2Triangle2TIResult {
        const result = defaultIntrTriangle2Triangle2TIResult();
        result.intersect =
            !separated(triangle0, triangle1) &&
            !separated(triangle1, triangle0);
        return result;
    }
}

// Find the convex polygon, segment or point of intersection of two
// triangles. The input triangles' vertices must be counterclockwise ordered.
export class IntrTriangle2Triangle2FI implements
    FIQuery<Triangle, Triangle, IntrTriangle2Triangle2FIResult> {

    find(triangle0: Triangle, triangle1: Triangle):
        IntrTriangle2Triangle2FIResult {
        const result = defaultIntrTriangle2Triangle2FIResult();

        // Start with triangle1 and clip against the edges of triangle0.
        let polygon: Vector[] = [
            triangle1.v[0].clone(), triangle1.v[1].clone(),
            triangle1.v[2].clone()
        ];

        const ppQuery = new IntrConvexPolygonHyperplaneFI();

        for (let i1 = 2, i0 = 0; i0 < 3; i1 = i0++) {
            // Create the clipping line for the current edge. The edge normal
            // N points inside the triangle.
            const P = triangle0.v[i0];
            const N = perp(sub(triangle0.v[i1], triangle0.v[i0]));
            const clippingLine = Hyperplane.fromNormalConstant(N, dot(N, P));

            // Do the clipping operation.
            const ppResult = ppQuery.find(polygon, clippingLine);
            if (ppResult.positivePolygon.length === 0) {
                // The current clipped polygon is outside triangle0.
                return result;
            }
            polygon = ppResult.positivePolygon;
        }

        result.intersection = polygon;
        return result;
    }
}
