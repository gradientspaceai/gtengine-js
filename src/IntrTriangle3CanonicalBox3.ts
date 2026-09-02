// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrTriangle3CanonicalBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The test-intersection query is based on the document
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection query clips the triangle against the faces of the
// canonical box.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The upstream
// TIQuery and FIQuery specializations become IntrTriangle3CanonicalBox3TI and
// IntrTriangle3CanonicalBox3FI; the private SeparatedBy* helpers become
// module-private functions. The C++ output reference parameters of
// GetTriangleProjection become a returned object literal.

import { CanonicalBox } from './CanonicalBox';
import type { Triangle } from './Triangle';
import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';
import { Hyperplane } from './Hyperplane';
import {
    IntrConvexPolygonHyperplaneConfiguration,
    IntrConvexPolygonHyperplaneFI
} from './IntrConvexPolygonHyperplane';
import { Vector, dot, sub } from './Vector';
import { cross } from './Vector3';

// The result of IntrTriangle3CanonicalBox3TI.test.
export interface IntrTriangle3CanonicalBox3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrTriangle3CanonicalBox3TIResult():
    IntrTriangle3CanonicalBox3TIResult {
    return { intersect: false };
}

// The result of IntrTriangle3CanonicalBox3FI.find.
export interface IntrTriangle3CanonicalBox3FIResult {
    // The portion of the triangle inside the box, as a convex polygon. It is
    // empty when the triangle does not intersect the box.
    //
    // Upstream documentation bug (corrected here, no code change): upstream
    // states the polygon has at most 7 vertices. The plane of the triangle
    // already cuts the box in a polygon of up to 6 edges, and each of the 3
    // triangle edges can contribute another, so the bound is 9. Octagons do
    // occur; see the test file. Nothing is preallocated, so the stale bound
    // never corrupted results.
    insidePolygon: Vector[];

    // The portions of the triangle outside the box, one convex polygon per
    // box face that clipped the triangle.
    outsidePolygons: Vector[][];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrTriangle3CanonicalBox3FIResult():
    IntrTriangle3CanonicalBox3FIResult {
    return { insidePolygon: [], outsidePolygons: [] };
}

// Test the direction of the triangle normal for separation.
function separatedByTriangleNormal(triangle: Triangle, box: CanonicalBox,
    edge: readonly Vector[]): boolean {
    const direction = cross(edge[0], edge[1]);
    const d = dot(direction, triangle.v[0]);

    const radius =
        Math.abs(box.extent.values[0] * direction.values[0]) +
        Math.abs(box.extent.values[1] * direction.values[1]) +
        Math.abs(box.extent.values[2] * direction.values[2]);

    return !(Math.abs(d) <= radius);
}

// Test the directions of the box face normals for separation.
function separatedByBoxFaceNormal(triangle: Triangle,
    box: CanonicalBox): boolean {
    for (let i = 0; i < 3; ++i) {
        const dot0 = triangle.v[0].values[i];
        const dot1 = triangle.v[1].values[i];
        const dot2 = triangle.v[2].values[i];
        let minDot = dot0, maxDot = dot0;
        if (dot1 < minDot) { minDot = dot1; } else if (dot1 > maxDot) { maxDot = dot1; }
        if (dot2 < minDot) { minDot = dot2; } else if (dot2 > maxDot) { maxDot = dot2; }
        const separated = (+box.extent.values[i] < minDot ||
            maxDot < -box.extent.values[i]);
        if (separated) {
            return true;
        }
    }
    return false;
}

// Compute the interval of projection of the triangle onto the specified
// direction.
function getTriangleProjection(direction: Vector, triangle: Triangle):
    { min: number, max: number } {
    const d: number[] = [
        dot(direction, triangle.v[0]),
        dot(direction, triangle.v[1]),
        dot(direction, triangle.v[2])
    ];

    let min = d[0];
    let max = min;

    if (d[1] < min) {
        min = d[1];
    } else if (d[1] > max) {
        max = d[1];
    }

    if (d[2] < min) {
        min = d[2];
    } else if (d[2] > max) {
        max = d[2];
    }

    return { min, max };
}

// Test the directions of the triangle-edge cross box-edge products for
// separation. The box edge directions are the standard Euclidean basis
// vectors, so the cross products are written out explicitly.
function separatedByTriangleEdgeCrossBoxEdge(triangle: Triangle,
    box: CanonicalBox, edge: readonly Vector[]): boolean {
    const e = [edge[0].values, edge[1].values, edge[2].values];
    const crossDirections: Vector[][] = [
        [
            Vector.fromArray([0, -e[0][2], +e[0][1]]),
            Vector.fromArray([0, -e[1][2], +e[1][1]]),
            Vector.fromArray([0, -e[2][2], +e[2][1]])
        ],
        [
            Vector.fromArray([+e[0][2], 0, -e[0][0]]),
            Vector.fromArray([+e[1][2], 0, -e[1][0]]),
            Vector.fromArray([+e[2][2], 0, -e[2][0]])
        ],
        [
            Vector.fromArray([-e[0][1], +e[0][0], 0]),
            Vector.fromArray([-e[1][1], +e[1][0], 0]),
            Vector.fromArray([-e[2][1], +e[2][0], 0])
        ]
    ];

    for (let i0 = 0; i0 < 3; ++i0) {
        for (let i1 = 0; i1 < 3; ++i1) {
            const direction = crossDirections[i0][i1];

            const { min, max } = getTriangleProjection(direction, triangle);

            const radius =
                Math.abs(box.extent.values[0] * direction.values[0]) +
                Math.abs(box.extent.values[1] * direction.values[1]) +
                Math.abs(box.extent.values[2] * direction.values[2]);

            const separated = (radius < min || max < -radius);
            if (separated) {
                return true;
            }
        }
    }
    return false;
}

// Test-intersection query for a triangle and a canonical box.
export class IntrTriangle3CanonicalBox3TI implements
    TIQuery<Triangle, CanonicalBox, IntrTriangle3CanonicalBox3TIResult> {

    test(triangle: Triangle, box: CanonicalBox):
        IntrTriangle3CanonicalBox3TIResult {
        const result = defaultIntrTriangle3CanonicalBox3TIResult();
        const edge: Vector[] = [new Vector(3), new Vector(3), new Vector(3)];

        // Test the direction of the triangle normal.
        edge[0] = sub(triangle.v[1], triangle.v[0]);
        edge[1] = sub(triangle.v[2], triangle.v[1]);
        if (separatedByTriangleNormal(triangle, box, edge)) {
            result.intersect = false;
            return result;
        }

        // Test the directions of the box faces.
        if (separatedByBoxFaceNormal(triangle, box)) {
            result.intersect = false;
            return result;
        }

        // Test the directions of the triangle-box edge cross products.
        edge[2] = sub(triangle.v[0], triangle.v[2]);
        result.intersect =
            !separatedByTriangleEdgeCrossBoxEdge(triangle, box, edge);
        return result;
    }
}

// Find-intersection query for a triangle and a canonical box.
export class IntrTriangle3CanonicalBox3FI implements
    FIQuery<Triangle, CanonicalBox, IntrTriangle3CanonicalBox3FIResult> {

    find(triangle: Triangle, box: CanonicalBox):
        IntrTriangle3CanonicalBox3FIResult {
        const result = defaultIntrTriangle3CanonicalBox3FIResult();
        const Cfg = IntrConvexPolygonHyperplaneConfiguration;

        // Start with the triangle and clip it against each face of the box.
        // (Upstream says the polygon has at most 7 vertices; see the Result
        // comments for why the bound is really 9.)
        result.insidePolygon = [
            triangle.v[0].clone(), triangle.v[1].clone(), triangle.v[2].clone()
        ];

        // Create planes for the box faces with normals that point inside the
        // box.
        const planes: Hyperplane[] = [
            Hyperplane.fromNormalConstant(negateUnit3(0), -box.extent.values[0]),
            Hyperplane.fromNormalConstant(negateUnit3(1), -box.extent.values[1]),
            Hyperplane.fromNormalConstant(negateUnit3(2), -box.extent.values[2]),
            Hyperplane.fromNormalConstant(Vector.unit(3, 0), -box.extent.values[0]),
            Hyperplane.fromNormalConstant(Vector.unit(3, 1), -box.extent.values[1]),
            Hyperplane.fromNormalConstant(Vector.unit(3, 2), -box.extent.values[2])
        ];

        const ppQuery = new IntrConvexPolygonHyperplaneFI();

        for (const plane of planes) {
            const ppResult = ppQuery.find(result.insidePolygon, plane);
            switch (ppResult.configuration) {
                case Cfg.SPLIT:
                    result.insidePolygon = ppResult.positivePolygon;
                    result.outsidePolygons.push(ppResult.negativePolygon);
                    break;
                case Cfg.POSITIVE_SIDE_VERTEX:
                case Cfg.POSITIVE_SIDE_EDGE:
                case Cfg.POSITIVE_SIDE_STRICT:
                    // The result.insidePolygon is already
                    // ppResult.positivePolygon, but to make it clear, assign
                    // it here.
                    result.insidePolygon = ppResult.positivePolygon;
                    break;
                case Cfg.NEGATIVE_SIDE_VERTEX:
                case Cfg.NEGATIVE_SIDE_EDGE:
                case Cfg.NEGATIVE_SIDE_STRICT:
                    result.insidePolygon = [];
                    result.outsidePolygons.push(ppResult.negativePolygon);
                    return result;
                case Cfg.CONTAINED:
                    // A triangle coplanar with a box face is processed as if
                    // it were inside the box.
                    result.insidePolygon = ppResult.intersection;
                    break;
                default:
                    result.insidePolygon = [];
                    result.outsidePolygons = [];
                    break;
            }
        }

        return result;
    }
}

// The negative of the i-th standard basis vector of 3-space.
function negateUnit3(i: number): Vector {
    const v = new Vector(3);
    v.values[i] = -1;
    return v;
}
