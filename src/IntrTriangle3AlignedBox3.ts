// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrTriangle3AlignedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The test-intersection query is based on the document
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection query clips the triangle against the faces of the
// box.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The TI query
// transforms to a canonical box and delegates to
// IntrTriangle3CanonicalBox3TI. The FI query clips with
// IntrConvexPolygonHyperplaneFI; the upstream 'PPQuery'/'PPResult' type
// aliases are unnecessary in the port.

import type { AlignedBox3 } from './AlignedBox';
import { CanonicalBox } from './CanonicalBox';
import type { FIQuery } from './FIQuery';
import { Hyperplane } from './Hyperplane';
import type { Plane3 } from './Hyperplane';
import {
    IntrConvexPolygonHyperplaneConfiguration, IntrConvexPolygonHyperplaneFI
} from './IntrConvexPolygonHyperplane';
import { IntrTriangle3CanonicalBox3TI } from './IntrTriangle3CanonicalBox3';
import { logAssert } from './Logger';
import { Triangle } from './Triangle';
import type { Triangle3 } from './Triangle';
import type { TIQuery } from './TIQuery';
import { Vector, add, dot, mul, negate, sub } from './Vector';

// The result of IntrTriangle3AlignedBox3TI.test.
export interface IntrTriangle3AlignedBox3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrTriangle3AlignedBox3TIResult():
    IntrTriangle3AlignedBox3TIResult {
    return { intersect: false };
}

// The result of IntrTriangle3AlignedBox3FI.find. The 'insidePolygon' is the
// portion of the triangle inside the box (empty when the triangle does not
// intersect the box). The 'outsidePolygons' are the portions of the triangle
// clipped away by the box faces.
export interface IntrTriangle3AlignedBox3FIResult {
    insidePolygon: Vector[];
    outsidePolygons: Vector[][];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrTriangle3AlignedBox3FIResult():
    IntrTriangle3AlignedBox3FIResult {
    return { insidePolygon: [], outsidePolygons: [] };
}

// Create the six planes for the box faces, with normals that point inside
// the box. This is shared by the aligned-box and oriented-box FI queries via
// the caller-supplied center, axes and extent.
export function intrTriangle3BoxFacePlanes(center: Vector,
    axis: readonly Vector[], extent: Vector): Plane3[] {
    const planes: Plane3[] = [];
    for (let i = 0; i < 3; ++i) {
        const normal = negate(axis[i]);
        planes.push(Hyperplane.fromNormalConstant(normal,
            dot(normal, center) - extent.values[i]));
    }
    for (let i = 0; i < 3; ++i) {
        const normal = axis[i].clone();
        planes.push(Hyperplane.fromNormalConstant(normal,
            dot(normal, center) - extent.values[i]));
    }
    return planes;
}

// Clip the triangle against the box-face planes. This is the body shared by
// the upstream aligned-box and oriented-box FIQuery::operator() functions.
export function intrTriangle3BoxClip(triangle: Triangle3,
    planes: readonly Plane3[]): IntrTriangle3AlignedBox3FIResult {
    const C = IntrConvexPolygonHyperplaneConfiguration;
    const result = defaultIntrTriangle3AlignedBox3FIResult();

    // Start with the triangle and clip it against each face of the box. The
    // largest number of vertices for the polygon of intersection is 7.
    result.insidePolygon = [
        triangle.v[0].clone(), triangle.v[1].clone(), triangle.v[2].clone()
    ];

    const ppQuery = new IntrConvexPolygonHyperplaneFI();
    for (const plane of planes) {
        const ppResult = ppQuery.find(result.insidePolygon, plane);
        switch (ppResult.configuration) {
            case C.SPLIT:
                result.insidePolygon = ppResult.positivePolygon;
                result.outsidePolygons.push(ppResult.negativePolygon);
                break;
            case C.POSITIVE_SIDE_VERTEX:
            case C.POSITIVE_SIDE_EDGE:
            case C.POSITIVE_SIDE_STRICT:
                // The result.insidePolygon is already
                // ppResult.positivePolygon, but to make it clear, assign it
                // here.
                result.insidePolygon = ppResult.positivePolygon;
                break;
            case C.NEGATIVE_SIDE_VERTEX:
            case C.NEGATIVE_SIDE_EDGE:
            case C.NEGATIVE_SIDE_STRICT:
                result.insidePolygon = [];
                result.outsidePolygons.push(ppResult.negativePolygon);
                return result;
            case C.CONTAINED:
                // A triangle coplanar with a box face will be processed as if
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

// Test-intersection query for a triangle and an axis-aligned box in 3D.
export class IntrTriangle3AlignedBox3TI implements
    TIQuery<Triangle3, AlignedBox3, IntrTriangle3AlignedBox3TIResult> {

    test(triangle: Triangle3, box: AlignedBox3):
        IntrTriangle3AlignedBox3TIResult {
        logAssert(triangle.dimension === 3 && box.dimension === 3,
            'IntrTriangle3AlignedBox3TI: mismatched sizes.');
        const result = defaultIntrTriangle3AlignedBox3TIResult();

        // Transform the aligned box to a canonical box. Transform the
        // vertices accordingly.
        const canonicalBox =
            CanonicalBox.fromExtent(mul(0.5, sub(box.max, box.min)));
        const alignedBoxCenter = mul(0.5, add(box.max, box.min));

        const transformedTriangle = new Triangle(3);
        for (let i = 0; i < 3; ++i) {
            transformedTriangle.v[i] = sub(triangle.v[i], alignedBoxCenter);
        }

        const query = new IntrTriangle3CanonicalBox3TI();
        result.intersect =
            query.test(transformedTriangle, canonicalBox).intersect;
        return result;
    }
}

// Find-intersection query for a triangle and an axis-aligned box in 3D.
export class IntrTriangle3AlignedBox3FI implements
    FIQuery<Triangle3, AlignedBox3, IntrTriangle3AlignedBox3FIResult> {

    find(triangle: Triangle3, box: AlignedBox3):
        IntrTriangle3AlignedBox3FIResult {
        logAssert(triangle.dimension === 3 && box.dimension === 3,
            'IntrTriangle3AlignedBox3FI: mismatched sizes.');

        const center = mul(0.5, add(box.max, box.min));
        const extent = mul(0.5, sub(box.max, box.min));
        const axis = [
            Vector.unit(3, 0), Vector.unit(3, 1), Vector.unit(3, 2)
        ];
        return intrTriangle3BoxClip(triangle,
            intrTriangle3BoxFacePlanes(center, axis, extent));
    }
}
