// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrTriangle3OrientedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The test-intersection query is based on the document
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection query clips the triangle against the faces of the
// oriented box.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream
// duplicates the face-plane construction and the clipping loop of
// IntrTriangle3AlignedBox3.h verbatim; the port reuses the exported
// 'intrTriangle3BoxFacePlanes' and 'intrTriangle3BoxClip' helpers of
// IntrTriangle3AlignedBox3.ts instead (behavior-preserving).

import type { FIQuery } from './FIQuery';
import { CanonicalBox } from './CanonicalBox';
import {
    intrTriangle3BoxClip, intrTriangle3BoxFacePlanes
} from './IntrTriangle3AlignedBox3';
import { IntrTriangle3CanonicalBox3TI } from './IntrTriangle3CanonicalBox3';
import { logAssert } from './Logger';
import type { OrientedBox3 } from './OrientedBox';
import { Triangle } from './Triangle';
import type { Triangle3 } from './Triangle';
import type { TIQuery } from './TIQuery';
import { Vector, dot, sub } from './Vector';

// The result of IntrTriangle3OrientedBox3TI.test.
export interface IntrTriangle3OrientedBox3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrTriangle3OrientedBox3TIResult():
    IntrTriangle3OrientedBox3TIResult {
    return { intersect: false };
}

// The result of IntrTriangle3OrientedBox3FI.find. The 'insidePolygon' is the
// portion of the triangle inside the box (empty when the triangle does not
// intersect the box). The 'outsidePolygons' are the portions of the triangle
// clipped away by the box faces.
export interface IntrTriangle3OrientedBox3FIResult {
    insidePolygon: Vector[];
    outsidePolygons: Vector[][];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrTriangle3OrientedBox3FIResult():
    IntrTriangle3OrientedBox3FIResult {
    return { insidePolygon: [], outsidePolygons: [] };
}

// Test-intersection query for a triangle and an oriented box in 3D.
export class IntrTriangle3OrientedBox3TI implements
    TIQuery<Triangle3, OrientedBox3, IntrTriangle3OrientedBox3TIResult> {

    test(triangle: Triangle3, box: OrientedBox3):
        IntrTriangle3OrientedBox3TIResult {
        logAssert(triangle.dimension === 3 && box.dimension === 3,
            'IntrTriangle3OrientedBox3TI: mismatched sizes.');
        const result = defaultIntrTriangle3OrientedBox3TIResult();

        // Transform the oriented box to a canonical box. Transform the
        // triangle vertices accordingly.
        const canonicalBox = CanonicalBox.fromExtent(box.extent);

        const transformedTriangle = new Triangle(3);
        for (let j = 0; j < 3; ++j) {
            const diff = sub(triangle.v[j], box.center);
            const v = new Vector(3);
            for (let i = 0; i < 3; ++i) {
                v.values[i] = dot(box.axis[i], diff);
            }
            transformedTriangle.v[j] = v;
        }

        // Execute the test-intersection query.
        const tcQuery = new IntrTriangle3CanonicalBox3TI();
        result.intersect = tcQuery.test(transformedTriangle, canonicalBox)
            .intersect;
        return result;
    }
}

// Find-intersection query for a triangle and an oriented box in 3D.
export class IntrTriangle3OrientedBox3FI implements
    FIQuery<Triangle3, OrientedBox3, IntrTriangle3OrientedBox3FIResult> {

    find(triangle: Triangle3, box: OrientedBox3):
        IntrTriangle3OrientedBox3FIResult {
        logAssert(triangle.dimension === 3 && box.dimension === 3,
            'IntrTriangle3OrientedBox3FI: mismatched sizes.');

        // Create planes for the box faces with normals that point inside the
        // box, then clip the triangle against them.
        return intrTriangle3BoxClip(triangle,
            intrTriangle3BoxFacePlanes(box.center, box.axis, box.extent));
    }
}
