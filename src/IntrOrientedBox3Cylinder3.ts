// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrOrientedBox3Cylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query is for finite cylinders. The cylinder and box are considered to
// be solids. The cylinder has center C, unit-length axis direction D, radius
// r and height h. The oriented box is converted to a canonical box after
// which a test-intersection query is performed on the finite cylinder and the
// canonical box. See the comments in IntrCanonicalBox3Cylinder3.ts for a
// brief description. The details are in
//   https://www.geometrictools.com/Documentation/IntersectionBoxCylinder.pdf
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, so the port has only IntrOrientedBox3Cylinder3TI.

import { CanonicalBox } from './CanonicalBox.js';
import { Cylinder3 } from './Cylinder3.js';
import { IntrCanonicalBox3Cylinder3TI } from './IntrCanonicalBox3Cylinder3.js';
import { Line } from './Line.js';
import { logAssert } from './Logger.js';
import type { OrientedBox3 } from './OrientedBox.js';
import { Vector, dot, sub } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrOrientedBox3Cylinder3TI.test.
export interface IntrOrientedBox3Cylinder3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrOrientedBox3Cylinder3TIResult {
    return { intersect: false };
}

// Test-intersection query for a solid oriented box and a solid finite
// cylinder in 3D.
export class IntrOrientedBox3Cylinder3TI implements
    TIQuery<OrientedBox3, Cylinder3, IntrOrientedBox3Cylinder3TIResult> {

    test(box: OrientedBox3, cylinder: Cylinder3):
        IntrOrientedBox3Cylinder3TIResult {
        logAssert(box.dimension === 3,
            'IntrOrientedBox3Cylinder3TI: mismatched sizes.');
        logAssert(cylinder.isFinite(),
            'Infinite cylinders are not yet supported.');

        // Convert the problem to one involving a finite cylinder and a
        // canonical box. This involves translating the box center to the
        // origin and then rotating the box axes to the standard coordinate
        // axes. The cylinder center must also be translated and rotated
        // accordingly.
        const cbox = CanonicalBox.fromExtent(box.extent);
        const diff = sub(cylinder.axis.origin, box.center);
        const origin = Vector.zero(3);
        const direction = Vector.zero(3);
        for (let i = 0; i < 3; ++i) {
            origin.values[i] = dot(box.axis[i], diff);
            direction.values[i] = dot(box.axis[i], cylinder.axis.direction);
        }
        const transformedCylinder = Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(origin, direction), cylinder.radius,
            cylinder.height);

        const bcQuery = new IntrCanonicalBox3Cylinder3TI();
        const bcResult = bcQuery.test(cbox, transformedCylinder);
        const result = defaultTIResult();
        result.intersect = bcResult.intersect;
        return result;
    }
}
