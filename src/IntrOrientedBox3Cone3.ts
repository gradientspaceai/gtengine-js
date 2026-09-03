// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrOrientedBox3Cone3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test for intersection of a box and a cone. The cone can be infinite
//   0 <= minHeight < maxHeight = infinity
// or finite (cone frustum)
//   0 <= minHeight < maxHeight < infinity.
// The algorithm is described in
// https://www.geometrictools.com/Documentation/IntersectionBoxCone.pdf
// and reports an intersection only when the intersection set has positive
// volume. For example, let the box be outside the cone. If the box is below
// the minHeight plane at the cone vertex and just touches the cone vertex, no
// intersection is reported. If the box is above the maxHeight plane and just
// touches the disk capping the cone, either at a single point, a line segment
// of points or a polygon of points, no intersection is reported.
//
// See IntrAlignedBox3Cone3.ts for the upstream TODO about infinite cones.
//
// Port notes: upstream derives TIQuery<Real, OrientedBox3, Cone3> from
// TIQuery<Real, AlignedBox3, Cone3> and adds no result members. The port
// aliases IntrAlignedBox3Cone3TIResult and holds an IntrAlignedBox3Cone3TI
// instance (composition instead of inheritance), which mirrors the upstream
// call 'TIQuery<Real, AlignedBox3, Cone3>::operator()'.

import { AlignedBox } from './AlignedBox.js';
import type { Cone3 } from './Cone.js';
import {
    IntrAlignedBox3Cone3TI, defaultIntrAlignedBox3Cone3TIResult
} from './IntrAlignedBox3Cone3.js';
import type { IntrAlignedBox3Cone3TIResult } from './IntrAlignedBox3Cone3.js';
import { logAssert } from './Logger.js';
import type { OrientedBox3 } from './OrientedBox.js';
import type { TIQuery } from './TIQuery.js';
import { Vector, add, dot, sub } from './Vector.js';

// The result of IntrOrientedBox3Cone3TI.test. Upstream derives its Result
// from the aligned-box result and adds no members.
export type IntrOrientedBox3Cone3TIResult = IntrAlignedBox3Cone3TIResult;

// The port of the upstream Result default constructor.
export function defaultIntrOrientedBox3Cone3TIResult():
    IntrOrientedBox3Cone3TIResult {
    return defaultIntrAlignedBox3Cone3TIResult();
}

// Test-intersection query for an oriented box and a solid cone in 3D.
export class IntrOrientedBox3Cone3TI implements
    TIQuery<OrientedBox3, Cone3, IntrOrientedBox3Cone3TIResult> {

    private mBCQuery: IntrAlignedBox3Cone3TI;

    constructor() {
        this.mBCQuery = new IntrAlignedBox3Cone3TI();
    }

    test(box: OrientedBox3, cone: Cone3): IntrOrientedBox3Cone3TIResult {
        logAssert(box.dimension === 3 && cone.dimension === 3,
            'IntrOrientedBox3Cone3TI: mismatched sizes.');

        // Transform the cone and box so that the cone vertex is at the origin
        // and the box is axis aligned. This allows us to call the aligned-box
        // query.
        const diff = sub(box.center, cone.ray.origin);
        const xfrmBoxCenter = Vector.fromArray([
            dot(box.axis[0], diff),
            dot(box.axis[1], diff),
            dot(box.axis[2], diff)
        ]);
        const xfrmBox = new AlignedBox(3);
        xfrmBox.min = sub(xfrmBoxCenter, box.extent);
        xfrmBox.max = add(xfrmBoxCenter, box.extent);

        const xfrmCone = cone.clone();
        for (let i = 0; i < 3; ++i) {
            xfrmCone.ray.origin.values[i] = 0;
            xfrmCone.ray.direction.values[i] =
                dot(box.axis[i], cone.ray.direction);
        }

        // Test for intersection between the aligned box and the cone.
        const bcResult = this.mBCQuery.test(xfrmBox, xfrmCone);
        const result = defaultIntrOrientedBox3Cone3TIResult();
        result.intersect = bcResult.intersect;
        return result;
    }
}
