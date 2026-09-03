// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrOrientedBox2Circle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The find-intersection query is based on the document
// https://www.geometrictools.com/Documentation/IntersectionMovingCircleRectangle.pdf
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The FIQuery is a
// dynamic (moving-objects) query whose 'operator()' takes four arguments, so
// it does not implement the two-argument FIQuery interface; the method is
// named 'find' with the four arguments, as in IntrAlignedBox2Circle2.ts.
//
// Upstream derives the FIQuery from FIQuery<AlignedBox2,Circle2> so that it
// can call that class's protected DoQuery after transforming the circle into
// box coordinates. In the port, the aligned-box DoQuery is a module-private
// function of IntrAlignedBox2Circle2.ts, so this file instead calls the
// public IntrAlignedBox2Circle2FI.find with an axis-aligned box centered at
// the origin (extents equal to the oriented-box extents) and a zero box
// velocity. That routine performs exactly the same steps upstream performs
// inline here (translate the circle into box coordinates, reflect into the
// first quadrant, DoQuery, undo the reflection, translate back by the zero
// box center), so the numerical results are identical; only the final
// rotation into world coordinates is applied here.

import { AlignedBox } from './AlignedBox.js';
import { DistPointOrientedBox } from './DistPointOrientedBox.js';
import { Hypersphere } from './Hypersphere.js';
import {
    IntrAlignedBox2Circle2FI,
    IntrAlignedBox2Circle2FIResultType
} from './IntrAlignedBox2Circle2.js';
import type { IntrAlignedBox2Circle2FIResult } from './IntrAlignedBox2Circle2.js';
import { logAssert } from './Logger.js';
import type { OrientedBox2 } from './OrientedBox.js';
import { Vector, add, dot, mul, negate, sub } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrOrientedBox2Circle2TI.test.
export interface IntrOrientedBox2Circle2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrOrientedBox2Circle2TIResult():
    IntrOrientedBox2Circle2TIResult {
    return { intersect: false };
}

// Upstream's derived query reuses the aligned-box Result verbatim.
export type IntrOrientedBox2Circle2FIResult = IntrAlignedBox2Circle2FIResult;

// The intersection query considers the box and circle to be solids; that is,
// the circle object includes the region inside the circular boundary and the
// box object includes the region inside the rectangular boundary. If the
// circle object and box object overlap, the objects intersect.
export class IntrOrientedBox2Circle2TI implements
    TIQuery<OrientedBox2, Hypersphere, IntrOrientedBox2Circle2TIResult> {

    test(box: OrientedBox2, circle: Hypersphere):
        IntrOrientedBox2Circle2TIResult {
        logAssert(box.dimension === 2 && circle.dimension === 2,
            'IntrOrientedBox2Circle2TI: mismatched sizes.');
        const pbQuery = new DistPointOrientedBox();
        const pbResult = pbQuery.compute(circle.center, box);
        return {
            intersect: (pbResult.sqrDistance <= circle.radius * circle.radius)
        };
    }
}

// Currently, only a dynamic query is supported. A static query will need to
// compute the intersection set of (solid) box and circle.
export class IntrOrientedBox2Circle2FI {
    find(box: OrientedBox2, boxVelocity: Vector, circle: Hypersphere,
        circleVelocity: Vector): IntrOrientedBox2Circle2FIResult {
        logAssert(box.dimension === 2 && boxVelocity.size === 2
            && circle.dimension === 2 && circleVelocity.size === 2,
            'IntrOrientedBox2Circle2FI: mismatched sizes.');

        // Transform the oriented box to an axis-aligned box centered at the
        // origin and transform the circle accordingly. Compute the velocity
        // of the circle relative to the box.
        const cdiff = sub(circle.center, box.center);
        const vdiff = sub(circleVelocity, boxVelocity);
        const C = Vector.zero(2);
        const V = Vector.zero(2);
        for (let i = 0; i < 2; ++i) {
            C.values[i] = dot(cdiff, box.axis[i]);
            V.values[i] = dot(vdiff, box.axis[i]);
        }

        const alignedBox = AlignedBox.fromMinMax(negate(box.extent),
            box.extent);
        const alignedCircle = Hypersphere.fromCenterRadius(C,
            circle.radius);
        const bcQuery = new IntrAlignedBox2Circle2FI();
        const result = bcQuery.find(alignedBox, Vector.zero(2),
            alignedCircle, V);

        if (result.intersectionType !== IntrAlignedBox2Circle2FIResultType.noContact) {
            // Transform back to the original coordinate system.
            result.contactPoint = add(box.center,
                add(mul(result.contactPoint.values[0], box.axis[0]),
                    mul(result.contactPoint.values[1], box.axis[1])));
        }
        return result;
    }
}
