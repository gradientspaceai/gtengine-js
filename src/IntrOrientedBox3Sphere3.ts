// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrOrientedBox3Sphere3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The find-intersection query is based on the document
// https://www.geometrictools.com/Documentation/IntersectionMovingSphereBox.pdf
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The FIQuery is a
// dynamic (moving-objects) query whose 'operator()' takes four arguments, so
// it does not implement the two-argument FIQuery interface; the method is
// named 'find' with the four arguments, as in IntrAlignedBox3Sphere3.ts.
//
// Upstream derives the FIQuery from FIQuery<AlignedBox3,Sphere3> so that it
// can call that class's protected DoQuery after transforming the sphere into
// box coordinates. In the port the aligned-box DoQuery is a module-private
// function of IntrAlignedBox3Sphere3.ts, so this file instead calls the
// public IntrAlignedBox3Sphere3FI.find with an axis-aligned box centered at
// the origin (extents equal to the oriented-box extents) and a zero box
// velocity. That routine performs the same DoQuery on the same inputs, and
// its extra separating-axis early-out only short-circuits configurations for
// which DoQuery reports no contact, so the results are identical; only the
// final rotation into world coordinates is applied here.

import { AlignedBox } from './AlignedBox';
import { DistPointOrientedBox } from './DistPointOrientedBox';
import { Hypersphere } from './Hypersphere';
import {
    IntrAlignedBox3Sphere3FI,
    IntrAlignedBox3Sphere3FIResultType
} from './IntrAlignedBox3Sphere3';
import type { IntrAlignedBox3Sphere3FIResult } from './IntrAlignedBox3Sphere3';
import { logAssert } from './Logger';
import type { OrientedBox3 } from './OrientedBox';
import { Vector, add, dot, mul, negate, sub } from './Vector';
import type { TIQuery } from './TIQuery';

// The result of IntrOrientedBox3Sphere3TI.test.
export interface IntrOrientedBox3Sphere3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrOrientedBox3Sphere3TIResult():
    IntrOrientedBox3Sphere3TIResult {
    return { intersect: false };
}

// The upstream derived FIQuery::Result adds no members.
export type IntrOrientedBox3Sphere3FIResult = IntrAlignedBox3Sphere3FIResult;

// The intersection query considers the box and sphere to be solids. For
// example, if the sphere is strictly inside the box (does not touch the box
// faces), the objects intersect.
export class IntrOrientedBox3Sphere3TI implements
    TIQuery<OrientedBox3, Hypersphere, IntrOrientedBox3Sphere3TIResult> {

    test(box: OrientedBox3, sphere: Hypersphere):
        IntrOrientedBox3Sphere3TIResult {
        logAssert(box.dimension === 3 && sphere.dimension === 3,
            'IntrOrientedBox3Sphere3TI: mismatched sizes.');
        const pbQuery = new DistPointOrientedBox();
        const pbResult = pbQuery.compute(sphere.center, box);
        return {
            intersect: (pbResult.sqrDistance <= sphere.radius * sphere.radius)
        };
    }
}

// Currently, only a dynamic query is supported. A static query will need to
// compute the intersection set of (solid) box and sphere.
export class IntrOrientedBox3Sphere3FI {
    find(box: OrientedBox3, boxVelocity: Vector, sphere: Hypersphere,
        sphereVelocity: Vector): IntrOrientedBox3Sphere3FIResult {
        logAssert(box.dimension === 3 && boxVelocity.size === 3
            && sphere.dimension === 3 && sphereVelocity.size === 3,
            'IntrOrientedBox3Sphere3FI: mismatched sizes.');

        // Transform the sphere and box so that the box center becomes the
        // origin and the box is axis aligned. Compute the velocity of the
        // sphere relative to the box.
        const cdiff = sub(sphere.center, box.center);
        const vdiff = sub(sphereVelocity, boxVelocity);
        const C = Vector.zero(3);
        const V = Vector.zero(3);
        for (let i = 0; i < 3; ++i) {
            C.values[i] = dot(cdiff, box.axis[i]);
            V.values[i] = dot(vdiff, box.axis[i]);
        }

        const alignedBox = AlignedBox.fromMinMax(negate(box.extent),
            box.extent);
        const alignedSphere = Hypersphere.fromCenterRadius(C, sphere.radius);
        const bsQuery = new IntrAlignedBox3Sphere3FI();
        const result = bsQuery.find(alignedBox, Vector.zero(3),
            alignedSphere, V);

        // Transform back to the original coordinate system.
        if (result.intersectionType !== IntrAlignedBox3Sphere3FIResultType.noContact) {
            const P = result.contactPoint;
            result.contactPoint = add(box.center,
                add(mul(P.values[0], box.axis[0]),
                    add(mul(P.values[1], box.axis[1]),
                        mul(P.values[2], box.axis[2]))));
        }
        return result;
    }
}
