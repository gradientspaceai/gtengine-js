// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrOrientedBox2Cone2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the box and cone to be solids.
//
// Define V = cone.ray.origin, D = cone.ray.direction, and cs = cone.cosAngle.
// Define C = box.center, U0 = box.axis[0], U1 = box.axis[1],
// e0 = box.extent[0], and e1 = box.extent[1]. A box point is
// P = C + x*U0 + y*U1 where |x| <= e0 and |y| <= e1. Define the function
//   F(P) = Dot(D, (P-V)/Length(P-V)) = F(x,y)
//        = (a0*x + a1*y + a2)/(x^2 + y^2 + 2*b0*x + 2*b1*y + b2)^{1/2}
// The function has an essential singularity when P = V. The box intersects
// the cone (with positive-area overlap) when at least one of the four box
// corners is strictly inside the cone. It is necessary that the numerator of
// F(P) be positive at such a corner. The interior of the solid cone is
// defined by the quadratic inequality
//   (Dot(D,P-V))^2 > |P-V|^2*(cone.cosAngle)^2
// This inequality is inexpensive to compute. In summary, overlap occurs when
// there is a box corner P for which
//   F(P) > 0 and (Dot(D,P-V))^2 > |P-V|^2*(cone.cosAngle)^2
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, which becomes IntrOrientedBox2Cone2TI.

import type { Cone } from './Cone.js';
import { IntrRay2OrientedBox2TI } from './IntrRay2OrientedBox2.js';
import { logAssert } from './Logger.js';
import type { OrientedBox2 } from './OrientedBox.js';
import { dot, sub } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrOrientedBox2Cone2TI.test.
export interface IntrOrientedBox2Cone2TIResult {
    // The value of 'intersect' is true when there is a box point that is
    // strictly inside the cone. If the box just touches the cone from the
    // outside, an intersection is not reported, which supports the common
    // operation of culling objects outside a cone.
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrOrientedBox2Cone2TIResult():
    IntrOrientedBox2Cone2TIResult {
    return { intersect: false };
}

// Test-intersection query for a solid oriented box and a solid cone in 2D.
export class IntrOrientedBox2Cone2TI implements
    TIQuery<OrientedBox2, Cone, IntrOrientedBox2Cone2TIResult> {

    test(box: OrientedBox2, cone: Cone): IntrOrientedBox2Cone2TIResult {
        logAssert(box.dimension === 2 && cone.dimension === 2,
            'IntrOrientedBox2Cone2TI: mismatched sizes.');
        const result = defaultIntrOrientedBox2Cone2TIResult();

        const rbQuery = new IntrRay2OrientedBox2TI();
        const rbResult = rbQuery.test(cone.ray, box);
        if (rbResult.intersect) {
            // The cone axis ray intersects the box, so the cone intersects
            // the box.
            result.intersect = true;
            return result;
        }

        const diff = sub(box.center, cone.ray.origin);
        const a0 = dot(cone.ray.direction, box.axis[0]);
        const a1 = dot(cone.ray.direction, box.axis[1]);
        const a2 = dot(cone.ray.direction, diff);
        const b0 = dot(box.axis[0], diff);
        const b1 = dot(box.axis[1], diff);
        const b2 = dot(diff, diff);
        const csSqr = cone.cosAngle * cone.cosAngle;

        for (let i1 = 0; i1 < 2; ++i1) {
            const sign1 = i1 * 2 - 1;
            const y = sign1 * box.extent.values[1];
            for (let i0 = 0; i0 < 2; ++i0) {
                const sign0 = i0 * 2 - 1;
                const x = sign0 * box.extent.values[0];
                const fNumerator = a0 * x + a1 * y + a2;
                if (fNumerator > 0) {
                    const dSqr = x * x + y * y + (b0 * x + b1 * y) * 2 + b2;
                    const nSqr = fNumerator * fNumerator;
                    if (nSqr > dSqr * csSqr) {
                        result.intersect = true;
                        return result;
                    }
                }
            }
        }

        result.intersect = false;
        return result;
    }
}
