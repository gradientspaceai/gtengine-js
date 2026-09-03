// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrAlignedBox2OrientedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the box to be a solid.
//
// The test-intersection query uses the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The set of potential separating directions includes the 2 edge normals of
// box0 and the 2 edge normals of box1. The integer 'separating' identifies
// the axis that reported separation; there may be more than one but only one
// is reported. The value is 0 when box0.axis[0] separates, 1 when
// box0.axis[1] separates, 2 when box1.axis[0] separates, or 3 when
// box1.axis[1] separates.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream provides
// only a TIQuery specialization for this pair of primitives.

import { AlignedBox } from './AlignedBox.js';
import { OrientedBox } from './OrientedBox.js';
import { sub, dot } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrAlignedBox2OrientedBox2TI.test.
export interface IntrAlignedBox2OrientedBox2TIResult {
    intersect: boolean;

    // Valid only when 'intersect' is false: 0 or 1 for a box0 axis, 2 or 3
    // for a box1 axis.
    separating: number;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrAlignedBox2OrientedBox2TIResult {
    return { intersect: false, separating: 0 };
}

export class IntrAlignedBox2OrientedBox2TI implements
    TIQuery<AlignedBox, OrientedBox, IntrAlignedBox2OrientedBox2TIResult> {

    test(box0: AlignedBox, box1: OrientedBox): IntrAlignedBox2OrientedBox2TIResult {
        const result = defaultTIResult();

        // Get the centered form of the aligned box. The axes are implicitly
        // A0[0] = (1,0) and A0[1] = (0,1).
        const { center: C0, extent: E0 } = box0.getCenteredForm();

        // Convenience variables.
        const C1 = box1.center;
        const A1 = box1.axis;
        const E1 = box1.extent;

        // Compute the difference of box centers.
        const D = sub(C1, C0);

        const absDot01: number[][] = [[0, 0], [0, 0]];
        let rSum = 0;

        // Test box0.axis[0] = (1,0).
        absDot01[0][0] = Math.abs(A1[0].values[0]);
        absDot01[0][1] = Math.abs(A1[1].values[0]);
        rSum = E0.values[0] + E1.values[0] * absDot01[0][0] +
            E1.values[1] * absDot01[0][1];
        if (Math.abs(D.values[0]) > rSum) {
            result.intersect = false;
            result.separating = 0;
            return result;
        }

        // Test axis box0.axis[1] = (0,1).
        absDot01[1][0] = Math.abs(A1[0].values[1]);
        absDot01[1][1] = Math.abs(A1[1].values[1]);
        rSum = E0.values[1] + E1.values[0] * absDot01[1][0] +
            E1.values[1] * absDot01[1][1];
        if (Math.abs(D.values[1]) > rSum) {
            result.intersect = false;
            result.separating = 1;
            return result;
        }

        // Test axis box1.axis[0].
        rSum = E1.values[0] + E0.values[0] * absDot01[0][0] +
            E0.values[1] * absDot01[1][0];
        if (Math.abs(dot(A1[0], D)) > rSum) {
            result.intersect = false;
            result.separating = 2;
            return result;
        }

        // Test axis box1.axis[1].
        rSum = E1.values[1] + E0.values[0] * absDot01[0][1] +
            E0.values[1] * absDot01[1][1];
        if (Math.abs(dot(A1[1], D)) > rSum) {
            result.intersect = false;
            result.separating = 3;
            return result;
        }

        result.intersect = true;
        return result;
    }
}
