// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrAlignedBox3Cylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query is for finite cylinders. The cylinder and box are considered to
// be solids. The cylinder has center C, unit-length axis direction D, radius
// r and height h. The aligned box is converted to a canonical box after which
// a test-intersection query is performed on the finite cylinder and the
// canonical box. See the comments in IntrCanonicalBox3Cylinder3.ts for a
// brief description. The details are in
//   https://www.geometrictools.com/Documentation/IntersectionBoxCylinder.pdf
//
// Port notes: upstream has only a TIQuery, so this file exports only the TI
// class. The result has no new members relative to the canonical-box query,
// so IntrAlignedBox3Cylinder3TIResult is a type alias per the Intr*
// precedent.

import type { AlignedBox } from './AlignedBox';
import { CanonicalBox } from './CanonicalBox';
import { Cylinder3 } from './Cylinder3';
import {
    IntrCanonicalBox3Cylinder3TI
} from './IntrCanonicalBox3Cylinder3';
import type {
    IntrCanonicalBox3Cylinder3TIResult
} from './IntrCanonicalBox3Cylinder3';
import { logAssert } from './Logger';
import { sub } from './Vector';
import type { TIQuery } from './TIQuery';

// The result of IntrAlignedBox3Cylinder3TI.test. Upstream adds no members to
// the canonical-box-versus-cylinder result.
export type IntrAlignedBox3Cylinder3TIResult =
    IntrCanonicalBox3Cylinder3TIResult;

// Test-intersection query for a solid aligned box and a solid finite cylinder
// in 3D.
export class IntrAlignedBox3Cylinder3TI implements
    TIQuery<AlignedBox, Cylinder3, IntrAlignedBox3Cylinder3TIResult> {

    private readonly bcQuery = new IntrCanonicalBox3Cylinder3TI();

    test(box: AlignedBox, cylinder: Cylinder3):
        IntrAlignedBox3Cylinder3TIResult {
        logAssert(box.dimension === 3,
            'IntrAlignedBox3Cylinder3TI: mismatched sizes.');
        logAssert(cylinder.isFinite(),
            'Infinite cylinders are not yet supported.');

        // Convert the problem to one involving a finite cylinder and a
        // canonical box. This involves translating the box center to the
        // origin. The cylinder center must also be translated.
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();
        const cbox = CanonicalBox.fromExtent(boxExtent);
        const translatedCylinder = cylinder.clone();
        translatedCylinder.axis.origin =
            sub(translatedCylinder.axis.origin, boxCenter);

        const bcResult = this.bcQuery.test(cbox, translatedCylinder);
        return { intersect: bcResult.intersect };
    }
}
